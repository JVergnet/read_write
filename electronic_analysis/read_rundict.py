# readRun.py projectPath
# plots the projected DOS of all the vaspRun of the specified project
# project/job/vasprun.xml

import gzip
import json
import logging
import os
import shutil
import traceback
import warnings
from multiprocessing import Pool, cpu_count

from pymatgen.apps.borg.hive import (SimpleVaspToComputedEntryDrone,
                                     VaspToComputedEntryDrone)
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.io.vasp.inputs import Incar
# from pymatgen.io.vasp.outputs import Oszicar  # Outcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

import rw_utils.platform_id as platform_id

import structure_analysis.structure_geometry_utils as cluster

import electronic_analysis.read_mag_props as nupdown

global PARAM

# DICT of parameters to set to true or false.
PARAM = {}

PARAM['logFile'] = True  # wether to print a log file of the results
PARAM['mainFolder'] = platform_id.setting_dir()
PARAM['verbose'] = 0

# DATA COLLECTION FUNCTIONS
# =========================================


def get_job_list(init_dir, file_system_choice=None):
    print("\n {} \n ".format(init_dir))
    if file_system_choice is None:
        file_system_choice = input(
            "[j]ob / [p]roject / [s]uper_project ?  :  ")

    if file_system_choice[0] == "p":
        # file_system = "p"
        subdir_list = [
            os.path.join(init_dir, o)
            for o in os.listdir(init_dir)
            if os.path.isdir(os.path.join(init_dir, o))]

    elif file_system_choice[0] == "s":
        # file_system = "s"
        subdir_list = []
        for folder in [
                os.path.join(init_dir, o)
                for o in os.listdir(init_dir)
                if os.path.isdir(os.path.join(init_dir, o))]:
            subdir_list += [
                os.path.join(folder, o)
                for o in os.listdir(folder)
                if os.path.isdir(os.path.join(folder, o))]

    else:
        # file_system = "j"
        subdir_list = [init_dir]
    return(subdir_list, file_system_choice)


class Rundict(ComputedStructureEntry):
    """Describe a VASP run AFTER computation

        Holds VASP inputs & outputs
        Compute and holds post-processing results
        """

    status_dict = {
        3:    "post-run : converged",
        1.5:  "post-run : corrupted / suppressed vasprun",
        2:    "post-run : not converged",
        1:    "pre-run",
        0:    "not a job folder"}

    def __init__(self, c_e, status, job_folder):

        self.status = status
        self.status_string = self.status_dict[status]
        self.job_folder = job_folder
        self.stacking = job_folder.split('/')[-2]
        self.str_id = job_folder.split('/')[-1].split("__")[-1]
        if c_e is not None:
            energy = 1000000 if status <= 1 else c_e.energy
            ComputedStructureEntry.__init__(self,
                                            c_e.structure, energy,
                                            correction=c_e.correction,
                                            parameters=c_e.parameters,
                                            data=c_e.data,
                                            entry_id=None)
            # define all attributes
            self.nb_cell = 1
            self.x_na = 0
            self.volume = None
            self.formula = None
            # self.name_tag = None
            self.spacegroup = None
            self.equivSiteList = None
            self.dOO_min = None
            self.dOO_min_indices = None
            # self.mag = None
            self.OO_pairs = None
            self.MMOO_quadruplets = None
            self.bader_done = False
            # self.structure_data = self.structure
            self._mag = None

        self.generate_tags()

    @property
    def nelect(self):
        nelect = self.parameters["incar"].get('NELECT', None)
        if nelect is None:
            print("NELECT not defined, returning 0")
            return 0

        return nelect

    @property
    def name_tag(self):
        return "{} : {} : {}".format(
            self.formula, self.stacking, self.str_id)

    @property
    def mag(self):
        "lazy parsing of the oszicar"
        if self._mag is None:
            self._mag = nupdown.get_mag_single(self)
        return self._mag

    def generate_tags(self):
        " incrementally generates tags & attributes depending on the convergence "
        if self.status == 0:
            return "{} Not a job folder".format(self.job_folder)
        if self.status > 0:  # pre-run : just the structure
            self.get_structure_tag()
            self.get_nametag()
        if self.status > 1:  # unconverged / corrupted run : just the energy
            self.energy_per_fu = self.energy / self.nb_cell
            self.structure_data = self.structure.copy()
        if self.status >= 3:  # converged run : DOS & co.
            print(self.name_tag, self.data["efermi"])
            # self.complete_dos = self.data["complete_dos"]
        return "tags generated"

    def get_nametag(self):
        "get nice formula of the current structure"
        self.nb_cell = get_nb_cell(self.structure)

        # Normalizing composition to get Nax My O2
        self.x_na, self.formula = get_xna_and_formula(
            self.structure, self.nb_cell)

    def get_structure_tag(self):
        "Structure and composition related tags "

        analyzer = SpacegroupAnalyzer(self.structure, symprec=0.1)
        self.spacegroup = analyzer.get_space_group_symbol()

        # get equivalent sites of the current structure in a list of dict
        self.equivSiteList = get_equiv_site_list(self.structure)

        (self.dOO_min, self.dOO_min_indices) = \
            cluster.get_min_OO_dist(self.structure)

        self.volume = self.structure.lattice.volume / self.nb_cell

    def get_MMOO_tags(self):
        """
        get oxygens pairs for quantification of A.R. distortion
                ACHTUNG ! : Cannot define layers in Na rich compositions
        """

        try:
            self.MMOO_quadruplets = \
                cluster.get_MMOO_quadruplets(self.structure)
            self.OO_pairs = [q["oxygen_pair"]
                             for q in self.MMOO_quadruplets]
        except Exception as ex:
            print(traceback.format_exc())
            print("could no define MMOO clusters because {}".format(ex))

    def as_dict(self):
        dict_of_run = super().as_dict()
        dict_of_run["@module"] = self.__class__.__module__
        dict_of_run["@class"] = self.__class__.__name__
        list_of_attr = ['status', 'status_string',
                        'job_folder', 'stacking']
        # even valid for invalid folders
        if self.status > 0:  # pre-run : just the structure
            list_of_attr += ['str_id', 'nb_cell', 'x_na', 'volume', 'formula']
        if self.status > 1:  # unconverged / corrupted run : just the energy
            list_of_attr += ["energy", "energy_per_fu"]
        if self.status >= 3 and self.bader_done:  # converged run : anything ?
            list_of_attr.append("structure_data")

        dict_of_run.update({key: getattr(self, key) for key in list_of_attr})
        return dict_of_run

    # def from_dict(cls, d):
    #     dec = MontyDecoder()  # /!\ not implemented !!
    #     return cls(dec.process_decoded(d["structure"]),
    #                d["energy"], d["correction"],
    #                parameters={k: dec.process_decoded(v)
    #                            for k, v in d.get("parameters", {}).items()},
    #                data={k: dec.process_decoded(v)
    #                      for k, v in d.get("data", {}).items()},
    #                entry_id=d.get("entry_id", None))

    @classmethod
    def from_structure(cls, structure, entry_id=None):
        status = 1  # pre-run status
        job_folder = ""
        c_e = ComputedStructureEntry(structure=structure,
                                     energy=100000000,
                                     entry_id=entry_id)
        return cls(c_e, status, job_folder)


def if_file_exist_gz(folder, non_comp_file):
    cwd = os.getcwd()
    os.chdir(folder)
    gz_file = "{}.gz".format(non_comp_file)
    if os.path.exists(gz_file):
        print("unzipping {} in {}".format(
            gz_file, folder))
        with gzip.open(gz_file, 'rb') as f_in:
            with open(non_comp_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(gz_file)
    os.chdir(cwd)
    return(os.path.exists(os.path.join(folder, non_comp_file)))


def collect_single_folder(job_folder, drone=None, vasprun_parsing_lvl=1):
    """
    Reading and parsing of VASP files in the job_folder folder
    vasprun_parsing_lvl > 0.5 : parse vasprun, else only parse vasp inputs
    returns a Rundict instance
    """
    # status :
    # 0 : not a job folder
    # 1 : pre-run
    # 2 : failed run
    # 3 : converged run

    # if no drone give, minimal parsing
    if drone is None:
        drone = drone = VaspToComputedEntryDrone(
            inc_structure=True,
            parameters=["incar"],
            data=['efermi', "converged"])

    status = 0

    with warnings.catch_warnings():
        # warnings.filterwarnings('once')
        if PARAM['verbose'] == 0:
            warnings.simplefilter("ignore")
        if vasprun_parsing_lvl > 0.5:
            if if_file_exist_gz(job_folder, "vasprun.xml"):
                # try  unzip the vasprun
                c_e = drone.assimilate(job_folder)
                if c_e is not None:
                    print(c_e.parameters['incar'])
                    status = 3 if c_e.data['converged'] else 2

        if status <= 0:
            # try  unzip the whole dir
            c_e = SimpleVaspToComputedEntryDrone(
                inc_structure=True).assimilate(job_folder)
            if c_e is not None:
                incar = Incar.from_file(
                    '{}/INCAR'.format(job_folder)).as_dict()
                c_e.parameters.update({"incar": incar})
                try:
                    status = 2 if c_e.energy < 10000000 else 1
                except TypeError as ex:
                    print("exception caught", ex)
                    # c_e.energy = 10000000
                    status = 1  # if there is an error, the run didn't converge

            # default energy of parser when ozscicar not read
            # to distinguish "pre-runs" from corrupted / suppressed vasprun

        if status > 0:
            for k in ['LDAUU', 'LDAUJ', 'LDAUL',
                      "@class", "@module"]:  # "MAGMOM"
                c_e.parameters['incar'].pop(k, None)
            for file_name in ["param", "data"]:
                file_path = os.path.join(
                    job_folder, '{}.json'.format(file_name))
                try:
                    with open(file_path) as file_object:
                        # load file data into object
                        param = json.load(file_object)
                        print(file_name, param)
                        c_e.parameters['custom'] = param
                except Exception:
                    # print(ex)
                    continue
                else:
                    break
            # print("finished gathering c_e data")
            # if hasattr(c_e, "run_type"):
            #     print(c_e.run_type)
            # print(c_e.parameters)
            # print(type(c_e.parameters["run_type"]))
            if not isinstance(c_e.parameters["run_type"], str):
                c_e.parameters["run_type"] = None

            # print("run type changed to None")
            # print(c_e)

    rundict = Rundict(c_e, status, job_folder)

    if PARAM["verbose"] > 0:
        print("sucessfully parsed : \n", rundict)

    return(rundict)


def collect_valid_runs(mainFolder, checkDiff=False,
                       vasprun_parsing_lvl=1, file_system_choice=None):
    """
    input : mainFolder ,{ param : value}
    output : [ {vaspRun : v , folder : F } ]

    COLLECT VALID VASPRUNS =================================================
    Walk trhough the folders to find valid vaspruns
    Convergence is tested by pymatgen built in fct
    If the converged structure if different from all the others in the list,
    it is added to the vaspRunList[]
    Create a list of all the subfolder of the run
    """
    global PARAM

    sub_dir_list, file_system_choice = get_job_list(
        mainFolder, file_system_choice=file_system_choice)

    if vasprun_parsing_lvl is None:
        try:
            parse_choice = input(
                "check vasprun level ?\n" +
                "\n".join([
                    "[f]ull vasprun parsing",
                    "[m]inimal vasprun parsing",
                    "[o]szicar parsing (energy & structure only)",
                    "[n]o parsing"]))
            if parse_choice[0] in ["f", "F"]:  # [f]ull parsing
                vasprun_parsing_lvl = 1
            elif parse_choice[0] in ["m", "M"]:  # [m]inimal parsing
                vasprun_parsing_lvl = 0.7
            elif parse_choice[0] in ["o", "O"]:  # [m]inimal parsing
                vasprun_parsing_lvl = 0.4
            elif parse_choice[0] in ["n", "N"]:  # [n]o parsing
                vasprun_parsing_lvl = 0
        except Exception as ex:
            print("Exception : {}\n default to full parsing".format(ex))
            vasprun_parsing_lvl = 1

    # print(str(subDirList)+"\n\n")
    drone = VaspToComputedEntryDrone(
        inc_structure=True,
        parameters=["incar"],
        data=['efermi', "complete_dos", "tdos", "converged"])

    with Pool(processes=cpu_count()) as parallel_runs:
        tmp_list = parallel_runs.starmap(
            collect_single_folder, [
                (run, drone, vasprun_parsing_lvl) for run in sub_dir_list])
        parallel_runs.close()
        parallel_runs.join()

    tmp_list.sort(key=lambda x: x.job_folder)
    # for d in tmp_list:
    #     print("{} / {} : {}".format(d.stacking, d.str_id,
    #                                 d.status_string))

    converged_runs = []

    if checkDiff:
        for run in [run for run in tmp_list if run.status > 0]:
            new_struct = True
            for existing_run in converged_runs:
                if run.structure.matches(
                        existing_run.structure):
                    print(
                        "Structure : {0} match with previous one :  {1}".format(
                            run.parameters['incar'].get(
                                'SYSTEM', "no_name"),
                            existing_run.parameters['incar'].get(
                                'SYSTEM', "no_name")))
                    new_struct = False
            if new_struct and run.status > 0:
                converged_runs.append(run)

    else:
        valid_runs = [run for run in tmp_list if run.status > 0]

    # os.chdir(valid_runs)
    print("\nnb of valid runs : {0}".format(len(valid_runs)))
    for run in valid_runs:
        print(run.name_tag, run.status_string)
    return(valid_runs)


# def write_log_file(rundict_list, log_folder):

#     log_file_name = platform_id.get_file_name(PARAM['mainFolder'], "log")
#     os.chdir(log_folder)
#     # log_file = open(logFileName,"a")
#     if input("print runs name & energies ?[Y]") == "Y":
#         print_stdout = True

#     print("\nSaving logFile in  : {0} \n".format(log_file_name))
#     with open(log_file_name, 'w') as log_file:
#         for rundict in rundict_list:
#             message = make_rundict_str(rundict)
#             if print_stdout:
#                 print(message)
#             log_file.write(message)


# def make_rundict_str(rundict):
#     message = "\n".join([
#         "NAME",
#         rundict.name_tag,
#         rundict.folder,
#         "ENERGY",
#         rundict.energy_tag,
#         "STRUCTURE",
#         rundict.structure_tag,
#     ])

#     if getattr(rundict, 'bader_tag', None) is not None:
#         message += "\n".join(["", "BADER", rundict.baderTag])
#     return message


def get_equiv_site_list(structure):
    """
    print a list of the equivalent sites of a structure
    with one representative position (mainIndex)
    the list on all the indices , the element on the site and the
    multiplicity
    """

    anal = SpacegroupAnalyzer(structure, symprec=0.05)
    sym_struct = anal.get_symmetrized_structure()
    sym_list = sym_struct.equivalent_indices
    equivalent_sites = []
    for k in range(0, len(sym_list), 1):
        # ACHTUNG ! : Position in the structure start at 0
        equivalent_sites.append({
            'mainIndex': sym_list[k][0],
            'indices': sym_list[k],
            'element': structure[sym_list[k][0]].species_string,
            'multiplicity': len(sym_list[k]),
            'charge': -1,
            'magnetization': -1})
    return equivalent_sites


def initialize(working_dir):
    global PARAM

    # Check if we are on frodon

    PARAM["run_on_cluster"] = platform_id.first_check_cluster()

    # get the path to the "project folder"
    # (i.e. the common parent directory of all the jobs to plot)
    PARAM['mainFolder'] = working_dir
    print("\nMain directory : \n{0} \n".format(PARAM['mainFolder']))

    print("\n preparing the run with the following parameters :\n{} \n\n".
          format(PARAM))

    log_file_name = platform_id.get_file_name(PARAM['mainFolder'], "log")
    logging.basicConfig(filename=log_file_name, level=logging.DEBUG)
    logging.info('Started')

    parameters = PARAM
    return parameters


def get_nb_cell(structure):
    compo_dict = structure.composition.get_el_amt_dict()
    nb_cell = 1
    # sulfide or oxide ?
    for spec in ['S', 'O']:
        if compo_dict.get(spec, 0) > 0:
            nb_cell = compo_dict[spec] / 2
            break
    return nb_cell


def get_xna_and_formula(structure, nb_cell=1):
    x_na = 0
    # Normalizing composition to get Nax My O2
    compo_dict = structure.composition.get_el_amt_dict()
    for key in compo_dict.keys():
        compo_dict[key] = round(compo_dict[key] / nb_cell, 2)

    # Na or Li insertion ?
    for spec in ['Na', 'Li']:
        if compo_dict.get(spec, 0) > 0:
            x_na = compo_dict[spec]
            break
    formula = Composition(compo_dict).formula
    return(x_na, formula)
