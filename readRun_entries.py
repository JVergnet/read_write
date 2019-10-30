# readRun.py projectPath
# plots the projected DOS of all the vaspRun of the specified project
# project/job/vasprun.xml
# from pymatgen.io.vasp import Vasprun
# from zipfile import ZipFile
import gzip
import json
# from pymatgen.electronic_structure.core import Spin, OrbitalType
# import pymatgen.analysis.diffraction.xrd as xrd
import logging
# from matplotlib.collections import LineCollection
import os
import shutil
import sys
# import importlib
import traceback
import warnings
from itertools import repeat
# import subprocess
from multiprocessing import Pool, cpu_count

from pymatgen.apps.borg.hive import (SimpleVaspToComputedEntryDrone,
                                     VaspToComputedEntryDrone)
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.io.vasp.inputs import Incar
from pymatgen.io.vasp.outputs import Oszicar  # Outcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

import structure_geometry_utils as cluster
import platform_id
# import lobster_coop as lob
# import generic_plot as generic_plot
import read_hull as hull

# from pymatgen.core.structure import Structure

# import matplotlib

# import matplotlib.pyplot as plt
# import numpy as np
# from operator import itemgetter
# import copy

# import platform


# import bailar_twist as bailar
# import readBader as bader
# import readO2 as O2
# import DOS_plot as DOS
# import energy_surface as PES
# import nupdown_scan as nupdown
# import launchDisordered as launch

global PARAM

# DICT of parameters to set to true or false.
PARAM = {}

PARAM['logFile'] = True  # wether to print a log file of the results
PARAM['mainFolder'] = platform_id.setting_dir()
PARAM['verbose'] = 0

# DATA COLLECTION FUNCTIONS
# =========================================


def get_job_list(main_folder, file_system_choice=None):
    print("\n {} \n ".format(main_folder))
    if file_system_choice is None:
        file_system_choice = input(
            "[j]ob / [p]roject / [s]uper_project ?  :  ")

    if file_system_choice[0] == "p":
        # file_system = "p"
        d = main_folder
        subdir_list = [
            os.path.join(d, o)
            for o in os.listdir(d)
            if os.path.isdir(os.path.join(d, o))]

    elif file_system_choice[0] == "s":
        # file_system = "s"
        d = main_folder
        subdir_list = []
        for folder in [
                os.path.join(d, o)
                for o in os.listdir(d)
                if os.path.isdir(os.path.join(d, o))]:
            subdir_list += [
                os.path.join(folder, o)
                for o in os.listdir(folder)
                if os.path.isdir(os.path.join(folder, o))]

    else:
        # file_system = "j"
        subdir_list = [main_folder]
    return(subdir_list, file_system_choice)


class Rundict(ComputedStructureEntry):
    "base class for vasprun computation"

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
            ComputedStructureEntry.__init__(self,
                                            c_e.structure, c_e.energy,
                                            correction=c_e.correction,
                                            parameters=c_e.parameters,
                                            data=c_e.data,
                                            entry_id=None)
        # define all attributes
        self.nb_cell = 1
        self.x_na = 0
        self.volume = None
        self.formula = None
        self.name_tag = None
        self.spacegroup = None
        self.equivSiteList = None
        self.dOO_min = None
        self.dOO_min_indices = None
        self.mag = None
        self.OO_pairs = None
        self.MMOO_quadruplets = None
        self.bader_done = False

        self.generate_tags()
#     if c_e is not None:

    @property
    def nelect(self):
        return self.parameters["incar"]['NELECT']

    @property
    def name_tag(self):
        return "{} : {} : {}".format(
            self.formula, self.stacking, self.str_id)

    @property
    def mag(self):
        try:
            mag_tot = Oszicar(os.path.join(
                self.job_folder, "OSZICAR")).ionic_steps[-1]["mag"]
            return mag_tot / len(self.structure.indices_from_symbol("Mn"))
        except Exception:
            print("could not define mag for", self.name_tag)
            return None

    def generate_tags(self):
        " progressively generates tags & attributes depending on the convergence "
        if self.status == 0:
            return("{} Not a job folder".format(self.job_folder))
        if self.status > 0:  # pre-run : just the structure
            self.get_structure_tag()
            self.get_nametag()
        if self.status > 1:  # unconverged / corrupted run : just the energy
            self.energy_per_fu = self.energy / self.nb_cell
            self.structure_data = self.structure.copy()
        if self.status >= 3:  # converged run : DOS & co.
            print(self.name_tag, self.data["efermi"])
            # self.complete_dos = self.data["complete_dos"]
        return("tags generated")

    def get_nametag(self):
        " nameTag : get name of the current structure in a string nameTag"
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
    # def coord(self, coord):
    #     return({"x_na": self.x_na, "doo_min": self.doo_min}[coord])

    # runDict['bandgap'] = runDict['vaspRun'].eigenvalue_band_properties[0]

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
        d = super().as_dict()
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d.update(dict(status=self.status,
                      status_string=self.status_string,
                      job_folder=self.job_folder,
                      stacking=self.stacking,
                      id=self.str_id,
                      # structure_data=self.structure_data.as_dict()
                      # mag=self.mag
                      ))
        return d

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
        return(cls(c_e, status, job_folder))


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


def collect_single_folder(job_folder, drone, vaspRun_parsing_lvl=1):
    # status :
    # 0 : not a job folder
    # 1 : pre-run
    # 2 : failed run
    # 3 : converged run

    status = 0

    with warnings.catch_warnings():
        # warnings.filterwarnings('once')
        if PARAM['verbose'] == 0:
            warnings.simplefilter("ignore")
        if vaspRun_parsing_lvl > 0.5:
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
                status = 2 if c_e.energy < 10000000 else 1

            # default energy of parser when ozscicar not read
            # to distinguish "pre-runs"
            #   from corrupted / suppressed vasprun
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
                except Exception as ex:
                    # print(ex)
                    continue
                else:
                    break
            # print("finished gathering c_e data")
            # if hasattr(c_e, "run_type"):
            #     print("WOWOWOWOWOWOWWOWOWOWOWOWWOWOWOWO")
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


def collect_valid_runs(
        mainFolder,
        checkDiff=False,
        vaspRun_parsing_lvl=1,
        file_system_choice=None):

    global PARAM

    # input : mainFolder ,{ param : value}
    # output : [ {vaspRun : v , folder : F } ]

    # COLLECT VALID VASPRUNS =================================================
    # Walk trhough the folders to find valid vaspruns
    # Convergence is tested by pymatgen built in fct
    # If the converged structure if different from all the others in the list,
    # it is added to the vaspRunList[]
    # Create a list of all the subfolder of the run

    sub_dir_list, file_system_choice = get_job_list(
        mainFolder, file_system_choice=file_system_choice)

    if vaspRun_parsing_lvl is None:
        try:
            parse_choice = input(
                "check vasprun level ?\n" +
                "\n".join([
                    "[f]ull vasprun parsing",
                    "[m]inimal vasprun parsing",
                    "[o]szicar parsing (energy & structure only)",
                    "[n]o parsing"]))
            if parse_choice[0] in ["f", "F"]:  # [f]ull parsing
                vaspRun_parsing_lvl = 1
            elif parse_choice[0] in ["m", "M"]:  # [m]inimal parsing
                vaspRun_parsing_lvl = 0.7
            elif parse_choice[0] in ["o", "O"]:  # [m]inimal parsing
                vaspRun_parsing_lvl = 0.4
            elif parse_choice[0] in ["n", "N"]:  # [n]o parsing
                vaspRun_parsing_lvl = 0
        except Exception as ex:
            print("Exception : {}\n default to full parsing".format(ex))
            vaspRun_parsing_lvl = 1

    # print(str(subDirList)+"\n\n")
    drone = VaspToComputedEntryDrone(
        inc_structure=True,
        parameters=["incar"],
        data=['efermi', "complete_dos", "tdos", "converged"])

    with Pool(processes=cpu_count()) as p:
        tmp_list = p.starmap(
            collect_single_folder, [
                (d, drone, vaspRun_parsing_lvl) for d in sub_dir_list])
        p.close()
        p.join()

    tmp_list.sort(key=lambda x: x.job_folder)
    # for d in tmp_list:
    #     print("{} / {} : {}".format(d.stacking, d.id,
    #                                 d.status_string))

    converged_runs = []

    if checkDiff:
        for d in [d for d in tmp_list if d.status > 0]:
            new_struct = True
            for existing_run in converged_runs:
                if d.structure.matches(
                        existing_run.structure):
                    print(
                        "Structure : {0} match with previous one :  {1}".format(
                            d.parameters['incar'].get('SYSTEM', "no_name"),
                            existing_run.parameters['incar'].get(
                                'SYSTEM', "no_name")))
                    new_struct = False
            if new_struct and d.status > 0:
                converged_runs.append(d)

    else:
        valid_runs = [d for d in tmp_list if d.status > 0]

    # os.chdir(valid_runs)
    print("\nnb of valid runs : {0}".format(len(valid_runs)))
    for run in valid_runs:
        print(run.name_tag, run.status_string)
    return(valid_runs)


def sort_run(rundict_list, sort_key="nelect"):
    # Sort the valid vasprun according to their energy
    if sort_key is None:
        sort_key = "energy_per_fu"

    class Sort:
        def nelect(self, run):
            return run.parameters['incar'].get("nelect", 0)

        def energy_per_fu(self, run):
            return run.energy_per_fu

    print("sorting by {}".format(sort_key))
    sort = Sort()
    sorted_run_list = sorted(rundict_list,
                             key=getattr(sort, sort_key))
    return(sorted_run_list)


# TAG GENERATION FUNCTIONS
# ==========================================================================


# def generate_tags(rundict_list, force=False, minimal=False):
#     global PARAM

#     # TAG GENERATION =========================================================
#     # Computation of  various measures on the structure, called Tags
#     # these tags are printed to the log file, to the stdout and on the graph

#     # input : [{vasprun : v , 'folder' : f}]
#     # output : [{vasprun : v , 'folder' : f,
#     # nameTag : "" , energyTag : " ", structureTag:"" , baderTag : "" }]

#     # if the tags are already generated , skip this function (unless it's
#     # forced)
#     if PARAM['generated_tags'] and not force:
#         # print(rundict_list)
#         return(rundict_list)

#     vasp_run_poll = []
#     with Pool(processes=cpu_count()) as p:
#         #        zip(rundict_list, repeat(minimal))
#         vasp_run_poll = p.starmap(
#             get_tag_single_run, zip(
#                 rundict_list, repeat(minimal)))
#         p.close()
#         p.join()
#     # print(vasp_run_poll)

#     try:
#         # [runDict for runDict in vasp_run_poll
#         # if runDict['x_na'] <= 1 and runDict['x_na'] >= 0]
#         sorted_entries = sorted(vasp_run_poll,
#                                 key=lambda x: (x.x_na, x.energy_per_atom))
#     except Exception as ex:
#         logging.warning(
#             "failed to sort by (x_na , Etot), fall back to (x_na, folder)")
#         sorted_entries = sorted(
#             vasp_run_poll, key=lambda x: (x.x_na, x.job_folder))

#     # Create a New "log" file with incremented name
#     write_log_file(sorted_entries, PARAM['mainFolder'])

#     PARAM['generated_tags'] = True
#     # print(sorted_entries)

#     return(sorted_entries)


def write_log_file(rundict_list, log_folder):

    log_file_name = get_file_name(PARAM['mainFolder'], "log")
    os.chdir(log_folder)
    # log_file = open(logFileName,"a")
    if input("print runs name & energies ?[Y]") == "Y":
        print_stdout = True

    print("\nSaving logFile in  : {0} \n".format(log_file_name))
    with open(log_file_name, 'w') as log_file:
        for rundict in rundict_list:
            message = make_rundict_str(rundict)
            if print_stdout:
                print(message)
            log_file.write(message)


def make_rundict_str(rundict):
    message = "\n".join([
        "NAME",
        rundict.name_tag,
        rundict.folder,
        "ENERGY",
        rundict.energy_tag,
        "STRUCTURE",
        rundict.structure_tag,
    ])

    if getattr(rundict, 'bader_tag', None) is not None:
        message += "\n".join(["", "BADER", rundict.baderTag])
    return message

    # log_file.close()


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
    return(equivalent_sites)


# def get_structure_tag(structure, structureAnalysis=False):

#     analyzer = SpacegroupAnalyzer(structure, symprec=0.1)
#     spacegroup = analyzer.get_space_group_symbol()
#     structureTag = "{}".format(spacegroup)

#     if structureAnalysis:
#         # calculate the metal environment of each alkali and add it to the
#         # structureTag
#         alkaliEnvt = cluster.alkaliBondCount(structure)
#         # this function returns 2 arrays :one with the number of bonds
#         # and the other normalized by the number of Na (the one we plot)
#         coordString = "|"
#         for coordNumber in alkaliEnvt:
#             coordString += "{0:.0f} |".format(coordNumber * 100)
#         structureTag += "\nNa coord :\n" + coordString

#     # calculate the metal environment of each metal and add it to the
#     # structureTag
#         metalBond = cluster.metalBondCount(structure)
#         # this function returns 2 arrays : one with the number of bonds
#         # and the other normalized by the number of Na (the one we plot)
#         coordString = "|"
#         for coordNumber in metalBond[0]:
#             coordString += "{0:.0f} |".format(coordNumber * 100)
#         structureTag += "\nMe bonds [AA,BB,AB] :\n" + coordString
#     return(structureTag)


def get_file_name(folder, name, ext=""):
    # create a unique filename by adding
    # incremental suffix to the name in argument
    k = 0
    basis = os.path.join(folder, name) + "_"
    file_name = basis + str(k)
    while os.path.exists(file_name + ext):
        k += 1
        file_name = basis + str(k)
        print(file_name)
    return(file_name)


def restrict_run_list(all_runs_input):
    """
    defines a restriction on the runs to analyze and plot
    if the given list is not locked, performs filtering on the list
    """
    selection_loop = True

    while selection_loop:

        # STACKING FAMILY FILTER ================
        # ========================================
        stacking_list = list(set([e.stacking for e in all_runs_input]))
        family_array = stacking_list + ["finished"]
        family_choice = []
        print("which stacking types to plot ? ",
              "\n {}".format(family_array),
              "\nNo choice = No filtering")
        while "finished" not in family_choice:
            try:
                family_choice.append(
                    family_array[int(input(
                        "add stacking to current choice ({}) ? : ".format(
                            family_choice)))])
            except Exception as ex:
                print("family choice terminated {}".format(ex))
                family_choice.append('finished')
        if len(family_choice) > 1:
            restricted_stack_runs = [r for r in all_runs_input
                                     if r.stacking in family_choice]
        else:
            print("No stacking filter")
            restricted_stack_runs = all_runs_input

        # CONVEX HULL FILTER ==================
        # ======================================

        try:
            if len(set([d.x_na for d in restricted_stack_runs])) == 1 and input(
                    "convex hull on the doo? [Y]") in ["Y", "y"]:
                restricted_stack_runs = hull.generate_hull_entries(
                    restricted_stack_runs, remove_extremes=None, coord="dOO_min")
            elif input("convex hull on x_na ? [Y]") in ["Y", "y"]:
                restricted_stack_runs = hull.generate_hull_entries(
                    restricted_stack_runs, remove_extremes=None)

            # 3 converged plots, 4 : all minimas (even off_hull) 5 : On hull
            # minima
            list_choice = None
            choice_dict = {
                "[a]ll": "all runs that have at least a POSCAR ",
                "[c]onverged": " all converged vasprun sorted by x_na then energy",
                "[m]inima": " structures of lowest energy for each x_na ",
                "[h]ull": " structures on the convex hull"
            }
            while list_choice is None:
                print("Avaliable filtering : "+" // ".join(choice_dict.keys))
                list_choice = input("filter structure list ? :")
                choice_dict = {"a": 1, "c": 3, "m": 4, "h": 5}
                if list_choice in choice_dict.keys():
                    sieve_lvl = choice_dict[list_choice]
                else:
                    print("\n".join(["{} : {}".format(*kv)
                                     for kv in choice_dict.items()]))
                    list_choice = None
            restricted_hull_runs = [
                d for d in restricted_stack_runs if d.status >= sieve_lvl]
        except Exception as ex:
            print("{} : no hull filtering".format(ex))
            print(traceback.format_exc())
            restricted_hull_runs = restricted_stack_runs

        print("number of selected runs : {}".format(len(restricted_hull_runs)))

        restricted_range_runs = restricted_hull_runs
        try:
            [x_min, x_max] = [0, 1]
            if len(set([d.x_na for d in restricted_stack_runs])) > 1 and \
               input("Change default x_na range [{},{}] ? [Y]".
                     format(x_min, x_max))[0] in ["Y", "y"]:

                x_min = float(input("Xmin : "))
                x_max = float(input("Xmax : "))
                restricted_range_runs = [r for r in restricted_hull_runs
                                         if (r.x_na >= x_min and r.x_na <= x_max)]
        except Exception as ex:
            print("{} : no range filtering".format(ex))

        print(
            "number of selected runs : {}".format(
                len(restricted_range_runs)))

        # INDIVIDUAL SELECTION ============
        # ==================================
        restricted_idv_runs = []
        if input("individual run  selection ? [Y]/[n] ") == "Y":
            for run in restricted_range_runs:
                if input(
                        "keep : {} ? Y/y".format(run.name_tag)) in ["Y", "y"]:
                    restricted_idv_runs.append(run)
        else:
            restricted_idv_runs = restricted_range_runs

        print(
            "nb of structures selected : {0} ".format(
                len(restricted_idv_runs)))

        if len(restricted_idv_runs) > 0:
            selection_loop = False
        else:
            print("""
================================
bad selection, please do it again
================================
            """)

    restricted_idv_runs = sort_run(restricted_idv_runs, sort_key=None)

    return(restricted_idv_runs)


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

    log_file_name = get_file_name(PARAM['mainFolder'], "log")
    logging.basicConfig(filename=log_file_name, level=logging.DEBUG)
    logging.info('Started')

    parameters = PARAM
    return(parameters)


def get_nb_cell(structure):
    compo_dict = structure.composition.get_el_amt_dict()
    nb_cell = 1
    # sulfide or oxide ?
    for spec in ['S', 'O']:
        if compo_dict.get(spec, 0) > 0:
            nb_cell = compo_dict[spec] / 2
            break
    return(nb_cell)


def get_xna_and_formula(structure, nb_cell=1):
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
