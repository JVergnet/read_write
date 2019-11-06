#!/usr/bin/env python3
# launchDisordered.py
# wait for the job in previousPath to pint a OK file before going
# Creates VASP imput files for vasp in the folder projectName/jobFformula


# import sys
import json
import os
# import subprocess
# import shutil
import time
from itertools import chain
from multiprocessing import Pool, cpu_count

import numpy as np
from pymatgen import Structure
# from pymatgen.io.vasp.inputs import Incar
# from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.sets import \
    MPRelaxSet  # , MVLRelax52Set  # , MPStaticSet

import createStructureList as create_list
import readRun_entries as read


class Job(MPRelaxSet):
    @property
    def structure(self):
        return self._structure

    @structure.setter
    def structure(self, value):
        self._structure = value

    def __init__(self, structure, entry_id,
                 user_param=None, user_incar=None,
                 job_folder="", **kwargs):
        # super().__init__(structure)
        self.explicit_jobpath = False
        self.entry_id = entry_id
        self.job_folder = job_folder
        # print(structure)

        self.user_kpoint = {'reciprocal_density': 100}
        self.user_param = user_param
        self.user_incar = user_incar if user_incar is not None \
            else default_incar()

        self.kwargs = kwargs
        for (k, v) in kwargs.items():
            setattr(self, k, v)

        # print(self.__dict__)
        self.standardize = False
        self.structure = structure

        self.nb_cell = 1
        self.x_na = 0
        self.formula = ""
        self.oldFolder = None
        self.job_name = self.set_job_name()

    # @structure.deleter
    # def structure(self):
    #     del self.structure

    def get_nametag(self):
        " name_tag : get name of the current structure in a string nameTag"
        self.nb_cell = read.get_nb_cell(self.structure)

        # Normalizing composition to get Nax My O2
        self.x_na, self.formula = read.get_xna_and_formula(
            self.structure, self.nb_cell)

    def set_job_name(self):
        self.get_nametag()
        self.job_name = self.formula + "__" + self.entry_id
        self.job_name = self.job_name.replace(' ', '-').replace('.', '')
        return self.job_name

    def set_job_folder(self, parent_path, explicit_jobpath=None):
        if explicit_jobpath is not None:
            self.explicit_jobpath = explicit_jobpath
        self.set_job_name()
        self.job_folder = parent_path if self.explicit_jobpath \
            else os.path.join(parent_path, self.job_name)
        return self.job_folder

    def write_data_input(self, parent_path=None, explicit_jobpath=None):
        if parent_path is not None:
            self.set_job_folder(parent_path, explicit_jobpath)
        print(self.kwargs)
        # self.kwargs['user_kpoints_settings'] = self.user_kpoint
        self.user_incar['SYSTEM'] = self.job_name
        # print(self.user_incar)

        if self.user_incar.get('MAGMOM', None) is not None:
            self.structure.add_spin_by_site(self.user_incar['MAGMOM'])
            self.user_incar.pop('MAGMOM')
        print(self.user_incar)
        super().__init__(self.structure,
                         user_kpoints_settings=self.user_kpoint,
                         force_gamma=True,
                         user_incar_settings=self.user_incar,
                         potcar_functional='PBE_54',
                         )
        # print("POTCAR FUNCTIONNAL :", self.potcar_functional)
        self.write_input(self.job_folder)

        if len(self.user_param) > 0:
            # we also dump the data of the run into the job_folder
            print(json.dumps(self.user_param))
            with open('{}/param.json'.format(self.job_folder), 'w') as outfile:
                json.dump(self.user_param, outfile)

    def copy(self):
        return(Job(self.structure.copy(), self.entry_id,
                   user_param=dict(self.user_param),
                   user_incar=dict(self.user_incar),
                   job_folder=self.job_folder))

    def copy_w_new_struct(self, structure, new_id="", id_mode="add"):
        id_str = ""
        if id_mode in ["add", "a", "A"]:
            id_str += self.entry_id
        # elif id_mode in ["replace", "r", "R"] and new_id != "":
        id_str += new_id
        return(Job(structure.copy(), id_str,
                   user_param=dict(self.user_param),
                   user_incar=dict(self.user_incar),
                   job_folder=self.job_folder))

    @classmethod
    def from_rundict(cls, rundict, new_folder=None):
        job = Job(rundict.structure, rundict.str_id,
                  user_param=dict(rundict.parameters.get('custom', {})),
                  user_incar=dict(rundict.parameters['incar']),
                  job_folder=rundict.job_folder)
        if new_folder is not None:
            job.oldFolder = job.job_folder
            job.set_job_folder(new_folder)
        return(job)


def default_incar():

    incar_default = {
        'SYSTEM': 'honeycomb',

        'ENCUT': 600, 'PREC': 'Normal',

        'EDIFF_PER_ATOM': 1e-03, 'EDIFFG': -0.01,

        'NELMDL': -12, 'NELM': 60, 'NELMIN': 4,

        # fichier DOS (impair pour fermi lvl)
        'NEDOS': 1001,
        # -5 pour les isolants (et simga=0)
        'ISMEAR': 0,
        # pametre de smearing( pour converger l'entropie)
        'SIGMA': 0.2,
        'IBRION': 2,

        'LCHARG': '.TRUE.', 'LAECHG': '.TRUE.',  # edition fichier charge
        'LASPH': '.TRUE.',                      # non isotropic d orbitals
        'ISIF': 3,                             # = 0 : juste les ions , =3 : cellule aussi
        'NSW': 200,                             # Nombre pas maximum relax ionique
        'ISYM': 0,                             # = 1 symetrie est gardee
        'NPAR': 4,                               # Parallelisation des bandes
        'LMAXMIX': 4,
        'IVDW': 12
    }
    return(incar_default)


def get_pristine_list(structure_files, cif_folder, is_poscar=False):
    os.chdir(cif_folder)

    # LIST OF PRISTINE_STRUC (fct of CIF_NAME list)
    # ===================================
    pristine_job_list = []
    for file_name in structure_files:
        if is_poscar:
            name = file_name.replace('_POSCAR', '')
        else:
            name = file_name.split(".")[-1]
        print("Structure raw file : ", file_name)
        pristine_struct = Structure.from_file(file_name)
        pristine_struct.make_supercell([1, 1, 3], to_unit_cell=True)
        # pristineStruct.replace_species({"Na": "Li"})
        pristine_job_list.append(Job(pristine_struct,
                                     name))
        print("OK\n")

    return(pristine_job_list)


def get_structure_list_wrapper(pristine_structure_list,
                               launch_choice):

    return(get_complete_structure_list(pristine_structure_list,
                                       launch_choice=launch_choice))


def get_complete_structure_list(
        pristine_job_list,
        launch_choice=None,
        parameter=["Na+"]):

    # LIST OF {PMG.STRUCTURES/ID} (fct of PRISTINE_STRUCT list)
    # ==================================
    # data structure : list of dict { 'structure' : pmg.struct , "id" : string
    # }
    choice_dict = {
        "s": "[s]imple",
        "d": "[d]esodiation_tree",
        "o": "[o]rder metallic",
        "l": "[l]ayer rotation",
        "r": "[r]eplacement",
        "b": "[b]ailar & trigonal distortions",
        "t": "bulk s[t]rain",
        "n": "[n]updown scan"
    }
    if launch_choice is None:
        print("Available functions to modifiy initial structures")
        print("\n".join(
            ["{}:{}".format(k, v)for (k, v) in choice_dict.items()]
        ))
        launch_choice = input("\nLaunch type (by letter) ?  ")

    job_list = []
    print("chosen structure treatment  : {}\n".format(
        choice_dict[launch_choice[0]]))

    for pristine_job in pristine_job_list:
        # pristine_job = read.get_nametag(pristine_job)
        if launch_choice[0] == "s":
            struct_list_tmp = [pristine_job]

        elif launch_choice[0] == "d":
            # Leaving  Na (in proportion) from starting structure (x=1)
            # (computing NaxMO2)
            struct_list_tmp = create_list.desodiation_tree(
                pristine_job,
                xmin=0,
                xmax=0.9,
                nb_compo=10,
                nb_struct_per_x=1,
                alkali="Na")

        elif launch_choice[0] == "l":
            struct_list_tmp = create_list.rotation_in_supercell(
                pristine_job,
                supercell_size=[1, 1, 1],
                angle_max=180, angle_step=15)

        elif launch_choice[0] == "o":
            struct_list_tmp = create_list.disorder_in_supercell(
                pristine_job,
                # {"Li+":0.67 ,"Ti2+" : 0.33},
                # substSpecies=None,
                supercell_size=[1, 1, 1],  # [2,2,2],
                number_of_struct=4)

        elif launch_choice[0] == "r":
            struct_list_tmp = create_list.replace(pristine_job,
                                                  # ["Li+" , "Ca+"]
                                                  substSpecies=parameter,
                                                  supercell_size=[1, 1, 1])
        elif launch_choice[0] == "b":
            # bailar = { "disto_min" : 0.07 , "disto_max" : 0.22, "nb_steps" : 6  }
            # trigo  = { "disto_min" : 0.13, "disto_max" : 0.27, "nb_steps" : 6 }
            bailar = {"disto_min": 0.3, "disto_max": 0.6, "nb_steps": 8}
            trigo = {"disto_min": 0.3, "disto_max": 0.6, "nb_steps": 8}
            struct_list_tmp = create_list.OO_distortion(
                pristine_job,
                bailar['disto_min'],
                bailar['disto_max'],
                bailar['nb_steps'],
                trigo['disto_min'],
                trigo['disto_max'],
                trigo['nb_steps'])
        elif launch_choice[0] == "t":
            struct_list_tmp = create_list.bulk_strain(pristine_job,
                                                      0.04,
                                                      4)
        elif launch_choice[0] == "n":
            struct_list_tmp = create_list.scan_nupdown(pristine_job, 8)

        # if launch_choice[0] != "s":
        #     for struct in struct_list_tmp:
        #         struct.entry_id = pristine_job.entry_id + struct.entry_id
        job_list += struct_list_tmp

    if launch_choice[0] not in ["b", "s", "t"]:
        job_list = create_list.remove_doubles(job_list)

    return(job_list)


def adjust_incar(job_list, incar_setting=None,
                 VdW=None, perturb=None):
    incar_out = incar_setting if incar_setting is not None else {}
    if VdW is None:
        VdW = True if input(
            "\n Van Der Waals ? Y / N   :   ") == "Y" else False

    if VdW:
        incar_out['IVDW'] = 12

    if perturb is None:
        try:
            perturb = float(
                input('Perturb the initial position of atoms ? in Angstrom '))
        except Exception as ex:
            print("No perturbation")
            perturb = 0

    for job in job_list:
        job.structure.perturb(perturb)
        job.user_incar.update(incar_out)

    return(True)


def generate_job_folders(final_job_list,
                         parent_folder, project_name,
                         selective_dynamic=None):

    # LIST OF JOB FOLDER (fct of VASP_INPUT_SET list)
    # ====================================================
    # Search and create a new directory for the project

    project_dir = read.get_file_name(parent_folder, project_name, ext="")

    print(project_dir)
    os.mkdir(project_dir)

    folder_list = []
    if selective_dynamic not in [None, "False"]:
        make_selective_dynamic(final_job_list, selective_dynamic)

        # if fukui is not None :
        #     inc = vasp_input.incar
        #     inc["NELECT"] =  vasp_input.nelect + fukui
        #     pos.write_file(folder_name+"/INCAR")

    for j, job in enumerate(final_job_list):
        job.write_data_input(project_dir)

        print("set written in {}".format(job.job_folder))

        folder_list += []

    return(folder_list)


def make_selective_dynamic(final_job_list, selective_dynamic):
    print("all sites are allowed to move")
    for job in final_job_list:
        print(job.job_name)
        mysd = np.ones([job.structure.num_sites, 3], bool)
        if selective_dynamic == "all_sites":
            allowed_sites = [n for n in range(job.structure.num_sites)]
        else:
            allowed_sites = []
            if selective_dynamic != "no_sites":
                for spec in selective_dynamic:
                    allowed_sites += job.structure.indices_from_symbol(
                        spec)
        indices_to_block = []

        while True:
            try:
                j = int(input("block site by index?"))
                indices_to_block.append(j)
                print("blocked indices : {}".format(indices_to_block))
            except Exception as ex:
                break

        for k in indices_to_block:
            mysd[k][1] = False
            mysd[k][0] = False
            mysd[k][2] = False
        # print(mysd)
        # pos = vasp_input.poscar
        # pos.selective_dynamics = mysd
        # pos.write_file(folder_name + "/POSCAR")
        job.structure.add_site_property("selective_dynamics", mysd)


def main():

    print("ALL RIGHT ! READY TO ROCK ! \n \n" + time.asctime() + "\n  \n")

    # incar_default = default_incar()

    # cif_list=["Na2MnO3_mp-769949_symmetrized.cif"]
    project_name = "new_project"
    parent_dir = "/home/jvergnet/frofro/Cu_phase/CuX1_2/"

    cif_folder = "/home/jvergnet/frofro/Cu_phase/Ni_Mn_Ti01__Cu01/first_scan_Na0/migrated/4_0"

    # cif_list_O3 = ["O3_POSCAR",
    #                "P3_POSCAR"]

    # cif_list_O2 = [  # "P2_EEL_POSCAR",
    #     # "P2_EF_POSCAR",
    #     "P2_EEL_POSCAR"
    #     # "P2P_EE_POSCAR",
    #     # "P2P_EF_POSCAR"
    #     # "MdB_L_POSCAR"
    # ]
    # # cif_list_tot = cif_list_O2 + cif_list_O3
    cwd = os.getcwd()
    if input("Work in CWD ? Y/N") == "Y":
        cif_folder = cwd
    else:
        cif_folder = str(input("Enter path to folder : \n"))

    print(" working in {}".format(cif_folder))

    file_list = [f for f in os.listdir(cif_folder) if (
        f.endswith('POSCAR') or f.endswith('cif'))]

    if len(file_list) > 0 and \
       input("found {} \n Work with these files ? Y/N".format(
           file_list)) in ["Y", "y"]:
        pristine_job_list = get_pristine_list(
            file_list, cif_folder, is_poscar=True)

    elif input("Browse existing runs in subfolder ? Y/N") == "Y":
        rundicts = read.collect_valid_runs(
            cif_folder, vaspRun_parsing_lvl=0)
        pristine_job_list = [Job(r.structure, r.str_id,
                                 user_incar=dict(r.parameters["incar"]))
                             for r in rundicts]
    else:
        pristine_job_list = []
    if len(pristine_job_list) == 0:
        print("You have not selected a file to work on !")
        return("Not OK")

        # ATTENTION !!!!
        # for entry in pristine_job_list :
        #     entry["structure"].remove_oxidation_states()
        #     entry["structure"].remove_species([Element("Na")])

    # for r in pristine_job_list:
    #     r['structure'].remove_site_property('selective_dynamics')

    job_list = get_complete_structure_list(
        pristine_job_list)

    adjust_incar(job_list, incar_setting={},
                 VdW=None, perturb=None)

    if input("create run folders in CWD ? Y/N") == "Y":
        parent_dir = cwd

    generate_job_folders(job_list,
                         parent_dir, project_name,
                         selective_dynamic=None)

    return("OK")


def substitution_desodiation_P2_P3():
    # Set the parameters ofthe run
    # parent_dir = "/home/jvergnet/frofro/honeycomb/"
    project_name = "substitued"
    # oxyde = True

    cif_folder = "/home/jvergnet/frofro/honeycomb/ref_POSCAR_MgMn"
    cif_list_O3 = ["O3_POSCAR",
                   "P3_POSCAR"]

    cif_list_O2 = [  # "P2_EEL_POSCAR",
        # "P2_EF_POSCAR",
        # "O2_L_POSCAR"
        # "P2P_EE_POSCAR",
        # "P2P_EF_POSCAR",
        "OP4_POSCAR"
    ]
    # cif_list_tot = cif_list_O2 + cif_list_O3

    # name_list = ["a", "a", "c"]

    cluster = True

    incar_default = default_incar()

    if cluster:
        parent_folder = "/home/jvergnet/frofro/honeycomb/"
    else:
        parent_folder = "/mnt/fbf83589-c1e0-460b-86cc-ec3cc4a64545/backup_cluster/backup_fro/honeycomb/"

    project_name = "Hull_all"

    # get the O2 /P2 / O3 with Mg
    Mg_pristine_stackings = get_pristine_list(
        cif_list_O2, cif_folder, is_poscar=True)
    # [ {'structure' : O2 , "id" : "O2 } ]

    projectDir = parent_folder + project_name
    i = 0
    while (os.path.exists(projectDir)):
        i += 1
        projectDir = parent_folder + project_name + str(i)

    os.mkdir(projectDir)
    os.chdir(projectDir)

    for specie in ["Zn"]:

        os.mkdir(specie)
        os.chdir(specie)
        specie_dir = projectDir + "/" + specie

        # get the pristine stackings for each specie
        # specie_pristine_stackings = get_complete_structure_list(Mg_pristine_stackings ,
        #                                                          launch_choice='r',
        #                                                         parameter  = [specie] )
        # specie_pristine_stackings = Mg_pristine_stackings

        specie_pristine_stackings = []
        for stacking in Mg_pristine_stackings:
            entry = dict(stacking)
            entry["structure"].replace_species(
                {Element("Mg"): Element(specie)})
            # entry["structure"].replace_species({ Element("Mn") : Element("Co")})
            specie_pristine_stackings.append(entry)

        for specie_stacking in specie_pristine_stackings:
            print(specie_stacking['id'])

            # desodiate each stacking and write create folders for each state
            final_job_list = get_complete_structure_list(
                [specie_stacking], launch_choice='d')

            # final_job_list = get_inputSet_list(desodiated_specie_stacking,
            #                                    incar_default,
            #                                    VdW=True, scan=False)
            adjust_incar(final_job_list, incar_setting={},
                         VdW=None, perturb=None)
            generate_job_folders(final_job_list, specie_dir,
                                 specie_stacking['id'])

        os.chdir(projectDir)

    return("OK")


def multi_distortion_desodiation_scan():
    cluster = True

    incar_static = default_incar()
    incar_static.update({
        'ENCUT': 500,
        'PREC': 'Normal',
        'EDIFF': 1e-04,
        'NELMDL': -15,
        'NELM': 400,
        # fichier DOS (impair pour fermi lvl)
        'NEDOS': 1001,
        # -5 pour les isolants (et simga=0)
        'ISMEAR': 0,
        # pametre de smearing( pour converger l'entropie)
        'SIGMA': 0.05,
        'IBRION': 2,
        'NSW': 300,                             # Nombre pas maximum relax ionique
        # on ne relaxe que la position des sodium, pas de la cellule
        'ISIF': 0,
        'LCHARG': 'FALSE', 'LAECHG': 'FALSE',  # edition fichier charge
        'LASPH': '.TRUE.',                      # non isotropic d orbitals
        'ISYM': 0,                             # = 1 symetrie est gardee
    })

    if cluster:
        parent_folder = "/home/jvergnet/frofro/honeycomb/"
    else:
        parent_folder = "/mnt/fbf83589-c1e0-460b-86cc-ec3cc4a64545/backup_cluster/backup_fro/"

    projectName = "prismatic_zoom"
    cif_folder = "/home/jvergnet/frofro/honeycomb/ref_POSCAR_MgMn/"
    cif_list = ["SYM_P2EEL_POSCAR"]
    # get the O2 /P2 / O3 with Mg
    pristine_stackings = get_pristine_list(
        cif_list, cif_folder, is_poscar=True)
    # [ {'structure' : O2 , "id" : "O2 } ]
    print(pristine_stackings)
    projectDir = read.get_file_name(parent_folder, projectName, ext="")

    os.mkdir(projectDir)
    os.chdir(projectDir)

    for specie_stacking in pristine_stackings:
        # stacking_dir = projectDir+"/"+specie_stacking['id']
        # print(stacking_dir)
        # os.mkdir(stacking_dir)
        # os.chdir(stacking_dir)
        # all_structures =  []

        ordered_stackings = [{'structure': specie_stacking['structure'],
                              'id': specie_stacking['id']}]
        # get_complete_structure_list([specie_stacking], launch_choice='s')
        for d in ordered_stackings:
            d['structure'].remove_species(["Na"])
        print([s["structure"] for s in ordered_stackings])
        with Pool(processes=cpu_count()) as p:
            list_of_struct_list = p.starmap(
                get_structure_list_wrapper, [
                    ([d], 'b') for d in ordered_stackings])
            p.close()
            p.join()

        final_job_list = list(chain(*list_of_struct_list))
        # for desodiated_structure in desodiated_stackings :
        #     print(desodiated_structure["id"])
        #     # desodiate each stacking and write create folders for each state
        #     all_structures += get_complete_structure_list([desodiated_structure],launch_choice='b')
        #     # for prismatic_disto_struct in prismatic_distorted_structures :
        #     #     print(prismatic_disto_struct["id"])
        #     #     all_structures += get_complete_structure_list([prismatic_disto_struct],launch_choice='t')

        # , name_format="sulfide")

        adjust_incar(final_job_list, incar_setting={},
                     VdW=None, perturb=None)

        generate_job_folders(final_job_list,
                             projectDir, specie_stacking['id'],
                             selective_dynamic=True)

    os.chdir(projectDir)

    return("OK")


if __name__ == '__main__':
    # substitution_desodiation_P2_P3 ()
    main()
    # multi_distortion_desodiation_scan()
