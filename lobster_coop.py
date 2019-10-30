#!/usr/bin/env python3

# launchDisordered.py
# wait for the job in previousPath to pint a OK file before going
# Creates VASP imput files for vasp in the folder projectName/jobFformula


import importlib
import os
import re
import shutil
import subprocess
import time

# import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from pymatgen.core.structure import Structure
# from operator import itemgetter
from pymatgen.electronic_structure.cohp import CompleteCohp
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp import Vasprun
# from pymatgen.electronic_structure.core import Spin, OrbitalType
from pymatgen.io.vasp.inputs import Incar
from scipy.integrate import simps

import cluster as cluster
import DOS_plot as DOS
import generic_plot
import platform_id
import readO2 as O2
import readRun_entries as read

# import sys


"""
    a) Générer un ficher de point K dans lequel la multiplicité de chaque point est 1.

    cp INCAR OLD_INCAR
    cp KPOINTS OLD_KPTS

    a – modifier l’INCAR -

        Ajouter : LSORBIT = TRUE
        Modifier : ISMEAR = 0
        Enlever:  NPAR = 4  // ISPIN = 1
        NELM = 100
        NSW = 0

    b – Lancer un calcul non – colinéaire avec VASP
        /home/sol/Vasp/Vasp5/vasp.5.3.5-openmpi-nc/vasp
    b1 – kill the job

    c-  if IBZKPT generated :
        Fichier IBZKPT est écrit
        (SI le Fichier KPOINT n’a pas déja la forme voulue !!)
        Cp IBZPKT KPOINTS

    Le fichier IBZKPTS doit contenir les points K, et la dernière colonne,
    correspondant à la multiplicité des points K doit être égale a 1.
"""
settingDir = platform_id.setting_dir()


def generate_IBZKPT(folder):

    shutil.copy2('INCAR', 'OLD_INCAR')
    incar = Incar.from_file("INCAR")
    os.remove("INCAR")

    print("old incar \n", incar)
    incar["ISPIN"] = 1  # spin less calculation
    incar.pop("MAGMOM", None)  # we don't need that anymore !
    incar["ISMEAR"] = 0  # No smearing
    incar["LSORBIT"] = "TRUE"  # Non colinear calculation
    incar["NSW"] = 0  # single point
    incar["ISYM"] = 0  # no symmetry for kpts
    incar.pop("NPAR", None)  # no need for parallel computing

    print("new NC  incar \n", incar)
    incar.write_file("INCAR")

    jobFileName = os.path.join(settingDir,
                               "job_scripts",
                               "vasp_job_single_nc")
    name = "NC_{}".format(abs(hash(folder)))
    job_string = 'sbatch -J {} --workdir {}  {}'.format(
        name, folder, jobFileName)
    print(job_string)
    subprocess.call([job_string], shell=True)

    # waiting for the ibzkpt file to come out of the vasp run
    count = 0
    while count < 4 and os.path.exists("IBZKPT") is False:
        print("waiting for IBZKPT to be computed since {} secs".format(count * 5))
        time.sleep(5)
        count += 1

    # kill the job as soon as the IBZKPT is out
    job_cancel_string = "scancel --name {}".format(name)
    print(job_cancel_string)
    subprocess.call([job_cancel_string], shell=True)

    shutil.copy2('INCAR', 'NC_INCAR')
    shutil.copy2('OLD_INCAR', 'INCAR')
    os.remove('OLD_INCAR')

    if os.path.exists("IBZKPT"):
        print("IBZKPT sucessfully generated")
        return(True)
    else:
        print("NO IBZKPT in this folder !! ")
        return(False)


def launch_COOP(folder, rerun=False):
    print("yoyo")

    # get the INCAR and KPOINTS ready for COOP calculation and launch it

    shutil.copy2('INCAR', 'OLD_INCAR')
    incar = Incar.from_file("INCAR")
    os.remove("INCAR")

#   print(incar)
    incar["ISPIN"] = 2  # BACK IN DA SPIN CALCULATION
    incar["ISYM"] = 0  # no symmetry for kpts
    incar["ISMEAR"] = 0  # gaussian smearing of e-
    incar["SIGMA"] = 0.01  # small smearing
    incar["NSW"] = 0  # single point
    incar["LWAVE"] = "TRUE"  # Print the wavecar
    incar["PREC"] = "High"
    incar["NELM"] = 200  # very accurate electronic relaxation
    if rerun:
        # restart COOP calculation from previous results with constant basis
        # set
        incar["ISTART"] = 2

#    print("new COOP  incar \n" , incar)

    incar.write_file("INCAR")

    jobFileName = os.path.join(settingDir,
                               "job_scripts",
                               "vasp_job_single")
    name = "COOP_{}".format(abs(hash(folder)))
    job_string = 'sbatch -J {} --workdir {}  {}'.format(
        name, folder, jobFileName)
    print(job_string)
    subprocess.call([job_string], shell=True)

    return(True)


def make_COOP_dir(folder):

    COOP_dir = "{}/COOP".format(folder)
    COOP_dir = read.get_file_name(folder, COOP_dir)
    os.mkdir(COOP_dir)
    shutil.copy2('{}/INCAR'.format(folder), '{}/INCAR'.format(COOP_dir))
    shutil.copy2('{}/KPOINTS'.format(folder), '{}/KPOINTS'.format(COOP_dir))
    shutil.copy2('{}/CONTCAR'.format(folder), '{}/POSCAR'.format(COOP_dir))
    shutil.copy2('{}/POTCAR'.format(folder), '{}/POTCAR'.format(COOP_dir))
    return(COOP_dir)


def get_oo_pairs_by_dist(structure, bond_lenght_max=2.75):
    oo_pairs = []
    o_indices = structure.indices_from_symbol("O")
    nb_O = len(o_indices)
    dist_M = structure.distance_matrix
    for i in range(nb_O - 1):
        for j in range(i + 1, nb_O):
            ind_a = o_indices[i]
            ind_b = o_indices[j]
            if dist_M[ind_a, ind_b] < bond_lenght_max:
                oo_pairs.append([ind_a + 1, ind_b + 1])
    return (oo_pairs)


def make_lobsterin(folder, bond="O_O"):

    lobster_param = ["COHPstartEnergy", "COHPendEnergy",
                     "includeorbitals", "cohpbetween", "gaussianSmearingWidth"]

    s = {}
    comment = {}

    comment["COHPstartEnergy"] = "! First, enter the energetic window in eV \
(relative to the Fermi level) :"
    s["COHPstartEnergy"] = [-10]
    s["COHPendEnergy"] = [10]

    comment["includeorbitals"] = "\n! Then, specify the types of valence orbitals:"
    s["includeorbitals"] = ["s p d"]

    comment["basisfunctions"] = "\n! You can also specify the basis functions per element manually, e.g. :"
    s["basisfunctions"] = ["O 2s 2p"]
    # "Ir 7s 6d" ! Sr sv potential
    # "O 2s 2p"
    if "cohpbetween" in lobster_param:
        comment["cohpbetween"] = "\n! Now define the pairs for which COHP analysis etc. should be done.\n\
! The atoms are numbered as per order in the POSCAR file."
        if bond == "O_O":
            atom_pair_list = cluster.get_OO_pairs(
                Structure.from_file("{}/POSCAR".format(folder)))
        elif bond == "M_O":
            atom_pair_list = cluster.get_MO_pairs(
                Structure.from_file("{}/POSCAR".format(folder)),
                metal_str="Mn")
        print(atom_pair_list)
        s["cohpbetween"] = [
            "atom{} atom{}".format(
                atom_pair[0] +
                1,
                atom_pair[1] +
                1) for atom_pair in atom_pair_list]

    # If you are interested in single orbital COHPs, you can get all the pairs
    # like s-s, s-p x, ..., p z-p z. Uncomment this line to switch it on :
    # cohpbetween atom 1 atom 2 orbitalwise
    # TO DO : add orbital wise generation
    comment["cohpGenerator"] = "\n! If you want to generate the COHPPairs automatically, use this to include \n\
! all pairs in a given distance range (in Angstrom, not in atomic units) :"
    if "cohpGenerator" in lobster_param:
        print(comment["cohpGenerator"])
        min_dist = eval(input("min dist : "))
        max_dist = eval(input("max dist : "))

        s["cohpGenerator"] = ["from {} to {} ".format(min_dist, max_dist)]

    s["gaussianSmearingWidth"] = [0.01]

    print(s)

    with open("{}/lobsterin".format(folder), 'w') as lobsterin:
        for key in lobster_param:
            print(key)
            line = "{}\n".format(comment.pop(key, ""))
            lobsterin.write(line)
            for value in s[key]:
                line = "{}  {}\n".format(key, value)
                lobsterin.write(line)


def COOP_from_folder(mainDir, bond="O_O"):
    COOP_folders = []
    complete_coop_list = []
    subDirList, fileSystem = read.get_job_list(mainDir)
    for d in subDirList:
        complete_coop_list.append(get_COOP_from_folder(d, bond=bond))

    return(complete_coop_list)


# def COOP_from_runListDict(vasprunlistDict, bond="O_O"):
#     for v in vasprunlistDict:
#         # this entry may be None
#         v.complete_coop = get_COOP_from_folder(v.job_folder, bond=bond)

#     return(vasprunlistDict)


def check_valid_folder(mainDir, bond="O_O"):
    COOP_folders = []
    # print(mainDir)
    coop_folder_list = [
        os.path.join(
            mainDir,
            f) for f in os.listdir(mainDir) if f.startswith("COOP_")]

    if len(coop_folder_list) > 0:
        # print(coop_folder_list)
        for folder in coop_folder_list:
            print(folder)
            os.chdir(folder)
            # for bond in ["O_O","M_O"] :
            print("BOND DIR : {}".format(bond))
            bond_dir = os.path.join(folder, bond)
            if not os.path.exists(bond_dir):
                os.mkdir(bond_dir)
                shutil.copy2(os.path.join(folder, "POSCAR"),
                             os.path.join(bond_dir, "POSCAR"))
            if os.path.exists(os.path.join(bond_dir, "COOPCAR.lobster")):
                print("lobster done in {}".format(bond_dir))
                COOP_folders.append(bond_dir)
            elif os.path.exists(os.path.join(folder, "vasprun.xml")) \
                    and platform_id.running_on_cluster():

                try:
                    print("attempt to parse vasprun")
                    converged = Vasprun("vasprun.xml").converged
                except BaseException:
                    converged = False
                    print("corrupted vasprun")

                if not converged:
                    print("invalid vaspRun")
                    if input("istart=1 and relaunch COOP ? Y/ N  : ") == "Y":
                        launch_COOP(folder, rerun=True)
                else:
                    print("vasprun converged")
                    if not os.path.exists(os.path.join(bond_dir, "lobsterin")):
                        if input("make lobsterin ? Y/N  : ") == "Y":
                            importlib.reload(cluster)
                            if not os.path.exists(os.path.join(bond_dir, "POSCAR")):
                                shutil.copy2(os.path.join(folder, "POSCAR"),
                                             os.path.join(bond_dir, "POSCAR"))
                            make_lobsterin(bond_dir, bond=bond)

                    if os.path.exists(os.path.join(bond_dir, "lobsterin")):
                        print("On cluster and lobsterin found")
                        if input("launch lobster ? Y/N  : ") == "Y":
                            # move lobsterin from the bond_dir to the main folder ,
                            # do the lobster
                            # and move everything back to the bond_dir folder
                            shutil.copy2(os.path.join(bond_dir, "lobsterin"),
                                         os.path.join(folder, "lobsterin"))
                            subprocess.call(['lobster_coop'], shell=True)
                            for f in ["lobsterin", "lobsterout"] +\
                                    [f for f in os.listdir(folder) if f.endswith(".lobster")]:
                                shutil.move(os.path.join(folder, f),
                                            os.path.join(bond_dir, f))
                            if os.path.exists(
                                os.path.join(
                                    bond_dir,
                                    "COOPCAR.lobster")):
                                print("lobster done")
                                COOP_folders.append(bond_dir)

    else:
        print("No COOP folder found !")
        if platform_id.running_on_cluster():
            if input("launch NC run to prepare COOP ? Y/N  : ") == "Y":
                prepare_COOP(mainDir)

    return(COOP_folders)


def get_EY_below_efermi(E, Ycoop):
    E_below, Y_below = [], []
    for e, y in zip(E, Ycoop):
        if e < 0:
            E_below.append(e)
            Y_below.append(y)
    if len(E_below) % 2 == 1:
        E_below.pop(0)
        Y_below.pop(0)
    return(E_below, Y_below)


def bond_order(E, Ycoop):
    E_below, Y_below = get_EY_below_efermi(E, Ycoop)
    bond_below = simps(Y_below, E_below)
    elec_below = simps(np.absolute(Y_below), E_below)
    bond_total = simps(Ycoop, E)
    elec_total = simps(np.absolute(Ycoop), E)

    print(
        "e below Ef = {}, bond order = {} \n  elec total = {}  total bond order = {}".format(
            elec_below,
            bond_below,
            elec_total,
            bond_total))
    return(elec_below, bond_below, elec_total, bond_total)


def get_COOP_from_folder(folder, bond="O_O"):

    folders = check_valid_folder(folder, bond=bond)
    if len(folders) == 0:
        print("No COOP in this folder")
        return(None)
    os.chdir(folders[0])
    complete_coop = CompleteCohp.from_file("LOBSTER", are_coops=True)
    print("cohp sucessfully read from {}".format(folders[0]))
    print("    /\/\/\      ")
    for key, bond_dict in complete_coop.bonds.items():
        print(key, bond_dict)
    # for key, bond_dict in complete_coop.bonds.items():
    #     print("    /\/\/\      ")
    #     print(key, bond_dict)
    #     # \D matches non digit characters
    #     [ind_a, ind_b] = [(int(re.sub(r"\D", "", ind_str)) - 1)
    #                       for ind_str in key.split("-")]
    #     bond_dict["atom_pair"] = [ind_a, ind_b]

    #     # print(complete_coop.bonds)

    return(complete_coop)


def sort_key_by_dist(COOP, bond_length_max=2.4):
    sorted_coop_list_of_tuple = sorted(
        COOP.bonds.items(), key=lambda x: x[1]['length'])
    # print(sorted_coop_list_of_tuple)
    valid_key_list = []
    for bond in sorted_coop_list_of_tuple:
        if bond[1]['length'] <= bond_length_max:
            valid_key_list.append(bond[0])
            print(bond[0], bond[1]['length'])

    nb_coop = len(valid_key_list)
    print("{} single coop for doo < {}".format(nb_coop, bond_length_max))
    return(valid_key_list)


def sort_key_by_oxidation(runDict, COOP, nb_sites=1):
    labile_O_list = np.array(O2.labile_O_indices(runDict,
                                                 nb_sites=nb_sites))

    # sorted_coop_list_of_tuple = sorted(COOP.bonds.items(),
    # key= lambda x: (x[1]["atom_pair"][0], x[1]["atom_pair"][1] ) )

    # print(sorted_coop_list_of_tuple)
    valid_key_list = []
    for bond_key, bond_dict in COOP.bonds.items():
        for O_index in labile_O_list:
            if O_index in bond_dict["atom_pair"]:
                valid_key_list.append(bond_key)
                print(" {} : {} A".format(
                    bond_dict["atom_pair"], bond_dict["length"]))

    nb_coop = len(valid_key_list)
    print("{} single coop in oxidation order".format(nb_coop))

    sorted_keys = sorted(valid_key_list,
                         key=lambda x: COOP.bonds[x]['length'])
    print(sorted_keys)

    return(sorted_keys)


def get_OO_pair_keys(runDict, COOP):

    OO_pair_list = runDict["OO_pairs"]
    valid_key_list = []
    for bond_key, bond_dict in COOP.bonds.items():
        for i_O1, i_O2 in OO_pair_list:
            if i_O1 in bond_dict["atom_pair"] and \
               i_O2 in bond_dict["atom_pair"]:
                valid_key_list.append(bond_key)
                # print(" {} : {} A".format(
                #   bond_dict["atom_pair"], bond_dict["length"]))

                # OO_pair_list.remove([O1,O2])

    nb_coop = len(valid_key_list)
    print("{} single coop of OO pairs".format(nb_coop))

    print(valid_key_list)

    return(valid_key_list)


# ===================================================================
# Highest level functions

def prepare_COOP(folder):

    COOP_dir = make_COOP_dir(folder)

    os.chdir(COOP_dir)

    if generate_IBZKPT(COOP_dir):
        shutil.copy2("IBZKPT", "KPOINTS")

    launch_COOP(COOP_dir)


def plot_all_coop_for_all_folder(mainDir):

    coop_list = COOP_from_folder(mainDir)
    for COOP in coop_list:
        plot_all_coop_for_one_structure(COOP, nb_coop_to_display=4)

    input("[ENTER] to close all figures")
    plt.close("all")


def plot_COOP_OO(RunDict_list, sorting="oxidation"):
    override_sorting = True
    try:
        bond_choice = input("[M]-O // [O]-O // [A]ll: Bond choice ?")
        if bond_choice in ["O", "o"]:
            bonds = ["O_O"]  # , "O_O" "M_O"
        elif bond_choice in ["M", "m"]:
            bonds = ["M_O"]
        elif bond_choice in ["A", "a"]:
            bonds = ["M_O", "O_O"]
    except Exception:
        print('bond = M-O')
        bonds = ["M_O"]

# ========= getting the COOP data
    for v in RunDict_list:
        v.coop_dict = {}
        for bond in bonds:
            v.coop_dict[bond] = get_COOP_from_folder(v.job_folder, bond=bond)

    plot_dos = True
    dos_ax = 1 if plot_dos else 0

    Emin = -3
    Emax = 3


# ========= plotting all bond_COOP per run
    for v in RunDict_list:
        valid_coops = {bond: v.coop_dict[bond]
                       for bond in bonds
                       if v.coop_dict[bond] is not None}
        nb_axes = len(valid_coops)
        if nb_axes == 0:
            print("no valid COOP runs for {}".format(v.id))
            return([])

        nb_axes += dos_ax

        fig = plt.figure("COOP", figsize=(20, 10))
        if nb_axes == 1:
            axes = [fig.add_subplot(111)]
        else:
            axes = fig.subplots(nb_axes, 1, sharex="col")
            fig.subplots_adjust(hspace=0)

        bottom_COOP, top_COOP = [0, 0]
        i = dos_ax
        for (bond, coop) in valid_coops.items():
            plot_coop(axes[i], coop, Emin, Emax,
                      pmg_efermi=v.data['efermi'],
                      sorting="all", override_sorting=True, title=bond)
            bottom_COOP = min(axes[i].get_ylim()[0], bottom_COOP)
            top_COOP = max(axes[i].get_ylim()[1], top_COOP)
            i += 1

        # for ax in axes[dos_ax:]:
        #     ax.set_ylim([bottom_COOP, top_COOP])

        if plot_dos:
            ax_dos = axes[0]
            pre_defined_colors = {"O": "xkcd:light red",
                                  "Mn": "xkcd:bluish",
                                  "Mg": "#20c073",
                                  "Na": "xkcd:dandelion"}
            DOS.plot_DOS_on_axe(ax_dos, v, Emin, Emax,
                                DOS_choice="perElt", spin_choice=0,  # 0 = Up & dpwn
                                pre_defined_colors=pre_defined_colors)
            ax_dos.legend()
            # (bottom_DOS, top_DOS) = ax_dos.get_ylim()
            # ax_dos.set_ylim([top_DOS * (bottom_COOP / top_COOP), top_DOS])

            # axes[-1].text(-0.05 , 0.5 , "Na {} \nEf={:.4f}".format(v["x_na"],v["vaspRun"].efermi),
            #      weight='bold',size='large',
            #     horizontalalignment='right', verticalalignment='center',
            #     multialignment='right', transform=axes[-1].transAxes)

        plt.show(block=False)
        generic_plot.save_fig(fig, "COOP_{}"
                              .format(v.job_folder.split("/")[-2]))

# # =========================
#     # for bond in bonds:
#     #     # plot each bond_coop of all runs
#     for bond in bonds:
#         fig = plt.figure("COOP by bonds", figsize=(20, 10))
#         if len(complete_coop_list) == 1:
#             axes = [fig.add_subplot(111)]
#         else:
#             axes = fig.subplots(len(complete_coop_list), 1, sharex="col")
#             fig.subplots_adjust(hspace=0)

#         runList = COOP_from_runListDict(rundict_list, bond=bond)
#         valid_coop_runs = [v for v in runList if v.complete_coop is not None]

    # def plot_coop_list(complete_coop_list, axes=None):
    #     if len(complete_coop_list) == 0:
    #         print("no valid COOP runs to show")
    #         return([])

    #     if axes is not None:

    #         # for i, v in enumerate(valid_coop_runs):
    #         #     COOP = v.complete_coop
    #         #     if sorting == "all" or override_sorting:
    #         #         key_list = [k for k in COOP.bonds.keys()]
    #         #     elif sorting == "distance":
    #         #         key_list = sort_key_by_dist(COOP, bond_length_max=3)[0]
    #         #     elif sorting == "oxidation":
    #         #         key_list = sort_key_by_oxidation(v, COOP, nb_sites=1)[0]
    #         #     elif sorting == "OO_pairs":
    #         #         key_list = get_OO_pair_keys(v, COOP)

    return(RunDict_list)


def plot_all_coop_for_one_structure(COOP, nb_coop_to_display=4):
    folder = COOP.bonds["folder"]
    key_list = sort_key_by_dist(COOP, bond_length_max=3)

    if len(key_list) == 0:
        print("no coop in in {}".format(folder))
        return(True)

    key_list_short = [key_list[i]
                      for i in range(min(nb_coop_to_display, len(key_list)))]
    nb_coop = len(key_list_short)
    fig = plt.figure(folder, figsize=(20, 10))
    fig.suptitle(
        "COOP for O-O pairs in structure \n{}".format(COOP.bonds["folder"]))
    axes = []

    for i, key in enumerate(key_list_short):
        if i == 0:
            axes.append(fig.add_subplot(nb_coop, 1, nb_coop - i))
        else:
            axes.append(
                fig.add_subplot(
                    nb_coop,
                    1,
                    nb_coop - i,
                    sharex=axes[0]))
            plt.setp(axes[i].get_xticklabels(), visible=False)
        plot_coop(axes[i], COOP, key)

    plt.show(block=False)
    read.save_fig(fig, "COOP_{}_{}"
                  .format(folder.split("/")[-3],
                          folder.split("/")[-2].split("-")[0]))


def plot_coop(ax, COOP, Emin, Emax,
              pmg_efermi=None,  title="",
              sorting="all", override_sorting=True):
    # although lobster shifts energies to align 0 and Efermi
    # the pymatgen import method re_shift it back to the vasp like format
    if sorting == "all" or override_sorting:
        key_list = [k for k in COOP.bonds.keys()]

        if len(key_list) == 0:
            print("no coop to plot")
            return(ax)
    pmg_fermi_first = True
    alternative_fermi_lvl = 0

    if pmg_fermi_first and pmg_efermi is not None:
        efermi = pmg_efermi
        alternative_fermi_lvl = COOP.efermi
    elif COOP.efermi != 0:
        efermi = COOP.efermi
        try:
            alternative_fermi_lvl = pmg_efermi
        except Exception as ex:
            alternative_fermi_lvl = 0
    elif pmg_efermi is not None:
        efermi = pmg_efermi
        alternative_fermi_lvl = COOP.efermi
    else:
        print("no valid efermi found !!")

        # print("recieved list : {}".format( key_list ) )
    E = COOP.all_cohps[key_list[0]].energies - efermi

    Ycoop = np.ufunc.reduce(
        np.add, [
            COOP.all_cohps[key].get_cohp(
                spin="up", integrated=False)[Spin.up] +
            COOP.all_cohps[key].get_cohp(
                spin="down", integrated=False)[Spin.down]
            for key in key_list])

    ax.plot(E, Ycoop, color='k', lw=2)

    ax.fill_between(
        E,
        Ycoop,
        0,
        where=Ycoop >= 0,
        facecolor='blue',
        interpolate=True,
        alpha=0.3)
    ax.fill_between(
        E,
        Ycoop,
        0,
        where=Ycoop <= 0,
        facecolor='red',
        interpolate=True,
        alpha=0.3)

    # ax.fill(E, Ycoop_bond, "b", E, Ycoop_antibond,  "r")
    ax.axhline(y=0, color='k', lw=0.5)

    if pmg_efermi is not None:
        pmg_fermi_lvl = pmg_efermi - efermi
        ax.axvline(x=pmg_fermi_lvl, color='k')

    indMin = 0
    indMax = 0
    for q, e in enumerate(E):
        if e < Emin:
            indMin = q
        if e < Emax:
            indMax = q

    ymax = max(Ycoop[indMin:indMax])
    ymin = min(Ycoop[indMin:indMax])

    ax.set_ylim([ymin * 1.1, ymax * 1.1])
    ax.set_xlim([Emin, Emax])
    ax.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        right=False,      # ticks along the bottom edge are off
        left=True,         # ticks along the top edge are off
        labelleft=True)  # labels along the bottom edge are off

    # avg_length = np.mean([COOP.bonds[key]['length'] for key in key_list])

    # ax.text(0.975 , 0.025 ,
    #         "d={:.4}\nEf={:.4}".format(avg_length, COOP.efermi),
    #         weight='bold',size='large',
    #         horizontalalignment='right', verticalalignment='bottom',
    #         multialignment='right', transform=ax.transAxes)

    elec_below, bond_below, elec_total, bond_total = bond_order(E, Ycoop)
    bond_order_txt = " {} \n pop : {:.3} / {:.3} = {:.3} \n\
    bond :\n {:.3} => {:.3} \n (tot bond {:.3})".format(
        title,
        elec_below, elec_total, elec_below / elec_total,
        bond_below, bond_below * 4 / elec_total, bond_total)

    ax.text(1.005, 0.5, bond_order_txt,
            size='small',
            horizontalalignment='left', verticalalignment='center',
            multialignment='left', transform=ax.transAxes)

    return([ymin, ymax])


if __name__ == '__main__':

    folder = os.getcwd()

    platform_id.first_check_cluster()

    print("CWD = {}".format(folder))
    if input("[R]ead or [W]rite COOP computation ?") == "R":
        plot_COOP_OO(folder)
    elif input("[W]rite COOP ?") == "W":
        prepare_COOP(folder)
