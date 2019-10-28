# readO2.py

# read and write O2 release enthalpy related quantities


import os
import subprocess
from multiprocessing import Pool, cpu_count
from operator import itemgetter

from pymatgen.io.vasp.outputs import Oszicar

import generic_plot as generic_plot
# import lobster_coop as lob
import launchDisordered as launch
# import read_hull as hull
import platform_id
# import matplotlib.pyplot as plt
import readRun_entries as read


def O2_computation(dosList, bader_done=False):
    vasp_run_poll = []
    with Pool(processes=cpu_count()) as p:
        vasp_run_poll = p.map(get_O2_tag_single, dosList)
        p.close()
        p.join()

    if platform_id.running_on_cluster() and bader_done and input(
            "generate O deficient structure ? Y / N  : ") == "Y":
        list_to_remove_O2 = [
            d for d in vasp_run_poll if d.get(
                "O2_release_enthalpy",
                None) is None]
        remove_O2(list_to_remove_O2, nb_sites=1, redundancy=False)

    valid_runs = [
        d for d in vasp_run_poll if d.get(
            "O2_release_enthalpy",
            None) is not None]
    print("{} valid runs".format(len(valid_runs)))
    if input("plot O2 release  evolution ? Y / N  : ") == "Y":
        generic_plot.plot_structure_value_evolution(
            valid_runs, prop_list=['O2_release_enthalpy'])

    return(vasp_run_poll)


def get_O2_tag_single(runDict):
    if runDict.get("O2_release_enthalpy", None) is None:
        O2_enthalpy = get_O2_release_enthalpy(runDict)
        if O2_enthalpy is not None:
            print("\n O2 release enthalpy = {}\n".format(O2_enthalpy))
            runDict['O2_release_enthalpy'] = O2_enthalpy

    return(runDict)


def remove_O2(vaspRundictList, nb_sites=1, redundancy=None):

    # for each structure given in argument,
    # remove the most oxidized oxygen
    # and create a VASP input folder in the job folder
    # option for redundancy to generate a larger supercell

    if redundancy is None:
        if input("redundancy ? Y/N")[0] == "Y":
            redundancy = True
        else:
            redundancy = False

    folder_list = []

    for run in vaspRundictList:
        structure_list = []
        parentDir = run['folder']
        projectName = "O_deficient"
        print("removing O from {}".format(parentDir))
        # print(run["equivSiteList"])

        for O_index in labile_O_indices(run, nb_sites):
            print("removed O index : {}".format(O_index))
            O_deficient_structure = run['vaspRun'].final_structure.copy()
            O_deficient_structure.remove_sites([O_index])
            structure_list.append({'structure': O_deficient_structure,
                                   "id": "rem_O{}".format(O_index)})
            if redundancy:
                big_struct = run['vaspRun'].final_structure.copy()
                big_struct.make_supercell([1, 2, 1])
                structure_list.append({'structure': big_struct,
                                       "id": "rem_O{}_big".format(O_index)})

        # print("O deficient structure list\n {}".format(structure_list))

        incar_default = launch.default_incar()

        vasp_input_set_list = launch.get_inputSet_list(structure_list,
                                                       incar_default,
                                                       VdW=False, scan_U=False,
                                                       oxyde=False)

        sub_folder_list = launch.generate_job_folders(vasp_input_set_list,
                                                      parentDir,
                                                      projectName)
        folder_list += [os.path.join(run['folder'], sub_folder)
                        for sub_folder in sub_folder_list]

    for folder in folder_list:
        settingDir = platform_id.setting_dir()
        jobFileName = os.path.join(
            settingDir, "job_scripts", "vasp_job_double")
        name = "degaz_{}".format(abs(hash(folder)))
        # workingDir = os.path.join(parentDir,folder)
        job_string = 'sbatch -J {} --workdir {}  {}'.format(
            name, folder, jobFileName)
        if input(
                "LAUNCH the following [Y/n] :\n{}\n".format(job_string)) == "Y":
            subprocess.call([job_string], shell=True)

    return(folder_list)


def get_O2_release_enthalpy(runDict):

    run_list = []

    # check the existence of the O2 deficient vasp result in the run folder
    # print(runDict)
    o_deficient_list = [
        os.path.join(
            runDict['folder'],
            f) for f in os.listdir(
            runDict['folder']) if f.startswith("O_deficient")]

    if len(o_deficient_list) == 0:
        print("no O_deficient folder")
        return(None)

    print(o_deficient_list)
    for project_folder in o_deficient_list:
        for job_folder in [
            os.path.join(
                project_folder,
                f) for f in os.listdir(project_folder)]:
            # print(job_folder)
            run = read.collect_single_folder(job_folder)
            if run.status >= 3:
                run_list.append(run)

    if len(run_list) == 0:
        print("no O_deficient converged run")
        return(None)

    #print("valid run lst {}".format(run_list))

    # get the run with the lowest energy
    run_list = sorted(run_list, key=lambda x: x['vaspRun'].final_energy)
    o_def_vasprun = run_list[0]["vaspRun"]
    # ============= IN CASE OF VDW CORRECTION =========================
    # TODO:
    E_no_vdw_O_def = Oszicar(
        run_list[0]['folder'] + "/OSZICAR").electronic_steps[-1][-1]["E"]
    # =================================================================
    E_plus_vdw_O_def = o_def_vasprun.final_energy

    # print("Energies 0 deficient : base {:.4f} // VDW {:.4f} (corr : {:.6f})".format(
    # E_no_vdw_O_def, E_plus_vdw_O_def, E_plus_vdw_O_def -  E_no_vdw_O_def  )
    # )

    #nbO_final = len(o_def_vasprun.final_structure.indices_from_symbol("O"))

    E_plus_vdw_normal = runDict['etot']
    E_no_vdw_normal = Oszicar(
        runDict['folder'] + "/OSZICAR").electronic_steps[-1][-1]["E"]

    # print("Energies : base {} // VDW {} (corr : {})".format(
    # E_no_vdw_normal, E_plus_vdw_normal, E_plus_vdw_normal -  E_no_vdw_normal
    # )  )

    #nbO_init = len(runDict['structure'].indices_from_symbol("O"))

    # H = E(normal) - E(O_deficient) - 1/2 E(O2)

    # runDict['ediff'] = runDict['etot'] / runDict["nb_cell"]
    # eV  -9.8897568 (no Wdv energy from dft O2 ) + 1.36 (ceder correction)
    EO2 = -8.52
    print("{}\n".format(runDict["nameTag"]))
    for VDW, E_normal, E_O_def in \
        [("with VdW correction", E_plus_vdw_normal, E_plus_vdw_O_def),
         ("without VdW correction", E_no_vdw_normal, E_no_vdw_O_def)]:

        H = E_O_def + 0.5 * EO2 - E_normal  # Normalized By atom of oxygen extracted

        print(
            "{} \nE_O_def   + 1/2.EO2 - E_normal  =  H \n{:.5f} + {:.5f} - 1/2.{:.5f} = {:.5f} " .format(
                VDW,
                E_O_def,
                EO2,
                E_normal,
                H))

        # H = 2*H # normalized by O2 molecules

    # return H without VdW correction
    return(H)


def labile_O_indices(dictStruct, nb_sites=1):
    # returns the indices of the N non-equivalent oxygen positions
    # with the lowest bader charge (most oxydized)

    O_indices = []

    O_site_list = sorted([s for s in dictStruct['equivSiteList']
                          if s['element'] == "O"],
                         key=itemgetter('charge'))

    for i in range(nb_sites):
        O_indices.append(O_site_list[i]['indices'][0])
    return(O_indices)
