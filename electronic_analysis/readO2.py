# readO2.py

"""read and write O2 release enthalpy related quantities"""


import os
import subprocess
from multiprocessing import Pool, cpu_count
from operator import itemgetter

from pymatgen.io.vasp.outputs import Oszicar

import write_job as launch
import electronic_analysis.rundict_utils as read
import rw_utils.generic_plot as generic_plot
import rw_utils.platform_id as platform_id


def O2_computation(rundict_list, bader_done=False):
    vasp_run_poll = []
    with Pool(processes=cpu_count()) as p:
        vasp_run_poll = p.map(get_O2_tag_single, rundict_list)
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
    """    for each structure given in argument,
    remove the most oxidized oxygen
    and create a VASP input folder in the job folder
    option for redundancy to generate a larger supercell"""

    if redundancy is None:
        if input("redundancy ? Y/N")[0] == "Y":
            redundancy = True
        else:
            redundancy = False

    folder_list = []

    for run in vaspRundictList:
        parent_path = os.path.join(run.job_folder, "O_deficient")
        print("removing O from {}".format(run.str_id))

        for O_index in labile_O_indices(run, nb_sites):
            print("removed O index : {}".format(O_index))
            job = launch.Job.from_rundict(run)
            job.structure.remove_sites([O_index])
            job.entry_id = "{}_rem_O{}".format(job.entry_id, O_index)
            folder_list.append(
                job.set_job_folder(parent_path, explicit_jobpath=False))
            job.write_data_input()

            if redundancy:
                job = launch.Job.from_rundict(run)
                job.structure.make_supercell([1, 2, 1])
                job.entry_id = "{}_big".format(job.entry_id)
                folder_list.append(
                    job.set_job_folder(parent_path, explicit_jobpath=False))
                job.write_data_input()

    for folder in folder_list:
        setting_dir = platform_id.setting_dir()
        job_file_name = os.path.join(
            setting_dir, "job_scripts", "vasp_job_double")
        name = "degaz_{}".format(abs(hash(folder)))
        job_string = 'sbatch -J {} --workdir {}  {}'.format(
            name, folder, job_file_name)
        if input(
                "LAUNCH the following [Y/n] :\n{}\n".format(job_string)) == "Y":
            subprocess.call([job_string], shell=True)

    return folder_list


def get_O2_release_enthalpy(rundict):
    """ compare energies of normal and O deficient runs

    check the existence of the O2 deficient vasp result in the run folder
    get the run with the lowest energy
    compare energies of normal and O deficient runs (w/ and wo/ WdW correction)
    return H without VdW correction"""

    run_list = []

    # check the existence of the O2 deficient vasp result in the run folder
    o_deficient_list = [
        os.path.join(
            rundict.job_folder,
            f) for f in os.listdir(
            rundict.job_folder) if f.startswith("O_deficient")]

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

    # print("valid run lst {}".format(run_list))

    # get the run with the lowest energy
    run_list = sorted(run_list, key=lambda x: x.energy)
    o_def_rundict = run_list[0]
    # ============= IN CASE OF VDW CORRECTION =========================
    E_plus_vdw_O_def = o_def_rundict.energy
    E_no_vdw_O_def = Oszicar(o_def_rundict.job_folder +
                             "/OSZICAR").electronic_steps[-1][-1]["E"]

    # print("Energies 0 deficient : base {:.4f} // VDW {:.4f} (corr : {:.6f})".format(
    # E_no_vdw_O_def, E_plus_vdw_O_def, E_plus_vdw_O_def -  E_no_vdw_O_def  )
    # )

    # nbO_final = len(o_def_vasprun.final_structure.indices_from_symbol("O"))

    E_plus_vdw_normal = rundict.energy
    E_no_vdw_normal = Oszicar(rundict.job_folder +
                              "/OSZICAR").electronic_steps[-1][-1]["E"]

    # print("Energies : base {} // VDW {} (corr : {})".format(
    # E_no_vdw_normal, E_plus_vdw_normal, E_plus_vdw_normal -  E_no_vdw_normal
    # )  )

    # nbO_init = len(runDict['structure'].indices_from_symbol("O"))

    # H = E(normal) - E(O_deficient) - 1/2 E(O2)

    # runDict['ediff'] = runDict.energy_per_fu / runDict["nb_cell"]
    # eV  -9.8897568 (no Wdv energy from dft O2 ) + 1.36 (ceder correction)
    EO2 = -8.52
    print("{}\n".format(rundict.name_tag))
    for VDW, E_normal, E_O_def in \
        [("with VdW correction", E_plus_vdw_normal, E_plus_vdw_O_def),
         ("without VdW correction", E_no_vdw_normal, E_no_vdw_O_def)]:

        H = E_O_def + 0.5 * EO2 - E_normal  # Normalized By atom of oxygen extracted

        print(
            "{} \nE_O_def   + 1/2.EO2 - E_normal  =  H \n{:.5f} + {:.5f} - 1/2.{:.5f} = {:.5f} " .format(
                VDW, E_O_def, EO2, E_normal, H))

        # H = 2*H # normalized by O2 molecules

    # return H without VdW correction
    return H


def labile_O_indices(rundict, nb_sites=1):
    """returns the indices of the N non-equivalent oxygen positions
    with the lowest bader charge (most oxydized)"""

    O_indices = []

    O_site_list = sorted([s for s in rundict.equivSiteList
                          if s['element'] == "O"],
                         key=itemgetter('charge'))

    for i in range(nb_sites):
        O_indices.append(O_site_list[i]['indices'][0])
    return(O_indices)
