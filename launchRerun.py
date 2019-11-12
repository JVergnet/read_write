#!/usr/bin/env python3
# launchRerun.py
# Creates VASP imput files for vasp in the folder projectName/jobFformula
# INPUT : the projectDirectory
# OUTPUT : creates rerun of the specified project with specified strategy
# (static, relax or non scf)


import os
import shutil
import sys

from pymatgen.io.vasp.inputs import Kpoints
# from pymatgen import Structure
from pymatgen.io.vasp.sets import MITRelaxSet  # , MPNonSCFSet  # ,MPStaticSet,

import launchDisordered as launch
# import platform_id
import readRun_entries as read
# from pymatgen.io.vasp.outputs import Vasprun
# from pymatgen.io.vasp.inputs import Incar
# from pymatgen.core.periodic_table import Element
# import platform
# import matplotlib
from drawKpoints import drawkpt
import filter_runs as filter_runs

# import read_hull as hull


# if platform.node() in 'bipbip.lsd.univ-montp2.fr':
#     matplotlib.use("Agg")


def less_precise_incar(struct):
    return(dict(ENCUT=600,
                PREC='Normal',
                EDIFFG=-1E-01,
                EDIFF=1E-04 * struct.num_sites,
                # far from minimum : conjugate gradient algorithm (2)
                IBRION=2,
                LCHARG="True",
                SIGMA=0.2,
                ISMEAR=0,
                LELF="False",
                ISYM=0))


def more_precise_incar(struct):
    return(dict(LCHARG="True",
                ENCUT=700,
                PREC='Accurate',
                EDIFFG=-1E-02,
                EDIFF=1E-06 * struct.num_sites,
                # Close to the local minimum : RMM-DIIS (1)
                IBRION=1,
                SIGMA=0.01,
                ISMEAR=0,
                LELF="True",
                ISYM=0))


def ultra_precise_incar():
    return(dict(NSW=100,
                ENCUT=650,
                PREC='Accurate',
                ALGO='Normal',
                EDIFFG=-1E-3,
                EDIFF=1E-8,  # * s.num_sites,
                # Close to the local minimum : RMM-DIIS (1)
                IBRION=1,
                # SIGMA=0.01,
                # ISMEAR=0,
                #  LELF  = "True",
                ISYM=0,
                NEDOS=10001,
                # EMIN=-8,
                # EMAX=8
                ))


def single_point_incar():
    return(dict(NSW=0,
                ENCUT=700,
                EDIFF_PER_ATOM=1e-06,
                EDIFFG=-0.001,
                ISMEAR=-5,
                SIGMA=0.05,
                NELM=300,
                PREC='Accurate'))


def prompt_rerun_type():
    rerun_type, incar_type = (None, None)

    rerun_dict = {'s': 'single_point',
                  'r': 'relaxation',
                  'i': 'identical',
                  'c': 'custom'}
    print("\n"+"".join(["[{}]: {}\n".format(*k)
                        for k in rerun_dict.items()]))

    rerun_choice = input(
        "rerun type ? : ")

    rerun_type = rerun_dict[rerun_choice[0]]
    print("rerun chosen type : {}\n".format(rerun_type))

    if rerun_type in ["single_point", "relaxation"]:

        if rerun_type == "single_point":
            incar_dict = {'s': 'static',
                          'd': 'DOS',
                          'n': 'non_SCF',
                          'f': 'fukui',
                          'o': 'poscar_only',
                          'p': "parcharg"}

        elif rerun_type == "relaxation":
            incar_dict = {'l': 'less_precise ',
                          'm': 'more_precise ',
                          'r': 'rebuild_from_scratch ',
                          'u': 'ultra_precise',
                          'p': 'poscar_only'}

        print("\n"+"".join(["[{}]: {}\n".format(*k)
                            for k in incar_dict.items()]))
        incar_choice = input(
            "incar type ? : ")

        incar_type = incar_dict[incar_choice[0]]
        print("incar chosen type : {}\n".format(incar_type))

    return(rerun_type, incar_type)


def main():
    """Set the parameters ofthe run"""
    try:
        main_dir = sys.argv[1]
    except IndexError:
        main_dir = os.getcwd()
        print("in current folder : {}".format(main_dir))

    rerun_type, incar_type = prompt_rerun_type()

    rerun_select = input(
        "Rerun [a]ll / [n]on-converged only / [c] converged-only ? : ")
    # do not parse vasprun if all job selected
    check_vasprun = 0 if rerun_select in ["a"] else 0.9

    try:
        file_system = input("[j]ob / [p]roject / [s]uper_project ?  :  ")[0]
    except Exception:
        file_system = "p"
    print("filesystem : {}".format(file_system))

    # Create a list of all the valid runs in the selected folders
    run_list = read.collect_valid_runs(main_dir, checkDiff=False,
                                       vasprun_parsing_lvl=check_vasprun,
                                       file_system_choice=file_system)

    converged_jobs = [d for d in run_list if d.status == 3]
    unconverged_jobs = [d for d in run_list if d.status < 3]

    # when reruning all jobs, they are all considered as unconverged
    if rerun_select in ["a", "n"]:
        rerun_list = unconverged_jobs

    elif rerun_select in ["c"]:
        rerun_list = converged_jobs

    rerun_list = filtering_runs(rerun_select, rerun_list)

    print("number of valid jobs to rerun : {}".format(len(rerun_list)))

    if len(rerun_list) == 0:
        print("no valid run")
        return 0

    print("selected runs : \n {}".format(
        [print(rundict.str_id) for rundict in rerun_list]))
    try:
        perturb = float(
            input('Perturb the initial position of atoms ? in Angstrom '))
    except Exception as ex:
        print("No perturbation")
        perturb = 0

    dirname = incar_type if incar_type is not None else rerun_type
    print("current dirname ={}".format(dirname))

    if incar_type == "fukui":
        fukui_nelec = float(
            input("nb elec for the fukui (>0: added, <0 : removed) ? "))
        print("fukui electrons : {}".format(fukui_nelec))

    elif rerun_type == "custom":
        try:
            tmpdir = str(input("Custom directory name ? :"))
            if len(tmpdir) > 0:
                dirname = tmpdir
            print(dirname)
        except Exception:
            print("error, default dirname to {}".format(dirname))

    if file_system in ["p", "s"]:
        dirname_path = read.get_file_name(main_dir, dirname)
    elif file_system in ["j"]:
        dirname_path = rerun_list[0].job_folder

    for rundict in rerun_list:

        # create Job from a RunDict
        job = launch.Job.from_rundict(rundict)
        # s = job.structure

        files_to_copy = []
        incar = {}
        if rerun_type == "identical":
            files_to_copy += ['INCAR', 'POTCAR', 'KPOINTS', 'CONTCAR']
            # quick and dirty copy
            print("identical set generated")

        if rerun_type == "relaxation":
            if incar_type == "poscar_only":
                pass
            elif incar_type == "rebuild_from_scratch":
                pass

            elif incar_type == "less_precise":
                incar.update(less_precise_incar(job.structure))
                print(" less precise set generated")

            elif incar_type == "more_precise":
                incar.update(more_precise_incar(job.structure))
                print("more precise set generated")

            elif incar_type == "ultra_precise":
                incar.update(ultra_precise_incar())
                print("ultra precise set generated")

        elif rerun_type == "custom":

            # incar['LDAUU']={'O': 6}

            # incar['ENCUT']= 600
            # incar['PREC'] = 'Normal'
            # incar['EDIFFG'] = -1E-02
            # incar['EDIFF'] = 1E-04
            # incar["NSW"]=100
            # ibrion far from min: conj. grad. algo (2)
            # incar["LCHARG"] = "True"
            # incar["SIGMA"] = 0.05
            # incar["ISMEAR"] = 0
            # incar["LELF"] = "False"
            # incar["ISYM"] = 0
            # incar_copy['LDAUU']={'O':2, 'Mn' : 3.9}
            # incar_copy['LDAUL']={'O':1, 'Mn' : 2}
            # incar_copy['SYSTEM'] = incar_copy['SYSTEM'].replace("cu", "ni")
            # incar_copy['SYSTEM'] = incar_copy['SYSTEM'].replace("Cu", "Ni")
            # incar_copy['ISIF'] = 3

            # incar_copy["NELMIN"] = 1
            # incar_copy['EDIFF'] = 1E-04
            # incar['LDAUU'] = {'Ti': 3.9}

            # incar["EMAX"] = 8
            # incar["EMIN"] = -5
            # incar["NEDOS"] = 5001

            # incar["IBRION"] = 1
            # incar["ALGO"] = "Normal"
            # incar["NCORE"] = 8
            # incar["KPAR"] = 2
            # incar["NUPDOWN"] = 0

            # incar['ICHARG'] = 0

            # HSE06
            # incar["NSW"] = 0
            # incar['LHFCALC'] = "TRUE"
            # incar['HFSCREEN'] = 0.2

            # PARCHG
            # incar["LWAVE"] = "True"

                # "ISMEAR": 0,
                # 'EDIFF': 1E-02,

                # "NELMDL": -10,
                # "PREC": "Accurate",
                # "ADDGRID": "True",
                # "IBRION": 1,
                # "NELM": 150,
                # "SIGMA": 0.1
            incar.update({
                "ALGO": "Normal",
                "EDIFF": 0.001,
                'ISMEAR': 0,
                "SIGMA": 0.1,
                "LCHARG": "True",
                "LAECHG": "True"
            })
            kpt = Kpoints.gamma_automatic(kpts=(1, 1, 1), shift=(0, 0, 0))
            job.user_kpoint = kpt
            print("yolo!!")

            print("MODIFIED PARAMETERS ========", incar)
            print("{} set generated".format(dirname))

        elif rerun_type == "single_point":
            if incar_type == "fukui":
                input_set = MITRelaxSet(job.structure)
                incar = rundict.parameters["incar"]
                incar["NELECT"] = input_set.nelect + fukui_nelec
                incar["NSW"] = 0
                print("fukui correction added :",
                      "\nNELECT read {} ==> wrote {}".format(
                          input_set.nelect, input_set.nelect + fukui_nelec))
            elif incar_type == "parcharg":
                efermi = rundict.data['efermi']
                print(efermi)
                incar["LPARD"] = "True"
                below_fermi = float(input("Emin (Efermi=0) ?"))
                above_fermi = float(input("Emax (Efermi=0) ?"))
                incar["EINT"] = "{} {}".format(efermi+below_fermi,
                                               efermi+above_fermi)
                dirname += "_{}_{}".format(below_fermi,
                                           above_fermi)
                files_to_copy.append("WAVECAR")

            elif incar_type in ["static", "DOS"]:
                incar.update(single_point_incar())
                if incar_type == "DOS":
                    incar['EMIN'] = -5
                    incar['EMAX'] = 5
                    incar["NEDOS"] = 2001
                    # folder = prev_folder + "/DOS"
                    # os.mkdir(folder)
                    kpt_settings = {'reciprocal_density': 1000}
                else:
                    kpt_settings = {'reciprocal_density': 300}

                job.user_kpoint = kpt_settings

            elif incar_type == "non_SCF":
                files_to_copy += ["CHGCAR", "CHG", "linear_KPOINTS"]
                incar.update({"IBRION": -1,
                              "LCHARG": False,
                              "LORBIT": 11,
                              "LWAVE": False,
                              "NSW": 0,
                              "ISYM": 0,
                              "ICHARG": 11,
                              "ISMEAR": 0,
                              "SIGMA": 0.01
                              })
                for k in ["NELMDL", "MAGMOM"]:
                    job.user_incar.pop(k, None)
                # job.set_job_folder(rerun_dir)
                kpt = drawkpt(rundict.structure)
                kpt.write_file(os.path.join(
                    job.old_folder, "linear_KPOINTS"))

        job.old_folder = job.job_folder
        if file_system == "j":
            job.set_job_folder(read.get_file_name(dirname_path, dirname),
                               explicit_jobpath=True)
        else:
            if file_system == "p":
                job.set_job_folder(dirname_path,
                                   explicit_jobpath=False)
            if file_system == "s":
                # print(rundict.stacking)
                job.set_job_folder(os.path.join(dirname_path,
                                                rundict.stacking),
                                   explicit_jobpath=False)

        if incar.get('EDIFF', None) is not None:
            incar['EDIFF'] = '{:0.1E}'.format(incar['EDIFF'])

        for k in ["MAGMOM", "EINT", "LPARD", "SIGMA"]:
            job.user_incar.pop(k, None)

        job.user_incar.update(incar)
        print("INCAR \n", incar)
        print("JOB INCAR\n", job.user_incar)

        job.structure.perturb(perturb)
        print("explicit jobpath", job.explicit_jobpath)

        if rerun_type == "identical":
            os.mkdir(job.job_folder)
        else:
            job.write_data_input()

        for f_name in files_to_copy:
            try:
                shutil.copy2('{0.old_folder}/{1}'.format(job, f_name),
                             '{0.job_folder}/{1}'.format(job, f_name))
            except Exception as ex:
                print("error when copying", f_name, ex)

    if input("[r]emove unconverged folders ? ") == "r":
        for rundict in unconverged_jobs:
            unconv_dir = rundict.old_folder
            # if input("remove {0} ? Y / N ".format(unconv_dir))=="Y" :
            shutil.rmtree(unconv_dir)
            print("{} deleted ".format(unconv_dir))


def filtering_runs(rerun_select, rerun_list):
    selected_runs = rerun_list
    if input("apply further selection on runs ? : Y / n ") == "Y":
        if rerun_select in ["c"]:
            sieve_lvl = filter_runs.select_sieve_level()
            selected_runs = filter_runs.hull_filtering(
                sieve_lvl, selected_runs)
            print("number of selected runs : {}".format(
                len(selected_runs)))

        if input("folder by folder ? : Y / n ") == "Y":
            selected_runs = filter_runs.idv_filtering(selected_runs)
            print("nb of structures : {0} ".format(len(rerun_list)))
    return rerun_list


if __name__ == '__main__':
    main()
