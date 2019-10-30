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
import platform_id
import readRun_entries as read
# from pymatgen.io.vasp.outputs import Vasprun
# from pymatgen.io.vasp.inputs import Incar
# from pymatgen.core.periodic_table import Element
# import platform
# import matplotlib
from drawKpoints import drawkpt

# import read_hull as hull


# if platform.node() in 'bipbip.lsd.univ-montp2.fr':
#     matplotlib.use("Agg")


def less_precise_incar(struct):
    return(dict(ENCUT=600,
                PREC='Normal',
                EDIFFG=-1E-01,
                EDIFF=1E-04 * struct.num_sites,
                NSW=100,
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
                NSW=100,
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
                          'p': 'poscar_only'}

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

    # Set the parameters ofthe run
    setting_dir = platform_id.setting_dir()

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
                                       vaspRun_parsing_lvl=check_vasprun,
                                       file_system_choice=file_system)

    converged_jobs = [d for d in run_list if d.status == 3]
    unconverged_jobs = [d for d in run_list if d.status < 3]

    # when reruning all jobs, they are all considered as unconverged
    if rerun_select in ["a", "n"]:
        rerun_list = unconverged_jobs

    elif rerun_select in ["c"]:
        rerun_list = converged_jobs

    if input("apply further selection on runs ? : Y / n ") == "Y":
        if rerun_select in ["c"] and input(
                "convex hull filtering ? : Y / n ") == "Y":
            rerun_list = read.restrict_run_list(rerun_list)
            print("selected runs : \n {}".format(
                [print(rundict.id) for rundict in rerun_list]))

        if input("folder by folder ? : Y / n ") == "Y":
            rerun_list_tmp = rerun_list
            rerun_list = []
            for run in rerun_list_tmp:
                if input("include {} : Y / n ".format(run.id)) == "Y":
                    rerun_list.append(run)
        print("nb of structures : {0} ".format(len(rerun_list)))

    print("number of valid jobs to rerun : {}".format(len(rerun_list)))

    if len(rerun_list) == 0:
        print("no valid run")
        return(0)

    print("selected runs : \n {}".format(
        [print(rundict.id) for rundict in rerun_list]))
    try:
        perturb = float(
            input('Perturb the initial position of atoms ? in Angstrom '))
    except Exception as ex:
        print("No perturbation")
        perturb = 0

    dirname = incar_type if incar_type is not None else rerun_type
    print("current dirname ={}".format(dirname))

    if incar_type == "fukui":
        fukui_nelec = 1 if (input("[+]1 or [-]1 elec ? ") == "+") else -1
        print("fukui n electrons : {}".format(fukui_nelec))

    elif rerun_type == "custom":
        try:
            tmpdir = str(input("Custom directory name ? :"))
            if len(tmpdir) > 0:
                dirname = tmpdir
            print(dirname)
        except Exception:
            print("error, default dirname to {}".format(dirname))
        parcharg = True if input("parcharg?") == "Y" else False

    if file_system in ["p", "s"]:
        dirname_path = read.get_file_name(main_dir, dirname)

    for rundict in rerun_list:

        # create Job from a RunDict
        job = launch.Job.from_rundict(rundict)
        # s = job.structure

        if file_system == "j":
            rerun_dir = read.get_file_name(rundict.job_folder, dirname)
            job.explicit_jobpath = True
        else:
            job.explicit_jobpath = False
            if file_system == "p":
                rerun_dir = dirname_path
            if file_system == "s":
                print(rundict.stacking)
                rerun_dir = os.path.join(dirname_path, rundict.stacking)
        job.oldFolder = job.job_folder
        job.set_job_folder(rerun_dir)

        if rerun_type == "identical":
            # quick and dirty copy
            os.makedirs(job.job_folder, exist_ok=True)
            for f_name in ['CONTCAR', 'POSCAR']:
                try:
                    shutil.copy2('{0.oldFolder}/{1}'.format(job, f_name),
                                 '{0.job_folder}/POSCAR'.format(job))
                    break
                except Exception as ex:
                    pass

            for f_name in ['INCAR', 'POTCAR', 'KPOINTS']:
                try:
                    shutil.copy2('{0.oldFolder}/{1}'.format(job, f_name),
                                 '{0.job_folder}/{1}'.format(job, f_name))
                except Exception as ex:
                    print(ex)
            print("identical set generated")

            continue

        incar = {}
        files_to_copy = []

        if rerun_type == "poscar_only":
            pass
            # incar_default = launch.default_incar()
            # incar_default["SYSTEM"] = runDict["job_name"]
            # s = runDict['structure']
            # s.remove_site_property("selective_dynamics")
            # inputSet = MITRelaxSet(s,
            #                        user_kpoints_settings={
            #                            'reciprocal_density': 100},
            #                        force_gamma=True,
            #                        user_incar_settings=incar_default)
            # inputSet.write_input(folder)

        elif rerun_type == "relaxation":

            if incar_type == "less_precise":
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
            # os.chdir(folder)
            # incar.write_file("INCAR")
            # settingCopy=dict(incar_setting)
            # incar.pop("MAGMOM",None)
            # settingCopy['SYSTEM']=jobName
            # incar_copy['LDAUU']={'O':2, 'Mn' : 3.9}
            # incar_copy['LDAUL']={'O':1, 'Mn' : 2}
            # incar_copy['SYSTEM'] = incar_copy['SYSTEM'].replace("cu", "ni")
            # incar_copy['SYSTEM'] = incar_copy['SYSTEM'].replace("Cu", "Ni")
            # incar_copy['ISIF'] = 3

            # incar_copy["NELMIN"] = 1
            # incar_copy['EDIFF'] = 1E-04
            # incar['LDAUU'] = {'Ti': 3.9}
            # print("NEW INCAR========", incar_copy)
            # incar["LELF"] = "False"
            # incar["NELM"] = 60
            # incar["EMAX"] = 8
            # incar["EMIN"] = -5
            # incar["NEDOS"] = 5001

            # incar["IBRION"] = 1
            # incar["ALGO"] = "Normal"
            # incar["NCORE"] = 8
            # incar["KPAR"] = 2
            # incar["NUPDOWN"] = 0
            # incar["NELM"] = 60
            # incar['PREC'] = "Accurate"
            # incar['ICHARG'] = 0

            # HSE06
            # incar["NSW"] = 0
            # incar['LHFCALC'] = "TRUE"
            # incar['HFSCREEN'] = 0.2

            # PARCHG
            # incar["LWAVE"] = "True"
            incar['EDIFF'] = 1E-02
            incar["ISMEAR"] = -5
            incar["LCHARG"] = "False"
            incar["LAECHG"] = "False"
            incar["NELMDL"] = -10
            incar["PREC"] = "Accurate"
            incar["ADDGRID"] = "True"
            incar["IBRION"] = 1
            incar["NELM"] = 150
            kpt = Kpoints.gamma_automatic(kpts=(3, 3, 3), shift=(0, 0, 0))
            job.user_kpoint = kpt
            print("yolo!!")

            if parcharg:
                efermi = rundict.data['efermi']
                print(efermi)
                incar["LPARD"] = "True"
                below_fermi = float(input("Emin (Efermi=0) ?"))
                above_fermi = float(input("Emax (Efermi=0) ?"))
                incar["EINT"] = "{} {}".format(efermi+below_fermi,
                                               efermi+above_fermi)
                # incar["EINT"] = "{} {}".format(efermi, efermi+0.5)

                print("MODIFIED PARAMETERS ========", incar)
                # settingCopy['LDAUU']={'Mn':U}
                # s.replace_species({Element("Cu"): Element("Ni")})
                #     job.write_vasp_input(

                files_to_copy.append("WAVECAR")

            # if fileSystem == "j":
            #     rerunDir = read.get_file_name(runDict.job_folder, dirname)
            # else:
            #     rerunDir = read.get_file_name(mainDir, dirname)
            #     if fileSystem == "s":
            #         rerunDir = os.path.join(rerunDir, runDict.stacking)
            # job.set_job_folder(rerunDir)
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

                # input_set = MITRelaxSet(
                #     s, force_gamma=True,
                #     user_kpoints_settings=kpt_settings,
                #     user_incar_settings=incar_copy)
                # print(input_set.incar)
                # folder_name = input_set.incar['SYSTEM'].replace(
                #     ' ', '-').replace('.', '_')
                # full_folder_name = os.path.join(new_parent_dir, folder_name)
                # input_set.write_input(full_folder_name)

                # print(full_folder_name, " static set generated")  #

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
                job.set_job_folder(rerun_dir)
                kpt = drawkpt(rundict.structure)
                kpt.write_file(os.path.join(
                    job.oldFolder, "linear_KPOINTS"))
                # inputSet = MPNonSCFSet.from_prev_calc(
                #     prev_folder,
                #     reciprocal_density=100,
                #     kpoints_line_density=20,
                #     user_incar_settings=incar)
                # print(folder, " non_SCF set generated")
                # inputSet.write_input(folder)

        if incar.get('EDIFF', None) is not None:
            incar['EDIFF'] = '{:0.1E}'.format(incar['EDIFF'])

        for k in ["MAGMOM", "EINT", "LPARD", "SIGMA"]:
            job.user_incar.pop(k, None)

        job.user_incar.update(incar)
        print("INCAR \n", incar)
        print("JOB INCAR\n", job.user_incar)

        job.structure.perturb(perturb)
        print("explicit jobpath", job.explicit_jobpath)
        job.write_data_input(rerun_dir)

        for f_name in files_to_copy:
            try:
                shutil.copy2('{0.oldFolder}/{1}'.format(job, f_name),
                             '{0.job_folder}/{1}'.format(job, f_name))
            except Exception as ex:
                print("error when copying", f_name, ex)

    if input("[r]emove unconverged folders ? ") == "r":
        for rundict in unconverged_jobs:
            unconv_dir = rundict.oldFolder
            # if input("remove {0} ? Y / N ".format(unconv_dir))=="Y" :
            shutil.rmtree(unconv_dir)
            print("{0} deleted ".format(unconv_dir))


if __name__ == '__main__':
    main()
