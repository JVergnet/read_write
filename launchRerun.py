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
from pymatgen.io.vasp.sets import MITRelaxSet, MPNonSCFSet  # ,MPStaticSet,

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


def less_precise_incar():
    return(dict(ENCUT=600,
                PREC='Normal',
                EDIFFG=-1E-01,
                EDIFF=1E-04 * s.num_sites,
                NSW=100,
                # far from minimum : conjugate gradient algorithm (2)
                IBRION=2,
                LCHARG="True",
                SIGMA=0.2,
                ISMEAR=0,
                LELF="False",
                ISYM=0))


def more_precise_incar():
    return(dict(LCHARG="True",
                ENCUT=700,
                PREC='Accurate',
                EDIFFG=-1E-02,
                EDIFF=1E-06 * s.num_sites,
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
    settingDir = platform_id.setting_dir()

    try:
        mainDir = sys.argv[1]
    except IndexError:
        mainDir = os.getcwd()
        print("in current folder : {}".format(mainDir))

    rerun_type, incar_type = prompt_rerun_type()

    rerun_select = input(
        "Rerun [a]ll / [n]on-converged only / [c] converged-only ? : ")
    # do not parse vasprun if all job selected
    check_vaspRun = 0 if rerun_select in ["a"] else 0.9

    try:
        fileSystem = input("[j]ob / [p]roject / [s]uper_project ?  :  ")[0]
    except Exception:
        fileSystem = "p"
    print("filesystem : {}".format(fileSystem))

    # Create a list of all the valid runs in the selected folders
    run_list = read.collect_valid_runs(mainDir, checkDiff=False,
                                       vaspRun_parsing_lvl=check_vaspRun,
                                       file_system_choice=fileSystem)

    converged_jobs = [d for d in run_list if d.status == 3]
    unconverged_jobs = [d for d in run_list if d.status < 3]

    # when reruning all jobs, they are all considered as unconverged
    if rerun_select in ["a", "n"]:
        rerunList = unconverged_jobs

    elif rerun_select in ["c"]:
        rerunList = converged_jobs

    if input("apply further selection on runs ? : Y / n ") == "Y":
        if rerun_select in ["c"] and input(
                "convex hull filtering ? : Y / n ") == "Y":
            rerunList = read.generate_tags(rerunList, minimal=True)
            rerunList = read.restrict_run_list(rerunList)
            print("selected runs : \n {}".format(
                [print(runDict.id) for runDict in rerunList]))

        if input("folder by folder ? : Y / n ") == "Y":
            rerun_list_tmp = rerunList
            rerunList = []
            for run in rerun_list_tmp:
                if input("include {} : Y / n ".format(run.id)) == "Y":
                    rerunList.append(run)
        print("nb of structures : {0} ".format(len(rerunList)))

    print("number of valid jobs to rerun : {}".format(len(rerunList)))

    if len(rerunList) == 0:
        print("no valid run")
        return(0)

    print("selected runs : \n {}".format(
        [print(runDict.id) for runDict in rerunList]))
    try:
        perturb = eval(
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

    if fileSystem in ["p", "s"]:
        dirname_path = read.get_file_name(mainDir, dirname)

    for runDict in rerunList:

        # create Job from a RunDict
        job = launch.Job.from_runDict(runDict)
        # s = job.structure

        if fileSystem == "j":
            rerunDir = read.get_file_name(runDict.jobFolder, dirname)
            job.explicit_jobpath = True
        else:
            job.explicit_jobpath = False
            if fileSystem == "p":
                rerunDir = dirname_path
            if fileSystem == "s":
                print(runDict.stacking)
                rerunDir = os.path.join(dirname_path, runDict.stacking)
        job.oldFolder = job.jobFolder
        job.set_jobFolder(rerunDir)

        if rerun_type == "identical":
            # quick and dirty copy
            os.makedirs(job.jobFolder, exist_ok=True)
            for F in ['CONTCAR', 'POSCAR']:
                try:
                    shutil.copy2('{0.oldFolder}/{1}'.format(job, F),
                                 '{0.jobFolder}/POSCAR'.format(job))
                    break
                except Exception as ex:
                    pass

            for F in ['INCAR', 'POTCAR', 'KPOINTS']:
                try:
                    shutil.copy2('{0.oldFolder}/{1}'.format(job, F),
                                 '{0.jobFolder}/{1}'.format(job, F))
                except Exception as ex:
                    print(ex)
            print("identical set generated")

            continue

        else:
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
                incar.update(less_precise_incar())
                print(" less precise set generated")

            elif incar_type == "more_precise":
                incar.update(more_precise_incar())
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
            incar["ISMEAR"] = -5
            incar["LCHARG"] = "False"
            incar["LAECHG"] = "False"
            incar["NELMDL"] = -10
            incar["PREC"] = "Accurate"
            incar["ADDGRID"] = "True"
            incar["IBRION"] = 2
            incar["NELM"] = 150
            kpt = Kpoints.gamma_automatic(kpts=(3, 3, 3), shift=(0, 0, 0))
            job.user_kpoint = kpt
            print("yolo!!")

            if parcharg:
                efermi = runDict.data['efermi']
                print(efermi)
                incar["LPARD"] = "True"
                below_fermi = eval(input("Emin (Efermi=0) ?"))
                above_fermi = eval(input("Emax (Efermi=0) ?"))
                incar["EINT"] = "{} {}".format(efermi+below_fermi,
                                               efermi+above_fermi)
                # incar["EINT"] = "{} {}".format(efermi, efermi+0.5)

                print("MODIFIED PARAMETERS ========", incar)
                # settingCopy['LDAUU']={'Mn':U}
                # s.replace_species({Element("Cu"): Element("Ni")})
                #     job.write_vasp_input(

                files_to_copy.append("WAVECAR")

            # if fileSystem == "j":
            #     rerunDir = read.get_file_name(runDict.jobFolder, dirname)
            # else:
            #     rerunDir = read.get_file_name(mainDir, dirname)
            #     if fileSystem == "s":
            #         rerunDir = os.path.join(rerunDir, runDict.stacking)
            # job.set_jobFolder(rerunDir)
            print("{} set generated".format(dirname))

        elif rerun_type == "single_point":
            if incar_type == "fukui":
                inputSet = MITRelaxSet(job.structure)
                incar = runDict.parameters["incar"]
                incar["NELECT"] = inputSet.nelect + fukui_nelec
                incar["NSW"] = 0
                print("fukui correction added :",
                      "\nNELECT read {} ==> wrote {}".format(
                          inputSet.nelect, inputSet.nelect + fukui_nelec))

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
                job.set_jobFolder(rerunDir)
                kpt = drawkpt(runDict.structure)
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
        job.write_data_input(rerunDir)

        for F in files_to_copy:
            try:
                shutil.copy2('{0.oldFolder}/{1}'.format(job, F),
                             '{0.jobFolder}/{1}'.format(job, F))
            except Exception as ex:
                print("error when copying", F, ex)

    if input("[r]emove unconverged folders ? ") == "r":
        for runDict in unconverged_jobs:
            unconv_dir = runDict.oldFolder
            # if input("remove {0} ? Y / N ".format(unconv_dir))=="Y" :
            shutil.rmtree(unconv_dir)
            print("{0} deleted ".format(unconv_dir))


if __name__ == '__main__':
    main()
