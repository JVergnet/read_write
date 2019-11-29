# launch_job.py
"""Convenience method to launch all the jobs in a prepared project folder
INPUT : path of the folder where the prepared vasp inputs are
"""

import subprocess
import os
import sys
from utils.platform_id import setting_dir


def main(init_dir):

    params = {}

    if len(sys.argv) > 1:
        project_folder = sys.argv[1]
    else:
        project_folder = os.getcwd()

    name = input("Name of the job ? : \n")
    print("Base name : {} \n".format(name))
    params['name'] = "-J {}".format(name)

    num_nodes = 1
    try:
        num = int(input('num nodes  / job ?'))
        if num > 0 and num <= 8:
            num_nodes = num
    except Exception:
        pass
    params['num_nodes'] = "-N {}".format(num_nodes)

    params['QOS'] = ""
    try:
        if input('QOS ? [BG] / [N]ormal') == "BG":
            params['QOS'] = "--qos=bg --requeue"
    except Exception:
        pass

    params['jobFileName'] = os.path.join(init_dir, "vasp_job")

    job_choice = input("job type : [s]ingle_run / [d]ouble_run ?\n")
    if job_choice == "d":
        params['jobFileName'] += "_double"
    else:
        params['jobFileName'] += "_single"

    print("\njobfile : ", params['jobFileName'], "\n")

    try:
        if input('gamma only ?') == "Y":
            params['jobFileName'] += "_gamma"
    except Exception:
        pass

    queue_strategy = input('queue strategy ? [O]ne by one / [A]rray \n')

    if queue_strategy == "A":

        params['jobFileName'] += "_array"
        # os.chdir(projectFolder)

        # Create a list of all the subfolder of the run
        params['workdir'] = "--workdir {}".format(project_folder)
        sub_dir = project_folder
        subdir_list = [os.path.join(sub_dir, o) for o in os.listdir(sub_dir)
                       if os.path.isdir(os.path.join(sub_dir, o))]
        print("Nb jobs to run : {}".format(len(subdir_list)))
        max_job = int(input("Nb of parallel job(max) ? [0-10] \n"))
        if max_job < 0:
            print(" <<< SKY IS THE LIMIT >>> My man !  !")
            max_job = 99
        elif max_job > 10:
            print("too many cores !")
            exit(1)

        array_size = len(subdir_list)
        params['array'] = "--array=1-{}%{}".format(array_size, max_job)
        subprocess.call(['sbatch {name} {num_nodes} {array} {workdir} {QOS} {jobFileName}'
                         .format(**params)],
                        shell=True)

    else:
        # Walk trhough the folders to find valid vasp inputs
        sub_dir = project_folder
        subdir_list = [os.path.join(sub_dir, o) for o in os.listdir(sub_dir)
                       if os.path.isdir(os.path.join(sub_dir, o))]
        print("Individually chose to run the following folders ? \n",
              "select [OK / NO ] \n")

        for job_folder in subdir_list:
            params['workdir'] = "--workdir {}".format(job_folder)

            run = False
            try:
                choice = input("launch {}? :".format(job_folder))
                print(choice)
                if choice in ["OK", "O", "o", "Ok"]:
                    run = True
            except Exception:
                print("invalid option")

            if run:
                subprocess.call([
                    'sbatch {name} {num_nodes} {workdir} {QOS} {jobFileName}'
                    .format(**params)], shell=True)

                print("\n")
            else:
                long_name = job_folder.split("/")[-1].split("__")
                short_name = long_name[0][:4] + "_" + long_name[-1][:3]
                print(
                    " No job launched for {0} in {1} :-( \n".format(short_name,
                                                                    params['jobFileName']))

        os.chdir(project_folder)


if __name__ == '__main__':

    SETTING_DIR = os.path.join(setting_dir(), "job_scripts")
    main(SETTING_DIR)
