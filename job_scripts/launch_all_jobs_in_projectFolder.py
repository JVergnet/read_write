# launch_all_jobs_in_projectFolder.py
# Convenient method to launch all the jobs in a prepared project folder
# INPUT : path of the folder where the prepared vasp inputs are
# output :

import subprocess
import os
import sys

def main(settingDir):
    
    
    params = {}

    if len(sys.argv) > 1:
        projectFolder = sys.argv[1]
    else:
        projectFolder = os.getcwd()


    name = input("Name of the job ? : \n")
    print("Base name : {} \n".format(name))
    params['name'] = "-J {}".format(name)

    num_nodes = 1
    try:
        num = eval(input('num nodes  / job ?'))
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

    params['jobFileName'] = os.path.join(settingDir, "vasp_job")

    jobChoice = input("job type : [s]ingle_run / [d]ouble_run ?\n")
    if jobChoice == "d":
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
        params['workdir'] = "--workdir {}".format(projectFolder)
        d = projectFolder
        subDirList = [os.path.join(d, o) for o in os.listdir(d)
                    if os.path.isdir(os.path.join(d, o))]
        print("Nb jobs to run : {}".format(len(subDirList)))
        maxJob = eval(input("Nb of parallel job(max) ? [0-10] \n"))
        if maxJob < 0:
            print(" <<< SKY IS THE LIMIT >>> My man !  !")
            maxJob = 99
        elif maxJob > 10:
            print("too many cores !")
            exit(1)

        array_size = len(subDirList)
        params['array'] = "--array=1-{}%{}".format(array_size, maxJob)
        subprocess.call(['sbatch {name} {num_nodes} {array} {workdir} {QOS} {jobFileName}'
                        .format(**params)],
                        shell=True)

    else:
        # Walk trhough the folders to find valid vasp inputs
        d = projectFolder
        subDirList = [os.path.join(d, o) for o in os.listdir(d)
                    if os.path.isdir(os.path.join(d, o))]
        print("Individually chose to run the following folders ? \n",
            "select [OK / NO ] \n")

        for i, jobFolder in enumerate(subDirList):
            params['workdir'] = "--workdir {}".format(jobFolder)

            run = False
            try:
                choice = input("launch {}? :".format(jobFolder))
                print(choice)
                if choice in ["OK", "O", "o", "Ok"]:
                    run = True
            except Exception:
                print("invalid option")

            if run:
                subprocess.call([
                    'sbatch {name} {num_nodes} {workdir} {QOS} {jobFileName}'
                    .format(**params)],
                    shell=True)

                print("\n")
            else:
                longName = jobFolder.split("/")[-1].split("__")
                shortName = longName[0][:4] + "_" + longName[-1][:3]
                print(
                    " No job launched for {0} in {1} :-( \n".format(shortName,
                                                                    params['jobFileName']))

        os.chdir(projectFolder)

if __name__ == '__main__':
    from utils.platform_id import setting_dir
    settingDir = os.path.join(setting_dir(), "job_scripts")
    main(settingDir)