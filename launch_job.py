# launch_job.py
"""Convenience method to launch all the jobs in a prepared project folder
INPUT : path of the folder where the prepared vasp inputs are

create a SlurJob to get Slurm Instructions
create a VaspExecutable to select Vasp Executable file
issue shell command accordingly
"""

import subprocess
import os
import sys
from rw_utils import platform_id


def main(project_folder):

    slurm_job = SlurmJob(project_folder)
    slurm_job.launch_job()


class SlurmJob(object):

    def __init__(self, project_folder):

        def define_values(self):
            self.project_folder = project_folder
            self.workdir = None
            self.array = None
            self.name = get_name_str()
            self.num_nodes = get_num_node_str()
            self.QOS = get_qos_str()

            self.vasp_exec = VaspExecutable()
            self.vasp_executable_path = self.vasp_exec.explicit_path()

        def get_name_str():
            try:
                name = input("Name of the job ? : \n")
                assert len(name) > 0, "invalid name"
                name = name.replace(" ", "_")
            except Exception:
                name = "vasp_job"
            print("Base name : {} \n".format(name))
            return "-J {}".format(name)

        def get_qos_str():
            qos_str = ""
            try:
                if input('QOS ? [BG] / [N]ormal : ')[0] in ["B", "b"]:
                    qos_str = "--qos=bg --requeue"
            except Exception:
                pass
            return qos_str

        def get_num_node_str():
            try:
                num = int(input('num nodes / job ? : '))
                assert num > 0 and num <= 8, "invalid node nb"
            except Exception:
                print("default to 1")
                num = 1
            return "-N {}".format(num)

        # ===========================
        # execution of all definitions
        define_values(self)

    def launch_job(self):
        def launch(self):
            if self.vasp_exec.array:
                launch_array_run(self)
            else:
                launch_individual_runs(self)
            os.chdir(self.project_folder)

        def launch_array_run(self):
            "Create a list of all the subfolder of the run"
            workdir = "--workdir {}".format(self.project_folder)
            array = get_array_str(self.project_folder)
            issue_shell_command(self, workdir, array)

        def get_array_str(project_folder):
            sub_dir = project_folder
            subdir_list = [os.path.join(sub_dir, o) for o in os.listdir(sub_dir)
                           if os.path.isdir(os.path.join(sub_dir, o))]
            print("Nb jobs to run : {}".format(len(subdir_list)))
            try:
                max_job = int(input("Nb of parallel job(max) ? [0-10] :  "))
                assert max_job < 10, "too many cores !"
                if max_job == -1:
                    print(" <<< SKY IS THE LIMIT >>> My man !  !")
                    max_job = 99
                assert max_job > 0, "Not enough cores !"
            except Exception as ex:
                print(ex, "default to 1")
                max_job = 1
            array_size = len(subdir_list)
            array_str = "--array=1-{}%{}".format(array_size, max_job)
            return array_str

        def launch_individual_runs(self):
            "Walk trhough the folders to find valid vasp inputs"

            print("Individually chose to run the following folders ? \n",
                  "select [OK / NO ] \n")
            subdir_list = [os.path.join(self.project_folder, o)
                           for o in os.listdir(self.project_folder)
                           if os.path.isdir(os.path.join(self.project_folder, o))]

            for job_folder in subdir_list:
                if ask_for_launch(job_folder):
                    workdir = "--workdir {}".format(job_folder)
                    issue_shell_command(self, workdir)
                else:
                    notify_no_launch(job_folder)

        def issue_shell_command(self, workdir, array=""):
            shell_command = \
                'sbatch {0.name} {0.num_nodes} {0.QOS} {1} {2}  {0.vasp_executable_path}'.format(
                    self, array, workdir)
            print(shell_command)
            # subprocess.call([shell_command], shell=True)

        def ask_for_launch(job_folder):
            launch = False
            try:
                choice = input("launch {}? :".format(job_folder))
                print(choice)
                if choice in ["OK", "O", "o", "Ok"]:
                    launch = True
            except Exception:
                print("invalid option")
            return launch

        def notify_no_launch(job_folder):
            long_name = job_folder.split("/")[-1].split("__")
            short_name = long_name[0][:4] + "_" + long_name[-1][:3]
            print(
                " No job launched for {0} :-( \n".format(short_name))

        # ==============
        # launch the job
        launch(self)


class VaspExecutable(object):

    def __init__(self, set_dir=None, base_job_name="vasp_job"):
        "defines the properties of the vasp job script"

        self.set_base_name(set_dir, base_job_name)

        # to avoid inconsistent exec name
        # the properties must be defined in this order
        self.nb_runs = "single"
        self.array = False
        self.gamma = False
        self.nc = False
        self.define_properties()

    def set_base_name(self, set_dir, base_job_name):
        if set_dir is not None:
            self.setting_dir = set_dir
        else:
            self.setting_dir = os.path.join(
                platform_id.setting_dir(), "job_scripts")
        self.base_job_name = base_job_name
        self.base_name = os.path.join(self.setting_dir, self.base_job_name)

    def define_properties(self):
        "define run properties in a consistent order"

        def def_all(self):
            "function called at the end of the definitions "
            get_nb_run(self)
            get_array_bool(self)

            if self.array and self.nb_runs == "single":
                get_gamma_bool(self)

            if not self.array and self.nb_runs == "single":
                get_nc_bool(self)

        def get_nb_run(self):
            try:
                if input("job type : [s]ingle_run / [d]ouble_run ? ")[0] in ["d", "D", "2"]:
                    self.nb_runs = "double"
            except Exception:
                pass

        def get_array_bool(self):

            print("queue strategy ? ")
            try:
                if input('[O]ne by one / [A]rray : ')[0] in ["A", "a"]:
                    self.array = True
            except Exception:
                pass

        def get_gamma_bool(self):
            """we only use gamma computation for
            - large nb of run (array)
            - without change in supercell (single run)"""
            try:
                if input('gamma only ? [Y]/n)')[0] == "Y":
                    self.gamma = True
            except Exception:
                pass

        def get_nc_bool(self):
            """we only use non-colinear computation for
            - unique run (no array)
            - no change in supercell (single run)"""
            try:
                if input('non-colinear ? (Y]/n)')[0] == "Y":
                    self.nc = True
            except Exception:
                pass

        # ============================================
        # calling all function above in specific order
        def_all(self)

    def explicit_path(self):

        str_list = [self.base_name, self.nb_runs]

        for prop in ['array', "gamma", "nc"]:
            if getattr(self, prop):
                str_list.append(prop)
        print("_".join(str_list))

        return "_".join(str_list)


if __name__ == '__main__':

    if len(sys.argv) > 1:
        PROJECT_FOLDER = sys.argv[1]
    else:
        PROJECT_FOLDER = os.getcwd()

    main(PROJECT_FOLDER)
