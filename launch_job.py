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
    """parametrize and launch one or several jobs
    via subprocess in SLURM"""

    def __init__(self, project_folder, name=None, num_nodes=None, QOS=None):

        self.project_folder = project_folder
        self.name_str = ""
        self.num_nodes_str = ""
        self.QOS_str = ""

        self._get_name_str(name)
        self._get_num_node_str(num_nodes)
        self._get_qos_str(QOS)

        self.vasp_exec = self._get_job_type()

        self.workdir = None
        self.array = None

    def _get_job_type(self):
        try:
            if input("do_bader in the array ? [Y]")[0] == "Y":
                return VaspExecutable.bader_exec()
        except Exception:
            pass
        return VaspExecutable()

    def _get_name_str(self, input_name):
        if input_name is not None:
            name = input_name
        else:
            try:
                name = input("Name of the job ? : \n")
                assert len(name) > 0, "invalid name"
                name = name.replace(" ", "_")
            except Exception:
                name = "vasp_job"

        print("Base name : {} \n".format(name))
        self.name_str = "-J {}".format(name)

    def _get_qos_str(self, input_qos_str):
        try:
            if input_qos_str is not None:
                qos_id = input_qos_str
            elif input('QOS ? [BG] / [N]ormal : ')[0] in ["B", "b"]:
                qos_id = "BG"
            else:
                qos_id = ""
        except Exception:
            qos_id = ""

        self.QOS_str = "--qos=bg --requeue" if qos_id == "BG" else ""

    def _get_num_node_str(self, num_node):
        try:
            if num_node is not None:
                num = num_node
            else:
                num = int(input('num nodes / job ? : '))
            assert num > 0 and num <= 8, "invalid node nb"
        except Exception:
            print("default to 1")
            num = 1
        self.num_nodes_str = "-N {}".format(num)

    def launch_job(self):
        if self.vasp_exec.array:
            self._launch_array_run()
        else:
            self._launch_individual_runs()
        os.chdir(self.project_folder)

    def _launch_array_run(self):
        "Create a list of all the subfolder of the run"

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

        workdir = "--workdir {}".format(self.project_folder)
        array = get_array_str(self.project_folder)
        self.issue_shell_command(workdir, array)

    def _launch_individual_runs(self):
        "Walk trhough the folders to find valid vasp inputs"

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

            print("Individually chose to run the following folders ? \n",
                  "select [OK / NO ] \n")
        subdir_list = [os.path.join(self.project_folder, o)
                       for o in os.listdir(self.project_folder)
                       if os.path.isdir(os.path.join(self.project_folder, o))]

        for job_folder in subdir_list:
            if ask_for_launch(job_folder):
                workdir = "--workdir {}".format(job_folder)
                self.issue_shell_command(workdir)
            else:
                notify_no_launch(job_folder)

    def issue_shell_command(self, workdir, array=""):
        shell_command = \
            'sbatch {0.name_str} {0.num_nodes_str} {0.QOS_str} {1} {2}  {3}'.format(
                self, array, workdir, self.vasp_exec.explicit_path())

        subprocess.call([shell_command], shell=True)

        # print(shell_command)


class VaspExecutable(object):
    """ define the parameters of a vaspjob and build de filename which corresponds"""

    def __init__(self, set_dir=None, base_job_name="vasp_job", **kw_args):
        """defines the properties of the vasp job script
        kw_args {nb_runs=int, array=bool, gamma=bool, nc=bool} """

        self.set_base_name(set_dir, base_job_name)

        # to avoid inconsistent exec name
        # the properties must be defined in this order
        self.nb_runs = None
        self.array = None
        self.gamma = None
        self.nc = None

        for kw in ["nb_runs", "array", "gamma", "nc"]:
            if kw_args.get(kw, None) is not None:
                setattr(self, kw, kw_args[kw])

        self.define_properties()

    def set_base_name(self, set_dir, base_job_name):
        if set_dir is not None:
            self.setting_dir = set_dir
        else:
            self.setting_dir = os.path.join(
                platform_id.setting_dir(), "job_scripts")
        self.base_job_name = base_job_name
        self.base_name = os.path.join(self.setting_dir, self.base_job_name)

    def define_properties(self,):
        """define run properties in a consistent order"""

        self._get_nb_run()
        self._get_array_bool()

        # we only use gamma computation for
        # - large nb of run (array)
        # - without change in supercell (single run)
        if self.array and self.nb_runs == "single":
            self._get_gamma_bool()

        # we only use non-colinear computation for
        # - unique run (no array)
        # - no change in supercell (single run)
        if not self.array and self.nb_runs == "single":
            self._get_nc_bool()

    def _get_nb_run(self):
        if self.nb_runs is not None:
            return

        self.nb_runs = "single"
        try:
            if input("job type : [s]ingle_run / [d]ouble_run ? ")[0] in ["d", "D", "2"]:
                self.nb_runs = "double"
        except Exception:
            pass

    def _get_array_bool(self):
        if self.array is not None:
            return

        self.array = False
        print("queue strategy ? ")
        try:
            if input('[O]ne by one / [A]rray : ')[0] in ["A", "a"]:
                self.array = True
        except Exception:
            pass

    def _get_gamma_bool(self):
        if self.gamma is not None:
            return

        self.gamma = False
        try:
            if input('gamma only ? [Y]/n)')[0] == "Y":
                self.gamma = True
        except Exception:
            pass

    def _get_nc_bool(self):
        if self.nc is not None:
            return
        self.nc = False
        try:
            if input('non-colinear ? (Y]/n)')[0] == "Y":
                self.nc = True
        except Exception:
            pass

    def explicit_path(self):
        """return explicit jobfile name with path and specified extension
        requires the job names to be properly formatted"""
        str_list = [self.base_name, self.nb_runs]

        for prop in ['array', "gamma", "nc"]:
            if getattr(self, prop):
                str_list.append(prop)
        print("_".join(str_list))

        return "_".join(str_list)

    @classmethod
    def bader_exec(cls):
        "return the name of the bader executable"
        return VaspExecutable(set_dir=None, base_job_name="bader_job",
                              nb_runs="single", array=True, gamma=False, nc=False)


if __name__ == '__main__':

    if len(sys.argv) > 1:
        PROJECT_FOLDER = sys.argv[1]
    else:
        PROJECT_FOLDER = os.getcwd()

    main(PROJECT_FOLDER)
