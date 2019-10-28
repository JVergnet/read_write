# platform_id.py

import os
import platform
from os.path import expanduser


def local_cluster_dir():

    home = expanduser("~")
    cluster_names = {'frodon.lsd.univ-montp2.fr': "frofro",
                     'bipbip.lsd.univ-montp2.fr': "bipbip",
                     'tornado': "tornado",
                     'vergnet-cse': "frofro"}
    cluster_name = cluster_names[platform.node()]
    return(os.path.join(home, cluster_name))


def setting_dir():
    return(os.path.join(local_cluster_dir(), "settings"))


def running_on_cluster():
    # Check if we are on frodon
    cluster_names = ['frodon.lsd.univ-montp2.fr',
                     'bipbip.lsd.univ-montp2.fr',
                     'tornado']
    if platform.node() in cluster_names:
        cluster = True
    else:
        cluster = False

    return(cluster)


def first_check_cluster():
    # Check if we are on frodon when the program is first called

    print("platform node : {}".format(platform.node))

    if running_on_cluster():
        print(
            """
=============================
ReadRun is executed on the cluster
so some analysis features will be turned off
but jobs may be launched from the script !
=============================
""")
        cluster = True
    else:
        print("not running on a cluster : more features but no job launch ! ")
        cluster = False

    return(cluster)
