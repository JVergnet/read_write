# energy_surface.py

import importlib
# import os
# import subprocess
# import sys
from operator import itemgetter

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # , axes3d
# from scipy.interpolate import griddata

import bailar_twist as bailar
import PES_complex_plot as PES_complex_plot
import read_hull as hull
import readRun_entries as read
import generic_plot


def plot_scatter_mesh_cell_shape(struct_list):
    colors = ['green', 'red', 'blue', 'black']
    marks = ['^', 'o', '*', '+']

    xtag = "disto"
    ytag = "x_na"
    XY = np.array([[d[xtag], d[ytag], d["eform"]] for d in struct_list])

    fig = plt.figure()
    ax = Axes3D(fig)

    # plot cell params  as a function of the distortion and x_na

    ZABC = np.array([list(d["structure"].lattice.lengths_and_angles[0])
                     for d in struct_list])

    for i, direction in enumerate(["a", "b", "c"]):
        color = colors[i]
        mark = marks[i]
        ax.scatter(XY[:, 0], XY[:, 1], ZABC[:, i], zdir='z',
                   c=color, marker=mark, label=direction, s=100)
        for x, y, z, energy in zip(
                XY[:, 0], XY[:, 1], ZABC[:, i], XY[:, 2]):  # plot each point + it's index as text above
            ax.text(
                x,
                y,
                z,
                "{:.3f}".format(energy),
                size=10,
                zorder=1,
                color='k')

    ax.set_xlabel(ytag)
    ax.set_ylabel(ytag)
    ax.set_ylabel("cell parameter")
    fig.legend()

    # plot  V of the cell as a function of the distortion and x_na
    fig2 = plt.figure()
    ax2 = Axes3D(fig2)

    ZV = np.array([d["structure"].lattice.volume for d in struct_list])

    ax2.scatter(XY[:, 0], XY[:, 1], ZV, zdir='z', c=color, marker=mark, s=150)
    ax2.set_xlabel(xtag)
    ax2.set_ylabel(ytag)
    ax2.set_zlabel("Volume")
#    fig2.legend()
    plt.show(block=False)
    generic_plot.save_fig(fig, "scatter surface", save_mode=None, folder=None)


def plot_scatter_surface(
        sorted_entries,
        xtag="disto",
        ytag="x_na",
        ztag="eform"):

    fig = plt.figure()
    ax = Axes3D(fig)

    for d in sorted_entries:
        d['disto'] = round(d["distortion"]["Mn"]["O:6"], 2)

    stacking_list = list(set([e["stacking"] for e in sorted_entries]))

    colors = ['green', 'red', 'blue', 'black']
    marks = ['^', 'o', '*', '+']

    for i, stacking in enumerate(stacking_list):
        color = colors[i]
        mark = marks[i]
        struct_list = [d for d in sorted_entries if d["stacking"] == stacking]

        XYZ = np.array([[d[xtag], d[ytag], d[ztag]] for d in struct_list])
        ax.scatter(XYZ[:, 0], XYZ[:, 1], XYZ[:, 2], zdir='z',
                   c=color, marker=mark, label=stacking, s=150)

    ax.set_xlabel(xtag)
    ax.set_ylabel(ytag)
    ax.set_zlabel(ztag)
    fig.legend()

    plt.show(block=False)

    generic_plot.save_fig(fig, "scatter surface", save_mode=None, folder=None)

    return(fig)


def generate_stacking_dict(sorted_entries, chem_env_done):
    stacking_list = list(set([e["stacking"] for e in sorted_entries]))

    stacking_dict = {}

    do_chem_env = False
    if chem_env_done and input("calculate (heavy) chem env ? [Y]/n ") == "Y":
        do_chem_env = True

    for stacking in stacking_list:
        # retrieve the computed single point vasp runs to build the PES
        struct_list = [d for d in sorted_entries if d["stacking"] == stacking]
        disto_mesh_folder = input(
            "Folder to find distortion mesh for stacking : {} ?? \n".format(stacking))
        # if
        # os.path.exists(os.path.join(disto_mesh_folder,"distortion_energy_array.npy")):
        distortion_mesh = None
        # try :
        #     # Read JSON file
        #     with open('distortion_mesh_data.json') as data_file:
        #         distortion_mesh = json.load(data_file)
        #         # distortion_mesh = np.load(os.path.join(disto_mesh_folder,"distortion_energy_array.npy"))
        #         print("distortion mesh data sucessfully loaded")
        # except Exception as ex :
        #     print("Could not load grid data file : {}".format(ex))

        if distortion_mesh is None or input(
                "force distortion mesh reloading ? Y / n ") == "Y":
            distortion_mesh = read.collect_valid_runs(disto_mesh_folder,
                                                      vasprun_parsing_lvl=0.5,
                                                      file_system_choice="p")
            # distortion_mesh = read.generate_tags(distortion_mesh, force=True)
            distortion_mesh = hull.generate_hull_entries(distortion_mesh)
            if do_chem_env:
                distortion_mesh = bailar.get_chem_env_tags(distortion_mesh)

                for d in distortion_mesh:
                    d['disto_prec'] = d["distortion"]["Mn"]["O:6"]

                distortion_mesh = sorted(
                    distortion_mesh, key=itemgetter(
                        'x_na', "disto_prec"))

            # # Write JSON file
            # with io.open('distortion_mesh_data.json', 'w') as outfile:
            #     json.dumps(data,outfile, indent=4, sort_keys=True,
            #                separators=(',', ': '), ensure_ascii=False)

        stacking_dict[stacking] = {
            "run_list": struct_list,
            "mesh": distortion_mesh}
    return(stacking_dict)


def plot_energy_surface_graphs(runList, chem_env_done):

    sorted_entries = [d for d in runList if d["status"] >= 3]
    print("\n==== ENERGY SURFACE GRAPHS ======\n")

    if chem_env_done:

        if input("plot scatter PES ? Y / N  : ") == "Y":
            plot_scatter_surface(sorted_entries)

        if input("plot shape evolution ? Y / N  : ") == "Y":
            plot_scatter_mesh_cell_shape(sorted_entries)

    if input("retrieve distortion mesh data ? Y / N  : ") == "Y":
        all_dicts = generate_stacking_dict(sorted_entries, chem_env_done)

        for stacking, stacking_dict in all_dicts.items():
            struct_list = stacking_dict["run_list"]
            distortion_mesh = stacking_dict["mesh"]
            if chem_env_done and input(
                    "plot shape evolution ? Y / N  : ") == "Y":
                plot_scatter_mesh_cell_shape(distortion_mesh)

            if input("start complex plot loop ? Y / N ") == "Y":
                finish = False
                # we do a loop to redo the plot without loosing all the mesh data that was loaded
                # allo both to debug and redraw a lot faster
                while not finish:
                    importlib.reload(PES_complex_plot)
                    try:
                        PES_complex_plot.plot_3d_angle_energy_from_struct_list(
                            struct_list, distortion_mesh)

                        if chem_env_done:
                            if input("plot PES 3D ? Y / N  : ") == "Y":
                                PES_complex_plot.plot_PES_3D(
                                    struct_list, distortion_mesh)

                            if input("plot PES projected ? Y / N  : ") == "Y":
                                PES_complex_plot.plot_PES_projected(
                                    struct_list, distortion_mesh)

                        if input(
                                "finished complex PES plot ? Y / N  : ") == "Y":
                            finish = True
                    except Exception as ex:
                        print(ex)

    return(True)
