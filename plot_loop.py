#!/usr/bin/env python3
import importlib
import os
import sys
import traceback

import matplotlib
import matplotlib.pyplot as plt

import bailar_twist as bailar
import DOS_plot as DOS
import energy_surface as PES
import generic_plot as generic_plot
import lobster_coop as lob
import nupdown_scan as nupdown
import platform_id
import read_hull as hull
import readBader as bader
import readO2 as O2
import readRun_entries as read

print(matplotlib.get_backend())
print(sys.version)


SETTING_DIR = platform_id.setting_dir()


def generate_plot(run_list_all, input_graph_type=None, allow_filtering=True):
    """
    allow to plot several quantities if graph_type option is set to None
    plot "none" to exit the loop
    """

    filter_loop = True

    # read.generate_tags(run_list_all, minimal=True)
    all_runs_input = run_list_all

    while filter_loop:

        # all_runs = all_runs_input

        restricted_runs = read.restrict_run_list(all_runs_input) \
            if allow_filtering else all_runs_input

        print("nb runs : {}".format(len(restricted_runs)))
        if len(restricted_runs) == 0:
            print("no enough runs ! ")
        else:
            print("============= PRE-TREATMENT OF DATA ====================")

            if input("generate bader tags ? : Y/n  : ") == "Y":
                bader_done = bader.get_bader_tags(restricted_runs)

            chem_env_done = False  # to avoid deprecated bailar twist functions switches
            # if not chem_env_done and input(
            #         "generate chem env tags ? : Y/n  : ") == "Y":
            #     restricted_runs = bailar.get_chem_env_tags(restricted_runs)
            #     chem_env_done = True

            print("============= PLOTTING ====================")
            general_plot_loop(
                restricted_runs,
                input_graph_type,
                chem_env_done,
                bader_done)

        if allow_filtering is False or input(
                "Continue plotting with different filter ? [Y/n]") != "Y":
            filter_loop = False

    return(True)


def general_plot_loop(restricted_runs, input_graph_type=None, chem_env_done=False, bader_done=False):
    " switch between quit // analysis "

    plot_loop = True
    while plot_loop:

        if input_graph_type is None:
            # if no input_graph_type is defined, ask user for graph_type
            graph_type = ask_graph_type()
        else:
            # exit after 1 loop if the function is called with an input_graph_type
            graph_type = input_graph_type
            plot_loop = False

        if graph_type is "QUIT":
            print("finished plotting \n")
            input("press any key to close all figures")
            plt.close("all")
            break

        chose_analysis_plot(
            graph_type,
            restricted_runs,
            chem_env_done,
            bader_done,)


def chose_analysis_plot(graph_type, restricted_runs, chem_env_done, bader_done):
    """
    switch between all analysis routines
    catch exceptions and reload analysis modules if necessary
    """
    reload_all = False
    try:
        if graph_type == 'Structure':
            bailar.plot_structure_graphs(
                restricted_runs, chem_env_done)

        elif graph_type == 'XRD':
            print("Not implemented yet")

        elif graph_type == 'Bader':
            if bader_done:
                bader.plot_charge_and_mag(restricted_runs)
            else:
                print("Generate bader tags before plotting !!")

        elif graph_type == 'COOP':
            sorting = "oxidation" if bader_done else "OO_pairs"
            restricted_runs = lob.plot_COOP_OO(
                restricted_runs, sorting=sorting)

        elif graph_type == "O2 release":
            restricted_runs = O2.O2_computation(
                restricted_runs, bader_done=bader_done)

        elif graph_type == "DOS":
            DOS.plot_DOS_graphs(restricted_runs)

        elif graph_type == 'hull':
            if len([d for d in restricted_runs if d.status >= 3]) >= 2:
                hull.plot_hull_graphs(restricted_runs)
            else:
                print("not enough converged runs to perform further analysis")

        elif graph_type == "energy surface":
            PES.plot_energy_surface_graphs(
                restricted_runs, chem_env_done)

        elif graph_type == "mag":
            nupdown.plot_mag_graphs(restricted_runs)

        elif graph_type == "reload":
            print("normal reload")
            reload_all = True

    except Exception as ex:
        print(traceback.format_exc())
        print(ex)
        reload_all = True

    if reload_all:
        plt.close("all")
        os.chdir(SETTING_DIR)
        try:
            importlib.reload(lob)
            importlib.reload(PES)
            importlib.reload(DOS)
            importlib.reload(hull)
            importlib.reload(bailar)
            importlib.reload(bader)
            importlib.reload(O2)
            importlib.reload(nupdown)
            importlib.reload(generic_plot)
            importlib.reload(read)
        except Exception:
            print(traceback.format_exc())
            print("PLEASE correct the issue and reload once more")
    else:
        plt.show(block=False)

    return reload_all


GRAPH_OPTION = [
    "Structure",
    "Bader",
    "XRD",
    "DOS",
    "hull",
    "COOP",
    "O2 release",
    "energy surface",
    "mag",
    "reload",
    "QUIT"]


def ask_graph_type(graph_type=None):
    "ask graph type is not given"
    graph_option_to_read = "Graph type option : \n" + \
        "\n".join(["|{} : {}|".format(i, option)
                   for (i, option) in enumerate(GRAPH_OPTION)])
    while graph_type is None:
        try:
            print(graph_option_to_read)
            graph_type = GRAPH_OPTION[
                int(input("\nEnter graph type number : \n"))]
            print("option choice : {}".format(graph_type))
        except Exception as ex:
            print("Bad Option {}".format(ex))
            if input(
                    "type C to [C]ontinue plotting or anything else to QUIT") != "C":
                graph_type = "QUIT"
    return graph_type


if __name__ == '__main__':

    read.initialize()

    vaspRunDictList = read.get_vasp_run_dict_list()

    nbRun = len(vaspRunDictList)
    print("number of valid runs " + str(nbRun))
    if nbRun == 0:
        print("no run to show, yow ! ")
        exit(1)

    generate_plot(vaspRunDictList)
