#!/usr/bin/env python3

import importlib
import os
import sys
import traceback

import matplotlib
import matplotlib.pyplot as plt

import electronic_analysis.energy_surface as PES
import electronic_analysis.read_bader as bader
import electronic_analysis.read_cohp as lob
import electronic_analysis.read_dos as DOS
import electronic_analysis.read_hull as hull
import electronic_analysis.read_mag_props as nupdown
import electronic_analysis.readO2 as O2
import electronic_analysis.readRun_entries as read

import structure_analysis.bailar_twist as bailar

import filtering_runs.filter_runs as filter_runs
import utils.generic_plot as generic_plot
import utils.platform_id as platform_id

print(matplotlib.get_backend())
print(sys.version)


SETTING_DIR = platform_id.setting_dir()


def filter_loop(run_list_all, input_graph_type=None, allow_filtering=True):
    """
    outer loop to select (and change) the dataset
    before entering the plot loop (with a fixed dataset)
    """

    continue_filter_loop = True

    while continue_filter_loop:
        x_coord = filter_runs.select_x_coord()
        if allow_filtering:
            restricted_runs = filter_runs.restrict_run_list(run_list_all,
                                                            x_coord)
        else:
            restricted_runs = run_list_all

        print("nb runs : {}".format(len(restricted_runs)))
        if len(restricted_runs) == 0:
            print("no enough runs ! ")
        else:
            plot_loop(restricted_runs,
                      input_graph_type=input_graph_type,
                      x_coord=x_coord)

        if allow_filtering is False or input(
                "Continue plotting with different filter ? [Y/n]") != "Y":
            continue_filter_loop = False

    return True


def plot_loop(restricted_runs, input_graph_type=None, x_coord="x_na"):
    " loop to draw several plot on a fixed dataset "

    looping_mode = True
    if input_graph_type is not None:
        # exit after 1 loop if the function is called with an input_graph_type
        looping_mode = False
        graph_type = input_graph_type

    print("============= PRE-TREATMENT OF DATA ====================")
    bader_done = False
    if input("generate bader tags ? : Y/n  : ") == "Y":
        bader_done = bader.get_bader_tags(restricted_runs)

    chem_env_done = False  # to avoid deprecated bailar twist functions switches
    # if not chem_env_done and input(
    #         "generate chem env tags ? : Y/n  : ") == "Y":
    #     restricted_runs = bailar.get_chem_env_tags(restricted_runs)
    #     chem_env_done = True

    print("============= PLOTTING ====================")
    while True:
        if looping_mode:
            graph_type = ask_graph_type()

        continue_loop = single_analysis_routine(restricted_runs,
                                                graph_type, x_coord,
                                                chem_env_done, bader_done,)

        # looping ONLY if initial input is undefined AND user want another graph
        if not (looping_mode and continue_loop):
            break


def single_analysis_routine(rundict_list, graph_type, x_coord, chem_env_done, bader_done):
    """
    chose & perform a single analysis routine
    catch exceptions and reload modules if necessary
    """
    reload_all = False
    continue_loop = True
    try:
        if graph_type == 'Structure':
            bailar.plot_structure_graphs(
                rundict_list, chem_env_done)

        elif graph_type == 'XRD':
            print("Not implemented yet")

        elif graph_type == 'Bader':
            if bader_done:
                bader.plot_charge_and_mag(rundict_list, coord=x_coord)
            else:
                print("Generate bader tags before plotting !!")

        elif graph_type == 'COOP':
            sorting = "oxidation" if bader_done else "OO_pairs"
            rundict_list = lob.plot_COOP_OO(
                rundict_list, sorting=sorting)

        elif graph_type == "O2 release":
            rundict_list = O2.O2_computation(
                rundict_list, bader_done=bader_done)

        elif graph_type == "DOS":
            DOS.plot_DOS_graphs(rundict_list)

        elif graph_type == 'hull':
            if len([d for d in rundict_list if d.status >= 3]) >= 2:
                hull.plot_hull_graphs(rundict_list, coord=x_coord)
            else:
                print("not enough converged runs to perform further analysis")

        elif graph_type == "energy surface":
            PES.plot_energy_surface_graphs(
                rundict_list, chem_env_done)

        elif graph_type == "mag":
            nupdown.plot_all_graphs(rundict_list)

        elif graph_type == "debug_infos":
            print("=== debug output ===")
            print_debug(rundict_list)

        elif graph_type == "reload":
            print("normal reload")
            reload_all = True

        elif graph_type is "QUIT":
            print("finished plotting \n")
            input("press any key to close all figures")
            plt.close("all")
            continue_loop = False

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

    return continue_loop


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
    "debug_infos",
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


def print_debug(rundict_list):
    "print all information of each run"
    for rundict in rundict_list:
        print("\n\n")
        print(rundict.as_dict())


def main():
    "main function"

    if len(sys.argv) > 1:
        working_dir = sys.argv[1]
    else:
        working_dir = os.getcwd()

    read.initialize(working_dir)

    rundict_list = read.collect_valid_runs(working_dir)

    print("number of valid runs " + str(len(rundict_list)))
    if len(rundict_list) == 0:
        print("no run to show, yow ! ")
        exit(1)

    filter_loop(rundict_list)


if __name__ == '__main__':
    main()
