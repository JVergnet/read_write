# generic_plot.py
import copy
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (AutoMinorLocator, FormatStrFormatter,
                               FuncFormatter, MaxNLocator, MultipleLocator,
                               StrMethodFormatter)

import platform_id
import read_hull as hull
import readRun_entries as read


def set_mpl_rc_params():
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
    mpl.rcParams['axes.labelsize'] = 25
    mpl.rcParams['xtick.labelsize'] = 20
    mpl.rcParams['ytick.labelsize'] = 20
    mpl.rc('legend', fontsize=20)


def set_mpl_rc_params_paper():
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
    mpl.rcParams['axes.labelsize'] = 10
    mpl.rcParams['xtick.labelsize'] = 8
    mpl.rcParams['ytick.labelsize'] = 8
    mpl.rc('legend', fontsize=8)


def minor_locator():
    return(AutoMinorLocator(n=2))


def major_locator():
    # , steps = [1,5]))
    return(MaxNLocator(nbins=5, steps=[1, 2, 5], min_n_ticks=3))
# def 2_digit_formatter():
#     return(FormatStrFormatter("%

# MultipleLocator(0.05)


def plot_site_value_evolution(
        sorted_entries,
        specie,
        value='charge',
        coord="x_na",
        plot_type=None,
        axe0=None):
    # get the value (charge or  magmom) of each site of a given specie
    # for each structure along the deintercalation path
    # plot them on a graph
    # each site individually and in average (accounting for site multiplicity)
    # plot_type = array of in representing plotting options in PLOT_TYPE_LIST
    # eg : plot_type = [0,1] for avg and persite plotting
    set_mpl_rc_params()

    PLOT_TYPE_LIST = ["per_site", "average", "min", "max", "sum"]

    if plot_type is None:
        plot_type = set()
        print(" \nPlotting type for {} of {} : \n {}\n "
              .format(value, specie, PLOT_TYPE_LIST))
        while True:
            try:
                plot_type.add(int(input(
                    "type plot type number or [Q]uit : ")))
                print("plot type choice  {}".format(
                    [PLOT_TYPE_LIST[n] for n in plot_type]))
            except Exception as ex:
                print("selection finished")
                continue

    plot_type_str = [PLOT_TYPE_LIST[n] for n in plot_type]

    if axe0 is None:
        fig = plt.figure("{} of {}".format(value, specie))
        # fig.suptitle("{} evolution of {}".format(value,specie),
        #              fontsize="large")
        axe = fig.add_subplot(1, 1, 1)
    else:
        axe = axe0

    stacking_list = list(set([e.stacking for e in sorted_entries]))
    for stack in stacking_list:
        rundict_list = copy.deepcopy(
            [s for s in sorted_entries if s.stacking == stack])

        X_all_sites = []
        Y_all_sites = []

        X = []
        Y_mean_list = []
        Y_min_list = []
        Y_max_list = []
        Y_sum_list = []
        if (len(rundict_list) > 1 and len(
                set([getattr(run, coord) for run in rundict_list])) == 1):
            coord = None

        for i, run in enumerate(rundict_list):
            try:
                nbSpecie = len(run.structure.indices_from_symbol(specie))
                # runDict["x_na"]= round(D["Na"]/nbCell , 2)
                # nbCell = run.nb_cell
                # print("nb cell : {}    nb specie  : {} " .format(nbCell , nbSpecie ) )
                if nbSpecie > 0:
                    x = i if coord is None else getattr(run, coord)
                    Ysum = 0
                    Y_min = 1e10
                    Y_max = -1e10
                    # get the value for each non equiv pt
                    for site in [
                            s for s in run.structure_data.sites
                            if s.specie.name == specie]:
                        Y_all_sites += [site.properties[value]]
                        X_all_sites += [x]
                        Ysum += site.properties[value]
                        Y_min = min(Y_min, site.properties[value])
                        Y_max = max(Y_max, site.properties[value])
                        # print(Ym)
                    Y_sum_list.append(Ysum / run.nb_cell)
                    Y_mean_list.append(Ysum / nbSpecie)
                    Y_min_list.append(Y_min)
                    Y_max_list.append(Y_max)
                    X.append(x)

                else:
                    print(" no {} in the structure !!".format(specie))
            except KeyError as ex:
                print(
                    "Missing property {} in structure {}".format(
                        ex, run.nameTag))

        if "per_site" in plot_type_str:
            axe.scatter(X_all_sites, Y_all_sites,  color="black",
                        s=12, label="per site {}".format(stack))
        if "average" in plot_type_str:
            axe.plot(X, Y_mean_list, linewidth=1.0, linestyle="--",
                     label="average {}".format(stack), color="red")
        if "min" in plot_type_str:
            axe.plot(X, Y_min_list, linewidth=1.0, linestyle="--",
                     label="minimum {}".format(stack))  # color="green",
        if "max" in plot_type_str:
            axe.plot(X, Y_max_list, linewidth=1.0, linestyle="--",
                     label="maximum {}".format(stack))  # color="green",

        if "sum" in plot_type_str:
            Y = np.array(Y_sum_list)
            Y = Y - min(Y)
            axe.plot(X, Y, "ko--", label="sum")
            axe.set_ylim([-0.05, 1])

    axe.legend()   # fontsize=18, fontname ="Arial")
    value_str = value
    if value == "charge":
        value_str = "Bader net population"
        # axe.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        axe.yaxis.set_major_locator(major_locator())
        axe.yaxis.set_minor_locator(minor_locator())
    axe.set_ylabel(value_str)  # , fontsize=18, fontname ="Arial" )

    if coord is not None:
        axe.set_xlabel(coord)  # ,fontsize=18, fontname ="Arial" )
        axe.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    else:
        # axe.xaxis.set_major_locator(MultipleLocator(1))
        labels = [r.nameTag for r in rundict_list]
        positions = [i for i in range(len(rundict_list))]
        # axe.set_xticklabels(labels)
        plt.xticks(positions, labels, rotation='horizontal')

    if axe0 is None:
        plt.show(block=False)
        return(fig)
    else:
        return(axe0)


def plot_structure_value_evolution(
        sorted_entries,
        prop_list=['dOO_min', 'bandgap'],
        legend=None,
        x_na_coords=None):
    stacking_list = list(set([e.stacking for e in sorted_entries]))
    # for e in sorted_entries :
    #     if e["stacking"] not in stacking_list :
    #         stacking_list.append(e["stacking"])
    set_mpl_rc_params()
    if legend is None:
        legend = ""
    else:
        legend = " of {}".format(legend)

    if x_na_coords is None:
        x_na_coords = "x_na" if len(
            set([run.x_na for run in sorted_entries])) > 1 else "name"
        print("automatic determination of coords : {}".format(x_na_coords))

    # prop_list = ['dOO_min' , 'bandgap']
    converged_entries = [d for d in sorted_entries if d.status >= 3]
    clean_entries = [d for d in sorted_entries if d.status == 4]
    hull_entries = [d for d in sorted_entries if d.status >= 5]

    for i, prop in enumerate(prop_list):
        plotTitle = prop + legend

        fig = plt.figure(plotTitle)
        fig.suptitle(plotTitle, fontsize="large")
        axe = fig.add_subplot(1, 1, 1)

        if x_na_coords == "name":
            # X_Y_name = np.array([ [ i,s.get(prop, 0),s["nameTag"] ] for (i,s) in enumerate(converged_entries)
            #                 if (s.get(prop, None) is not None ) ] )
            #print(X , Y)
            #axe.yaxis.set_major_formatter(FuncFormatter(lambda x, loc: "{:.2f}".format(x) ) )
            #axe.xaxis.set_major_formatter(FuncFormatter(lambda x, loc: "{:.2f}".format(x) ) )
            # axe.plot(X_Y_name[:,0],X_Y_name[:,1],)
            axe.plot([getattr(s, prop) for s in converged_entries
                      if hasattr(s, prop)], "o")
            labels = [s.nameTag
                      for s in converged_entries
                      if hasattr(s, prop)]
            plt.xticks([r for r in range(len(labels))],
                       labels, rotation='horizontal')
            #axe.set_yticklabels(['{:.2f}'.format(x) for x in axe.get_yticks()])
            #y_fmt = FormatStrFormatter('%2.2f')
            axe.yaxis.set_major_formatter(StrMethodFormatter("{x:.3f}"))

            # print(X_Y_name[:,1])

        else:
            for stack in stacking_list:
                stack_list = copy.deepcopy(
                    [s for s in sorted_entries if s.stacking == stack])
                if x_na_coords == "x_na":
                    stack_list = hull.generate_hull_entries(stack_list)
                    X_Y = np.array([[s.x_na, getattr(s, prop)] for s in stack_list
                                    if (s.status >= 4 and hasattr(s, prop))])
                    # print(X , Y)
                    axe.plot(X_Y[:, 0], X_Y[:, 1], "o--",
                             label="minimal {}".format(stack))

                X_Y = np.array([[getattr(s, x_na_coords), getattr(s, prop)]
                                for s in stack_list
                                if (s.status == 3 and hasattr(s, prop))])
                # print(X , Y)
                try:
                    axe.plot(X_Y[:, 0], X_Y[:, 1], "x",
                             label="non-min {}".format(stack))
                except Exception as ex:
                    print(ex)
            if x_na_coords == "x_na":
                for entries, color, lab in [
                        (hull_entries, "black", "on hull"),
                        (clean_entries, "red", "off-hull minimum")]:
                    try:
                        X_Y = np.array([[s.x_na, getattr(s, prop)]
                                        for s in entries
                                        if hasattr(s, prop)])
                        X = X_Y[:, 0]
                        Y = X_Y[:, 1]
                        axe.scatter(X, Y,
                                    s=180,
                                    facecolors='white',
                                    edgecolor=color,
                                    label=lab)
                    except Exception as ex:
                        print(
                            "exception while plotting {} circles : {}".format(
                                lab, ex))

            axe.set_xlabel(x_na_coords)  # ,fontsize=18, fontname ="Arial" )
            axe.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            axe.legend()
        print(prop.replace("_", " "))
        axe.set_ylabel(prop.replace("_", " "))

#        plt.show(block=False)
#        save_fig(fig,plotTitle)
    return(fig)


def save_fig(fig, plotTitle, save_mode=None, folder=None):
    figName = None
    if save_mode is None:
        save_mode = input(
            "Type [s]ave to save {0} ([d]efault name) \n".format(plotTitle))

    if folder is None:
        folder = os.path.join(platform_id.local_cluster_dir(), "figures")

    os.chdir(folder)
    print("in", os.getcwd())

    try:
        if save_mode[0] == "s":
            figName = input("Type the file name (SVG format) : ")

        elif save_mode[0] == "d":

            figName = read.get_file_name(
                folder, plotTitle.replace(
                    ' ', '_'), ext=".svg")
            print("default name : {}".format(figName))
            #print(os.getcwd(), param['mainFolder'])
    except IndexError:
        print("option not recognized, figure not saved")

    if figName is not None:
        fig.savefig("{}.svg".format(figName), bbox_inches='tight')
        print("{} saved in {}".format(figName.split("/")[-1], folder))

    return(figName)
