# read_hull.py

# import logging

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.ticker import AutoMinorLocator

# from operator import itemgetter
# import readRun_entries as read
import generic_plot as generic_plot


def generate_hull_entries(run_list, remove_extremes=False, coord="xNa"):

    # sort the vasRunDictList entries (a list of dict) in 3 categories
    # add formation energy ("eform" key) to each vasprun_dict

    # input ; [{vaspRun : v , param : value }
    # output : [sorted_entries , clean_entries , hull_entries]

    # clean_entries : selected  entries of lowest energy for each xNa
    # hull_entries : clean_entries that lie on the convex hull

    # Sort the entries according to 1) their Na proportion and 2) their energy

    converged_entries = [d for d in run_list if d.status >= 3]
    sorted_entries = sorted(converged_entries,
                            key=lambda x: (
                                getattr(x, coord),
                                x.energy_per_FU))

    if remove_extremes is None and input(
            "REMOVE  extreme composition x=0 x=1 ?") == "Y":
        remove_extremes = True

    if remove_extremes:
        sorted_entries = [x for x in sorted_entries
                          if x.xNa < 1 and x.xNa > 0]

    # Clean Entries
    # =================================================
    # Select the structure of lowest energy for each Na
    # xNa = []
    # voltage = []

    clean_entries = []
    current_Na = -1
    for entry in sorted_entries:
        if getattr(entry, coord) > current_Na:
            clean_entries.append(entry)
            current_Na = getattr(entry, coord)
    #       print(entry)
    # print(clean_entries)
    # [print([r.xNa for r in L])
    #  for L in [sorted_entries, clean_entries]]
    # compute the formation energy in meV
    # with respect to endmembers (x->0 & x->1)
    [[xNa0, Ex0], [xNa1, Ex1]] = [[getattr(clean_entries[i], attr)
                                   for attr in [coord, "energy_per_FU"]]
                                  for i in [0, -1]]
    print([[xNa0, Ex0], [xNa1, Ex1]])
    for entry in sorted_entries:
        x = getattr(entry, coord)
        entry.eform = (entry.energy_per_FU
                       - (x - xNa0) / (xNa1 - xNa0) * Ex1
                       - (xNa1 - x) / (xNa1 - xNa0) * Ex0) * 1000
    # print("sorted entries {0} \n\n clean entries {1}\n\n"
    #      .format(sorted_entries,clean_entries))

    # Hull Entries
    # =============================================
    # Select structures that lie on the convex hull

    """
from scipy.spatial import ConvexHull
hull = ConvexHull(points)
hull_pts = points[hull.vertices,0], points[hull.vertices,1]
"""

    current_index = 0
    hull_entries = [clean_entries[current_index]]

    while current_index < len(clean_entries) - 1:
        min_slope = +1e10
        next_index = current_index + 1
        # print("current index : {}".format(current_index))

        for j in range(current_index + 1, len(clean_entries), 1):
            [dNa, delta_E] = [getattr(clean_entries[j], attr) -
                              getattr(clean_entries[current_index], attr)
                              for attr in [coord, 'eform']]
            slope = delta_E/dNa
            # print(j)
            if slope < min_slope:
                min_slope = slope
                next_index = j
                # print("new index : {} , slope : {}".format(j, min_slope))

        # seg_list.append([current_index,next_index])
        current_index = next_index
        hull_entries.append(clean_entries[current_index])

    for entry in clean_entries:
        entry.status = 4

    for entry in hull_entries:
        entry.status = 5

    return(sorted_entries)


# HULL RELATED GRAPHS
# =============================================================================


def plot_hull_graphs(sorted_entries, **kwargs):
    "kwargs : 'coord', 'simple_energy_graph' ,'hull','voltage' = True/False"
    # print(hull_entries)
    converged_entries = [d for d in sorted_entries if d.status >= 3]
    if "coord" in kwargs.keys():
        coord = kwargs["coord"]
    elif len(set([d.xNa for d in sorted_entries])) > 1:
        coord = "xNa"
    elif input("convex hull on the doo? [Y]") in ["Y", "y"]:
        coord = "dOO_min"
    else:
        coord = None

    if "simple_energy_graph" in kwargs.keys():
        simple_energy_graph = kwargs["simple_energy_graph"]
    elif input("plot simple energy graph ? [Y/N] ") == "Y":
        simple_energy_graph = True
    else:
        simple_energy_graph = False
    if simple_energy_graph:
        sorted_converged_entries = sorted(
            converged_entries, key=lambda x: x.energy_per_FU)
        generic_plot.plot_structure_value_evolution(
            sorted_converged_entries,
            prop_list=['energy_per_FU'],
            legend=None,
            xNa_coords=None)
        for s in sorted_converged_entries:
            print("{} {}/n".format(s.id, s.energy_per_FU))

    if coord is not None:
        print("\n==== HULL GRAPHS ======\n")
        if "hull" in kwargs.keys():
            hull = kwargs["hull"]
        else:
            hull = input("plot convex_hull ? [Y/N] ") == "Y"

        if hull:
            sorted_entries = generate_hull_entries(sorted_entries)
            converged_entries = [d for d in sorted_entries if d.status >= 3]
            hull_entries = [d for d in sorted_entries if d.status >= 5]
            plot_convex_hull(converged_entries, coord=coord)

            if coord == "xNa":
                if "voltage" in kwargs.keys():
                    voltage = kwargs["hull"]
                else:
                    voltage = (input("plot voltage curve ? [Y/N] ") == "Y")

                if voltage:
                    plot_voltage_curve(hull_entries)


def plot_convex_hull(sorted_entries, coord='xNa'):

    clean_entries = [d for d in sorted_entries if d.status >= 4]
    hull_entries = [d for d in sorted_entries if d.status >= 5]

    stacking_list = set([e.stacking for e in sorted_entries])

    print("stacking list : ", stacking_list)

    # [(1,0,0,1),(0,0,1,1),(0,1,0,1)]
    colors = ['green', 'red', 'blue', 'black']

    color_dict = {}
    for i, stacking in enumerate(stacking_list):
        color_dict[stacking] = colors[i]
    print(color_dict)
    specific_color_dict = {"P3": "#b700ffff",
                           "O3": "#00baffff"}
    color_dict.update(specific_color_dict)
    plotTitle = "Convex Hull"

    fig = plt.figure(plotTitle)

    # fig.suptitle(hull_entries[-1]["folder"]+"\n"
    #              +hull_entries[-1]["formula"]+"\n"
    #              +plotTitle,
    #              fontsize="large")
    # plotting
    dotSize = 40

    axe = fig.add_subplot(1, 1, 1)

    # Plot error bar on energies of structures with lowest energy at given xNa
    # symmetric or not depending if they are (or not) on the convex hull

    if True:  # input("add error bars for ambient temp ? [Y/N] ")=="Y" :
        thermal_error = 25  # in meV
        XE1 = np.array([[getattr(e, coord), e.eform] for e in clean_entries])
        # print(XE1)
        X_clean = XE1[:, 0]
        E_clean = XE1[:, 1]
        yerr_down = [(thermal_error if (e.status < 4) else 0)
                     for e in clean_entries]
        yerr_up = thermal_error * np.ones_like(E_clean)
        # print([yerr_down, yerr_up])
        axe.errorbar(
            X_clean, E_clean,
            yerr=[yerr_down, yerr_up],
            fmt='',
            linestyle="None")

    # # Plot all entries (hollow circles) [sorted_entries]
    for entry in sorted_entries:
        try:
            print(entry.eform)
        except Exception as ex:
            print("Exception raised for {} : {}".format(entry.nameTag, ex))
    for stacking in stacking_list:
        print("current stacking", stacking)
        X_E = np.array([[getattr(entry, coord), entry.eform]
                        for entry in sorted_entries
                        if entry.stacking == stacking])
        X = X_E[:, 0]
        E = X_E[:, 1]
        # print(E)
        axe.scatter(
            X, E, s=dotSize,
            facecolors='none',
            edgecolor=color_dict[stacking],
            label=stacking)

    # Plot entries with lowest energies (filled circles) [clean_entries]

    X_E_C = [[getattr(entry, coord), entry.eform, color_dict[entry.stacking]]
             for entry in clean_entries]
    X = [x[0] for x in X_E_C]
    E = [x[1] for x in X_E_C]
    C = [x[2] for x in X_E_C]
    # print(E)
    axe.scatter(X, E, s=dotSize, color=C, edgecolor='face')

    # link entries on the convex hull (black line) [hull_entries]

    X_E = np.array([[getattr(entry, coord), entry.eform]
                    for entry in hull_entries])
    X = X_E[:, 0]
    E = X_E[:, 1]
    axe.plot(X, E, "k-")

    axe.legend()
    axe.set_ylabel('E total (meV)')
    axe.set_xlabel('Na content'
                   )
    # X = [getattr(entry, coord) for entry in sorted_entries]
    axe.invert_xaxis()
    # axe.tick_params(axis='x', which='minor', bottom=True)
    axe.xaxis.set_minor_locator(AutoMinorLocator(n=2))
    # xlim(max(X), min(X))

    plt.show(block=False)
    generic_plot.save_fig(fig, plotTitle)

    return(fig)


def plot_voltage_curve(hull_entries):
    # compute OCV between stable intermediate compositions
    # see ref below :
    # https://www.nature.com/articles/npjcompumats20162

    voltage_list = []
    eNa = -2  # -1.47
    for i in range(0, len(hull_entries) - 1):
        dNa = hull_entries[i + 1].xNa - hull_entries[i].xNa
        E1 = hull_entries[i].energy_per_FU
        E2 = hull_entries[i + 1].energy_per_FU
        # print("E1={0} \n E2={1} \n dNa = {2} \nV = -(E2-E1-eNa) = {3}"
        #     .format(E1,E2,dNa,-((E2-E1)/dNa-eNa)))

        voltage_list.append(-((E2 - E1) / dNa - eNa))

    voltage_col = LineCollection([[(hull_entries[i].xNa, voltage_list[i]),
                                   (hull_entries[i + 1].xNa, voltage_list[i])]
                                  for i in range(0, len(hull_entries) - 1)])

    plotTitle = "Predicted Voltage Curve"

    fig = plt.figure(plotTitle)

    # fig.suptitle(hull_entries[-1]['folder']+"\n"
    #              +hull_entries[-1]['formula']+"\n"
    #              +plotTitle,
    #              fontsize="large")

    axe = fig.add_subplot(1, 1, 1)
    axe.add_collection(voltage_col)
    axe.set_ylim([min(voltage_list) - 1, max(voltage_list) + 1])

    # axe.legend()

    axe.set_ylabel('V equilibrium')
    axe.set_xlabel('Na content')

    plt.show(block=False)

    generic_plot.save_fig(fig, plotTitle)

    return(fig)
