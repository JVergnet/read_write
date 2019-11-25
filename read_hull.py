# read_hull.py

# import logging

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.ticker import AutoMinorLocator

# from operator import itemgetter
# import readRun_entries as read
import generic_plot as generic_plot


def generate_hull_entries(run_list, remove_extremes=False, coord="x_na"):
    """
     sort the vasRunDictList entries (a list of dict) in 3 categories
    add formation energy ("eform" key) to each vasprun_dict

    input ; [{vaspRun : v , param : value }
    output : [sorted_entries , clean_entries , hull_entries]

    clean_entries : selected  entries of lowest energy for each x_na
    hull_entries : clean_entries that lie on the convex hull

    Sort the entries according to 1) their Na proportion and 2) their energy
    """
    print("computing hull in 2D : {} and energy".format(coord))
    converged_entries = [d for d in run_list if d.status >= 3]
    sorted_entries = sorted(converged_entries,
                            key=lambda x: (
                                getattr(x, coord),
                                x.energy_per_fu))

    # Clean Entries
    # =================================================
    # Select the structure of lowest energy for each Na
    # x_na = []
    # voltage = []

    clean_entries = []
    current_Na = -1
    for entry in sorted_entries:
        if getattr(entry, coord) > current_Na:
            clean_entries.append(entry)
            current_Na = getattr(entry, coord)
    #       print(entry)
    # print(clean_entries)
    # [print([r.x_na for r in L])
    #  for L in [sorted_entries, clean_entries]]
    # compute the formation energy in meV
    # with respect to endmembers (x->0 & x->1)
    [[x_na0, Ex0], [x_na1, Ex1]] = [[getattr(clean_entries[i], attr)
                                     for attr in [coord, "energy_per_fu"]]
                                    for i in [0, -1]]
    print([[x_na0, Ex0], [x_na1, Ex1]])
    for entry in sorted_entries:
        x = getattr(entry, coord)
        entry.eform = (entry.energy_per_fu
                       - (x - x_na0) / (x_na1 - x_na0) * Ex1
                       - (x_na1 - x) / (x_na1 - x_na0) * Ex0) * 1000
    # print("sorted entries {0} \n\n clean entries {1}\n\n"
    #      .format(sorted_entries,clean_entries))

    # Hull Entries
    # =============================================
    # Select structures that lie on the convex hull

    # from scipy.spatial import ConvexHull
    # hull = ConvexHull(points)
    # hull_pts = points[hull.vertices,0], points[hull.vertices,1]

    current_index = 0
    hull_entries = [clean_entries[current_index]]

    while current_index < len(clean_entries) - 1:
        min_slope = +1e10
        next_index = current_index + 1
        # print("current index : {}".format(current_index))

        for j in range(current_index + 1, len(clean_entries), 1):
            [d_x, d_e] = [getattr(clean_entries[j], attr) -
                          getattr(clean_entries[current_index], attr)
                          for attr in [coord, 'eform']]
            slope = d_e/d_x
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


def plot_hull_graphs(sorted_entries, coord="x_na", **kwargs):
    "kwargs : 'coord', 'simple_energy_graph' ,'hull','voltage' = True/False"
    # print(hull_entries)
    converged_entries = [d for d in sorted_entries if d.status >= 3]
    if len(set([getattr(d, coord) for d in sorted_entries])) <= 1:
        print("not enough runs with distinct {}".format(coord))
        coord = None

    if "simple_energy_graph" in kwargs.keys():
        simple_energy_graph = kwargs["simple_energy_graph"]
    else:
        simple_energy_graph = True if \
            input("plot simple energy graph ? [Y/N]")[0] == "Y" \
            else False

    if simple_energy_graph:
        sorted_converged_entries = sorted(
            converged_entries, key=lambda x: x.energy_per_fu)
        generic_plot.plot_structure_value_evolution(
            sorted_converged_entries,
            prop_list=['energy_per_fu'],
            legend=None,
            x_na_coords=None)
        for s in sorted_converged_entries:
            print("{} {}/n".format(s.str_id, s.energy_per_fu))

    if coord is not None:
        print("\n==== HULL GRAPHS ======\n")
        if "hull" in kwargs.keys():
            hull = kwargs["hull"]
        else:
            string = "plot convex_hull on {}? [Y/N] ".format(coord)
            hull = bool(input(string) == "Y")

        if hull:
            sorted_entries = generate_hull_entries(sorted_entries)
            converged_entries = [d for d in sorted_entries if d.status >= 3]
            hull_entries = [d for d in sorted_entries if d.status >= 5]
            plot_convex_hull(converged_entries, coord=coord)

            if coord == "x_na":
                if "voltage" in kwargs.keys():
                    voltage = kwargs["hull"]
                else:
                    voltage = (input("plot voltage curve ? [Y/N] ") == "Y")

                if voltage:
                    plot_voltage_curve(hull_entries)


def plot_convex_hull(sorted_entries, coord='x_na'):

    clean_entries = [d for d in sorted_entries if d.status >= 4]
    hull_entries = [d for d in sorted_entries if d.status >= 5]

    stacking_list = {e.stacking for e in sorted_entries}

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
    plot_title = "Convex Hull"

    fig = plt.figure(plot_title)

    # fig.suptitle(hull_entries[-1]["folder"]+"\n"
    #              +hull_entries[-1]["formula"]+"\n"
    #              +plotTitle,
    #              fontsize="large")
    # plotting
    dot_size = 40

    axe = fig.add_subplot(1, 1, 1)

    # Plot error bar on energies of structures with lowest energy at given x_na
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
        x_e = np.array([[getattr(entry, coord), entry.eform]
                        for entry in sorted_entries
                        if entry.stacking == stacking])
        axe.scatter(
            x_e[:, 0], x_e[:, 1], s=dot_size,
            facecolors='none',
            edgecolor=color_dict[stacking],
            label=stacking)

    # Plot entries with lowest energies (filled circles) [clean_entries]

    x_e_c = np.array([[getattr(entry, coord),
                       entry.eform,
                       color_dict[entry.stacking]]
                      for entry in clean_entries])
    # print(E)
    axe.scatter(x_e_c[:, 0], x_e_c[:, 1], s=dot_size,
                color=x_e_c[:, 2], edgecolor='face')

    # link entries on the convex hull (black line) [hull_entries]
    x_e = np.array([[getattr(entry, coord), entry.eform]
                    for entry in hull_entries])
    axe.plot(x_e[:, 0], x_e[:, 1], "k-")

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
    generic_plot.save_fig(fig, plot_title)

    return(fig)


def plot_voltage_curve(hull_entries):
    """
    compute OCV between stable intermediate compositions
    see ref : https://www.nature.com/articles/npjcompumats20162
    """

    e_na = -2  # -1.47

    voltages = []
    voltage_coord = []
    for run_1, run_2 in zip(hull_entries[:-1], hull_entries[1:]):
        d_na = run_2.x_na - run_1.x_na
        d_e = run_2.energy_per_fu - run_1.energy_per_fu
        v_1_2 = -(d_e / d_na - e_na)
        voltages.append(v_1_2)
        voltage_coord.append([(run_1.x_na, v_1_2), (run_2.x_na, v_1_2)])

    plot_title = "Predicted Voltage Curve"
    fig = plt.figure(plot_title)
    axe = fig.add_subplot(1, 1, 1)
    axe.add_collection(LineCollection(voltages))
    axe.set_ylabel('V equilibrium')
    axe.set_xlabel('Na content')
    axe.set_ylim([min(voltages) - 1, max(voltages) + 1])

    # axe.legend()

    plt.show(block=False)

    generic_plot.save_fig(fig, plot_title)

    return fig
