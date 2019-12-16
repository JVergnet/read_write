# read_hull.py

# import logging

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.ticker import AutoMinorLocator

# from operator import itemgetter
# import rundict_utils as read
import rw_utils.generic_plot as generic_plot
from filtering_runs import filter_runs

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
        try:
            choice = input("plot simple energy graph ? [Y/N]")
            assert len(choice) > 0, "default to false"
            simple_energy_graph = bool(choice[0] == "Y")
        except Exception:
            simple_energy_graph = False

    if simple_energy_graph:
        sorted_converged_entries = sorted(
            converged_entries, key=lambda x: x.energy_per_fu)
        generic_plot.plot_structure_value_evolution(
            sorted_converged_entries,
            prop_list=['energy_per_fu'],
            legend=None,
            x_na_coords=coord)
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
            sorted_entries = filter_runs.generate_hull_tags(
                sorted_entries, coord=coord)
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

    #  int => stacking string
    stacking_list = list({e.stacking for e in sorted_entries})
    print("stacking set : ", stacking_list)

    # stacking string => int
    stacking_index = {s: i for i, s in enumerate(stacking_list)}
    print("stacking index : ", stacking_index)

    colors = ['green', 'red', 'blue', 'black']

    # stacking string => color
    color_dict = {}
    for i, stacking in enumerate(stacking_list):
        color_dict[stacking] = colors[i]
    print(color_dict)
    color_dict.update({"P3": "#b700ffff",
                       "O3": "#00baffff"})

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
        x_e_1 = np.array([[getattr(e, coord), e.eform] for e in clean_entries])
        # print(XE1)
        x_clean = x_e_1[:, 0]
        e_clean = x_e_1[:, 1]
        yerr_down = [(thermal_error if (e.status < 4) else 0)
                     for e in clean_entries]
        yerr_up = thermal_error * np.ones_like(e_clean)
        # print([yerr_down, yerr_up])
        axe.errorbar(
            x_clean, e_clean,
            yerr=[yerr_down, yerr_up],
            fmt='',
            linestyle="None")

    # # Plot all entries (hollow circles) [sorted_entries]
    for entry in sorted_entries:
        try:
            print(entry.eform)
        except Exception as ex:
            print("Exception raised for {} : {}".format(entry.name_tag, ex))
    for stacking in stacking_list:
        print("current stacking", stacking)
        x_e = np.array([[getattr(entry, coord), entry.eform]
                        for entry in sorted_entries
                        if entry.stacking == stacking])

        print(x_e)
        axe.scatter(
            x_e[:, 0], x_e[:, 1], s=dot_size,
            facecolors='none',
            edgecolor=color_dict[stacking],
            label=stacking)

    # Plot entries with lowest energies (filled circles) [clean_entries]

    x_e_c = np.array([[getattr(entry, coord),
                       entry.eform,
                       stacking_index[entry.stacking]]
                      for entry in clean_entries])
    print(x_e_c)
    axe.scatter(x_e_c[:, 0], x_e_c[:, 1], s=dot_size,
                color=[color_dict[stacking_list[int(index)]]
                       for index in x_e_c[:, 2]],
                edgecolor='face')

    # link entries on the convex hull (black line) [hull_entries]
    x_e = np.array([[getattr(entry, coord), entry.eform]
                    for entry in hull_entries])
    print(x_e)
    axe.plot(x_e[:, 0], x_e[:, 1], "k-")

    axe.legend()
    axe.set_ylabel('E total (meV)')
    axe.set_xlabel('Na content'
                   )
    # X = [getattr(entry, coord) for entry in sorted_entries]
    if coord == "x_na":
        axe.invert_xaxis()
    # axe.tick_params(axis='x', which='minor', bottom=True)
    axe.xaxis.set_minor_locator(AutoMinorLocator(n=2))
    # xlim(max(X), min(X))

    plt.show(block=False)
    generic_plot.save_fig(fig, plot_title)

    return fig


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
