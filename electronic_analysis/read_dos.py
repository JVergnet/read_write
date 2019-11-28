"DOS PLOTTING FONCTIONS"

import matplotlib.pyplot as plt
import numpy as np
from pymatgen.electronic_structure.core import OrbitalType, Spin

import utils.generic_plot as generic_plot

# import pymatgen.analysis.diffraction.xrd as xrd


# DOS PLOTTING FONCTIONS
# ============================================================================


def DOS_per_site(runDict, spin):

    spec_count = {}
    structure = runDict.structure

    # each type of atom has a differnt colormap
    # so that non-equivalent sites of same element
    # have similar  colors

    # first get the number of sites per element
    for site_dict in runDict.equivSiteList:
        elt_name = site_dict['element']
        if spec_count.get(elt_name, None) is None:
            spec_count[elt_name] = 0
        spec_count[elt_name] += 1

    print(spec_count)

    cmap_colors = ['Blues', 'Greens', 'Purples', 'Reds', 'Greys', 'Oranges']
    cmaps = [plt.get_cmap(s) for s in cmap_colors]
    cmap_list = {}
    for elt_name, cmap in zip(spec_count.keys(), cmaps):
        elt_colors = np.linspace(0.8, 0.45, spec_count[elt_name]) \
            if spec_count[elt_name] > 1 else [0.7]
        cmap_list[elt_name] = cmap(elt_colors)
    densities = {}
    colors = {}
    for site_dict in runDict.equivSiteList:
        elt_name = site_dict['element']
        leg = "{0}:{1} (x{2})".format(elt_name,
                                      site_dict['mainIndex'],
                                      site_dict['multiplicity'])
        site_dos = runDict.data["complete_dos"].get_site_dos(
            structure[site_dict['mainIndex']])
        value = site_dos.densities[spin] * site_dict['multiplicity']
        densities[leg] = value
        spec_count[elt_name] -= 1
        colors[leg] = cmap_list[elt_name][spec_count[elt_name]]

    return(densities, colors)


def DOS_on_index(rundict, spin, index_list):
    densities = {}
    struct = rundict.structure
    colors = {}
    cmap = plt.get_cmap('nipy_spectral')
    color_list = cmap(np.linspace(0.5, 0.1, len(index_list)))
    for index in index_list:
        leg = "atom {} ({})".format(index, struct[index].species_string)
        site_dos = rundict.data["complete_dos"].get_site_dos(
            struct[index]).densities[spin]
        densities[leg] = site_dos
        colors[leg] = color_list[index]

    return (densities, colors)


def DOS_spd(rundict, spin):
    # dict of {orbital:DOS}
    spd_dos = rundict.data["complete_dos"].get_spd_dos()
    densities = {"s": spd_dos[OrbitalType.s].densities[spin],
                 "p": spd_dos[OrbitalType.p].densities[spin],
                 "d": spd_dos[OrbitalType.d].densities[spin]}

    colors = {"s": (1, 0, 0, 1), "p": (0, 0, 1, 1), "d": (0, 1, 0, 1)}  # R,B,G
    return(densities, colors)


def DOS_per_elt(rundict, spin, pre_defined_colors=None):
    elt_dos = rundict.data["complete_dos"].get_element_dos()
    densities = {}
    element_list = elt_dos.keys()
    # print([ e.symbol for e in element_list])
    cmap = plt.get_cmap('nipy_spectral')
    colors = {}

    color_list = [cmap(i) for i in np.linspace(0.1, 0.85, len(element_list))]
    for i, element in enumerate(element_list):
        # print(element)
        densities[element] = elt_dos[element].densities[spin]
        if pre_defined_colors is not None \
           and element.symbol in pre_defined_colors.keys():
            colors[element] = pre_defined_colors[element.symbol]
        else:
            colors[element] = color_list[i]
    return(densities, colors)


def DOS_total(rundict, spin):
    densities = {"total dos": rundict.data["tdos"].densities[spin]}
    colors = {"total dos": (0, 0, 0, 1)}  # black
    return(densities, colors)


def plot_DOS_on_axe(axe, rundict, Emin, Emax,
                    DOS_choice="total",
                    spin_choice=1,
                    doo_projected=False,
                    pre_defined_colors=None,
                    N_move_mean=1):

    # v = runDict['vaspRun']
    [ymin, ymax] = [0, 0]
    e_values = rundict.data["tdos"].energies - rundict.data["efermi"]
    ind_min = 0
    ind_max = 0
    for q, e in enumerate(e_values):
        if e < Emin:
            ind_min = q
        if e < Emax:
            ind_max = q
    # XRD / DOS PLOTTING =====================================================
    # Generation of the graph XRD or DOS of the current structure
    # if param['graph_type'] == 'DOS' :
    # Each DOS plotting defines
    # a dict of { label : [Y values]} (densities)
    # and the corresponding colormap (colors)

    densities_list = []

    for s in [Spin.up, Spin.down]:
       # projection on S, P and D orbitals
        if DOS_choice == "spd":
            (spin_densities, colors) = DOS_spd(rundict, s)

        # Projection on each weighted non-equivalent site
        elif DOS_choice == "perSite":
            (spin_densities, colors) = DOS_per_site(rundict, s)
        # Projection on elements
        elif DOS_choice == "perElt":
            (spin_densities, colors) = DOS_per_elt(rundict, s,
                                                   pre_defined_colors=pre_defined_colors)

        # total DOS (default)
        else:
            (spin_densities, colors) = DOS_total(rundict, s)

        # list of dict {label : Y values of dos }
        # spin_array=[[Spin.up],[Spin.down],[Spin.up,Spin.down]]

        # plot projected orbitals on O in peroxo like bonds
        if doo_projected:
            (index_spin_densities, colors_site) = DOS_on_index(
                rundict, s, rundict.dOO_min_indices)
            spin_densities.update(index_spin_densities)
            colors.update(colors_site)

        densities_list.append(spin_densities)
    # print("densities : \n {} \n colors : \n {}"
    #       .format(densities_list,colors ))

    def move_mean(x, n_input):
        if n_input == 1:
            return(x)
        else:
            half = np.ceil(n_input) // 2
            n = half * 2 + 1
            cumsum = np.cumsum(np.insert(x, 0, 0))
            s = (cumsum[n:] - cumsum[:-n]) / float(n)
            s = np.insert(s, len(s), x[-1]*np.ones(half))
            s = np.insert(s, 0, x[0]*np.ones(half))
            return(s)

    if spin_choice == 0:  # spin up & down
        for n, key in enumerate(densities_list[0].keys()):
            dos_up = move_mean(densities_list[0][key], N_move_mean)
            axe.plot(e_values, dos_up, color=colors[key], label=key, lw=2)
            dos_down = -move_mean(densities_list[1][key], N_move_mean)
            axe.plot(e_values, dos_down, color=colors[key], lw=2)

            ymax = max(ymax, max(dos_up[ind_min:ind_max]))
            ymin = min(ymin, min(dos_down[ind_min:ind_max]))

    else:
        if spin_choice is 1:  # sum of up + down
            a = 1
        elif spin_choice is 2:  # DIFF  of up - down
            a = -1

        for n, key in enumerate(densities_list[0].keys()):
            dos = move_mean(
                densities_list[0][key] + a * densities_list[1][key],
                N_move_mean)

            ymax = max(ymax, max(dos[ind_min:ind_max]))
            ymin = min(ymin, min(dos[ind_min:ind_max]))

            axe.plot(e_values, dos, color=colors[key], label=key, lw=2)

    axe.axvline(x=0)  # the fermi lvl
    axe.set_ylim([ymin * 1.1, ymax * 1.1])
    axe.set_xlim([Emin, Emax])
    axe.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        right=False,      # ticks along the bottom edge are off
        left=False,         # ticks along the top edge are off
        labelleft=False)  # labels along the bottom edge are off

    if DOS_choice == 'perSite':
        axe.legend(loc="upper right", ncol=2,
                   bbox_to_anchor=(1, 1),
                   bbox_transform=axe.transAxes)

    # # generation of XRD patterns
    # elif param['graph_type'] == 'XRD' :

    #     XRDCalc=xrd.XRDCalculator()
    #     XRDCalc.get_xrd_plot(structure,annotate_peaks=False, ax=axe, with_labels=False)

    # Print the name and energy Tags on the left of the graph
    legend = "{} Na{}\n({})".format(
        rundict.stacking, rundict.x_na, rundict.str_id)
    axe.text(-0.03, 0.5, legend,
             horizontalalignment='right', verticalalignment='center',
             multialignment='right', transform=axe.transAxes)
    # Print the structure tag  on the right of the graph
    axe.text(1.03, 0.5, rundict.spacegroup,
             horizontalalignment='left', verticalalignment='center',
             multialignment='left', transform=axe.transAxes,
             fontsize=8)


def plot_DOS(rundict_list,
             spin_choice=None,
             DOS_choice=None,
             Erange=None,
             N_move_mean=None):
    "plot DOS or XRD for each structure in the vasprundictList"

    nbRun = len(rundict_list)

    # OPTIONS
    # ========================

    print("Plotting DOS")

    # Spin choice
    if spin_choice is None:
        spin_choice = 0
        try:
            spin_choice = int(
                input("0 = up & down , 1 = summ, 2 = diff  \nSpin type ? "))
        except Exception:
            print("spin to default : spin up & down")

    # graph type choice
    dos_options = ["spd", "perSite", "perElt", "total"]
    while DOS_choice is None:
        try:
            DOS_choice = dos_options[int(input(
                "DOS type option : {0} \nEnter DOS type number : \n"
                .format(dos_options)))]
        except Exception as ex:
            print("Bad Option", ex)

    doo_projected = False
    # if input("Plot min doo projection? Y / N  : ")=="Y" :
    #     doo_projected = True

    # energy range
    if Erange is not None:
        [Emin, Emax] = Erange
    else:
        [Emin, Emax] = [-3, 3]
        if input("change default energy range [{},{}] ? Y / N "
                 .format(Emin, Emax)) == "Y":
            Emin = float(input("Emin : "))
            Emax = float(input("Emax : "))
    if N_move_mean is None:
        try:
            N_move_mean = int(
                int(input("running mean window ? (default to None)")))
        except Exception as ex:
            N_move_mean = 1

    # PLOTTING
    # =====================================
    #  initialization of the plotting figure

    plotTitle = "DOS " + DOS_choice + " spin " + \
        ["up & down", "sum", "diff"][spin_choice]

    fig = plt.figure(plotTitle, figsize=(20, 10))

    fig.suptitle(rundict_list[-1].job_folder + "\n" +
                 rundict_list[-1].formula + "\n"
                 + plotTitle,
                 fontsize="large")

    fig = plt.figure(plotTitle, figsize=(20, 20))
    # nbRun = len(rundict_list)
    if nbRun == 1:
        axes = [fig.add_subplot(111)]
    else:
        axes = fig.subplots(nbRun, 1, sharex="col")
        fig.subplots_adjust(hspace=0)

    # defining common colors
    elements = []
    for run in rundict_list:
        elements += run.structure.composition.get_el_amt_dict().keys()
    elements = set(elements)
    print("elements : {}".format(elements))

    cmap = plt.get_cmap('tab10')
    pre_defined_colors = {}
    for element, x_elem in \
            zip(elements, np.linspace(0.04, 0.95, len(elements))):
        pre_defined_colors[element] = cmap(x_elem)

    # pre_defined_colors = {"O": "xkcd:light red",
    #                       "Mn": "xkcd:vibrant purple",
    #                       "Mg": "xkcd:yellow orange",
    #                       "Na": "xkcd:vibrant green"}

    for i, runDict in enumerate(rundict_list):
        plot_DOS_on_axe(axes[(nbRun - 1) - i],
                        runDict,
                        Emin,
                        Emax,
                        DOS_choice=DOS_choice,
                        spin_choice=spin_choice,
                        doo_projected=doo_projected,
                        pre_defined_colors=pre_defined_colors,
                        N_move_mean=N_move_mean)
    # If all the graphs have the same legend (ie for all options except "dos per site")
    # print only one legend at the top left side of the figure
    if DOS_choice != 'perSite':
        axes[nbRun - 1].legend(loc="lower left", bbox_to_anchor=(1, 0.6),
                               bbox_transform=axes[nbRun - 1].transAxes)

    # plt.show(block=False)
    # read.save_fig(fig,plotTitle)

    return(fig)


def plot_DOS_graphs(sorted_entries):
    print("\n==== ELECTRONS IN RECIPROCAL SPACE (DOS) ======\n")

    if input("plot DOS  of some structures ? Y / N  : ") == "Y":
        # param['graph_type']="DOS"
        fig = plot_DOS(sorted_entries, spin_choice=None, DOS_choice=None)
        plt.show(block=False)
        generic_plot.save_fig(fig, "DOS", save_mode=None, folder=None)

    # if input("plot bandgap evolution ? Y / N  : ") == "Y":
    #     generic_plot.plot_structure_value_evolution(
    #         sorted_entries, prop_list=['bandgap'])
