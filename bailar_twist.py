import importlib
import itertools
# import launchDisordered as launch
# import structure_geometry_utils as cluster
# import math
from multiprocessing import Pool, cpu_count
# import subprocess
# import sys
# import os
from operator import itemgetter

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pymatgen.util.coord as util_coord
from matplotlib.ticker import FormatStrFormatter
# from pymatgen.core.structure import Structure
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import \
    LocalGeometryFinder

import DOS_plot as DOS
import generic_plot as generic_plot
import readBader as bader
# import lobster_coop as lob
import readRun_entries as read

try:
    import PES_complex_plot as PES_plot
except Exception as ex:
    print(ex)

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
mpl.rcParams['axes.labelsize'] = 17
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rc('legend', fontsize=15)


global species
species = ["Na", "Mn", "Mg"]


global coords
coords = ['T:6', 'O:6']


def get_chem_env_tags(runList):
    vasp_run_poll = []
    # for d in runList :
    #     vasp_run_poll.append(get_chem_env_tag_single(d))

    with Pool(processes=cpu_count()) as p:
        vasp_run_poll = p.map(get_chem_env_tag_single, runList)
        p.close()
        p.join()

    return(vasp_run_poll)


def get_chem_env_tag_single(dict_car):

    print("\n\nnext structure : {}".format(dict_car.folder))

    if dict_car.data.get("distortion", None) is not None:
        print("already done, skipping this one ! ")
        return(dict_car)

    S = dict_car.structure
    L = LocalGeometryFinder(
        permutations_safe_override=False,
        plane_ordering_override=False,
        debug_level=None,
        plane_safe_permutations=False,
        only_symbols=[
            "O:6",
            "T:6"])

    L.setup_structure(S)
    C = L.compute_structure_environments(
        max_cn=6, min_cn=6, additional_conditions=[
            1, 4], only_symbols=[
            "O:6", "T:6"])

    print("\n{} : OP disto ".format(dict_car["folder"]))

    dict_car.data["distortion"] = {}

    for specie in species:
        dict_car.data["distortion"][specie] = {}
        # print("\nspecie {} ! ".format(specie))
        sites = S.indices_from_symbol(specie)
        if len(sites) == 0:
            continue

        for geom in coords:
            n = 0
            csm = 0
            for site in sites:
                try:
                    csm += C.get_csm(site, geom)['symmetry_measure']
                    n += 1
                except Exception as ex:
                    print(
                        "{} geom for {} at site {} couldn't be computed ".format(
                            geom, specie, site))
            print("number of valid sites : {}".format(n))
            if n > 0:
                dict_car.data["distortion"][specie][geom] = csm / n

    return(dict_car)


def plot_CSM(struct_list, stacking_name=None, axe0=None):

    if axe0 is None:
        if stacking_name is None:
            stacking_name = struct_list[0].stacking
        fig = plt.figure("dist against x {}".format(stacking_name))
        axe = fig.add_subplot(111)
    else:
        axe = axe0

    print("x_na :{}".format([dict_car.x_na for dict_car in struct_list]))
    for specie in species:
        for coord in coords:
            XY = np.array([[dict_car.x_na,
                            dict_car.data["distortion"][specie].get(coord, 0)]
                           for dict_car in struct_list
                           if dict_car.data["distortion"][specie].get(coord, None)
                           is not None])
            X = XY[:, 0]
            Y = XY[:, 1]
            legend = "{} : dist from {}".format(specie, coord)
            axe.plot(X, Y, label=legend, marker="*")
            print("{}, {} :{}".format(specie, coord, Y))

    if axe0 is None:
        fig.legend()
        return(fig)
    else:
        return(axe)
    # plt.show(block = False )
    # read.save_fig(fig,"disto_{}".format(stacking),save_mode=None, folder=None)


def plot_bailar_path(
        struct_list,
        stacking_name=None,
        axe0=None,
        add_disto_path=True,
        annotate=True):

    if axe0 is None:
        if stacking_name is None:
            stacking_name = struct_list[0]["stacking"]
        fig = plt.figure("bailar path {}".format(stacking_name))
        axe = fig.add_subplot(111)
    else:
        axe = axe0

    for specie in species:
        XY = np.array([[dict_car.data["distortion"][specie].get(c, 0)
                        for c in coords]
                       for dict_car in struct_list
                       if None not in [dict_car.data["distortion"][specie].get(c, None)
                                       for c in coords]])
        X = XY[:, 0]
        Y = XY[:, 1]

        # [ [dict_car["distortion"][specie].get(coord, 0) for dict_car in struct_list ] for coord in coords]
        axe.scatter(X, Y, label="{} ({})".format(specie, stacking_name))
        if annotate:
            tags = [dict_car["x_na"] for dict_car in struct_list]
            for label, x, y in zip(tags, X, Y):
                axe.annotate(
                    label,
                    xy=(x, y), xytext=(20, 20),
                    textcoords='offset points', ha='right', va='bottom',
                    bbox=dict(boxstyle='round,pad=0.1',
                              fc='yellow', alpha=0.2),
                    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

    if add_disto_path:
        bailarX = np.linspace(0, 17.3, num=50)
        bailarY = np.square(4.16 - np.sqrt(bailarX))
        axe.plot(bailarX, bailarY, label="bailar path")

        TP_X = np.linspace(0, 5, num=20)
        YcompTP = 0.75 * TP_X + 16.40
        YelongTP = 0.97 * TP_X + 16.88
        # axe.plot(TP_X,YcompTP, label="compressed TP")
        axe.plot(TP_X, YelongTP, label="elongated TP")

        JT_X = np.linspace(17, 22, num=20)
        JT_Y = 1.23 * JT_X - 20.52
        axe.plot(JT_X, JT_Y, label="Jahn-Teller")

        ETAP_X = np.linspace(16.5, 22, num=20)
        ETAP_Y = 2.03 * ETAP_X - 33
        axe.plot(ETAP_X, ETAP_Y, label="compressed TAP")

    axe.set_xlabel('$CSM(D_{3H})$')
    axe.set_ylabel('$CSM(O_{H})$')

    if axe0 is None:
        fig.legend()
        return(fig)
    else:
        return(axe)
    # plt.show(block = False )

    # read.save_fig(fig,"bailar_{}".format(stacking),save_mode=None,
    # folder=None


def plot_cell_param(struct_list, ax0=None):
    colors = ['green', 'red', 'blue', 'black']
    marks = ['^', 'o', '*', '+']

    if ax0 is not None:
        axe = ax0
    else:
        fig = plt.figure()
        axe = fig.add_subplot(111)

    # plot cell params  as a function of the distortion and x_na
    X = np.array([d.x_na for d in struct_list])

    YABC = np.array([list(d.structure.lattice.lengths_and_angles[0])
                     for d in struct_list])

    for i, direction in enumerate(["a", "b", "c"]):
        color = colors[i]
        mark = marks[i]
        axe.plot(X, YABC[:, i], c=color, marker=mark, label=direction, ls="--")

    axe.set_xlabel("$X_{Na}$")
    axe.set_ylabel("Super-Cell parameter\n(A)")

    if ax0 is None:
        fig.legend()
        plt.show(block=False)
        generic_plot.save_fig(fig, "cell parameters",
                              save_mode=None, folder=None)
        return(fig)
    else:
        return(axe)


def plot_cell_volume(struct_list, ax0=None):

    # plot  V of the cell as a function of the distortion and x_na
    if ax0 is not None:
        axe = ax0
    else:
        fig = generic_plot.plot_structure_value_evolution(struct_list,
                                                          prop_list=['volume'],
                                                          x_na_coords=False)
        # fig = plt.figure()
    #     axe = fig.add_subplot(111)

    # X = np.array([d["x_na"] for d in struct_list])
    # YV = np.array([d["structure"].lattice.volume / d["nb_cell"]
    #                for d in struct_list])

    # axe.plot(X, YV, c="black", marker="o", ls="--")  # label="Cell Volume" ,

    # axe.set_xlabel("$X_{Na}$")
    # axe.set_ylabel("Cell volume\n(A³/F.U.)")

    if ax0 is None:
        fig.legend()
        plt.show(block=False)
        generic_plot.save_fig(fig, "volume", save_mode=None, folder=None)
        return(fig)
    else:
        return(axe)


def plot_OO_compression(struct_list, ax0=None):

    Y_all = []
    X_all = []

    for d in struct_list:
        struct = d.structure
        Y = [abs(util_coord.pbc_diff(struct[oo[0]].frac_coords,
                                     struct[oo[1]].frac_coords)[2])
             * struct.lattice.c for oo in d.OO_pairs]

        X = np.ones(len(d.OO_pairs)) * d.x_na

        Y_all.append(Y)
        X_all.append(X)

    # print(Y_all)
    Y_all = np.vstack(Y_all)
    X_all = np.vstack(X_all)
    # print(Y_all)
    print(Y_all.shape, X_all.shape)
    Y_min_list = np.amin(Y_all, axis=1)
    Y_max_list = np.amax(Y_all, axis=1)
    Y_mean_list = np.sum(Y_all, axis=1) * (1 / len(d.OO_pairs))
    X_mean = np.array([d.x_na for d in struct_list])
    print(Y_mean_list.shape, X_mean.shape)
    # plot  V of the cell as a function of the distortion and x_na
    if ax0 is not None:
        axe = ax0
    else:
        fig = plt.figure()
        axe = fig.add_subplot(111)
    if True:
        #    if "per_site" in plot_type :
        axe.scatter(X_all.flatten(), Y_all.flatten(), color="blue",
                    s=10)  # , label = "non-equivalent sites")
    # if "average" in plot_type :
        axe.plot(X_mean, Y_mean_list, "ro-")  # ,  label = "weighted mean")
    # if "min" in plot_type :
        # axe.plot(X_mean ,Y_min_list ,"go--",  label = "minimum")
    # if "max" in plot_type :
        # axe.plot(X_mean ,Y_max_list ,"go--",  label = "maximum")

    axe.set_xlabel("$X_{Na}$")
    axe.set_ylabel("short O-O distance\nalong C direction (A)")
    if ax0 is None:
        fig.legend()
        plt.show(block=False)
        generic_plot.save_fig(fig, "compression", save_mode=None, folder=None)
        return(fig)
    else:
        return(axe)


def new_dihedral(p):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 = b1 / np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def measure_quadruplet_angles(struct, MMOO):
    print(MMOO)

    O1, O2 = MMOO['oxygen_pair']
    XO1 = struct.lattice.get_cartesian_coords(struct[O1].frac_coords)

    d, jimage = struct[O1].distance_and_image(struct[O2], jimage=None)
    XO2 = struct.lattice.get_cartesian_coords(
        struct[O2].frac_coords + np.array(jimage))

    bailar_angle = np.linalg.norm(XO2[0:2] - XO1[0:2])
    trigonal_length = abs(XO1[2] - XO2[2])
    # print( bailar_angle, trigonal_length)
    return(bailar_angle, trigonal_length)


def measure_quadruplet_angles_old(struct, MMOO):
    print(MMOO)

    O1, O2 = MMOO['oxygen_pair']
    M1, M2 = MMOO['metal_pair']
    XM, XO1, XO2 = [get_coord(struct, O1, i) for i in [M1, O1, O2]]
    # define two "virtual points" above and below the metal (en the upper and
    # lower Oxygen planes)
    XA1, XA2 = [np.array([XM[0], XM[1], Xi[2]]) for Xi in [XO1, XO2]]
    # print(XA1,XA2)
    bailar_angle = abs(new_dihedral([XO1, XA1, XA2, XO2]))
    trigonal_length = abs(XO1[2] - XO2[2])
    # print( bailar_angle, trigonal_length)
    return(bailar_angle, trigonal_length)


def get_coord(struct, O1, index):
    d, jimage = struct[O1].distance_and_image(struct[index], jimage=None)
    return(struct.lattice.get_cartesian_coords(
        struct[index].frac_coords + np.array(jimage)))


def measure_quadruplet_angles_old_old(struct, MMOO):

    # abc = np.array(struct.lattice.abc)
    bailar_angles = []
    trigonal_angles = []
    # print(MMOO)
    O1, O2 = MMOO['oxygen_pair']
    M1, M2 = MMOO['metal_pair']

    XM1, XM2, XO1, XO2 = [get_coord(struct, O1, Xi) for Xi in [M1, M2, O1, O2]]
    # print("XM1",XM1,"XM2", XM2,"XO1",XO1, "XO2", XO2)

    XOm = (XO1 + XO2) / 2
    # print(XOm)

    for XM in [XM1, XM2]:
        XO_flat_z = np.array([[x[0], x[1], XOm[2]] for x in [XO1, XO2]])
        # print(XO_flat_z)
        v = XO_flat_z - XM
        # print(v)
        cosine_angle = np.dot(v[0], v[1]) / \
            (np.linalg.norm(v[0]) * np.linalg.norm(v[1]))
        angle = np.arccos(cosine_angle)
        bailar_angles.append(np.degrees(angle))

        XO_flat_xy = np.array([[XOm[0], XOm[1], x[2]] for x in [XO1, XO2]])
        v = XO_flat_xy - XM
        cosine_angle = np.dot(v[0], v[1]) / \
            (np.linalg.norm(v[0]) * np.linalg.norm(v[1]))
        angle = np.arccos(cosine_angle)
        trigonal_angles.append(np.degrees(angle))

    bailar_angle = np.mean(bailar_angles)
    trigonal_angle = np.mean(trigonal_angles)
    # print("bailar {:.1f} trigo  {:.1f}".format(bailar_angle, trigonal_angle) )
    return(bailar_angle, trigonal_angle)


def plot_OO_angles(struct_list, axes0=None):

    if axes0 is not None:
        axes = axes0
    else:
        fig = plt.figure()
        axes = fig.subplots(2, 1, sharex="col")
        fig.subplots_adjust(hspace=0)

    Y_bailar_all = []
    Y_trigo_all = []
    X_all = []

    for d in struct_list:
        struct = d.structure
        Y = np.array([measure_quadruplet_angles(struct, MMOO)
                      for MMOO in d.MMOO_quadruplets])
        Y_bailar = Y[:, 0]
        Y_trigo = Y[:, 1]
        X = np.ones(len(d.MMOO_quadruplets)) * d.x_na

        Y_bailar_all.append(Y_bailar)
        Y_trigo_all.append(Y_trigo)
        X_all.append(X)

    axes[0].set_xlabel("$X_{Na}$")
    axes[0].set_ylabel("bailar angle (in °)")
    axes[1].set_ylabel(r"trigonal compression (in $\AA$)")
    for i, Y_all in enumerate([Y_bailar_all, Y_trigo_all]):
        # print(Y_all)
        Y_all = np.vstack(Y_all)
        X_all = np.vstack(X_all)
        # print(Y_all)
        print(Y_all.shape, X_all.shape)
        # Y_min_list = np.amin(Y_all, axis=1)
        # Y_max_list = np.amax(Y_all, axis=1)
        Y_mean_list = np.sum(Y_all, axis=1) * (1 / len(d.MMOO_quadruplets))
        X_mean = np.array([d.x_na for d in struct_list])
        print(Y_mean_list.shape, X_mean.shape)
        # plot  V of the cell as a function of the distortion and x_na

        if True:
            #    if "per_site" in plot_type :
            axes[i].scatter(X_all.flatten(), Y_all.flatten(), color="blue",
                            s=10)  # , label = "non-equivalent sites")
            #   if "average" in plot_type :
            # ,  label = "weighted mean")
            axes[i].plot(X_mean, Y_mean_list, "ro-")
            #  if "min" in plot_type :
            # axe.plot(X_mean ,Y_min_list ,"go--",  label = "minimum")
            # if "max" in plot_type :
            # axe.plot(X_mean ,Y_max_list ,"go--",  label = "maximum")
    column_name = ["bailar angle", "trigo_angle"]
    row_name = ["{} Na{}\n({})".format(v.stacking, v.x_na, v.str_id)
                for v in struct_list]
    column_content_list = [Y_bailar_all, Y_trigo_all]
    plot_histo_charge(column_name, row_name, column_content_list,
                      plotTitle="Distortion angles")

    if axes0 is None:
        fig.legend()
        plt.show(block=False)
        generic_plot.save_fig(fig, "OO angles", save_mode=None, folder=None)
        return(fig)
    else:
        return(axes)


def get_AB_dist(s, A, B, max_length=1000):
    AB_dist = []
    M = s.distance_matrix

    A_index = np.array(list(s.indices_from_symbol(A)))
    nb_A = len(A_index)
    if A != B:
        B_index = np.array(list(s.indices_from_symbol(B)))
        nb_B = len(B_index)
        for i in range(nb_B):
            for j in range(nb_A):
                if M[A_index[j], B_index[i]] <= max_length:
                    AB_dist.append(M[A_index[j], B_index[i]])
    else:
        for i in range(nb_A):
            for j in range(i + 1, nb_A):
                if M[A_index[j], A_index[i]] <= max_length:
                    AB_dist.append(M[A_index[j], A_index[i]])

    # MO_dist = sorted(AB_dist)
    return(AB_dist)


def plot_AB_dist(vasprun_dict_list, A_list, B_list, max_length, raw=False):
    (AB_name_list, id_list, AB_dist_list) = get_AB_dist_histo(
        vasprun_dict_list, A_list, B_list, max_length, raw=False)
    fig = plot_histo_charge(AB_name_list, id_list, AB_dist_list)
    return(fig)


def plot_angles_histo(
        vasprun_dict_list,
        A_list,
        B_list,
        max_length,
        raw=False):

    (AB_name_list, id_list, AB_dist_list) = get_AB_dist_histo(
        vasprun_dict_list, A_list, B_list, max_length, raw=False)
    fig = plot_histo_charge(AB_name_list, id_list, AB_dist_list)
    return(fig)


def get_AB_dist_histo(
        vasprun_dict_list,
        A_list,
        B_list,
        max_length,
        raw=False):
    # nb_run = len(vasprun_dict_list)
    # print("nb run : ",nb_run)

    # gather the length of the specified bonds for each structure

    id_list = []
    s_list = []
    AB_dist_list = []
    AB_name_list = []

    for i, v in enumerate(vasprun_dict_list):
        id_list.append("{} Na{}\n({})".format(
            v.stacking, v.x_na, v.str_id))
        s_list.append(v.structure)

    for A in A_list:  # "Mn,Mg,O"
        for B in B_list:
            AB_dist_list.append([get_AB_dist(s, A, B, max_length)
                                 for s in s_list])
            AB_name_list.append("{}-{}".format(A, B))
    return(AB_name_list, id_list, AB_dist_list)


def plot_histo_charge(column_name, row_name, column_content_list,
                      plotTitle="Interatomic distances"):
    # content = np.array(column_content_list).T
    # print("data shape : {}".format(content.shape))
    nb_column = len(column_content_list)
    nb_row = len(column_content_list[0])
    print(nb_column, nb_row)
    # plot histograms of bond lengths for each structure (with identical scale)
    fig = plt.figure(plotTitle)  # ,figsize=(20,20))
    axs_init = fig.subplots(nb_row, nb_column, sharex="col")
    fig.subplots_adjust(hspace=0)
    try:
        print(axs_init)
        print(axs_init[0, 0])
        print(axs_init[0][0])
    except Exception as ex:
        print(ex)
    # ax[row][column]
    # 11 12
    # 21 22
    # ax[ [11,12] [21, 22] ]
    if nb_column == 1 and nb_row == 1:
        # ax[ [11] ]
        axs = [[axs_init]]
    elif nb_row == 1:
        # ax[ [11,12] ]
        axs = [axs_init]
    elif nb_column == 1:
        # ax[ [11],[21] ]
        axs = [[ax] for ax in axs_init]
    else:
        axs = axs_init
    try:
        print(axs)
        print("[0][0]", axs[0][0])
        print("[0,0]", axs[0, 0])
    except Exception as ex:
        print(ex)

    for k in range(nb_column):
        axe = axs[0][k]
        axe.text(0.5, 1.1, column_name[k],
                 size='large', transform=axe.transAxes)
        column_content = list(
            itertools.chain.from_iterable(column_content_list[k]))
        # print(column_name[k])
        # print(column_content)
        if len(column_content) > 0:
            x_lim = [min(column_content), max(column_content)]
            for i in range(nb_row):
                dist = column_content_list[k][i]
                if len(dist) != 0:
                    axe = axs[i][k]
                    a, b, c = axe.hist(dist, 30, range=x_lim,
                                       density=True, stacked=True)

    for i in range(nb_row):
        axe = axs[i][0]
        axe.text(-0.15, 0.5, row_name[i],
                 horizontalalignment='right', verticalalignment='center',
                 multialignment='center', transform=axe.transAxes)

    return(fig)


def plot_structure_graphs(runList, chem_env_done):

    # chem_env_done = False
    # sorted_entries = [d for d in run_list if d["status"] >= filter_lvl]
    print("\n==== STRUCTURE GRAPHS ======\n")
    if input(
            "plot cell parameters volume and trigonal compression ? Y / N  : ") == "Y":
        fig = plt.figure()
        axes = fig.subplots(3, 1, sharex="col")
        fig.subplots_adjust(hspace=0)
        try:
            plot_cell_param(runList, ax0=axes[0])
            plot_cell_volume(runList, ax0=axes[1])
            plot_OO_compression(runList, ax0=axes[2])
        except Exception as ex:
            print(ex)
        fig.legend()
        # plt.pause(0.0001)
        # plt.show(block=False)

    if input("plot distances histograms during discharge ? Y / N  : ") == "Y":
        elements = set()
        for run in runList:
            elements.update(
                set(run.structure.composition.get_el_amt_dict().keys()))
        A_list = [e for e in elements]  # if e is not "Na"]
        for chalco in ["O", "S", "F"]:
            if chalco in elements:
                B_list = [chalco]
                break
        max_length = 4
        fig_dist = plot_AB_dist(runList, A_list, B_list, max_length)

        # plt.pause(0.0001)
        # plt.show(block=False)
        # read.save_fig(fig_dist,"Atomic Distances evolution")

    if input("plot angle variation during discharge ? Y / N  : ") == "Y":
        fig_angle = plot_OO_angles(runList)

    if input("plot angles/energy relation  for converged runs  ? Y / N  : ") == "Y":
        converged_runs = [d for d in runList if d.status >= 3]
        importlib.reload(PES_plot)
        fig_done = PES_plot.plot_3d_angle_energy_from_struct_list(
            converged_runs, struct_mesh=None)
        print(fig_done)

    if input("plot properties as function of min dOO ? Y / N  : ") == "Y":
        runList = [r for r in runList if r.status >= 3]
        # np.arange(0.005,0.1,0.005) : # optimal value determined to be a nice
        # figure
        for binsize in [0.065]:

            DM = np.array([[d.dOO_min, d.mag] for d in runList])

            bins = np.arange(DM[:, 0].min(), DM[:, 0].max(), binsize)

            # finding the most stable structure per bin
            stable_struct_list = []
            for x in bins:
                run_in_bin = [d for d in runList if (
                    d.dOO_min >= x and d.dOO_min < (x + binsize))]
                if len(run_in_bin) > 0:
                    run_in_bin = sorted(
                        run_in_bin, key=lambda x: x.energy_per_fu)
                    print("doo : {:.2f} to {:.2f} : \n \
                    most stable : {} Ediff = {:.2f} among {}".format(x,
                                                                     x + binsize,
                                                                     run_in_bin[0]["nameTag"],
                                                                     run_in_bin[0]["ediff"] -
                                                                     run_in_bin[-1]["ediff"],
                                                                     len(run_in_bin)))

                    stable_struct_list.append(run_in_bin[0])

            DM_bin = np.array([[d.dOO_min, d.mag]
                               for d in stable_struct_list])

            fig = plt.figure("doo and mag : binzize of {}".format(binsize))
            fig.tight_layout()

            ax_all = fig.add_subplot(222)
            ax_all.scatter(DM[:, 0], DM[:, 1])
            ax_all.set_xlabel("d$_{oo}$")
            ax_all.set_ylabel("S$_{Z}$")
            ax_all.scatter(DM_bin[:, 0], DM_bin[:, 1], s=180,
                           facecolors='none', edgecolor="orange")  # , label=lab)
            for n in bins:
                ax_all.axvline(x=n)

            ax_hist = fig.add_subplot(224)
            ax_hist.set_xlim(ax_all.get_xlim())
            hist_bins = np.append(bins, bins[-1] + binsize)
            # including the last bin not created by np arrange
            ax_hist.hist(DM[:, 0], bins=hist_bins)

            ax_bin = fig.add_subplot(121)
            ax_bin.scatter(DM_bin[:, 0], DM_bin[:, 1])
            ax_bin.set_xlabel("d$_{oo}$")
            ax_bin.set_ylabel("S$_{Z}$")

            fig_simple = plt.figure("doo and mag")
            fig_simple.tight_layout()
            ax_simple = fig_simple.add_subplot(111)
            ax_simple.scatter(DM_bin[:, 0], DM_bin[:, 1])
            ax_simple.set_xlabel("d$_{oo}$")
            ax_simple.set_ylabel("S$_{Z}$")
            ax_simple.yaxis.set_major_locator(generic_plot.major_locator())
            ax_simple.yaxis.set_minor_locator(generic_plot.minor_locator())
            ax_simple.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            bader.plot_charge_and_mag(stable_struct_list, detailled=None)
            DOS.plot_DOS(stable_struct_list,
                         spin_choice=None, DOS_choice=None)
            plt.show(block=False)

            generic_plot.save_fig(fig_simple, "doo and mag",
                                  save_mode=None, folder=None)

    if chem_env_done:

        if input("plot CSM ? Y / N  : ") == "Y":
            plot_CSM(runList, stacking_name=None, axe0=None)

        if input("plot distortion map ? Y / N  : ") == "Y":
            plot_bailar_path(runList, stacking_name=None,
                             axe0=None, add_disto_path=True, annotate=True)

    return(True)


def plot_doo_mag(runList):

    plt.figure()

    return(True)
