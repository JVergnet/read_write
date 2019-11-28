import os
import sys

import matplotlib.patches as mpatches
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
# from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.electronic_structure.core import Spin  # , OrbitalType
from pymatgen.io.vasp import Vasprun

from structure_analysis.drawKpoints import drawkpt

# import matplotlib
# matplotlib.use('TkAgg')


def rgbline(ax, k, e, red, green, blue, alpha=1.):
    # creation of segments based on
    # http://nbviewer.ipython.org/urls/raw.github.com/dpsanders/matplotlib-examples/master/colorline.ipynb
    # transpose Kpt and energy list to list of 2elt lists (ie pts)
    pts = np.array([k, e]).T.reshape(-1, 1, 2)
    # concatenate pts i to pts i+1
    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)
    nseg = len(k) - 1
    r = [0.5 * (red[i] + red[i + 1]) for i in range(nseg)]
    g = [0.5 * (green[i] + green[i + 1]) for i in range(nseg)]
    b = [0.5 * (blue[i] + blue[i + 1]) for i in range(nseg)]
    a = np.ones(nseg, np.float) * alpha
    # add list of segments with given color & alpha to the ax
    lc = LineCollection(seg, colors=list(zip(r, g, b, a)), linewidth=2)
    ax.add_collection(lc)


def plotband(bandV, kpt_path="KPOINTS"):
    bands = bandV.get_band_structure(kpt_path, line_mode=True)
    pbands = bands.get_projection_on_elements()
    struct = bandV.final_structure

    emin = 0
    emax = 0
    for spin in bands.bands.keys():
        for b in range(bands.nb_bands):
            emin = min(emin, min(bands.bands[spin][b]))
            emax = max(emax, max(bands.bands[spin][b]))

    # emin = max(-10, emin)
    # emax = min(10, emax)

    emin = emin - bands.efermi
    emax = emax - bands.efermi

    contrib = np.zeros((bands.nb_bands, len(bands.kpoints), 3))
    # for b in range(bands.nb_bands):
    #     for k in range(len(bands.kpoints)):
    #         sc = pbands[Spin.up][b][k][name]["s"]**2
    #         pc = pbands[Spin.up][b][k][name]["p"]**2
    #         dc = pbands[Spin.up][b][k][name]["d"]**2
    #         tot = sc + pc + dc
    #         if tot != 0.0:
    #             contrib[b, k, 0] = sc / tot
    #             contrib[b, k, 1] = pc / tot
    #             contrib[b, k, 2] = dc / tot

    species = ['Mn', 'O']

    for b in range(bands.nb_bands):
        for k in range(len(bands.kpoints)):
            specContrib = []
            for specName in species:
                specContrib.append(
                    sum([
                        pbands[s][b][k][specName]**2
                        for s in [Spin.up, Spin.down]
                    ])
                )
            tot = sum(specContrib)
            if tot != 0.0:
                for i in range(len(species)):
                    # print(i)
                    contrib[b, k, i] = specContrib[i] / tot
                for i in range(len(species), 3):
                    # print(i)
                    contrib[b, k, 2] = 0.5

    red_patch = mpatches.Patch(color='red', label=species[0])
    green_patch = mpatches.Patch(color='green', label=species[1])

    fig, ax = plt.subplots()

    fig.suptitle("Band Structure of \n {0}".format(
        struct.get_primitive_structure().formula))

    fig.legend(handles=[red_patch, green_patch],
               labels=[species[0], species[1]])

    # plot bands using rgb mapping
    for b in range(bands.nb_bands):
        rgbline(ax,
                range(len(bands.kpoints)),
                [e - bands.efermi for e in bands.bands[Spin.up][b]],
                contrib[b, :, 0],
                contrib[b, :, 1],
                contrib[b, :, 2])

    index = 0
    branch_list = []
    while index < len(bands.kpoints):
        branch = bands.get_branch(index)[0]
        branch_list.append(branch)
        index = branch['end_index'] + 2

    # print(branch_list)

    point_list = []
    for i, branch in enumerate(branch_list):
        pt_dict = [branch['start_index'], branch['name'].split("-")[0]]
        # print(pt_dict,i)
        if i == 0 or (pt_dict[1] != point_list[-1][1]):
            point_list.append(pt_dict)
        point_list.append([branch['end_index'], branch['name'].split("-")[-1]])

    print(point_list)

    for point in point_list:
        p = ax.vlines(point[0], emin, emax, "k")

    Xticks = ax.set_xticks([p[0] for p in point_list])
    Xtickslabel = ax.set_xticklabels([p[1] for p in point_list])

    drawkpt(struct)

    # fig.show()
    # if input("press [s]ave to save band Structure as it is \n") == "s":
    #     fig.savefig("bands_pet_elt.png", bbox_inches='tight')

    # input("press any key to close all figures")
    # plt.close("all")
    return(fig)


if __name__ == '__main__':
    job_folder = sys.argv[1]
    os.chdir(job_folder)
    print("parsing vasprun in {0}".format(job_folder))
    bandV = Vasprun("vasprun.xml", parse_projected_eigen=True)
    plotband(bandV)
