#!/usr/bin/env python
# drawKpoints.py
# -*- coding=utf-8 -*-

from __future__ import print_function

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pymatgen as mg
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath

"""
Create a KPOINTS file for a band structure calculation. This script use
methods of pymatgen in order to compute and select high symetry lines
in the first brillouin zone.
SYNTAX
        makeKpoints.py [OPTIONS] [STRUCTURE FILE]
STRUCTURE FILE
        must contain a structure. For example a POSCAR file.
OPTIONS
        -d ndiv
                ndiv is a integer corresponding to the number of
                k-points needed along each symetry line
"""

__author__ = "Germain Salvato-Vallverdu"
__email__ = "germain.vallverdu@univ-pau.fr"
__licence__ = "GPL"
__date__ = "April 2014"






def drawkpt(struct, ndiv="10"):

    # symmetry information
    struct_sym = SpacegroupAnalyzer(struct)
    print("\nLattice details:")
    print("----------------")
    print("lattice type : {0}".format(struct_sym.get_lattice_type()))
    print(
        "space group  : {0} ({1})".format(
            struct_sym.get_space_group_symbol(),
            struct_sym.get_space_group_number()))

    # Compute first brillouin zone

    ibz = HighSymmKpath(struct)
    print("ibz type     : {0}".format(ibz.name))
    [x, y, z] = list(map(list, zip(*ibz.get_kpoints()[0])))

    fig = plt.figure("Brillouin Zone and High Symm Pts")
    ax = fig.gca(projection='3d')
    ax.plot(x, y, z)

    for i, name in enumerate(ibz.get_kpoints()[1]):
        if name != '':
            #print(" name {0} : [{1},{2},{3}]".format(name, x[i],y[i],z[i]))
            ax.text(x[i], y[i], z[i], '%s' % (name),
                    color='k', size="15")

    new_lat = ibz.prim_rec
    bz_array = new_lat.get_wigner_seitz_cell()
    bz_faces = Poly3DCollection(bz_array)

    bz_faces.set_edgecolor('k')
    bz_faces.set_facecolor((0, 1, 1, 0.4))
    ax.add_collection3d(bz_faces)

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    fig.suptitle(
        "Brillouin Zone and K_Path of \n {0}".format(
            struct.get_primitive_structure().formula))

    fig.show()

    # if input("press [s]ave to save brillouin zone figure as it is \n") == "s":
    #     fig.savefig("BZ_KPath.svg", bbox_inches='tight')

    # plt.close(fig)

    # print specific kpoints in the first brillouin zone
    print("\nList of high symmetry k-points:")
    print("-------------------------------")
    for key, val in ibz.kpath["kpoints"].items():
        print("%8s %s" % (key, str(val)))

    # suggested path for the band structure
    print("\nSuggested paths in first brillouin zone:")
    print("----------------------------------------")
    for i, path in enumerate(ibz.kpath["path"]):
        print("   %2d:" % (i + 1), " -> ".join(path))

    # write the KPOINTS file
    print("\nWrite file KPOINTS")
    kpt = Kpoints.automatic_linemode(ndiv, ibz)

    # if input("write kpt in cwd ?") == "Y":
    #     kpt.write_file(os.path.join(
    #                    folder, "linear_KPOINTS"))
    return(kpt)


if __name__ == '__main__':

    # default args
    fstruct = "POSCAR"
    ndiv = 20

    # read args
    if len(sys.argv) > 1:
        if len(sys.argv) == 2:
            fstruct = sys.argv[1]

    if os.path.exists(fstruct):
        struct = mg.Structure.from_file(fstruct)
    else:
        print("file {0} does not exist in {1}".format(fstruct, os.getcwd()))
        exit(1)

    drawkpt(struct)
