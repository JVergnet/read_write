# cluster.py
from operator import itemgetter

# import matplotlib.pyplot as plt
import numpy as np
from pymatgen.analysis.local_env import MinimumDistanceNN
# from pymatgen import Structure
from pymatgen.core.operations import SymmOp
from pymatgen.transformations.standard_transformations import \
    RemoveSpeciesTransformation

# import math
# import os
# import readRun_entries as read


def get_min_OO_dist(s):
    o_index = np.array(list(s.indices_from_symbol("O")))
    nb_O = len(o_index)
    min_dist = 1000
    M = s.distance_matrix
    indices = None
    for i in range(nb_O):
        for j in range(i + 1, nb_O):
            #print(M[o_index[i],o_index[j]] , min_dist)
            if min_dist > M[o_index[i], o_index[j]]:
                min_dist = M[o_index[i], o_index[j]]
                indices = [i, j]
    return(min_dist, indices)


def get_Mn_in_sublattice(struct):

    layer = get_layer_indices_by_number(struct, ["Na"], layer_index=1)

    all_Mn = struct.indices_from_symbol("Mn")
    all_Mg = struct.indices_from_symbol("Mg")
    Mn_layer = list(set(all_Mn) & set(layer))
    Mg_layer = list(set(all_Mg) & set(layer))

    print("Mn in layer : {} ".format(Mn_layer))

    lattice_1 = []
    lattice_2 = []
    lattice_3 = Mg_layer

    new_Mn_1 = [Mn_layer[0]]
    while len(new_Mn_1) > 0:
        Mn = new_Mn_1.pop()
        print("central Mn : {}".format(Mn))
        lattice_1.append(Mn)
        for next_M in get_metal_next_neighbor(struct, Mn):
            if next_M in Mn_layer and next_M not in lattice_2 + lattice_3:
                print("next M :  : {}".format(next_M))
                lattice_2.append(next_M)
                new_Mn_1 += [
                    M for M in get_metal_next_neighbor(
                        struct,
                        next_M) if M not in lattice_1 +
                    new_Mn_1 and M in Mn_layer]

    return(lattice_1, lattice_2)


def get_MO_pairs(struct, metal_str="Mn"):
    print("getting {} - O pairs".format(metal_str))
    atom_pair_list = []
    mini = MinimumDistanceNN(tol=0.3)
    for center_index in struct.indices_from_symbol(metal_str):
        neighbor_O = mini.get_nn_info(struct, center_index)
        atom_pair_list += [[center_index, O['site_index']] for O in neighbor_O
                           if O['site'].species_string in ["O", "S"]]
    return(atom_pair_list)


def get_metal_next_neighbor(struct, center_index):
    next_metal_neighbors = []

    # gather the indices of all metal next neighbors
    # metals around neighbor of neighbor Os
    mini = MinimumDistanceNN(tol=0.3)
    neighbor_O = mini.get_nn_info(struct, center_index)
    for O_index in [O['site_index'] for O in neighbor_O]:
        next_metal_neighbors += [M['site_index']
                                 for M in mini.get_nn_info(struct, O_index)
                                 if (M['site'].specie.symbol not in ["Na", "Li", "O"])]

    # remove the doubles and the metal center
    next_metal_neighbors = set(next_metal_neighbors)
    next_metal_neighbors.discard(center_index)

    return(next_metal_neighbors)


def get_Mn_Mn_pairs(struct):

    # returns all the Mn-Mn pairs between Mn from the [0,0,0] image and any other image
    # output : list of (index_Mn_A (in 0,0,0) index_Mn_B , image of Mn_B )
    Mn_Mn_pairs = []
    Mn_left = set(struct.indices_from_symbol("Mn"))
    for current_Mn in struct.indices_from_symbol("Mn"):
        neighbor_Mn_list = [
            n for n in get_metal_next_neighbor(
                struct, current_Mn) if n in Mn_left]

        Mn_Mn_pairs += [[current_Mn, neighbor_Mn]
                        for neighbor_Mn in neighbor_Mn_list]

        #print("found Mn Mn pair : {} - {} ".format( current_Mn,neighbor_Mn ))
        Mn_left.discard(current_Mn)

    return(list(Mn_Mn_pairs))


def get_OO_pairs(struct):

    AB_list = get_Mn_Mn_pairs(struct)
    print("nb of pairs to test : {}".format(len(AB_list)))
    O_found = set()
    OO_pairs = set()

    for A, B in AB_list:
        #print("A : {} , B : {}".format(A,B))
        # getting the common O in the neighboorhods of the two metals

        mini = MinimumDistanceNN(tol=0.35)

        O_A_indices = set([O['site_index']
                           for O in mini.get_nn_info(struct, A)])
        #print("found O around A : {}".format( O_A_indices) )

        O_B_indices = set([O['site_index']
                           for O in mini.get_nn_info(struct, B)])
        #print("found O around B : {}".format( O_B_indices) )

        O_pair = O_A_indices & O_B_indices
        #print("found O around A still untreated : {}".format(O_A_indices) )

        nb_O = len(O_pair)

        if nb_O > 2:
            O_pair -= O_found

            O_pair_dict = [O for O in mini.get_nn_info(
                struct, B) if O['site_index'] in O_pair]
            # print("found O pair : {}".format([ (dictA["image"],dictA["site_index"]) for dictA in O_pair]))

            # O_pair = [ O for O in O_pair if O["site_index"] not in O_found]

            print("common O between atoms {} are not a pair : {}".format(
                [A, B], O_pair))
            #print("B site {} and image : {} ".format(struct[B],np.array(site_B["image"]) )  )
            for O in O_pair_dict:
                #jimg = np.array([site_B["image"][i] for i in [0,1] ]+[0] )
                dA, imgA = O["site"].distance_and_image(struct[B], jimage=None)

                dB, imgB = O["site"].distance_and_image(struct[A], jimage=None)
                # print(imgA,imgB)

                O["hash"] = hash((tuple(imgA), tuple(imgB)))
                #print("O site : {}, dist : {}".format(O["site_index"],dA+dB))

            O_pair = [O["site_index"]
                      for O in sorted(O_pair, key=itemgetter('hash'))]
            print("sorted O list : {}".format(
                [dictA["site_index"] for dictA in O_pair]))

        elif nb_O < 2:
            print("Not enough O in the pair !!")
            continue

        while len(O_pair) > 0:
            O1 = O_pair.pop()
            O2 = O_pair.pop()
            O_found.add(O1)
            O_found.add(O2)
            OO_pairs.add((O1, O2))

    print("Valid OO pairs : {}".format(len(OO_pairs)))
    return(list(OO_pairs))


def get_MMOO_quadruplets(struct):

    AB_list = get_Mn_Mn_pairs(struct)
    print("nb of pairs to test : {}".format(len(AB_list)))

    quadruplets = []

    for A, B in AB_list:
        mini = MinimumDistanceNN(tol=0.35)

        O_A_indices = [O['site_index'] for O in mini.get_nn_info(struct, A)]
        O_pair = [O['site_index'] for O in mini.get_nn_info(struct, B)
                  if O["site_index"] in O_A_indices]

        # O_pair = O_A_indices & O_B_indices

        nb_O = len(O_pair)

        if nb_O == 2:
            quadruplets.append({"metal_pair": [A, B], "oxygen_pair": O_pair})

        elif nb_O > 2:

            O_pair_dict = [O for O in mini.get_nn_info(
                struct, B) if O['site_index'] in O_pair]

            print("common O between atoms {} are not a pair : {}".format(
                [A, B], O_pair))

            for O in O_pair_dict:

                dA, imgA = O["site"].distance_and_image(struct[B], jimage=None)

                dB, imgB = O["site"].distance_and_image(struct[A], jimage=None)
                # print(imgA,imgB)

                O["hash"] = hash((tuple(imgA), tuple(imgB)))
                #print("O site : {}, dist : {}".format(O["site_index"],dA+dB))

            O_pair = [O["site_index"]
                      for O in sorted(O_pair, key=itemgetter('hash'))]
            print("sorted O list : {}".format(
                [dictA["site_index"] for dictA in O_pair]))
            while len(O_pair) > 0:
                O1 = O_pair.pop()
                O2 = O_pair.pop()
                quadruplets.append(
                    {"metal_pair": [A, B], "oxygen_pair": [O1, O2]})

        else:
            print("Not enough O in the pair !!")

    print("Valid MMOO quadruplets : {}".format(len(quadruplets)))
    return(quadruplets)


def contract_OO_pairs(s, O_pair, bailar_disto=0, trigonal_disto=0):
    if trigonal_disto == 0 and bailar_disto == 0:
        print("OO distortion impossible :  Bailar = trigonal = 0")
        return(s)

    # disto : amount of distortion : 0 ==> no distortion, 1 ==> d(OO) = 0
    # bailar along A and B and trigonal along C
    disto_list = np.array([bailar_disto, bailar_disto, trigonal_disto])

    # compressing the O_O distance
    O1, O2 = O_pair
    struct = s.copy()

    #abc = np.array(struct.lattice.abc)

    d, jimage = struct[O1].distance_and_image(struct[O2], jimage=None)

    X1 = struct[O1].coords
    X2 = struct.lattice.get_cartesian_coords(
        struct[O2].frac_coords + np.array(jimage))
    Xo_list = [X1, X2]

    Xm = [(Xi + Xj) / 2 for (Xi, Xj) in zip(X1, X2)]

    # for both of the two oxygens :
    for i, Xo in enumerate(Xo_list):
        # contraction towards the barycenter weighted by "disto" in that
        # direction
        translation = np.array([(disto) * (Xj - Xi)
                                for (Xi, Xj, disto) in zip(Xo, Xm, disto_list)])

        struct.translate_sites([O_pair[i]], translation,
                               frac_coords=False, to_unit_cell=True)

    return(struct)


def metalBondCount(struc):

    # Counts the metal metal bounds in a structure
    # argument pymatgen.core.structure
    # return [ normalizedBonds, normalizedAA, normalizedBB]

    # normalizedBond = 3 elt list =[ nb AA bonds, nb BB bonds, nb  AB bonds ] / nb of metal
    # normalizedAA = 7 element list : the ith elt is the number of A with i A-neighbor
    # (ie if there are 1/4 of the Mg surrounded by 3 Mg, normalizedAA[3] = 0.25)
    # normalizedBB = idem for B

    # Logically, if nbA = nbB, normalizedAA = normalizedBB

    A = "Mg"
    B = "Mn"
    countAA = 0
    countAB = 0
    countBB = 0

    neighborAA = [0 for i in range(7)]
    neighborBB = [0 for i in range(7)]

    A_indices = struc.indices_from_symbol(A)
    B_indices = struc.indices_from_symbol(B)
    nbA = len(A_indices)
    nbB = len(B_indices)
    nbTot = nbA + nbB
    print("{0} {1}, {2} {3} : {4}".format(nbA, A, nbB, B, nbTot))

    sCopy = RemoveSpeciesTransformation(
        ["Na", "O"]).apply_transformation(struc)

#    struc.remove_species('O')
    mini = MinimumDistanceNN(tol=0.3)

    for i, site in enumerate(sCopy):
        # print(str(i))
        neigList = mini.get_nn_info(sCopy, i)
        #print("{0} closest neighboors \n".format(len(neigList)))
        coord = 0

        for neighbor in neigList:
            index = neighbor['site_index']
            name = neighbor['site'].species_string
 #           print(" ( {0} at {1} ) ".format(name,index))

            if site.species_string == A:
                if name == A:
                    countAA += 1
                    coord += 1
                if name == B:
                    countAB += 1
            if site.species_string == B:
                if name == A:
                    countAB += 1
                if name == B:
                    countBB += 1
                    coord += 1

        if site.species_string == A:
            neighborAA[coord] += 1
        elif site.species_string == B:
            neighborBB[coord] += 1

    bonds = [countAA // 2, countBB // 2, countAB // 2]
    normalizedBonds = [i / nbTot for i in bonds]
    print("[{3}{3}, {4}{4}, {3}{4}] = {0} normalized by {1} : {2} \n".format(
        bonds, nbTot, normalizedBonds, A, B))
    if nbA > 0:
        normalizedAA = [i / nbA for i in neighborAA]
        print(
            "[Nombre de {0} avec [1-6] {0} autour] = {1} normalized by {2} : {3} \n".format(
                A,
                neighborAA,
                nbA,
                normalizedAA))
    else:
        normalizedAA = neighborAA
    if nbB > 0:
        normalizedBB = [i / nbB for i in neighborBB]
        print(
            "[Nombre de {0} avec [1-6] {0} autour] = {1} normalized by {2} : {3} \n".format(
                B,
                neighborBB,
                nbB,
                normalizedBB))
    else:
        normalizedBB = neighborBB
    return [normalizedBonds, normalizedAA, normalizedBB]


def alkaliBondCount(struc):

    # Counts the Na environment in a structure
    #argument : pymatgen.core.structure
    # return : 7 element list : the i th element is the propotion of Na with i
    # A neighbors

    A = "Mg"
    B = "Mn"
    C = "Na"

    neighborC = [0 for i in range(0, 7, 1)]
    C_indices = struc.indices_from_symbol(C)
    nbC = len(C_indices)

    #print("{0} {1}".format(nbC,C))

    sCopy = RemoveSpeciesTransformation(["O"]).apply_transformation(struc)
    # print(sCopy)
    mini = MinimumDistanceNN(tol=0.3)

    for i, site in enumerate(sCopy):
        # print(str(i))
        # print(site)
        siteName = site.species_string
        if siteName == C:
            neigList = mini.get_nn_info(sCopy, i)
            #print("{0} closest neighboors \n".format(len(neigList)))
            coordA = 0
            coordB = 0
            for neighbor in neigList:
                index = neighbor['site_index']
                neighborName = neighbor['site'].species_string
                #print(" ( {0} at {1} ) ".format(neighborName,index))
                if neighborName == A:
                    coordA += 1
                if neighborName == B:
                    coordB += 1

            neighborC[coordA] += 1

    #print("neighborC list :" , neighborC)
    if nbC == 0:
        normalizedNeighborC = [0 for i in neighborC]
    else:
        normalizedNeighborC = [i / nbC for i in neighborC]
        print("coordination {0} : {1} \nnormalized by {2} : {3} \n".format(
            C, neighborC, nbC, normalizedNeighborC))

    return normalizedNeighborC


def metal_envt(struc):

    # Counts the Me environment in a structure
    #argument : pymatgen.core.structure
    # return : 2elt list of (7 element list) : the i th element is the
    # propotion of A with i A neighbors

    A = "Mg"
    B = "Mn"
    C = "Na"

    neighborA = [0 for i in range(0, 7, 1)]
    neighborB = [0 for i in range(0, 7, 1)]
    A_indices = struc.indices_from_symbol(C)
    nbC = len(C_indices)

    #print("{0} {1}".format(nbC,C))

    sCopy = RemoveSpeciesTransformation(["O"]).apply_transformation(struc)
    # print(sCopy)
    mini = MinimumDistanceNN(tol=0.3)

    for i, site in enumerate(sCopy):
        # print(str(i))
        # print(site)
        siteName = site.species_string
        if siteName == C:
            neigList = mini.get_nn_info(sCopy, i)
            #print("{0} closest neighboors \n".format(len(neigList)))
            coordA = 0
            coordB = 0
            for neighbor in neigList:
                index = neighbor['site_index']
                neighborName = neighbor['site'].species_string
                #print(" ( {0} at {1} ) ".format(neighborName,index))
                if neighborName == A:
                    coordA += 1
                if neighborName == B:
                    coordB += 1

            neighborC[coordA] += 1

    #print("neighborC list :" , neighborC)
    if nbC == 0:
        normalizedNeighborC = [0 for i in neighborC]
    else:
        normalizedNeighborC = [i / nbC for i in neighborC]
        print("coordination {0} : {1} \nnormalized by {2} : {3} \n".format(
            C, neighborC, nbC, normalizedNeighborC))

    return normalizedNeighborC

    return good_Mn_list


def get_layer_indices(struct, separator_specie, refIndex=-1):

    # fonction to gather the indexes of all the atoms in layers separated by Na atoms
    # input : the layered structure, the index of one non-Na reference atom
    # returns : an ordered list of indexes of all the atoms in the layer of
    # the reference atom

    A = "Mn"
    B = "Mg"
    C = "O"

    # if no reference is given, takes the first atom of A as a reference
    if refIndex == -1:
        layerIndices = [struct.indices_from_symbol(C)[0]]
    else:
        layerIndices = [refIndex]

    mini = MinimumDistanceNN(tol=0.5)

    # algorithm : 2 buffer list (prev and new atomList) and a definitive list
    # adds all neighbor of the old buffer  atoms to a  new buffer list
    # if they are non-Na and non-already counted in the definitive List
    # Stops when the new buffer is empty (all non-Na atoms next to the layer
    # have been counted)

    prevAtomList = layerIndices
    finish = False
    while not finish:
        newAtomList = []
        for prevAtomIndex in prevAtomList:
            neigList = mini.get_nn_info(struct, prevAtomIndex)
            for neighbor in neigList:
                index = neighbor['site_index']
                name = neighbor['site'].species_string
                # print(name)
                if not (
                        name in separator_specie + [
                            struct[prevAtomIndex].species_string]) and not (
                        index in layerIndices):
                    layerIndices.append(index)
        if len(newAtomList) > 0:
            previousAtomList = newAtomList
            layerIndices = layerIndices + newAtomList
        else:
            finish = True

    layerIndices.sort()

    return layerIndices


def find_symmop(struct, angle=0, translation=(0, 0, 0), center=None, is_frac=True):

    # SYMMETRY PARAMETERS
    # ====================================================
    # angle=45
    # translation=(0,0,0)

    a = struct.lattice.get_cartesian_coords([0, 1, 0])
    b = struct.lattice.get_cartesian_coords([1, 0, 0])
    cStar = np.cross(a, b)

    if center is None:
        origin = struct.lattice.get_cartesian_coords([0.5, 0.5, 0.5])
    else:
        origin = center

    axis = cStar
    trans = tuple(struct.lattice.get_cartesian_coords(list(translation))) \
        if is_frac else translation

    if angle != 0:
        symOperation = SymmOp.from_origin_axis_angle(origin,
                                                     axis,
                                                     angle,
                                                     angle_in_radians=False)
    if angle == 0:
        symOperation = SymmOp.from_axis_angle_and_translation(
            axis, angle, angle_in_radians=False, translation_vec=trans)

    return(symOperation)


def find_translated_sites(struct, layerIndices, symOperation):

    # STRUCTURE MODIFICATION
    # ================================================
    # removal of the specified atoms
    remainingIndices = list(set(range(len(struct.sites))) - set(layerIndices))

    structLayerOnly = struct.copy()
    structLayerOnly.remove_sites(remainingIndices)
    structLayerOnly.apply_operation(symOperation)
    structLayerOnly.make_supercell([2, 2, 1])

    lattice1 = struct.lattice
    lattice2 = structLayerOnly.lattice
    [X, Y, Z] = lattice2.matrix
    [A1, B1, C1] = lattice1.abc
    C1 = 3 * C1

    trueSites = []  # {"specie": "Na" , "coordCart" : [X,Y,Z]}
    remainingSites = []

    # finding the right modulo !
    for currentSite in structLayerOnly.sites:
        siteDict = {
            "specie": currentSite.species_and_occu,
            "coordFrac2": currentSite.frac_coords,
            "coordCart": currentSite.coords}
        # print("\nFrac2 : {0} \n".format(siteDict["coordFrac2"]))
        coordFrac1 = lattice1.get_fractional_coords(siteDict["coordCart"])
        [a, b, c] = coordFrac1
        # print("Frac1 : {0} \n".format(coordFrac1))
        if a <= 1 and a > 0 and b <= 1 and b > 0:
            siteDict["coordFrac1"] = coordFrac1
            trueSites.append(siteDict)
        else:
            remainingSites.append(siteDict)

    print("Sites", len(structLayerOnly.sites))
    print("trueSites", len(trueSites))
    print("remainingSites", len(remainingSites))

    # remainingSites=structLayerOnly.sites

    i = 1
    while len(remainingSites) > 0 and i < 3:
        badSites = remainingSites
        remainingSites = []
        for siteDict in badSites:
            done = False
            #print("startLoop for {0}".format(siteDict["specie"]))
            coordCart = siteDict["coordCart"]
            for n in range(-i, +i, 1):
                for m in range(-i, +i, 1):
                    for l in range(-3, 3, 1):
                        if not done:
                            [Aslide, Bslide, Cslide] = [cart + n * x + m * y +
                                                        l * z for (cart, x, y, z) in zip(coordCart, X, Y, Z)]
                            # print("cart",Aslide,Bslide,Cslide)
                            coordFrac1 = lattice1.get_fractional_coords(
                                [Aslide, Bslide, Cslide])
                            [a, b, c] = coordFrac1
                            # print("frac1",a,b,c,"\n")
                            err = 0.000001

                            # print(n,m,l,i)
                            if a <= (1 + err * A1) and a > (0 - err * A1) and b <= (1 + err * B1) and b > (
                                    0 - err * B1) and c <= (0.5 + err * C1) and c > (-0.5 - err * C1):
                                print(
                                    "n:{0} A: {1:.3f} || m:{2} B:{3:.3f} ||  l:{4} C:{5:.3f}|| i: {6}\n".format(
                                        n, Aslide, m, Bslide, l, Cslide, i))
                                print("Frac1 : {0} \n \n".format(coordFrac1))
                                print("OK !!")
                                siteDict["coordFrac1"] = coordFrac1
                                trueSites.append(siteDict)
                                done = True
            if not done:
                #print("did not converge and i={0}".format(i))
                remainingSites.append(siteDict)
        i += 1

    return(trueSites)


def replace_shifted_sites(struct, layerIndices, trueSites):
    structWithoutLayer = struct.copy()
    structWithoutLayer.remove_sites(layerIndices)

    for siteDict in trueSites:
        try:
            structWithoutLayer.append(
                siteDict["specie"],
                siteDict["coordFrac1"],
                coords_are_cartesian=False,
                validate_proximity=True)
        except ValueError:
            print("# this site is to close !")

    # structWithoutLayer.to(fmt="poscar",filename="layer_POSCAR")

    return(structWithoutLayer)


def shift_indices_2(struct, layerIndices, angle, translation):
    symOperation = find_symmop(struct, angle, translation, center=None)
    trueSites = find_translated_sites(struct, layerIndices, symOperation)
    shifted_struct = replace_shifted_sites(struct, layerIndices, trueSites)

    return(shifted_struct)


def get_layer_indices_by_number(struct, specie_list, layer_index=-1):
    i = 0
    past_indices = []
    indices_to_exclude = []
    for specie in specie_list:
        indices_to_exclude += struct.indices_from_symbol(specie)
    remainingIndices = set(range(len(struct.sites))) - set(indices_to_exclude)
    print("layer index :", layer_index)

    while len(remainingIndices) > 0 and i <= layer_index:
        atom_height = [struct[n].c for n in remainingIndices]
        print(["i:{} h:{:.3f}".format(n, h)
               for h, n in sorted(zip(atom_height, remainingIndices))][0::5])
        remainingIndices_sorted = [x for _, x in sorted(
            zip(atom_height, remainingIndices))]
        n = list(remainingIndices_sorted)[0]
        print("index to start with : {}".format(n))
        current_layer = get_layer_indices(struct, specie_list, refIndex=n)
        remainingIndices = remainingIndices - set(current_layer)
        #print("i : {} , n : {}\ncurrent layer : {} \nRemaining indices :  {}".format(i,n,current_layer, remainingIndices))
        i += 1

    if i <= layer_index:
        print("this layer does not exist, returning last layer")

    return(current_layer)


def p2_to_O2_old(P2_struct):

    all_layer = get_layer_indices(P2_struct, ["Na"])
    # O_global=P2_struct.indices_from_symbol("O")
    #O_layerlist = (set(O_global) & set(all_layer))
    O2_struct = shift_indices(P2_struct, all_layer,
                              angle=0, translation=(1 / 3, 1 / 3, 0))
    O2_struct.to(fmt="poscar", filename="O2_POSCAR")
    return(O2_struct)


def p2_to_O2(P2_struct):

    MgO_layer = get_layer_indices_by_number(P2_struct, ["Na"], 1)
    for i in MgO_layer:
        print("{} {}".format(i, P2_struct[i].species_string))

    O_global = P2_struct.indices_from_symbol("O")
    Na_layer = list(set(get_layer_indices_by_number(
        P2_struct, ["Mg", "Mn"], 2)) - set(O_global))

    ind_to_move = Na_layer + MgO_layer

    for i in Na_layer:
        print("{} {}".format(i, P2_struct[i].species_string))

    symOperation = find_symmop(P2_struct, angle=0, translation=(1 / 6, 0, 0))
    trueSites = find_translated_sites(P2_struct, MgO_layer, symOperation)
    trueSites += find_translated_sites(P2_struct, Na_layer, symOperation)

    O2_struct = replace_shifted_sites(P2_struct, ind_to_move, trueSites)

    O2_struct.to(fmt="poscar", filename="O2_POSCAR")
    return(O2_struct)


def P2_to_OP4(P2_struct):
    # doubler la maille p2 en hauteur
    large_P2 = P2_struct.copy()
    large_P2.make_supercell([1, 1, 2])

    # selectionner deux layers de MGO consécutifs
    MgO_layer = get_layer_indices_by_number(large_P2, ["Na"], 1)
    MgO_layer += get_layer_indices_by_number(large_P2, ["Na"], 2)
    for i in MgO_layer:
        print("{} {}".format(i, large_P2[i].species_string))

    # selectionner deux layers de Na consécutifs
    O_global = large_P2.indices_from_symbol("O")
    Na_layer = list(set(get_layer_indices_by_number(
        large_P2, ["Mg", "Mn"], 2)) - set(O_global))
    Na_layer += list(set(get_layer_indices_by_number(large_P2,
                                                     ["Mg", "Mn"], 1)) - set(O_global))

    Na_layer_octa = list(set(get_layer_indices_by_number(
        large_P2, ["Mg", "Mn"], 4)) - set(O_global))

    for i in Na_layer:
        print("{} {}".format(i, large_P2[i].species_string))

    symOperation = find_symmop(large_P2, angle=0, translation=(1 / 6, 0, 0))
    trueSites = find_translated_sites(large_P2, MgO_layer, symOperation)
    trueSites += find_translated_sites(large_P2, Na_layer, symOperation)

    # should be enhanced to keep symmetry but nevermind
    symOperation_octa = find_symmop(
        large_P2, angle=0, translation=(
            1 / 6, 0, 0))
    trueSites += find_translated_sites(large_P2,
                                       Na_layer_octa, symOperation_octa)

    ind_to_move = MgO_layer + Na_layer + Na_layer_octa
    O2_struct = replace_shifted_sites(large_P2, ind_to_move, trueSites)

    O2_struct.to(fmt="poscar", filename="OP4_POSCAR")
    return(O2_struct)


def O3_to_P3(O3_struct, layers_3=True):
    large_O3 = O3_struct.copy()
    if layers_3:
        large_O3.make_supercell([1, 1, 3])
    test_folder = "/home/jvergnet/Téléchargements/structures/test_struct/O3_to_P3"
    size = "large" if layers_3 else "small"
    large_O3.to(
        fmt="poscar", filename="{}/O3_ {}_POSCAR".format(test_folder, size))

    print('ok')
    Na_global = large_O3.indices_from_symbol("Na")

    Me_layer_1 = get_layer_indices_by_number(large_O3, ["Na"], 1)
    layer_1 = Me_layer_1
    layer_2 = []

    if layers_3:
        Na_O_layer_1 = get_layer_indices_by_number(
            large_O3, ["Ni", "Ti", "Cu", "Mn"], 1)
        Na_layer_1 = list(set(Na_global) & set(Na_O_layer_1))
        layer_1 += Na_layer_1

        Me_layer_2 = get_layer_indices_by_number(large_O3, ["Na"], 2)
        Na_O_layer_2 = get_layer_indices_by_number(
            large_O3, ["Ni", "Ti", "Cu", "Mn"], 2)
        Na_layer_2 = list(set(Na_global) & set(Na_O_layer_2))
        layer_2 = Me_layer_2 + Na_layer_2

        # Na_O_layer_3 = get_layer_indices_by_number (large_O3 , ["Ni","Ti","Cu","Mn"], 2)
        # Na_layer_3 = list(set(Na_global) & set(Na_O_layer_3))

        print("Me_layer_1:{} \n + Na_layer_1{}\n".format(Me_layer_1, Na_layer_1))
        print("Me_layer_2:{} \n + Na_layer_2{}".format(Me_layer_2, Na_layer_2))

    trans_a = 0.141   # 0.1666
    trans_b = 0.065   # 0.1666
    print(trans_a, trans_b)

    symOperation = find_symmop(
        large_O3, angle=0, translation=(trans_a, trans_b, 0))
    trueSites = find_translated_sites(large_O3, layer_1, symOperation)

    if layers_3:
        symOperation_2 = find_symmop(
            large_O3, angle=0, translation=(2 * trans_a, 2 * trans_b, 0))
        trueSites += find_translated_sites(large_O3,
                                           Me_layer_2, symOperation_2)

        symOperation_3 = find_symmop(
            large_O3, angle=0, translation=(
                2 * trans_a, 2 * trans_b + 0.02, 0))
        trueSites += find_translated_sites(large_O3,
                                           Na_layer_2, symOperation_2)

    P3_struct = replace_shifted_sites(large_O3, layer_1 + layer_2, trueSites)

    P3_struct.to(fmt="poscar", filename="P3_new_POSCAR")

    return(P3_struct)


def shifted_sites_after_prismatic_disto(struct, metal_index, angle):
    mini = MinimumDistanceNN(tol=0.3)
    print("rotating around index : ", metal_index)
    neighbor_O = mini.get_nn_info(struct, metal_index)

    [X, Y, Z] = struct[metal_index].coords
    print("metal coords : ", [X, Y, Z])

    trueSites = []  # list of dict with index and shifted coordinates
    O_above = [O['site_index'] for O in neighbor_O if O['site'].z - Z > 0]
    O_below = [O['site_index'] for O in neighbor_O if O['site'].z - Z < 0]

    trueSites += find_translated_sites(struct, O_above,
                                       angle, translation=(0, 0, 0), center=[X, Y, Z])
    trueSites += find_translated_sites(struct, O_below, -
                                       angle, translation=(0, 0, 0), center=[X, Y, Z])

    return(trueSites)


def O3_to_O1(O3_struct, layers_3=False):
    # translate the layers Na(2) and TM
    large_O3 = O3_struct.copy()
    # large_O3.make_supercell([1, 1, 2])
    test_folder = "/home/jvergnet/Téléchargements/structures/test_struct/O3_to_O1"
    # size = "small"
    # large_O3.to(
    #     fmt="poscar", filename="{}/O3_ {}_POSCAR".format(test_folder, size))

    print('ok')
    Na_global = large_O3.indices_from_symbol("Na")

    Me_layer_1 = get_layer_indices_by_number(large_O3, ["Na"], 1)
    layer_1 = Me_layer_1
    layer_2 = []

    if layers_3:
        Na_O_layer_1 = get_layer_indices_by_number(
            large_O3, ["Ni", "Ti", "Cu", "Mn"], 1)
        Na_layer_1 = list(set(Na_global) & set(Na_O_layer_1))
        layer_1 += Na_layer_1

        Me_layer_2 = get_layer_indices_by_number(large_O3, ["Na"], 2)
        Na_O_layer_2 = get_layer_indices_by_number(
            large_O3, ["Ni", "Ti", "Cu", "Mn"], 2)
        Na_layer_2 = list(set(Na_global) & set(Na_O_layer_2))
        layer_2 = Me_layer_2 + Na_layer_2

        # Na_O_layer_3 = get_layer_indices_by_number (large_O3 , ["Ni","Ti","Cu","Mn"], 2)
        # Na_layer_3 = list(set(Na_global) & set(Na_O_layer_3))

        print("Me_layer_1:{} \n + Na_layer_1{}\n".format(Me_layer_1, Na_layer_1))
        print("Me_layer_2:{} \n + Na_layer_2{}".format(Me_layer_2, Na_layer_2))

    trans_a = 0.141   # 0.1666
    trans_b = 0.065   # 0.1666
    print(trans_a, trans_b)

    symOperation = find_symmop(
        large_O3, angle=0, translation=(trans_a, trans_b, 0))
    trueSites = find_translated_sites(large_O3, layer_1, symOperation)

    if layers_3:
        symOperation_2 = find_symmop(
            large_O3, angle=0, translation=(2 * trans_a, 2 * trans_b, 0))
        trueSites += find_translated_sites(large_O3,
                                           Me_layer_2, symOperation_2)

        symOperation_3 = find_symmop(
            large_O3, angle=0, translation=(
                2 * trans_a, 2 * trans_b + 0.02, 0))
        trueSites += find_translated_sites(large_O3,
                                           Na_layer_2, symOperation_2)

    P3_struct = replace_shifted_sites(large_O3, layer_1 + layer_2, trueSites)

    P3_struct.to(fmt="poscar", filename="P3_new_POSCAR")

    return(P3_struct)

#    shifted_struct = replace_shifted_sites(struct,O_above+O_below,trueSites)


# def shift_indices(struct,layerIndices, angle,translation) :

#     structWithoutLayer=struct.copy()
#     structWithoutLayer.remove_sites(layerIndices)

#     #rotation of one layer around the axis perpendicular to a and b (ie c*)
#     #that passes through the middle of the cell
#     #and translation in the ab plane
#     #input : structure, list of indices,
#     # angle in degrees, translation vector in fractional coordinates
#     #output : structure with rotated/translated layer

# # SYMMETRY PARAMETERS ====================================================
#     # angle=45
#     # translation=(0,0,0)
#     a=struct.lattice.get_cartesian_coords([0,1,0])
#     b=struct.lattice.get_cartesian_coords([1,0,0])
#     cStar=np.cross(a,b)
#     origin=struct.lattice.get_cartesian_coords([0.5,0.5,0.5])
#     axis=cStar
#     trans=tuple(struct.lattice.get_cartesian_coords(list(translation)))
#     symRot=SymmOp.from_origin_axis_angle(origin,
#                                          axis,
#                                          angle,
#                                          angle_in_radians=False)

#     symTrans =SymmOp.from_axis_angle_and_translation(axis,
#                                                    angle,
#                                                    angle_in_radians=False,
#                                                     translation_vec=trans)

#     # STRUCTURE MODIFICATION  ==========================================
#     # removal of the specified atoms
#     remainingIndices=list(set(range(0,len(struct.sites),1))-set(layerIndices))


#     structLayerOnly=struct.copy()
#     structLayerOnly.remove_sites(remainingIndices)
#     structLayerOnly.apply_operation(symTrans)
#     structLayerOnly.make_supercell([2,2,1])


#     lattice1=struct.lattice
#     lattice2=structLayerOnly.lattice
#     [X,Y,Z]=lattice2.matrix
#     [A1,B1,C1]=lattice1.abc
#     C1=3*C1


#     trueSites=[] # {"specie": "Na" , "coordCart" : [X,Y,Z]}
#     remainingSites=[]

#     # finding the right modulo !
#     for currentSite in structLayerOnly.sites :
#         siteDict={"specie": currentSite.species_and_occu , "coordFrac2" : currentSite.frac_coords , "coordCart" : currentSite.coords }
#         # print("\nFrac2 : {0} \n".format(siteDict["coordFrac2"]))
#         coordFrac1=lattice1.get_fractional_coords(siteDict["coordCart"])
#         [a,b,c]=coordFrac1
#         # print("Frac1 : {0} \n".format(coordFrac1))
#         if a<= 1 and a>0 and b <= 1 and b>0 :
#             siteDict["coordFrac1"]=coordFrac1
#             trueSites.append(siteDict)
#         else :
#             remainingSites.append(siteDict)

#     # print("Sites",len(structLayerOnly.sites))
#     # print("trueSites",len(trueSites))
#     # print("remainingSites",len(remainingSites))

#     # remainingSites=structLayerOnly.sites

#     i=1
#     while len( remainingSites )>0 and i < 3 :
#         badSites=remainingSites
#         remainingSites = []
#         for siteDict in badSites :
#             done=False
#             #print("startLoop for {0}".format(siteDict["specie"]))
#             coordCart=siteDict["coordCart"]
#             for n in range(-i,+i,1):
#                 for m in range(-i,+i,1) :
#                     for l in range(-3,3,1) :
#                         if done == False :
#                             [Aslide,Bslide,Cslide]=[cart+n*x+m*y+l*z for (cart,x,y,z) in zip(coordCart,X,Y,Z)]
#                             # print("cart",Aslide,Bslide,Cslide)
#                             coordFrac1=lattice1.get_fractional_coords([Aslide,Bslide,Cslide])
#                             [a,b,c]=coordFrac1
#                             # print("frac1",a,b,c,"\n")
#                             err=0.000001

#                             # print(n,m,l,i)
#                             if a<=(1+err*A1) and a>(0-err*A1) and b <= (1+err*B1) and b>(0-err*B1) and c<=(0.5+err*C1) and c>(-0.5-err*C1):
#                                 print("n:{0} A: {1:.3f} || m:{2} B:{3:.3f} ||  l:{4} C:{5:.3f}|| i: {6}\n".format(n,Aslide,m,Bslide,l,Cslide,i))
#                                 print("Frac1 : {0} \n \n".format(coordFrac1))
#                                 print("OK !!")
#                                 siteDict["coordFrac1"]=coordFrac1
#                                 trueSites.append(siteDict)
#                                 done=True
#             if done==False :
#                 #print("did not converge and i={0}".format(i))
#                 remainingSites.append(siteDict)
#         i+=1

#     for siteDict in trueSites :
#         try:
#             structWithoutLayer.append(siteDict["specie"],siteDict["coordFrac1"], coords_are_cartesian=False, validate_proximity=True)
#         except ValueError :
#             print("# this site is to close !")

#     #structWithoutLayer.to(fmt="poscar",filename="layer_POSCAR")

#     return structWithoutLayer

# def old_rotate_translate_layer(struct,angle=0,translation=(0,0,0)) :

#     #rotation of one layer around the axis perpendicular to a and b (ie c*)
#     #that passes through the middle of the cell
#     #and translation in the ab plane
#     #input : structure, angle in degrees, translation vector in fractional coordinates
#     #output : structure with rotated/translated layer

#     # SYMMETRY PARAMETERS ==============================================
#     a=struct.lattice.get_cartesian_coords([0,1,0])
#     b=struct.lattice.get_cartesian_coords([1,0,0])
#     cStar=np.cross(a,b)

#     origin=struct.lattice.get_cartesian_coords([0.5,0.5,0.5])
#     axis=cStar

#     symRot=SymmOp.from_axis_angle_and_translation(axis,
#                                                    angle,
#                                                    angle_in_radians=False,
#                                                    translation_vec=translation)


#     # Structural modification ==========================================
#     # removal of the layer atoms (metal + oxygens)
#     layerIndices=cluster.get_layer_indices(struct)
#     structWithout=struct.copy()
#     structWithout.remove_sites(layerIndices)
#     structWithout.to(fmt="poscar",filename="without_POSCAR")

#     # insertion of the rotated atoms with the same index
#     for siteIndex in layerIndices :
#         initialSite=struct[siteIndex]
#         structWithout.insert(siteIndex,initialSite.species_and_occu,symRot.operate(initialSite.coords), coords_are_cartesian=True, validate_proximity=False, properties=None)

#     #structWithout.to(fmt="poscar",filename="with_POSCAR")
#     return structWithout

# def old_rotate_layer(struct,angle=0,translation=(0,0,0)) :

#     #rotation of one layer around the axis perpendicular to a and b (ie c*)
#     #that passes through the middle of the cell
#     #and translation in the ab plane
#     #input : structure, angle in degrees, translation vector in fractional coordinates
#     #output : structure with rotated/translated layer

# # SYMMETRY PARAMETERS ====================================================
#     # angle=15
#     # translation=(0,0,0)
#     a=struct.lattice.get_cartesian_coords([0,1,0])
#     b=struct.lattice.get_cartesian_coords([1,0,0])
#     cStar=np.cross(a,b)

#     origin=struct.lattice.get_cartesian_coords([0.5,0.5,0])
#     axis=cStar
#     symRot=SymmOp.from_origin_axis_angle(origin,
#                                   axis,
#                                   angle,
#                                   angle_in_radians=False)
#     # symTrans =SymmOp.from_axis_angle_and_translation(axis,
#     #                                                angle=0,
#     #                                                angle_in_radians=False,
#     #                                                translation_vec=translation)

#     layerIndices=get_layer_indices(struct)
#     translateStruct=struct.copy()
#     for i in layerIndices :
#         initialSite=struct[i]
#         translaVectos=symRot.operate(initialSite.coords)-initialSite.coords
# #        translaVectos=struct.lattice.get_fractional_coords(translaVectos)

#         translateStruct.translate_sites(indices=[i],vector=translaVectos,frac_coords=False, to_unit_cell=True)
#     return translateStruct
