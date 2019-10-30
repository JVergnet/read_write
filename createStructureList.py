# createStructureList.py
from copy import deepcopy

import numpy as np
import pymatgen.analysis.magnetism.analyzer as mag_anal
from pymatgen.io.vasp.sets import MITRelaxSet, MPStaticSet
from pymatgen.transformations import site_transformations as st
from pymatgen.transformations import standard_transformations as sd
from pymatgen.transformations.advanced_transformations import \
    EnumerateStructureTransformation

import cluster
import launchDisordered as launch

# METAL DISORDER
# ----------------------------------------


def disorder_in_supercell(pristineStruct,
                          substSpecies={},
                          supercell_size=[1, 1, 1],
                          number_of_struct=10):

    ini = pristineStruct.copy()

    # substitued structures =====================================
    # substitued_1 = sd.SubstitutionTransformation(
    #     {"Na":  {"Li": 0.5, "Fe": 0.5}}).apply_transformation(ini)
    substitued_2 = sd.SubstitutionTransformation(
        {"Ti":  {"Ti": 7/8, "Fe": 1/8}}).apply_transformation(ini)
    ini = substitued_2
    # GENERATE OXIDATION STATES ================================

    try:
        oxi = sd.AutoOxiStateDecorationTransformation().apply_transformation(ini)
    except ValueError:
        print("Forced to set oxidation by hand")
        oxi.add_oxidation_state_by_element({"Na": 1,
                                            "Mg": 2,
                                            "Mn": 3.5,
                                            "O": -2,
                                            "Li": 1,
                                            "Ca": 2,
                                            "Fe": 3,
                                            "Zn": 4,
                                            "Co": 2.5,
                                            "Ti": 4,
                                            "S": -2,
                                            'P': 5})
    print(oxi)
    # transform into 2,2,2 supercell
    [x, y, z] = supercell_size
    megaOxi = sd.SupercellTransformation.from_scaling_factors(
        x, y, z).apply_transformation(oxi)
    # dictPri = {'structure': megaOxi, 'id': "pristine"}
    print(megaOxi)
    # chosenList.append(dictPri) #uncomment to compute pistine !

    #
    # generation of a list of ordered structures from the disordered one
    chosenList = []
    try:
        type_of_algo = input(
            "type of algo : [E]numlib (more efficient but weaker) // [P]ymatgen built-in (slower but safer)")[0]
    except BaseException:
        type_of_algo = "E"
        print("trying E then P sequence")
    try:
        if type_of_algo == "E":
            enum = EnumerateStructureTransformation(
                min_cell_size=1,
                max_cell_size=None,
                symm_prec=0.1,
                refine_structure=False,
                enum_precision_parameter=0.1,
                check_ordered_symmetry=True,
                max_disordered_sites=100,
                sort_criteria='ewald',
                timeout=None)

            enumList = enum.apply_transformation(
                megaOxi, return_ranked_list=number_of_struct)
    except Exception as ex:
        print(
            "error during enumlib method : {} \n fall back to [P]ymatgen buit-in method",
            ex)
        type_of_algo = "P"
    try:
        if type_of_algo == "P":
            enumList = sd.OrderDisorderedStructureTransformation().apply_transformation(
                megaOxi, return_ranked_list=number_of_struct)
    except Exception as ex:
        print(
            "error during [P]ymatgen buit-in method {} \n You're on your own beyond this point...",
            ex)

    for i, chosenStruct in enumerate(enumList):
        chosenList.append(
            {'structure': chosenStruct['structure'], "id": str(i + 1)})
    print(chosenList)

    return chosenList

# INTERCALANT DISORDER
# ---------------------------------------------------


def sodium_disorder(
        pristine_job,
        number_of_struct=10,
        fracRemoved=-1,
        nb_removed=-1,
        alkali="Na"):

    oxi = pristine_job.structure.copy()
    try:
        oxi = sd.AutoOxiStateDecorationTransformation().apply_transformation(oxi)
    except Exception as ex:
        print("Forced to set oxidation by hand")
        not_oxidizeable = {
            "Na": 1,
            "Li": 1,
            "Mg": 2,
            "Cu": 2,
            "Ca": 2,
            "Zn": 2,
            "Ti": 4}
        oxidizeable = {"Ni": 2, "Co": 2, "Fe": 3, "Mn": 4}
        anions = {"O": -2, "S": -2}
        actual_oxidation_states = not_oxidizeable
        actual_oxidation_states.update(oxidizeable)
        actual_oxidation_states.update(anions)
        print(actual_oxidation_states)
        oxi.add_oxidation_state_by_element(actual_oxidation_states)

    print("alkali :", alkali)
    try:
        indices = pristine_job.structure.indices_from_symbol(alkali)
    except Exception as ex:
        print("error : {} \n Probably no {} in the structure".format(ex, alkali))
        return [pristine_job.copy_w_new_struct(oxi, new_id="_0")]

    print("nb_removed{} , fracRemoved {}".format(nb_removed, fracRemoved))
    if nb_removed == -1 and fracRemoved >= 0:
        fraction_removed = len(indices) * fracRemoved
        nb_removed = int(round(fraction_removed, 0))
        print("removing {} e- ({:.3})".format(nb_removed, fraction_removed))
    elif fracRemoved == -1 and nb_removed == -1:
        print("Number of deintercalated atom not defined ")
        exit(1)

    enum = st.PartialRemoveSitesTransformation(
        [list(indices)], [[fracRemoved]], algo=3)
    # print(" frac removed ({0}) x indices ({1}) = nb removed ({2})".format(fracRemoved, len(indices), nb_removed))
    enumList = enum._fast_ordering(oxi,
                                  num_remove_dict={indices: nb_removed},
                                  num_to_return=number_of_struct)

    chosenList = [pristine_job.copy_w_new_struct(
        chosenStruct['structure'], new_id="_{}".format(i))
        for i, chosenStruct in enumerate(enumList)]
    return chosenList


def remove_doubles(tmpList, change_id=False):
    # comparing all the tmp structure found to the already found structures
    k = 0
    nextList = []
    # print("tmp list lenght : {0}".format(len(tmpList)))
    for tmpJob in tmpList:
        new = True
        for nextStruct in nextList:
            if tmpJob.structure.matches(nextStruct.structure,
                                        # ignored_species=["O","Mn","Mg"],
                                        scale=False):
                # print("Structure : {0} match with previous one : {1}"
                #       .format(tmpJob['id'],nextStruct['id']) )
                new = False
        if new:
            if change_id:
                tmpJob.entry_id += "-{0}".format(k)
            # print("new struc {0}".format(tmpJob.entry_id))
            nextList.append(tmpJob)
            k += 1
    print("TMP {0} => NEXT {1}\n".format(len(tmpList), len(nextList)))
    return(nextList)


def desodiation_tree(
        pristine_job,
        xmin,
        xmax,
        nb_compo=2,
        nb_struct_per_x=1,
        alkali="Na"):

    # xmin_rem = 1 - xmin
    # xmax_rem = 1 - xmax

    pristineStruct = pristine_job.structure
    abs_frac_left = np.append([pristine_job.x_na],
                              np.linspace(xmax, xmin,
                                          nb_compo,
                                          endpoint=True))

    rel_frac_left = [(abs_frac_left[i + 1] / abs_frac_left[i])
                     for i in range(0, len(abs_frac_left) - 1, 1)]
    rel_frac_rem = [1 - r for r in rel_frac_left]

    print("absolute fraction left {} \n".format(abs_frac_left))
    print("Relative fraction left {} \n".format(rel_frac_left))
    # print("absolute fraction removed {} \n".format(rel_frac_rem))
    print("Relative fraction Removed {} \n".format(rel_frac_rem))

    prevList = [pristine_job.copy_w_new_struct(
        pristineStruct, new_id="-P")]
    chosenList = prevList

    for i, x in enumerate(rel_frac_rem):
        print("i : {0} for structure :{1} "
              .format(i, prevList[0].structure.composition.formula))
        tmpList = []
        for dictStruct in prevList:
            tmpList += sodium_disorder(dictStruct,
                                       number_of_struct=nb_struct_per_x,
                                       fracRemoved=x, alkali=alkali)

        prevList = remove_doubles(tmpList, change_id=True)
        chosenList += prevList

    return chosenList


def desodiation_incremental(pristineStruct, xmin=0.5, xmax=1):

    return([])


# OXYGEN DISORDER
# -----------------------------------

def OO_distortion(pristine,
                  bailar_disto_min=0, bailar_disto_max=0, bailar_nb_steps=0,
                  trigo_disto_min=0, trigo_disto_max=0, trigo_nb_steps=0):

    bailar_disto_list = [bailar_disto_max] if bailar_nb_steps == 0 else [
        bailar_disto_min +
        i *
        (
            (bailar_disto_max -
             bailar_disto_min) /
            bailar_nb_steps) for i in range(
            0,
            bailar_nb_steps +
            1)]
    trigo_disto_list = [trigo_disto_max] if trigo_nb_steps == 0 else [
        trigo_disto_min +
        i *
        (
            (trigo_disto_max -
             trigo_disto_min) /
            trigo_nb_steps) for i in range(
            0,
            trigo_nb_steps +
            1)]

    pairs = cluster.get_OO_pairs(pristine)
    print("nb pairs : ", len(pairs))
    print("pairs", pairs)
    distorted_pristine_list = []
    for bailar_disto in bailar_disto_list:
        for trigo_disto in trigo_disto_list:
            s_copy = pristine.copy()
            if trigo_disto != 0 or bailar_disto != 0:
                for pair in pairs:
                    s_copy = cluster.contract_OO_pairs(
                        s_copy, pair, bailar_disto=bailar_disto, trigonal_disto=trigo_disto)
            # s_copy.to(fmt="POSCAR", filename="test_{}_POSCAR".format(i))
            # C = CifWriter(p, symprec=0.5, write_magmoms=False)
            # C.write_file("disto_{}_P2_full.cif".format(n))
            distorted_pristine_list.append(
                {"structure": s_copy, "id": "_Dbailar_{:.4f}_Dtrigo_{:.4f}".format(bailar_disto, trigo_disto)})
            print(
                "Disto OO : bailar : {:.4f} + trigo : {:.4f}".format(bailar_disto, trigo_disto))

    return(distorted_pristine_list)


def trigonal_distortion(pristine, disto_max, nb_steps):
    pairs = cluster.get_OO_pairs(pristine)
    print("nb pairs : ", len(pairs))
    print("pairs", pairs)
    distorted_pristine_list = []
    for i in range(0, nb_steps + 1):
        s_copy = pristine.copy()
        n = i * (disto_max / nb_steps)
        for pair in pairs:

            s_copy = cluster.contract_OO_pairs(
                s_copy, pair, bailar_disto=False, trigonal_disto=True)
        distorted_pristine_list.append(
            {"structure": s_copy, "id": "_Dtrigo_{:.4f}".format(n)})
        print("Dtrigo_{}".format(n))
    return(distorted_pristine_list)




def rotation_in_supercell(initialStruct,
                          supercell_size=[2, 2, 1],
                          angle_max=180,
                          angle_step=15):

    [x, y, z] = supercell_size
    superCell = sd.SupercellTransformation.from_scaling_factors(
        x, y, z).apply_transformation(initialStruct)  #
    # uncomment to compute pistine !

    current_angle = 0
    currentStruct = superCell

    chosenList = []
    while True:
        print(currentStruct)
        chosenList.append({'structure': currentStruct,
                           'id': "{0}deg".format(current_angle)})
        current_angle += angle_step
        currentStruct = cluster.rotate_oxygens_in_metal_layer(
            superCell, angle=current_angle, translation=(0, 0, 0))
        if current_angle > angle_max:
            break

    return chosenList

 # (currentStructure,settingCopy, 0, 6, 1)


def scan_nupdown(job, step):
    mag = mag_anal.CollinearMagneticStructureAnalyzer(
        job.structure,
        overwrite_magmom_mode='replace_all',
        round_magmoms=False,
        detect_valences=False,
        make_primitive=False,
        default_magmoms={
            "Mn": 3},
        threshold=0.1)

    mag_list = mag.magmoms
    job_list = []
    for N in range(0, step + 1):
        # varying from 2 to 4 = 1+1.5*(N/step)
        mag_list_div = (mag_list / 3) * (2 + 2 * (N / step))
        NUPDOWN = np.sum(mag_list_div)
        new_job = job.copy()
        new_job.structure.add_site_property('magmom', mag_list_div)
        new_job.incar['SYSTEM'] += "_SZ_{:.2f}".format(
            NUPDOWN / len(new_job.structure.indices_from_symbol("Mn")))
        new_job.incar["NUPDOWN"] = NUPDOWN
        print("{}\n\n".format(new_job.incar["SYSTEM"]))
        job_list.append(new_job)

    return job_list


def scan_U(structure, incarSetings, Umin, Umax, step):
    vaspInputSetList = []
    for U in range(Umin, Umax + 1, step):
        settingCopy = dict(incarSetings)
        print(U)
        # print(settingCopy['LDAUU'])
        # settingCopy['LDAUU']={'Mn':0}
        settingCopy['LDAUU'] = {'Mn': U}
        settingCopy['SYSTEM'] += "U{0}".format(U)
        inputSet = MITRelaxSet(
            structure,
            user_kpoints_settings={
                'reciprocal_density': 50},
            force_gamma=True,
            user_incar_settings=settingCopy)
        vaspInputSetList.append(inputSet)

    return vaspInputSetList


{
    'LDAUJ': {
        'F': {
            'Ag': 0,
            'Co': 0,
            'Cr': 0,
            'Cu': 0,
            'Fe': 0,
            'Mn': 0,
            'Mo': 0,
            'Nb': 0,
            'Ni': 0,
            'Re': 0,
            'Ta': 0,
            'V': 0,
            'W': 0},
        'O': {
            'Ag': 0,
            'Co': 0,
            'Cr': 0,
            'Cu': 0,
            'Fe': 0,
            'Mn': 0,
            'Mo': 0,
            'Nb': 0,
            'Ni': 0,
            'Re': 0,
            'Ta': 0,
            'V': 0,
            'W': 0},
        'S': {
            'Fe': 0,
            'Mn': 0}},
    'LDAUL': {
        'F': {
            'Ag': 2,
            'Co': 2,
            'Cr': 2,
            'Cu': 2,
            'Fe': 2,
            'Mn': 2,
            'Mo': 2,
            'Nb': 2,
            'Ni': 2,
            'Re': 2,
            'Ta': 2,
            'V': 2,
            'W': 2},
        'O': {
            'Ag': 2,
            'Co': 2,
            'Cr': 2,
            'Cu': 2,
            'Fe': 2,
            'Mn': 2,
            'Mo': 2,
            'Nb': 2,
            'Ni': 2,
            'Re': 2,
            'Ta': 2,
            'V': 2,
            'W': 2},
        'S': {
            'Fe': 2,
            'Mn': 2.5}},
    'LDAUTYPE': 2,
    'LDAUU': {
        'F': {
            'Ag': 1.5,
            'Co': 3.4,
            'Cr': 3.5,
            'Cu': 4,
            'Fe': 4.0,
            'Mn': 3.9,
            'Mo': 4.38,
            'Nb': 1.5,
            'Ni': 6,
            'Re': 2,
            'Ta': 2,
            'V': 3.1,
            'W': 4.0},
        'O': {
            'Ag': 1.5,
            'Co': 3.4,
            'Cr': 3.5,
            'Cu': 4,
            'Fe': 4.0,
            'Mn': 3.9,
            'Mo': 4.38,
            'Nb': 1.5,
            'Ni': 6,
            'Re': 2,
            'Ta': 2,
            'V': 3.1,
            'W': 4.0},
        'S': {
            'Fe': 1.9,
            'Mn': 2.5}}}


def replace(pristineStruct, substSpecies=["Li", "Fe", "Ca"],
            supercell_size=[1, 1, 1]):
    struct_dict_list = []
    for subst_specie in substSpecies:
        substitued_struct = sd.SubstitutionTransformation(
            {"Mg": subst_specie}).apply_transformation(pristineStruct)
        struct_dict_list.append({'structure': substitued_struct,
                                 'id': "{}".format(subst_specie)})
    return struct_dict_list


def selective_dynamic(vaspInputset, indices_to_move):
    sd_poscar = vaspInputset.poscar
    mysd = np.zeros([sum(sd_poscar.natoms), 3], bool)
    for i in indices_to_move:
        mysd[i][1] = True
        mysd[i][0] = True
        mysd[i][2] = True
    sd_poscar.selective_dynamics = mysd
    return (vaspInputset)


def bulk_strain(job, max_strain, steps):
    job_list = []

    # s = job.structure
    # str_id = job.entry_id
    # data = job.data
    for strain in np.linspace(-max_strain, max_strain,
                              2*int(steps)+1, endpoint=True):
        print("strain : {}".format(strain))
        new_job = job.copy()
        print(new_job)
        new_job.structure.apply_strain(strain)
        new_job.entry_id += "strain_{:.4f}".format(strain)
        new_job.data.update({"strain": strain})
        new_job.incar['ISIF'] = 2
        job_list.append(new_job)

    return(job_list)
