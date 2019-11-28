# readBader.py
"""
Analyze CHGCAR using Henkelman algo (if not already done)
Parse rundict.job_folder/BADER/*/acf.dat (* = CHARGE or SPIN))
Compare bader pop vs potcar pop to get oxidation state
Plot bader pop & oxidation state (site-wise or structure average)
"""
import os
import subprocess
from multiprocessing import Pool, cpu_count

import matplotlib.pyplot as plt
import numpy as np
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.io.vasp.inputs import Poscar, Potcar


# import readRun_entries as read
import run_utils.platform_id as platform_id
import run_utils.generic_plot as g_plot


def get_bader_tags(rundict_list):
    """
    fetch the bader values of each run via parallel threads
    then feed the data within rundict.structure_data
    for readable acess later on
    """
    # nb : poll map cannot handle (pickle) big objects
    # so we poll (small) arrays instead of (big) Rundict instances

    poll_bader_array = []
    try:
        with Pool(processes=cpu_count()) as parallel_threads:
            poll_bader_array = parallel_threads.map(
                get_bader_data_single, rundict_list)
            parallel_threads.close()
            parallel_threads.join()

    except Exception as ex:
        print("EXCEPTION", ex)

    print("nb of bader sucess : {} / {}".format(
        [poll[0] for poll in poll_bader_array].count(True),
        len(rundict_list)))

    for (rundict, (bader_success, bader_array)) in zip(rundict_list, poll_bader_array):
        if not bader_success:
            continue
        if rundict.bader_done:
            continue
        struct = rundict.structure_data
        for j, prop in enumerate(
                ['charge', 'vol_chg', 'magnetization', 'vol_mag']):
            struct.add_site_property(prop, [bader_array[i][j]
                                            for i in range(struct.num_sites)])
        rundict.bader_done = True

    print("bader tag after polling")
    print(rundict_list[0].structure_data[0].properties)

    return rundict_list


def get_bader_data_single(rundict):
    """
    Check existence of bader files / performs the calculation if no file is found
    read the data
    feed it to the "structure data" attribute of the corresponding rundict instance
    return :
        bool : sucess (True) or failure (False) of bader computation
        array : bader_data_array (only if necessary)
        bader_array[i] = (bader_charges,charge_vol,bader_magmom,magmom_vol)
    """

    if rundict.bader_done:
        print("bader already done for {}".format(rundict.name_tag))
        return(True, None)

    valid_bader_files = generate_bader_files(rundict.job_folder)

    if not valid_bader_files:
        print("Bader computation aborted for {}".format(rundict.name_tag))
        return(False, None)

    bader_data_array = read_bader_files(rundict.job_folder)
    return(True, bader_data_array)


def generate_bader_files(folder, remove_mag=True):
    """
    Run bader analysis on the vasprun directories
    only keep .dat files
    if the /BADER folder aleardy exists, we skip the computation
    Require the alias "do_bader" to be set in the ~/.bashrc of the user
    return True if the .dat files are created or found
    """
    sucessful_bader = False
    try:
        print("in do_bader for folder : \n {}".format(folder))

        if not os.path.exists("{}/BADER/CHARGE/ACF.dat".format(folder)) and \
           not os.path.exists("{}/BADER/SPIN/ACF.dat".format(folder)):
            print(folder, "BADER data not found !")
            if os.path.exists("{}/CHGCAR".format(folder)):
                print(" Starting bader in {} ".format(folder))
                os.chdir(folder)
                proc = subprocess.Popen(['do_bader'], shell=True)
                print("process launched")
                proc.wait()
            else:
                print(folder, "no CHGCAR found")

        if os.path.exists("{}/BADER/CHARGE/ACF.dat".format(folder)) and \
           os.path.exists("{}/BADER/SPIN/ACF.dat".format(folder)):
            sucessful_bader = True

            if remove_mag:
                print("removing useless bader files  : ")
                bader_folder = os.path.join(folder, "BADER")
                data_folder_list = [bader_folder] + \
                    [os.path.join(bader_folder, data)
                     for data in ["SPIN", "CHARGE"]]

                for data_folder in data_folder_list:
                    for bader_file in [os.path.join(data_folder, o)
                                       for o in os.listdir(data_folder)]:
                        if not os.path.isdir(bader_file) \
                                and bader_file.split('.')[-1] != "dat":
                            print("removing", bader_file.split("/")[-1])
                            os.remove(bader_file)

                print(
                    "sucessfully removed previous bader files (CHG,CHGCAR,AECCAR) \n\n")

    except Exception as ex:
        print("do_bader failed in {} : {} : ".format(folder, ex))

    return sucessful_bader


def read_bader_files(folder):
    """
    Read the bader (charge + magmom) calculation done in a folder
    (Performs the bader calculation if the file is not present
    """

    # bader_charges = []
    # charge_vol = []
    # bader_magmom = []
    # magmom_vol = []

    os.chdir(folder)
    if os.path.exists("{}/BADER/CHARGE/ACF.dat".format(folder)):
        print('bader charge : OK ')
        bad_chg = np.loadtxt("{}/BADER/CHARGE/ACF.dat".format(folder),
                             comments=['#', '---', 'VACUUM', 'NUMBER'],
                             usecols=[4, 6], )
        # with open('BADER/CHARGE/ACF.dat') as f_charge:
        #     lines = [line.rstrip() for line in f_charge]
        #     for l in lines[2: -4]:
        #         # print(l)
        #         bader_charges.append(float(l.split()[4]))
        #         charge_vol.append(float(l.split()[6]))
    if os.path.exists("{}/BADER/SPIN/ACF.dat".format(folder)):
        print('bader spin : OK ')
        bad_mag = np.loadtxt('{}/BADER/SPIN/ACF.dat'.format(folder),
                             comments=['#', '---', 'VACUUM', 'NUMBER'],
                             usecols=[4, 6], )
        # with open('BADER/SPIN/ACF.dat') as f_spin:
        #     lines = [line.rstrip() for line in f_spin]
        #     for l in lines[2: -4]:
        #         # print(l)
        #         bader_magmom.append(float(l.split()[4]))
        #         magmom_vol.append(float(l.split()[6]))
    # bader_array = [
    #     n for n in zip(
    #         bader_charges,
    #         charge_vol,
    #         bader_magmom,
    #         magmom_vol)]
    bader_array = np.hstack((bad_chg, bad_mag))
    return bader_array


def get_charge(self, atom_index):
    """
    Convenience method to get the charge on a particular atom.

    Args:
        atom_index:
            Index of atom.

    Returns:
        Charge associated with atom from the Bader analysis.
    """
    return self.data[atom_index]["charge"]


def get_charge_transfer(atom_index, structure, potcar, poscar):
    """
    Returns the charge transferred for a particular atom. Requires POTCAR
    to be supplied.

    Args:
        atom_index:
            Index of atom.

    Returns:
        Charge transfer associated with atom from the Bader analysis.
        Given by final charge on atom - nelectrons in POTCAR for
        associated atom.
    """
    if potcar is None:
        raise ValueError("POTCAR must be supplied in order to calculate "
                         "charge transfer!")
    potcar_indices = []
    # natoms = [1,3] for TiS3
    for i, nb_same_specie in enumerate(poscar.natoms):
        potcar_indices += [i] * nb_same_specie  # [1,2,2,2] for TiS3
        # s.composition.element_composition[site.specie]

    # atom 3 (S) : potcar[2]
    nelect = potcar[potcar_indices[atom_index]].nelectrons

    return structure[atom_index].properties["charge"] - nelect


def get_oxidation_state_decorated_structure(structure, potcar, poscar):
    """
    Returns an oxidation state decorated structure.

    Returns:
        Returns an oxidation state decorated structure. Requires POTCAR
        to be supplied.
    """
    # structure = structure.copy()
    charges = [-get_charge_transfer(i, structure, potcar, poscar)
               for i in range(len(structure))]
    structure.add_oxidation_state_by_site(charges)
    structure.add_site_property("bader_oxi", charges)
    return structure


def get_madelung_tag_single(rundict, force=False, verbose=0, oxi_int=False):
    struct_data = rundict.structure_data
    if not force and struct_data[0].properties.get("E_mad", None) is not None:
        print("skipping")
        return True
    poscar, potcar = [p.from_file(rundict.job_folder+f) for (p, f) in
                      [(Poscar, "/POSCAR"), (Potcar, "/POTCAR")]]
    s_dict = {}
    if verbose > 0:
        print("poscar/potcar parsed")
    if oxi_int:
        oxi = rundict.structure.copy()
        try:
            oxi.add_oxidation_state_by_guess()
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
        ew_sum_int = EwaldSummation(oxi)
        oxi.add_site_property("E_mad_int",
                              [ew_sum_int.get_site_energy(i)
                               for i in range(oxi.num_sites)])
        s_dict['oxi_int'] = oxi
        struct_data.add_site_property("E_mad_int",
                                      [ew_sum_int.get_site_energy(i)
                                       for i in range(oxi.num_sites)])
        if verbose > 0:
            print("integer oxidation state computed")

    if struct_data[0].properties.get("charge", None) is None:
        print("No bader charges => no bader madelung")
        return False
        # Test if bader charges were computed correctly

    get_oxidation_state_decorated_structure(struct_data,
                                            potcar, poscar)
    if verbose > 0:
        print("bader oxidation state computed")

    ew_sum = EwaldSummation(rundict.structure_data)
    # madelung_array = []
    struct_data.add_site_property("E_mad",
                                  [ew_sum.get_site_energy(i)
                                   for i in range(struct_data.num_sites)])
    s_dict["oxi_bader"] = struct_data
    print(" ".join(s_dict.values()))
    return True


# ===== PLOTTING ======
# =====================
def plot_charge_and_mag(rundict_list, detailled=None, **kwargs):
    """
    plot charge and magmom for each specie of a starting structure
    then save them with proper name
    """

    print("\n==== ELECTRONS IN REAL SPACE (BADER) ======\n")

    coord = kwargs["coord"] if "coord" in kwargs.keys() else "x_na"
    print(rundict_list[0].structure_data)
    if detailled is None:
        detailled = 0
        try:
            print("""
plotting charge and magmom
0 = total charge only
1 = site charges
2 = site chg and mag,
3 = chg & mag volumes 
( < 0 to pass)\n
            """)
            detailled = int(input("Level of verbosity ?   : "))
        except BaseException:
            print("default to 0")

    figures = {}
    elements = set()
    for run in rundict_list:
        elements.update(
            set(run.structure.composition.get_el_amt_dict().keys()))

    if detailled == 0:
        for specie in elements:
            fig_name = "total_charge_on_{}".format(specie)
            figures[fig_name] = plt.figure(fig_name)

            axe = figures[fig_name].add_subplot(1, 1, 1)
            g_plot.plot_site_value_evolution(
                rundict_list,
                specie,
                value="charge",
                coord=coord,
                plot_type=[4],
                axe0=axe)
            # plot_type = 0:sites 1:avg 2:min 3:max 4:sum/nbCell

    if detailled > 0:
        prop_list = ["charge"]
        if detailled > 1:
            prop_list.append('magnetization')
        if detailled > 2:
            prop_list += ['vol_chg', 'vol_mag']
        for specie in elements:
            for prop in prop_list:
                figures["{}_of_{}".format(prop, specie)] = g_plot.plot_site_value_evolution(
                    rundict_list, specie, value=prop, coord=coord, plot_type=[0, 1])
                # plot_type = 0:sites 1:avg 2:min 3:max 4:sum

    try:
        print('Type[s]ave to save all bader figures \n')
        if input(": ") == "s":
            prefix = input(" prefix ? : ")
            for f_name in figures.keys():
                current_fig_name = platform_id.get_file_name(
                    os.getcwd(), "{}_{}".format(prefix, f_name), ext=".svg")
                figures[f_name].savefig(
                    "{}.svg".format(current_fig_name),
                    bbox_inches='tight')
    except Exception:
        pass
    return True
