# filter_runs.py


import traceback


def generate_hull_tags(run_list, coord="x_na"):
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
    current_x = -1000000
    for entry in sorted_entries:
        if getattr(entry, coord) > current_x:
            clean_entries.append(entry)
            current_x = getattr(entry, coord)

    for entry in clean_entries:
        entry.status = 4

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

    for entry in hull_entries:
        entry.status = 5

    return(sorted_entries)


def restrict_run_list(all_runs_input, x_coord):
    """
    defines a restriction on the runs
    if the given list is not locked, performs filtering on the list
    """
    selection_loop = True
    while selection_loop:

        # STACKING FAMILY FILTER ================
        restricted_stack_runs = stacking_filter(all_runs_input)
        print("number of selected runs : {}".format(
            len(restricted_stack_runs)))

        # CONVEX HULL FILTER ==================
        sieve_lvl = select_sieve_level()

        restricted_hull_runs = hull_filtering(
            sieve_lvl, restricted_stack_runs, x_coord)
        print("number of selected runs : {}".format(
            len(restricted_hull_runs)))

        # RANGE FILTERING ============

        restricted_range_runs = range_filtering(restricted_hull_runs, x_coord)

        # INDIVIDUAL SELECTION ============
        restricted_idv_runs = idv_filtering(restricted_range_runs)

        print(
            "nb of structures selected : {0} ".format(
                len(restricted_idv_runs)))

        if len(restricted_idv_runs) > 0:
            selection_loop = False
        else:
            print("================================\n",
                  " Bad selection, please do it again\n",
                  "================================"
                  )

    restricted_idv_runs = sort_run(restricted_idv_runs, sort_key=x_coord)

    return restricted_idv_runs


def sort_run(rundict_list, sort_key="nelect"):
    """Sort the valid vasprun according to a defined keyword

        valid sort_key = [nelect, energy_per_fu, x_na, d_oo]"""
    if sort_key is None:
        sort_key = "energy_per_fu"

    class Sort:
        def nelect(self, run):
            return run.parameters['incar'].get("nelect", 0)

        def energy_per_fu(self, run):
            return run.energy_per_fu

        def x_na(self, run):
            return run.x_na

        def d_oo(self, run):
            return run.d_oo

    print("sorting by {}".format(sort_key))
    sort = Sort()
    sorted_run_list = sorted(rundict_list,
                             key=getattr(sort, sort_key))
    return sorted_run_list


def idv_filtering(restricted_range_runs):
    restricted_idv_runs = []
    if input("individual run  selection ? [Y]/[n] ") == "Y":
        for run in restricted_range_runs:
            if input(
                    "keep : {} ? Y/y".format(run.name_tag)) in ["Y", "y"]:
                restricted_idv_runs.append(run)
    else:
        restricted_idv_runs = restricted_range_runs
    return restricted_idv_runs


def range_filtering(restricted_hull_runs, x_coord="x_na"):
    " restrict the runs according to the values of x_coord"
    restricted_range_runs = restricted_hull_runs
    try:
        x_coord_values = {getattr(d, x_coord) for d in restricted_hull_runs}
        if len(x_coord_values) > 1:
            [x_min, x_max] = [min(x_coord_values), max(x_coord_values)]
            if input("Change current {} range [{},{}] ? [Y]".
                     format(x_coord, x_min, x_max))[0] in ["Y", "y"]:

                x_min = float(input("X min : "))
                x_max = float(input("X max : "))
                restricted_range_runs = [r for r in restricted_hull_runs
                                         if (getattr(r, x_coord) >= x_min
                                             and getattr(r, x_coord) <= x_max)]
    except Exception as ex:
        print("{} : no range filtering".format(ex))
    print("number of selected runs : {}".format(
        len(restricted_range_runs)))
    return restricted_range_runs


def hull_filtering(sieve_lvl, restricted_stack_runs, x_coord="x_na"):
    """
    filter rundicts based en convergence
    hull filtering uses a user defined x_axis against energy
    """
    selected_runs = []
    if sieve_lvl <= 3:
        selected_runs = [
            d for d in restricted_stack_runs if d.status >= sieve_lvl]

    if sieve_lvl > 3:

        if len({getattr(d, x_coord) for d in restricted_stack_runs}) == 1:
            print("not enough runs for that coord !",
                  "no filtering performed")
            selected_runs = restricted_stack_runs
        else:
            restricted_stack_runs = generate_hull_tags(
                restricted_stack_runs,
                coord=x_coord)
            selected_runs = [
                d for d in restricted_stack_runs if d.status >= sieve_lvl]
    return selected_runs


def select_x_coord():
    print("select your x axis : ",
          " // ".join(["[x]_na", "[d]_OO", "[n]elect"]))
    try:
        list_choice = str(input("filter structure list ? :"))
        assert (len(list_choice) > 0), "Invalid Input"
    except AssertionError as msg:
        print(msg, "default to x_na")
        list_choice = "x"
    x_coord = {"x": "x_na", "d": "d_oo", "n": "nelect"}[list_choice]
    return x_coord


def stacking_filter(all_runs_input):
    "STACKING FAMILY FILTER "

    stacking_set = {e.stacking for e in all_runs_input}
    if len(stacking_set) < 2:
        return all_runs_input
    family_array = stacking_set.add("finished")
    print(family_array)
    family_choice = []
    print("which stacking types to plot ? ",
          "\n".join(iter(family_array)),
          "\nNo choice = No filtering")
    while "finished" not in family_choice:
        try:
            family_choice.append(
                family_array[int(input(
                    "add stacking to current choice ({}) ? : ".format(
                        family_choice)))])
        except Exception as ex:
            print("family choice terminated {}".format(ex))
            family_choice.append('finished')
    if len(family_choice) > 1:
        restricted_stack_runs = [r for r in all_runs_input
                                 if r.stacking in family_choice]
    else:
        print("No stacking filter")
        restricted_stack_runs = all_runs_input
    return restricted_stack_runs


def select_sieve_level():

    # 3 converged plots, 4 : all minimas (even off_hull) 5 : On hull
    # minima
    list_choice = None
    choice_dict = {
        "[a]ll": "all runs that have at least a POSCAR ",
        "[c]onverged": " all converged vasprun sorted by x_na then energy",
        "[m]inima": " structures of lowest energy for each x_na ",
        "[h]ull": " structures on the convex hull"
    }
    try:
        while list_choice is None:
            print("Avaliable filtering : "+" // ".join(choice_dict.keys()))
            list_choice = input("filter structure list ? :")
            choice_dict = {"a": 1, "c": 3, "m": 4, "h": 5}
            if list_choice in choice_dict.keys():
                sieve_lvl = choice_dict[list_choice]
            else:
                print("\n".join(["{} : {}".format(*kv)
                                 for kv in choice_dict.items()]))
                list_choice = None
    except Exception as ex:
        print("{} : no hull filtering".format(ex))
        print(traceback.format_exc())
        sieve_lvl = 1

    return sieve_lvl
