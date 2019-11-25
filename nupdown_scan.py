# nupdown_scan.py
""" Plots 2D heatmaps (svg) and 3d surfaces (html) from abstract run data"""
# from multiprocessing import Pool, cpu_count

import os

# import matplotlib.cm as cm
# from matplotlib.colors import LogNorm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
# generric libraries
import numpy as np
# pymatgen libraries
# from pymatgen.io.vasp.outputs import Oszicar
from scipy.interpolate import griddata

# personnal libraries
# import readRun_entries as read
import generic_plot as g_plot
import platform_id

# from skimage import measure


try:
    # plotting libraries
    # import plotly.offline as po
    # import plotly.plotly as py
    import plotly.graph_objects as go
    import plotly.io as p_io
    # from plotly import tools
    plotly_available = True
except Exception as ex:
    print(ex)
    plotly_available = False

# === MAIN FUNCTION ===========


def plot_all_graphs(rundict_list, nupdown=None, nelect=None):
    g_plot.set_mpl_rc_params()
    plot_choice = {}
    for param, param_str in zip([nupdown, nelect], ['nupdown', 'nelect']):
        if param is None:
            try:
                input_plot_choice = str(
                    input("plot {} graphs ? [Y]/n".format(param_str)))
                assert len(input_plot_choice) > 0, "incorrect length"
                plot_choice[param_str] = bool(input_plot_choice[0] == "Y")
            except AssertionError:
                print('default to False')
                plot_choice[param_str] = False
        else:
            plot_choice[param_str] = param

    if plot_choice['nupdown']:
        selected_runs = prepare_nupdown_data(rundict_list)
        if input("Plot mag 3d surface ? [Y]") in ["Y", "y"]:
            colorcode = "O_mag" if input("colorcode ? O_mag ") == "Y" \
                else "energy"
            plot_mag_surface(selected_runs, colorcode=colorcode)

        if input("Plot mag heatmap ? [Y]") in ["Y", "y"]:
            selected_runs = [s for s in selected_runs if s["x_na"] > 0.28]
            plot_mag_heatmap(selected_runs)

    if plot_choice['nelect']:
        if input("Plot nelect heatmap ? [Y]") in ["Y", "y"]:
            plot_nelect_heatmap(rundict_list)


# === GETTING DATA =============


def prepare_nupdown_data(rundict_list):

    sorted_runs = sorted(rundict_list,
                         key=lambda k: (k.x_na, k.mag, k.energy_per_fu))
    selected_runs = [sorted_runs[0]]
    for run in sorted_runs:
        if run.mag != selected_runs[-1].mag or run.x_na != selected_runs[-1].x_na:
            selected_runs.append(run)
    for run in selected_runs:
        print(run.x_na, run.mag, run.stacking)
    return selected_runs

# === FORMATTING DATA =========


def build_x_y_e(rundict_list, attr_x, attr_y):
    """
    build a X / Y / E array to be plotted
    define a vertical slice (constant X) : all structure with same NELECT (x axis)
    scale to zero energy the bottom of each slice
    """
    for attr_x_value in {getattr(d, attr_x) for d in rundict_list}:
        struct_list_slice = [
            d for d in rundict_list if getattr(d, attr_x) == attr_x_value]
        e_min = min([d.energy_per_fu for d in struct_list_slice])
        print(" {}={} :  {} runs, E min = {:.2f} ".
              format(attr_x, attr_x_value, len(struct_list_slice), e_min))

        for rundict in struct_list_slice:
            rundict.e_valley = (rundict.energy_per_fu -
                                e_min)*1000   # in meV

    x_y_e = np.array(
        [[getattr(d, attr_x), getattr(d, attr_y), d.e_valley] for d in rundict_list])
    return x_y_e


def build_x_y_c(rundict_list, attr_x, attr_y,
                specie="O", attr_z="charge", attr_z_func="min"):
    """
    build a X / Y / CHARGE array to be plotted
    (specie, attr_z, value_type, cmap)) = ("O", "charge", "min", "Reds")
    """
    x_y_c = []
    for rundict in rundict_list:
        value_of_struct = 0
        nb_species = len(rundict.structure.indices_from_symbol(specie))

        if nb_species == 0:
            print(" no {} in the structure !!".format(specie))
            value_of_struct = np.nan
            continue
        print("{} {} in the struct {}".format(
            nb_species, specie, rundict.name_tag))
        # get the value for each non equiv pt
        try:
            sites_values = [s.properties[attr_z]
                            for s in rundict.structure_data.sites
                            if s.specie.name == specie]
        except KeyError:
            print("no {} in {}".format(attr_z, rundict.name_tag))
        else:
            # perform reduction of data (min/max/mean/std) following attr_z_func
            print("retrieved {} of {} ".format(attr_z, specie))
            value_of_struct = getattr(np, attr_z_func)(sites_values)
            print("computed {} of {} of {} ".format(
                attr_z_func, attr_z, specie))
            x_y_c.append([getattr(rundict, attr_x),
                          getattr(rundict, attr_y),
                          value_of_struct])
    print("all structures have been treated")
    all_charge_array = np.array(x_y_c)
    print("All charge array : shape {0.shape},  size : {0.size}".format(
        all_charge_array))

    assert all_charge_array.size > 0, "not enough charge data, Bader data probably missing"

    return all_charge_array


def interp_xye(x_y_z):
    "interpolation of the grid of computed points"

    def make_linspace(x_1d_vect):
        return np.linspace(min(x_1d_vect), max(x_1d_vect), num=100)

    grid_x, grid_y = np.meshgrid(make_linspace(x_y_z[:, 0]),
                                 make_linspace(x_y_z[:, 1]))
    # print( XYE[:,2])
    interp_e = griddata(x_y_z[:, 0:2], x_y_z[:, 2],
                        (grid_x, grid_y), method='cubic')
    min_e = np.nan_to_num(interp_e).min()
    max_e = np.nan_to_num(interp_e).max()
    # interp_e = interp_e - min_e
    print("min E {:.4f} , max E {:.4f}".format(min_e, max_e))
    return grid_x, grid_y, interp_e


# === PLOTTING DATA ===========

def plot_mag_heatmap(selected_runs):
    """
    specific caller to plot_abstract_heatmap
    draw x_na vs mag energy heatmap and other bader heatmaps
    """
    landscapes_to_draw = [("O", "charge", "min", "Reds"),
                          ("Mn", "charge", "max", "Greens")]
    attr_x = "x_na"
    attr_y = "mag"

    plot_abstract_heatmap(selected_runs,
                          attr_x=attr_x, attr_y=attr_y,
                          bader_landscapes=landscapes_to_draw)
    # if specie == "O" and value == "charge" and value_type == "min":
    # min_value = 6.6
    # max_value = 7.2
    # bounds = np.arange(min_value, max_value, 0.05)
    # print(bounds)
    # norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
    # axe_img = axes[-1].contourf(grid_x,
    #                             grid_y,
    #                             surface_color,
    #                             bounds,
    #                             cmap=cmap,
    #                             norm=norm,
    #                             extend='both')
    # axe_img.cmap.set_under('white')
    # axe_img.cmap.set_over('black')

    # elif specie == "Mn" and value == "charge" and value_type == "max":
    # min_value = 5.25
    # max_value = 5.38
    # bounds = np.arange(min_value, max_value, 0.01)
    # norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
    # axe_img = axes[-1].contourf(grid_x,
    #                             grid_y,
    #                             surface_color,
    #                             bounds,
    #                             cmap=cmap,
    #                             norm=norm,
    #                             extend='both')
    # axe_img.cmap.set_under('white')
    # axe_img.cmap.set_over('black')
    # axes[-1].set_xlabel('X$_{Na}$')
    # axes[-1].set_ylabel('S$_{Z}$')


def plot_nelect_heatmap(input_rundict_list):
    """
    specific caller to plot_abstract_heatmap
    draw nelect vs doo energy heatmap and other bader heatmaps
    """
    attr_x = "nelect"
    attr_y = "doo"
    landscapes_to_draw = [("S", "charge", "mean", "Reds"),
                          ("H", "charge", "mean", "Blues")]

    rundict_list = ([d for d in input_rundict_list if (
        d.status >= 3)])
    for rundict in rundict_list:
        o_x = rundict.structure.indices_from_symbol("O")
        if len(o_x) == 0:
            # print("sulfide !")
            o_x = rundict.structure.indices_from_symbol("S")
        rundict.doo = rundict.structure.distance_matrix[o_x[0], o_x[1]]

    plot_abstract_heatmap(input_rundict_list,
                          attr_x=attr_x, attr_y=attr_y,
                          bader_landscapes=landscapes_to_draw)


def plot_abstract_heatmap(rundict_list, attr_x="nelect", attr_y="doo", bader_landscapes=None):
    """
    generic method to draw energy heatmap and other bader heatmap if required
    """
    plot_data = []

    x_y_e = build_x_y_e(rundict_list, attr_x, attr_y)
    plot_data.append({"title": 'landscape_energy',
                      "data": x_y_e,
                      "cmap": 'coolwarm'})

    if bader_landscapes is not None:
        short_2_long = {"std": "Standard Deviation",
                        "min": "Minimum",
                        "max": "Maximum",
                        "mean": "Average"}
        for(specie, attr_z, attr_z_func, cmap) in bader_landscapes:
            try:
                x_y_c = build_x_y_c(rundict_list, attr_x, attr_y,
                                    specie=specie, attr_z=attr_z,
                                    attr_z_func=attr_z_func)
            except AssertionError as ex:
                # the assert in "build_x_y_c" check that bader data is present
                # if raised, we do not plot the axe
                print(ex)
            else:
                plot_title = 'Landscape_{}_of_{}_{}'.format(
                    short_2_long[attr_z_func], specie, attr_z)
                plot_data.append(
                    {"title": plot_title, "data": x_y_c, "cmap": cmap})

    nb_plots = len(plot_data)
    fig = plt.figure("2D_landscapes", figsize=(7*nb_plots, 7))
    axe = fig.add_subplot(1, nb_plots, 1)
    data_dict = plot_data[0]
    plot_energy_landscape(data_dict["data"], fig, axe,
                          attr_x=attr_x, attr_y=attr_y,
                          plot_title=data_dict["title"],
                          cmap=data_dict["cmap"])

    for (j, data_dict) in enumerate(plot_data[1:]):
        axe = fig.add_subplot(1, nb_plots, j+2)
        plot_charge_landscape(data_dict["data"], fig, axe,
                              attr_x=attr_x, attr_y=attr_y,
                              plot_title=data_dict["title"],
                              cmap=data_dict["cmap"])

    fig.tight_layout()
    plt.show(block=False)


def plot_energy_landscape(x_y_e, fig, axe, attr_x="nelect", attr_y="doo", plot_title="", cmap='coolwarm'):
    " plot energy heatmap using abstract attributes for x and y"

    # interpolation of the grid of computed points
    grid_x, grid_y, interp_e = interp_xye(x_y_e)
    # surface_color = interp_e
    # colorbar_title = "Energy"

    bounds = np.linspace(0, 3000, num=16, endpoint=True)
    interp_nice = interp_e[~np.isnan(interp_e)]
    bounds = np.linspace(
        interp_nice.min(), interp_nice.max()/2, num=41, endpoint=True)

    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
    e_img = axe.contourf(grid_x,
                         grid_y,
                         interp_e,
                         bounds,
                         cmap=cmap,
                         norm=norm,
                         extend='max')
    # axe_img.cmap.set_under('white')
    e_img.cmap.set_over('xkcd:brick red')  # color for E above max defined

    axe.scatter(x_y_e[:, 0], x_y_e[:, 1], c="black", marker="+",
                s=30, label="computed structures",
                alpha=0.3)
    # cbar = \
    fig.colorbar(e_img, ax=axe)
    #  ,  norm=colors.PowerNorm(gamma=1./3.) ) //,
    # axes[-1].imshow(interp_E, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)),
    #                             cmap=cm.jet) #, norm=LogNorm())

    # contours = plt.contour(E_img, levels=[25],colors='r' )
    # plt.clabel(E_img, [25], inline=True, fmt=["25 meV"], fontsize=10)
    # Add the contour line levels to the colorbar
    # cbar.add_lines(contours)

    axe.set_xlabel(attr_x)
    axe.set_ylabel(attr_y)
    axe.title.set_text('Energy landscape')


def plot_charge_landscape(x_y_c, fig, axe,
                          attr_x="nelect", attr_y="doo",
                          plot_title="", cmap='Reds'):
    " plot charge heatmap using abstract attributes for x and y"

    # may raise an AssertionError catched in plot_abstract_heatmap

    # interpolation of the grid of computed points
    grid_x, grid_y, interp_c = interp_xye(x_y_c)

    # surface_color = interp_e
    # colorbar_title = "Energy"

    # removing all NaN values from the array
    interp_nice = interp_c[~np.isnan(interp_c)]
    min_e = interp_nice.min()
    max_e = interp_nice.max()
    bounds = np.linspace(min_e, max_e, num=21, endpoint=True)
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
    c_img = axe.contourf(grid_x,
                         grid_y,
                         interp_c,
                         bounds,
                         cmap=cmap,
                         norm=norm,
                         extend='max')
    c_img.cmap.set_under('white')
    c_img.cmap.set_over('black')  # color for E above max defined

    axe.scatter(x_y_c[:, 0], x_y_c[:, 1], c="black", marker="+",
                s=30, label="computed structures",
                alpha=0.3)
    # cbar =  \
    fig.colorbar(c_img, ax=axe)
    #  ,  norm=colors.PowerNorm(gamma=1./3.) ) //,
    # axes[-1].imshow(interp_E, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)),
    #                             cmap=cm.jet) #, norm=LogNorm())

    # contours = plt.contour(E_img, levels=[25],colors='r' )
    # plt.clabel(E_img, [25], inline=True, fmt=["25 meV"], fontsize=10)
    # Add the contour line levels to the colorbar
    # cbar.add_lines(contours)

    axe.set_xlabel(attr_x)
    axe.set_ylabel(attr_y)

    axe.title.set_text(plot_title)
    return True


# === PLOTLY REQUIRED
def plot_mag_surface(rundict_list, colorcode="O_mag", attr_x="x_na", attr_y="mag"):

    if not plotly_available:
        print("plotly could not be imported, \n 3D plots not available")
        return False

    x_y_e = build_x_y_e(rundict_list, attr_x, attr_y)

    # interpolation of the grid of computed points
    grid_x, grid_y, interp_e = interp_xye(x_y_e)

    colorbar_title = "Energy"
    file_name = 'Sz_landscape_energy.html'
    # if colorcode == "O_mag":
    specie = "O"
    attr_z = "charge"
    attr_z_func = "min"
    # cmap = 'Reds'

    try:
        x_y_c = build_x_y_c(rundict_list, attr_x, attr_y,
                            specie=specie, attr_z=attr_z,
                            attr_z_func=attr_z_func)
        print("generated x_y_c")
        grid_x, grid_y, interp_c = interp_xye(x_y_c)
        print("interpolation done")
        surface_color = interp_c
        colorbar_title = "{} {} of {}".format(attr_z_func, attr_z, specie)
    except KeyError as ex:
        print("key {} not found".format(ex),
              "fallback surface color = energy")
        surface_color = interp_e
        colorbar_title = "Energy"

    # defining text on hover
    textz = [
        ["".join([
            'x_na: {:0.2f}'.format(grid_x[i][j]),
            '<br>Sz: {:0.2f}'.format(grid_y[i][j]),
            '<br>E: {:0.5f}'.format(interp_e[i][j]),
            '<br>{}: {:0.5f}'.format(
                colorbar_title, surface_color[i][j])
        ])

            for j in range(grid_x.shape[1])
        ]

        for i in range(grid_x.shape[0])
    ]
    # main surface with colors
    data_surface = go.Surface(
        x=grid_x,
        y=grid_y,
        z=interp_e,
        surfacecolor=surface_color,
        colorscale='Jet',
        reversescale=True,
        cmax=7,
        cmin=6.3,
        cmid=6.8,
        colorbar={"title": colorbar_title, "titleside": 'top'},

        text=textz,
        hoverinfo='text',
        contours=dict(
            x={"show": False, "highlight": False},
            y={"show": False, "highlight": False},
            z={"show": True,
               "start": 25,
               "size": 100,
               "color": 'rgb(111,151,255)',
                "width": 1,
                "highlight": False}
        ),
        lighting=dict(
            ambient=0.8,  # 0.8,
            specular=0.1,  # 0.05,
            diffuse=0.7,  # 0.8,
            roughness=0.3  # 0.5
        )

    )
    data = [data_surface]

    orange_dots = False
    if orange_dots:
        # orange dots to find computed points
        data_scatter = go.Scatter3d(
            x=x_y_e[:, 0],
            y=x_y_e[:, 1],
            z=x_y_e[:, 2],
            mode='markers',
            marker=dict(
                size=5,
                symbol="diamond"
            )
        )
        data.append(data_scatter)

    axis = dict(
        title='x Axis',
        showbackground=False,
        backgroundcolor="rgb(230, 230,230)",
        gridcolor="rgb(255, 255, 255)",
        zerolinecolor="rgb(255, 255, 255)",
    )

    scene = dict(
        xaxis=dict(axis), yaxis=dict(axis), zaxis=dict(axis),
        aspectratio=dict(x=2, y=1.3, z=1)
    )

    scene['xaxis']['title'] = "Na content"
    scene['yaxis']['title'] = "constrained NUPDOWN <br>Sz / UF"
    scene['zaxis']['title'] = "Energy <br> (meV)"

    layout = go.Layout(
        title='Constrained magmom energy landscape',
        autosize=True,
        scene=scene
    )

    fig = go.Figure(data=data, layout=layout)

    camera = dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=-0.5),
        eye=dict(x=1.6, y=-0.1, z=0.3)
    )
    fig.update_layout(scene_camera=camera)
    file_name = 'Sz_landscape_{}.html'.format(colorbar_title)
    path_to_fig = os.path.join(platform_id.figure_dir(), file_name)
    p_io.write_html(fig, path_to_fig, auto_open=False)
    print("wrote fig at : {}".format(path_to_fig))

    for ext in ["svg", "pdf"]:
        file_name = 'Sz_landscape_{}.{}'.format(colorbar_title, ext)
        path_to_fig = os.path.join(platform_id.figure_dir(), file_name)
        fig.write_image(path_to_fig)
        print("wrote fig at : {}".format(path_to_fig))

    # po.plot(fig, filename=file_name)


# DEPRECATED ===========
# ======================

    # def plot_mag_heatmap(input_rundict_list):

    #     landscapes_to_draw = [("O", "charge", "min", "Reds"),
    #                           ("Mn", "charge", "max", "Greens")]
    #     # ("Zn", "charge","min","Purples")]
    #     nb_plots = len(landscapes_to_draw) + 1
    #     axes = []

    #     rundict_list = ([d for d in input_rundict_list if (
    #         d.status >= 3 and d.mag >= 2 and d.mag <= 3.7)])

    #     for x_na in set([d.x_na for d in rundict_list]):
    #         struct_list_slice = [d for d in rundict_list if d.x_na == x_na]
    #         e_min = min([d.energy_per_fu for d in struct_list_slice])
    #         print(
    #             " x_na ={} :  {} runs, E min = {:.2f} ".format(
    #                 x_na,
    #                 len(struct_list_slice),
    #                 e_min))
    #         for d in struct_list_slice:
    #             d.e_valley = 1000 * (d.energy_per_fu - e_min)

    #     x_y_e = np.array([[d.x_na, d.mag, d.e_valley] for d in rundict_list])
    #     # XYE =  np.flipud( XYE)
    #     # interpolation of the grid of computed points
    #     x = np.linspace(min(x_y_e[:, 0]), max(x_y_e[:, 0]), num=100)
    #     y = np.linspace(min(x_y_e[:, 1]), max(x_y_e[:, 1]), num=100)

    #     grid_x, grid_y = np.meshgrid(x, y)
    #     # print( XYE[:,2])
    #     interp_e = griddata(x_y_e[:, 0:2], x_y_e[:, 2],
    #                         (grid_x, grid_y), method='cubic')
    #     min_e = np.nan_to_num(interp_e).min()
    #     interp_e = interp_e - min_e
    #     surface_color = interp_e
    #     colorbar_title = "Energy"
    #     file_name = 'Sz_landscape_energy'

    #     # print(interp_E)
    #     print("min E", interp_e.min(), "max E", interp_e.max())
    #     fig = plt.figure(file_name, figsize=(15, 5))
    #     axes.append(fig.add_subplot(1, nb_plots, 1))
    #     bounds = np.array([0, 10, 25, 50, 100, 150, 200, 300, 400])
    #     norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
    #     e_img = axes[-1].contourf(grid_x,
    #                               grid_y,
    #                               interp_e,
    #                               bounds,
    #                               cmap='coolwarm',
    #                               norm=norm,
    #                               extend='max')
    #     # axe_img.cmap.set_under('white')
    #     e_img.cmap.set_over('black')

    #     # CS3.cmap.set_under('yellow')
    #     # CS3.cmap.set_over('cyan')
    #     #
    #     cbar = fig.colorbar(e_img, ax=axes[-1])
    #     #  ,  norm=colors.PowerNorm(gamma=1./3.) ) //,
    #     # axes[-1].imshow(interp_E, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)),
    #     #                             cmap=cm.jet) #, norm=LogNorm())

    #     # contours = plt.contour(E_img, levels=[25],colors='r' )
    #     # plt.clabel(E_img, [25], inline=True, fmt=["25 meV"], fontsize=10)
    #     # Add the contour line levels to the colorbar
    #     # cbar.add_lines(contours)
    #     axes[-1].set_xlabel('X$_{Na}$')
    #     axes[-1].set_ylabel('S$_{Z}$')
    #     axes[-1].title.set_text('Energy landscape')

    #     short_2_long = {
    #         "dev": "Standard Deviation",
    #         "min": "Minimum",
    #         "max": "Maximum",
    #     }
    #     print(short_2_long)

    #     for(j, (specie, value, value_type, cmap)) in enumerate(landscapes_to_draw):
    #         try:
    #             ax_nb = j + 2

    #             # specie = "O"
    #             # value = "charge"
    #             # value_type = "min"
    #             value_list = np.zeros(len(rundict_list))
    #             for i, run in enumerate(rundict_list):
    #                 nbSpecie = len(run.structure.indices_from_symbol(specie))
    #                 if nbSpecie == 0:
    #                     print(" no {} in the structure !!".format(specie))
    #                     value_list[i] = np.nan
    #                     break
    #                 else:
    #                     # get the value for each non equiv pt

    #                     sites_values = [s.properties[value]
    #                                     for s in run.structure_data.sites
    #                                     if s.specie.name == specie]

    #                     if value_type == "mean":
    #                         value_list[i] = np.sum(sites_values) / nbSpecie
    #                     elif value_type == "min":
    #                         value_list[i] = min(sites_values)
    #                     elif value_type == "max":
    #                         value_list[i] = max(sites_values)
    #                     elif value_type == "dev":  # variance
    #                         value_list[i] = np.sqrt(np.sum(np.square(
    #                             sites_values - np.average(sites_values))
    #                         ) / nbSpecie)

    #             if np.isnan(value_list).any() is False:
    #                 surface_color = griddata(
    #                     x_y_e[:, 0:2], value_list, (grid_x, grid_y), method='cubic')
    #                 # print(surface_color)
    #                 # colorbar_title = "{} {} of {}".format(value_type, value,specie)
    #                 file_name = 'Landscape_{}_of_{}_{}'.format(
    #                     short_2_long[value_type], specie, value)

    #             axes.append(fig.add_subplot(1, nb_plots, ax_nb))
    #             # define preselected contour values for wanted plots
    #             surface_color_norm = surface_color[~np.isnan(surface_color)]
    #             if specie == "O" and value == "charge" and value_type == "min":
    #                 min_value = 6.6
    #                 max_value = 7.2
    #                 bounds = np.arange(min_value, max_value, 0.05)
    #                 print(bounds)
    #                 norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
    #                 axe_img = axes[-1].contourf(grid_x,
    #                                             grid_y,
    #                                             surface_color,
    #                                             bounds,
    #                                             cmap=cmap,
    #                                             norm=norm,
    #                                             extend='both')
    #                 axe_img.cmap.set_under('white')
    #                 axe_img.cmap.set_over('black')

    #             elif specie == "Mn" and value == "charge" and value_type == "max":
    #                 min_value = 5.25
    #                 max_value = 5.38
    #                 bounds = np.arange(min_value, max_value, 0.01)
    #                 norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
    #                 axe_img = axes[-1].contourf(grid_x,
    #                                             grid_y,
    #                                             surface_color,
    #                                             bounds,
    #                                             cmap=cmap,
    #                                             norm=norm,
    #                                             extend='both')
    #                 axe_img.cmap.set_under('white')
    #                 axe_img.cmap.set_over('black')

    #             else:
    #                 axe_img = axes[-1].contourf(grid_x,
    #                                             grid_y, surface_color, cmap=cmap)

    #             # axe.imshow(surface_color,
    #                       extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)) )

    #             cbar = fig.colorbar(axe_img, ax=axes[-1])

    #             axes[-1].set_xlabel('X$_{Na}$')
    #             axes[-1].set_ylabel('S$_{Z}$')
    #             axes[-1].title.set_text('{} of {} {}'.format(
    #                 short_2_long[value_type], specie, value))
    #         except KeyError as key_err:
    #             print(
    #                 "No value for the key {}".format(key_err))
    #             print(": you should probably compute bader charges ! ")

    #     fig.tight_layout()
    #     plt.show(block=False)
