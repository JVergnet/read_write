# nupdown_scan.py

from multiprocessing import Pool, cpu_count

# import matplotlib.cm as cm
# from matplotlib.colors import LogNorm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
# generric libraries
import numpy as np
# pymatgen libraries
from pymatgen.io.vasp.outputs import Oszicar
from scipy.interpolate import griddata

# personnal libraries
# import readRun_entries as read
import generic_plot as generic_plot

# from skimage import measure


try:
    # plotting libraries
    import plotly.offline as po
    # import plotly.plotly as py
    import plotly.graph_objs as go
    # from plotly import tools
    plotly_available = True
except Exception as ex:
    print(ex)
    plotly_available = False




def get_mag_tag_list(vaspRunDictList):

    vasp_run_poll = []
    with Pool(processes=cpu_count()) as p:
        vasp_run_poll = p.map(get_mag_tag_single, vaspRunDictList)
        p.close()
        p.join()

    return(vasp_run_poll)


def get_mag_tag_single(runDict):
    if not hasattr(runDict, "mag"):
        runDict.mag = Oszicar(runDict.jobFolder + "/OSZICAR"
                              ).ionic_steps[-1]["mag"] /\
            runDict.nb_cell
    return(runDict)


def plot_mag_graphs(runDict_list):
    generic_plot.set_mpl_rc_params()
    nupdown = False
    if nupdown:
        runDict_list = get_mag_tag_list(runDict_list)
        sorted_runs = sorted(
            runDict_list,
            key=lambda k: (
                k.xNa,
                k.mag,
                k.energy_per_FU))
        selected_runs = [sorted_runs[0]]
        for r in sorted_runs:
            if r.mag != selected_runs[-1].mag or r.xNa != selected_runs[-1].xNa:
                selected_runs.append(r)
        for r in selected_runs:
            print(r.xNa, r.mag, r.stacking)

        # run_index = {}

        if input("Plot mag 3d surface ? [Y]") in ["Y", "y"]:
            colorcode = "O_mag" if input(
                "colorcode ? O_mag ") == "Y" else "energy"
            plot_mag_surface(selected_runs, colorcode=colorcode)

        if input("Plot mag heatmap ? [Y]") in ["Y", "y"]:
            selected_runs = [s for s in selected_runs if s["xNa"] > 0.28]
            plot_mag_heatmap(selected_runs)

    if input("Plot nelect heatmap ? [Y]") in ["Y", "y"]:
        plot_nelect_heatmap(runDict_list)


def plot_nelect_heatmap(input_runDict_list):

    nb_plots = 1
    axes = []

    runDict_list = ([d for d in input_runDict_list if (
        d.status >= 3)])
    for runDict in runDict_list:
        runDict.nelect = runDict.parameters["incar"]['NELECT']
        o_1, o_2 = runDict.structure.indices_from_symbol("O")
        runDict.doo = runDict.structure.distance_matrix[o_1, o_2]

    # define a vertical slice : all structure with same NELECT (x axis)
    for nelect in set([d.nelect for d in runDict_list]):
        struct_list_slice = [d for d in runDict_list if d.nelect == nelect]
        E_min = min([d.energy_per_FU for d in struct_list_slice])
        print(" nelect={} :  {} runs, E min = {:.2f} ".
              format(nelect, len(struct_list_slice), E_min))

        for d in struct_list_slice:
            d.e_valley = (d.energy_per_FU - E_min)  # in meV

    XYE = np.array([[d.nelect, d.doo, d.e_valley] for d in runDict_list])
    # XYE =  np.flipud( XYE)

    # interpolation of the grid of computed points
    x = np.linspace(min(XYE[:, 0]), max(XYE[:, 0]), num=100)
    y = np.linspace(min(XYE[:, 1]), max(XYE[:, 1]), num=100)

    grid_x, grid_y = np.meshgrid(x, y)
    #print( XYE[:,2])
    interp_E = griddata(XYE[:, 0:2], XYE[:, 2],
                        (grid_x, grid_y), method='linear')
    min_E = np.nan_to_num(interp_E).min()
    max_E = np.nan_to_num(interp_E).max()
    interp_E = interp_E - min_E
    surface_color = interp_E
    colorbar_title = "Energy"
    file_name = 'OO_landscape_energy'

    # print(interp_E)
    print("min E {} , max E {}".format(min_E, max_E))
    fig = plt.figure(file_name, figsize=(11, 9))
    axes.append(fig.add_subplot(1, nb_plots, 1))
    bounds = np.linspace(0, 5, num=21, endpoint=True)
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
    E_img = axes[-1].contourf(grid_x,
                              grid_y,
                              interp_E,
                              bounds,
                              cmap='coolwarm',
                              norm=norm,
                              extend='max')
    # axe_img.cmap.set_under('white')
    E_img.cmap.set_over('black')

    # CS3.cmap.set_under('yellow')
    # CS3.cmap.set_over('cyan')
    #
    cbar = fig.colorbar(E_img, ax=axes[-1])
    #  ,  norm=colors.PowerNorm(gamma=1./3.) ) //,
    # axes[-1].imshow(interp_E, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)),
    #                             cmap=cm.jet) #, norm=LogNorm())

    # contours = plt.contour(E_img, levels=[25],colors='r' )
    # plt.clabel(E_img, [25], inline=True, fmt=["25 meV"], fontsize=10)
    # Add the contour line levels to the colorbar
    # cbar.add_lines(contours)
    axes[-1].set_xlabel('$N_{electrons}$')
    axes[-1].set_ylabel('$d_{OO}$')
    axes[-1].title.set_text('Energy landscape')
    fig.tight_layout()
    plt.show(block=False)


def plot_mag_heatmap(input_runDict_list):

    landscapes_to_draw = [("O", "charge", "min", "Reds"),
                          ("Mn", "charge", "max", "Greens")]
    # ("Zn", "charge","min","Purples")]
    nb_plots = len(landscapes_to_draw) + 1
    axes = []

    runDict_list = ([d for d in input_runDict_list if (
        d.status >= 3 and d.mag >= 2 and d.mag <= 3.7)])

    for xNa in set([d.xNa for d in runDict_list]):
        struct_list_slice = [d for d in runDict_list if d.xNa == xNa]
        E_min = min([d.energy_per_FU for d in struct_list_slice])
        print(
            " xNa ={} :  {} runs, E min = {:.2f} ".format(
                xNa,
                len(struct_list_slice),
                E_min))
        for d in struct_list_slice:
            d.e_valley = 1000 * (d.energy_per_FU - E_min)

    XYE = np.array([[d.xNa, d.mag, d.e_valley] for d in runDict_list])
    # XYE =  np.flipud( XYE)
    # interpolation of the grid of computed points
    x = np.linspace(min(XYE[:, 0]), max(XYE[:, 0]), num=100)
    y = np.linspace(min(XYE[:, 1]), max(XYE[:, 1]), num=100)

    grid_x, grid_y = np.meshgrid(x, y)
    #print( XYE[:,2])
    interp_E = griddata(XYE[:, 0:2], XYE[:, 2],
                        (grid_x, grid_y), method='cubic')
    min_E = np.nan_to_num(interp_E).min()
    interp_E = interp_E - min_E
    surface_color = interp_E
    colorbar_title = "Energy"
    file_name = 'Sz_landscape_energy'

    # print(interp_E)
    print("min E", interp_E.min(), "max E", interp_E.max())
    fig = plt.figure(file_name, figsize=(15, 5))
    axes.append(fig.add_subplot(1, nb_plots, 1))
    bounds = np.array([0, 10, 25, 50, 100, 150, 200, 300, 400])
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
    E_img = axes[-1].contourf(grid_x,
                              grid_y,
                              interp_E,
                              bounds,
                              cmap='coolwarm',
                              norm=norm,
                              extend='max')
    # axe_img.cmap.set_under('white')
    E_img.cmap.set_over('black')

    # CS3.cmap.set_under('yellow')
    # CS3.cmap.set_over('cyan')
    #
    cbar = fig.colorbar(E_img, ax=axes[-1])
    #  ,  norm=colors.PowerNorm(gamma=1./3.) ) //,
    # axes[-1].imshow(interp_E, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)),
    #                             cmap=cm.jet) #, norm=LogNorm())

    # contours = plt.contour(E_img, levels=[25],colors='r' )
    # plt.clabel(E_img, [25], inline=True, fmt=["25 meV"], fontsize=10)
    # Add the contour line levels to the colorbar
    # cbar.add_lines(contours)
    axes[-1].set_xlabel('X$_{Na}$')
    axes[-1].set_ylabel('S$_{Z}$')
    axes[-1].title.set_text('Energy landscape')

    short_2_long = {
        "dev": "Standard Deviation",
        "min": "Minimum",
        "max": "Maximum",
    }
    print(short_2_long)

    for(i, (specie, value, value_type, cmap)) in enumerate(landscapes_to_draw):
        try:
            ax_nb = i + 2

            # specie = "O"
            # value = "charge"
            # value_type = "min"
            value_list = np.zeros(len(runDict_list))
            for i, run in enumerate(runDict_list):
                nbSpecie = len(run.structure.indices_from_symbol(specie))
                if nbSpecie == 0:
                    print(" no {} in the structure !!".format(specie))
                    value_list[i] = np.nan
                    break
                else:
                    # get the value for each non equiv pt

                    sites_values = [s.properties[value]
                                    for s in run.structure.sites
                                    if s.specie.name == specie]

                    if value_type == "mean":
                        value_list[i] = np.sum(sites_values) / nbSpecie
                    elif value_type == "min":
                        value_list[i] = min(sites_values)
                    elif value_type == "max":
                        value_list[i] = max(sites_values)
                    elif value_type == "dev":  # variance
                        value_list[i] = np.sqrt(np.sum(np.square(
                            sites_values - np.average(sites_values))
                        ) / nbSpecie)

            if np.isnan(value_list).any() == False:
                surface_color = griddata(
                    XYE[:, 0:2], value_list, (grid_x, grid_y), method='cubic')
                # print(surface_color)
                #colorbar_title = "{} {} of {}".format(value_type, value,specie)
                file_name = 'Landscape_{}_of_{}_{}'.format(
                    short_2_long[value_type], specie, value)

            axes.append(fig.add_subplot(1, nb_plots, ax_nb))
            # define preselected contour values for wanted plots
            surface_color_norm = surface_color[~np.isnan(surface_color)]
            if specie == "O" and value == "charge" and value_type == "min":
                min_value = 6.6
                max_value = 7.2
                bounds = np.arange(min_value, max_value, 0.05)
                print(bounds)
                norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
                axe_img = axes[-1].contourf(grid_x,
                                            grid_y,
                                            surface_color,
                                            bounds,
                                            cmap=cmap,
                                            norm=norm,
                                            extend='both')
                axe_img.cmap.set_under('white')
                axe_img.cmap.set_over('black')

            elif specie == "Mn" and value == "charge" and value_type == "max":
                min_value = 5.25
                max_value = 5.38
                bounds = np.arange(min_value, max_value, 0.01)
                norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
                axe_img = axes[-1].contourf(grid_x,
                                            grid_y,
                                            surface_color,
                                            bounds,
                                            cmap=cmap,
                                            norm=norm,
                                            extend='both')
                axe_img.cmap.set_under('white')
                axe_img.cmap.set_over('black')

            else:
                axe_img = axes[-1].contourf(grid_x,
                                            grid_y, surface_color, cmap=cmap)

            #axe.imshow(surface_color,  extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)) )

            cbar = fig.colorbar(axe_img, ax=axes[-1])

            axes[-1].set_xlabel('X$_{Na}$')
            axes[-1].set_ylabel('S$_{Z}$')
            axes[-1].title.set_text('{} of {} {}'.format(
                short_2_long[value_type], specie, value))
        except KeyError as key_err:
            print(
                "No value for the key {} : you should probably compute bader charges ! ".format(key_err))

    fig.tight_layout()
    plt.show(block=False)


def plot_mag_surface(runDict_list, colorcode="O_mag"):
    if not plotly_available:
        print("plotly could not be imported, \n 3D plots not available")
        return(False)

    mag_list = runDict_list
    for xNa in set([d.xNa for d in mag_list]):
        struct_list_slice = [d for d in mag_list if d.xNa == xNa]
        E_min = min([d.energy_per_FU for d in struct_list_slice])
        print(
            " xNa ={} :  {} runs, E min = {:.2f} ".format(
                xNa, len(struct_list_slice), E_min))
        for d in struct_list_slice:
            d.e_valley = 1000 * (d.energy_per_FU - E_min)

    XYE = np.array([[d.xNa, d.mag, d.e_valley]
                    for d in mag_list if d.status >= 3])

    # interpolation of the grid of computed points
    x = np.linspace(min(XYE[:, 0]), max(XYE[:, 0]), num=80)
    y = np.linspace(min(XYE[:, 1]), max(XYE[:, 1]), num=50)

    grid_x, grid_y = np.meshgrid(x, y)

    interp_E = griddata(XYE[:, 0:2], XYE[:, 2],
                        (grid_x, grid_y), method='cubic')

    surface_color = interp_E
    colorbar_title = "Energy"
    file_name = 'Sz_landscape_energy.html'
    if colorcode == "O_mag":
        specie = "O"
        value = "charge"
        value_type = "min"
        value_list = np.zeros(len(mag_list))
        for i, run in enumerate(mag_list):
            nbSpecie = len(run.structure.indices_from_symbol(specie))
            if nbSpecie == 0:
                print(" no {} in the structure !!".format(specie))
                value_list[i] = np.nan
                break
            else:
                # get the value for each non equiv pt

                if value_type == "mean":
                    sites_values = [[site[value] * site['multiplicity']]
                                    for site in run['equivSiteList']
                                    if site['element'] == specie]
                    value_list[i] = np.sum(sites_values) / nbSpecie
                if value_type == "min":
                    Y_min = min([site[value]
                                 for site in run['equivSiteList']
                                 if site['element'] == specie])
                    value_list[i] = Y_min
                if value_type == "max":
                    Y_max = max([site[value]
                                 for site in run['equivSiteList']
                                 if site['element'] == specie])
                    value_list[i] = Y_max
        if not np.isnan(value_list).any():  # if no Nan value in the vector
            surface_color = griddata(
                XYE[:, 0:2], value_list, (grid_x, grid_y), method='cubic')
            colorbar_title = "{} {} of {}".format(value_type, value, specie)
            file_name = 'Sz_landscape_{}_{}_of_{}.html'.format(
                value_type, value, specie)

    textz = [
        [
            'xNa: ' +
            '{:0.2f}'.format(
                grid_x[i][j]) +
            '<br>Sz: ' +
            '{:0.2f}'.format(
                grid_y[i][j]) +
            '<br>E: ' +
            '{:0.5f}'.format(
                interp_E[i][j]) +
            '<br>{}:'.format(colorbar_title) +
            '{:0.5f}'.format(
                surface_color[i][j]) for j in range(
                grid_x.shape[1])] for i in range(
            grid_x.shape[0])]

    data_surface = go.Surface(
        x=x,
        y=y,
        z=interp_E,
        surfacecolor=surface_color,
        colorscale='Jet',
        colorbar=dict(
            title=colorbar_title,
            titleside='top'
        ),
        text=textz,
        hoverinfo='text',

    )

    data_scatter = go.Scatter3d(
        x=XYE[:, 0],
        y=XYE[:, 1],
        z=XYE[:, 2],
        mode='markers',
        marker=dict(
            size=5,
            symbol="diamond-tall"
        )
    )

    axis = dict(
        title='x Axis',
        showbackground=True,
        backgroundcolor="rgb(230, 230,230)",
        gridcolor="rgb(255, 255, 255)",
        zerolinecolor="rgb(255, 255, 255)",
    )

    data = [data_surface, data_scatter]

    scene = dict(
        xaxis=dict(axis), yaxis=dict(axis), zaxis=dict(axis),
        cameraposition=[[0.2, 0.5, 0.5, 0.2], [0, 0, 0], 4.8],
        aspectratio=dict(x=1, y=1, z=1)
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

    po.plot(fig, filename=file_name)
