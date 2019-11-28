# PES_complex_plot.py

# standard libraries
# from operator import itemgetter
# import os
# import sys
# import subprocess
import importlib

# import matplotlib as mpl
import matplotlib.pyplot as plt
# scientific libraries
import numpy as np
from scipy.interpolate import griddata

# personnal libraries
# import readRun_entries as read
import bailar_twist as bailar
# import read_hull as hull
import nupdown_scan as nupdown

try:
    # plotting libraries
    import plotly.offline as po
    # import plotly.plotly as py
    import plotly.graph_objs as go
    from plotly import tools
    from skimage import measure
except Exception as ex:
    print(ex)

# from mpl_toolkits.mplot3d import Axes3D, axes3d


def plot_3d_angle_energy_from_struct_list(struct_list, struct_mesh=None):
    if input("get mag tag ? Y/n") == "Y":
        importlib.reload(nupdown)

        prop = "mag"
    else:
        prop = "ediff"

    XYZE_relax_mean, XYZE_relax_all = get_OO_angles_prop(
        struct_list, prop=prop)
    # wether to include inddual OO pairs of just the average
    all_points = XYZE_relax_all if input(
        "include ALL points (instead of mean value) ? Y / n") == "Y" else None
    # equals None if no mesh is provided (plot relaxed path only )
    XYZE_mesh, XYZE_mesh_all = get_OO_angles_prop(
        struct_mesh, prop=prop) if struct_mesh is not None else (None, None)

    fig_done = plot_3D_angles_PES(XYZE_relax_mean, struct_list,
                                  XYZE_mesh=XYZE_mesh, XYZE_relax_all=XYZE_relax_all,)

    return(fig_done)


def get_OO_angles_prop(struct_list, prop="ediff"):

    x_y_e = nupdown.build_x_y_e(struct_list, "x_na", "mag")

    Z_bailar_all = []  # list (nb structures)  of lists (nb of angles)
    Y_trigo_all = []
    X_all = []
    E_all = []

    for i, d in enumerate(struct_list):
        struct = d["structure"]
        YZ = np.array([bailar.measure_quadruplet_angles(struct, MMOO)
                       for MMOO in d["MMOO_quadruplets"]])
        Z_bailar_all.append(YZ[:, 0])  # Z_bailar
        Y_trigo_all.append(YZ[:, 1])  # Y_trigo
        X_all.append(np.ones_like(YZ[:, 1]) * x_y_e[:, 0][i])
        E_all.append(np.ones_like(YZ[:, 1]) * x_y_e[:, 2][i])

    Z_bailar_all = np.vstack(Z_bailar_all)
    Y_trigo_all = np.vstack(Y_trigo_all)
    X_all = np.vstack(X_all)
    E_all = np.vstack(E_all)

    XYZE_all_list = [list(t) for t in zip(X_all.flatten(),
                                          Y_trigo_all.flatten(),
                                          Z_bailar_all.flatten(),
                                          E_all.flatten())]
    XYZE_all = np.vstack(XYZE_all_list)
    print(X_all.shape, Y_trigo_all.shape, Z_bailar_all.shape, E_all.shape)

    # list (nb of structures) of float (average angle)
    Z_bailar_mean = np.mean(Z_bailar_all, axis=1)
    Y_trigo_mean = np.mean(Y_trigo_all, axis=1)
    X_mean = np.mean(X_all, axis=1)
    E_mean = np.mean(E_all, axis=1)
    print(X_mean.shape, Y_trigo_mean.shape, Z_bailar_mean.shape, E_mean.shape)

    XYZE_mean_list = [list(t) for t in zip(
        X_mean, Y_trigo_mean, Z_bailar_mean, E_mean)]
    print("length of 4D coordinates (list) : {}".format(len(XYZE_mean_list)))
    XYZE_mean = np.vstack(XYZE_mean_list)
    print("shape of 4D coordinates (XYZ + energy) : {}".format(XYZE_mean.shape))
    return(XYZE_mean, XYZE_all)


def plot_3D_angles_PES(
        XYZE_relax,
        struct_list,
        XYZE_mesh=None,
        XYZE_relax_all=None):

    # preparing the plot :
    # ============================
    # PLOTLY 3D BUILD
    # ===========================
    axis = dict(
        title='x Axis',
        showbackground=True,
        backgroundcolor="rgb(230, 230,230)",
        gridcolor="rgb(255, 255, 255)",
        zerolinecolor="rgb(255, 255, 255)",
    )

    scene = dict(
        xaxis=dict(axis), yaxis=dict(axis), zaxis=dict(axis),
        cameraposition=[[0.2, 0.5, 0.5, 0.2], [0, 0, 0], 4.8],
        aspectratio=dict(x=1, y=1, z=1)
    )

    scene['xaxis']['title'] = "Na content"
    scene['yaxis']['title'] = "vertical O-O distance (in $\AA$)"
    scene['zaxis']['title'] = "bailar angle (in °)"

    data_relax_list = [get_scatterplot_data(
        XYZE_relax, midpoint=100, numtick=6)]

    if XYZE_relax_all is not None:
        MMOO_list = []
        for d in struct_list:
            vesta_offset = len(d["structure"]) - \
                len(d["structure"].indices_from_symbol("O"))
            for MMOO in d["MMOO_quadruplets"]:
                MMOO_list. append("{}".format(
                    np.array(MMOO['oxygen_pair']) - vesta_offset))

        #O1,O2 = MMOO['oxygen_pair']
        hover_txt = ["E =  {:.3f} meV \n MMOO = {}".format(e, MMOO)
                     for e, MMOO in zip(XYZE_relax_all[:, 3], MMOO_list)]
        print(hover_txt)
        data_relax_list.append(
            get_scatterplot_data(
                XYZE_relax_all,
                midpoint=100,
                numtick=6,
                secondary_data=True,
                hover_txt=hover_txt))

    if input("plot relax data only ? Y / n") == "Y":
        plot_relax_data_only(data_relax_list, scene)

    if XYZE_mesh is not None:
        # ============
        # Mesh data has been provided
        # ============
        data_mesh = get_scatterplot_data(XYZE_mesh, midpoint=100, numtick=6)
        if input("plot mesh data only ? Y / n") == "Y":
            plot_mesh_data_only(data_mesh, scene)

        nb_grid_pts = 100
        interp_E = interpolate_3D_PEV(XYZE_mesh, nb_grid_pts)
        if input("interpolate single isosurface ? Y / n") == "Y":
            plot_single_isosurf(XYZE_mesh, interp_E,
                                data_relax_list, scene, nb_grid_pts=100)

        if input("slider for  multiple  isosurface ? Y / n") == "Y":
            plot_slider_isosurf(
                XYZE_mesh,
                interp_E,
                data_relax_list,
                scene,
                nb_grid_pts=100,
                energy_values=[
                    10,
                    25,
                    40,
                    50,
                    75,
                    100])

        if input("plot slice in the mesh ? Y / n") == "Y":
            plot_slice_in_volume(XYZE_mesh, nb_grid_pts, scene)

    return(True)


def plot_relax_data_only(data_relax_list, scene):

    layout = go.Layout(
        margin=dict(l=0, r=0, b=0, t=0),
        height=600,
        width=800
    )
    fig = go.Figure(data=data_relax_list, layout=layout)
    fig['layout']['scene'].update(scene)
    fig['layout'].update(
        title='Energy of relaxed structures normalized by lowest energy at each x_na')
    po.plot(fig, filename='energy_angles_relax_only.html')
    return(fig)


def plot_mesh_data_only(data_mesh, scene):
    data = [data_mesh]
    layout = go.Layout(
        margin=dict(l=0, r=0, b=0, t=0),
        height=600,
        width=800
    )
    fig = go.Figure(data=data, layout=layout)
    fig['layout']['scene'].update(scene)
    fig['layout'].update(
        title='Energy of frozen structures normalized by lowest energy at each x_na')
    po.plot(fig, filename='energy_angles_mesh_only.html')
    return(fig)


def plot_slice_in_volume(XYZE, nb_grid_pts, scene):

    # define the points in the slice
    min_max_list = [[min(XYZE[:, i]), max(XYZE[:, i])] for i in range(3)]
    x = np.linspace(min_max_list[0][0], min_max_list[0][1], nb_grid_pts)
    z = np.linspace(min_max_list[2][0], min_max_list[2][1], nb_grid_pts)

    x, z = np.meshgrid(x, z)
#    alpha=np.pi/10
    # y=np.ones(x.shape)*2.0 # slice of constant trigonal distortion
    # y=-(1-x)*np.tan(alpha)+2.2
    y = -(0.66 - x) * 0.1 + 2.05
    print("grid generated")
    # interpolate the mesh data on the slice
    surfacecolor = griddata(XYZE[:, 0:3], XYZE[:, 3], (x, y, z))
    print("interpolation done")
    log_surfacecolor = np.log10(surfacecolor + 1)
    midpoint = 100
    numtick = 4
    ticks = np.concatenate([np.linspace(min(surfacecolor.flatten()),
                                        midpoint,
                                        num=numtick,
                                        endpoint=False),
                            np.linspace(midpoint,
                                        max(surfacecolor.flatten()),
                                        num=numtick)])
    print("ticks generated : {}".format(ticks))
#    sminz, smaxz= [np.min(surfacecolor),np.max(surfacecolor)]

    # d = dict(cmax=max(log_surfacecolor),
    #         cmin=min(log_surfacecolor),
    #         colorscale='Jet',   # choose a colorscale
    #         colorbar=dict(
    #             title='Energies',
    #             titleside = 'top',
    #             tickmode = 'array',
    #             tickvals = np.log10(ticks+1),
    #             ticktext = ["{:.0f}".format(t) for t in ticks],
    #             ticks = 'outside'
    #         )
    # )

    # plot the slice on a 3d fig
    surf = go.Surface(x=x,  # https://plot.ly/python/reference/#surface
                      y=y,
                      z=z,
                      surfacecolor=log_surfacecolor,
                      colorscale='Jet',
                      colorbar=dict(
                          title='Energies',
                          titleside='top',
                          tickmode='array',
                          tickvals=np.log10(ticks + 1),
                          ticktext=["{:.0f}".format(t) for t in ticks],
                          ticks='outside'
                      )



                      )
    print("surface generated")
    # ,
    #         colorscale=colorscale,
    #         showscale=showscale)
    data = [surf]
    layout = go.Layout(
        margin=dict(l=0, r=0, b=0, t=0),
        height=600,
        width=800
    )
    print("layout generated")
    fig = go.Figure(data=data, layout=layout)
    fig['layout']['scene'].update(scene)
    fig['layout'].update(title='slice in the nergy landscape')
    po.plot(fig, filename='Slice_across_energy_landscape.html')

    return(True)


def plot_single_isosurf(XYZE_mesh, interp_E, data_relax_list,
                        scene, nb_grid_pts=100, iso_value=None):

    if iso_value is None:
        iso_value = 50
        try:
            input_iso = float(
                int(input("Energy of the isosurface ([integer] in mEV) : ")))
            if input_iso > 0 and input_iso < max(XYZE_mesh[:, 3]):
                iso_value = input_iso
            else:
                print(
                    "invalid energy value for isosurface, set back to default = 50 mEv")
        except Exception as ex:
            print(ex)
            print("invalid energy value for isosurface, set back to default = 50 mEv")

    data_interp = get_isosurface_data(
        interp_E, iso_value, XYZE_mesh, nb_grid_pts)

    fig = tools.make_subplots(rows=1, cols=2,
                              specs=[[{'is_3d': True}, {'is_3d': True}]])
    fig['layout']['scene1'].update(scene)
    fig['layout']['scene2'].update(scene)

    # only raw mesh points on left fig
    data_mesh = {}
    data_mesh['scene'] = 'scene1'
    fig.append_trace(data_mesh, 1, 1)

    # both isosurface and relax points on right fig
    data_interp['scene'] = 'scene2'
    fig.append_trace(data_interp, 1, 2)
    for data_relax in data_relax_list:
        data_relax['scene'] = 'scene2'
        fig.append_trace(data_relax, 1, 2)

    po.plot(fig, filename='4D_PES_isosurf_{:.0f}.html'.format(iso_value))
    #del data_relax_all['scene']
    #del data_relax['scene']
    return(fig)


def plot_slider_isosurf(
    XYZE_mesh,
    interp_E,
    data_relax_list,
    scene,
    nb_grid_pts=100,
    energy_values=[
        10,
        25,
        40,
        50,
        75,
        100]):

    print("let's go for the slider effect !")
    figure = {
        'data': [],
        'layout': {},
        'frames': [],
        'config': {'scrollzoom': True}
    }

    # fill in most of layout
    figure['layout']['scene'] = dict(scene)
    figure['layout']['hovermode'] = 'closest'
    figure['layout']['plot_bgcolor'] = 'rgb(223, 232, 243)'

    sliders_dict = {
        'active': 0,
        'yanchor': 'top',
        'xanchor': 'left',
        'currentvalue': {
            'font': {'size': 20},
            'prefix': 'energy (mEV) : ',
            'visible': True,
            'xanchor': 'right'
        },
        'transition': {'duration': 300, 'easing': 'cubic-in-out'},
        'pad': {'b': 10, 't': 50},
        'len': 0.9,
        'x': 0.1,
        'y': 0,
        'steps': []
    }

    for energy in energy_values:
        frame = {'data': [], 'name': str(energy)}
        data_interp = get_isosurface_data(
            interp_E, energy, XYZE_mesh, nb_grid_pts)
        if len(figure['data']) == 0:
            figure['data'].append(data_interp)
            figure['data'] += data_relax_list
        frame['data'].append(data_interp)
        frame['data'] += data_relax_list
        figure['frames'].append(frame)
        slider_step = {'args': [[energy],
                                {'frame': {'duration': 300, 'redraw': False},
                                 'mode': 'immediate',
                                 'transition': {'duration': 300}}],
                       'label': energy,
                       'method': 'animate'}

        sliders_dict['steps'].append(slider_step)

    figure['layout']['sliders'] = [sliders_dict]
    po.plot(figure, filename='4D_PES_isosurf_slider.html')
    return(figure)


def interpolate_3D_PEV(XYZE, nb_grid_pts):
    # interpolation of the grid of computed points
    grid_x, grid_y, grid_z = np.mgrid[min(XYZE[:, 0]):max(XYZE[:, 0]):(nb_grid_pts *
                                                                       np.complex(0, 1)), min(XYZE[:, 1]):max(XYZE[:, 1]):(nb_grid_pts *
                                                                                                                           np.complex(0, 1)), min(XYZE[:, 2]):max(XYZE[:, 2]):(nb_grid_pts *
                                                                                                                                                                               np.complex(0, 1))]
    #print(np.shape(grid_x),np.shape(grid_y), np.shape(grid_z) )
    # print(grid_x) #,np.shape(grid_y), np.shape(grid_z) )
    interp_E = griddata(XYZE[:, 0:3], XYZE[:, 3], (grid_x, grid_y, grid_z))
    #print("interp E shape : {}".format(interp_E.shape))
    #print("interp E value : {}".format(interp_E))
    # print("interp E  contains NaN values ? \nif yes, they lie at : \n  : {}".
    #      format(np.argwhere(np.isnan(interp_E)) ) )
    return(interp_E)


def get_scatterplot_data(
        XYZE,
        midpoint=100,
        numtick=6,
        secondary_data=False,
        hover_txt=None):
    #midpoint = 100
    ticks = np.concatenate([np.linspace(min(XYZE[:, 3]), midpoint, num=numtick, endpoint=False),
                            np.linspace(midpoint, max(XYZE[:, 3]), num=numtick)])
    if hover_txt is None:
        hover_txt = ["Sz =  {:.3f} meV".format(e) for e in XYZE[:, 3]]
    print("hover txt len : {}".format(len(hover_txt)))
    data = go.Scatter3d(
        x=XYZE[:, 0],
        y=XYZE[:, 1],
        z=XYZE[:, 2],
        mode='markers',
        text=hover_txt,
        marker=dict(
            size=5 if not secondary_data else 2,
            symbol="circle" if not secondary_data else "x",
            cmax=max(np.log10(XYZE[:, 3] + 1)),
            cmin=min(np.log10(XYZE[:, 3] + 1)),
            # set color to an array/list of desired values
            color=np.log10(XYZE[:, 3] + 1),
            colorscale='Jet',   # choose a colorscale
            opacity=0.8,
            colorbar=dict(
                title='Sz',
                titleside='top',
                tickmode='array',
                tickvals=np.log10(ticks + 1),
                ticktext=["{:.0f}".format(t) for t in ticks],
                ticks='outside'
            ),
        )
    )
    return(data)


def get_isosurface_data(interp_E, iso_value, XYZE, nb_grid_pts):

    spacing = tuple([(max(XYZE[:, 0]) - min(XYZE[:, 0])) / nb_grid_pts,
                     (max(XYZE[:, 1]) - min(XYZE[:, 1])) / nb_grid_pts,
                     (max(XYZE[:, 2]) - min(XYZE[:, 2])) / nb_grid_pts])
    print(spacing)

    vertices, faces, normals, values = measure.marching_cubes(
        interp_E, iso_value, spacing=spacing, step_size=1)
    print("marching cubes done")
    print(type(vertices))
    print(vertices)
    x, y, z = [vertices[:, i] + min(XYZE[:, i]) for i in range(3)]
    I, J, K = zip(*faces)
    for name, vect in zip(("x", "y", "z"), (x, y, z)):
        if np.isnan(np.dot(vect, vect)):
            print("{} contains NaN values\nPositions : {}".format(
                name, np.argwhere(np.isnan(vect))))

    triangles = go.Mesh3d(x=x,
                          y=y,
                          z=z,
                          i=I,
                          j=J,
                          k=K,
                          name='',
                          # alphahull=5,
                          opacity=0.4,
                          color='#00FFFF'
                          )

    return(triangles)


def interpolate_PES(computed_mesh):

    XYZ = np.array([[d["disto_prec"], d["x_na"], d['eform']] for d in computed_mesh
                    if d["status"] >= 3 and d["disto_prec"] > 0.39])

    # interpolation of the grid of computed points
    y = np.unique(XYZ[:, 1])
    print(y)
    X, Y = [np.reshape(XYZ[:, i], (len(y), -1)) for i in [0, 1]]
    print(X, Y)
    # x,y = [ np.unique(XYZ[:,i]) for i in [0,1] ]
    # print(x , y )
    # X, Y = np.meshgrid(x, y)
    z = XYZ[:, 2]
    Z = np.array(z).reshape(len(y), -1)
    print(Z)
    points = XYZ[:, 0:2]
    print(np.shape(X), np.shape(Y), np.shape(Z))
    grid_x, grid_y = np.mgrid[min(X[0, :]):max(
        X[0, :]):40j, min(Y[:, 0]):max(Y[:, 0]):80j]
    print(np.shape(grid_x), np.shape(grid_y))
    interp_Z = griddata(points, z, (grid_x, grid_y), method='cubic')
    # print(interp_Z)
    return(grid_x, grid_y, interp_Z)


def plot_PES_3D(struct_list, distortion_mesh):

    # input : stacking_dict[stacking] = {"run_list" : struct_list , "mesh" : sorted_distortion_mesh }

    stacking = struct_list[0]["stacking"]

    fig = plt.figure(stacking + " Energy surface")
    ax = Axes3D(fig)

    # plot the relaxed runs
    for filter_lvl, color, legend, size in [
            (5, "red", "relaxed : on hull", 200), (4, "orange", "relaxed : off-hull", 100)]:

        XYZ = np.array([[d["disto"], d["x_na"], d['eform']]
                        for d in struct_list if d["status"] == filter_lvl])
        ax.scatter(XYZ[:, 0], XYZ[:, 1], XYZ[:, 2], zdir='z',
                   c=color, marker="*", label=legend, s=size)

    # plot the grid of distorted Single Point runs
    for filter_lvl, color, legend, size in [(5, "red", "mesh : on hull", 200),
                                            (4, "blue", "mesh : off-hull ", 100),
                                            (3, "blue", "Calculated mesh point", 10)]:
        XYZ = np.array([[d["disto"], d["x_na"], d['eform']]
                        for d in distortion_mesh if d["status"] == filter_lvl])
        ax.scatter(XYZ[:, 0], XYZ[:, 1], XYZ[:, 2], zdir='z',
                   c=color, marker="o", label=legend, s=size)

    # interpolation of the grid of points
    grid_x, grid_y, interp_Z = interpolate_PES(distortion_mesh)
    surf_interp = ax.plot_surface(grid_x, grid_y, interp_Z, alpha=0.2)
    wirefram_interp = ax.plot_wireframe(
        grid_x,
        grid_y,
        interp_Z,
        rstride=3,
        cstride=3,
        linewidth=0.4,
        label="interpolation")

    ax.set_xlabel('Mn dist from octa')
    ax.set_ylabel('$X_{Na}$')
    ax.set_zlabel('E (meV)')
    fig.legend()
    plt.show(block=False)
    # plot the minima of energy at each Na (calculated)
    # ie : the minimal path "by hand"


def plot_PES_projected(struct_list, distortion_mesh):
    stacking = struct_list[0]["stacking"]

    # input : stacking_dict[stacking] = {"run_list" : struct_list ,
    # "mesh" : sorted_distortion_mesh }

    grid_x, grid_y, interp_Z = interpolate_PES(distortion_mesh)
    print(np.shape(grid_x), np.shape(grid_y), np.shape(interp_Z))

    Y = grid_y[0, :]  # x_na
    x = grid_x[:, 0]  # distortion
    print(x, Y)
    # for j in range(len(Y)) :
    #     print(interp_Z[:,j])
    # plot a valley of energy (by scaling the bottom of the valley to zero )
    rescaled_Z = np.array([interp_Z[:, j] - min(interp_Z[:, j])
                           for j in range(len(Y))])

    # plot the path at the bottom of the valley (i.e. the X of minimal Z for
    # each Y )
    argX = [np.argmin(interp_Z[:, j]) for j in range(len(Y))]
    X = grid_x[argX, 0]

    print("x boundaries {} {} ".format(min(x), max(x)))
    print("Y boundaries {} {} ".format(min(Y), max(Y)))
    fig = plt.figure(stacking + " Energy valley rescaled projected ")
    ax = fig.add_subplot(111)

    heatmap = plt.imshow(rescaled_Z[2:-1,
                                    :],
                         extent=(min(x) - 0.1,
                                 max(x),
                                 min(Y),
                                 max(Y)),
                         origin='lower',
                         cmap=plt.get_cmap("jet"))
    # , cmap=cmap, norm=norm, boundaries=bounds, ticks=[0, 5, 10])
    plt.colorbar(heatmap)

    ax.plot(X, Y, color="red", label=stacking + " minimum energy path")
    # plot the relaxed runs
    for filter_lvl, color, legend, size in [
            (5, "red", "relaxed : on hull", 100), (4, "orangered", "relaxed : off-hull", 50)]:

        XY = np.array([[d["disto"], d["x_na"]]
                       for d in struct_list if d["status"] == filter_lvl])
        ax.scatter(XY[:, 0], XY[:, 1], c=color,
                   marker="*", label=legend, s=size)

    ax.set_xlabel('Mn dist from octa')
    ax.set_ylabel('$X_{Na}$')
    ax.set_ylim(0, 1)
    #plt.gcf().set_size_inches(6, 6)
    fig.legend()
    plt.show(block=False)

    if input("close all graphs now ?? Y/N ") == "Y":
        plt.close("all")
