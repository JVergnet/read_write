# Read_Run_pages.py

import copy
import itertools
# import os
import tkinter as tk
from functools import partial
from tkinter import filedialog, ttk  # , messagebox

import matplotlib
# import numpy as np
# from matplotlib import pyplot as plt
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)
# from matplotlib.figure import Figure

from electronic_analysis import read_dos as DOS
from electronic_analysis import read_hull as hull
from electronic_analysis import read_bader as bader
# import readRun_entries as readRun
from utils import generic_plot
from filtering_runs import filter_runs

matplotlib.use("TkAgg")


LARGE_FONT = ("Verdana", 12)


def analysis_pages():
    return([DOS_page, Hull_page, Bader_page, COOP_page, O2_page])


class Blank_page(tk.Frame):
    "base class for every page"

    def __init__(self, tab_control, controller):
        # a Frame is contained in its tab_control
        ttk.Frame.__init__(self, tab_control)
        self.tab_control = tab_control
        print(self.graph_type)
        self.controller = controller
        frame_title = tk.Frame(self)
        frame_title.pack(fill=tk.BOTH, expand=True)
        label = tk.Label(
            frame_title,
            text="{} page".format(
                self.graph_type),
            font=LARGE_FONT)
        label.pack(fill=tk.BOTH, expand=True)  # print "North"
        #label.pack(pady = 10, padx = 10)

    # def create_buttons(self, frame_names):
    #     self.buttons = {}
    #     self.valid_frame_names = [n for n in frame_names if self.graph_type != n ]
    #     for i,page_name in enumerate(self.valid_frame_names):
    #         action = partial(self.controller.show_frame, page_name)
    #         self.buttons[page_name] = ttk.Button(self, text = page_name, command=action)
    # self.buttons[page_name].grid(row = 6 , column = i) #, sticky = tk.S)


class Plot_page(Blank_page):
    "base class for plot & analysis pages"

    def __init__(self, tab_control, controller, filtered_run_list):
        Blank_page.__init__(self, tab_control, controller)
        self.figure_widget = None
        self.toolbar = None
        self.filtered_run_list = filtered_run_list

    def display_plot(self, fig, frame):
        # plot
        # try:
        #     self.figure_widget.destroy()
        #     self.toolbar.destroy()
        # except :
        #     pass
        self.canvas = FigureCanvasTkAgg(fig, frame)
        self.canvas.draw()
        self.figure_widget = self.canvas.get_tk_widget()
        self.figure_widget.pack(
            fill=tk.BOTH,
            expand=1,
            padx=5,
            pady=5,
            side=tk.TOP)

        self.toolbar = NavigationToolbar2Tk(self.canvas, frame)
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(
            fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        def on_key_press(event):
            print("you pressed {}".format(event.key))
            key_press_handler(event, self.canvas, self.toolbar)

        self.canvas.mpl_connect("key_press_event", on_key_press)

    def save_fig(self, fig):
        try:
            from pathlib import Path
            home = str(Path.home())
            print("pathlib : fome Path : {}".format(home))
        except Exception as ex:
            print(ex)
            home = "/home/jvergnet/"

        self.fig_file = filedialog.asksaveasfile(
            initialdir=home, title="Select File", filetypes=[
                ('PNG', ".png"), ("SVG", ".svg"), ("PDF", ".pdf")])
        try:
            print(self.fig_file, self.fig_file.name)
            fig.savefig("{}".format(self.fig_file.name), bbox_inches='tight')
        except Exception as ex:
            print(ex)


class DOS_page(Plot_page):

    graph_type = "DOS"

    def __init__(self, tab_control, controller, filtered_run_list):
        Plot_page.__init__(self, tab_control, controller, filtered_run_list)
        # self.figure_widget = None
        # self.toolbar = None

        def clicked_button_select():
            self.hull_spin = {
                "up & down": 0, "sum": 1, "diff": 2}[
                combo_spin.get()]
            self.hull_proj = combo_proj.get()
            self.Erange = [emin_var.get(), emax_var.get()]
            # choice_dict = {"a" : 1, "c" : 3 , "m" : 4 , "h" : 5 }
            # self.sieve_lvl = choice_dict[self.hull_proj[0]]
            text_select.config(
                text="Selected parameters : \nSpin : {} \nProj {}\nEnergy range :[{} , {} ]".format(
                    self.hull_spin, self.hull_proj, self.Erange[0], self.Erange[1]))
            #text_result.config(text="Initial runs : {}".format(len(filtered_run_list)))

        def clicked_button_PLOT():
            print(
                filtered_run_list,
                self.hull_spin,
                self.hull_proj,
                self.Erange)
            try:
                self.fig.clf()
            except BaseException:
                pass
            self.fig = DOS.plot_DOS(
                filtered_run_list,
                spin_choice=self.hull_spin,
                DOS_choice=self.hull_proj,
                Erange=self.Erange)
            # plt.show(block=False)
            self.display_plot(self.fig, frame_bottom)
            text_result.config(
                text="PLOTTED runs : {} \n spin {} proj {}\nEnergy range :[{} , {} ]".format(
                    len(filtered_run_list),
                    self.hull_spin,
                    self.hull_proj,
                    self.Erange[0],
                    self.Erange[1]))

        def clicked_button_SAVE():
            self.save_fig(self.fig)

        def clicked_button_help():
            help_txt = """
HELP ON DOS !!
"""
            tk.messagebox.showinfo('Proj help', help_txt)

        # FRAME TOP  =========================================================
        frame_top = tk.Frame(self)
        frame_top.pack(fill=tk.BOTH, expand=True)

        # === SELECTION : spin value
        frame_spin = tk.Frame(frame_top)
        frame_spin.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)

        text_spin = tk.Label(frame_spin, text="Spin selection")
        text_spin.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        combo_spin = ttk.Combobox(frame_spin)
        combo_spin['values'] = ("up & down", "sum", "diff")
        combo_spin.current(0)
        combo_spin.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        self.hull_spin = combo_spin.get()

        # === SELECTION : PROJECTION value
        frame_proj = tk.Frame(frame_top)
        frame_proj.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)

        text_proj = tk.Label(frame_proj, text="Projection selection ?")
        text_proj.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        combo_proj = ttk.Combobox(frame_proj)
        combo_proj['values'] = ("spd", "perSite", "perElt", "total")
        combo_proj.current(2)
        combo_proj.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        self.hull_proj = combo_proj.get()

        # === SELECTION : ENERGY RANGE OF DOS
        frame_erange = tk.Frame(frame_top)
        frame_erange.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)

        # == emin frame
        frame_emin = tk.Frame(frame_erange)
        frame_emin.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)

        emin_var = tk.IntVar()
        emin_var.set(-3)
        text_emin = tk.Label(frame_emin, text="Emin ?")
        text_emin.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)
        comboroll_emin = tk.Spinbox(
            frame_emin,
            from_=-10,
            to=0,
            textvariable=emin_var,
            width=5)
        comboroll_emin.pack(
            fill=tk.BOTH,
            expand=1,
            padx=5,
            pady=5,
            side=tk.TOP)

        # == emin frame
        frame_emax = tk.Frame(frame_erange)
        frame_emax.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)

        emax_var = tk.IntVar()
        emax_var.set(3)
        text_emax = tk.Label(frame_emax, text="Emax ?")
        text_emax.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)
        comboroll_emax = tk.Spinbox(
            frame_emax,
            from_=0,
            to=10,
            textvariable=emax_var,
            width=5)
        comboroll_emax.pack(
            fill=tk.BOTH,
            expand=1,
            padx=5,
            pady=5,
            side=tk.TOP)

        # FRAME MIDDLE  =======================================================
        frame_mid = tk.Frame(self)
        frame_mid.pack(fill=tk.BOTH, expand=True)

        # === SELECT
        frame_select = tk.Frame(frame_mid)
        frame_select.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)

        button_select = ttk.Button(frame_select, text="SELECT parameters",
                                   command=clicked_button_select)
        button_select.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        text_select = tk.Label(frame_select, text="DOS parameters")
        text_select.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)
        text_select.config(
            text="Selected parameters : \n spin {} proj {}".format(
                self.hull_spin, self.hull_proj))

        # === APPLY
        frame_apply = tk.Frame(frame_mid)
        frame_apply.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)

        button_PLOT = ttk.Button(frame_apply, text="PLOT",
                                 command=lambda: clicked_button_PLOT())
        button_PLOT.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        text_result = tk.Label(
            frame_apply,
            text="Filtered runs : {}".format(
                len(filtered_run_list)))
        text_result.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        # === HELP
        frame_help_save = tk.Frame(frame_mid)
        frame_help_save.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)
        button_help = ttk.Button(frame_help_save, text="HELP",
                                 command=clicked_button_help)
        button_help.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        button_SAVE = ttk.Button(frame_help_save, text="SAVE",
                                 command=lambda: self.save_fig(self.fig))
        button_SAVE.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        # FRAME BOTTOM ========================================================
        frame_bottom = tk.Frame(self)
        frame_bottom.pack(fill=tk.BOTH, expand=True)


class Hull_page(Plot_page):

    graph_type = "Hull"

    def __init__(self, tab_control, controller, filtered_run_list):
        Plot_page.__init__(self, tab_control, controller, filtered_run_list)
        # if len( [d for d in filtered_run_list if d["status"] >= 3 ]) >= 1 :
        self.setup_hull_page(tab_control, controller, filtered_run_list)

    def setup_hull_page(self, tab_control, controller, filtered_run_list):

        def clicked_button_plot(hull_plot_type):

            try:
                fig = self.figures[hull_plot_type]
                fig.clf()
            except KeyError as ex:
                print(ex)

            if hull_plot_type == "simple":
                fig = generic_plot.plot_structure_value_evolution(
                    filtered_run_list, prop_list=['ediff'], legend=None)

            elif hull_plot_type == "hull":
                coord = "x_na"
                fig = hull.plot_convex_hull(filtered_run_list, coord=coord)

            elif hull_plot_type == "voltage":
                hull_entries = copy.deepcopy(filtered_run_list)
                hull_entries = filter_runs.generate_hull_tags(hull_entries)
                hull_entries = [s for s in hull_entries if s["status"] == 5]
                fig = hull.plot_voltage_curve(hull_entries)

            self.figures[hull_plot_type] = fig

            self.display_plot(fig, self.frame_fig[plot_type])

        def clicked_button_save(hull_plot_type):
            self.save_fig(self.figures[hull_plot_type])

        def clicked_button_help():
            help_txt = """
HELP ON HULL !!
"""
            tk.messagebox.showinfo('Proj help', help_txt)

        # FRAME TOP  =========================================================
        frame_top = tk.Frame(self)
        frame_top.pack(fill=tk.BOTH, expand=True)

        # === SELECTION : X coord
        frame_X_coord = tk.Frame(frame_top)
        frame_X_coord.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)

        text_X_coord = tk.Label(frame_X_coord, text="X coord selection")
        text_X_coord.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        combo_X_coord = ttk.Combobox(frame_X_coord)
        combo_X_coord['values'] = ("ediff", "dOO_min")
        combo_X_coord.current(0)
        combo_X_coord.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        self.coord = combo_X_coord.get()

        # === SELECTION : x_na range
        frame_x_range = tk.Frame(frame_top)
        frame_x_range.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)

        # == x_min frame
        frame_x_min = tk.Frame(frame_x_range)
        frame_x_min.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)

        x_min_var = tk.IntVar()
        x_min_var.set(-3)
        text_x_min = tk.Label(frame_x_min, text="X_Min ?")
        text_x_min.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)
        comboroll_x_min = tk.Spinbox(
            frame_x_min,
            from_=-10,
            to=0,
            textvariable=x_min_var,
            width=5)
        comboroll_x_min.pack(
            fill=tk.BOTH,
            expand=1,
            padx=5,
            pady=5,
            side=tk.TOP)

        # == x_max frame
        frame_x_max = tk.Frame(frame_x_range)
        frame_x_max.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)

        x_max_var = tk.IntVar()
        x_max_var.set(3)
        text_x_max = tk.Label(frame_x_max, text="X_Max ?")
        text_x_max.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)
        comboroll_x_max = tk.Spinbox(
            frame_x_max,
            from_=0,
            to=10,
            textvariable=x_max_var,
            width=5)
        comboroll_x_max.pack(
            fill=tk.BOTH,
            expand=1,
            padx=5,
            pady=5,
            side=tk.TOP)

        # === HELP
        frame_help_save = tk.Frame(frame_top)
        frame_help_save.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)
        button_help = ttk.Button(frame_help_save, text="HELP",
                                 command=clicked_button_help)
        button_help.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        # FRAME MIDDLE  =======================================================
        frame_mid = tk.Frame(self)
        frame_mid.pack(fill=tk.BOTH, expand=True)

        # FRAME BOTTOM ========================================================
        frame_bottom = tk.Frame(self)
        frame_bottom.pack(fill=tk.BOTH, expand=True)

        frame_button = {}
        button_plot = {}
        button_save = {}
        self.frame_fig = {}
        self.figures = {}

        for plot_type in ["simple", "hull", "voltage"]:

            frame_button[plot_type] = tk.Frame(frame_mid)
            frame_button[plot_type].pack(
                fill=tk.BOTH, expand=True, side=tk.LEFT)

            plot_action = partial(clicked_button_plot, plot_type)
            button_plot[plot_type] = ttk.Button(
                frame_button[plot_type],
                text="PLOT {}".format(plot_type),
                command=plot_action)
            button_plot[plot_type].pack(
                fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

            save_action = partial(clicked_button_save, plot_type)
            button_save[plot_type] = ttk.Button(
                frame_button[plot_type],
                text="SAVE {}".format(plot_type),
                command=save_action)
            button_save[plot_type].pack(
                fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

            self.frame_fig[plot_type] = tk.Frame(frame_bottom)
            self.frame_fig[plot_type].pack(
                fill=tk.BOTH, expand=True, side=tk.LEFT)


class Bader_page(Plot_page):

    graph_type = "Bader"

    def __init__(self, tab_control, controller, filtered_run_list):
        Plot_page.__init__(self, tab_control, controller, filtered_run_list)

        def clicked_button_do_bader():

            self.filtered_run_list = bader.get_bader_tags(
                self.filtered_run_list)

        def clicked_button_plot():
            selected_qty = [q for q in quantities if qty_check[q].get() == 1]
            selected_elt = [e for e in elements if elt_check[e].get() == 1]

            self.frame_fig = {}
            self.figures = {}
            self.frame_rows = {}

            frame_elt_label = tk.Frame(frame_bottom)
            frame_elt_label.pack(fill=tk.BOTH, expand=True, side=tk.TOP)

            for j, elt in enumerate(selected_elt):
                text_elt_label = tk.Label(frame_elt_label, text=elt)
                text_elt_label.pack(
                    fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.LEFT)

            for i, qty in enumerate(selected_qty):
                self.frame_fig[qty] = {}
                self.figures[qty] = {}

                self.frame_rows[qty] = tk.Frame(frame_bottom)
                self.frame_rows[qty].pack(
                    fill=tk.BOTH, expand=True, side=tk.TOP)

                text_qty_label = tk.Label(self.frame_rows[qty], text=qty)
                text_qty_label.pack(
                    fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.LEFT)

                for j, elt in enumerate(selected_elt):

                    frame = tk.Frame(self.frame_rows[qty])
                    frame.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)

                    # self.figures["{}{}".format(elt,qty)] =  ttk.Button(self.frame_rows[qty], text = "PLOT{}".format(elt))
                    # self.figures["{}{}".format(elt,qty)].pack(fill=tk.BOTH, expand=True,side=tk.LEFT)

                    try:
                        self.figures[qty][elt].clf()
                        print("figure cleaned : {} of {}".format(elt, qty))
                    except BaseException:
                        pass

                    fig = generic_plot.plot_site_value_evolution(
                        self.filtered_run_list,
                        elt,
                        value=qty,
                        coord="x_na",
                        plot_type=[0, 1])

                    self.display_plot(fig, frame)

                    self.frame_fig[qty][elt] = frame
                    self.figures[qty][elt] = fig

                    print("{}{}".format(elt, qty))
                    input("click me !")

        def clicked_button_SAVE():
            for f in self.figures():
                self.save_fig(self.fig)

        # FRAME TOP  =  SELECTION : bader quantity ===========================
        frame_top = tk.Frame(self)
        frame_top.pack(fill=tk.BOTH, expand=True)

        qty_check = {}
        quantities = ["charge", "magnetization", 'vol_chg', 'vol_mag']
        for qty in quantities:
            qty_check[qty] = tk.IntVar()
            tk.Checkbutton(
                frame_top,
                text=qty,
                variable=qty_check[qty]).pack(
                fill=tk.Y,
                expand=True,
                side=tk.LEFT)

        # FRAME MID  = SELECTION : ELEMENTS  =============
        frame_middle = tk.Frame(self)
        frame_middle.pack(fill=tk.BOTH, expand=True)

        elt_check = {}
        irreg_list_of_elt = [[str for str in s["structure"].composition.get_el_amt_dict(
        ).keys()] for s in filtered_run_list]
        elements = set(itertools.chain(*irreg_list_of_elt))
        for elt in elements:
            elt_check[elt] = tk.IntVar()
            tk.Checkbutton(
                frame_middle,
                text=elt,
                variable=elt_check[elt]).pack(
                fill=tk.Y,
                expand=True,
                side=tk.LEFT)
        # PLOT BUTTON ================
        frame_mid = tk.Frame(self)
        frame_mid.pack(fill=tk.BOTH, expand=True)

        button_do_bader = ttk.Button(
            frame_mid,
            text="get bader tag",
            command=clicked_button_do_bader)
        button_do_bader.pack(fill=tk.BOTH, expand=True, side=tk.TOP)

        button_plot = ttk.Button(
            frame_mid,
            text="PLOT",
            command=clicked_button_plot)
        button_plot.pack(fill=tk.BOTH, expand=True, side=tk.TOP)
        # FRAME BOTTOM ========================================================
        frame_bottom = tk.Frame(self)
        frame_bottom.pack(fill=tk.BOTH, expand=True)


#        plot_action = partial(clicked_button_plot, plot_type)


class COOP_page(Plot_page):

    graph_type = "COOP"

    def __init__(self, tab_control, controller, filtered_run_list):
        Plot_page.__init__(self, tab_control, controller, filtered_run_list)


class O2_page(Plot_page):

    graph_type = "O2 release"

    def __init__(self, tab_control, controller, filtered_run_list):
        Plot_page.__init__(self, tab_control, controller, filtered_run_list)
