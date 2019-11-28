#!/usr/bin/env python3
# Read_Run_app.py

import importlib
import os
import tkinter as tk
import traceback
# Imports
import types
from functools import partial
from tkinter import filedialog, messagebox, ttk

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)
from matplotlib.figure import Figure

import electronic_analysis.read_hull as hull
import electronic_analysis.readRun_entries as readRun
import Read_Run_pages as pages

matplotlib.use("TkAgg")


# import read_dos as DOS


# Initialisation of the different parameters and lists
# V = [-1.,0.01,1.]
# I = [3.,4.,5.]

# n = 4. #Nb of electrons exchanged
# F = 96485 #C/mol, Faraday's constant
# R = 8.314 #J.k-1.mol-1, perfect gaz constant
# T = 298.15 #Operating temperature
# f = F/(R*T) #Reduced Faraday's constant
# alpha = 0.5 #Electron transfer coefficient
# A = 1 #Area of the Electrode

LARGE_FONT = ("Verdana", 12)

valid_run_list = False
initial_run_list = []
filtered_run_list = []

# Functions


def get_bool(prompt):
    while True:
        answer = ""
        try:
            answer = input("{} [Y/N] :".format(prompt))
            if answer[0] in ["Y", "y", "T", "t"]:
                return True
            elif answer[0] in ["N", "n", "F", "f"]:
                return False
        except KeyError:
            pass
        print("Invalid input please enter T,t,Y,y or F,f,N,n !")


def Export():
    "Creates a plotly plot on the internet"
    #plotly.tools.set_credentials_file(username='TanguyP', api_key='Cy0FkAUAjVjlMVSrqzrh')
    print("File exported")


def Help():
    "Call for help"
    tk.messagebox.showinfo(
        'Looking for help ?',
        "Ask God or Jean-Marie, they'll know")


def plot():
    "Plot something"
    print("wow it plots ! ")

# Lecture de dossier

# a = get_bool("Do you want to read runs   ?")

# if True :
#     Name = input("What is the name of the file ?")
#     Path = Name
#     print(Path)

# with open("LSV Screening @64°C_C01.txt", "r") as f:
 #   for line in f:
  #      myLine = line.strip()
   #     Pot, PowPot, Cur, PowCur = [float(elt) for elt in myLine.split("E"," ","E")]
    #    Potentiel.append(Pot), Current.append(Cur)
    #   print(Pot)

#print(Pot, Current)


# import importlib


def reload_package(package):
    assert(hasattr(package, "__package__"))
    fn = package.__file__
    fn_dir = os.path.dirname(fn) + os.sep
    module_visit = {fn}
    del fn

    def reload_recursive_ex(module):
        importlib.reload(module)

        for module_child in vars(module).values():
            if isinstance(module_child, types.ModuleType):
                fn_child = getattr(module_child, "__file__", None)
                if (fn_child is not None) and fn_child.startswith(fn_dir):
                    if fn_child not in module_visit:
                        # print("reloading:", fn_child, "from", module)
                        module_visit.add(fn_child)
                        reload_recursive_ex(module_child)

    return reload_recursive_ex(package)

# Interface

# Initialisation


class ReadRun_app(tk.Tk):

    def __init__(self, *args, **kwargs):

        tk.Tk.__init__(self, *args, **kwargs)

        # tk.Tk.iconbitmap(self, default = "CDF.ico")
        tk.Tk.wm_title(self, "Vasp read run")

        # the container is a Frame that holds other frames
        self.tab_control = ttk.Notebook(self)
        self.tab_control.pack(side="top", fill="both", expand=True)
        self.tab_control.grid_rowconfigure(0, weight=1)
        self.tab_control.grid_columnconfigure(0, weight=1)

        self.frames = {}
        self.start_frames = {}
        self.plot_frames = {}

        frame = Start_page(self.tab_control, self)
        name = frame.graph_type
        self.frames[name] = frame
        self.start_frames[name] = frame
        frame.grid(row=0, column=0, sticky="nsew")
        self.tab_control.add(frame, text=name)

        self.show_frame("start")  # Start_page)

    def show_frame(self, name):
        frame = self.frames[name]
        frame.tkraise()

    def create_plot_frames(self, filtered_run_list):
        for F in pages.analysis_pages():
            frame = F(self.tab_control, self, filtered_run_list)
            frame.pack(fill=tk.BOTH, expand=True)
            name = frame.graph_type
            self.frames[name] = frame
            self.plot_frames[name] = frame
            self.tab_control.add(frame, text=name)

        self.show_frame("DOS")

    def destroy_plot_frames(self):
        for f in self.plot_frames.values():
            self.tab_control.forget(f)
            f.pack_forget()
            f.destroy()
        print("pages destroyed")

    def reload_plot_frames(self, filtered_run_list):
        try:
            reload_package(pages)
            print("pages lib reloaded")
        except Exception as ex:
            print(traceback.format_exc())
            print(ex)

        self.destroy_plot_frames()

        self.create_plot_frames(filtered_run_list)
        print("pages created")

# class Blank_page(tk.Frame):

#     def __init__(self, tab_control, controller, **kwargs):
#         ttk.Frame.__init__(self,tab_control) # a Frame is contained in its tab_control
#         self.tab_control = tab_control
#         print(self.graph_type)
#         self.controller = controller
#         frame_title = tk.Frame(self)
#         frame_title.pack(fill=tk.BOTH, expand=True)
#         label = tk.Label(frame_title, text = "{} page".format(self.graph_type), font = LARGE_FONT)
#         label.pack(fill=tk.BOTH, expand=True) #print "North"
#         #label.pack(pady = 10, padx = 10)

#     # def create_buttons(self, frame_names):
#     #     self.buttons = {}
#     #     self.valid_frame_names = [n for n in frame_names if self.graph_type != n ]
#     #     for i,page_name in enumerate(self.valid_frame_names):
#     #         action = partial(self.controller.show_frame, page_name)
#     #         self.buttons[page_name] = ttk.Button(self, text = page_name, command=action)
#     #         self.buttons[page_name].grid(row = 6 , column = i) #, sticky = tk.S)


class Start_page(pages.Blank_page):

    graph_type = "start"

    def __init__(self, tab_control, controller, **kwargs):
        self.path = ""

        def valid_path_test(path):
            collect_tag_button.config(state="disabled")
            if os.path.exists(path):
                text = "CORRECT PATH :\n {}".format(path)
                collect_run_button.config(state="normal")
                text_collect.config(text="Valid folder, ready to run")
            else:
                text = "INVALID PATH :\n {}".format(path)
                collect_run_button.config(state="disabled")
                text_collect.config(text="Require valid path")
            text_chosed_path.configure(text=text)

        def clicked_path_button():
            self.path = path_entry.get()
            valid_path_test(self.path)

        def clicked_browse_button():
            self.path = filedialog.askdirectory()
            valid_path_test(self.path)

        def clicked_collect_button():
            self.job_type = job_type_var.get()
            print("job type", self.job_type)
            self.vasprundictList = readRun.collect_valid_runs(
                self.path, file_system_choice=self.job_type)
            print(self.vasprundictList)
            global initial_run_list
            initial_run_list = self.vasprundictList
            nb_runs = len(self.vasprundictList)
            text_collect.configure(
                text="Job type : {}\n{} converged runs".format(
                    self.job_type, nb_runs))
            if nb_runs > 0:
                collect_tag_button.config(state="normal")
                text_tag.config(text="GENERATE tag(s)")
            else:
                collect_tag_button.config(state="disabled")
                text_tag.config(text="Not enough valid runs !")

        def clicked_tag_button():
            # readRun.generate_tags(self.vasprundictList)
            self.vasprundictList = None
            global valid_run_list
            valid_run_list = True
            global initial_run_list
            initial_run_list = self.vasprundictList
            text_tag.configure(
                text="Tag DONE !!\n {} runs tagged".format(
                    len(initial_run_list)))
            button_reload.config(state="normal")

            if len(initial_run_list) == 1:
                global filtered_run_list
                filtered_run_list = self.vasprundictList
                try:
                    controller.destroy_plot_frames()
                except Exception as ex:
                    print((ex))

                controller.create_plot_frames(filtered_run_list)

            else:
                frame = Filter_page(tab_control, controller)
                frame.pack(fill=tk.BOTH, expand=True)
                name = frame.graph_type
                controller.frames[name] = frame
                tab_control.add(frame, text=name)
                controller.show_frame("filter")

        pages.Blank_page.__init__(self, tab_control, controller, **kwargs)
        # FRAME TOP =========================================================
        frame_select_folder = tk.Frame(self)
        frame_select_folder.pack(fill=tk.BOTH, expand=True)
        # TEXT : folder for the run
        text_path = tk.Label(
            frame_select_folder,
            text="ENTER path to the run root")
        text_path.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.LEFT)

        self.pathVar = tk.StringVar(
            None, "/home/jvergnet/frofro/honeycomb/Hull_all/Mg/P2_EE/Na17")
        path_entry = tk.Entry(frame_select_folder, textvariable=self.pathVar)
        path_entry.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.LEFT)

        # Select and display info on path
        path_button = ttk.Button(
            frame_select_folder,
            text="SELECT this folder",
            command=clicked_path_button)
        path_button.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.LEFT)

        browse_button = ttk.Button(
            frame_select_folder,
            text="BROWSE a folder",
            command=clicked_browse_button)
        browse_button.pack(
            fill=tk.BOTH,
            expand=1,
            padx=5,
            pady=5,
            side=tk.LEFT)

        # FRAME MID =========================================================
        # frame_selected_folder = tk.Frame(self)
        # frame_selected_folder.pack(fill=tk.BOTH, expand=True)

        text_chosed_path = tk.Label(self, text="No chosen path yet")
        text_chosed_path.pack(
            fill=tk.BOTH,
            expand=1,
            padx=5,
            pady=5,
            side=tk.TOP)

        # FRAME BOTTOM ========================================================
        frame_read_run = tk.Frame(self)
        frame_read_run.pack(fill=tk.BOTH, expand=True)

        # === radio button job / folder / superfolder
        frame_radio = tk.Frame(frame_read_run)
        frame_radio.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        job_type_var = tk.StringVar(None, "p")
        rad_dict = {}
        for i, name, value in [
                (0, "project", "p"), (1, "super project", "s"), (2, "job", "j")]:
            rad_dict[name] = tk.Radiobutton(
                frame_radio, text=name, value=value, variable=job_type_var)
            rad_dict[name].pack(
                fill=tk.BOTH,
                expand=1,
                padx=5,
                pady=5,
                side=tk.TOP)

        # === Read Run
        frame_collect = tk.Frame(frame_read_run)
        frame_collect.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        # BUTTON : Collect run
        collect_run_button = ttk.Button(
            frame_collect,
            text="READ run in this folder",
            command=clicked_collect_button,
            state="disabled")
        collect_run_button.pack(
            fill=tk.BOTH,
            expand=1,
            padx=5,
            pady=5,
            side=tk.TOP)

        text_collect = tk.Label(frame_collect, text="Require valid path")
        text_collect.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        # === Tag Run
        frame_tag = tk.Frame(frame_read_run)
        frame_tag.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # BUTTON : generate Tag
        collect_tag_button = ttk.Button(
            frame_tag,
            text="TAG the collected run(s)",
            command=clicked_tag_button,
            state="disabled")
        collect_tag_button.pack(
            fill=tk.BOTH,
            expand=1,
            padx=5,
            pady=5,
            side=tk.TOP)
        text_tag = tk.Label(frame_tag, text="Require at least 1 valid run")
        text_tag.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        button_reload = ttk.Button(
            self,
            text="RELOAD",
            command=lambda: controller.reload_plot_frames(filtered_run_list),
            state="disabled")
        button_reload.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)


class Filter_page(pages.Blank_page):

    graph_type = "filter"

    def __init__(self, tab_control, controller, **kwargs):
        pages.Blank_page.__init__(self, tab_control, controller, **kwargs)

        def clicked_button_select():
            self.hull_value = combo_value.get()
            self.hull_filter = combo_filter.get()
            choice_dict = {"a": 1, "c": 3, "m": 4, "h": 5}
            self.sieve_lvl = choice_dict[self.hull_filter[0]]
            text_apply.config(
                text="Selected parameters : \n value {} filter {}".format(
                    self.hull_value, self.hull_filter))
            text_result.config(
                text="Initial runs : {}".format(
                    len(initial_run_list)))

        def clicked_button_apply():
            global filtered_run_list
            restricted_stack_runs = hull.generate_hull_entries(
                initial_run_list, remove_extremes=None, coord=self.hull_value)
            restricted_hull_runs = [
                d for d in restricted_stack_runs if d["status"] >= self.sieve_lvl]
            text_result.config(
                text="Initial runs : {} \n Valid runs left {}".format(
                    len(initial_run_list),
                    len(restricted_hull_runs)))
            if len(restricted_hull_runs) != 0:
                filtered_run_list = restricted_hull_runs
                try:
                    controller.destroy_plot_frames()
                except Exception as ex:
                    print((ex))
                controller.create_plot_frames(filtered_run_list)

        def clicked_button_help():
            help_txt = """
[a]ll : all runs that have at least a POSCAR
[c]onverged : all converged vaspruns
[m]inima : structures of lowest energy for each x_na
[h]ull : structures on the convex hull
"""
            tk.messagebox.showinfo('Filter help', help_txt)
        # FRAME TOP  =========================================================
        frame_top = tk.Frame(self)
        frame_top.pack(fill=tk.BOTH, expand=True)

        # === SELECTION : hull value
        frame_select = tk.Frame(frame_top)
        frame_select.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)

        text_value = tk.Label(frame_select, text="Variable to minimize")
        text_value.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        combo_value = ttk.Combobox(frame_select)
        combo_value['values'] = ("x_na", "min OO distance")
        combo_value.current(0)
        combo_value.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)
        self.hull_value = combo_value.get()

        # === SELECTION : hull filter
        frame_filter = tk.Frame(frame_top)
        frame_filter.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)

        text_filter = tk.Label(frame_filter, text="Type of filter")
        text_filter.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        combo_filter = ttk.Combobox(frame_filter)
        combo_filter['values'] = ("all", "converged", "minima", "hull")
        combo_filter.current(1)
        combo_filter.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        self.hull_filter = combo_filter.get()
        # === APPLY & HELP
        frame_choice = tk.Frame(frame_top)
        frame_choice.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)

        button_help = ttk.Button(frame_choice, text="HELP",
                                 command=clicked_button_help)
        button_help.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        button_select = ttk.Button(
            frame_choice,
            text="SELECT filter parameters",
            command=clicked_button_select)
        button_select.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.TOP)

        # FRAME TOP  =========================================================
        frame_bottom = tk.Frame(self)
        frame_bottom.pack(fill=tk.BOTH, expand=True)

        text_apply = tk.Label(
            frame_bottom,
            text="Filter parameters :\n  value : N.A. \n param : N.A. ")
        text_apply.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.LEFT)

        button_apply = ttk.Button(frame_bottom, text="APPLY filter",
                                  command=clicked_button_apply)
        button_apply.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.LEFT)

        text_result = tk.Label(
            frame_bottom,
            text="initial runs : {} \n Filter results : N.A.".format(
                len(initial_run_list)))
        text_result.pack(fill=tk.BOTH, expand=1, padx=5, pady=5, side=tk.LEFT)


# class DOS_page(Blank_page):

#     graph_type = "DOS"

#     def __init__(self, tab_control, controller, **kwargs):
#         Blank_page.__init__(self,tab_control, controller, **kwargs)
#         self.figure_widget = None
#         self.toolbar = None

#         def clicked_button_select():
#             self.hull_spin ={"up & down":0,"sum":1, "diff":2}[combo_spin.get()]
#             self.hull_proj = combo_proj.get()
#             self.Erange = [ emin_var.get(), emax_var.get()]
#             # choice_dict = {"a" : 1, "c" : 3 , "m" : 4 , "h" : 5 }
#             # self.sieve_lvl = choice_dict[self.hull_proj[0]]
#             text_select.config(text="Selected parameters : \nSpin : {} \nProj {}\nEnergy range :[{} , {} ]".format(
#                 self.hull_spin, self.hull_proj, self.Erange[0], self.Erange[1]) )
#             #text_result.config(text="Initial runs : {}".format(len(filtered_run_list)))

#         def clicked_button_apply() :

#             # plot
#             try:
#                 self.fig.clf()
#                 self.figure_widget.destroy()
#                 self.toolbar.destroy()
#             except :
#                 pass

#             self.fig = DOS.plot_DOS(filtered_run_list, spin_choice = self.hull_spin,  DOS_choice=self.hull_proj, Erange=self.Erange)

#             self.canvas = FigureCanvasTkAgg(self.fig, self)
#             self.canvas.draw()
#             self.figure_widget = self.canvas.get_tk_widget()
#             self.figure_widget.pack(fill=tk.BOTH, expand=1,padx=5, pady=5,side=tk.TOP)

#             self.toolbar = NavigationToolbar2Tk(self.canvas, self)
#             self.toolbar.update()
#             self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=1,padx=5, pady=5,side=tk.TOP)

#             def on_key_press(event):
#                 print("you pressed {}".format(event.key))
#                 key_press_handler(event, self.canvas, self.toolbar)


#             self.canvas.mpl_connect("key_press_event", on_key_press)

#             text_result.config(text="PLOTTED runs : {} \n spin {} proj {}\nEnergy range :[{} , {} ]".format(
# len(filtered_run_list),self.hull_spin, self.hull_proj, self.Erange[0],
# self.Erange[1]) )

#         def clicked_button_save() :
#             try :
#                 from pathlib import Path
#                 home = str(Path.home())
#                 print("pathlib ok")
#             except Exception as ex :
#                 print(ex)
#                 home =  "/home/jvergnet/"

#             self.fig_file = filedialog.asksaveasfile(initialdir =home ,
#                                                      title = "Select File",
#                                                      filetypes=[('PNG', ".png"),("SVG",".svg"),("PDF",".pdf")])
# #            fmt = "pdf","svg","png","ps"] :
#             try :
#                 print(self.fig_file, self.fig_file.name)
#                 self.fig.savefig("{}".format(self.fig_file.name) , bbox_inches='tight')
#             except Exception as ex :
#                 print(ex)

#         def clicked_button_help() :
#             help_txt ="""
# HELP ON DOS !!
# """
#             tk.messagebox.showinfo('Proj help',help_txt)

#         # FRAME TOP  ===================================================
#         frame_top = tk.Frame(self)
#         frame_top.pack(fill=tk.BOTH, expand=True)

#         #=== SELECTION : spin value
#         frame_spin = tk.Frame(frame_top)
#         frame_spin.pack(fill=tk.BOTH, expand=True,side=tk.LEFT)

#         text_spin  = tk.Label(frame_spin, text = "Spin selection")
#         text_spin.pack(fill=tk.BOTH, expand=1,padx=5, pady=5,side=tk.TOP)

#         combo_spin = ttk.Combobox(frame_spin)
#         combo_spin['values']= ("up & down","sum","diff")
#         combo_spin.current(0)
#         combo_spin.pack(fill=tk.BOTH, expand=1,padx=5, pady=5,side=tk.TOP)

#         self.hull_spin = combo_spin.get()

#         #=== SELECTION : PROJECTION value
#         frame_proj = tk.Frame(frame_top)
#         frame_proj.pack(fill=tk.BOTH, expand=True,side=tk.LEFT)

#         text_proj  = tk.Label(frame_proj, text = "Projection selection ?")
#         text_proj.pack(fill=tk.BOTH, expand=1,padx=5, pady=5,side=tk.TOP)

#         combo_proj = ttk.Combobox(frame_proj)
#         combo_proj['values']= ( "spd","perSite","perElt","total" )
#         combo_proj.current(2)
#         combo_proj.pack(fill=tk.BOTH, expand=1,padx=5, pady=5,side=tk.TOP)

#         self.hull_proj = combo_proj.get()

#         #=== SELECTION : ENERGY RANGE OF DOS
#         frame_erange = tk.Frame(frame_top)
#         frame_erange.pack(fill=tk.BOTH, expand=True,side=tk.LEFT)

#         # == emin frame
#         frame_emin = tk.Frame(frame_erange)
#         frame_emin.pack(fill=tk.BOTH, expand=True,side=tk.LEFT)

#         emin_var = tk.IntVar()
#         emin_var.set(-3)
#         text_emin  = tk.Label(frame_emin, text = "Emin ?")
#         text_emin.pack(fill=tk.BOTH, expand=1,padx=5, pady=5,side=tk.TOP)
#         comboroll_emin = tk.Spinbox(frame_emin, from_=-10, to=0,textvariable=emin_var, width=5)
#         comboroll_emin.pack(fill=tk.BOTH, expand=1,padx=5, pady=5,side=tk.TOP)

#         # == emin frame
#         frame_emax = tk.Frame(frame_erange)
#         frame_emax.pack(fill=tk.BOTH, expand=True,side=tk.LEFT)

#         emax_var = tk.IntVar()
#         emax_var.set(3)
#         text_emax  = tk.Label(frame_emax, text = "Emax ?")
#         text_emax.pack(fill=tk.BOTH, expand=1,padx=5, pady=5,side=tk.TOP)
#         comboroll_emax = tk.Spinbox(frame_emax, from_=0, to=10,textvariable=emax_var, width=5)
#         comboroll_emax.pack(fill=tk.BOTH, expand=1,padx=5, pady=5,side=tk.TOP)

#         # FRAME MIDDLE  ================================================
#         frame_mid = tk.Frame(self)
#         frame_mid.pack(fill=tk.BOTH, expand=True)

#         # === SELECT
#         frame_select = tk.Frame(frame_mid)
#         frame_select.pack(fill=tk.BOTH, expand=True,side=tk.LEFT)

#         button_select = ttk.Button(frame_select, text = "SELECT parameters",
#                                         command=clicked_button_select)
#         button_select.pack(fill=tk.BOTH, expand=1,padx=5, pady=5,side=tk.TOP)


#         text_select  = tk.Label(frame_select, text = "DOS parameters")
#         text_select.pack(fill=tk.BOTH, expand=1,padx=5, pady=5,side=tk.TOP)
#         text_select.config(text="Selected parameters : \n spin {} proj {}".format(self.hull_spin, self.hull_proj) )

#         # === APPLY
#         frame_apply = tk.Frame(frame_mid)
#         frame_apply.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)

#         button_apply = ttk.Button(frame_apply, text = "PLOT",
#                                         command=clicked_button_apply)
#         button_apply.pack(fill=tk.BOTH, expand=1,padx=5, pady=5,side=tk.TOP)

#         text_result  = tk.Label(frame_apply, text = "Filtered runs : {}".format(len(filtered_run_list)) )
#         text_result.pack(fill=tk.BOTH, expand=1,padx=5, pady=5,side=tk.TOP)

#         # === HELP
#         frame_help_save = tk.Frame(frame_mid)
#         frame_help_save.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)
#         button_help = ttk.Button(frame_help_save, text = "HELP",
#                                         command=clicked_button_help)
#         button_help.pack(fill=tk.BOTH, expand=1,padx=5, pady=5,side=tk.TOP)

#         button_save = ttk.Button(frame_help_save, text = "SAVE",
#                                         command=clicked_button_save)
#         button_save.pack(fill=tk.BOTH, expand=1,padx=5, pady=5,side=tk.TOP)


# class Hull_page(Blank_page):

#     graph_type = "Hull"

#     def __init__(self, tab_control, controller, **kwargs):
#         Blank_page.__init__(self,tab_control, controller, **kwargs)


app = ReadRun_app()

# ###Buttons
# button1 = tk.Button(app, text = "Export plots", command = plot)
# button1.pack()

# Menu Bar
menubar = tk.Menu(app)

menu_file = tk.Menu(menubar, tearoff=0)
menu_file.add_command(label="Export on Plotly", command=Export)
menu_file.add_command(label="Export on Window", command=plot)
menu_file.add_separator()
menu_file.add_command(label="Quit", command=app.quit)
menubar.add_cascade(label="File", menu=menu_file)


menu_navig = tk.Menu(menubar, tearoff=0)

for p in app.frames.keys():
    action = partial(app.show_frame, p)
    menu_navig.add_command(label=p,
                           command=action)  # Start_page) )
menubar.add_cascade(label="Navigation", menu=menu_navig)

menu_help = tk.Menu(menubar, tearoff=0)
menu_help.add_command(label="Help !", command=Help)
menubar.add_cascade(label="Help", menu=menu_help)


app.config(menu=menubar)

app.mainloop()
