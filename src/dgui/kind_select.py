import tkinter as tk
from tkinter import ttk

import systemsDucque as SD
from dgui.grid_geometry import Geometry as G
import dgui.kind_build as KB
import dgui.kind_rand as KR
import dgui.kind_transmute as KT
import dgui.kind_transmutelinker as KTL
import dgui.kind_nbase as KN
import dgui.kind_transmutenbase as KTN
#import dgui.kind_xyzpdb as KX # disabled module

#  +--------------------------------------------------+
#  |                   SELECT                         |
#  +--------------------------------------------------+
class SelectApp(tk.Tk):

    def __init__(self, title : str): #master
        # baseline stuff
        super().__init__() #master
        self.title("Ducque : " + title)
        self.geometry(G.window_size_SELECT)

        # Set Parent Frame
        self.content = ttk.Frame(self)

        # set labels
        self.set_labels()

        # set buttons
        self.set_buttons()

        # set option menu
        self.set_omenu()

        # place all the widgets
        self.place_widgets()

    def start(self):
        self.mainloop()

    def set_labels(self): 
        self.label_title = ttk.Label(self.content, text="Choose one of the following modules : ")

        self.label_empty = ttk.Label(self.content, text="")

    def set_buttons(self):
        self.Launch = ttk.Button(self.content, text="Launch", command=self.on_destroy)
        self.Launch.configure(width=G.SELECT_launch)

    def on_destroy(self):
        self.module_kind = self.module_choices.get()
        SD.print_launch(self.module_kind)
        self.destroy()
        self.module = self.gui_module(self.module_kind)

    def start_module(self, App, module):
        self.app = App(module)
        self.mainloop()         # ... RUN ! DUH DUUUH DUNDUNDUDUNUDDDUUNNN


    def gui_module(self, module):
        if module == "build" :
            self.start_module(KB.BuildApp, module)

        if module == "nbase" :
            self.start_module(KN.NbaseApp, module)

        if module == "randomise" :
            self.start_module(KR.RandomiseApp, module)

        if module == "transmute" :
            self.start_module(KT.TransmuteApp, module)

        if module == "tlinker" :
            self.start_module(KTL.TransmuteLinkerApp, module)

        if module == "tbase" :
            self.start_module(KTN.TransmuteTbaseApp, module)

        # This module has been disabled for the GUI
        # if module == "xyz_pdb" :
        #     self.start_module(KX.FormatPdbApp, module)


    def set_omenu(self):
        self.module_choices = tk.StringVar()
        module_opts = ["build",
                       "nbase",
                       "transmute",
                       "tlinker",
                       "tbase",
                       "randomise",
                     # "xyz_pdb" # This module has been disabled for the GUI
                        ]

        self.module_choices.set("build") # default value
        self.omenu_module = tk.OptionMenu(self.content, self.module_choices, *module_opts)
        self.omenu_module.configure(width=G.SELECT_optmenu)

    def place_widgets(self): 

        self.content.grid(column=0,row=0)
        # labels
        self.label_title.grid(column=1, row=1, columnspan=2)

        # option menu
        self.omenu_module.grid(column=1, row=2, columnspan=2)

        # button
        self.Launch.grid(column=1, row=3, columnspan=2)
