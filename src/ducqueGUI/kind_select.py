import tkinter as tk
from tkinter import ttk

from ducqueGUI.grid_geometry import GridGeometry
import ducqueGUI.kind_build as KB
import ducqueGUI.kind_rand as KR
import ducqueGUI.kind_xyzpdb as KX
import ducqueGUI.kind_transmute as KT

#  +--------------------------------------------------+
#  |                   SELECT                         |
#  +--------------------------------------------------+
class SelectApp(tk.Tk):

    def __init__(self, title : str): #master
        # baseline stuff
        super().__init__() #master
        self.title("Ducque : " + title)
        self.geometry(GridGeometry.window_size)

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

    def set_buttons(self):
        self.yeppers = ttk.Button(self.content, text="yeppers", command=self.on_destroy)

    def on_destroy(self):
        self.module_kind = self.module_choices.get()
        print("Initiating " + self.module_kind + " module ...")
        self.destroy()
        self.module = self.gui_module(self.module_kind)

    def start_module(self, App, module):
        self.app = App(module)
        self.mainloop()         # ... RUN ! DUH DUUUH DUNDUNDUDUNUDDDUUNNN

    def gui_module(self, module):
        if module == "build" :
            self.start_module(KB.BuildApp, module)

        if module == "randomise" :
            self.start_module(KR.RandomiseApp, module)

        if module == "transmute" :
            self.start_module(KT.TransmuteApp, module)


        if module == "xyz_pdb" :
            self.start_module(KX.FormatPdbApp, module)




    def set_omenu(self):
        self.module_choices = tk.StringVar()
        module_opts = ["build", "transmute", "randomise", "xyz_pdb"]
        self.module_choices.set("build") # default value
        self.omenu_module = tk.OptionMenu(self.content, self.module_choices, *module_opts)

    def place_widgets(self): 

        self.content.grid(column=0,row=0)
        # labels
        self.label_title.grid(column=1, row=1, columnspan=2)

        # option menu
        self.omenu_module.grid(column=1, row=2, columnspan=2)

        # button
        self.yeppers.grid(column=1, row=3)
