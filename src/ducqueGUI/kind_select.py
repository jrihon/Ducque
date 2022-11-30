import tkinter as tk
from tkinter import ttk

from ducqueGUI.grid_geometry import GridGeometry as GG

#  +--------------------------------------------------+
#  |                    BUILD                         |
#  +--------------------------------------------------+
class SelectApp(tk.Tk):

    def __init__(self, title):
        # baseline stuff
        super().__init__()
        self.title("Ducque : " + title)
        self.geometry(GG.window_size)

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

    def set_labels(self): 
        self.label_title = ttk.Label(self.content, textvariable=self.title)

    def set_buttons(self):
        self.yeppers = ttk.Button(self.content, text="yeppers", command=self.start_module)

#    def start_module(self):
#        self.destroy()
#        module = self.module_choices.get()
#        gui_window(module)

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
