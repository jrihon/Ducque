import tkinter as tk
from tkinter import ttk

from shutil import which   # Run Ducque
from subprocess import run # Run Ducque

from builder.builder_library import backbone_codex # import possibilities to build complementary strand
from ducqueGUI.grid_geometry import GridGeometry as GG

#  +--------------------------------------------------+
#  |                    RANDOMISE                     |
#  +--------------------------------------------------+

class TransmuteApp(tk.Tk):

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

        # set entries
        self.set_entries()

        # set entries
        self.set_checkbutton()

        # place all the widgets
        self.place_widgets()

