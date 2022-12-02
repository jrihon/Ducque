import tkinter as tk
from tkinter import filedialog
from tkinter import ttk

from os import getcwd
#from shutil import which   # Run Ducque
#from subprocess import run # Run Ducque

from ducqueGUI.grid_geometry import GridGeometry as GG

#  +--------------------------------------------------+
#  |                    XYZ 2 PDB                     |
#  +--------------------------------------------------+

class FormatPdbApp(tk.Tk):

    def __init__(self, title):
        # baseline stuff
        super().__init__()
        self.title("Ducque : " + title)
        self.geometry(GG.window_size)
        self.cwd = getcwd()

        # Set Parent Frame
        self.content = ttk.Frame(self)

        # set labels
        self.set_labels()


        # set buttons
        self.set_buttons()

        # set entries
#        self.set_entries()

        # set checkbutton
#        self.set_checkbutton()

        # place all the widgets
        self.place_widgets()

    def set_labels(self):

        self.xyztoptext = ttk.Label(self.content, text="Convert from a `xyz` to a `pdb` format.")

    def set_flag_labels(self): 

        self.str_fname = "--xyz"
        self.str_resname = "--residue"
        self.str_atomlist = "--atomname_list"

        self.label_fname = ttk.Label(self.content, text=self.str_fname)
        self.label_resname = ttk.Label(self.content, text=self.str_resname)
        self.label_atomlist = ttk.Label(self.content, text=self.str_atomlist)

        self.label_fname.grid(column=0, row=3)
        self.label_resname.grid(column=0, row=4)
        self.label_atomlist.grid(column=0, row=5)

    def set_flag_entries(self): 

        self.ent_fname = tk.StringVar()
        self.ent_resname = tk.StringVar()

        self.entry_fname = ttk.Entry(self.content, textvariable=self.ent_fname)
        self.entry_resname = ttk.Entry(self.content, textvariable=self.ent_resname)

        self.entry_fname.grid(column=1, row=3)
        self.entry_resname.grid(column=1, row=4)

    def file_dialog(self):

        filetypes = (("All files", "*.*"), ("xyz-files", "*.xyz"))

        select_files = filedialog.askopenfilenames(
                                        title="Open files : " + self.cwd,
                                        initialdir= self.cwd,
                                        filetypes=filetypes
                                        )
        try :
            file_queried = select_files[0]
        except IndexError:
            print("No file selected. Please try again")
            return

        with open(file_queried, "r") as xyz :
            file_content = xyz.readlines()[2:]

        self.elements = [x.split()[0] for x in file_content]
        self.len_xyzfile = len(self.elements)

        # Generate the labels and entries
        self.geometry("769x300") # reset geometry of window

        width = 120
        self.elem_label = ", ".join([str(x + 1) + "." + self.elements[x] for x in range(self.len_xyzfile)])
        self.atom_labels = ttk.Label(self.content, text=self.elem_label)
        self.atom_labels.config(font=('TkDefaultFont', 9))
        self.atom_labels.grid(column=1, row=6, columnspan=3)


        self.atom_str = tk.StringVar()
        self.atom_str.set(self.elem_label)
        self.atom_entry = ttk.Entry(self.content, textvariable=self.atom_str, width=width)
        self.atom_entry.grid(column=1, row=5, columnspan=3)

        self.set_flag_labels()
        self.set_flag_entries()


    def print_entries(self):
        pass

    def set_buttons(self):
        self.btn_printentries = ttk.Button(self.content, text="Print!", command=self.print_entries)
            

    def place_widgets(self):

        self.content.grid(column=0,row=0)

        self.open_files_btn = ttk.Button(self.content, text="OPEN FILES", command=self.file_dialog)
        self.open_files_btn.grid(column=1, row=2)

        self.xyztoptext.grid(column=1, columnspan=2, row=1)

