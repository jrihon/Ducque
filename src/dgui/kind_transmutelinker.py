import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from dgui.grid_geometry import Geometry as G

from os import getcwd 
from os.path import isfile
from shutil import which   # Check if Ducque is on the $PATH
from subprocess import run # Run Ducque

import systemsDucque as SD
from ducquelib.library import TABLE_BACKBONE, TABLE_LINKER

#  +--------------------------------------------------+
#  |                    TRANSMUTE                     |
#  +--------------------------------------------------+

class TransmuteLinkerApp(tk.Tk):

    def __init__(self, title):
        # baseline stuff
        super().__init__()
        self.title("Ducque : " + title)
        self.geometry(G.window_size_TRANSMUTE)
        self.cwd = getcwd()
        self.padding = {"padx" : 3, "pady" : G.pady}

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

        # set optionmenu
        self.set_optionmenu()

        # place all the widgets
        self.place_widgets()

    def set_labels(self):
        self.label_pdbfname = ttk.Label(self.content, text="--pdb")
        self.label_chemistry = ttk.Label(self.content, text="--chemistry")
        self.label_moiety = ttk.Label(self.content, text="--moiety")
        self.label_conformation = ttk.Label(self.content, text="--conformation")
        self.label_bondangles = ttk.Label(self.content, text="--bondangles")
        self.label_dihedrals = ttk.Label(self.content, text="--dihedrals")

        self.label_ang1 = ttk.Label(self.content, text="angle 1")
        self.label_dihr1 = ttk.Label(self.content, text="dihedral 1")
        self.label_dihr2 = ttk.Label(self.content, text="dihedral 2")


    def set_buttons(self):
        self.btn_write = ttk.Button(self.content, text="Write input file", command=self.write_inputfile, width=20)
        self.btn_transmute = ttk.Button(self.content, text="Transmute!", command=self.transmute_input, width=20)


#    def toggle_zeta(self):
#        if self.int_zeta.get() == 1 : 
#            self.ent_ang_z.config(state="enabled")
#            self.ent_dihr_z.config(state="enabled")
#        elif self.int_zeta.get() == 0  and self.int_nu.get() == 1: 
#            SD.angle_exclusivity()
#
#            self.ent_ang_z.config(state="enabled")
#            self.ent_dihr_z.config(state="enabled")
#            self.int_zeta.set(1)
#        else :
#            self.ent_ang_z.config(state="disabled")
#            self.ent_dihr_z.config(state="disabled")
#
#    def toggle_nu(self):
#        if self.int_nu.get() == 1  and self.int_zeta.get() == 1: 
#            self.ent_ang_n.config(state="enabled")
#            self.ent_dihr_n.config(state="enabled")
#        elif self.int_nu.get() == 1  and self.int_zeta.get() == 0: 
#            SD.angle_exclusivity()
#
#            self.ent_ang_n.config(state="disabled")
#            self.ent_dihr_n.config(state="disabled")
#            self.int_nu.set(0)
#        else :
#            self.ent_ang_n.config(state="disabled")
#            self.ent_dihr_n.config(state="disabled")
#
#
    def set_checkbutton(self):
#        self.int_zeta = tk.IntVar()
#        self.int_nu = tk.IntVar()
#        self.chkbtn_zeta = ttk.Checkbutton(self.content, text="Zeta", variable=self.int_zeta, command=self.toggle_zeta,
#                                                onvalue=1, offvalue=0)
#        self.chkbtn_nu = ttk.Checkbutton(self.content, text="Eta", variable=self.int_nu, command=self.toggle_nu,
#                                                onvalue=1, offvalue=0)
#
#        self.int_zeta.set(1)
#        self.ent_ang_n.config(state="disabled")
#        self.ent_dihr_n.config(state="disabled")

        self.int_overwrite = tk.IntVar()
        self.int_overwrite.set(0)
        self.chkbtn_overwrite = ttk.Checkbutton(self.content, variable=self.int_overwrite, text="Overwrite", onvalue=1, offvalue=0)

        self.int_build = tk.IntVar()
        self.chkbtn_build = ttk.Checkbutton(self.content, variable=self.int_build, command=self.build_toggle, text="Build module", onvalue=1, offvalue=0)

    def set_entries(self):
        self.str_pdbfname = tk.StringVar()
#        self.str_nucleobase = tk.StringVar()
        
        self.entr_pdbfname = ttk.Entry(self.content, textvariable=self.str_pdbfname, width=42)
#        self.entr_nucleobase = ttk.Entry(self.content, textvariable=self.str_nucleobase)

        # moiety entry
        self.str_moiety = tk.StringVar()
        self.entr_moiety = ttk.Entry(self.content, textvariable=self.str_moiety)
        self.str_moiety.set("linker")
        self.entr_moiety.config(state="disabled")

        # bondangles
        self.str_angle1 = tk.StringVar() 
        self.ent_angle1 = ttk.Entry(self.content, textvariable=self.str_angle1) # Angle_1
    
        # Dihedrals
        self.str_dihr1 = tk.StringVar() # Dihedral_1
        self.str_dihr2 = tk.StringVar() # Dihedral_2

        self.ent_dihr1 = ttk.Entry(self.content, textvariable=self.str_dihr1) # Dihedral_1
        self.ent_dihr2 = ttk.Entry(self.content, textvariable=self.str_dihr2) # Dihedral_2


    def set_optionmenu(self):
        # chemistries
        self.choice_chem = tk.StringVar()
        opt_chems = list(TABLE_BACKBONE.keys())
        self.opt_chems = ["..."] + sorted([x for x in opt_chems if x != "PHOSPHATE"], key=str.casefold) # replace the phosphate key with the `...` key and sort
        self.omenu_chem = ttk.OptionMenu(self.content, self.choice_chem, *self.opt_chems)
        self.omenu_chem.configure(width=15)

        # conformations
        self.choice_conf = tk.StringVar()
        self.opt_confs = ["...", "NONE", "R", "S"]
        self.omenu_confs = ttk.OptionMenu(self.content, self.choice_conf, *self.opt_confs)
        self.omenu_confs.configure(width=15)


    def file_dialog(self):

        filetypes = (("All files", "*.*"), ("input-files", "*.tinp"))

        select_files = filedialog.askopenfilenames(
                                        title="Import files : " + self.cwd,
                                        initialdir= self.cwd,
                                        filetypes=filetypes
                                        )

        # This is why Python will never live up to the standards of Rust
        # This is an unreadable block of code but it works and I dislike it

        try :   # see if self.file_queried exists
            self.file_queried 
        except IndexError :
            SD.print_empty_query("IMPORT INPUT FILE")
            return
        except :   # if it does not exist yet, create it with whatever has been queried in the file dialog
            try :
                self.file_queried = select_files[0]  # try to parse from the file dialog
            except IndexError:                       # if the list does not exist, return early
                SD.print_empty_query("IMPORT INPUT FILE")
                return
        else :  # if self.file_queried already exists, parse from the given list or do nothing if nothing has been queried
            try :
                self.file_queried = select_files[0]  # try to parse from the list
            except : pass




        # Setting strings to the variables.
        # If this errors out, it will be caught when Transmute runs
        with open(self.file_queried, "r") as inputfile :
            file_content = inputfile.readlines()

        for line in file_content :
            flag , inp = line.split(" ", maxsplit=1)

            if flag == "--pdb" : 
                self.entr_pdbfname.configure()
                self.str_pdbfname.set(inp.strip())

            if flag == "--chemistry" :
                opt = inp.strip()
                if opt in list(TABLE_BACKBONE.keys()) : self.choice_chem.set(inp.strip())
                else : SD.print_invalid_argument(opt, "`--chemistry`")

            if flag == "--conformation" :
                opt = inp.strip()
                if opt in ["NONE", "R", "S"] : self.choice_conf.set(inp.strip())
                else : SD.print_invalid_argument(opt, "`--conformation`")

#            if flag == "--moiety" :
#                if inp.strip() == "nucleoside": 
#                    print("--moiety `nucleoside` correctly prompted.")

            if flag == "--bondangles" : 
                angs = list(map(lambda x: x.strip(), inp.split(",")))
                if len(angs) != 1:
                    SD.print_insuf_amount("--bondangles")
                    return
                self.str_angle1.set(angs[0])

            if flag == "--dihedrals" :
                dihrs = list(map(lambda x: x.strip(), inp.split(",")))
                if len(dihrs) != 2 : 
                    SD.print_insuf_amount("--dihedrals")
                    return
                self.str_dihr1.set(dihrs[0])
                self.str_dihr2.set(dihrs[1])


    def build_toggle(self):

        if self.int_build.get() == 1:
            self.set_build_button()

        else : 
            self.build_files_btn.destroy()
            self.build_run_btn.destroy()
            try :
                self.buildlabel
            except : pass
            else : self.buildlabel.destroy()

    def set_build_button(self):
        # set filedialog for build option
        self.build_files_btn = ttk.Button(self.content, text="Import build file", command=self.builder_dialog, width=19)
        self.build_files_btn.grid(column=2, row=12, **self.padding)

        self.build_run_btn = ttk.Button(self.content, text="Build!", command=self.build_command, width=19)
        self.build_run_btn.grid(column=3, row=12, **self.padding)

    def build_command(self):

        try : 
            self.buildinput
        except : 
            SD.print_empty_query("IMPORT BUILD FILE")
            return

        # At this point, this would not be necessary, but better safe than sorry
        if not which("Ducque"): 
            SD.print_cant_find_Ducque()

        run(["Ducque", "--build", self.buildinput])

    def builder_dialog(self):

        if self.int_build.get() == 1:
            build_filetypes = (("All files", "*.*"), ("input-files", "*.binp"))

            select_files = filedialog.askopenfilenames(
                                            title="Import files : " + self.cwd, # OR the $DUCQUE/transmute directory
                                            initialdir= self.cwd,               
                                            filetypes=build_filetypes
                                            )

            try :
                file_queried = select_files[0]
            except IndexError:
                SD.print_empty_query("IMPORT BUILD (`.binp`) FILE")
                return

            try :
                self.buildlabel
            except :
                self.buildinput = file_queried
                self.buildlabel = ttk.Label(self.content, text= "Input file :  " + self.buildinput)
                self.buildlabel.grid(column=4, row=12, columnspan=7)
            else :
                self.buildlabel.destroy()
                self.buildinput = file_queried
                self.buildlabel = ttk.Label(self.content, text= "Input file :  " + self.buildinput)
                self.buildlabel.grid(column=4, row=12, columnspan=7)


        if self.int_build.get() == 0:
            pass

    def place_widgets(self):

        self.content.grid(column=0, row=0)

        self.label_words = ttk.Label(self.content, text="Convert your `pdb` to an input for Ducque")
        self.label_words.grid(column=1, row=1, columnspan=2)

        # set filedialog
        self.open_files_btn = ttk.Button(self.content, text="Import input file", command=self.file_dialog)
        self.open_files_btn.grid(column=2, row=2)

        # set labels
        self.label_pdbfname.grid(column=0, row=3, **self.padding)
        self.label_chemistry.grid(column=0, row=4, **self.padding)
        self.label_conformation.grid(column=0, row=5, **self.padding)
        self.label_moiety.grid(column=0, row=6, **self.padding)

        # bond angles
        self.label_ang1.grid(column=1, row=7)
        self.label_bondangles.grid(column=0, row=8, **self.padding)
        self.ent_angle1.grid(column=1, row=8, **self.padding)

        # dihedrals
        self.label_dihr1.grid(column=1, row=9)
        self.label_dihr2.grid(column=2, row=9)

        self.label_dihedrals.grid(column=0, row=10, **self.padding)
        self.ent_dihr1.grid(column=1, row=10, **self.padding)
        self.ent_dihr2.grid(column=2, row=10, **self.padding)

        # set entries
        self.entr_pdbfname.grid(column=1, row=3, columnspan=2, **self.padding)
        self.entr_moiety.grid(column=1, row=6, **self.padding)

        # set optionmenu
        self.omenu_chem.grid(column=1, row=4, **self.padding)
        self.omenu_confs.grid(column=1, row=5, **self.padding)

        # set buttons
        self.btn_write.grid(column=7, row=11, **self.padding)
        self.btn_transmute.grid(column=7, row=12, **self.padding)
        self.chkbtn_overwrite.grid(column=8, row=12, **self.padding)

        self.chkbtn_build.grid(column=1, row=13, **self.padding)


    def write_inputfile(self):

        import re   # regex library

        def format_angle(name_angle : str, angle : str):

            if len(angle.strip()) == 0 :
                SD.print_empty_query(name_angle)
                return ""

            a = re.sub(",", ".", angle) # replace any commas with points for each angle entry
            try :
                float(a)
            except ValueError: 
                SD.print_conversion_err(name_angle, angle)
                return ""
            else :
                return a
        # Handle empty inputs for these fields
        for string in [self.str_pdbfname, self.choice_chem, self.choice_conf] :
            if len(string.get()) == 0 or string.get() == "..." : 
                SD.print_empty_query("--pdb, --chemistry or --conformation")
                return

#        # Sort of handle numerical inputs
        list_ang = [ format_angle("angle(angle_1)", self.str_angle1.get())]
        list_dih = [ format_angle("dihedral(dihedral_1)", self.str_dihr1.get()),
                     format_angle("dihedral(dihedral_2)", self.str_dihr2.get())]

        # Check if any angles are empty, if so, return early
        if "" in list_ang or "" in list_dih:
            SD.print_empty_query("Angles or Dihedrals")
            return

        # This means that when we import a file to be read in by the GUI, that we will have it start out in the current directory
        # OR the $DUCQUE/transmute directory
        # Add the related residue name to this, from 
        # Linking the nucleoside to the correct linker moiety 
        try: 
            linker = TABLE_LINKER[self.choice_chem.get().upper()]
        except: 
            SD.print_invalid_key(self.choice_chem.get().upper(), TABLE_LINKER)
            return

        input_fname = self.choice_chem.get().lower() + "_" + linker.lower() + ".tinp"


        if self.choice_conf.get() != "NONE": # meaning it is either R or S
            input_fname = self.choice_conf.get().lower() + input_fname

        if not self.int_overwrite.get() == 1 and isfile(input_fname) :
            SD.print_no_overwrite(input_fname, getcwd())
            return

            
        # This file has to be generated in the $DUCQUEHOME/transmute/* folder and not in the current one
        with open(input_fname , "w") as fileto :
            fileto.write("--pdb " + self.str_pdbfname.get() + "\n"
                        "--chemistry " + self.choice_chem.get() + " \n"  
                        "--conformation " + self.choice_conf.get() + " \n"
                        "--moiety " + self.str_moiety.get() + " \n"
                        "--bondangles " + ", ".join(list_ang) + " \n"
                        "--dihedrals " + ", ".join(list_dih) + " \n"
                    )

        SD.print_writing(input_fname)


    def transmute_input(self):

        # Check if chemistry is queried
        chemchoice = self.choice_chem.get().lower()
        if len(chemchoice) == 0 or chemchoice == "...":
            SD.print_empty_query("--chemistry")
            return

        # filename for the transmute file
        try: 
            linker = TABLE_LINKER[self.choice_chem.get().upper()]
        except: 
            SD.print_invalid_key(self.choice_chem.get().upper(), TABLE_LINKER)
            return

        json_fname = chemchoice + "_" + linker.lower() + ".json"
        input_fname = chemchoice + "_" + linker.lower() + ".tinp"

        if self.choice_conf.get() != "NONE": # meaning it is either R or S
            input_fname = self.choice_conf.get().lower() + input_fname
            json_fname = self.choice_conf.get().lower() + json_fname

        if not which("Ducque"):  # At this point, this would not be necessary, but better safe than sorry
            SD.print_cant_find_Ducque()

        # Check if the transmute input file is present on the current working directory
        if not isfile(input_fname):
            SD.print_filenotfound(input_fname)
            return

        # Check if the json is not already present in the DUCQUEHOME/json/
        json_dir = SD.return_DUCQUEHOME() + "json/"
        if not isfile(json_dir + json_fname): 
            run(["Ducque", "--transmute", input_fname])

        else :
            # check if the overwrite button is in which state 
            if self.int_overwrite.get() == 0 :
                SD.print_no_overwrite(json_fname, json_dir)
                return
            if self.int_overwrite.get() == 1 :
                run(["Ducque", "--transmute", input_fname])
