import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from dgui.grid_geometry import Geometry as G

from os import getcwd 
from os.path import isfile
from shutil import which   # Check if Ducque is on the $PATH
from subprocess import run # Run Ducque

import systemsDucque as SD
#from builder.builder_library import TABLE_BACKBONE 
#from transmute.transmute_library import TABLE_NUCLEOBASE, TABLE_BACKBONE
from ducquelib.library import TABLE_NUCLEOBASE, TABLE_BACKBONE

#  +--------------------------------------------------+
#  |                    TRANSMUTE                     |
#  +--------------------------------------------------+

class TransmuteApp(tk.Tk):

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
        self.label_conformation = ttk.Label(self.content, text="--conformation")
        self.label_moiety = ttk.Label(self.content, text="--moiety")
        self.label_bondangles = ttk.Label(self.content, text="--bondangles")
        self.label_dihedrals = ttk.Label(self.content, text="--dihedrals")
        self.label_nucleobase = ttk.Label(self.content, text="--nucleobase")

        self.label_alpha = ttk.Label(self.content, text="Alpha")
        self.label_beta = ttk.Label(self.content, text="Beta")
        self.label_gamma = ttk.Label(self.content, text="Gamma")
        self.label_delta = ttk.Label(self.content, text="Delta")
        self.label_epsilon = ttk.Label(self.content, text="Epsilon")
        self.label_zeta = ttk.Label(self.content, text="Zeta")
        self.label_nu = ttk.Label(self.content, text="Eta")
        self.label_chi = ttk.Label(self.content, text="Chi")


    def set_buttons(self):
        self.btn_write = ttk.Button(self.content, text="Write input file", command=self.write_inputfile, width=20)
        self.btn_transmute = ttk.Button(self.content, text="Transmute!", command=self.transmute_input, width=20)


    def toggle_zeta(self):
        if self.int_zeta.get() == 1 : 
            self.ent_ang_z.config(state="enabled")
            self.ent_dihr_z.config(state="enabled")
        elif self.int_zeta.get() == 0  and self.int_nu.get() == 1: 
            SD.angle_exclusivity()

            self.ent_ang_z.config(state="enabled")
            self.ent_dihr_z.config(state="enabled")
            self.int_zeta.set(1)
        else :
            self.ent_ang_z.config(state="disabled")
            self.ent_dihr_z.config(state="disabled")

    def toggle_nu(self):
        if self.int_nu.get() == 1  and self.int_zeta.get() == 1: 
            self.ent_ang_n.config(state="enabled")
            self.ent_dihr_n.config(state="enabled")
        elif self.int_nu.get() == 1  and self.int_zeta.get() == 0: 
            SD.angle_exclusivity()

            self.ent_ang_n.config(state="disabled")
            self.ent_dihr_n.config(state="disabled")
            self.int_nu.set(0)
        else :
            self.ent_ang_n.config(state="disabled")
            self.ent_dihr_n.config(state="disabled")


    def set_checkbutton(self):
        self.int_zeta = tk.IntVar()
        self.int_nu = tk.IntVar()
        self.chkbtn_zeta = ttk.Checkbutton(self.content, text="Zeta", variable=self.int_zeta, command=self.toggle_zeta,
                                                onvalue=1, offvalue=0)
        self.chkbtn_nu = ttk.Checkbutton(self.content, text="Eta", variable=self.int_nu, command=self.toggle_nu,
                                                onvalue=1, offvalue=0)

        self.int_zeta.set(1)
        self.ent_ang_n.config(state="disabled")
        self.ent_dihr_n.config(state="disabled")

        self.int_overwrite = tk.IntVar()
        self.int_overwrite.set(0)
        self.chkbtn_overwrite = ttk.Checkbutton(self.content, variable=self.int_overwrite, text="Overwrite", onvalue=1, offvalue=0)

        self.int_build = tk.IntVar()
        self.chkbtn_build = ttk.Checkbutton(self.content, variable=self.int_build, command=self.build_toggle, text="Build module", onvalue=1, offvalue=0)

    def set_entries(self):
        self.str_pdbfname = tk.StringVar()
        self.str_conformation = tk.StringVar()
        self.str_nucleobase = tk.StringVar()
        
        self.entr_pdbfname = ttk.Entry(self.content, textvariable=self.str_pdbfname, width=42)
        self.entr_conformation = ttk.Entry(self.content, textvariable=self.str_conformation)
        self.entr_nucleobase = ttk.Entry(self.content, textvariable=self.str_nucleobase)

        # moiety entry
        self.str_moiety = tk.StringVar()
        self.entr_moiety = ttk.Entry(self.content, textvariable=self.str_moiety)
        self.str_moiety.set("nucleoside")
        self.entr_moiety.config(state="disabled")

        # bondangles
        self.str_ang_a = tk.StringVar() # alpha
        self.str_ang_b = tk.StringVar() # beta
        self.str_ang_g = tk.StringVar() # gamma
        self.str_ang_d = tk.StringVar() # delta
        self.str_ang_e = tk.StringVar() # epsilon
        self.str_ang_z = tk.StringVar() # zeta
        self.str_ang_n = tk.StringVar() # nu
        self.str_ang_x = tk.StringVar() # chi

        self.ent_ang_a = ttk.Entry(self.content, textvariable=self.str_ang_a) # alpha
        self.ent_ang_b = ttk.Entry(self.content, textvariable=self.str_ang_b) # beta
        self.ent_ang_g = ttk.Entry(self.content, textvariable=self.str_ang_g) # gamma
        self.ent_ang_d = ttk.Entry(self.content, textvariable=self.str_ang_d) # delta
        self.ent_ang_e = ttk.Entry(self.content, textvariable=self.str_ang_e) # epsilon
        self.ent_ang_z = ttk.Entry(self.content, textvariable=self.str_ang_z) # zeta
        self.ent_ang_n = ttk.Entry(self.content, textvariable=self.str_ang_n) # nu
        self.ent_ang_x = ttk.Entry(self.content, textvariable=self.str_ang_x) # chi
    
        # Dihedrals
        self.str_dihr_a = tk.StringVar() # alpha
        self.str_dihr_b = tk.StringVar() # beta
        self.str_dihr_g = tk.StringVar() # gamma
        self.str_dihr_d = tk.StringVar() # delta
        self.str_dihr_e = tk.StringVar() # epsilon
        self.str_dihr_z = tk.StringVar() # zeta
        self.str_dihr_n = tk.StringVar() # nu
        self.str_dihr_x = tk.StringVar() # chi

        self.ent_dihr_a = ttk.Entry(self.content, textvariable=self.str_dihr_a) # alpha
        self.ent_dihr_b = ttk.Entry(self.content, textvariable=self.str_dihr_b) # beta
        self.ent_dihr_g = ttk.Entry(self.content, textvariable=self.str_dihr_g) # gamma
        self.ent_dihr_d = ttk.Entry(self.content, textvariable=self.str_dihr_d) # delta
        self.ent_dihr_e = ttk.Entry(self.content, textvariable=self.str_dihr_e) # epsilon
        self.ent_dihr_z = ttk.Entry(self.content, textvariable=self.str_dihr_z) # zeta
        self.ent_dihr_n = ttk.Entry(self.content, textvariable=self.str_dihr_n) # nu
        self.ent_dihr_x = ttk.Entry(self.content, textvariable=self.str_dihr_x) # chi


    def set_optionmenu(self):
        # chemistries
        self.choice_chem = tk.StringVar()
        opt_chems = list(TABLE_BACKBONE.keys())
        self.opt_chems = ["..."] + sorted([x for x in opt_chems if x != "PHOSPHATE"], key=str.casefold) # replace the phosphate key with the `...` key and sort
        self.omenu_chem = ttk.OptionMenu(self.content, self.choice_chem, *self.opt_chems)
        self.omenu_chem.configure(width=15)

        # moiety
#        self.choice_moi = tk.StringVar()
#        self.opt_moi = ["...", "nucleoside", "linker"]
#        self.omenu_moi = ttk.OptionMenu(self.content, self.choice_moi, *self.opt_moi)
#        self.choice_moi.set("nucleoside")
#        self.omenu_moi.configure(width=15)


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
                self.str_conformation.set(inp.strip())
            
            if flag == "--moiety" :
                if inp.strip() == "nucleoside": 
                    print("--moiety `nucleoside` correctly prompted.")

            if flag == "--nucleobase" :
                self.str_nucleobase.set(inp.strip()) 

            if flag == "--bondangles" : 
                angs = list(map(lambda x: x.strip(), inp.split(",")))
                if len(angs) <= 5 or len(angs) >= 9 : 
                    SD.print_insuf_amount("--bondangles")
                    return
                if len(angs) == 6 :
                    self.int_zeta.set(0)
                    self.ent_ang_z.config(state="disabled")
                    self.ent_dihr_z.config(state="disabled")
                    self.int_nu.set(0)
                    self.ent_ang_n.config(state="disabled")
                    self.ent_dihr_n.config(state="disabled")

                    self.str_ang_a.set(angs[0])
                    self.str_ang_b.set(angs[1])
                    self.str_ang_g.set(angs[2])
                    self.str_ang_d.set(angs[3])
                    self.str_ang_e.set(angs[4])
                    self.str_ang_x.set(angs[5])
                if len(angs) == 7 :

                    self.int_zeta.set(1)
                    self.ent_ang_z.config(state="enabled")
                    self.ent_dihr_z.config(state="enabled")
                    self.int_nu.set(0)
                    self.ent_ang_n.config(state="disabled")
                    self.ent_dihr_n.config(state="disabled")

                    self.str_ang_a.set(angs[0])
                    self.str_ang_b.set(angs[1])
                    self.str_ang_g.set(angs[2])
                    self.str_ang_d.set(angs[3])
                    self.str_ang_e.set(angs[4])
                    self.str_ang_z.set(angs[5])
                    self.str_ang_x.set(angs[6])
                if len(angs) == 8 :
                    self.int_zeta.set(1)
                    self.ent_ang_z.config(state="enabled")
                    self.ent_dihr_z.config(state="enabled")
                    self.int_nu.set(1)
                    self.ent_ang_n.config(state="enabled")
                    self.ent_dihr_n.config(state="enabled")

                    self.str_ang_a.set(angs[0])
                    self.str_ang_b.set(angs[1])
                    self.str_ang_g.set(angs[2])
                    self.str_ang_d.set(angs[3])
                    self.str_ang_e.set(angs[4])
                    self.str_ang_z.set(angs[5])
                    self.str_ang_n.set(angs[6])
                    self.str_ang_x.set(angs[7])

            if flag == "--dihedrals" :
                dihrs = list(map(lambda x: x.strip(), inp.split(",")))
                if len(dihrs) <= 5 or len(dihrs) >= 9 : 
                    SD.print_insuf_amount("--dihedrals")
                    return
                #
                if len(dihrs) == 6 :
                    self.int_zeta.set(0)
                    self.ent_ang_z.config(state="disabled")
                    self.ent_dihr_z.config(state="disabled")
                    self.int_nu.set(0)
                    self.ent_ang_n.config(state="disabled")
                    self.ent_dihr_n.config(state="disabled")

                    self.str_dihr_a.set(dihrs[0])
                    self.str_dihr_b.set(dihrs[1])
                    self.str_dihr_g.set(dihrs[2])
                    self.str_dihr_d.set(dihrs[3])
                    self.str_dihr_e.set(dihrs[4])
                    self.str_dihr_x.set(dihrs[5])
                #
                if len(dihrs) == 7 :
                    self.int_zeta.set(1)
                    self.ent_ang_z.config(state="enabled")
                    self.ent_dihr_z.config(state="enabled")
                    self.int_nu.set(0)
                    self.ent_ang_n.config(state="disabled")
                    self.ent_dihr_n.config(state="disabled")

                    self.str_dihr_a.set(dihrs[0])
                    self.str_dihr_b.set(dihrs[1])
                    self.str_dihr_g.set(dihrs[2])
                    self.str_dihr_d.set(dihrs[3])
                    self.str_dihr_e.set(dihrs[4])
                    self.str_dihr_z.set(dihrs[5])
                    self.str_dihr_x.set(dihrs[6])
                #
                if len(dihrs) == 8 :
                    self.int_zeta.set(1)
                    self.ent_ang_z.config(state="enabled")
                    self.ent_dihr_z.config(state="enabled")
                    self.int_nu.set(1)
                    self.ent_ang_n.config(state="enabled")
                    self.ent_dihr_n.config(state="enabled")

                    self.str_dihr_a.set(dihrs[0])
                    self.str_dihr_b.set(dihrs[1])
                    self.str_dihr_g.set(dihrs[2])
                    self.str_dihr_d.set(dihrs[3])
                    self.str_dihr_e.set(dihrs[4])
                    self.str_dihr_z.set(dihrs[5])
                    self.str_dihr_n.set(dihrs[6])
                    self.str_dihr_x.set(dihrs[7])


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
                                            initialdir= self.cwd,               # OR the $DUCQUE/transmute directory
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
        self.label_nucleobase.grid(column=2, row=5, **self.padding)

        # alphabet
        self.label_alpha.grid(column=1, row=7)
        self.label_beta.grid(column=2, row=7)
        self.label_gamma.grid(column=3, row=7)
        self.label_delta.grid(column=4, row=7)
        self.label_epsilon.grid(column=5, row=7)
        self.chkbtn_zeta.grid(column=6, row=7)
        self.chkbtn_nu.grid(column=7, row=7)
        self.label_chi.grid(column=8, row=7)

        self.label_bondangles.grid(column=0, row=8, **self.padding)
        self.label_dihedrals.grid(column=0, row=9, **self.padding)

        # bondangles
        self.ent_ang_a.grid(column=1, row=8, **self.padding)
        self.ent_ang_b.grid(column=2, row=8, **self.padding)
        self.ent_ang_g.grid(column=3, row=8, **self.padding)
        self.ent_ang_d.grid(column=4, row=8, **self.padding)
        self.ent_ang_e.grid(column=5, row=8, **self.padding)
        self.ent_ang_z.grid(column=6, row=8, **self.padding)
        self.ent_ang_n.grid(column=7, row=8, **self.padding)
        self.ent_ang_x.grid(column=8, row=8, **self.padding)

        # dihedrals
        self.ent_dihr_a.grid(column=1, row=9, **self.padding)
        self.ent_dihr_b.grid(column=2, row=9, **self.padding)
        self.ent_dihr_g.grid(column=3, row=9, **self.padding)
        self.ent_dihr_d.grid(column=4, row=9, **self.padding)
        self.ent_dihr_e.grid(column=5, row=9, **self.padding)
        self.ent_dihr_z.grid(column=6, row=9, **self.padding)
        self.ent_dihr_n.grid(column=7, row=9, **self.padding)
        self.ent_dihr_x.grid(column=8, row=9, **self.padding)

        # set entries
        self.entr_pdbfname.grid(column=1, row=3, columnspan=2, **self.padding)
        self.entr_conformation.grid(column=1, row=5, **self.padding)
        self.entr_nucleobase.grid(column=3, row=5, **self.padding)

        # set optionmenu
        self.omenu_chem.grid(column=1, row=4, **self.padding)
        self.entr_moiety.grid(column=1, row=6, **self.padding)

        # set buttons
        self.btn_write.grid(column=7, row=10, **self.padding)
        self.btn_transmute.grid(column=7, row=11, **self.padding)
        self.chkbtn_overwrite.grid(column=8, row=11, **self.padding)

        self.chkbtn_build.grid(column=1, row=12, **self.padding)


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

        # get valid key from the TABLE_NUCLEOBASE
        if self.str_nucleobase.get() == "":
            SD.print_empty_query("--nucleobase"); return

        try :
            B = self.str_nucleobase.get().upper()
            TABLE_NUCLEOBASE[B]
        except :
            B = self.str_nucleobase.get().upper()
            SD.print_invalid_key(B, TABLE_NUCLEOBASE); return
        else :
            B = self.str_nucleobase.get().upper()
            self.base = TABLE_NUCLEOBASE[B]



        # Handle empty inputs for these fields
        for string in [self.str_pdbfname, self.str_conformation, self.choice_chem, self.str_moiety] :
            if len(string.get()) == 0 or string.get() == "..." : 
                SD.print_empty_query("--pdb, --chemistry, --conformation or --nucleobase")
                return

        # Sort of handle numerical inputs
        list_ang = [ format_angle("angle(alpha)", self.str_ang_a.get()),
                     format_angle("angle(beta)", self.str_ang_b.get()),
                     format_angle("angle(gamma)", self.str_ang_g.get()),
                     format_angle("angle(delta)", self.str_ang_d.get()),
                     format_angle("angle(epsilon)", self.str_ang_e.get()),
                    ]
        list_dih = [ format_angle("dihedral(alpha)", self.str_dihr_a.get()),
                     format_angle("dihedral(beta)", self.str_dihr_b.get()),
                     format_angle("dihedral(gamma)", self.str_dihr_g.get()),
                     format_angle("dihedral(delta)", self.str_dihr_d.get()),
                     format_angle("dihedral(epsilon)", self.str_dihr_e.get()),
                    ]
        if self.int_zeta.get() == 1 :
            list_ang.append(format_angle("angle(zeta)", self.str_ang_z.get()))
            list_dih.append(format_angle("dihedral(zeta)", self.str_dihr_z.get()))
        if self.int_nu.get() == 1 :
            list_ang.append(format_angle("angle(nu)", self.str_ang_n.get()))
            list_dih.append(format_angle("dihedral(nu)", self.str_dihr_n.get()))

        list_ang.append(format_angle("angle(chi)", self.str_ang_x.get()))
        list_dih.append(format_angle("dihedral(chi)", self.str_dihr_x.get()))

        # Check if any angles are empty, if so, return early
        if "" in list_ang or "" in list_dih:
            SD.print_empty_query("Angles or Dihedrals")
            return

        # This means that when we import a file to be read in by the GUI, that we will have it start out in the current directory
        # OR the $DUCQUE/transmute directory
        input_fname = self.choice_chem.get().lower() + self.str_nucleobase.get().lower() +  "_" + self.str_conformation.get().lower() + ".tinp"


        if not self.int_overwrite.get() == 1 and isfile(input_fname) :
            SD.print_no_overwrite(input_fname, getcwd())
            return

            
        # This file has to be generated in the $DUCQUEHOME/transmute/* folder and not in the current one
        with open(input_fname , "w") as fileto :
            fileto.write("--pdb " + self.str_pdbfname.get() + "\n"
                        "--chemistry " + self.choice_chem.get() + " \n"  
                        "--conformation " + self.str_conformation.get() + " \n"
                        "--moiety " + self.str_moiety.get() + " \n"
                        "--bondangles " + ", ".join(list_ang) + " \n"
                        "--dihedrals " + ", ".join(list_dih) + " \n"
                    )

        SD.print_writing(input_fname)


    def transmute_input(self):


#        # get key from the nucleobase
#        try :
#            self.str_nucleobase.get()
#        except :
#            SD.print_empty_query("--nucleobase"); return
#
#        try :
#            B = self.str_nucleobase.get().upper()
#            TABLE_NUCLEOBASE[B]
#        except :
#            B = self.str_nucleobase.get().upper()
#            SD.print_invalid_chemistry(B); return
#        else :
#            B = self.str_nucleobase.get().upper()
#            self.base = TABLE_NUCLEOBASE[B]

        for string in [self.choice_chem , self.str_conformation, self.str_nucleobase]:
            s = string.get().lower()
            if len(s) == 0 :
                SD.print_empty_query("--chemistry, --conformation or --nucleobase")
                return

        # filename for the transmute file
        json_fname = self.choice_chem.get().lower() + "_" + self.str_nucleobase.get().lower() + "_" + self.str_conformation.get().lower() + ".json"
        input_fname = self.choice_chem.get().lower() + self.str_nucleobase.get().lower() +  "_" + self.str_conformation.get().lower() + ".tinp"


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
