import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

from os import getcwd
#from shutil import which   # Run Ducque
#from subprocess import run # Run Ducque

from builder.builder_library import backbone_codex # import possibilities to build complementary strand
#from ducqueGUI.grid_geometry import GridGeometry as GG

#  +--------------------------------------------------+
#  |                    TRANSMUTE                     |
#  +--------------------------------------------------+

class TransmuteApp(tk.Tk):

    def __init__(self, title):
        # baseline stuff
        super().__init__()
        self.title("Ducque : " + title)
        self.geometry("1500x512")
        self.cwd = getcwd()

        # Set Parent Frame
        self.content = ttk.Frame(self)

        # set labels
        self.set_labels()

        # set buttons
#        self.set_buttons()

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

        self.label_alpha = ttk.Label(self.content, text="α")
        self.label_beta = ttk.Label(self.content, text="β")
        self.label_gamma = ttk.Label(self.content, text="γ")
        self.label_delta = ttk.Label(self.content, text="δ")
        self.label_epsilon = ttk.Label(self.content, text="ε")
        self.label_zeta = ttk.Label(self.content, text="ζ")
        self.label_nu = ttk.Label(self.content, text="η")
        self.label_chi = ttk.Label(self.content, text="χ")

    def toggle_zeta(self):
        if self.int_zeta.get() == 1 : 
            self.ent_ang_z.config(state="enabled")
            self.ent_dihr_z.config(state="enabled")
        else :
            self.ent_ang_z.config(state="disabled")
            self.ent_dihr_z.config(state="disabled")

    def toggle_nu(self):
        if self.int_nu.get() == 1 : 
            self.ent_ang_n.config(state="enabled")
            self.ent_dihr_n.config(state="enabled")
        else :
            self.ent_ang_n.config(state="disabled")
            self.ent_dihr_n.config(state="disabled")



    def set_checkbutton(self):
        self.int_zeta = tk.IntVar()
        self.int_nu = tk.IntVar()
        self.chkbtn_zeta = ttk.Checkbutton(self.content, text="ζ", variable=self.int_zeta, command=self.toggle_zeta,
                                                onvalue=1, offvalue=0)
        self.chkbtn_nu = ttk.Checkbutton(self.content, text="η", variable=self.int_nu, command=self.toggle_nu,
                                                onvalue=1, offvalue=0)

        self.int_zeta.set(1)
        self.ent_ang_n.config(state="disabled")
        self.ent_dihr_n.config(state="disabled")

    def set_entries(self):
        self.str_pdbfname = tk.StringVar()
        self.str_conformation = tk.StringVar()
        
        self.entr_pdbfname = ttk.Entry(self.content, textvariable=self.str_pdbfname)
        self.entr_conformation = ttk.Entry(self.content, textvariable=self.str_conformation)

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
    
#    def populate_angles(self, angles_list):
#
#        ln = len(angles_list)

    def set_optionmenu(self):
        # chemistries
        self.choice_chem = tk.StringVar()
        opt_chems = list(backbone_codex.keys())
        self.opt_chems = ["..."] + sorted([x for x in opt_chems if x != "Phosphate"], key=str.casefold) # replace the phosphate key with the `...` key and sort
        self.omenu_chem = ttk.OptionMenu(self.content, self.choice_chem, *self.opt_chems)

        # moiety
        self.choice_moi = tk.StringVar()
        self.opt_moi = ["...", "nucleoside", "linker"]
        self.omenu_moi= ttk.OptionMenu(self.content, self.choice_moi, *self.opt_moi)

    def file_dialog(self):

        filetypes = (("All files", "*.*"), ("input-files", "*.in*"))

        select_files = filedialog.askopenfilenames(
                                        title="Import files : " + self.cwd,
                                        initialdir= self.cwd,
                                        filetypes=filetypes
                                        )

        try :
            file_queried = select_files[0]
        except IndexError:
            print("No file selected. Please try again")
            return

        with open(file_queried, "r") as inputfile :
            file_content = inputfile.readlines()

        for line in file_content :
            flag , inp = line.split(" ", maxsplit=1)

            if flag == "--pdb" : 
                self.entr_pdbfname.configure(width=len(inp))
                self.str_pdbfname.set(inp.strip())

            if flag == "--id" :
                if inp.strip() in list(backbone_codex.keys()) : self.choice_chem.set(inp.strip())
                else : print(f"{inp} is not an available type for `--id`. ")

            if flag == "--conformation" :
                self.str_conformation.set(inp.strip())
            
            if flag == "--moiety" :
                if inp.strip() in ["nucleoside", "linker"]: self.choice_moi.set(inp.strip())
                else : print(f"{inp} is not an available type for `--moiety`. ")

            if flag == "--bondangles" : pass
            if flag == "--dihedrals" : pass


#        self.set_flag_labels()
#        self.set_flag_entries(file_queried)

    def place_widgets(self):

        self.content.grid(column=0, row=0)

        self.label_words = ttk.Label(self.content, text="Convert your `pdb` to an input for Ducque")
        self.label_words.grid(column=1, row=1, columnspan=2)
        # set labels
        self.label_pdbfname.grid(column=0, row=2)
        self.label_chemistry.grid(column=0, row=3)
        self.label_conformation.grid(column=0, row=4)
        self.label_moiety.grid(column=0, row=5)

        # alphabet
        self.label_alpha.grid(column=1, row=6)
        self.label_beta.grid(column=2, row=6)
        self.label_gamma.grid(column=3, row=6)
        self.label_delta.grid(column=4, row=6)
        self.label_epsilon.grid(column=5, row=6)
        self.chkbtn_zeta.grid(column=6, row=6)
        self.chkbtn_nu.grid(column=7, row=6)
        self.label_chi.grid(column=8, row=6)

        self.label_bondangles.grid(column=0, row=7)
        self.label_dihedrals.grid(column=0, row=8)

        # bondangles
        self.ent_ang_a.grid(column=1, row=7)
        self.ent_ang_b.grid(column=2, row=7)
        self.ent_ang_g.grid(column=3, row=7)
        self.ent_ang_d.grid(column=4, row=7)
        self.ent_ang_e.grid(column=5, row=7)
        self.ent_ang_z.grid(column=6, row=7)
        self.ent_ang_n.grid(column=7, row=7)
        self.ent_ang_x.grid(column=8, row=7)

        # dihedrals
        self.ent_dihr_a.grid(column=1, row=8)
        self.ent_dihr_b.grid(column=2, row=8)
        self.ent_dihr_g.grid(column=3, row=8)
        self.ent_dihr_d.grid(column=4, row=8)
        self.ent_dihr_e.grid(column=5, row=8)
        self.ent_dihr_z.grid(column=6, row=8)
        self.ent_dihr_n.grid(column=7, row=8)
        self.ent_dihr_x.grid(column=8, row=8)

        # set entries
        self.entr_pdbfname.grid(column=1, row=2)
        self.entr_conformation.grid(column=1, row=4)

        # set optionmenu
        self.omenu_chem.grid(column=1, row=3)
        self.omenu_moi.grid(column=1, row=5)

        # set filedialog
        self.open_files_btn = ttk.Button(self.content, text="Import input file", command=self.file_dialog)
        self.open_files_btn.grid(column=2, row=10)
