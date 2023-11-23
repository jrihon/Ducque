import tkinter as tk
from tkinter import ttk

from os.path import isfile
from shutil import which   # Run Ducque
from subprocess import run # Run Ducque

from builder.builder_library import TABLE_BACKBONE # import possibilities to build complementary strand
from dgui.grid_geometry import Geometry as G

import systemsDucque as SD
import dgui.kind_build as KB

#  +--------------------------------------------------+
#  |                    RANDOMISE                     |
#  +--------------------------------------------------+

class RandomiseApp(tk.Tk):

    def __init__(self, title):
        # baseline stuff
        super().__init__()
        self.title("Ducque : " + title)
        self.geometry(G.window_size_RANDOMISE)
        self.padding = {"padx" : G.padx, "pady" : G.pady}

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
        
        #set optionmenus
        self.set_optionmenu_chemistry()
        self.set_optionmenu_complement()

        # place all the widgets
        self.place_widgets()


    def set_labels(self):

        E = "e"   # align labels to east
        # labels
        self.label_title = ttk.Label(self.content, text="Randomise a sequence of your choice")
        self.lab_chm = ttk.Label(self.content, text="--chemistry : ", anchor=E, width=12)

    def set_checkbutton(self):

        # Checkbutton for the possibility of changing the name of the input files
        self.fname_btn = tk.IntVar()
        self.checkbtn_fname = ttk.Checkbutton(self.content, text="--filename : ", command=self.toggle_fname_checkbutton, variable=self.fname_btn,
                                                onvalue=1, offvalue=0, width=G.BUILD_label)

        self.fname_str = tk.StringVar()
        self.entry_fname = ttk.Entry(self.content, textvariable=self.fname_str, state='disabled', width=G.BUILD_label)


        # Checkbutton for complement, as this is an optional value
        self.com_int = tk.IntVar()
        self.checkbtn_com = ttk.Checkbutton(self.content, variable=self.com_int, width=G.BUILD_label, text="--complement : ",
                                                            command=self.toggle_compl_checkbutton, onvalue=1, offvalue=0)

        # Toggles for sequence and for length
        self.seq_toggle_int = tk.IntVar()
        self.len_toggle_int = tk.IntVar()
        self.checkbtn_sequence = ttk.Checkbutton(self.content, variable=self.seq_toggle_int, command=self.toggle_sequence,
                onvalue=1, offvalue=0, text="--sequence : ", width=G.BUILD_label)
        self.checkbtn_length = ttk.Checkbutton(self.content, variable=self.len_toggle_int, command=self.toggle_length,
                onvalue=1, offvalue=0, text="--length : ", width=G.BUILD_label)

    def set_optionmenu_chemistry(self):

        self.chem_choices = tk.StringVar()
        chemlist = self.reveal_chemistry_keys()
        chemlist.remove("HOMO")
        self.chem_choices.set("DNA") # default value
        self.omenu_chem = tk.OptionMenu(self.content, self.chem_choices, *chemlist)
        self.omenu_chem.configure(width=16)

    def set_optionmenu_complement(self):

        self.compl_choices = tk.StringVar()
        compllist = self.reveal_chemistry_keys()
        self.compl_choices.set("HOMO") # default value
        self.omenu_compl = tk.OptionMenu(self.content, self.compl_choices, *compllist)
        self.omenu_compl.configure(width=16)
        self.omenu_compl.configure(state="disabled")

    def toggle_sequence(self):

        if self.seq_toggle_int.get() == 1 :
            self.len_toggle_int.set(0)
            self.len_str.set('')
            self.entry_len.config(state="disabled")
            self.entry_seq.config(state="enabled")

        if self.seq_toggle_int.get() == 0  and self.len_toggle_int.get() == 0 :
            self.seq_toggle_int.set(1)
            print("Either `--sequence` or `--length` have to be active at any time.")
            return

    def toggle_length(self):

        if self.len_toggle_int.get() == 1 :
            self.seq_toggle_int.set(0)
            self.seq_str.set('')
            self.entry_seq.config(state="disabled")
            self.entry_len.config(state="enabled")

        if self.seq_toggle_int.get() == 0  and self.len_toggle_int.get() == 0 :
            self.len_toggle_int.set(1)
            print("Either `--sequence` or `--length` have to be active at any time.")
            return



    def toggle_fname_checkbutton(self):
        """ the fname_btn variable only has two states, which are 1 and 0 """
        if self.fname_btn.get() == 1 : 
            self.entry_fname.config(state="enabled")
        else : 
            self.entry_fname.delete(0, 'end')
            self.entry_fname.config(state="disabled")


    def toggle_compl_checkbutton(self):
        if self.com_int.get() == 1 : 
            self.omenu_compl.config(state="active")

        else :
            self.omenu_compl.config(state="disabled")


    def set_and_place_optionmenu_complement(self):
        """ Set the list for all the possible chemistries for the `--complement` flag"""

        self.compl_choices = tk.StringVar()
        compllist = self.reveal_chemistry_keys()
        self.compl_choices.set("HOMO") # default value
        self.omenu_compl = tk.OptionMenu(self.content, self.compl_choices, *compllist)
        self.omenu_compl.configure(width=16)

        self.omenu_compl.grid(column=1, row=5) # place widget

            
    def set_and_place_optionmenu_chemistry(self):
        """ Set the list for all the possible chemistries for the `--complement` flag"""

        self.chem_choices = tk.StringVar()
        chemlist = self.reveal_chemistry_keys()
        chemlist.remove("HOMO")
        self.chem_choices.set("DNA") # default value
        self.omenu_chem = tk.OptionMenu(self.content, self.chem_choices, *chemlist)
        self.omenu_chem.configure(width=16)

        self.omenu_chem.grid(column=1, row=2) # place widget

    def reveal_chemistry_keys(self):

        chemistries = list(TABLE_BACKBONE.keys())
        chemistries[chemistries.index("PHOSPHATE")] = "HOMO" # replace the phosphate key with the `homoduplex` key
        return chemistries
#
    def set_buttons(self):
        # buttons
        self.btn_write = ttk.Button(self.content, text="write to file!", command=self.write_inputfile)
        self.btn_rand = ttk.Button(self.content, text="Randomise", command=self.randomise_set_of_inputs)
        self.btn_go2build = ttk.Button(self.content, text="Go to Build module", command=self.go_to_build)

    def set_entries(self):
        # set sequence
        self.seq_str = tk.StringVar()
        self.entry_seq = ttk.Entry(self.content, textvariable=self.seq_str, state="disabled")
        # set length
        self.len_str = tk.StringVar()
        self.entry_len = ttk.Entry(self.content, textvariable=self.len_str, state="enabled")

    def go_to_build(self):
        self.destroy()
        self.app = KB.BuildApp("build")
        self.mainloop() 

    def place_widgets(self):

        self.content.grid(column=0,row=0)
        # labels
        self.label_title.grid(column=1, row=1, columnspan=2)

        self.lab_chm.grid(column=0, row=2, sticky=tk.W)

        # buttons
        self.btn_write.grid(column=2, row=8, **self.padding )
        self.btn_rand.grid(column=2, row=9, **self.padding )
        self.btn_go2build.grid(column=0, row=9, **self.padding)

        # entries
        self.entry_len.grid(column=1, row=3, **self.padding )
        self.entry_seq.grid(column=1, row=4, **self.padding )
        self.entry_fname.grid(column=1, row=6, **self.padding )

        # option menu
        self.omenu_chem.grid(column=1, row=2, **self.padding )
        self.omenu_compl.grid(column=1, row=5, **self.padding )

        # checkbutton
        self.checkbtn_length.grid(column=0, row=3, sticky=tk.W)
        self.checkbtn_sequence.grid(column=0, row=4, sticky=tk.W)
        self.checkbtn_com.grid(column=0, row=5, sticky=tk.W)
        self.checkbtn_fname.grid(column=0, row=6, sticky=tk.W)

        # set default
        self.len_toggle_int.set(1)


    def write_inputfile(self):

        # Filename
        if self.fname_btn.get() == 1 :
            fname = self.entry_fname.get()
        else : 
            fname = "randomised_sequence"

        # `--sequence` and `--length` are exclusive flags
        with open("./" + fname + ".rinp", "w") as fileto :


            ## Chemistry
            fileto.write("--chemistry " + self.chem_choices.get() )

            ## Sequence
            if self.seq_toggle_int.get() == 1 :
                fileto.write("\n--sequence " + self.seq_str.get() )

            ## Length
            if self.len_toggle_int.get() == 1 :
                fileto.write("\n--length " + self.len_str.get() )

            ## Complement
            if self.com_int.get() == 1 :
                fileto.write("\n--complement " + self.compl_choices.get().upper() )

            fileto.write("\n")

        SD.print_writing(fname + ".rinp")


    def randomise_set_of_inputs(self):

        # If else clause to check if it exists
        if self.fname_btn.get() == 1 :
            fname = self.entry_fname.get()
        else : # what to do if it does not exist
            fname = "randomised_sequence"

        # At this point, this would not be necessary, but better safe than sorry
        if not which("Ducque"): 
            SD.print_cant_find_Ducque()

        try :
            isfile(fname)
        except :
            SD.print_empty_query("randomise inputfile")
        else :
            run(["Ducque", "--randomise", fname + ".rinp"])
