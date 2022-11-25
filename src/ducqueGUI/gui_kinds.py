import tkinter as tk
from tkinter import ttk

from shutil import which   # Run Ducque
from subprocess import run # Run Ducque

from builder.builder_library import backbone_codex # import possibilities to build complementary strand




#  +--------------------------------------------------+
#  |                    BUILD                         |
#  +--------------------------------------------------+
class BuildApp(tk.Tk):

    def __init__(self, title):
        # baseline stuff
        super().__init__()
        self.title("Ducque : " + title)
        self.geometry('512x300')

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

    def set_labels(self):

        E = "e"   # align labels to east
        # labels
        self.label_title = ttk.Label(self.content, text="Build a sequence of your choice")
        self.lab_seq = ttk.Label(self.content, text="--sequence", anchor=E, width=12)
        self.lab_com = ttk.Label(self.content, text="--complement", anchor=E, width=12)
        self.lab_out = ttk.Label(self.content, text="--out", anchor=E, width=12)

    def set_checkbutton(self):

        # Checkbutton for the possibility of changing the name of the input files
        self.fname_btn = tk.IntVar()
        self.checkbtn_fname = ttk.Checkbutton(self.content, text="filename : ", command=self.toggle_fname_checkbutton, variable=self.fname_btn,
                                                onvalue=1, offvalue=0)

        self.fname_str = tk.StringVar()
        self.entry_fname = ttk.Entry(self.content, textvariable=self.fname_str, state='disabled')
        self.entry_fname.grid(column=1, row=5)

        # Checkbutton to toggle a list of possible chemistries from the `--complement`
        self.chem_btn = tk.IntVar()
        self.checkbtn_chemistries = ttk.Checkbutton(self.content, variable=self.chem_btn, command=self.set_and_place_optionmenu,
                                                onvalue=1, offvalue=0)

    def toggle_fname_checkbutton(self):
        """ the fname_btn variable only has two states, which are 1 and 0 """
        if self.fname_btn.get() == 1 : 
            self.entry_fname.config(state="enabled")
        else : 
            self.entry_fname.delete(0, 'end')
            self.entry_fname.config(state="disabled")

    def set_and_place_optionmenu(self):
        """ Set the list for all the possible chemistries for the `--complement` flag"""

        if self.chem_btn.get() == 1 :
            self.chem_choices = tk.StringVar()
            chemlist = self.reveal_chemistry_keys()
            self.chem_choices.set("homo") # default value
            self.omenu_chem = tk.OptionMenu(self.content, self.chem_choices, *chemlist)
            self.com_str.set('')
            self.entry_com.destroy()
            self.omenu_chem.grid(column=1, row=3)
        else :
            self.omenu_chem.destroy()
            self.com_str = tk.StringVar()
            self.entry_com = ttk.Entry(self.content, textvariable=self.com_str)
            self.entry_com.grid(column=1, row=3)
            

    def reveal_chemistry_keys(self):

        chemistries = list(backbone_codex.keys())
        chemistries[chemistries.index("Phosphate")] = "homo" # replace the phosphate key with the `homoduplex` key
        return chemistries

    def set_buttons(self):
        # buttons
        self.btn_write = ttk.Button(self.content, text="write to file!", command=self.write_inputfile)
        self.btn_build = ttk.Button(self.content, text="Build", command=self.build_structure)

    def set_entries(self):
        # set sequence
        self.seq_str = tk.StringVar()
        self.entry_seq = ttk.Entry(self.content, textvariable=self.seq_str)
        # set complement
        self.com_str = tk.StringVar()
        self.entry_com = ttk.Entry(self.content, textvariable=self.com_str)
        # set output file
        self.out_str = tk.StringVar()
        self.entry_out = ttk.Entry(self.content, textvariable=self.out_str)


    def place_widgets(self):

        self.content.grid(column=0,row=0)
        # labels
        self.label_title.grid(column=1, row=1, columnspan=2)

        self.lab_seq.grid(column=0, row=2)
        self.lab_com.grid(column=0, row=3)
        self.lab_out.grid(column=0, row=4)

        # buttons
        self.btn_write.grid(column=2, row=6)
        self.btn_build.grid(column=2, row=7)

        # entries
        self.entry_seq.grid(column=1, row=2)
        self.entry_com.grid(column=1, row=3)
        self.entry_out.grid(column=1, row=4)

        # checkbutton
        self.checkbtn_fname.grid(column=0, row=5)
        self.checkbtn_chemistries.grid(column=2, row=3)


    def write_inputfile(self):

        # If else clause to check if it exists
        if self.fname_btn.get() == 1 :
            fname = self.entry_fname.get()
        else : # what to do if it does not exist
            fname = self.out_str.get()

        # try except clause are stupid, because python is a dumb language
        # Rust for the win! 
        if self.chem_btn.get() == 1 :
            complement = self.chem_choices.get()
        else :
            complement = self.com_str.get()


        with open("./" + fname + ".in", "w") as fileto :
            fileto.write("--sequence " + self.seq_str.get() +
                        "\n--complement " + complement + 
                        "\n--out " + self.out_str.get() + "\n"
                    )

        print("File written to : " + fname + ".in")


    def build_structure(self):

        # If else clause to check if it exists
        if self.fname_btn.get() == 1 :
            fname = self.entry_fname.get()
        else : # what to do if it does not exist
            fname = self.out_str.get()

        # At this point, this would not be necessary, but better safe than sorry
        if not which("Ducque"): 
            print("Ducque not found in the $PATH. Please add `Ducque` to the search path.\n")

        run(["Ducque", "--build", fname + ".in"])






#  +--------------------------------------------------+
#  |                    RANDOMISE                     |
#  +--------------------------------------------------+

class RandomiseApp(tk.Tk):

    def __init__(self, title):
        # baseline stuff
        super().__init__()
        self.title("Ducque : " + title)
        self.geometry('512x300')

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

    def set_labels(self):

        E = "e"   # align labels to east
        # labels
        self.label_title = ttk.Label(self.content, text="Randomise a sequence of your choice")
#        self.lab_len = ttk.Label(self.content, text="--length", anchor=E, width=12)
#        self.lab_seq = ttk.Label(self.content, text="--sequence", anchor=E, width=12)
        self.lab_chm = ttk.Label(self.content, text="--chemistry", anchor=E, width=12)

    def set_checkbutton(self):

        # Checkbutton for the possibility of changing the name of the input files
        self.fname_btn = tk.IntVar()
        self.checkbtn_fname = ttk.Checkbutton(self.content, text="filename : ", command=self.toggle_fname_checkbutton, variable=self.fname_btn,
                                                onvalue=1, offvalue=0)

        self.fname_str = tk.StringVar()
        self.entry_fname = ttk.Entry(self.content, textvariable=self.fname_str, state='disabled')
        self.entry_fname.grid(column=1, row=6)

        # Checkbutton to toggle a list of possible chemistries from the `--complement`
        self.compl_btn_toggle = tk.IntVar()
        self.checkbtn_compl_chems = ttk.Checkbutton(self.content, variable=self.compl_btn_toggle, command=self.set_and_place_optionmenu,
                                                onvalue=1, offvalue=0)
        self.checkbtn_compl_chems.config(state="disabled")

        # Checkbutton for complement, as this is an optional value
        self.com_toggle_int = tk.IntVar()
        self.checkbtn_com_toggle = ttk.Checkbutton(self.content, variable=self.com_toggle_int, width=12, text="--complement",
                                                            command=self.toggle_compl_checkbutton, onvalue=1, offvalue=0)

        # Toggles for sequence and for length
        self.seq_toggle_int = tk.IntVar()
        self.len_toggle_int = tk.IntVar()
        self.checkbtn_sequence = ttk.Checkbutton(self.content, variable=self.seq_toggle_int, command=self.toggle_sequence,
                                                onvalue=1, offvalue=0, text="--sequence")
        self.checkbtn_length = ttk.Checkbutton(self.content, variable=self.len_toggle_int, command=self.toggle_length,
                                                onvalue=1, offvalue=0, text="--length")


    def toggle_sequence(self):

        if self.seq_toggle_int.get() == 1 :
            self.len_toggle_int.set(0)
            self.len_str.set('')
            self.entry_len.config(state="disabled")
            self.entry_seq.config(state="enabled")

    def toggle_length(self):
        if self.len_toggle_int.get() == 1 :
            self.seq_toggle_int.set(0)
            self.seq_str.set('')
            self.entry_seq.config(state="disabled")
            self.entry_len.config(state="enabled")


    def toggle_fname_checkbutton(self):
        """ the fname_btn variable only has two states, which are 1 and 0 """
        if self.fname_btn.get() == 1 : 
            self.entry_fname.config(state="enabled")
        else : 
            self.entry_fname.delete(0, 'end')
            self.entry_fname.config(state="disabled")


    def toggle_compl_checkbutton(self):
        """ the complement checkbutton variable only has two states, which are 1 and 0 """
        if self.com_toggle_int.get() == 1 : 
            self.com_str.set('')

            if self.entry_com.winfo_exists() == 0 :
                self.entry_com = ttk.Entry(self.content, textvariable=self.com_str)

            self.entry_com.config(state="enabled")              # entry for complement enabled
            self.checkbtn_compl_chems.config(state="enabled")   # checkbutton for the optionmenu enabled
        else : 

            if self.omenu_chem.winfo_exists() == 1 :
                self.compl_btn_toggle.set(0)
                self.chem_choices.set('')
                self.omenu_chem.destroy()

            self.com_str.set('')
            self.entry_com = ttk.Entry(self.content, textvariable=self.com_str, state="disabled") # entry for complement disabled
            self.entry_com.grid(column=1, row=5)
            self.checkbtn_compl_chems.config(state="disabled")  # checkbutton for the optionmenu disabled

    def set_and_place_optionmenu(self):
        """ Set the list for all the possible chemistries for the `--complement` flag"""

        if self.compl_btn_toggle.get() == 1 :
            self.chem_choices = tk.StringVar()
            chemlist = self.reveal_chemistry_keys()
            self.chem_choices.set("homo") # default value
            self.omenu_chem = tk.OptionMenu(self.content, self.chem_choices, *chemlist)

            self.com_str.set('')
            if self.entry_com.winfo_exists() == 1 :
                self.entry_com.config(state="disabled")
                self.entry_com.destroy()        # for some reason does not always destroy ...

            self.omenu_chem.grid(column=1, row=5)
        else :
            self.chem_choices.set('')
            self.omenu_chem.destroy()

            self.entry_com = ttk.Entry(self.content, textvariable=self.com_str, state="enabled")
            self.entry_com.grid(column=1, row=5)
            

    def reveal_chemistry_keys(self):

        chemistries = list(backbone_codex.keys())
        chemistries[chemistries.index("Phosphate")] = "homo" # replace the phosphate key with the `homoduplex` key
        return chemistries
#
    def set_buttons(self):
        # buttons
        self.btn_write = ttk.Button(self.content, text="write to file!", command=self.write_inputfile)
        self.btn_rand = ttk.Button(self.content, text="Build", command=self.randomise_set_of_inputs)

    def set_entries(self):
        # set sequence
        self.seq_str = tk.StringVar()
        self.entry_seq = ttk.Entry(self.content, textvariable=self.seq_str, state="disabled")
        # set complement
        self.com_str = tk.StringVar()
        self.entry_com = ttk.Entry(self.content, textvariable=self.com_str, state="disabled")
        # set length
        self.len_str = tk.StringVar()
        self.entry_len = ttk.Entry(self.content, textvariable=self.len_str, state="enabled")
        # set chemistry
        self.chm_str = tk.StringVar()
        self.entry_chm = ttk.Entry(self.content, textvariable=self.chm_str)


    def place_widgets(self):

        self.content.grid(column=0,row=0)
        # labels
        self.label_title.grid(column=1, row=1, columnspan=2)

        self.lab_chm.grid(column=0, row=2)
#        self.lab_len.grid(column=0, row=3)
#        self.lab_seq.grid(column=0, row=4)
#        self.lab_com.grid(column=0, row=5)
        self.checkbtn_com_toggle.grid(column=0, row=5)

        # buttons
        self.btn_write.grid(column=2, row=7)
        self.btn_rand.grid(column=2, row=8)

        # entries
        self.entry_chm.grid(column=1, row=2)
        self.entry_len.grid(column=1, row=3)
        self.entry_seq.grid(column=1, row=4)
        self.entry_com.grid(column=1, row=5)

        # checkbutton
        self.checkbtn_fname.grid(column=0, row=6)
        self.checkbtn_compl_chems.grid(column=2, row=5)
        self.checkbtn_length.grid(column=0, row=3)
        self.checkbtn_sequence.grid(column=0, row=4)

        # set default
        self.len_toggle_int.set(1)


    def write_inputfile(self):

        # If else clause to check if it exists
        if self.fname_btn.get() == 1 :
            fname = self.entry_fname.get()
        else : # what to do if it does not exist
            fname = "randomised_sequence"

        # `--sequence` and `--length` are exclusive flags
        with open("./" + fname + ".in", "w") as fileto :
            fileto.write("--chemistry " + self.chm_str.get() )

            if self.seq_toggle_int.get() == 1 :
                fileto.write("\n--sequence " + self.seq_str.get() )

            if self.len_toggle_int.get() == 1 :
                fileto.write("\n--length " + self.len_str.get() )

                fileto.write("\n--complement " + self.com_str.get()

                        )

        print("File written to : " + fname + ".in")


    def randomise_set_of_inputs(self):

        # If else clause to check if it exists
        if self.fname_btn.get() == 1 :
            fname = self.entry_fname.get()
        else : # what to do if it does not exist
            fname = "randomised_sequence"

        # At this point, this would not be necessary, but better safe than sorry
        if not which("Ducque"): 
            print("Ducque not found in the $PATH. Please add `Ducque` to the search path.\n")

        run(["Ducque", "--randomise", fname + ".in"])
