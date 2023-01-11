import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

from shutil import which   # Run Ducque
from subprocess import run # Run Ducque
from os import getcwd

from builder.builder_library import backbone_codex # import possibilities to build complementary strand
from dgui.grid_geometry import Geometry as G
import systemsDucque as SD

#  +--------------------------------------------------+
#  |                    BUILD                         |
#  +--------------------------------------------------+
class BuildApp(tk.Tk):

    def __init__(self, title):
        # baseline stuff
        super().__init__()
        self.title("Ducque : " + title)
        self.geometry(G.window_size_BUILD)
        self.cwd = getcwd()

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

        # place all the widgets
        self.place_widgets()

    def set_labels(self):

        # labels
        self.label_title = ttk.Label(self.content, text="Build a sequence of your choice")
        self.lab_seq = ttk.Label(self.content, text="--sequence", anchor=G.E, width=G.BUILD_label)
        self.lab_com = ttk.Label(self.content, text="--complement", anchor=G.E, width=G.BUILD_label)
        self.lab_out = ttk.Label(self.content, text="--out", anchor=G.E, width=G.BUILD_label)

    def set_checkbutton(self):

        # Checkbutton for the possibility of changing the name of the input files
        self.fname_btn = tk.IntVar()
        self.checkbtn_fname = ttk.Checkbutton(self.content, text="filename : ", command=self.toggle_fname_checkbutton, variable=self.fname_btn,
                                                onvalue=1, offvalue=0)

        self.fname_str = tk.StringVar()
        self.entry_fname = ttk.Entry(self.content, textvariable=self.fname_str, state='disabled')
        self.entry_fname.grid(column=1, row=6)

        # Checkbutton to toggle a list of possible chemistries from the `--complement`
        self.chem_btn = tk.IntVar()
        self.checkbtn_chemistries = ttk.Checkbutton(self.content, text="toggle list", variable=self.chem_btn, command=self.set_and_place_optionmenu,
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
            self.omenu_chem.configure(width=15)
            self.com_str.set('')
            self.entry_com.destroy()
            self.omenu_chem.grid(column=1, row=4)
        else :
            self.omenu_chem.destroy()
            self.com_str = tk.StringVar()
            self.entry_com = ttk.Entry(self.content, textvariable=self.com_str)
            self.entry_com.grid(column=1, row=4, **self.padding)
            

    def reveal_chemistry_keys(self):

        chemistries = list(backbone_codex.keys())
        chemistries[chemistries.index("Phosphate")] = "homo" # replace the phosphate key with the `homoduplex` key
        return chemistries

    def populate_entries(self):

        filetypes = (("All files", "*.*"), ("input-files", "*.in*"))

        select_files = filedialog.askopenfilenames(
                                        title="Import files : " + self.cwd,
                                        initialdir= self.cwd,
                                        filetypes=filetypes
                                        )

        try :
            self.file_queried = select_files[0]
        except IndexError:
            SD.print_empty_querty()
            return

        with open(self.file_queried, "r") as inputfile :
            file_content = inputfile.readlines()

        sizeline = 0
        for line in file_content :
            flag , inp = line.split(" ", maxsplit=1)

            if flag == "--sequence" : 
                self.entry_seq.configure(width=len(inp))
                self.seq_str.set(inp.strip())

                sizeline = 15 * len(inp)
                if sizeline > 470 :
                    sizeline = 470

            if flag == "--complement" :
                self.com_str.set(inp.strip())
                self.entry_com.configure(width=G.BUILD_label)

            if flag == "--out" :
                self.out_str.set(inp.strip())
                self.entry_out.configure(width=G.BUILD_label)

        self.geometry(str(sizeline + 512) + "x300")
#        self.geometry("975x300")

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

        # file dialog
        self.open_files_btn = ttk.Button(self.content, text="Import input file", command=self.populate_entries)
        self.open_files_btn.grid(column=2, row=2)

        self.lab_seq.grid(column=0, row=3)
        self.lab_com.grid(column=0, row=4)
        self.lab_out.grid(column=0, row=5)

        # buttons
        self.btn_write.grid(column=2, row=7, **self.padding)
        self.btn_build.grid(column=2, row=8, **self.padding)

        # entries
        self.entry_seq.grid(column=1, row=3, **self.padding)
        self.entry_com.grid(column=1, row=4, **self.padding)
        self.entry_out.grid(column=1, row=5, **self.padding)

        # checkbutton
        self.checkbtn_fname.grid(column=0, row=6, **self.padding)
        self.checkbtn_chemistries.grid(column=2, row=4)



    def write_inputfile(self):

        # If else clause if we imported a file
        try : 
            self.file_queried
        except : 
            if self.fname_btn.get() == 1 :
                self.outputfname = self.entry_fname.get()
            else :
                print("Consider adding a filename.")
                return
        else :
            self.outputfname = self.file_queried.split(".")[0]

        
        if self.out_str.get().strip() == "" :
            SD.print_empty_querty()
            return

        if self.fname_btn.get() == 1 :
            self.outputfname = getcwd() + "/" + self.entry_fname.get()

        # try except clause are stupid, because python is a dumb language
        # Rust for the win! 
        if self.chem_btn.get() == 1 :
            complement = self.chem_choices.get()
        else :
            complement = self.com_str.get()


        with open(self.outputfname + ".binp", "w") as fileto :
            fileto.write("--sequence " + self.seq_str.get() +
                        "\n--complement " + complement + 
                        "\n--out " + self.out_str.get() + "\n"
                    )

        SD.print_writing(self.outputfname + ".binp")


    def build_structure(self):

        # At this point, this would not be necessary, but better safe than sorry
        if not which("Ducque"): 
            SD.print_cant_find_Ducque()

        run(["Ducque", "--build", self.outputfname + ".binp"])
