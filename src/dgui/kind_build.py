import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

from shutil import which   # Run Ducque
from subprocess import run # Run Ducque
from os import getcwd
from os.path import basename

from builder.builder_library import TABLE_BACKBONE # import possibilities to build complementary strand
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

        # set optionmenu
        self.set_optionmenu()

        # set entries
        self.set_checkbutton()

        # place all the widgets
        self.place_widgets()

    def set_labels(self):

        # labels
        self.label_title = ttk.Label(self.content, text="Build a sequence of your choice")
        self.lab_seq = ttk.Label(self.content, text="--sequence", anchor=G.E, width=G.BUILD_label)
        self.lab_com = ttk.Label(self.content, text="--complement", anchor=G.E, width=G.BUILD_label)
        self.lab_out = ttk.Label(self.content, text="--pdbname", anchor=G.E, width=G.BUILD_label)

    def set_checkbutton(self):

        # Checkbutton for the possibility of changing the name of the input files
        self.fname_btn = tk.IntVar()
        self.checkbtn_fname = ttk.Checkbutton(self.content, text="filename : ", command=self.toggle_fname_checkbutton, variable=self.fname_btn,
                                                onvalue=1, offvalue=0)

        self.fname_str = tk.StringVar()
        self.entry_fname = ttk.Entry(self.content, textvariable=self.fname_str, state='disabled')


    def toggle_fname_checkbutton(self):
        """ the fname_btn variable only has two states, which are 1 and 0 """
        if self.fname_btn.get() == 1 : 
            self.entry_fname.config(state="enabled")
        else : 
            self.entry_fname.delete(0, 'end')
            self.entry_fname.config(state="disabled")

    def set_optionmenu(self):
        """ Set the list for all the possible chemistries for the `--complement` flag"""

        self.chem_choices = tk.StringVar()
        chemlist = self.reveal_chemistry_keys()
        self.chem_choices.set("HOMO") # default value
        self.omenu_chem = tk.OptionMenu(self.content, self.chem_choices, *chemlist)
        self.omenu_chem.configure(width=16)
            

    def reveal_chemistry_keys(self):

        chemistries = list(TABLE_BACKBONE.keys())
        chemistries[chemistries.index("PHOSPHATE")] = "HOMO" # replace the phosphate key with the `homoduplex` key
        return chemistries

    def populate_entries(self):

        filetypes = (("All files", "*.*"), ("input-files", "*.binp"))

        select_files = filedialog.askopenfilenames(
                                        title="Import files : " + self.cwd,
                                        initialdir= self.cwd,
                                        filetypes=filetypes
                                        )

        try :
            self.file_queried = select_files[0]
        except IndexError:
            SD.print_empty_query("IMPORT INPUT FILE")
            return

        with open(self.file_queried, "r") as inputfile :
            file_content = inputfile.readlines()

        # Set filename as the file queried when importing 
        self.fname_btn.set(1)
        self.fname_str.set(basename(self.file_queried))
        self.entry_fname.config(state="enabled")

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
                chemlist = self.reveal_chemistry_keys()
                if inp.strip() not in chemlist: 
                    self.chem_choices.set("")
                    SD.print_invalid_chemistry(inp.strip())
                else :
                    self.chem_choices.set(inp.strip())

            if flag == "--pdbname" :
                self.pdb_str.set(inp.strip())
                self.entry_pdb.configure(width=G.BUILD_label)

        self.geometry(str(sizeline + 512) + "x300")


    def set_buttons(self):
        # buttons
        self.btn_write = ttk.Button(self.content, text="Write to file!", command=self.write_inputfile)
        self.btn_build = ttk.Button(self.content, text="Build", command=self.build_structure)

    def set_entries(self):
        # set sequence
        self.seq_str = tk.StringVar()
        self.entry_seq = ttk.Entry(self.content, textvariable=self.seq_str)
        # set output file
        self.pdb_str = tk.StringVar()
        self.entry_pdb = ttk.Entry(self.content, textvariable=self.pdb_str)


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
        self.entry_pdb.grid(column=1, row=5, **self.padding)
        self.entry_fname.grid(column=1, row=6, **self.padding)

        # optionmenu
        self.omenu_chem.grid(column=1, row=4)

        # checkbutton
        self.checkbtn_fname.grid(column=0, row=6, **self.padding)



    def write_inputfile(self):
        
        if self.pdb_str.get().strip() == "" :
            SD.print_empty_query("--pdb")
            return

        if self.fname_btn.get() == 1 :
            self.outputfname = getcwd() + "/" + self.entry_fname.get()
        else :  # if no file has been queried and the button is off, consider adding a filename
            try : 
                self.file_queried
            except : 
                if self.fname_btn.get() == 1 :
                    self.outputfname = self.entry_fname.get()
                else :
                    print(f"[EMPTY QUERY]     : Consider adding a filename")
                    return

        # try except clause are stupid, because python is a dumb language
        complement = self.chem_choices.get()


        if not self.outputfname.endswith(".binp"): 
            self.outputfname += ".binp"

        with open(self.outputfname, "w") as fileto :
            fileto.write("--sequence " + self.seq_str.get() +
                        "\n--complement " + complement + 
                        "\n--pdbname " + self.pdb_str.get() + "\n"
                    )

        try : 
            self.file_queried
        except : 
            SD.print_writing(self.outputfname)
        else : 
            if self.file_queried == self.outputfname :
                SD.print_writing(self.outputfname + ". Overwritten ..")
            else : 
                SD.print_writing(self.outputfname)


    def build_structure(self):

        try : 
            self.outputfname
        except : 
            SD.print_empty_query("IMPORT INPUT FILE")
            return 

        # At this point, this would not be necessary, but better safe than sorry
        if not which("Ducque"): 
            SD.print_cant_find_Ducque()

        run(["Ducque", "--build", self.outputfname])
