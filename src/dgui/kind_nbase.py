import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

from shutil import which   # Run Ducque
from subprocess import run # Run Ducque
from os import getcwd
from os.path import isfile

from ducquelib.library import TABLE_NUCLEOBASE_MODS # import possibilities of nucleobase modifications
from dgui.grid_geometry import Geometry as G
import systemsDucque as SD

#  +--------------------------------------------------+
#  |                    BUILD                         |
#  +--------------------------------------------------+
class NbaseApp(tk.Tk):

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

        # Count amount of modification in a given file : 
        self.amountOfMods = 0

    def set_labels(self):

        # labels
        self.label_title = ttk.Label(self.content, text="Modify nucleobases of your choice")
        self.lab_pdb = ttk.Label(self.content, text="--pdb", anchor=G.E, width=G.BUILD_label)
        self.lab_nbase = ttk.Label(self.content, text="--nucleobase", anchor=G.E, width=G.BUILD_label)
        self.lab_pos = ttk.Label(self.content, text="-position", anchor=G.E, width=G.BUILD_label)
        self.lab_mod = ttk.Label(self.content, text="-mod", anchor=G.E, width=G.BUILD_label)
        self.lab_resname = ttk.Label(self.content, text="-resname", anchor=G.E, width=G.BUILD_label)

        self.lab_from = ttk.Label(self.content, text="from", anchor=G.C, width=G.BUILD_label)
        self.lab_to = ttk.Label(self.content, text="to", anchor=G.C, width=G.BUILD_label)

        self.lab_residueNumber = ttk.Label(self.content, text="residue number", anchor=G.C, width=G.BUILD_label)
        self.lab_chainLetter = ttk.Label(self.content, text="chain", anchor=G.C, width=G.BUILD_label)
#        self.lab_reorient = ttk.Label(self.content, text="-reorient", anchor=G.E, width=G.BUILD_label)

        self.lab_line = ttk.Label(self.content, text="---------------------------------", width=G.BUILD_label)

        self.lab_output_fname = ttk.Label(self.content, text="filename", anchor=G.E, width=G.BUILD_label)

        self.lab_empty = ttk.Label(self.content, text="    ", anchor=G.E, width=G.BUILD_label)

        self.lab_counter = ttk.Label(self.content, text="Mod. Counter :", anchor=G.E, width=G.BUILD_label)
        self.lab_counter_value = ttk.Label(self.content, text="0", anchor=G.C, width=G.BUILD_label)

    def set_checkbutton(self):

        # Checkbutton for the possibility of changing the name of the input files
        self.btn_reorient = tk.IntVar()
        self.str_reorient = tk.StringVar()
        self.btn_reorient.set(0)

        self.lab_reorient = ttk.Label(self.content, text="False")

        self.lab_reorient.grid(column=1, row=10, **self.padding)

        self.checkbtn_reorient = ttk.Checkbutton(self.content, text="-reorient", command=self.toggle_reorient_checkbutton, variable=self.btn_reorient, onvalue=1, offvalue=0)


    def toggle_reorient_checkbutton(self):

        if self.btn_reorient.get() == 1 : 
            self.lab_reorient.destroy()
            self.lab_reorient = ttk.Label(self.content, text="True")
            self.lab_reorient.grid(column=1, row=10, **self.padding)
        else : 
            self.lab_reorient.destroy()
            self.lab_reorient = ttk.Label(self.content, text="False")
            self.lab_reorient.grid(column=1, row=10, **self.padding)


    def set_buttons(self):
        # buttons
        self.btn_fname = ttk.Button(self.content, text="Prompt pdb file!", command=self.prompt_pdbfile)
        self.btn_write = ttk.Button(self.content, text="Write to file!", command=self.write_inputfile)
        self.btn_modify = ttk.Button(self.content, text="Modify", command=self.modify_nbases)
        self.btn_ninp  = ttk.Button(self.content, text="Prompt ninp file!", command=self.prompt_ninpfile)

    def set_entries(self):
        # set input/ouput file
        self.str_pdb = tk.StringVar()
        self.entry_pdb = ttk.Entry(self.content, textvariable=self.str_pdb)

        self.str_outputFname = tk.StringVar()
        self.entry_outputFname = ttk.Entry(self.content, textvariable=self.str_outputFname)

        # set position file
        self.str_residueNumber = tk.StringVar()
        self.str_chainLetter = tk.StringVar()
        self.entry_residueNumber = ttk.Entry(self.content, textvariable=self.str_residueNumber)
        self.entry_chainLetter = ttk.Entry(self.content, textvariable=self.str_chainLetter)

        # set residue name
        self.str_residueName = tk.StringVar()
        self.entry_residueName = ttk.Entry(self.content, textvariable=self.str_residueName)

    def set_optionmenu(self):
        """ Set the list for all the possible chemistries for the `--complement` flag"""

        self.from_modification = tk.StringVar()
        self.to_modification = tk.StringVar()
        from_list = self.reveal_nbases_keys()
        to_list = self.reveal_nbases_keys()

        # Set two 
        self.omenu_from = tk.OptionMenu(self.content, self.from_modification, *from_list)
        self.omenu_to = tk.OptionMenu(self.content, self.to_modification, *to_list)

        self.from_modification.set(from_list[0])
        self.to_modification.set(to_list[0])

        self.omenu_from.configure(width=16)
        self.omenu_to.configure(width=16)
            
    def reveal_nbases_keys(self):

        return list(TABLE_NUCLEOBASE_MODS.keys())


    def prompt_pdbfile(self):

        filetypes = (("All files", "*.*"), ("input-files", "*.pdb"))

        select_files = filedialog.askopenfilenames(
                                        title="Import files : " + self.cwd,
                                        initialdir= self.cwd,
                                        filetypes=filetypes
                                        )

        try :
            self.file_queried_pdb = select_files[0]
        except IndexError:
            SD.print_empty_query("IMPORT PDB FILE")
            return

        self.str_pdb.set(self.file_queried_pdb)


    def prompt_ninpfile(self):

        filetypes = (("All files", "*.*"), ("input-files", "*.ninp"))

        select_files = filedialog.askopenfilenames(
                                        title="Import files : " + self.cwd,
                                        initialdir= self.cwd,
                                        filetypes=filetypes
                                        )

        try :
            self.file_queried_ninp = select_files[0]
        except IndexError:
            SD.print_empty_query("IMPORT NINP FILE")
            return

        # Check if there is a `--pdb` flag in the file already
        foundPdbFlag = False
        with open(self.file_queried_ninp, "r") as content : 
            lines = content.readlines()
            for line in lines : 

                # Match on a line that starts with the `--pdb` flag
                if line.startswith("--pdb") : 

                    try : 
                        _ = line.split()[1]
                    except IndexError:
                        SD.print_empty_query("FILE CONTAINS EMPTY QUERY FOR `--pdb")
                        return
                    else : 
                        containedPdbEntry = line.split()[1]
                        matchesPdbEntry = self.str_pdb.get() == containedPdbEntry
                    
                    # Check if the currently instanced pdb entry matches the one in the file 
                    hasPdbEntry = self.str_pdb.get() != ""

                    # If the file has a `--pdb` entry, check if it matches with the current one
                    if matchesPdbEntry : 
#                        print("The prompted `.ninp` file matches the `--pdb` flag!")
                        self.str_pdb.set(containedPdbEntry)
                        self.str_outputFname.set(self.file_queried_ninp)
                        foundPdbFlag = True
                        break
                    else : 
                        # If we found a `--pdb` in the file and the entry is empty
                        if len(containedPdbEntry) > 0 : 
                            self.str_pdb.set(containedPdbEntry)
                            self.str_outputFname.set(self.file_queried_ninp)
                            foundPdbFlag = True
                        else : 
                            SD.print_mismatch_flag(self.file_queried_ninp, "--pdb")
                            return

                    # if --pdb flag is set nonetheless
                    if hasPdbEntry : 
                        self.str_pdb.set(containedPdbEntry)
                        self.str_outputFname.set(self.file_queried_ninp)
                        foundPdbFlag = True
                        break 

            if not foundPdbFlag : 
                SD.print_invalid_file(self.file_queried_ninp)
                return


#                if line.startswith("--pdb") : 
#                    
#                    try : 
#                        _ = line.split()[1]
#                    except IndexError:
#                        SD.print_empty_query("FILE CONTAINS EMPTY QUERY FOR `--pdb")
#                        return
#                    else : 
#                        containedPdbEntry = line.split()[1]
#                        matchesPdbEntry = self.str_pdb.get() == containedPdbEntry
#
#                    # if --pdb flag is empty
#                    if self.str_pdb.get() == "" : 
#                        self.str_pdb.set(containedPdbEntry)
#                        self.str_outputFname.set(self.file_queried_ninp)
#                        break 
#
#                    # Check if the currently instanced pdb entry matches the one in the file 
#                    hasPdbEntry = self.str_pdb.get() != ""
#
#                    # If the file has a `--pdb` entry, check if it matches with the current one
#                    if not hasPdbEntry and not matchesPdbEntry : 
#                    if hasPdbEntry and matchesPdbEntry : 
#                        print("The prompted `.ninp` file matches the `--pdb` flag!")
#                    else : 
#                        SD.print_mismatch_flag(self.file_queried_ninp, "--pdb")
#                        return
#

            # Count the instances of the word `--nucleobase`
            amountOfMods = 0
            for ljne in lines : 

                if ljne.startswith("--nucleobase") : 
                    amountOfMods += 1

            self.amountOfMods = amountOfMods

            self.lab_counter_value.destroy()
            self.lab_counter_value = ttk.Label(self.content, text=str(self.amountOfMods))
            self.lab_counter_value.grid(column=1, row=13, **self.padding)

        self.str_outputFname.set(self.file_queried_ninp)


    def place_widgets(self):

        self.content.grid(column=0,row=0)
        # labels
        self.label_title.grid(column=1, row=1, columnspan=2)

        # file dialog
        self.lab_line.grid(column=1, row=4, columnspan=2)

        self.lab_pdb.grid(column=0, row=3)
        self.lab_nbase.grid(column=0, row=4)
        self.lab_pos.grid(column=0, row=6)
        self.lab_mod.grid(column=0, row=8)
        self.lab_resname.grid(column=0, row=9)
        self.lab_empty.grid(column=0, row=11)
        self.lab_output_fname.grid(column=0, row=12)


        # buttons
        self.btn_fname.grid(column=2, row=3, **self.padding)
        self.btn_write.grid(column=3, row=12, **self.padding)
        self.btn_modify.grid(column=3, row=13, **self.padding)
        self.btn_ninp.grid(column=2, row=12, **self.padding)

        # nucleobase counter 
        self.lab_counter.grid(column=0, row=13, **self.padding)
        self.lab_counter_value.grid(column=1, row=13, **self.padding)

        # entries
        self.entry_pdb.grid(column=1, row=3, **self.padding)
        self.entry_residueName.grid(column=1, row=9, **self.padding)

        # position arguments
        self.lab_residueNumber.grid(column=1, row=5)
        self.lab_chainLetter.grid(column=2, row=5)
        self.entry_residueNumber.grid(column=1, row=6)
        self.entry_chainLetter.grid(column=2, row=6)
        self.entry_outputFname.grid(column=1, row=12)
#        self.entry_fname.grid(column=1, row=6, **self.padding)

        # optionmenu from and to modification
        self.lab_from.grid(column=1, row=7)
        self.lab_to.grid(column=2, row=7)
        self.omenu_from.grid(column=1, row=8)
        self.omenu_to.grid(column=2, row=8)

        # checkbutton
        self.checkbtn_reorient.grid(column=0, row=10, **self.padding)
#        self.lab_reorient.grid(column=1, row=9, **self.padding)

    def write_inputfile(self):

        # Check pdb filename
        if not isfile(self.str_pdb.get()) : 
            SD.print_filenotfound(self.str_pdb.get())
            return

        if self.str_pdb.get().strip() == "" :
            SD.print_empty_query("--pdb")
            return

        # Output filename
        if self.str_outputFname.get().strip() == "" :
            SD.print_empty_query("filename")
            return

        filenameOutput = self.str_outputFname.get()
        if not filenameOutput.endswith(".ninp"):
            filenameOutput += ".ninp"

        # Residue name
        if len(self.str_residueName.get()) > 3 : 
            SD.print_invalid_query("-resname")
            return

        # TODO : 
        # empty position 
        if self.str_residueNumber.get() == "" : 
            SD.print_empty_query("residue number")
            return
        
        try : 
            _ = int(self.str_residueNumber.get())
        except : 
            SD.print_invalid_query("residue number")
            return
        
        # TODO : 
        # Ok so something is off 
        # I want to get the following 
        #
        # If the existing ninp file contains the flag pdb, we should : 
        # -> Check if it matches the existing entry
        # -> If the existing entry is empty, then fill it out with what is in the file
        # -> If there is not --pdb entry in the given ninp file, then it is invalid!
        #
        # => Stepping away for now, making a GUI is hella dumb
        #
        #
        #
        #
        #
        #
        #
        #
        #
        #

        # If the outputfile already exists
        if isfile(filenameOutput) : 

            # Check if there is a `--pdb` flag in the file already
            with open(filenameOutput, "r") as content : 
                lines = content.readlines()

            foundPdbFlag = False
            for line in lines : 

                # Match on a line that starts with the `--pdb` flag
                if line.startswith("--pdb") : 

                    try : 
                        _ = line.split()[1]
                    except IndexError:
                        SD.print_empty_query("FILE CONTAINS EMPTY QUERY FOR `--pdb")
                        return
                    else : 
                        containedPdbEntry = line.split()[1]
                        matchesPdbEntry = self.str_pdb.get() == containedPdbEntry
                    
                    # Check if the currently instanced pdb entry matches the one in the file 
                    hasPdbEntry = self.str_pdb.get() != ""

                    # If the file has a `--pdb` entry, check if it matches with the current one
                    if matchesPdbEntry : 
#                        print("The prompted `.ninp` file matches the `--pdb` flag!")
                        self.str_pdb.set(containedPdbEntry)
                        self.str_outputFname.set(filenameOutput)
                        foundPdbFlag = True
                        break
                    else : 
                        SD.print_mismatch_flag(self.file_queried_ninp, "--pdb")
                        return

                    # if --pdb flag is set nonetheless
                    if hasPdbEntry : 
                        self.str_pdb.set(containedPdbEntry)
                        self.str_outputFname.set(self.file_queried_ninp)
                        foundPdbFlag = True
                        break 

            if not foundPdbFlag : 
                SD.print_invalid_file(filenameOutput)
                return


#            # Count the instances of the word `--nucleobase`
#            for ljne in lines : 
#
#                if ljne.startswith("--nucleobase") : 
#                    self.amountOfMods += 1
#
#            self.lab_counter_value.destroy()
#            self.lab_counter_value = ttk.Label(self.content, text=str(self.amountOfMods))
#            self.lab_counter_value.grid(column=1, row=13, **self.padding)

            # Write out to a file
            with open(filenameOutput, "r") as content : 
                lines = content.readlines()

            with open(filenameOutput, "w") as fileto : 

                for lkne in lines : 
                    fileto.write(lkne)

                # Write nucleobase line
                fileto.write(f"\n--nucleobase\n")
                # write position
                if self.str_chainLetter.get() == "": 
                    fileto.write(f"-position {self.str_residueNumber.get()}\n")
                else : 
                    fileto.write(f"-position {self.str_residueNumber.get()}.{self.str_chainLetter.get()}\n")
                # write from to modification
                fileto.write(f"-mod {self.from_modification.get()},{self.to_modification.get()}\n")
                # write new residue name
                fileto.write(f"-resname {self.str_residueName.get()}\n")
                # write reorient
                if self.btn_reorient.get() == 1 : 
                    fileto.write("-reorient\n")

            # Getting here means a succesful write has occured
            self.amountOfMods += 1
            self.lab_counter_value.destroy()
            self.lab_counter_value = ttk.Label(self.content, text=str(self.amountOfMods))
            self.lab_counter_value.grid(column=1, row=13, **self.padding)

        else : 
        
            with open(filenameOutput, "w") as fileto : 

                fileto.write(f"--pdb {self.str_pdb.get()}\n \n")
                
                # Write nucleobase line
                fileto.write("\n--nucleobase\n")
                # write position
                if self.str_chainLetter.get() == "": 
                    fileto.write(f"-position {self.str_residueNumber.get()}\n")
                else : 
                    fileto.write(f"-position {self.str_residueNumber.get()}.{self.str_chainLetter.get()}\n")
                # write from to modification
                fileto.write(f"-mod {self.from_modification.get()},{self.to_modification.get()}\n")
                # write new residue name
                fileto.write(f"-resname {self.str_residueName.get()}\n")
                # write reorient
                if self.btn_reorient.get() == 1 : 
                    fileto.write("-reorient\n")

            # Getting here means a succesful write has occured
            self.amountOfMods += 1
            self.lab_counter_value.destroy()
            self.lab_counter_value = ttk.Label(self.content, text=str(self.amountOfMods))
            self.lab_counter_value.grid(column=1, row=13, **self.padding)

        SD.print_writing(filenameOutput)

    def modify_nbases(self):

        try : 
            self.str_outputFname.get()
        except : 
            SD.print_empty_query("NO OUTPUT FILE GIVEN")
            return 

        # At this point, this would not be necessary, but better safe than sorry
        if not which("Ducque"): 
            SD.print_cant_find_Ducque()

        run(["Ducque", "--nbase", self.str_outputFname.get()])
