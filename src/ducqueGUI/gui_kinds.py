import tkinter as tk
from tkinter import ttk


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
        self.fname_btn = tk.IntVar()
        self.checkbtn_fname = ttk.Checkbutton(self.content, text="filename : ", command=self.toggle_checkbutton, variable=self.fname_btn,
                                                onvalue=1, offvalue=0)

        self.fname_str = tk.StringVar()
        self.entry_fname = ttk.Entry(self.content, textvariable=self.fname_str, state='disabled')
        self.entry_fname.grid(column=1, row=5)

    def toggle_checkbutton(self):
        """ the fname_btn variable only has two states, which are 1 and 0 """
        if self.fname_btn.get() == 1 : 
            self.entry_fname.config(state="enabled")
        else : 
            self.entry_fname.delete(0, 'end')
            self.entry_fname.config(state="disabled")

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


    def write_inputfile(self):

        if len(self.entry_fname.get()) > 0 :
            fname = self.entry_fname.get()
        else : # what to do if it does not exist
            fname = self.out_str.get()

        with open("./" + fname + ".in", "w") as fileto :
            fileto.write("--sequence " + self.seq_str.get() +
                        "\n--complement " + self.com_str.get() + 
                        "\n--out " + self.out_str.get() + "\n"
                    )

        print("File written to : " + fname + ".in")


    def build_structure(self):
        from shutil import which
        from subprocess import run

        if not which("Ducque"): 
            print("Ducque not found in the $PATH. Please add `Ducque` to the search path.\n")

        run(["Ducque", "--build", self.out_str.get() + ".in"])
