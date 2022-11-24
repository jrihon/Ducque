from ducqueGUI.gui_kinds import BuildApp


def run(title):

    buildapp = BuildApp(title)
    buildapp.mainloop()         # ... RUN DUH DUUUH DUNDUNDUDUNUDDDUUNNN




def gui_window(cli_argument):
    if cli_argument == "build" :
        run(cli_argument)

#    if cli_argument == "transmute" :
#        run(cli_argument)
#
#    if cli_argument == "randomise" :
#        run(cli_argument)
#
#    if cli_argument == "xyz_pdb" :
#        run(cli_argument)
#
#



#class App(tk.Frame):
#    def __init__(self, master):
#        super().__init__(master)
#        self.pack()
#
#        self.entrythingy = tk.Entry()
#        self.entrythingy.pack()
#
#        # Create the application variable.
#        self.contents = tk.StringVar()
#        # Set it to some value.
#        self.contents.set("this is a variable")
#        # Tell the entry widget to watch this variable.
#        self.entrythingy["textvariable"] = self.contents
#
#        # Define a callback for when the user hits return.
#        # It prints the current value of the variable.
#        self.entrythingy.bind('<Key-Return>',
#                             self.print_contents)
#
#    def print_contents(self, event):
#        print("Hi. The current entry content is:",
#              self.contents.get())
#
