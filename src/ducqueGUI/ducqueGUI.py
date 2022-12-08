from ducqueGUI.kind_select import SelectApp

from ducqueGUI.kind_rand import RandomiseApp
from ducqueGUI.kind_build import BuildApp
from ducqueGUI.kind_xyzpdb import FormatPdbApp
from ducqueGUI.kind_transmute import TransmuteApp


def run(App, title : str):
    """ Run the tkinter GUI Application """
    app = App(title)
    app.mainloop()         # ... RUN ! DUH DUUUH DUNDUNDUDUNUDDDUUNNN




def gui_window(cli_argument):
    if cli_argument == "build" :
        run(BuildApp, cli_argument)

    if cli_argument == "randomise" :
        run(RandomiseApp, cli_argument)

    if cli_argument == "transmute" :
        run(TransmuteApp, cli_argument)

    if cli_argument == "xyz_pdb" :
        run(FormatPdbApp, cli_argument)

def select_window(cli_argument : str):

    if cli_argument == "NO_FLAG" :
        SelectApp("Choose Module").start()

    else : 
        gui_window(cli_argument)
