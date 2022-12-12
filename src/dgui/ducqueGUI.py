from dgui.kind_select import SelectApp

from dgui.kind_rand import RandomiseApp
from dgui.kind_build import BuildApp
from dgui.kind_xyzpdb import FormatPdbApp
from dgui.kind_transmute import TransmuteApp


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
