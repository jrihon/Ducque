from tk import * 


# Build an initial window that supports all other windows.

# So we can create and destroy windows from this class without having to close and open the app again
class TitleScreen(tk.Tk):


    def __init__(self):
        super().__init__()





# A class to make set presets for all grid geometries
class GridGeometry:

    window_size = "512x300"
