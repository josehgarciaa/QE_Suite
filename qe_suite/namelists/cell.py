from .namelist import Namelist


class Cell(Namelist):

    options = dict

    def __init__(self):
        self.set_parameters(
            {"cell_dynamics": None, "press": None, "wmass": None,
             "cell_factor": None, "press_conv_thr": None, "cell_dofree": None}
        )
        self.set_name("&CELL");
