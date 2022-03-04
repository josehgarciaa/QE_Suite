from .control import Control
from .system  import System
from .electrons import Electrons
from .ions import Ions
from .cell import Cell
from .fcp import FCP
class Handler:

    def __init__(self, calculation, structure, electronic_state):
        self.__dict__.update(
            {"control": Control(calculation="scf") , 
             "system": System(), 
             "electrons": Electrons(), 
             "ions": Ions(), 
             "cell": Cell(), 
             "fpc": FCP(), 
             })

    def __str__(self):
        out = "";
        for key, namelist in self.__dict__.items():
            out += str(namelist);
        return out


