from .control import Control
from .system import System
from .electrons import Electrons
from .ions import Ions
from .cell import Cell
from .fcp import FCP


class Handler:

    def __init__(self, calculation, structure, electronic_state):
        self.parameters = {
            "control": Control(calculation="scf"),
            "system": System(),
            "electrons": Electrons(),
            "ions": Ions(),
            "cell": Cell(),
            "fpc": FCP()};
        self.__dict__.update(self.parameters); #passed by reference!

    def set(self, **kwargs ):
        attributes = self.__dict__;
        for k, v in kwargs.items():
            print(k,v)
            if (k not in attributes) or  (v is  None):
                raise ValueError("key:",k," does not exists or ",  v, "is None")
            attributes.update({k:v})
            
    def __str__(self):
        out = ""
        for key, value in self.parameters.items():
            out += str(value)
        return out
