from .control import Control
from .system import System
from .electrons import Electrons
from .ions import Ions
from .cell import Cell
from .fcp import FCP
import qe_suite.io as qe_io

class Handler:

    def __init__(self, calculation, structure, electronic_state):
        self.parameters = {
            "control": Control(calculation="scf"),
            "system": System(),
            "electrons": Electrons(),
            "ions": Ions(),
            "cell": Cell(),
            "fpc": FCP()};
        self.__dict__.update(self.parameters); #passed 

    def set_parameters(self, parameters):
        self.parameters = parameters;

    def update_parameters(self, parameters):
        for key in self.parameters.keys():
            if key in parameters:
                self.parameters.update({key:parameters[key]})

    def get_parameters(self):
        return self.parameters;

    def set(self, **kwargs ):
        attributes = self.__dict__;
        for k, v in kwargs.items():
            if (k not in attributes) or  (v is  None):
                raise ValueError("key:",k," does not exists or ",  v, "is None")
            attributes.update({k:v})

    def __str__(self):
        attributes = (self.__dict__)
        self.update_parameters(self.__dict__);        
        out = ""
        for key, value in self.parameters.items():
            out += str(value)
        return out
