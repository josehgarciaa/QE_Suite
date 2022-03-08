from email.headerregistry import ContentDispositionHeader
from .atomic_positions import AtomicPositions
from .atomic_species import AtomicSpecies
from .cell_parameters import CellParameters
from .k_points import KPoints
from .constraints import Constraints
from .occupations import Occupations
import qe_suite.io as qe_io


class Handler:

    def __init__(self, calculation, structure, electronic_state):
        pass
        self.parameters = {
            "atomic_species": AtomicSpecies(),
            "atomic_positions": AtomicPositions(),
            "k_points": KPoints(),
            "cell_parameters": CellParameters(),
            "constraints": Constraints(),
            "occupations": Constraints()
        }
        self.__dict__.update(self.parameters)

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
        #Update parameters from the initialize parameters in the class 
        attributes = (self.__dict__)
        self.update_parameters(attributes);

        out = "";
        for key, card in self.get_parameters().items():
            if card.value is not None: #Only initialize cards are shown
                out+= str(card)
        return out
