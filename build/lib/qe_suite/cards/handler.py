from .atomic_positions import AtomicPositions
from .atomic_species import AtomicSpecies
from .cell_parameters import CellParameters
from .k_points import KPoints
from .constraints import Constraints
from .occupations import Occupations


class Handler:

    def __init__(self, calculation, structure, electronic_state):
        pass
        self.parameters = {
                        "atomic_species": AtomicSpecies(),
                        "atomic_positions": AtomicPositions(),
                        "k_points": KPoints(),
                        "cell_parameters": CellParameters(),
            #            "constraints": Constraints(),
            #            "occupations": Occupations()
        }
        self.__dict__.update(self.parameters)

    def __str__(self):
        out = ""
        for key, value in self.parameters.items():
            out += str(value)
        return out
