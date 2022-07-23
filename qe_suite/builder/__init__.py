from .structure import Structure as Structure
from .pw_input import PWInput as PWInput
from .calculation import SCF, Bands
from .system import System
from .wannier90_input import Wannier90Input as Wannier90Input

"""MODULE generator of input files for Quantum Espresso.

    Parameters:
    ------------
    name:
        A  name to identify the system. 
    structure:
        An instance of a :py:class:`~qe_suite.Structure`
"""
