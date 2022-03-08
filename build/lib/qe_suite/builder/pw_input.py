from re import M
from .. import namelists
from .. import cards


class PWInput():
    """
    A generator of input files for Quantum Espresso.

    Parameters:

        name:
            A  name to identify the system. 
        structure:
            An instance of a :py:class:`~qe_suite.builder.Structure`
        Structure:
            An instance of a :py:class:`~qe_suite.builder.Calculation`
    
    Methods:
    """

    def __init__(self, name="qe_suite", structure=None, calculation=None, electronic_state=None):

        self.elec_state = None
        self.structure = structure
        self.namelists = namelists.Handler( self.structure, self.structure, self.elec_state)
        self.cards     = cards.Handler( self.structure, self.structure, self.elec_state)
        
        if calculation is not None:
            if calculation.valid():
                self.namelists.set( control = calculation.get_control() );
                self.cards.set( 
                                atomic_species= calculation.get_atomic_species(),
                                atomic_positions = calculation.get_atomic_positions(),
                                k_points = calculation.get_k_points(),
                                cell_parameters  = calculation.get_cell_parameters()
                                )


    def set_state(self, two_dimensional=False, insultator=False, magnetic=False):
        """
        Set state function 
        """
        pass

    def set_structure(self, structure):
        """
        The  cite function 
        """
        pass

    def set_calculation(self, calc):
        """
        The  cite function 
        """
        pass

    def write(self, ofname="out"):

        with open(ofname, 'w') as f:
            f.write(str(self))
        return 0



    def cite(self):
        """
        
        The  cite function  

        """
        pass

    def __str__(self):
        out = str(self.namelists)+str(self.cards);
        return out