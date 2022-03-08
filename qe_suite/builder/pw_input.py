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

    def __init__(self, name="qe_suite", system=None, calculation=None, electronic_state=None):
      
        self.elec_state = None
        self.structure = None
        self.namelists = namelists.Handler( self.structure, self.structure, self.elec_state)
        self.cards     = cards.Handler( self.structure, self.structure, self.elec_state)
        
        if calculation is not None:
            if calculation.valid():
                self.namelists.set( control = calculation.get_control() );
                self.cards.set( k_points = calculation.get_k_points() )

        if system is not None:
            if system.valid():
                self.namelists.set( system = system.get_system() );
                self.cards.set( 
                                atomic_species  = system.get_atomic_species(),
                                atomic_positions= system.get_atomic_positions(),
                                cell_parameters = system.get_cell_parameters()
                                )



    def set_state(self, two_dimensional=False, insultator=False, magnetic=False):
        """
        Set state function 
        """
        pass

    def set_structure(self, structure):
        print("I am setting the structure")
        return self;

    def set_calculation(self, calc):
        """
        The  cite function 
        """
        pass

    def write(self, ofname="qe_suite.out"):

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