from re import M
from .. import namelists
from .. import cards
from ..builder import calculation as calc_type

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
                self.cards.set( k_points = calculation.get_k_points() );
                if isinstance( calculation, calc_type.Relaxation ):
                    self.namelists.set( ions = calculation.get_ions() );
                    self.namelists.set( cell = calculation.get_cell() );

        if system is not None:
            if system.valid():
                self.namelists.set( system = system.get_system() );
                if system.get_atomic_species() is not None:
                    self.cards.set( atomic_species = system.get_atomic_species() )
                if system.get_atomic_positions() is not None:
                    self.cards.set( atomic_positions=system.get_atomic_positions() )
                if system.get_cell_parameters() is not None:
                    self.cards.set( cell_parameters= system.get_cell_parameters() )

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
        return ofname



    def cite(self):
        """
        
        The  cite function  

        """
        pass

    def __str__(self):
        out = str(self.namelists)+str(self.cards);
        return out