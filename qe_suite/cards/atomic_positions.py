from .card import Card

class AtomicPositions(Card):
    """The definition of the atomic positions used in thecsimulation. 

      Attributes:
      ------------
        This class implement attributes that follow the naming convenction and with
        the same description as given in the Quantum Espresso Card: `\&ATOMIC_POSITIONS <https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm1311`_
        Unless otherwise specifed, the attributes are to `None` and won't be included in the input_file which 
        will force Quantum Espresso to use its default values.

        When used within the python API, the atomic position should be provided as a list of tuples:
         [ (X_0,x_0,y_0,z_0), (X_1,x_1,y_1,z_1), ...,(X_n,x_n,y_n,z_n)] where X_i the label that identifies the position of
         the ith atom that defines the crystalline structure and x_i,y_i,z_i its cartesian position. n here should be equal to number of atoms
         nat 

      Parameters:
      ------------
      option:
        A string defining the units of the atomic position: (bohr, angstrom, crystal, crystal_sg)
      structure:
        An instance of a :py:class:`~qe_suite.builder.Structure`

      Example
      ------------
      >>> # qe_input.cards.atomic_positions.option	= "bohr";
    """
    def __init__(self):
        self.set_option("angstrom");
        self.set_value( value=None, type=dict );
        self.set_name( "ATOMIC_POSITIONS" );
        
    def set_atomic_positions(self, structure , coords="crystal"):
        self.karg=coords;
        xyz = [ (s,x,y,z)  for s, (x,y,z) in zip( structure.get_chemical_symbols(), structure.get_scaled_positions()) ];
        self.xyz = xyz
        
        self.options[self.key] ="";
        for s,x,y,z in xyz:     
            self.options[self.key]+="{} {} {} {} \n".format(s,x,y,z);
 
    def __str__(self):
      out=self.header();
      for v in self.get_value():
        out+="{}\t {} {} {}\n".format(*v)
      return out





    
