from ase import Atom
from ase import Atoms


class AtomicSpecies:
    """The definition of the atomic species used in thecsimulation. 

      Attributes:
      ------------
        This class implement attributes that follow the naming convenction and with
        the same description as given in the Quantum Espresso Card: `\&ATOMIC_SPECIES <https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm1301`_
        Unless otherwise specifed, the attributes are to `None` and won't be included in the input_file which 
        will force Quantum Espresso to use its default values.

        When used within the python API, the atomic specie should be provided as dictionary with key= X and 
        value = (Mass_X, PseudoPot_X), where X the label that identifies the atomic species, Mass_X
        the atomic mass of the specie, and PseudoPot_X the name of the pseudopotential 

      Examples:
      ------------
      For instance, we could define a set of atomic species as in different ways

      >>> qe_input.cards.atomic_species = dict( {"C":(16,"C.upf") }

      Parameters:
      ------------
      structure:
          An instance of a :py:class:`~qe_suite.builder.Structure`

      Example
      ------------
      >>> # qe_input.control.etot_conv_thr	= 1e04;
    """

    def __init__(self, structure):
        self.atomic_species = None

    def set_atomic_species(self, species):
        self.species = species

    def get_atomic_species(self):
        return self.species

    def __str__(self):
        return "\n"+self.key+"\n" + self.options[self.key]

