from .. import namelists


class QEInput():
    """A generator of input files for Quantum Espresso.

    Parameters:
    ------------
    name:
        A  name to identify the system. 
    structure:
        An instance of a :py:class:`~qe_suite.Structure`
    Structure:
        An instance of a :py:class:`~qe_suite.Calculation`
    """

    def __init__(self, name="qe_suite", structure=None, calculation=None, electronic_state=None):

        self.elec_state = None
        self.structure = structure
        self.namelists = namelists.Handler(self.structure, self.structure, self.elec_state)

    def set_state(self, two_dimensional=False, insultator=False, magnetic=False):
        pass

    def set_structure(self, structure):
        pass

    def set_calculation(self, calc):
        pass

    def write(self, ofname="out"):
        print("Writing the following information")
        print(self.namelists)
        with open(ofname, 'w') as f:
            f.write(str(self.namelists));
        return 0;

        pass
