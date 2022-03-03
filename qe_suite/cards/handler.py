from . import atomic_positions, atomic_species, kpoints, cell_params, constraints 

class Handler:

    def __init__(self, calculation, structure, electronic_state):
        print("I am handling the cards using", calculation, structure, electronic_state)
#        self.c   = control.handler();
#        self.s   = system.handler();
#        self.e   = electrons.handler();
#        self.ions= ions.handler();
#        self.cell= cell.handler();
#        self.fcp = fcp.handler();