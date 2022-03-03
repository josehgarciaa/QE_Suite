from .control import Control
#, system, electrons, ions, cell, fcp


class Handler:

    def __init__(self, calculation, structure, electronic_state):
        self.c   = Control(calculation="scf");
#        self.s   = system.handler();
#        self.e   = electrons.handler();
#        self.ions= ions.handler();
#        self.cell= cell.handler();
#        self.fcp = fcp.handler();

