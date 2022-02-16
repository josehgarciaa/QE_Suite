from ase import Atom
from ase import Atoms

class handler:
    """
    The species class storages and handle everything relate to different atomic species in the simulation.

    Attributes:

    """
    def __init__(self):
        self.kptype  = "automatic";
        self.kpoints = [20,20,1]
        self.shifts  = [1,1,1];
        self.options = dict();
        self.set_kpoints( ((20,20,1),(1,1,1)), kptype = "automatic")

    def set_kpoints(self, kpoints, kptype = "automatic"):
        s = self;
        s.kptype = kptype;
        if kptype == "automatic":
            s.kpoints, s.shifts = kpoints
            s.options["K_POINTS"]=s.kptype+"\n";
            for k in s.kpoints:
                s.options["K_POINTS"]+=str(k)+" ";
            for k in s.shifts:
                s.options["K_POINTS"]+=str(k)+" ";

    def text(self):
        key = "\nK_POINTS ";
        return key + self.options["K_POINTS"];

    def print(self):
        print(self.text());




    
