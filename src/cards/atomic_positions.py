from ase import Atom
from ase import Atoms

class handler:
    """
    The species class storages and handle everything relate to different atomic species in the simulation.

    Attributes:

    """
    def __init__(self):
        self.key  ="ATOMIC_POSITIONS";
        self.karg ="crystal";
        self.xyz= None;
        self.options = dict();
        self.options[self.key]="";
        
    def set_atomic_positions(self, structure , coords="crystal"):
        self.karg=coords;
        xyz = [ (s,x,y,z)  for s, (x,y,z) in zip( structure.get_chemical_symbols(), structure.get_scaled_positions()) ];
        self.xyz = xyz
        
        self.options[self.key] ="";
        for s,x,y,z in xyz:     
            self.options[self.key]+="{} {} {} {} \n".format(s,x,y,z);
 
    def text(self):
        return "\n"+self.key+" "+self.karg+"\n" + self.options[self.key];

    def print(self):
        print(self.text());





    
