from ase import Atom
from ase import Atoms

class handler:
    """
    The species class storages and handle everything relate to different atomic species in the simulation.

    Attributes:

    """
    def __init__(self):
        self.key ="ATOMIC_SPECIES";
        self.species = None;
        self.options = dict();
        self.options[self.key]=None;
        
    def set_atomic_species(self, species):
        self.species = species;
        self.options[self.key]="";
        for s,m,ps in self.species:
            self.options[self.key]+="{} {} {} \n".format(s,m,ps);

    def get_atomic_species(self):
        return self.species;

    def text(self):
        return "\n"+self.key+"\n" + self.options[self.key];

    def print(self):
        print(self.text());