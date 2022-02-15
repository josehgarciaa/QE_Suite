import qe_io 
import seekpath 
import spglib as spg
from ase import Atoms
class handler:
    
    options = dict();
    
    def __init__(self):
        self.options["ibrav"] = 0;

    def text(self):
        out = "&SYSTEM\n";
        for k,v in self.options.items(): 
            out+= k+"="+qe_io.format(v)+"\n";
        out += "/";
        return out;

    def print(self):
        print(self.text());
