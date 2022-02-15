import qe_io 

class handler:
    
    options = dict();
    
    def __init__(self):
        self.options["ion_positions"] = 'default';

    def text(self):
        out = "&IONS\n";
        for k,v in self.options.items(): 
            out+= k+"="+qe_io.format(v)+"\n";
        out += "/";
        return out;

    def print(self):
        print(self.text());
