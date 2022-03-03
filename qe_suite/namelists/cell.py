import qe_suite.io as qe_io 

class handler:
    
    options = dict();
    
    def __init__(self):
        self.options["CASE"] = 'none';

    def text(self):
        out = "&CELL\n";
        for k,v in self.options.items(): 
            out+= k+"="+qe_io.format(v)+"\n";
        out += "/";
        return out;

    def print(self):
        print(self.text());
