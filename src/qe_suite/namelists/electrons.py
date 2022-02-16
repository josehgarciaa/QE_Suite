import qe_suite.io as qe_io 

class handler:
    
    options = dict();
    
    def __init__(self):
        self.options["diagonalization"] = 'david';
        self.options["conv_thr"] =   4e-10;
        self.options["electron_maxstep"] = 200;
        self.options["mixing_beta"] = 4e-1;

    def text(self):
        out = "&ELECTRONS\n";
        for k,v in self.options.items(): 
            out+= k+"="+qe_io.format(v)+"\n";
        out += "/\n";
        return out;

    def print(self):
        print(self.text());
