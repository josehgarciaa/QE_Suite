import qe_io 

class handler:
    
    options = dict();
    
    def __init__(self):
        self.options["diagonalization"] = 'david';
        self.options["conv_thr"] =   1.2e-9;
        self.options["electron_maxstep"] = 200;
        self.options["mixing_beta"] = 200;

    def text(self):
        out = "&ELECTRONS\n";
        for k,v in self.options.items(): 
            out+= k+"="+qe_io.format(v)+"\n";
        out += "/";
        return out;

    def print(self):
        print(self.text());
