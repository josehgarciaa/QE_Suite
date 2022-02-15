import qe_io 

class handler:
    
    options = dict();
    
    def __init__(self):
        self.options["fcp_mu"] = 0;

    def text(self):
        out = "&FCP\n";
        for k,v in self.options.items(): 
            out+= k+"="+qe_io.format(v)+"\n";
        out += "/";
        return out;

    def print(self):
        print(self.text());
