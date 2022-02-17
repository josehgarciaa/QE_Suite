import qe_suite.io as qe_io 
class handler:
    
    options = dict();
    
    def __init__(self):
        self.options["ibrav"] = 0;
        self.options["ntyp"] = 0;
        self.options["nat"] = 0;
        self.options["occupations"] = 'smearing'
        self.options["smearing"] = 'cold'
        self.options["degauss"] = 1.46997236e-02
        self.celldm = None;

    def set_bravais_lattice(self, ibrav, celldm):
        self.options["ibrav"] = ibrav;
        self.celldm = celldm;

    def set_cutoff(self, ecutwfc, ecutrho=None):
        if ecutrho is None:
            ecutrho= 4*ecutwfc;
        self.options["ecutwfc"] = ecutwfc;
        self.options["ecutrho"] = ecutrho;

    def text(self):
        opt = self.options;
        out = "&SYSTEM\n";
        for k,v in opt.items():
            if( k !="ibrav" ): 
                out+= k+"="+qe_io.format(v)+"\n";

        k ="ibrav";
        if( opt[k] ==0 ):
            out+= k+"="+qe_io.format(opt[k])+"\n";
        else:
            out+= k+"="+qe_io.format(opt[k])+"\n";
            for i,c in enumerate(self.celldm):
                if c!=0:
                    out+= "celldm("+qe_io.format(1+i)+")="+qe_io.format(c)+"\n";

        out += "/";
        return out;

    def print(self):
        print(self.text());
