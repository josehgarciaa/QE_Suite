import qe_io 

class handler:
  options = dict();
  def __init__(self):
    self.options["calculation"]="scf";
    self.options["outdir"]= "./qesuite_out" ;
    self.options["etot_conv_thr"]=1e-4;
    self.options["forc_conv_thr"]=1e-3;
    self.options["pseudo_dir"]="./";
    self.options["lfcp"] = False;
  def text(self):
    out = "&CONTROL\n";
    for k,v in self.options.items(): 
      out+= k+"="+qe_io.format(v)+"\n";
    out += "/";
    return out;

  def print(self):
    print(self.text());
