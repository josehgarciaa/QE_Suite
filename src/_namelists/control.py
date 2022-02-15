
class handler:
  options = dict();
  def __init__(self):
    self.options["calculation"]="scf";
    self.options["nstep"]=1;
    self.options["outdir"]= "./qesuite_out" ;
    self.options["etot_conv_thr"]=1e-4;
    self.options["forc_conv_thr"]=1e-3;
    self.options["pseudo_dir"]="./";
    