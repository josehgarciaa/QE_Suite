import qe_suite.io as qe_io
import qe_suite.io.wannier90 as wann_io
from qe_suite.parse.xml import WannierInput
import subprocess
import numpy as np


class Wannier90Input():
    """
    A generator of input files for Wannier90.

    Parameters:

        name:
            A  name to identify the system. 
        structure:
            An instance of a :py:class:`~qe_suite.builder.Structure`
        Structure:
            An instance of a :py:class:`~qe_suite.builder.Calculation`
    
    Methods:
    """

    def __init__(self, num_wann,mp_grid, use_spin=False, seedname="qe_suite.wannier", xml=None):
        self.seedname =  seedname ;
        self.set_parameters({"num_wann": None, "num_bands": None,"unit_cell_cart": None,"atoms_frac":None,"mp_grid":None,
                             "kpoints":None, "gamma_only":None,"exclude_bands":None,"spinors":None, "projections":None,
                             "select_projections":None,"auto_projections":None,"spin":None,"translate_home_cell":None,"write_xyz":None, "write_hr":None,
                             "write_rmn":None,"write_tb":None,"hr_cutoff":None, "dist_cutoff":None, "use_ws_distance":None,"fermi_energy":None, 
                             "guiding_centres":None, "ws_distance_tol":None, "ws_search_size":None, "write_u_matrices":None, "dis_win_min":None,
                             "dis_win_max":None, "dis_froz_min":None, "dis_froz_max":None, "dis_num_iter":None, 
                             "dis_mix_ratio":None, "dis_conv_tol":None, "dis_conv_window":None, "dis_spheres_num":None,
                             "dis_spheres_first_wann":None, "dis_spheres":None, "num_iter" :None, "num_cg_steps":None, "conv_window"
                             "conv_tol":None, "use_bloch_phases":None, "site_symmetry":None, "symmetrize_eps":None});


        #User defined variables
        self.set(num_wann = num_wann)
        self.set(mp_grid  = mp_grid); 
        self.set(kpoints  = self.get_kpoints() );

        #Variables read from the previous calculation
        winp = WannierInput(xml=xml);
        self.set(unit_cell_cart = winp.get_cell() );
        self.set(atoms_frac = winp.get_fractional_atomic_positions() );
        self.set(num_bands = winp.get_num_bands() );
        self.set(spinors=winp.get_spin_state()['spin']);
        self.set(fermi_energy=winp.get_fermi_energy())

        #Additional variables
        self.set(write_xyz = True);
        self.set(translate_home_cell=True);
        self.set(write_tb=True);
        self.set( projections = "random\n");

        #PW2Wannier90
        self.pw2wann90_params={
            "seedname":seedname, "outdir":winp.get_outdir(),
            "prefix":winp.get_prefix(),"write_dmn":True,"read_sym":True,
            "write_spn":winp.get_spin_state()['spin'],
            "write_uHu":False, "write_uIu":True,"uIu_formatted":True
            };

    #METHODS
    def use_autoprojections(self,mu=None, sigma=0.1, entanglement="erfc"):
        self.set(auto_projections = True);
        self.set( projections = "");
        self.set(guiding_centres=False);

        if mu is None:
            mu = self.fermi_energy;

        entanglement_option = {"isolated",'erfc',"gaussian"};
        #assert (entanglement in entanglement_option);

        self.pw2wann90_params.update(
            {"scdm_proj":True,"scdm_entanglement":'erfc',
             "scdm_mu":mu, "scdm_sigma":sigma
             });

    def use_bloch_phases(self):
        self.set(use_bloch_phases = True);
        self.set(auto_projections = False);
        self.write_dmn =False;
    
    def get_kpoints(self, qe_format =False):
        kx,ky,kz = [ np.linspace(0, 1, n,endpoint=False) for n in self.mp_grid ];
        kxv, kyv, kzv = [ R.flatten() for R in np.meshgrid(kx,ky,kz, indexing='ij') ];

        if qe_format:
            return np.transpose([kxv,kyv,kzv,0*kzv+1.0]);
        return np.transpose([kxv,kyv,kzv]);  

    def set(self, **kwargs ):
        for k, v in kwargs.items():
            if (k not in self.parameters) or  (v is  None):
                raise ValueError("key:",k," does not exists or ",  v, "is None")
            self.parameters.update({k:v})
        self.__dict__.update(self.parameters);
        
    def set_parameters(self, parameters):
        self.parameters = parameters;

    def update_parameters(self, parameters):
        for key in self.parameters.keys():
            if key in parameters:
                self.parameters.update({key:parameters[key]})

    def get_parameters(self):
        return self.parameters;

    def __str__(self):
        #Update parameters from the initialize parameters in the class 
        self.update_parameters(self.__dict__);
        #Check is there is at leas a  None, which means the namelist is not set
        out = ""
        if not all(v is None for v in self.get_parameters().values()):
            for k, v in self.get_parameters().items():
                if v is not None:
                    out += wann_io.key_value_format(k,v)+"\n"
            out += "\n"
            return out
        #If all parameters are none it means the namelist is not present
        return ""

    def use_projection_scheme(self,projections):
        self.set(guiding_centres=True);




    def write_wannier90(self):
        fname =self.seedname+".win";
        try:
            with open(self.seedname+".win", "w") as f:
                    f.write( self.__str__() );
        except:
            raise FileNotFoundError
        return fname



    def run(self, w90_path="./", logfile=None, postprocessing=False):
        inpfile = self.write_wannier90();
        try:
            with open(inpfile) as f:
                pass
        except:
            raise FileNotFoundError
        if logfile is None:
            logfile = inpfile+".log"

        exec = "wannier90.x ";        
        if postprocessing:
            exec+="-pp ";

        proc= subprocess.run([w90_path+exec+inpfile+" |tee "+logfile ], shell=True);


    def pw2wannier90(self, qe_path="./", logfile=None):
        inpfile = self.seedname+".pw2wann90";
        try:
            with open(inpfile, "w") as f:
                f.write("&inputpp\n")
                for k, v in self.pw2wann90_params.items():
                    f.write("\t"+qe_io.key_format(k)+" = "+qe_io.format(v)+"\n");
                f.write("/\n")
        except:
            raise FileNotFoundError
        if logfile is None:
            logfile = inpfile+".log"
        proc= subprocess.run([qe_path+"pw2wannier90.x -inp "+inpfile+" |tee "+logfile ], shell=True);

class Projections():

    def __init__(self, projections=None, use_spin = False):
        self.projections = None;

    def __str__(self) -> str:
        out = "";
        for p in self.projections:
            pos,l, mr, xaxis, zaxis, zona, radial, spin, quant_dir = p;
            out+="f={},{},{}:l={}:mr={}:".format(*pos,l,mr)+\
                 "zaxis={},{},{}:xaxis={},{},{}".format(*zaxis,*xaxis)+\
                 "radial={}:,zona:{}".format(radial,zona);
            if self.use_spin:
                out+="({})[{},{},{}]".format(spin,*quant_dir);
            out+="\n";
         
    def set_projection(self,pos, l, mr, xaxis=(1,0,0), zaxis=(0,0,1), zona=1.0, radial=1.0, spin=None, quant_dir=None):
        self.projections+= (pos,l, mr, xaxis, zaxis, zona, radial, spin, quant_dir);


