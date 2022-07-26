import qe_suite.io.wannier90 as wann_io
from qe_suite.parse.xml import WannierInput

class Projections():

    def __init__(self, projections=None;, use_spin = False):
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

    def __init__(self, num_wann, name="qe_suite.wannier", xml=None):
        self.name = name
        self.set_parameters({"num_wann": None, "num_bands": None,"unit_cell_cart": None,"atoms_frac":None,"mp_grid":None,
                             "kpoints":None, "gamma_only":None,"exclude_bands":None,"spinors":None,
                             "set_projections":None,"spin":None,"translate_home_cell":None,"write_xyz":None, "write_hr":None,
                             "write_rmn":None,"write_tb":None,"hr_cutoff":None, "dist_cutoff":None, "use_ws_distance":None, 
                             "ws_distance_tol":None, "ws_search_size":None, "write_u_matrices":None, "dis_win_min":None,
                             "dis_win_max":None, "dis_froz_min":None, "dis_froz_max":None, "dis_num_iter":None, 
                             "dis_mix_ratio":None, "dis_conv_tol":None, "dis_conv_window":None, "dis_spheres_num":None,
                             "dis_spheres_first_wann":None, "dis_spheres":None, "num_iter" :None, "num_cg_steps":None, "conv_window"
                             "conv_tol":None, "use_bloch_phases":None, "site_symmetry":None, "symmetrize_eps":None});

        self.write_tb  = True;
        self.write_xyz = True;

        #Variables read from the pwinput
        #self.set(mp_grid=);

        winp = WannierInput(xml=xml);
        self.set(num_wann = num_wann)
        self.set(unit_cell_cart = winp.get_cell() );
        self.set(atoms_frac = winp.get_fractional_atomic_positions() );
        self.set(kpoints =  winp.get_kpoints() );
        self.set(num_bands = winp.get_num_bands() );
        self.set(spinors=winp.get_spin_state()['spin']);
        self.set(write_xyz = True);
        self.set(translate_home_cell=True);
        self.set(auto_projections = True);
        self.set(write_tb=True);

    def use_bloch_phases(self):
        self.set(use_bloch_phases = True);
        self.set(auto_projections = False);


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
            out += "/\n"
            return out
        #If all parameters are none it means the namelist is not present
        return ""

    def use_projection_scheme(self,projections):
        self.set(guiding_centres=True);
