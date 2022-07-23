import qe_suite.io.wannier90 as wann_io

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

    def __init__(self, num_wann, pwinput, name="qe_suite.wannier"):
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

        self.num_wann=num_wann;
        self.write_tb  = True;
        self.write_xyz = True;

        #Variables read from the pwinput
        #self.set(num_bands=);
        #self.set(unit_cell_cart=);
        #self.set(atoms_frac=);
        #self.set(mp_grid=);
        #self.set(kpoints=);
        #self.set(spinors=);


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
        if not all(v is None for v in self.get_parameters().values()):
            out = self.name+"\n";
            for k, v in self.get_parameters().items():
                if v is not None and k != "name":
                    out += qe_io.key_format(k)+"="+qe_io.format(v)+"\n"
            out += "/\n"
            return out
        #If all parameters are none it means the namelist is not present
        return ""