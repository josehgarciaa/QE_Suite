import qe_suite.io as qe_io 

class Electrons:

    electron_maxstep= None; scf_must_converge= None;
    conv_thr= None; adaptive_thr=None; conv_thr_init=None;
    conv_thr_multi= None; mixing_mode=None;
    mixing_beta= None; mixing_ndim= None;
    mixing_fixed_ns= None; diagonalization= None;
    diago_thr_init= None; diago_cg_maxiter= None;
    diago_david_ndim= None; diago_full_acc= None;
    diago_rmm_ndim= None; diago_rmm_conv= None;
    efield= None; efield_cart= None; efield_phase= None;
    startingpot= None; startingwfc= None; tqr= None;
    real_space= None;
    
    """The electronic parameters used in the simulation. 

        Attributes:
        ------------
        This class implement attributes that follow the naming convenction and with
        the same description as given in the namelist Quantum Espresso in `\&ELECTRONS <https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm820>`_
        Unless otherwise specifed, the attributes are to `None` and won't be included in the input_file which 
        will force Quantum Espresso to use its default values.

        Examples of attribute initialization:
        ------------
        For instance, we could define the electronic convergencescin the &ELECTRONS namelist in different ways

        >>> # qe_input.electrons.conv_thr = 400"

        Parameters:
        ------------
        None:
            

        Example
        ------------
        >>> # qe_input.control.etot_conv_thr	= 1e04;
    """    

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
