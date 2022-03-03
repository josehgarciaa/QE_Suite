import qe_suite.io as qe_io 
class System:

    #Typical  required structural parameters
    ibrav= None; celldm= None;
    A, B, C, cosAB, cosAC, cosBC = [None]*6;
    nat= None; ntyp=None; nbnd=None;

    #Typical  required simulation parameters parameters
    ecutwfc=None; ecutrho= None; 
    ecutfock = None;	

    #Initial electronic and magnetic configuration
    tot_charge= None; starting_charge= None;
    tot_magnetization= None; starting_magnetization= None;
    
    # Symmetries and FFT  
    nr1, nr2, nr3 = [None]*3;
    nr1s, nr2s, nr3s = [None]*3;
    nosym =None; nosym_evc = None;
    noinv =None; no_t_rev =None;
    force_symmorphic = None; use_all_frac = None;

    #Occupation parameters
    occupations =None; one_atom_occupations = None;
    starting_spin_angle = None;	
    degauss = None; smearing = None;
    nspins = None; noncolin=None;
    
    #hybrid functional definitions
    ecfixed = None; qcutz= None; q2sigma=None;
    input_dft = None; ace = None; exx_fraction =None;
    screening_parameter = None; exxdiv_treatment =None
    x_gamma_extrapolation = None; ecutvcut =None;
    nqx1, nqx2, nqx3 = [None]*3;
    localization_thr = None;

    #DFT+U parameters
    lda_plus_u = None; lda_plus_u_kind = None;
    Hubbard_U = [None,]; Hubbard_J0= [None,]; Hubbard_V = [ [None,],];
    Hubbard_alpha = [None,]; Hubbard_beta= [None,]; Hubbard_J = [ [None,],];
    starting_ns_eigenvalue = [None,]; U_projection_type = None;
    Hubbard_parameters = None; 

    #DMFT parameters
    dmft = None; dmft_prefix= None;
    ensemble_energies = None; 
    
    #Electric field parameter
    edir =None;
    emaxpos = None; eopreg=None; eamp = None;
    zgate = None; relaxz = None; block=None;
    block_1=None; block_2 = None; block_height=None;

    #Non-collinear calculation parameters
    angle1= None; angle2 = None; lforcet = None;
    constrained_magnetization = None;
    fixed_magnetization = None;
    report = None;

    #Spin-orbit coupling parameters
    Lambda = None;
    lspinorb = None;	

    # Isolated calculation parameters
    assume_isolated = None; esm_bc = None;	
    esm_w = None; esm_efield=None;
    esm_nfit = None;
    
    # Grand Canonical Self consistent calculations parameters
    lgcscf	= None;
    gcscf_mu = None; gcscf_conv_thr = None; gcscf_beta = None;	
    
    # Van der Waal parameters
    vdw_corr = None;
    london_s6= None; london_c6 = None; london_rvdw = None;
    london_rcut = None;
    dftd3_version = None; dftd3_threebody = None; ts_vdw_econv_thr = None;
    ts_vdw_isolated = None;
    xdm_a1 = None; xdm_a2 = None;
    
    # Symmetry considetation parameters
    space_group = None;
    uniqueb = None;         
    origin_choice = None;     
    rhombohedral = None;
    
    """The system parameters used in the simulation. 

        Attributes:
        ------------
        This class implement attributes that follow the naming convenction and with
        the same description as given in the namelist Quantum Espresso in `\&SYSTEM <https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm199>`_
        with the exception of lambda which was changed to Lambda
        Unless otherwise specifed, the attributes are to `None` and won't be included in the input_file which 
        will force Quantum Espresso to use its default values.

        Examples of attribute initialization:
        ------------
        For instance, we could define the wave function cutoff in the &SYSTEM namelist in different ways

        >>> # qe_input.system.ecutwfc = 400"

        Parameters:
        ------------
        structure:
            An instance of a :py:class:`~qe_suite.builder.Structure`

        Example
        ------------
        >>> # qe_input.control.etot_conv_thr	= 1e04;
    """    
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
                if c is not None:
                    out+= "celldm("+qe_io.format(1+i)+")="+qe_io.format(c)+"\n";

        out += "/";
        return out;

    def print(self):
        print(self.text());
