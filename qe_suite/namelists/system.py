from .namelist import Namelist


class System(Namelist):
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
        super().__init__()

        self.set_parameters({"ibrav": None, "celldm_1": None,"celldm_2": None, "celldm_3": None,
                              "celldm_4": None,"celldm_5": None,"celldm_6": None,
                              "A": None, "B": None, "C": None, "cosAB": None, "cosAC": None, "cosBC": None,
                              "nat": None, "ntyp": None, "nbnd": None, "ecutwfc": None, "ecutrho": None, "ecutfock": None,
                              "tot_charge": None, "starting_charge": None, "tot_magnetization": None, "starting_magnetization": None,
                              "nr1": None, "nr2": None, "nr3": None, "nr1s": None, "nr2s": None, "nr3s": None,
                              "nosym": None, "nosym_evc": None, "noinv": None, "no_t_rev": None, "force_symmorphic": None, "use_all_frac": None,
                              "occupations": None, "one_atom_occupations": None, "starting_spin_angle": None, "degauss": None, "smearing": None,
                              "nspins": None, "noncolin": None, "ecfixed": None, "qcutz": None, "q2sigma": None,
                              "input_dft": None, "ace": None, "exx_fraction": None, "screening_parameter": None, "exxdiv_treatment": None,
                              "x_gamma_extrapolation": None, "ecutvcut": None, "nqx1": None, "nqx2": None, "nqx3": None, "localization_thr": None,
                              "lda_plus_u": None, "lda_plus_u_kind": None, "Hubbard_U": None, "Hubbard_J0": None, "Hubbard_V": None,
                              "Hubbard_alpha": None, "Hubbard_beta": None, "Hubbard_J": None, "starting_ns_eigenvalue": None,
                              "U_projection_type": None, "Hubbard_parameters": None,
                              "dmft": None, "dmft_prefix": None, "ensemble_energies": None,
                              "edir": None, "emaxpos": None, "eopreg": None, "eamp": None, "zgate": None,
                              "relaxz": None, "block": None, "block_1": None, "block_2": None, "block_height": None,
                              "angle1": None, "angle2": None, "lforcet": None, "constrained_magnetization": None,
                              "fixed_magnetization": None, "report": None, "Lambda": None, "lspinorb": None,
                              "assume_isolated": None, "esm_bc": None, "esm_w": None, "esm_efield": None, "esm_nfit": None,
                              "lgcscf": None, "gcscf_mu": None, "gcscf_conv_thr": None, "gcscf_beta": None, "vdw_corr": None,
                              "london_s6": None, "london_c6": None, "london_rvdw": None, "london_rcut": None, "dftd3_version": None,
                              "dftd3_threebody": None, "ts_vdw_econv_thr": None, "ts_vdw_isolated": None, "xdm_a1": None, "xdm_a2": None,
                              "space_group": None, "uniqueb": None, "origin_choice": None, "rhombohedral": None});
        self.set_name("&SYSTEM");
        self.set(ibrav=0)        

    def set_bravais_lattice(self, ibrav, celldm):
        self.options["ibrav"] = ibrav
        self.celldm = celldm

    def set_cutoff(self, ecutwfc, ecutrho=None):
        if ecutrho is None:
            ecutrho = 4*ecutwfc
        self.options["ecutwfc"] = ecutwfc
        self.options["ecutrho"] = ecutrho

