from .namelist import Namelist


class Ions(Namelist):

    """The ions parameters used in the simulation. 

        Attributes:
        ------------
        This class implement attributes that follow the naming convenction and with
        the same description as given in the namelist Quantum Espresso in `\&IONS <https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm932>`_

        Unless otherwise specifed, the attributes are to `None` and won't be included in the input_file which 
        will force Quantum Espresso to use its default values.

        Examples of attribute initialization:
        ------------
        For instance, we could define the wave function cutoff in the &IONS namelist in different ways

        >>> # qe_input.ions.? = "

        Parameters:
        ------------
        structure:
            An instance of a :py:class:`~qe_suite.builder.Structure`
        calculation:
            An instance of a :py:class:`~qe_suite.builder.Calculation`

        Example
        ------------
        >>> # qe_input.control.etot_conv_thr	= 1e04;
    """

    options = dict()

    def __init__(self):
        self.__dict__.update(
            {"ion_positions": None, "ion_velocities": None, "ion_dynamics": None,
                "pot_extrapolation": None, "wfc_extrapolation": None, "remove_rigid_rot": None,
                "ion_temperature": None, "tempw": None, "tolp": None, "delta_t": None,
                "nraise": None, "refold_pos": None, "upscale": None, "bfgs_ndim": None,
                "trust_radius_max": None, "trust_radius_min": None, "trust_radius_ini": None,
                "w_1": None, "w_2": None, "fire_alpha_init": None, "fire_falpha": None,
                "fire_nmin": None, "fire_f_inc": None, "fire_f_dec": None, "fire_dtmax": None, "dtmax": None}
        )
        self.set_namelist_name("&IONS");
