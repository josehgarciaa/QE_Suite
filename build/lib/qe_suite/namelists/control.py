#from ..builder.calculation import Target
from .namelist import Namelist
#default_control = { "outdir": "./qesuite_out", "pseudo_dir":"./", "tprnfor":True,"tstress":True, "verbosity":'high'};
#performance = default_control.update( {"etot_conv_thr":1e-4, "forc_conv_thr":1e-3 });
#efficiency  = default_control.update( {"etot_conv_thr":1e-4, "forc_conv_thr":1e-3 });

# option = { Target.default: default_control,
#           Target.efficiency:performance,
#           Target.performance:performance }


class Control(Namelist):

    """The control parameters that defines a simulation. 

      Attributes:
      ------------
        This class implement attributes that follow the naming convenction and with
        the same description as given in the namelist Quantum Espresso in `\&CONTROL <https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm32>`_
        Unless otherwise specifed, the attributes are to `None` and won't be included in the input_file which 
        will force Quantum Espresso to use its default values.

      Examples of attribute initialization:
      ------------
      For instance, we could define the calculation in the &CONTROL namelist in different ways

      >>> # qe_input.control.calculation = "scf";

      Similarly, the convergence threshold on total energy could be defined as

      >>> # qe_input.control.etot_conv_thr	= 1e04;

      Parameters:
      ------------
      calculation:
          An instance of a :py:class:`~qe_suite.builder.Calculation`

      Example
      ------------
      >>> # qe_input.control.etot_conv_thr	= 1e04;
    """

    def __init__(self, calculation=None):
        self.__dict__.update(
            {"calculation": None, "title": None, "verbosity": None,
             "restart_mode": None, "nstep": None, "iprint": None,
             "tstress": None, "tprnfor": None, "dt": None, "outdir": None,
             "wfcdir": None, "prefix": None, "max_seconds": None,
             "etot_conv_thr": None, "forc_conv_thr": None, "disk_io": None,
             "pseudo_dir": None, "tefield": None, "dipfield": None, "lelfield": None,
             "nberrycyc": None, "lorbm": None, "lberry": None, "gdir": None,
             "nppstr": None, "lfcp": None, "gate": None})
        self.set_namelist_name("&CONTROL");
