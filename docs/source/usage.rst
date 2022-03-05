Documentation
----------------

The core  of |QESuite| is the Class :py:class:`~qe_suite.builder.PWInput` which serves as a container of an inputfile. 

PWInput is structured in namelists and cards following the coventions in `pw.x manual <https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm3>`_, 
changing certain inputs into python objects as explained in DICTIONARY DOCUMENTATION. Therefore, there is a unique mapping
from the :py:class:`~qe_suite.builder.PWInput` object to an inputfile. An example of this is given below

.. doctest::
  
  >>> import qe_suite.builder.PWInput as PWInput

  >>> pw_input = PWInput();
   
  >>> control = pw_input.namelists.control;
  >>> control.calculation = 'scf'
  >>> control.etot_conv_thr =   2.0e-05
  >>> control.forc_conv_thr =   1.e-04
  >>> control.outdir = './out/'
  >>> control.prefix = 'aiida'
  >>> control.pseudo_dir = './pseudo/'
  >>> control.tprnfor = True
  >>> control.tstress = True
  >>> control.verbosity = 'high'
  >>> 
  >>> system = pw_input.namelists.system;
  >>> system.degauss =   1.46997e-02
  >>> system.ecutrho =   3.60000e+02
  >>> system.ecutwfc =   4.50000e+01
  >>> system.ibrav = 0
  >>> system.nat = 2
  >>> system.nosym = False
  >>> system.ntyp = 1
  >>> system.occupations = 'smearing'
  >>> system.smearing = 'cold'

  >>> electrons = pw_input.namelists.electrons;
  >>> electrons.conv_thr =   4.0e-10
  >>> electrons.electron_maxstep = 80
  >>> electrons.mixing_beta =   4.0e-01

  >>> cards = pw_input.cards;
  >>> cards.atomic_species = {"C": (12.0107,C.pbesol-n-kjpaw_psl.1.0.0.UPF) };
  >>> cards.atomic_postions.option = "crystal"
  >>> cards.atomic_postions =  [ ("C",0.6667,0.3333,0.5), 
  >>>                            ("C",0.3333,0.6667,0.5) ] 
  >>> cards.k_points.option = "automatic"
  >>> cards.k_points =  [15,15,3,0,0,0]
  >>> cards.cell_parameters.option = "angstrom"
  >>> cards.cell_parameters.option = [ ( 2.467, 0.000, 0.0), 
  >>>                                  (-1.234, 2.137, 0.0), 
  >>>                                  ( 0.000, 0.000,15.0)];

  >>> pw_input.print();
  &CONTROL
    calculation = 'scf'
    etot_conv_thr =   2.0d-05
    forc_conv_thr =   1.0d-04
    outdir = './out/'
    prefix = 'aiida'
    pseudo_dir = './pseudo/'
    tprnfor = .true.
    tstress = .true.
    verbosity = 'high'
  /
  &SYSTEM
    degauss =   1.46997d-02
    ecutrho =   3.60000d+02
    ecutwfc =   4.50000d+01
    ibrav = 0
    nat = 2
    nosym = .false.
    ntyp = 1
    occupations = 'smearing'
    smearing = 'cold'
  /
  &ELECTRONS
    conv_thr =   4.0d-10
    electron_maxstep = 80
    mixing_beta =   4.0d-01
  /
  ATOMIC_SPECIES
  C      12.0107 C.pbesol-n-kjpaw_psl.1.0.0.UPF
  ATOMIC_POSITIONS crystal
  C            0.6667       0.3333       0.5 
  C            0.3333       0.6667       0.5 
  K_POINTS automatic
  15 15 3 0 0 0
  CELL_PARAMETERS angstrom
        2.467       0.000       0.0
      -1.234       2.1367       0.0
        0.000       0.000      15.0

We must higlight however that the previous form is not recommended since no check is performed on the inputs and there are other
:py:class:`~qe_suite.builder.PWInput` methods that will provided better readibility


Modules
=======

.. automodule:: qe_suite.builder
  :members:  
  :undoc-members:  

Classes
_______

.. autoclass:: qe_suite.builder.PWInput
  :members:  
  :undoc-members:  

.. autoclass:: qe_suite.builder.Structure
  :members:  
  :undoc-members:  



