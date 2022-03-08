# Package directory
import os
import sys
sys.path.insert(0, os.path.abspath('../../'))
sys.path.insert(0, os.path.abspath('.'))



from qe_suite.builder import PWInput
from qe_suite.builder import SCF


scf = SCF()
scf.set_atomic_species( {"C": (12.0107, "C.pbesol-n-kjpaw_psl.1.0.0.UPF")} )
scf.set_k_points("automatic", [15, 15, 3, 0, 0, 0])
scf.set_atomic_positions("crystal",
                              [("C", 0.6667, 0.3333, 0.5),
                               ("C", 0.3333, 0.6667, 0.5)])
scf.set_cell_parameters("angstrom",
                         [(2.467, 0.000, 0.0),
                          (-1.234, 2.137, 0.0),
                          (0.000, 0.000, 15.0)])

pw_input = PWInput(calculation = scf );



calc = SCF()
