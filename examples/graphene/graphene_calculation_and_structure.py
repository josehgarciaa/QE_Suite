# Package directory
import os
import sys
sys.path.insert(0, os.path.abspath('../../'))
sys.path.insert(0, os.path.abspath('.'))

from math import sqrt

from qe_suite.builder import PWInput, System, SCF



#Define calculation
scf = SCF();

graphene = System();

graphene.set_atomic_species( {"C": (12.0107, "C.pbesol-n-kjpaw_psl.1.0.0.UPF")} )
graphene.set_atomic_positions("crystal",
                              [("C", 0.6667, 0.3333, 0.5),
                               ("C", 0.3333, 0.6667, 0.5)])
graphene.set_cell_parameters("angstrom",
                         [(2.467, 0.000, 0.0),
                          (-1.234, 2.137, 0.0),
                          (0.000, 0.000, 15.0)]);
graphene.use_structure_as_symmetries();
graphene.refine_symmetries();



scf = SCF().set_k_points("automatic", [15, 15, 3, 0, 0, 0])

pw_input = PWInput(calculation = scf, system=graphene  );
pw_input.write("calculation.out")