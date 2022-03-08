# Package directory
import os
import sys
sys.path.insert(0, os.path.abspath('../../'))
sys.path.insert(0, os.path.abspath('.'))

from math import sqrt

from qe_suite.builder import Structure, System, PWInput
from qe_suite.builder.calculation import SCF, Bands



#Define the system
hexagonal= Structure( cell = [(2.467, 0.000, 0.0), (-1.234, 2.137, 0.0), (0.000, 0.000, 15.0)],
                      fractional_positions = [ (0.6667, 0.3333, 0.5), (0.3333, 0.6667, 0.5) ],
                      atomic_symbols = ["C", "C"]
                    );
graphene = System();
scf = SCF();
graphene.set_atomic_species( {"C": (12.0107, "C.pbesol-n-kjpaw_psl.1.0.0.UPF")} )
graphene.set_structure(hexagonal)

#Create a band structure calculation
scf = SCF().set_k_points("automatic", [15, 15, 3, 0, 0, 0]).set_pseudopot_dir("../SSSP")
pw_input = PWInput(calculation = scf, system=graphene  );
output_file = "qe_suite.scf.out";
pw_input.write(output_file)
import subprocess
proc= subprocess.run(["pw.x -inp "+output_file], shell=True);

#Creates a band structure calculation
bands = Bands(scf=scf).set_band_path( hexagonal.symmetrized_band_path() );
pw_input = PWInput(calculation = bands, system=graphene  );
output_file = "qe_suite.bands.out";
pw_input.write(output_file)
import subprocess
proc= subprocess.run(["pw.x -inp "+output_file], shell=True);

