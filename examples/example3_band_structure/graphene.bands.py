# Package directory
import os
from re import M
import sys
sys.path.insert(0, os.path.abspath('../../'))
sys.path.insert(0, os.path.abspath('.'))

# Package directory
from math import sqrt
import qe_suite as qes
from qe_suite.builder import Structure, System, PWInput
from qe_suite.builder.calculation import SCF, Bands



#First we define the structure of the system
hexagonal= Structure( cell = [(2.467, 0.000, 0.0), (-1.234, 2.137, 0.0), (0.000, 0.000, 15.0)],
                      fractional_positions = [ (0.6667, 0.3333, 0.5), (0.3333, 0.6667, 0.5) ],
                      atomic_symbols = ["C", "C"]
                    );
hexagonal.set_as_2D();

#That structure is used to define the system to be study
graphene = System();
scf = SCF();
graphene.set_atomic_species( {"C": (12.0107, "C.pbesol-n-kjpaw_psl.1.0.0.UPF")} )
graphene.set_structure(hexagonal)

#Then we create a self consistent calculation to determine the density matrix
scf = SCF().set_pseudopot_dir("../SSSP")
scf.set_k_points("automatic", hexagonal.get_kpoints(type="automatic"))
pw_input = PWInput(calculation = scf, system=graphene  );
qes.run_pw(inpfile=pw_input.write("qe_suite.scf.inp") );

#Finally we analyze the system by performing a band structure calculation
bands = Bands(scf=scf).set_band_path( hexagonal.symmetrized_band_path() );
pw_input = PWInput(calculation = bands, system=graphene.set_numbands(5) );
qes.run_pw(inpfile=pw_input.write("qe_suite.bands.inp") );


