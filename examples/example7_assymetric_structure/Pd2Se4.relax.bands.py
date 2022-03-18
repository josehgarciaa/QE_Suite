# Package directory
import os
from re import M
import sys
sys.path.insert(0, os.path.abspath('../../'))
sys.path.insert(0, os.path.abspath('.'))

from math import sqrt
import qe_suite as qes
from qe_suite.builder import Structure, System, PWInput
from qe_suite.builder.calculation import Relaxation, Bands

#First we define the structure of the system
monoclinic = Structure().load_from(xyz_file="Pd2Se4.xyz", cell_file="Pd2Se4.uc",symmetrize=True);
monoclinic.set_as_2D();

#That structure is used to define the system to be study
Pd2Se4 = System();
Pd2Se4.set_atomic_species( ["Pd","Se"], use_SSSP= "../SSSP/SSSP_1.1.2_PBEsol_precision.json" );
Pd2Se4.set_structure(monoclinic, use_symmetries=True);

#Then we create a self consistent calculation to determine the density matrix
relax = Relaxation().set_pseudopot_dir("../SSSP").set_cell_do_free('ibrav');
relax.set_k_points("automatic", monoclinic.get_kpoints(type="automatic"))

pw_input = PWInput(calculation = relax, system=Pd2Se4  );
inpfile=pw_input.write("qe_suite.relax.inp")
qes.run_pw(inpfile=pw_input.write("qe_suite.relax.inp"),shell=False );

#Finally we analyze the system by performing a band structure calculation
bands = Bands(scf=relax).set_band_path( monoclinic.symmetrized_band_path() );
pw_input = PWInput(calculation = bands, system=Pd2Se4 );
inpfile=pw_input.write("qe_suite.bands.inp")
qes.run_pw(inpfile=pw_input.write("qe_suite.bands.inp") );


import qe_suite.parse as parse
import matplotlib.pyplot as plt

bs = parse.BandStructure(xml="qe_suite/qe_suite.xml");
hsp =monoclinic.get_symmetrized_high_symmetry_points();

for band in bs.bands(shift_Efermi=True):
    plt.plot(*band);
plt.xticks( *bs.xticks_high_symmetry_points( hsp= hsp ) );
plt.ylabel( bs.get_ylabel() );
plt.savefig('band_structure.pdf');  
