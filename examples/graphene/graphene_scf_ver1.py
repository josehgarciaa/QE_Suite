import os
import sys
#Package directory
sys.path.insert(0, os.path.abspath('../../'))
sys.path.insert(0, os.path.abspath('.'))

from qe_suite.builder import PWInput
pw_input = PWInput()

control = pw_input.namelists.control
control.calculation = 'scf'
control.etot_conv_thr = 2.0e-05
control.forc_conv_thr = 1.e-04
control.outdir = './out/'
control.prefix = 'aiida'
control.pseudo_dir = './pseudo/'
control.tprnfor = True
control.tstress = True
control.verbosity = 'high'

system = pw_input.namelists.system
system.degauss = 1.46997e-02
system.ecutrho = 3.60000e+02
system.ecutwfc = 4.50000e+01
system.ibrav = 0
system.nat = 2
system.nosym = False
system.ntyp = 1
system.occupations = 'smearing'
system.smearing = 'cold'

electrons = pw_input.namelists.electrons
electrons.conv_thr = 4.0e-10
electrons.electron_maxstep = 80
electrons.mixing_beta = 4.0e-01

cards = pw_input.cards
cards.atomic_species.value = {"C": (12.0107, "C.pbesol-n-kjpaw_psl.1.0.0.UPF")}
cards.atomic_positions.option = "crystal"
cards.atomic_positions.value = [("C", 0.6667, 0.3333, 0.5),
                         ("C", 0.3333, 0.6667, 0.5)]
cards.k_points.option = "automatic"
cards.k_points.value = [15, 15, 3, 0, 0, 0]
cards.cell_parameters.option = "angstrom"
cards.cell_parameters.value = [(2.467, 0.000, 0.0),
                                (-1.234, 2.137, 0.0),
                                (0.000, 0.000, 15.0)]


print(str(pw_input))

pw_input.write("out")