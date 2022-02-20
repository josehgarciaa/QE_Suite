import qe_suite.qe_suite as qes
import qe_suite.io as qe_io
import qe_suite.vdw as vdw



a_syst = "C2"
xyz  = qe_io.load_xyz(a_syst+".xyz");
cell = qe_io.load_cell(a_syst+".uc");
a_qes_h =  qes.generate_from_xyz(xyz=xyz, cell=cell, two_dimensional=True);

b_syst = "Pd2Se4"
xyz  = qe_io.load_xyz(b_syst+".xyz");
cell = qe_io.load_cell(b_syst+".uc");
b_qes_h =  qes.generate_from_xyz(xyz=xyz, cell=cell, two_dimensional=True);


a_structure = a_qes_h.get_structure();
b_structure = b_qes_h.get_structure();

min_scatms, min_diff,min_ds = vdw.get_vdw_cell( a_structure, b_structure, max_strain=0.1, strain_cell="a", max_size=30 );

import numpy as np
from ase import Atoms


sccell_a,sccell_b =min_scatms
structure = a_structure;
sc_cell= sccell_a; 


