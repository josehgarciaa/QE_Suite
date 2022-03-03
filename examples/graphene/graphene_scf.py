
import math as m

import qe_suite.builder as builder;


#import qe_suite.qe_suite as qes
#import qe_suite.io as qe_io


system = "graphene"


#Define Graphene structure
a= 2.467  # lattice constant in nm which is the default unit. 
cell = [ [a,0,0], [-a/2, a*m.sqrt(3)/2,0], [0,0,5*a] ];
fractional_positions = [ [0,0,0.0], [2/3,1/3,0.0] ]
atomic_symbols = ['C', 'C']
structure = builder.Structure(cell, fractional_positions, atomic_symbols);
#Evaluate the symmetries of the structure and symmetrize it
print("The space group is: ", structure.hm_symbol() )

#Create a configuration file based on graphene's structure
qe_input = builder.QEInput(structure= structure);

#xyz  = qe_io.load_xyz(syst+".xyz");
#cell = qe_io.load_cell(syst+".uc");

#qes_handler =  qes.generate_from_xyz(xyz=xyz, cell=cell, two_dimensional=True);

#qes_handler.set_calculation("scf");
#qes_handler.use_symmetries();
#qes_handler.use_SSSP(functional="PBEsol", target="precision", path="../SSSP");

#qes_handler.write_input_file("QEsuite.scf");



