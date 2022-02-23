import qe_suite.qe_suite as qes
import qe_suite.vdw as vdw
import qe_suite.io as qeio
from qe_suite.io import  load_xyz,load_cell 
from qe_suite.qe_suite import  generate_from_xyz
import ase
from ase import Atoms
import numpy as np

systs = ( "C2", "Pd2Se4" );


#load the geometry of the files in two structures
structures = [];
for s in systs:
    xyz =load_xyz(s+".xyz");
    cell=load_cell(s+".uc");
    positions=  [(x,y,z) for s,(x,y,z) in xyz];
    symbols  =  [s for s,(x,y,z) in xyz];
    structure=  Atoms( positions=(positions), symbols=symbols, cell=cell);
    structures.append(structure);

#Find the most similar supercell structures within a strain tolerance
opt_scells, opt_diff,opt_strain = vdw.get_vdw_cell( *structures, max_strain=0.3, max_size=30 );

if( opt_diff>1e-1):
    print("The resulting cell is not optimal")
else:
    print(opt_scells, opt_diff,opt_strain)


sc_structures= [ vdw.expand_supercell(struct,sc) for struct,sc in zip(structures,opt_scells) ];

#Construct an stacked structure
vdw_cell =sc_structures[0].get_cell();
print("vdw_cell",vdw_cell)
vdw_structure = ase.build.stack(*sc_structures, axis=2, maxstrain=0.1, distance=2);

#fix height
cell   = vdw_structure.get_cell();
cell[2]= vdw_cell[2]; 
vdw_structure.set_cell(cell);

#get xyz info
positions= vdw_structure.get_positions();
symbols  = vdw_structure.get_chemical_symbols();
xyz = [ (x,list(p)) for x,p in zip(symbols,positions)];

qeio.write_xyz(xyz,"vdw.xyz")
qes_handler = generate_from_xyz(xyz=xyz, cell=cell, two_dimensional=True) 
qes_handler.set_calculation("scf");
qes_handler.use_SSSP(functional="PBEsol", target="precision", path="../SSSP");
qes_handler.write_input_file("QEsuite.scf");



