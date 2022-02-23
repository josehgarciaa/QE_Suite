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

icells = [s.get_cell() for s in structures];
max_cell = np.diag((5,5,1));
max_area = np.linalg.det( max_cell.dot(icells[0])[:2,:2]);
cells,params,min_diff=  vdw.get_vdw_cell( *structures, max_strain=0.1, max_area=max_area ) ; 
print(cells[0].area(2),cells[1].area(2))
print(cells[0],cells[1])
print("diff",min_diff)
print("params",params)

transfs =[ vdw.get_cell_transformation(icells[0],cells[0]), np.round( vdw.get_cell_transformation(icells[1],cells[1]) ) ];
sc_structures =[None,None];
sc_structures[0]=ase.build.make_supercell(structures[0],transfs[0], wrap=True, tol=1e-05)
sc_structures[1]=ase.build.make_supercell(structures[1],transfs[1], wrap=True, tol=1e-05)
rcell = sc_structures[0].get_cell();
scell = sc_structures[1].get_cell()
strain= vdw.get_cell_transformation(scell,rcell)
sc_structures[1].set_cell(rcell, scale_atoms=True)
vdw_structure = ase.build.stack(*sc_structures, axis=2, maxstrain=1e-3, distance = 2.0);

##fix height
cell   = vdw_structure.get_cell();
cell[2]= icells[0][2]; 
vdw_structure.set_cell(cell);


#get xyz info
positions= vdw_structure.get_positions();
symbols  = vdw_structure.get_chemical_symbols();
xyz = [ (x,list(p)) for x,p in zip(symbols,positions)];

qeio.write_xyz(xyz,"vdw.xyz")
qes_handler = generate_from_xyz(xyz=xyz, cell=cell, two_dimensional=True) 
#qes_handler.use_symmetries();
qes_handler.set_calculation("scf");
qes_handler.use_SSSP(functional="PBEsol", target="precision", path="../SSSP");
qes_handler.write_input_file("QEsuite.scf");



