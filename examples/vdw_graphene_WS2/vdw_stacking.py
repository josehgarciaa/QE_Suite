import qe_suite.qe_suite as qes
import qe_suite.vdw as vdw
import qe_suite.io as qeio
from qe_suite.io import  load_xyz,load_cell 
from qe_suite.qe_suite import  generate_from_xyz
from ase import Atoms
import numpy as np

systs = ( "C2", "WS2" );


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
aicell,bicell =  [s.get_cell() for s in structures ];


def get_vdw_cell( a_structure, b_structure, max_strain=0.2, strain_cell="b", max_size=10, max_area=100 ):

    if strain_cell=="a":
        #invert since the algorithm will always strain b
        a_structure, b_structure= b_structure, a_structure

    from scipy.optimize import differential_evolution
    a_cell,b_cell = a_structure.get_cell(),b_structure.get_cell();
    def diff(x):
        ds= x[0]
        S = np.diag([1+ds,1+ds,1]);
        closest_cells = vdw.get_closest_cells( a_cell, S.dot(b_cell), max_size=max_size);
        if closest_cells is None:
            return np.inf;
        ab_scells, diff= closest_cells
        return diff
    bounds = [(-max_strain,max_strain)]
    res = differential_evolution(diff, bounds,  polish=True );
    opt_strain = 1+res.x[0];
    S = np.diag([opt_strain,opt_strain,1]);

    closest_cells = vdw.get_closest_cells( a_cell, S.dot(b_cell), max_size=max_size);
    opt_ab_scells, opt_diff= closest_cells
    opt_a,opt_b = list(map(vdw.convert_to_cell,opt_ab_scells));

    if strain_cell == "a":
        #invert the resulting cell since the algorithm assumed strained b
        opt_a,opt_b = opt_b,opt_a;

    return opt_ab_scells, opt_diff,opt_strain


opt_ab_scells, opt_diff,opt_strain = get_vdw_cell( *structures, max_strain=0.3, strain_cell="a", max_size=30 );

if( opt_diff!=0):
    print("The resulting cell is not optimal")


#icell =vdw.convert_to_cell(min_scell[0])
#a,b,c,bc,ac,ab= icell.cellpar();
#icell =vdw.convert_to_cell(np.round(icell/a,4));
#print("par",a, icell.cellpar() )


#icell =vdw.convert_to_cell(min_scell[1])
#a,b,c,bc,ac,ab= icell.cellpar();
#icell =vdw.convert_to_cell(np.round(icell/a,4));
#print("par",a, icell.cellpar() )

#sc_structures= [ vdw.expand_supercell(struct,sc) for struct,sc in zip(structures,min_scell) ];

#Construct an stacked structure
#vdw_cell =sc_structures[0].get_cell();
#print("vdw_cell",vdw_cell)
#vdw_structure = ase.build.stack(*sc_structures, axis=2, maxstrain=0.1, distance=2);

#fix height
#cell   = vdw_structure.get_cell();
#cell[2]= vdw_cell[2]; 
#vdw_structure.set_cell(cell);

#get xyz info
#positions= vdw_structure.get_positions();
#symbols  = vdw_structure.get_chemical_symbols();
#xyz = [ (x,list(p)) for x,p in zip(symbols,positions)];

#qeio.write_xyz(xyz,"vdw.xyz")
#qes_handler = generate_from_xyz(xyz=xyz, cell=cell, two_dimensional=True) 

#qes_handler.set_calculation("scf");
#qes_handler.use_SSSP(functional="PBEsol", target="precision", path="../SSSP");
#qes_handler.write_input_file("QEsuite.scf");


