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

def atoms_in_cell(positions,cell):
    try:
        pos_iter = iter(positions);
        scal_pos = cell.scaled_positions(positions);
        allowed_pos =  np.all( (scal_pos < [1,1,1])*(scal_pos >= [0,0,0]),axis=1 )
        return scal_pos[allowed_pos];
    except TypeError as te:
        print( (positions), 'is not iterable')


r0 =0;
L  =1;
cell = a_structure.get_cell() ;
spos = a_structure.get_scaled_positions() ; 
sccell_a,sccell_b = min_scatms;
scal     = np.linalg.norm(sccell_a, axis=0)/np.linalg.norm(sccell_b, axis=0);
sccell_a*= scal;

nmin,nmax=-10,10;
for n0 in range(nmin,nmax):
    for n1 in range(nmin,nmax):
        positions = cell.cartesian_positions ( spos + [n0,n1,0] );
        cell_pos = atoms_in_cell(positions,sccell_a);
        if len(cell_pos)!= 0:
            print(n0,n1, atoms_in_cell(positions,sccell_a) ) 
        

