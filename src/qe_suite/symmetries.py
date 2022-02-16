import seekpath as seek
import spglib as spg
import numpy as np

def get_brav_params( system ):
    ibrav = 0;
    spgnum=0;
    celldm= np.zeros(6);

    structure = (system.get_cell(),
                 system.get_scaled_positions(),
                 system.get_atomic_numbers() 
                 );
    structure = spg.standardize_cell(structure, symprec=1e-2);
    kp_path = seek.get_path( structure , symprec=1e-5);
    cell,spos,anum = structure;
    brav_lat= kp_path["bravais_lattice"];
    spgnum  = kp_path["spacegroup_number"];

    system.set_cell(cell);
    system.set_scaled_positions(spos);
    system.set_atomic_numbers(anum); 
    cell = system.get_cell();

    celldm[0] = np.linalg.norm(cell[0] );
    if brav_lat== 'hP': #monoclinic
        celldm[2] = np.linalg.norm(cell[2] )/celldm[0];
        ibrav= 4;

    if brav_lat== 'mP': #monoclinic
        celldm[0] = np.linalg.norm(cell[0] );
        celldm[1] = np.linalg.norm(cell[1] )/celldm[0];
        celldm[2] = np.linalg.norm(cell[2] )/celldm[0];
        celldm[3] = np.dot( cell[0], cell[2])/celldm[0]/celldm[2];
        ibrav=-12;

    return system, ibrav, spgnum, celldm;