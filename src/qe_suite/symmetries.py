import seekpath as seek
import spglib as spg
import numpy as np


def system_to_structure(system):
    structure = (system.get_cell(),
                 system.get_scaled_positions(),
                 system.get_atomic_numbers() 
                 );
    return structure;

def get_brav_params( system ):
    ibrav = 0;
    spgnum=0;
    celldm= np.zeros(6);
    structure = system_to_structure(system);
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

def get_band_path( system ):
    structure= system_to_structure(system);
    structure= spg.standardize_cell(structure, symprec=1e-2);
    kp_path  = seek.get_path( structure , symprec=1e-5);

    #Get the paths compatibles with the periodic boundary conditions
    coords= dict();    
    for k,v in kp_path['point_coords'].items(): 
        # If the coordinate is not zero or zero and periodic is True, then should be true
        # If the coordinate is not zero and periodic is False, then should be False
        # If the coordinate is  zero and periodic is False, then should be True
        keep_coord = np.all( [ ( x==0.0 or per ) for x, per in zip(v, system.pbc) ] );
        if keep_coord:
            coords[k]=v;

    path  = kp_path["path"];
    labels= coords.keys();
    path  = [ (l1,l2)  for l1,l2 in path if (l1 in labels) and (l2 in labels) ]

    density = np.ones(len(path), dtype=int)*20;
    kpath  = { (l1,l2):np.linspace(coords[l1],coords[l2],n,endpoint=False )  for (l1,l2),n in zip(path,density) }
    return kpath;