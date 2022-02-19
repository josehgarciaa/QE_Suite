from matplotlib.pyplot import sca
import seekpath as seek
import spglib as spg
import numpy as np
import qe_suite.constants as const

from  ase.spacegroup import symmetrize

def get_2D_orientation( system ):
    cell = np.array(system.get_cell());
    largest_vector_pos = np.argmax( np.linalg.norm(cell, axis=1) );
    v1,v2 = [x for i,x in enumerate(cell) if i !=largest_vector_pos];
    if largest_vector_pos==1:
        v1,v2= v2,v1; #To preserve the handness
    plane_dir= np.cross(v1,v2);
    return plane_dir/np.linalg.norm(plane_dir);

def is_standard_2D_orientation( system ):
    orientation = get_2D_orientation( system );
    return tuple( np.round(orientation) ) == (0,0,1);

def standard_2D_angles(system):
    up   = get_2D_orientation( system );
    x,y,z= up;
    rp   = np.sqrt( x**2 + y**2);
    theta= np.arctan( rp/z) if z!=0 else np.pi/2;

    sy   =  1 if y>=0 else -1;
    if x>0:
        phi = np.arctan( y/x );
    elif x <0:
        phi = np.arctan( y/x ) + sy*np.pi;
    else:
        phi = sy*np.pi/2; 

    return phi,theta;

def standard_2D_transformation(system):
    phi,theta = standard_2D_angles(system);
    c,s = np.cos(phi), np.sin(phi);
    Rz1 = np.array([ [ c,-s, 0],
                    [ s, c, 0],
                    [ 0, 0, 1] ]);
    c,s = np.cos(theta), np.sin(theta);
    Ry = np.array([ [ c, 0, s],
                    [ 0, 1, 0],
                    [-s, 0, c]]);
    c,s = np.cos(phi-np.pi/2), np.sin(phi-np.pi/2);
    Rz2 = np.array([ [ c,-s, 0],
                    [ s, c, 0],
                    [ 0, 0, 1] ]);

    R    = np.dot(Rz2,np.dot(Ry,Rz1));
    return R;

def standard_2D_form(system):
    print("The cell primitive cell was transformed")
    trans= standard_2D_transformation(system);
    cell = system.get_cell();
    spos = system.get_scaled_positions();
    system.set_cell(np.dot(trans,np.dot( cell, trans.T)), scale_atoms=True);
    print(np.array(cell),"-->", np.array(system.get_cell()) ); 
    if not is_standard_2D_orientation( system ):
        print("The cell was not properly set in the stanarda 2D form")
    return system;

def system_to_structure(system):
    structure = (system.get_cell(),
                 system.get_scaled_positions(),
                 system.get_atomic_numbers() 
                 );
    return structure;

def get_brav_params( system ):

    symm  = symmetrize.refine_symmetry(system);
    kp_path = seek.get_path( system_to_structure(system) , symprec=1e-5);
    brav_lat= kp_path["bravais_lattice"];
    spgnum  = kp_path["spacegroup_number"];
    ibrav = 0;
    spgnum=0;
    celldm= np.zeros(6);

    #For 2D system use orient the plane vector along z
    if ( system.pbc == (True,True,False) ).all(): 
        if not is_standard_2D_orientation( system ):
            system = standard_2D_form(system);

    #Express the cell in bohr
    cell = system.get_cell()*const.Ang2Bohr;
    celldm[0] = np.linalg.norm(cell[0] );

    #Use the bravai lattice information to define ibrav
    if brav_lat== 'hP': #hexagonal
        celldm[2] = np.linalg.norm(cell[2] )/celldm[0];
        ibrav= 4;

    if brav_lat== 'mP': #monoclinic
        celldm[1] = np.linalg.norm(cell[1] )/celldm[0];
        celldm[2] = np.linalg.norm(cell[2] )/celldm[0];
        celldm[3] = np.dot( cell[0], cell[1])/celldm[0]/celldm[2];
        ibrav= 12;

    return system, ibrav, spgnum, celldm;

def get_band_path( system ):

    structure = system_to_structure(system);
    kp_path= seek.get_path( structure , symprec=1e-5);

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
    kpath  = { (l1,l2):np.linspace(coords[l1],coords[l2],n,endpoint=False )  for (l1,l2),n in zip(path,density) };
    return kpath;
