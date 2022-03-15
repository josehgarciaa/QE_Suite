import seekpath as seek
import spglib as spg
import numpy as np

def informations(structure, symprec):
    return spg.get_symmetry_dataset(structure, symprec);

def band_path( structure, symprec ):
    kp_path= seek.get_path( structure , symprec=1e-3);

    #Get the paths compatibles with the periodic boundary conditions
    coords= kp_path['point_coords'];    
    path  = kp_path["path"];
    return coords,path;

