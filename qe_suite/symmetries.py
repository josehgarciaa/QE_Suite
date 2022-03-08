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
    labels= coords.keys();
    path  = [ (l1,l2)  for l1,l2 in path if (l1 in labels) and (l2 in labels) ]
    density = np.ones(len(path), dtype=int)*20;
    kpath  = { (l1,l2):np.linspace(coords[l1],coords[l2],n,endpoint=False )  for (l1,l2),n in zip(path,density) };
    return kpath;
