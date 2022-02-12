import spglib as sp
import seekpath

import format as f
import json
def load_xyz(fname):
    with open(fname) as file:
        data = f.remove_double(" ", file.read().replace("\t"," ") );
        data = f.remove_empty( data.split("\n") );
        natm = int( data.pop(0) );
        data = [ f.remove_empty(x.split(" ")) for x in data ];
        assert natm == len(data);
        data = [ (s,(float(x),float(y),float(z))) for s,x,y,z in data];


    return data

def load_lattice(fname):
    with open(fname) as file:
        data = f.remove_double(" ", file.read().replace("\t"," ") );
        data = f.remove_empty( data.split("\n") );
        data = [ list(map(float,f.remove_empty( x.split(" ")))) for x in data ];
    return data


fname = "Pd2Se4.xyz";
xyz_data= load_xyz(fname);
symbols, positions = zip(*xyz_data);

fname = "Pd2Se4.uc";
lattice = load_lattice(fname);


import json
import numpy as np
with open('Pd2Se4.json') as json_file:
    data = json.load(json_file)
data = data[ str(*data["ids"]) ];


lattice = None;
for k,v in data["cell"]["array"].items():
    if k == "__ndarray__":
       shape, dtype, values = v;
       lattice = np.array(values, dtype=dtype).reshape(shape)

tofrac = np.linalg.inv(lattice);
positions = None;
for k,v in data["positions"].items():
    if k == "__ndarray__":
       shape, dtype, values = v;
       positions = np.array(values, dtype=dtype).reshape(shape);
       positions = positions.dot(tofrac);

numbers = None;
for k,v in data["numbers"].items():
    if k == "__ndarray__":
       shape, dtype, values = v;
       numbers = np.array(values, dtype=dtype).reshape(shape)

cell = (lattice, positions, numbers);
spacegroup = sp.get_spacegroup(cell, symprec=1e-5)

symm_and_path = seekpath.get_path(cell, with_time_reversal=True, recipe="hpkot",  symprec=1e-5);

path_point_coords = symm_and_path["point_coords"];
path =  symm_and_path["path"];


path_point_coords = { k:v for k,v in path_point_coords.items() if v[2]==0 };

points = path_point_coords.keys();
path   = [ (x,y) for x,y in path if x in points and y in points ]

print(path)