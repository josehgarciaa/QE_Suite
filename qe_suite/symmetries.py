import seekpath as seek
import spglib as spg
import numpy as np
import json

#retrieve information from the current path
import inspect
fname = __name__.split(".")[-1]+".py"
curr_path = inspect.getfile(inspect.currentframe()).replace(fname,"");
spg_hall_numbers =None;#obtained form http://cci.lbl.gov/sginfo/itvb_2001_table_a1427_hall_symbols.html
with open(curr_path+'spg_hall_numbers.json') as fp:
    spg_hall_numbers= json.loads(fp.read());
spgnum2HM ={ v["Number"].split(":")[0]:v["Hermann-Mauguin"] for k,v in spg_hall_numbers.items() }


def informations(structure, symprec):
    symm_dataset= spg.get_symmetry_dataset(structure, symprec);
    symm_dataset.update(spg_hall_numbers[str(symm_dataset["hall_number"])]);
    assert symm_dataset["number"]!=symm_dataset["Number"], print("incompatibility between spglib and hall_numbers");    
    return symm_dataset;

def band_path( structure, symprec ):
    kp_path= seek.get_path( structure , symprec=1e-3);

    #Get the paths compatibles with the periodic boundary conditions
    coords= kp_path['point_coords'];    
    path  = kp_path["path"];
    return coords,path;

def get_crystal_from_spgnum(spgnum):
    if 195 <=spgnum <=230:
        return "cubic";
    if 168 <=spgnum <=194:
        return "hexagonal";
    if 143 <=spgnum <=167:
        return "trigonal"
    if 75 <=spgnum <= 142:
        return "tetragonal"
    if 16 <=spgnum <= 74:
        return "orthorhombic"
    if 3 <=spgnum <= 15:
        return "monoclinic"
    if 1 <=spgnum <= 2:
        return "triclinic"
    return None;

