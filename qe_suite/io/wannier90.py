import qe_suite.format as f
import numpy as np

SIG_DIGITS = 6

def key_value_format(k,v):
    if not isinstance(k, str):
        raise ValueError("Error value in key_format")
    k = k.lower();

    if k.strip() == "unit_cell_cart":
        out = "begin unit_cell_cart\n"
        for x in v:
            x = np.round(x,SIG_DIGITS)
            out += "{} {} {}\n".format(*np.round(x,SIG_DIGITS));
        out+="end unit_cell_cart\n"
        return out

    if k.strip() == "atoms_frac":
        out = "begin atoms_frac\n"
        for x in v:
            out += "{} {} {} {}\n".format(x[0],*np.round(x[1],SIG_DIGITS));
        out+="end atoms_frac\n"
        return out

    if k.strip() == "kpoints":
        out = "begin kpoints\n"
        for x in v:
            out += "{} {} {}\n".format(*np.round(x,SIG_DIGITS));
        out+="end kpoints\n"
        return out
    return k+"="+str(v)



def key_format(x):
        
    if not isinstance(x, str):
        raise ValueError("Error value in key_format")

    return x.upper();

def format(x):
        
    if isinstance(x, str):
        return "\'"+x+"\'";

    if isinstance(x, bool):
        return "true" if x else "false";

    try: 
        x = str(x)
    except ValueError:
        print(x," is not a valid string") 

    return str(x);

