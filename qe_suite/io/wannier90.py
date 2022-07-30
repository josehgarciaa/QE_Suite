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
            out += "{:.6f} {:.6f} {:.6f}\n".format(*x);
        out+="end unit_cell_cart\n"
        return out

    if k.strip() == "atoms_frac":
        out = "begin atoms_frac\n"
        for x in v:
            out += "{} {:.6f} {:.6f} {:.6f}\n".format(x[0],*x[1]);
        out+="end atoms_frac\n"
        return out

    if k.strip() == "kpoints":
        out = "begin kpoints\n"
        for x in v:
            out += "{:.8f} {:.8f} {:.8f}\n".format(*x);
        out+="end kpoints\n"
        return out

    if k.strip() == "projections":
        out = "\nbegin projections\n"
        if v is not "":
            out += v;
        out+="end projections\n"
        return out


    if k.strip() == "mp_grid":
        out = k.strip();
        out += "{} {} {}\n".format(*v);
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

