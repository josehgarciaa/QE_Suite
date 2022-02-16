import qe_suite.format as f


def format(x):
    
    if isinstance(x, str):
        return "\'"+x+"\'";

    if isinstance(x, bool):
        return ".TRUE." if x else ".FALSE.";

    return str(x);



def load_xyz(fname):
    with open(fname) as file:
        data = f.remove_double(" ", file.read().replace("\t"," ") );
        data = f.remove_empty( data.split("\n") );
        natm = int( data.pop(0) );
        data = [ f.remove_empty(x.split(" ")) for x in data ];
        assert natm == len(data);
        data = [ (s,(float(x),float(y),float(z))) for s,x,y,z in data];


    return data

def load_cell(fname):
    with open(fname) as file:
        data = f.remove_double(" ", file.read().replace("\t"," ") );
        data = f.remove_empty( data.split("\n") );
        data = [ list(map(float,f.remove_empty( x.split(" ")))) for x in data ];
    return data