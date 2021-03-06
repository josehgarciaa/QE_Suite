import qe_suite.format as f

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


def write_xyz(xyz, ofname):
    with open(ofname, "w") as f:
        f.write( str(len(xyz))+"\n\n")
        for s,p in xyz:
            f.write("{} {} {} {} \n".format(s,*p) )
    return ;
