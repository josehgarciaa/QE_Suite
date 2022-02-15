import format as f

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



import namelists as nls
import cards 
class qe_handler:
    namelists = dict(); 

def load_inputf(inputf):
    """
    A function that reads a quantum espresso input file and returns a qehandler object. 
    Args:
        inputf (str): Pathname to the location of the file.
    return An instance of class qe_handler.  
    """

    #Check if file exists and open it.
    #  Not onde

    with open(inputf) as  file:
        data = f.remove_double(" ",file.read().replace("\t"," ")); #Remove extra spaces
        nl_blocks = { nl: nls.get_namelist_block(nl, data) for nl in nls.namelists };

        c= cards.cards_fromtxt(data)["ATOMIC_SPECIES"];
        print( c )

             #All namelist start with a &
#            if "&" in line and (line in qe_namelist):
#                namelist = line;
#                qe_namelist[namelist] = dict();
                
                #read options
#                option  = next(file).strip()
#                while option!= "/":
#                    key,value = list(map(str.strip, option.split('=')));
#                    qe_namelist[namelist][key]=value;
#                    option  = next(file).strip();
    return 0;
