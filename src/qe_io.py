import format 



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
        data = format.remove_double(" ",file.read().replace("\t"," ")); #Remove extra spaces
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

def write_QEnamelist( qe_namelist, ofile= "test"):
    with open( ofile, 'w') as file:
        for namespace, options in qe_namelist.items():
            if options != None:    
                file.write(namespace+"\n")
                for option, value in options.items():
                    file.write(" "+option+"="+str(value)+"\n")
                file.write("/ \n")


def read_QEnamelist(inputf):
    qe_namelist = {"&CONTROL":None, "&SYSTEM":None, "&ELECTRONS":None, "&IONS":None, "&CELL":None}

    with open(inputf) as  file:
        for line in file:
            line = line.strip();
            
            #All namelist start with a &
            if "&" in line and (line in qe_namelist):
                namelist = line;
                qe_namelist[namelist] = dict();
                
                #read options
                option  = next(file).strip()
                while option!= "/":
                    key,value = list(map(str.strip, option.split('=')));
                    qe_namelist[namelist][key]=value;
                    option  = next(file).strip();
    return qe_namelist

def write_atomic_specties(atm_spec,inputf):
    with open(inputf, "a") as file:
        file.write("ATOMIC_SPECIES\n")
        for a, m, pseudo in atm_spec:
            file.write(a+" "+m+" "+pseudo+"\n")

def read_atomic_species(inputf):
    atm_spec = [];
    with open(inputf) as f:
        for line in f:
            if "atomic_species" in line.lower():
                items = next(f).strip().split(" ");
                items = [x for x in items if x!=""]
                while len(items)== 3:
                    atm_spec.append(items);
                    items = next(f).strip().split(" ");
                    items = [x for x in items if x!=""]

                return atm_spec;
    print("atomic positions not found. Returning 0");
    return 0;

def write_atomic_positions(atm_pos,inputf):
    with open(inputf, "a") as file:
        file.write("ATOMIC_POSITIONS crystal\n")
        for p in atm_pos:
            file.write("{} {} {} {}\n".format(*p))

def read_atomic_positions(inputf):
    atm_pos = [];
    with open(inputf) as f:
        for line in f:
            if "atomic_positions" in line.lower():
                print(line)
                items = next(f).strip().split(" ");
                items = [x for x in items if x!=""]

                while len(items)== 4:
                    atm_pos.append(items);
                    items = next(f).strip().split(" ");
                    items = [x for x in items if x!=""]
                return atm_pos;
    print("atomic positions not found. Returning 0");
    return 0;

def read_kpoints(inputf):
    with open(inputf) as f:
        for line in f:
            if ("K_POINTS" in line) and ("automatic" in line):
                kpoints = next(f)
                kpoints = np.array(list(map(int,kpoints.strip().split(" "))));
                return kpoints;
        print("the K_POINTS atomatic flat not found. Returning zero")
        return 0;
        
def write_kpoints(kpoints,inputf):
    with open(inputf, "a") as file:
        file.write("K_POINTS automatic\n")
        file.write("{} {} {} {} {} {}\n".format(*kpoints));
        
def mod_kpoints(kpoint, dkp,shift):
    kp = kpoint[:3];
    ks = kpoint[3:]*0 +shift; 

    if (kp!=1).all():
        print("Gamma point calculation not supported, change your kpoints beyond 1")
        return kpoint;
    
    kp[kp!=1]+= dkp;
    return np.array([*kp,*ks]);

def best_kpoints( kpoints, conv_kps ):
    conv_kps = np.array(conv_kps);
    idx = np.argmin(conv_kps[:,0]);
    niter = conv_kps[idx][0];
    ciks = conv_kps[idx][1];
    cink = conv_kps[idx][2];
    ink = kpoints[:3];
    ink = [ ( cink if nk!= 1 else 1) for nk in ink ];
    iks = np.array([1,1,1]); iks.fill(ciks)
    return [*ink,*iks];

def read_finalenergy(inputf):
    energy = None;
    with open(inputf) as  file:
        for line in file:
            if "!    total energy" in line:
                key, value = line.split("=");
                energy = float( value.replace("Ry","").strip() ); 
    return energy;

def read_finalocurrence(label, inputf):
    obs = None;
    with open(inputf) as  file:
        for line in file:
            if label in line:
                value = line.split("=")[1];
                obs = value.strip().split(" ")[0];
                obs = float( obs ); 
    return obs;


def is_converged( array, tol=1e-5):
    if len(array)>=2:
        rel = np.abs( (array[-1]-array[-2])/array[-2]);
        if rel<= tol:
            return True;
    return False;

def final_diff(array):
    array = np.array(array);
    if len(array)>=2:
        return np.abs( (array[-1]-array[-2]) );
    return 0;

def plot_conv(x,y, outputf):
    x=np.array(x);
    y=np.array(y);
    num_iter = len(y);
    x=x[1:num_iter];

    rel_chg = np.abs(np.diff(y)/y[:-1]);
    plt.plot(x,rel_chg,  "-o" ); 
    plt.gca().set_yscale('log')
    plt.savefig(outputf);
    return plt.gcf();


