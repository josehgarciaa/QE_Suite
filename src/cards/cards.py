import format as f

import _cards.specie as specie
class kpoints:
    """
    K_POINTS { tpiba | automatic | crystal | gamma | tpiba_b | crystal_b | tpiba_c | crystal_c }
    if (gamma)
    nothing to read
    if (automatic)
    nk1, nk2, nk3, k1, k2, k3
    if (not automatic)
    nks
    xk_x, xk_y, xk_z,  wk
    if (tpipa_b or crystal_b in a 'bands' calculation) see Doc/brillouin_zones.pdf
    """
    x=0;

class cell:
    """
    [ CELL_PARAMETERS { alat | bohr | angstrom }
    v1(1) v1(2) v1(3)
    v2(1) v2(2) v2(3)
    v3(1) v3(2) v3(3) ]
    """
    x=0;

class positions:
    x=0;

class occupations:
    x=0;

class constraints:
    x=0;

class velocities:
    x=0;

class forces:
    x=0;




def species_fromtxt(data):
    s= specie.specie;
    return [s(*l.split(" ")) for l in f.remove_empty(data.split("\n") ) ]
    
def positions_fromtxt(data):
    positions = f.remove_empty(data.split("\n") );
    coord= positions.pop(0);
    return [p.split(" ") for p in positions ];

def kpoints_fromtxt(data):
    kpts   = f.remove_empty(data.split("\n") );
    kp_type= kpts.pop(0);
    return [kp.split(" ") for kp in kpts ];

def cell_fromtxt(data):
    vectors = f.remove_empty(data.split("\n") );
    coord= vectors.pop(0);
    return [v.split(" ") for v in vectors ];

def occupations_fromtxt(data):
    return None;

def constraints_fromtxt(data):
    return None;

def velocities_fromtxt(data):
    return None;

def forces_fromtxt(data):
    return None;

def cards_fromtxt(data):

    cards = { "ATOMIC_SPECIES":species_fromtxt,
            "ATOMIC_POSITIONS":positions_fromtxt,
            "K_POINTS":kpoints_fromtxt,
            "CELL_PARAMETERS":cell_fromtxt,
            "OCCUPATIONS":occupations_fromtxt,
            "CONSTRAINTS":constraints_fromtxt,
            "ATOMIC_VELOCITIES":velocities_fromtxt,
            "ATOMIC_FORCES":forces_fromtxt
            };

    cards_begpos =sorted([ (data.find(k), k) for k in cards.keys() if k in data]);
    begpos, cardnames = list(zip(*cards_begpos));
    begpos = list(begpos);
    endpos = list(begpos)[1:] + [-1] ;
    cards  = { c:cards[c]( data[i+len(c):f] ) for c,i,f in zip(cardnames, begpos, endpos) } ;
    return cards;


