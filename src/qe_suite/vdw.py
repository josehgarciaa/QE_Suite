import numpy as np
from ase.build.niggli import niggli_reduce_cell
import ase 
from ase import Atoms

from fractions import Fraction


def convert_to_cell(cell):
    return ase.cell.Cell( cell );

def transform_cell(T,cell ):
    return ase.cell.Cell( np.dot(T,cell) );
 
def niggli_cell_2D(cell, eps=1e-2, loop_max=100, transformation=False): #Acta Cryst. (1976). A32, 297 
    icell = np.copy(cell)
    def G_metrix(cell):
        a,b,c= cell;
        A = np.dot(a,a);
        B = np.dot(b,b);
        Y = 2*a.dot(b) ;
        return A,B,Y;

    itnum = 0;
    nTrans= np.eye(3);
    A,B,Y = G_metrix(cell);
    is_reduced = (A<=B+eps and np.abs(Y) <= A +eps and np.abs(Y) <= B+eps  and Y<=0 );
    while(not is_reduced):
        A,B,Y = G_metrix(cell);
        itnum+=1;
        if A>B+eps:
            M = [ [0,1,0],[1,0,0], [0,0,1] ];
            cell  = np.dot(M,cell);
            nTrans= np.dot(M,nTrans);
        if Y>0:
            M = [ [1,0,0],[0,-1,0], [0,0,1] ];
            cell  = np.dot(M,cell);
            nTrans= np.dot(M,nTrans);
        if np.abs(Y) > (A + eps):   
            M = [ [1,0,0],[-np.sign(Y),1,0], [0,0,1] ];
            cell  = np.dot(M,cell);
            nTrans= np.dot(M,nTrans);
        if np.abs(Y) > (B + eps):   
            M = [ [1,-np.sign(Y),0],[0,1,0], [0,0,1] ];
            cell  = np.dot(M,cell);
            nTrans= np.dot(M,nTrans);
        is_reduced = (A<=B+eps and np.abs(Y) <= A +eps and np.abs(Y) <= B+eps  and Y<=0 );
        if itnum > loop_max:
            print("Cell not reduced, returning None");
            return None

    if transformation:
        return cell, nTrans
    return cell;
    



def get_cells_rational_ratio(a_cell,b_cell, max_size=10):
    #The algorithm assumes a large and a small celll
    l_cell  = ase.cell.Cell(a_cell);
    s_cell  = ase.cell.Cell(b_cell);
    inverted= False;
    #Check if the assignation is correc
    if l_cell.area(2)< s_cell.area(2):
        l_cell,s_cell = s_cell,l_cell;
        inverted =True;
    #Compute the raatio of the areas and its approximation
    r_ratio  = s_cell.area(2)/l_cell.area(2);
    a_ratio  = Fraction(r_ratio).limit_denominator(max_size);
    if a_ratio.denominator ==0 or a_ratio.numerator ==0:
        return (0,0),r_ratio,np.inf;
 
    if inverted:
        diff = 1.0 - r_ratio/a_ratio;
        ratio= 1/r_ratio
        return (a_ratio.denominator,a_ratio.numerator), ratio,diff;

    diff = 1 - a_ratio/r_ratio;
    return (a_ratio.numerator,a_ratio.denominator),r_ratio,diff;

def divisors(n):
    return [ i for i in range(1,n+1) if n%i==0 ];

def get_compatible_supercell_transformations(n):
    #Citation: Journal of Applied Physics 55, 378 (1984); doi: 10.1063/1.333084
    #We are building a matrix that construct a super cell with area equal to n times
    #he unit cell
    #    | k  j 0|
    # A =| 0  m 0|
    #    | 0  0 1|
    # that satisfy the constrains
    # k*m = n. Since we want n/m to be an integer, the ms are the divisors of n.
    # k,m > 0
    # 0<= j <= m-1 
    ms = divisors(n);
    ks = [n//m for m in ms];
    js_m= [[j for j in range(m)] for m in ms ];
    Us = [ [[k,j,0],[0,m,0],[0,0,1]] for m,k,js in zip(ms,ks,js_m) for j in js] 
    return np.array(Us, dtype=int);

def get_closest_cells( a_cell, b_cell, max_size=10, tol=1e-2):
    a_cell,b_cell = niggli_cell_2D(a_cell),niggli_cell_2D(b_cell)
    (a_n,b_n) ,ratio,err =  get_cells_rational_ratio( a_cell,b_cell,max_size=max_size );
    if a_n == 0 or b_n ==0:
        return None;

    #Compute the set of possible transformations, its correspinding cells
    a_trans,b_trans  = [ get_compatible_supercell_transformations(n) for n in (a_n,b_n) ]
    a_scells,b_scells= [ np.dot(trans,cell) for trans,cell in zip( (a_trans,b_trans),(a_cell, b_cell))]

    #Compute the maximally reduced niggli cells, make it 2D and try to get the unique transformations
    a_ncells,b_ncells= [ np.array(list(map(niggli_cell_2D,scell))) for scell in (a_scells,b_scells)];
    a_ncells,b_ncells= [ np.unique( np.round(scell,3), axis=0)  for scell in (a_ncells,b_ncells) ]

    #Construct a difference matrix and get the minimum index
    diff_matrix =  np.array([ [ np.linalg.norm( (a_ncell - b_ncell)[:2,:2]) for b_ncell in b_ncells] for a_ncell in a_ncells])
    a_idx, b_idx= np.unravel_index(np.argmin(diff_matrix, axis=None), diff_matrix.shape)

    return (a_ncells[a_idx],b_ncells[b_idx]),diff_matrix[a_idx, b_idx];
    

def get_vdw_cell( a_structure, b_structure, max_strain=0.2, strain_cell="a", max_size=10, max_area=100 ):
    #Brute force scanning of strain
    min_diff,min_ds,min_ab_scells = np.inf,np.inf,None;
    for ds in np.linspace(-max_strain,max_strain,100):
        a_cell,b_cell = a_structure.get_cell(),b_structure.get_cell();

        strain = np.diag([1+ds,1+ds,1]);
        if strain_cell == "a":
            a_cell = strain.dot(a_cell);
        if strain_cell == "b":
            b_cell = strain.dot(b_cell);

        closest_cells = get_closest_cells( a_cell, b_cell, max_size=max_size);
        if closest_cells is None:
            return None;        
        ab_scells, diff= closest_cells;
        if diff <= min_diff and abs(ds)< abs(min_ds):
             min_ab_scells,min_diff,min_ds = ab_scells,diff,ds;

    min_ab_scell = list(map(cell.Cell,min_ab_scells));
    return min_ab_scell, min_diff,min_ds

def match_cells(a_structure, b_structure):
    a_cell = a_structure.get_cell() ;
    b_cell = b_structure.get_cell() ;

    norm   = np.linalg.norm;
    invert = a_cell[2][2]< b_cell[2][2];
    if invert:
        transf = np.diag( norm(b_cell, axis=1)/norm(a_cell, axis=1) );
        a_structure.set_cell(transf.dot(a_cell));
    else:
        transf = np.diag( norm(a_cell, axis=1)/norm(b_cell, axis=1) );
        b_structure.set_cell(transf.dot(b_cell));

    return (a_structure, b_structure), transf;

def get_cell_transformation(ocell,fcell):
    return ocell.scaled_positions(np.array(fcell))

def get_atoms_in_cell(structure,cell,shift=[0,0,0]):
    icell    = structure.get_cell(); 
    numbers  = structure.get_atomic_numbers();
    scal_pos = structure.get_scaled_positions();
    cart_pos = icell.cartesian_positions( scal_pos + shift );
    scal_pos = cell.scaled_positions ( cart_pos  ) ;
    allowed  =  np.all( (scal_pos < [1,1,1])*(scal_pos >= [0,0,0]),axis=1 );
    return list(scal_pos[allowed]), list(numbers[allowed]);

def expand_supercell(structure, sc_cell, nmax=10 ):
    sc_spos,sc_anum = [], []
    pns = np.arange(0 , nmax);
    nns = np.arange(-nmax,0)[::-1];
    for nxs,nys in ( (pns,pns),(pns,nns),(nns,pns),(nns,nns)):
        for nx in nxs:
            natoms= len(sc_spos);
            for ny in nys:
                shift  = [nx,ny,0];
                scell_pos,scell_num= get_atoms_in_cell(structure,sc_cell,shift= shift);
                if len(scell_pos)==0:
                    break;
                sc_spos.append(scell_pos);
                sc_anum.append(scell_num);
            if natoms == len(sc_spos):
                break;
    sc_anum = [n for nums in sc_anum for n in nums];
    sc_spos = [p for spos in sc_spos for p in spos];
    #sort by atomic number from high to low
    indexes = np.argsort(sc_anum)[::-1];
    sc_anum = np.array(sc_anum)[indexes];
    sc_spos = np.array(sc_spos)[indexes];

    return Atoms( numbers=sc_anum, scaled_positions=sc_spos, cell=sc_cell, pbc=structure.pbc );
