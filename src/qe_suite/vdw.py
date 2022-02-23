
from numpy.linalg import norm
import numpy as np
from ase.build.niggli import niggli_reduce_cell
import ase 
from ase import Atoms

from fractions import Fraction
import warnings
warnings.simplefilter("error")

def cell_param(cell, axis=None):
    norms = np.linalg.norm(cell, axis=1 );
    ucell = np.diag(1/norms).dot(cell);
    inners= np.dot(ucell, ucell.T);
    angles= np.arccos(inners[ np.mask_indices(3, np.triu, 1) ]) #Get only upper elements    
    if axis is None:
        return np.array( ( *norms, *angles ) );
    else:
        return np.array(( *norms[norms!= norms[axis]],angles[axis-2]));


def cell_diff(cell_a, cell_b, axis=None):
#    d1 = np.abs(1- cell_param(cell_a, axis=axis)/cell_param(cell_b, axis=axis) );
#    d2 = np.abs(1- cell_param(cell_b, axis=axis)/cell_param(cell_a, axis=axis) );
#    rel_diff = np.linalg.norm( (d1+d2) )
    rel_diff =np.linalg.norm( cell_a - cell_b)
    return rel_diff



def convert_to_cell(cell):
    return ase.cell.Cell( cell );

def transform_cell(T,cell ):
    return ase.cell.Cell( np.dot(T,cell) );
 
def niggli_cell_2D(cell, eps=1e-1, loop_max=100, transformation=False): #Acta Cryst. (1976). A32, 297 
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
        if np.abs(Y) > A + eps:   
            M = [ [1,0,0],[-np.sign(Y),1,0], [0,0,1] ];
            cell  = np.dot(M,cell);
            nTrans= np.dot(M,nTrans);
        if np.abs(Y) > B + eps:   
            M = [ [1,-np.sign(Y),0],[0,1,0], [0,0,1] ];
            cell  = np.dot(M,cell);
            nTrans= np.dot(M,nTrans);
        is_reduced = (A<=B+eps and np.abs(Y) <= A +eps and np.abs(Y) <= B+eps  and Y<=0 );
        if itnum > loop_max:
            print("Cell not reduced, returning None");
            return None

    if transformation:
        return ase.cell.Cell(cell), nTrans
    return ase.cell.Cell(cell);
    

def best_unit_fraction_approximation( ratio):
    max_n= int( np.ceil(1/ratio) ); tol_n= 2;
    min_n = 1 if max_n-tol_n<1 else max_n-tol_n;
    max_n = max_n+tol_n;
    ns   = np.arange( min_n ,max_n + 1);
    n    = ns[ np.argmin( np.abs(ratio -1/ns)) ];
    return n;



def divisors(n):
    return [ i for i in range(1,n+1) if n%i==0 ];

def get_compatible_supercell_transformations( initial_area, target_area ):
    assert initial_area<= target_area, "Initial area should be smaller tan target area"
    r_ratio  = initial_area/target_area;
    n = best_unit_fraction_approximation( r_ratio)

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

def get_closest_cells( a_cell, b_cell, max_area=None, tol=1e-2):
    a_cell,b_cell = niggli_cell_2D(a_cell),niggli_cell_2D(b_cell)
    a_area,b_area = [a.area(2) for a in (a_cell,b_cell) ];
    if max_area is None:
        max_area = a_area; 

    #Compute the set of possible transformations and return None if any of the cells does not have
    #a transformation
    a_trans,b_trans  = [ get_compatible_supercell_transformations(  initial_area=area, 
                                                                    target_area=max_area)   
                                                                    for area in (a_area,b_area)]
    if a_trans is None or b_trans is None:
        return None;

    #Compute the supercell cells and  reduce them to niggli form
    a_scells,b_scells= [ np.dot(trans,cell) for trans,cell in zip( (a_trans,b_trans),(a_cell, b_cell)) ]

    #Compute the maximally reduced niggli cells, make it 2D and try to get the unique transformations
    a_ncells,b_ncells= [ [ncell for ncell in map(niggli_cell_2D,scell) if ncell is not None ] for scell in (a_scells,b_scells)];

    a_ncells,b_ncells= [ np.unique( np.round(scell,3), axis=0)  for scell in (a_ncells,b_ncells) ]

    #Construct a difference matrix and get the minimum index
#    diff_matrix =  np.array([ [ np.linalg.norm( (a_ncell-b_ncell)[:2,:2] ) for b_ncell in b_ncells] for a_ncell in a_ncells])
    diff_matrix =  np.array([ [ cell_diff(a_ncell,b_ncell, axis=2) for b_ncell in b_ncells] for a_ncell in a_ncells])
    a_idx, b_idx= np.unravel_index(np.argmin(diff_matrix, axis=None), diff_matrix.shape)

    rel_diff = diff_matrix[a_idx, b_idx];
    return (a_ncells[a_idx],b_ncells[b_idx]),rel_diff;



def get_compatible_supercell_transformations_old(n):
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

    

def get_vdw_cell( a_structure, b_structure, max_strain=0.01, strain_cell="b", max_area=None ):
#    options = (popsize=100, mutation=(0.7, 1), recombination=0.5,workers=-1);
    if strain_cell=="a":
        #invert since the algorithm will always strain b
        a_structure, b_structure= b_structure, a_structure

    from scipy.optimize import differential_evolution
    a_cell,b_cell = a_structure.get_cell(),b_structure.get_cell();
    min_area = (1+max_strain)**2*np.max( [a_cell.area(2),b_cell.area(2)]);
    def optimize_cell(params, get_cells=False):
        rel_strain,area = params;
        strain  = np.diag([1+rel_strain,1+rel_strain,1]);
        closest_cells = get_closest_cells( a_cell, strain.dot(b_cell), area);
        if closest_cells is None:
            return np.inf;        
        if get_cells:
            return closest_cells;
        ab_scells, diff= closest_cells
        return diff;
    bounds = [(-max_strain,max_strain), (min_area,max_area)];
    opt_params = differential_evolution(optimize_cell, bounds, init="sobol");
    (opt_a,opt_b),min_diff = optimize_cell(opt_params.x, get_cells=True);

    if strain_cell == "a":
        #invert the resulting cell since the algorithm assumed strained b
        opt_a,opt_b = opt_b,opt_a;

    return (convert_to_cell(opt_a),convert_to_cell(opt_b)), opt_params.x,min_diff

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
