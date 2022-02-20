from ase.build.niggli import niggli_reduce_cell
from ase import cell
from fractions import Fraction

def get_cells_rational_ratio(a_cell,b_cell, max_size=10):
    #The algorithm assumes a large and a small celll
    l_cell  = cell.Cell(a_cell);
    s_cell  = cell.Cell(b_cell);
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
def get_compatible_supercell_transformations(n):
    ms = divisors(n);
    ks = [n//m for m in ms];
    js_m= [[j for j in range(m)] for m in ms ];
    Us = [ [[k,j,0],[0,m,0],[0,0,1]] for m,k,js in zip(ms,ks,js_m) for j in js] 
    return np.array(Us);

def get_closest_cells( a_cell, b_cell, max_size=10):
    (a_n,b_n) ,ratio,err =  get_cells_rational_ratio( a_cell,b_cell,max_size=max_size);
    if a_n == 0 or b_n ==0:
        return None;

    niggli_cells =[];
    for cell in (a_cell,b_cell):
        comp_cells = get_compatible_supercell_transformations(a_n).dot(a_cell);
        niggli_cells.append([ niggli_reduce_cell(sc)[0] for sc in comp_cells]);
    
    amin,bmin=None,None 
    min_diff = np.inf;
    nig_a, nig_b = np.array(niggli_cells);
    for ia,na in enumerate(nig_a):
        for ib,nb in enumerate(nig_b):
            diff = np.linalg.norm( (na-nb)[:2,:2])
            if diff <min_diff:
                min_diff= diff;
                amin,bmin = (ia,ib)
    return (nig_a[amin],nig_b[bmin]),min_diff;


def get_vdw_cell( a_atoms, b_atoms, max_strain=0.2, strain_cell="a", max_size=10 ):
    a_icell,b_icell = a_atoms,b_atoms;
    min_diff,min_ds,min_scatms = np.inf,np.inf,None;
    
    for ds in np.linspace(-max_strain,max_strain,100):
        a_cell,b_cell = a_icell,b_icell;

        strain = np.diag([1+ds,1+ds,1]);
        if strain_cell == "a":
            a_cell = strain.dot(a_icell);
        if strain_cell == "b":
            b_cell = strain.dot(b_icell);

        closest_cells = get_closest_cells( a_cell, b_cell, max_size=max_size);
        if closest_cells is not None:
            sc_atoms, diff= closest_cells;
            if diff <= min_diff and abs(ds)< abs(min_ds):
                min_scatms, min_diff,min_ds = sc_atoms,diff,ds;

    return min_scatms, min_diff,min_ds

