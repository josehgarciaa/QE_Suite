import numpy as np
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
    return np.array(Us);

def get_closest_cells( a_cell, b_cell, max_size=10):


    (a_n,b_n) ,ratio,err =  get_cells_rational_ratio( a_cell,b_cell,max_size=max_size );
    if a_n == 0 or b_n ==0:
        return None;

    niggli_cells =[];
    for n,cell in zip( (a_n,b_n), (a_cell,b_cell) ):
        comp_cells = get_compatible_supercell_transformations(n).dot(cell);
        niggli_cells.append([ niggli_reduce_cell(sc)[0] for sc in comp_cells]);

    amin,bmin=None,None 
    min_diff = np.inf;
    nig_a, nig_b = [ np.array(x,dtype=float) for x in niggli_cells];
    for ia,na in enumerate(nig_a):
        for ib,nb in enumerate(nig_b):
            diff = np.linalg.norm( (na-nb)[:2,:2]);
            if diff < min_diff:
                min_diff = diff;
                amin,bmin= (ia,ib)
    return (nig_a[amin],nig_b[bmin]),min_diff;


def get_vdw_cell( a_structure, b_structure, max_strain=0.2, strain_cell="a", max_size=10 ):
    a_icell,b_icell = a_structure.get_cell(),b_structure.get_cell();
    min_diff,min_ds,min_scatms = np.inf,np.inf,None;
    
    #Brute force scanning of strain
    for ds in np.linspace(-max_strain,max_strain,100):
        a_cell,b_cell = a_icell,b_icell;

        strain = np.diag([1+ds,1+ds,1]);
        if strain_cell == "a":
            a_cell = strain.dot(a_cell);
        if strain_cell == "b":
            b_cell = strain.dot(b_cell);

        closest_cells = get_closest_cells( a_cell, b_cell, max_size=max_size);
        if closest_cells is None:
            return None;
        
        sc_atoms, diff= closest_cells;
        if diff <= min_diff and abs(ds)< abs(min_ds):
            min_scatms, min_diff,min_ds = sc_atoms,diff,ds;

    a_cell , b_cell = min_scatms
    return (cell.Cell(a_cell),cell.Cell(b_cell)), min_diff,min_ds

def match_cells(a_structure, b_structure):
    a_cell = a_structure.get_cell() ;
    b_cell = a_structure.get_cell() ;
    scal   = np.linalg.norm(a_cell, axis=0)/np.linalg.norm(b_cell, axis=0);
    b_structure.set_cell(b_cell*scal);
    return a_structure, b_structure;


def get_atoms_in_cell(structure,cell,shift=[0,0,0]):
    icell    = structure.get_cell(); 
    numbers  = structure.get_atomic_numbers();
    scal_pos = structure.get_scaled_positions();
    cart_pos = icell.cartesian_positions( scal_pos + shift );
    scal_pos = cell.scaled_positions ( cart_pos  ) ;
    allowed  =  np.all( (scal_pos < [1,1,1])*(scal_pos >= [0,0,0]),axis=1 );
    return list(scal_pos[allowed]), list(numbers[allowed]);

def create_expand_supercell(structure, sc_cell ):
    sc_spos,sc_anum = [], []
    zero= np.array([0,0,0],dtype=int);  
    xvec= np.copy(zero);
    while( True ):
        natoms= len(sc_spos);
        yvec  = np.copy(zero);
        while( True  ):
            scell_pos,scell_num= get_atoms_in_cell(structure,sc_cell,shift= xvec + yvec);
            if len(scell_pos)==0:
                break;
            sc_spos.append(scell_pos);
            sc_anum.append(scell_num);
            yvec += [0,1,0];
        if natoms == len(sc_spos):
            break;
        xvec += [1,0,0];

    sc_anum = [n for nums in sc_anum for n in nums];
    sc_spos = [p for spos in sc_spos for p in spos];
    return Atoms( numbers=sc_anum, scaled_positions=sc_spos, cell=sc_cell, pbc=structure.pbc );
