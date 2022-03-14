import numpy as np
from ase import Atom, Atoms
from .. import symmetries as symm

class Structure(Atoms):
    """An atomic structure.

    The atomic structure defines the type of atoms and their position within a unit cell.
    
    Parameters:
    ------------
    cell:
        array of shape (3,3) containing the lattice vectors 
        in row_major form, i.e, cell[0] corresponds to the first lattice vector
    fractional_positions: 
        array of shape (n,3) containing the three-dimensional coordinates of the n atoms
        in fractional coordinates
    atomic_symbols: 
        array of shape (n,1) containing the atomic symbols of the n atoms in the same order as 
        fractional_positions

    Examples:

    >>> import math as m
    >>> from  qe_suite.builder import Structure 
    >>>
    >>> a= 0.246 ;  # lattice constant in nm which is the default unit. 
    >>> cell = [ [a,0,0], [-a/2,a*m.sqrt(3)/2,0], [0,0,10*a]] ;
    >>> fractional_positions = [ [0,0,0], [2/3,1/3,0] ] ;
    >>> atomic_symbols = ['C', 'C'] ;
    >>> s = Structure(cell, fractional_positions, atomic_symbols) ;
    
    """
    def __init__(self, cell, fractional_positions, atomic_symbols):
        
        #Check if inputs are arrays
        natm  = len(atomic_symbols);
        arrays = [cell,fractional_positions, atomic_symbols];
        shapes= ( (3,3), (natm,3), (natm,1) );
        for i,(a,s) in enumerate(zip(arrays,shapes)):
            try:
                arrays[i] = np.array(a);
            except TypeError:
                print("The input", a, " in Structure does not have the proper type ")
                raise
            arrays[i].reshape(s);

        self.periodicity= (True, True, True)            
        #To be replaced by my own functions
        super().__init__(symbols = atomic_symbols,
                         scaled_positions=fractional_positions,
                         cell = cell,
                         pbc = self.periodicity);
        self.symm_dataset  = None;
    


    def get_atomic_positions(self):
        return ( "crystal",[ (s,*x) for s,x in zip(self.get_chemical_symbols(), self.get_scaled_positions())] );


    def get_cell_parameters(self):
        return ("angstrom", list(self.get_cell()) );

    def get_reciprocal_vectors(self):
        return np.linalg.inv(self.get_cell()).T*2*np.pi;



    #COSAS DE SIEMTRIA
    def get_symmetry_informations(self, symprec=1e-2):
        symbols= self.get_chemical_symbols();
        sym2num= { s:i for i,s in enumerate(set(self.symbols)) };
        num2sym= { i:s for i,s in enumerate(set(self.symbols)) };
        numbers= [ sym2num[s] for s in self.symbols];
        #Get the symmetrized structure
        symm_structure = (self.get_cell(), self.get_scaled_positions(), numbers);

        symm_dataset = symm.informations(symm_structure, symprec=1e-2);
        symm_dataset["std_symbols"] = [ num2sym[x] for x in symm_dataset["std_types"] ];

        return symm_dataset;      

    def get_periodicity(self):
        return self.get_pbc();

    def set_periodicity(self,periodicity):
            self.set_pbc(periodicity)
            return self

    def set_as_2D(self):
        self.set_periodicity((True, True,False))
        return self

    def symmetrize(self, symprec=1e-2):
        symm_dataset = self.get_symmetry_informations(symprec=1e-2);
        print("The structure was symmetrized to the spacegroup:",symm_dataset["international"])

        #Use it for the lattice
        self.set_cell( symm_dataset["std_lattice"] )
        self.set_chemical_symbols( symm_dataset[ "std_symbols"] )
        self.set_scaled_positions( symm_dataset["std_positions"] )

        return self;    


    def get_kpoints(self,type="automatic", kp_distance=0.15, shifts=[1,1,1]):

        rec_vec_lengths = np.linalg.norm( self.get_reciprocal_vectors(),axis=1 );
        rec_vec_divs    = rec_vec_lengths/kp_distance
        kpoints = np.ceil(rec_vec_divs ).astype(int)
        not_periodic = np.logical_not(self.get_periodicity());
        kpoints[  not_periodic] = 1;

        return [*kpoints,*shifts]


    def symmetrized_band_path(self, symprec=1e-2):
        symbols= self.get_chemical_symbols();
        sym2num= { s:i for i,s in enumerate(set(self.symbols)) };
        numbers= [ sym2num[s] for s in self.symbols];
        self.symmetrize();
        symm_structure = (self.get_cell(), self.get_scaled_positions(), numbers);
        kpath = symm.band_path( symm_structure, symprec )

        path_norm = np.linalg.norm( np.array(list(kpath.values())),axis=1);
        print("kpath")

        return kpath

    

