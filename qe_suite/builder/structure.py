import numpy as np
from ase import Atom, Atoms
from .. import symmetries as symm
from  .. import parse 
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
    def __init__(self, cell=None, fractional_positions=None, atomic_symbols=None, symmetrize=False):
        if all((cell, fractional_positions, atomic_symbols)):
            super().__init__(symbols = atomic_symbols,
                             scaled_positions=fractional_positions,
                             cell = cell)
        self.symprec = 1e-2; 
        self.was_symmetrized = False;
        if symmetrize:
            self.symmetrize();


    def set_atomic_symbols(self, atomic_symbols):
        self.set_atomic_numbers(atomic_symbols);
        return self;

    def get_atomic_positions(self):
        return ( "crystal",[ (s,*x) for s,x in zip(self.get_chemical_symbols(), self.get_scaled_positions())] );


    def get_cell_parameters(self):
        return ("angstrom", list(self.get_cell()) );

    def get_reciprocal_vectors(self):
        return np.linalg.inv(self.get_cell()).T*2*np.pi;


    def get_periodicity(self):
        return self.get_pbc();

    def set_periodicity(self,periodicity):
            self.set_pbc(periodicity)
            return self

    def set_as_2D(self):
        self.set_periodicity((True, True,False))
        return self

    def get_crystallographic_constants(self):
        cell = np.array(self.get_cell());
        norms=np.linalg.norm(cell,axis=1)
        ncell= cell.dot(np.diag(1/norms))
        ncell_params = ncell.dot(ncell.T);
        A,B,C = norms;
        cosAB,cosAC,cosBC = ncell_params[np.triu_indices(3,k=1)];
        return A,B,C,cosAB,cosAC,cosBC;

    def get_kpoints(self,type="automatic", kp_distance=0.25, shifts=[1,1,1]):

        rec_vec_lengths = np.linalg.norm( self.get_reciprocal_vectors(),axis=1 );
        rec_vec_divs    = rec_vec_lengths/kp_distance
        kpoints = np.ceil(rec_vec_divs ).astype(int)
        not_periodic = np.logical_not(self.get_periodicity());
        kpoints[  not_periodic] = 1;

        return [*kpoints,*shifts]

    #COSAS DE SIEMTRIA

    def _get_spglib_structure(self):
        symbols= self.get_chemical_symbols();
        sym2num= { s:i for i,s in enumerate(set(self.symbols)) };
        numbers= [ sym2num[s] for s in self.symbols];
        return (self.get_cell(), self.get_scaled_positions(), numbers);

    def get_symmetry_informations(self):
        #Get the symmetrized structure
        symm_structure = self._get_spglib_structure();
        symm_dataset = symm.informations(symm_structure, symprec=self.symprec);
        num2sym= { i:s for i,s in enumerate(set(self.get_chemical_symbols())) };
        symm_dataset["std_symbols"] = [ num2sym[x] for x in symm_dataset["std_types"] ];
        return symm_dataset;      

    def symmetrize(self):
        symm_structure = self._get_spglib_structure();
        symm_dataset = self.get_symmetry_informations();
        print("The structure was symmetrized to the spacegroup:",symm_dataset["international"])

        #Use it for the lattice
        print(symm_dataset["std_lattice"] )
        self.set_cell( symm_dataset["std_lattice"] )
        self.set_chemical_symbols( symm_dataset[ "std_symbols"] )
        self.set_scaled_positions( symm_dataset["std_positions"] )
        self.was_symmetrized =True;   

        return self;    

    def get_symmetrized_high_symmetry_points(self):
        hsymm_points,path = symm.band_path( self._get_spglib_structure(), self.symprec );
        labels = np.array(list(hsymm_points.keys()))
        coords = np.array(list(hsymm_points.values()));
        #When a system is non periodic in a given direction, any
        #displacement along that direction is forbiden. Therefore
        #we compare the coordinates with zero to identify those
        #compatible with non_periodic directions, i.e, the Gamma
        #point [0,0,0] is always compatible with any non periodic system
        non_pbc_comp_coords = np.isclose(coords,[0,0,0]);

        #We negate the periodicty flags to obtain non periodic directions if any
        #and compare those directions with the compatible coordinates with non_pbc directions
        non_pbc   =  np.logical_not(self.get_periodicity());
        comp_index=  np.all(non_pbc_comp_coords[:,non_pbc],axis=1);
        #The outcomes of this process are all high symmetry points
        #compatible with the periodicity of the system
        hsymm_points = dict(zip(labels[comp_index], coords[comp_index]));
        return hsymm_points;

    def symmetrized_band_path(self, density=None):
        full_hsymm_points,path = symm.band_path( self._get_spglib_structure(), self.symprec );
        hsymm_points = self.get_symmetrized_high_symmetry_points();
        labels = hsymm_points.keys()
        path  = [ (l1,l2)  for l1,l2 in path if (l1 in labels) and (l2 in labels) ]

        density = np.ones(len(path), dtype=int)*20;
        kpath  = { (l1,l2):np.linspace(hsymm_points[l1],hsymm_points[l2],n,endpoint=False )  for (l1,l2),n in zip(path,density) };

        return kpath

    def load_from(self,xyz_file, cell_file, symmetrize=False):
        xyz = parse.load_xyz(xyz_file);
        cell= parse.load_cell(cell_file);
        atomic_symbols, atomic_positions = list(zip(*xyz));
        fractional_positions =  list(np.dot(atomic_positions,np.linalg.inv(cell).T ))
        self.__init__(cell=cell, 
                      fractional_positions=fractional_positions, 
                      atomic_symbols=atomic_symbols,\
                      symmetrize=symmetrize);
        return self
