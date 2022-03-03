import numpy as np
from ase import Atom, Atoms
import spglib as spg

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
            
        #To be replaced by my own functions
        super().__init__(symbols = atomic_symbols,
                         scaled_positions=fractional_positions,
                         cell = cell);

        self.symm_dataset  = None;
        
    def _spglib_structure(self):
        structure = (self.get_cell(),
                     self.get_scaled_positions(),
                     self.get_atomic_numbers() 
                    );
        return structure;

    def symmetries(self, symprec=1e-2):
        """Determine the structure's symmetry operations.

        Returns:
        ------------
        Dictionary:
            The dictionary consist of three keys \"rotation\", \"translations\", and \"equivalent_atoms\".
            The first two correspond to combined rotation and translation operator that generate one of the symmetries
            of the crystal, while the third item indicates the inequivalent atoms that cannot be build using symmetries. 
            
            This function returns the same output `spglib.get_symmetry <https://spglib.github.io/spglib/python-spglib.html#get-symmetry>`_
    """
        return spg.get_symmetry(self._spglib_structure(), symprec=symprec);

    def spacegroup(self, symprec=1e-2):
        """Return the structure's spacegroup as a string.
        """
        print("STADANRD CELL",spg.standardize_cell(self._spglib_structure(), symprec=symprec)
)
        return spg.get_spacegroup(self._spglib_structure(), symprec=symprec)

    def hall_number(self, symprec=1e-2):
        """Return the structure's hall_number as a string.
        """
        if self.symm_dataset is None:
            self.symm_dataset = spg.get_symmetry_dataset(self._spglib_structure(), symprec=symprec);
        return self.symm_dataset["hall_number"]


    def hm_symbol(self, symprec=1e-2):
        """Return the structure's  (full) Hermann-Mauguin symbol.
        """
        if self.symm_dataset is None:
            self.symm_dataset = spg.get_symmetry_dataset(self._spglib_structure(), symprec=symprec);
        
        spacegroup =spg.get_spacegroup_type( self.hall_number(self) );
        return spacegroup["international_full"]
