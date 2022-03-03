
#import json
#import numpy as np
#import qe_suite.symmetries as symmetries
#import qe_suite.io as qe_io
#from qe_suite.namelists import control, system, electrons, ions, cell, fcp
#from qe_suite.cards import atomic_positions, atomic_species, kpoints, cell_params, constraints 

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



class Calculation():
    x = 0;


class QEInput():
    """A generator of input files for Quantum Espresso.

    Parameters:
    ------------
    name:
        A  name to identify the system. 
    structure:
        An instance of a :py:class:`~qe_suite.Structure`
    Structure:
        An instance of a :py:class:`~qe_suite.Calculation`
    """

    def __init__(self, name="qe_suite", structure = None, calculation = None, electronic_state = None):
        self.structure = structure;
        "I am here"

class QEInput_():
    """A generator of input files for Quantum Espresso.

    
    """

    def __init__(self, structure = None, calculation = None, electronic_state = None):
        #Default initialization of namelists
        self.c   = control.handler();
        self.s   = system.handler();
        self.e   = electrons.handler();
        self.ions= ions.handler();
        self.cell= cell.handler();
        self.fcp = fcp.handler();
        #Default initialization of cards
        self.ap    = atomic_positions.handler();
        self.ae    = atomic_species.handler();
        self.kpts  = kpoints.handler();
        self.cp    = cell_params.handler();
        self.constr= constraints.handler();
        
        self.structure = None;
        self.pbc = (True, True, True);

    def set_state(self,two_dimensional = False, insultator=False, magnetic=False ):
        if insultator:
            self.s.options["occupations"]="fixed";

    def set_structure(self, structure ):

        self.structure = structure;
        symbols = structure.get_chemical_symbols();
        species = [ (s, Atom(s).mass, s+".UPF") for s in set(symbols)]
        self.s.options["nat"]  = len(symbols);
        self.s.options["ntyp"]= len(species);

        if ( self.structure.pbc == (True,True,False) ).all() :
            self.s.options["assume_isolated"]='2D';
            self.kpts.set_kpoints( self.kpts.get_kpoints(), pbc=self.structure.pbc );

        self.ae.set_atomic_species(species);
        self.ap.set_atomic_positions(self.structure);

        return True;

    def get_structure(self):
        return self.structure;

    def set_calculation(self, calc  ):

        self.c.options["calculation"]= calc;
        return 0;

    def use_symmetries(self):
        structure, ibrav, spgnum, celldm = symmetries.get_brav_params( self.structure );
        self.s.set_bravais_lattice(ibrav, celldm);
        self.set_structure(structure);

        if self.c.options["calculation"]=="bands":
            kpoints = symmetries.get_band_path( self.structure );
            self.kpts.set_kpoints( kptype= "crystal_b", kpoints= kpoints);

        return 0;

    def write(self, ofname="out" ):
        text = self.c.text();
        text += "\n"+self.s.text();
        text += "\n"+self.e.text();

        calc = self.c.options["calculation"];
        if calc =="relax" or calc =="md" or calc =="vc-relax"or calc =="vc-md":
            text += "\n"+self.ions.text();

        if self.c.options["lfcp"]:
            text += "\n"+self.fcp.text();

                
        text += self.ae.text();
        text += self.ap.text();
        text += self.kpts.text();

        with open(ofname, 'w') as f:
            f.write(text);
        return 0;


    def use_SSSP(self, functional="PBEsol", target="precision", path="./SSSP", use_cutoff=True):
        """
        Use the standard Solid State Pseudo Potential (SSSP). The type define whereas you want to use efficiency or efficiency 
        or precision. 
        If you use this option please read the acknowledgment instrunctions in :
        https://www.materialscloud.org/discover/sssp/table/efficiency

        Args:
            type ( string) : 
            path ( string) : 
        """

        self.c.options["pseudo_dir"]=path+"/";

        with open(path+"/versions.yaml") as f:
            version = f.read().split(":")[-1].replace("\'","").replace(" ","").replace("\n","");

        dbfname = path+"/SSSP_"+version+"_"+functional+"_"+target+".json";
        with open(dbfname) as f:
            sp_info = json.loads(f.read());

        ae = self.ae;
        species = [ (s,m ,sp_info[s]["filename"]) for s,m,sp in ae.get_atomic_species() ]
        ae.set_atomic_species(species);
        
        
        if use_cutoff:
            cutoffs=[0,0];
            for s,m,sp in ae.get_atomic_species():
                sp = sp_info[s];
                keys = ( "cutoff_wfc", "cutoff_rho");
                cutoffs = [v if sp[k]<v else sp[k] for k,v in zip( keys,cutoffs )];
                self.s.set_cutoff(*cutoffs);

        return 0;



def generate_from_xyz(xyz, cell, two_dimensional = False, magnetic=False ):
#    Generate a QESuite handler using the structural informaiton
#    in the form of a xyz and cell file. 
#
#    Args:
#        cell (string): An array containing the lattice vectors in arbitrary units, where cell[0] the first vector.
#        xyz (string): A list of tuples (s,x,y,z), whe s detones atomic specie located at the (x,y,z) position in the same units as cell.
#        magnetic (bool): A label to indicates if the system is magnetic. Default False.
#        preiodic (bool): A tuple that indicates the periodicity directions. Default (True, True, True ).
#    """
    #Convert the xyz tuple into an atom object
    structure    = Atoms( [ Atom(*a) for a in xyz ] ) ;
    structure.set_cell(cell);

    if two_dimensional:
        structure.set_pbc( (True, True, False) );

    #Initialize the QESuite handler
    qes_h = handler();
    qes_h.set_structure( structure = structure, two_dimensional=two_dimensional );

    return qes_h;


