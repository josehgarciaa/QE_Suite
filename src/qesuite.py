import qe_io
import seekpath 
import spglib as spg
import numpy as np
from ase import Atom
from ase import Atoms


import _namelists.control as control


#c = control.handler();
#print( c.options.keys() )
#q = qes.handler();

#print( q.c.options.keys() )

#import json
#fname ="SSSP_1.1.2_PBEsol_precision.json";
#with open(fname) as f:
#    pseudo_potentials = json.loads(f.read());


class handler():
    def __init__(self):
        x=0
        #self.c = control.handler();

    def set_structure(self, structure  ):
        """
        Set the atomic position and lattice vectors

        Args:
            structure ( tuple) : A tuple containing (cell, positions, numbers)n array containing the lattice vectors in arbitrary units, where cell[0] the first vector.

        """
        return True;

    def set_calculation(self, calc  ):
        """
        Defines the calculation to be performed. This flag will initialize different options required for the calculation
        using the default values

        Args:
            calc ( string) : 
        """
        return 0;

    def use_SSSP(self, type="efficiency", path="."):
        """
        Use the standard Solid State Pseudo Potential (SSSP). The type define whereas you want to use efficiency or efficiency 
        or precision. 
        If you use this option please read the acknowledgment instrunctions in :
        https://www.materialscloud.org/discover/sssp/table/efficiency

        Args:
            type ( string) : 
            path ( string) : 
        """
        return 0;

    def use_symmetries(self):


        
        #def get_ibrav(cell, brav_lat ):
        #    if brav_lat== 'mP': #monoclinic
        #        cdm1 = np.linalg.norm(cell[0] );
        #        cdm2 = np.linalg.norm(cell[1] )/cdm1;
        #        cdm3 = np.linalg.norm(cell[2] )/cdm1;
        #        cdm4 = np.dot( cell[0], cell[1])/cdm1/cdm2;
        #        ibrav= 12


        """
        Analyze the symmetries in the crystal strturess and force the cell to respect it
        """
        return 0;

    def write_input_file(self, ofname ):
        """
        Write the input file for QESpresso
        """
        print("Hola mundo")        

        return 0;





def generate_from_xyz(xyz, cell, magnetic=False, periodic=(True, True,True) ):
    """
    Generate a QESuite handler using the structural informaiton
    in the form of a xyz and cell file. 

    Args:
        cell (string): An array containing the lattice vectors in arbitrary units, where cell[0] the first vector.
        xyz (string): A list of tuples (s,x,y,z), whe s detones atomic specie located at the (x,y,z) position in the same units as cell.
        magnetic (bool): A label to indicates if the system is magnetic. Default False.
        preiodic (bool): A tuple that indicates the periodicity directions. Default (True, True, True ).
    """
    #Convert the xyz tuple into an atom object
    atoms    = Atoms( [ Atom(*a) for a in xyz ], cell=cell ) ;
    structure= (atoms.get_cell(),
                atoms.get_scaled_positions(),
                atoms.get_atomic_numbers() );
    structure= spg.standardize_cell(structure, symprec=1e-3);

    #Initialize the QESuite handler
    qes_h = handler();
    qes_h.set_structure( structure = structure );

    return qes_h;

