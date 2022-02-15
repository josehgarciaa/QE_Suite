import qe_io
import namelists.control
import namelists.system
import namelists.electrons
import namelists.ions
import namelists.cell
import namelists.fcp
import cards.atomic_species
import cards.kpoints
import cards.cell_params
import cards.constraints

from ase import Atom
from ase import Atoms


class handler():
    def __init__(self):
        #Default initialization of namelists
        self.c   = namelists.control.handler();
        self.s   =  namelists.system.handler();
        self.e   = namelists.electrons.handler();
        self.ions= namelists.ions.handler();
        self.cell= namelists.cell.handler();
        self.fcp = namelists.fcp.handler();
        #Default initialization of cards
        self.ae    = cards.atomic_species.handler();
        self.kpts  = cards.kpoints.handler();
        self.cp    = cards.cell_params.handler();
        self.constr= cards.constraints.handler();
        #Default initialization of cards

    def set_structure(self, structure  ):
        """
        Set the atomic position and lattice vectors

        Args:
            structure ( Atoms) : A tuple containing (cell, positions, numbers)n array containing the lattice vectors in arbitrary units, where cell[0] the first vector.

        """
        symbols = structure.get_chemical_symbols();
        unique_symbols = [];
        for s in symbols:
            if s not in unique_symbols:
                unique_symbols.append(s);

        self.s.options["nat"]  = len(symbols);
        self.s.options["ntype"]= len(unique_symbols);

#        self.s.set_atoms(structure);
#        scal_pos= structure.get_scaled_positions();
#        symbols = structure.get_chemical_symbols();
#        xyz = [ (s,x,y,z) for s,(x,y,z) in zip(symbols,scal_pos)]
#        self.s.set_atomic_positions(xyz);
        return True;

    def set_calculation(self, calc  ):
        """
        Defines the calculation to be performed. This flag will initialize different options required for the calculation
        using the default values

        Args:
            calc ( string) : 
        """
        self.c.options["calculation"]= calc;
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
        text = self.c.text();
        text += "\n"+self.s.text();
        text += "\n"+self.e.text();

        calc = self.c.options["calculation"];
        if calc =="relax" or calc =="md" or calc =="vc-relax"or calc =="vc-md":
            text += "\n"+self.ions.text();

        if self.c.options["lfcp"]:
            text += "\n"+self.fcp.text();



        print(text)
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



def generate_from_xyz(xyz, cell, magnetic=False, pbc=(True, True,True) ):
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
    structure    = Atoms( [ Atom(*a) for a in xyz ], cell=cell, pbc=pbc ) ;
    #Initialize the QESuite handler
    qes_h = handler();
    qes_h.set_structure( structure = structure );

    return qes_h;



#c = control.handler();
#print( c.options.keys() )
#q = qes.handler();

#print( q.c.options.keys() )

#import json
#fname ="SSSP_1.1.2_PBEsol_precision.json";
#with open(fname) as f:
#    pseudo_potentials = json.loads(f.read());

