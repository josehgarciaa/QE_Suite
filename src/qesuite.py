import qe_io
import json
import symmetries
import namelists.control
import namelists.system
import namelists.electrons
import namelists.ions
import namelists.cell
import namelists.fcp
import cards.atomic_positions
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
        self.s   = namelists.system.handler();
        self.e   = namelists.electrons.handler();
        self.ions= namelists.ions.handler();
        self.cell= namelists.cell.handler();
        self.fcp = namelists.fcp.handler();
        #Default initialization of cards
        self.ap    = cards.atomic_positions.handler();
        self.ae    = cards.atomic_species.handler();
        self.kpts  = cards.kpoints.handler();
        self.cp    = cards.cell_params.handler();
        self.constr= cards.constraints.handler();
        #Default initialization of cards
        self.structure = None;

    def set_structure(self, structure  ):
        """
        Set the atomic position and lattice vectors

        Args:
            structure ( Atoms) : A tuple containing (cell, positions, numbers)n array containing the lattice vectors in arbitrary units, where cell[0] the first vector.

        """
        symbols = structure.get_chemical_symbols();
        species = [ (s, Atom(s).mass, s+".UPF") for s in set(symbols)]
        self.s.options["nat"]  = len(symbols);
        self.s.options["ntype"]= len(species);

        self.ae.set_atomic_species(species);
        self.ap.set_atomic_positions(structure);
        self.structure = structure;
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
        system, ibrav, spgnum, celldm = symmetries.get_brav_params( self.structure );
        self.s.set_bravais_lattice(ibrav, celldm);
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

                
        text += self.ae.text();
        text += self.ap.text();
        text += self.kpts.text();

        print(text)
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

        with open(path+"/versions.yaml") as f:
            version = f.read().split(":")[-1].replace("\'","").replace(" ","").replace("\n","");

        dbfname = path+"/SSSP_"+version+"_"+functional+"_"+target+".json";
        with open(dbfname) as f:
            sp_info = json.loads(f.read());

        ae = self.ae;
        species = [ (s,m ,path+"/"+sp_info[s]["filename"]) for s,m,sp in ae.get_atomic_species() ]
        ae.set_atomic_species(species);

        if use_cutoff:
            cutoffs=[0,0];
            for s,m,sp in ae.get_atomic_species():
                sp = sp_info[s];
                keys = ( "cutoff_wfc", "cutoff_rho");
                cutoffs = [v if sp[k]<v else sp[k] for k,v in zip( keys,cutoffs )];
                self.s.set_cutoff(*cutoffs);

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
    structure    = Atoms( [ Atom(*a) for a in xyz ] ) ;
    structure.set_cell(cell);
    structure.set_pbc(pbc);
    #Initialize the QESuite handler
    qes_h = handler();
    qes_h.set_structure( structure = structure );

    return qes_h;
