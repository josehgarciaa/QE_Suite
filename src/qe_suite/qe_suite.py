
import json
import numpy as np
import qe_suite.symmetries as symmetries
import qe_suite.io as qe_io
from ase import Atom, Atoms
from qe_suite.namelists import control, system, electrons, ions, cell, fcp
from qe_suite.cards import atomic_positions, atomic_species, kpoints, cell_params, constraints 

class handler():
    def __init__(self):
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

    def set_structure(self, structure, two_dimensional = False  ):
        """
        Set the atomic position and lattice vectors

        Args:
            structure ( Atoms) : A tuple containing (cell, positions, numbers)n array containing the lattice vectors in arbitrary units, where cell[0] the first vector.

        """
        self.structure = structure;
        symbols = structure.get_chemical_symbols();
        species = [ (s, Atom(s).mass, s+".UPF") for s in set(symbols)]
        self.s.options["nat"]  = len(symbols);
        self.s.options["ntyp"]= len(species);

        if ( self.structure.pbc == (True,True,False) ).all() :
            self.s.options["assume_isolated"]='2D';
            #Define the heigh using the guidelines with QE
            zs     = [ z for (x,y,z) in structure.get_positions() ];
            delta_z= np.max( [24, 7.0 + np.max(zs)- np.min(zs)] ); #Angstrong
            cell   = structure.get_cell(); 
            cell[2]= [0,0,delta_z];
            structure.set_cell(cell);

            #Define the kpoints
            self.kpts.set_kpoints( self.kpts.get_kpoints(), pbc=structure.pbc );

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
        structure, ibrav, spgnum, celldm = symmetries.get_brav_params( self.structure );
        self.s.set_bravais_lattice(ibrav, celldm);
        self.set_structure(structure);
        
        if self.c.options["calculation"]=="bands":
            kpoints = symmetries.get_band_path( self.structure );
            self.kpts.set_kpoints( kptype= "crystal_b", kpoints= kpoints);
            

        return 0;

    def write_input_file(self, ofname="out" ):
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

    if two_dimensional:
        structure.set_pbc( (True, True, False) );

    #Initialize the QESuite handler
    qes_h = handler();
    qes_h.set_structure( structure = structure, two_dimensional=two_dimensional );

    return qes_h;
