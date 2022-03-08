#import json
#import numpy as np
#import qe_suite.symmetries as symmetries
#import qe_suite.io as qe_io

from . import _builder, namelists, cards 

class QEInput():
    """A generator of input files for Quantum Espresso.

    Parameters:
    ------------
    name:
        A  name to identify the system. 
    structure:
        An instance of a :py:class:`~qe_suite.Structure`
    """

    def __init__(self, name="qe_suite", structure = None, calculation = None, electronic_state = None):
        print("I am here importing a module")

        self.electronic_state = None;
        self.structure = structure;
        self.namelist_h= namelists.handler("scf","structure", "electrons");
        self.cards_h   = cards.handler("scf","structure", "electrons");

    def set_state(self,two_dimensional = False, insultator=False, magnetic=False ):
        pass

    def set_structure(self, structure ):
        pass

    def set_calculation(self, calc  ):
        pass

    def write(self, ofname="out" ):
        pass


class QEInput_():
    """A generator of input files for |QuantumEspreso|.
    
    
    """

        

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


