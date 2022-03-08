import numpy as np
from ..cards import  atomic_species, atomic_positions, cell_parameters
from ..namelists import system
class System:

    def __init__(self, name="", structure = None):
        self.k_points         = None;
        self.structure        = None;
        self.atomic_species   = None;
        self.cell_parameters  = None;
        self.atomic_positions = None;

        self.system = system.System();

    def get_system(self):
            return self.system;

    def set_structure(self, structure):

        #MOVE TO VALID
#        symbols = list(set( structure.get_chemical_symbols() ));
#        symbols = ["K", "C"]
#        defined_pseudos = self.get_atomic_species().get_value().keys();

#        found =  np.all([ np.any([ s==p for p in defined_pseudos])  for s in symbols ]);
#        if not found:
#            raise 

        #    print("Falta un pseudo",notFound)
#        found = ( ( p==p for p in defined_pseudos ) for s in symbols );
#        print(found)

        self.set_atomic_positions(*structure.get_atomic_positions() )
        self.set_cell_parameters(*structure.get_cell_parameters() );

        return self;


    def set_atomic_species(self, value):       
        if self.atomic_species is None:
            self.atomic_species = atomic_species.AtomicSpecies();

        self.atomic_species.set("", value);
        self.system.set( ntyp=len(value) )
        return self;

    def get_atomic_species(self):
            return self.atomic_species

    def set_atomic_positions(self, option, value):
        if self.atomic_positions is None:
            self.atomic_positions = atomic_positions.AtomicPositions()
        self.atomic_positions.set(option,value);
        self.system.set( nat=len(value) )
        return self;

    def get_atomic_positions(self) :
        return self.atomic_positions;

    def set_cell_parameters(self, option, value):
        if self.cell_parameters is None:
            self.cell_parameters = cell_parameters.CellParameters()
        self.cell_parameters.set(option,value);

        lat_const = np.min(np.linalg.norm(value, axis=1));
        if self.system.ibrav!=0:    
           self.system.set( A=lat_const )
        return self;

    def get_cell_parameters(self):
        return self.cell_parameters;


    def valid(self):
        return True;


    def set_atomic_species_from(self,library = "SSSP"):
        print(" I will generate the pseudopotentials alone")
        return self;


    def use_structure_as_symmetries(self):
        print("This flag will determine the symmetries of the structure")

    def refiene_symmetries(self,):
        print("Get the propert symmetries")


    def symmetry_based_bandpath(self,):
        print("Get the band path")


