from .structure import Structure
from ..namelists import control
from ..cards import atomic_species, atomic_positions, cell_parameters, k_points
class Calculation():
    
    def __init__(self) -> None:
        self.control          = control.Control();
        self.atomic_species   = None;
        self.atomic_positions = None;
        self.k_points         = None;
        self.cell_parameters  = None;
    def get_control(self):
        return self.control;

    def get_atomic_positions(self):
        return self.atomic_positions;

    def set_atomic_species(self, value):
        if self.atomic_species is None:
            self.atomic_species = atomic_species.AtomicSpecies();
        self.atomic_species.set("",value);
        
    def set_atomic_positions(self, option, value):
        if self.atomic_positions is None:
            self.atomic_positions = atomic_positions.AtomicPositions()
        self.atomic_positions.set(option,value);
        
    def set_k_points(self, option, value):
        if self.k_points is None:
            self.k_points = k_points.KPoints()
        self.k_points.set(option,value);
        
    def set_cell_parameters(self, option, value):
        if self.cell_parameters is None:
            self.cell_parameters = cell_parameters.CellParameters()
        self.cell_parameters.set(option,value);

    
        
class SCF(Calculation):

    def __init__(self, atomic_positions=None, k_point = None, target = "default", pseudo_dir='.' ) -> None:
        super().__init__()
        self.control.calculation = 'scf'
        self.control.outdir = './out/'
        self.control.prefix = 'qe_suite'
        self.control.pseudo_dir = pseudo_dir
        self.control.tprnfor = True
        self.control.tstress = True

        if target == "precission":
            self.control.etot_conv_thr = 1e-6;
            self.control.etot_conv_thr = 1e-3;


        if target == "efficiency":
            self.control.etot_conv_thr = 1e-5;
            self.control.etot_conv_thr = 1e-4;

        if target == "default":
            self.control.etot_conv_thr = 1e-4;
            self.control.etot_conv_thr = 1e-3;

        self.atomic_positions = None;        
        



class NSCF(Calculation):

    def __init__(self, scf = None , startingpot='file' ) -> None:
        super().__init__()
        self.control = scf.get_control();
        self.control.calculation = 'nscf'
        self.control.startingpot = startingpot
        print("Check if startingpot/charge-density.xml exists")



#types = [ 'bands', 'relax','md', 'vc-relax', 'vc-md'];


class NSCF(Calculation):

    def __init__(self, scf = None , startingpot='file' ) -> None:
        super().__init__()
        self.control = scf.get_control();
        self.control.calculation = 'nscf'
        self.control.startingpot = startingpot
        print("Check if startingpot/charge-density.xml exists")
