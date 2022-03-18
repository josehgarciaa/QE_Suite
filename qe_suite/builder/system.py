import numpy as np
import json
import qe_suite.constants as qesc
from ..cards import  atomic_species, atomic_positions, cell_parameters
from ..namelists import system
class System:

    def __init__(self, name="", structure = None):
        self.k_points         = None;
        self.structure        = None;
        self.atomic_species   = None;
        self.cell_parameters  = None;
        self.atomic_positions = None;
        self.periodicity      = (True,True,True);
        self.system = system.System();
        self.set_energy_cutoff_wfc(ecutwfc=60, dual=4);

    def get_system(self):
            return self.system;

    def set_structure(self, structure, use_symmetries=False):

        if use_symmetries:
            ibrav = structure.get_symmetry_informations()["ibrav"];
            self.system.set(ibrav=ibrav )
            if ibrav != 0:
                A, B, C, cosAB, cosAC, cosBC  = structure.get_crystallographic_constants();
                self.system.set(A=A, B=B, C=C, cosAB=cosAB, cosAC=cosAC, cosBC=cosBC )

        self.set_atomic_positions(*structure.get_atomic_positions() )
        self.set_cell_parameters(*structure.get_cell_parameters() );

        return self;

    def set_atomic_species(self, species, use_SSSP=None):    

        if use_SSSP is not None:
            try:
                f= open(use_SSSP);
            except:
                raise FileNotFoundError;
            psinfo = json.loads(f.read());
            f.close();
            
            max_cutoff_wfc = np.max( [ psinfo[specie]["cutoff_wfc"] for specie in species] )
            max_cutoff_rho = np.max( [ psinfo[specie]["cutoff_rho"] for specie in species] )
            self.set_energy_cutoff_wfc(max_cutoff_wfc);
            self.set_energy_cutoff_rho(max_cutoff_rho);

            species = { specie:(qesc.atomic_masses[specie],psinfo[specie]["filename"]) for specie in species }         

        if self.atomic_species is None:
            self.atomic_species = atomic_species.AtomicSpecies();

        self.atomic_species.set("", species);
        self.system.set( ntyp=len(species) )
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

    def get_ibrav(self):
        return self.system.ibrav;

    def set_cell_parameters(self, option, value):
        if self.cell_parameters is None and self.get_ibrav==0:
            self.cell_parameters = cell_parameters.CellParameters()
            self.cell_parameters.set(option,value);
        return self;

    def get_cell_parameters(self):
        return self.cell_parameters;

    def set_numbands(self, nbnd ):
        self.system.set_numbands(nbnd);
        return self;

    def valid(self):
        return True;

    def set_energy_cutoff_wfc(self, ecutwfc, dual=4):
        self.system.set(ecutwfc=ecutwfc);
        if self.system.ecutrho is None:
            self.system.set(ecutrho=dual*ecutwfc);
        return self;

    def set_energy_cutoff_rho(self, ecutrho, dual=4):
        self.system.set(ecutrho=ecutrho);
        return self;

    def get_energy_cutoff_wfc(self):
        return self.system.ecutwfc;

    def get_energy_cutoff_rho(self):
        return self.system.ecutrho;
