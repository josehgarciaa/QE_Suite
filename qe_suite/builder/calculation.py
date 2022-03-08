from inspect import signature
from xmlrpc.client import Boolean

import qe_suite
from .structure import Structure
from ..namelists import control
from ..cards import   k_points
class Calculation():
    
    def __init__(self) -> None:
        self.control          = control.Control();
        self.k_points         = None;

    def set_control(self, control):
        return self.control;

    def get_control(self):
        return self.control;

    def set_k_points(self, option, value):
        if self.k_points is None:
            self.k_points = k_points.KPoints()
        self.k_points.set(option,value);
        return self;

    def get_k_points(self):
        return self.k_points;

        
class SCF(Calculation):

    def __init__(self, target = "default", pseudo_dir='.' ) -> None:
        super().__init__()
        self.control.calculation = 'scf'
        self.control.outdir = './out/'
        self.control.prefix = 'qe_suite'
        self.control.pseudo_dir = pseudo_dir
        self.control.tprnfor = True
        self.control.tstress = True

        if target == "precission":
            self.control.etot_conv_thr = 1e-6;
            self.control.forc_conv_thr = 1e-3;


        if target == "efficiency":
            self.control.etot_conv_thr = 1e-5;
            self.control.forc_conv_thr = 1e-4;

        if target == "default":
            self.control.etot_conv_thr = 1e-4;
            self.control.forc_conv_thr = 1e-3;


    def valid(self) -> Boolean:
        return True;



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
