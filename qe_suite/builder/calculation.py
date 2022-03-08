from .structure import Structure
from ..namelists import control
from ..cards import   k_points
class Calculation():
    
    def __init__(self) -> None:
        self.control          = control.Control();
        self.k_points         = None;

    def set_pseudopot_dir(self, psdir):
        self.control.pseudo_dir = psdir
        return self;

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
        self.control.outdir = './qe_suite/'
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


    def valid(self):
        return True;


from pathlib import Path
import errno
class NSCF(Calculation):

    def __init__(self, scf = None , startingpot='file' ) -> None:
        super().__init__()
        self.control = scf.get_control();
        self.control.calculation = 'nscf'
        self.k_points         = None;



        prefix = self.control.prefix;
        outdir = self.control.outdir;
        xml_0 = outdir+prefix+".save/data-file-schema.xml";
        xml_1 = outdir+prefix+".xml";
        if ( not Path(xml_0).is_file()) and ( not Path(xml_1).is_file()):
            print("A NSCF calculation requires a valid xml file either at",xml_0, "or", xml_1)
            raise FileNotFoundError

        def valid(self) :
            return True;

#types = [ 'bands', 'relax','md', 'vc-relax', 'vc-md'];
from pathlib import Path
import errno

class Bands(Calculation):

    def __init__(self, scf = None , startingpot='file' ) -> None:
        super().__init__()
        self.control = scf.get_control();
        self.control.calculation = "bands";

        prefix = self.control.prefix;
        outdir = self.control.outdir;
        xml_0 = outdir+prefix+".save/data-file-schema.xml";
        xml_1 = outdir+prefix+".xml";
        if ( not Path(xml_0).is_file()) and ( not Path(xml_1).is_file()):
            print("A Bands calculation requires a valid xml file either at",xml_0, "or", xml_1)
            raise FileNotFoundError


    def set_band_path(self, bandpath):
        kpoints = [];
        for kps in bandpath.values():
            kpoints+= list(kps)
        self.k_points = k_points.KPoints()
        self.k_points.set("crystal_b", kpoints)
        return self

    def valid(self) :
        return True;