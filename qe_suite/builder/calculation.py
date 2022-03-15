from .structure import Structure
from ..namelists import control,ions,cell
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



class Relaxation(Calculation):
        
    def __init__(self ) -> None:
        super().__init__()
#        self.system  = system;
        self.control.set(calculation="vc-relax");
        self.ions = ions.Ions();
        self.cell = cell.Cell();

        self.set_cell_pressure_threshold( threshold=0.1);
        self.set_ions_dynamics( dynamic= 'bfgs');
        self.set_cell_dynamics( dynamic= 'bfgs');
        self.set_cell_do_free( freedom="all");


    def set_cell_pressure_threshold(self, threshold=0.1):
        self.cell.set(press_conv_thr = threshold);
        return self;


    def set_ions_dynamics(self, dynamic):
        if dynamic not in ('bfgs','damp'):
            print("not proper dynamics=",dynamic, "in ions_dynamics function")
            raise ValueError;
        self.ions.set(ion_dynamics = dynamic);
        return self;

    def set_cell_dynamics(self, dynamic):
        if dynamic not in ('bfgs','damp'):
            print("not proper dynamics=",dynamic, "in cell_dynamics function")
            raise ValueError;
        self.cell.set(cell_dynamics = dynamic);
        return self;


    def set_cell_do_free(self, freedom):
        allowed_freedoms= ('all','ibrav','x','y','z','xy','xz','yz','xyz','shape','volume','2Dxy','2Dshape','epitaxial_ab','epitaxial_ac', 'epitaxial_bc')
        if freedom not in allowed_freedoms:
            print("not proper freedom=",freedom, "in set_cell_do_free function")
            raise ValueError;
        self.cell.set(cell_dofree = freedom);
        return self;

    def get_ions(self):
        return self.ions;

    def get_cell(self):
        return self.cell;


    def valid(self) :
        return True;





