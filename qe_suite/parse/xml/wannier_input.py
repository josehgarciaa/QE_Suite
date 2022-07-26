from logging import exception
import qe_suite.format as f
import qe_suite.constants as c
import numpy as np
import xml.etree.ElementTree as ET

def strvec2list(x):
    return list(map(float,x.split(" ")));    

class WannierInput():

    def __init__(self, xml) -> None:
        self.xml    = None;
        self.tree   = None;
        self.outdir = None;
        self.prefix = None;
        self.cell   = None;
        self.atomic_positions = None;
        self.spin_state = None;
        self.kpoints = None;
        self.num_bands = None;

        self.set_xml_tree(xml=xml);

    def set_xml_tree(self,xml=None):
        if xml is not None:
            self.xml=xml;
        try:
            tree = ET.parse(self.xml);
        except:
            raise FileNotFoundError
        self.tree=tree;
        return self;

    def get_prefix(self):
        if self.prefix is None:
            root = self.tree.getroot()
            self.prefix = root.find("input/control_variables/prefix").text ;
        return self.prefix;

    def get_outdir(self):
        if self.outdir is None:
            root = self.tree.getroot()
            self.outdir = root.find("input/control_variables/outdir").text ;
        return self.outdir;


    def get_cell(self):
        if self.cell is None:
            root= self.tree.getroot()
            cell  = root.find("output/atomic_structure/cell");
            self.cell  = [ strvec2list( cell.find("a"+str(i)).text ) for i in (1,2,3) ];           
        return self.cell;

    def get_atomic_positions(self):
        if self.atomic_positions is None:
            root= self.tree.getroot()
            atom_pos = root.find("output/atomic_structure/atomic_positions");
            atom_pos = [ (x.attrib["name"],strvec2list(x.text) ) for x in atom_pos.iter("atom") ];
            self.atomic_positions = atom_pos;
        return self.atomic_positions;
    
    def get_fractional_atomic_positions(self):
        atm_pos= self.get_atomic_positions();
        to_frac= np.linalg.inv( self.get_cell() );
        return [ ( k,to_frac.dot(v) ) for k,v in atm_pos ];

    def get_kpoints(self):
        if self.kpoints is None:
            root = self.tree.getroot()
            kpoints = root.find("input/k_points_IBZ");
            kpoints = [strvec2list(x.text) for x in kpoints.iter("k_point")]
            self.kpoints = kpoints;
        return self.kpoints;

    def get_num_bands(self):
        if self.num_bands is None:
            root = self.tree.getroot()
            self.num_bands = int( root.find("input/bands/nbnd").text );
        return self.num_bands;

    def get_spin_state(self):
        if self.spin_state is None:
            root= self.tree.getroot()
            spin_state = root.find("input/spin").iter();
            self.spin_state = { x.tag: False if x.text=="false" else True for x in spin_state } ;
        return self.spin_state;