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
        self.cell   = None;
        self.atomic_positions = None;
        self.spin_state = None;
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
    
    def get_spin_state(self):
        if self.spin_state is None:
            root= self.tree.getroot()
            spin_state = root.find("input/spin").iter();
            self.spin_state = { x.tag: False if x.text=="false" else True for x in spin_state } ;
        return self.spin_state;