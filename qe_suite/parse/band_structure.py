from logging import exception
import qe_suite.format as f
import qe_suite.constants as c
import numpy as np
import xml.etree.ElementTree as ET

class BandStructure():

    def __init__(self, xml) -> None:
        self.xml=None;
        self.tree= None;
        self.kpoints= None;
        self.nelec  = None;
        self.ylabel = "E (eV)"
        self.set_xml_tree(xml=xml);
        self.set_kpoints();


    def set_xml_tree(self,xml=None):
        if xml is not None:
            self.xml=xml;
        try:
            tree = ET.parse(self.xml);
        except:
            raise FileNotFoundError
        self.tree=tree;
        return self;


    def xticks_high_symmetry_points(self, hsp ):
        root   = self.tree.getroot()
        rec_lat= root.find("output/basis_set/reciprocal_lattice");
        rec_vec= np.array([ rec_lat.find(b).text.split(" ") for b in ("b1","b2","b3") ], dtype=float);

        labels=  list(hsp.keys());
        kpts  = np.dot( list(hsp.values()),rec_vec);

        hsp_idx = [ np.argmin(np.linalg.norm(self.get_kpoints() -k, axis=1))  for k in kpts]
        xticks = self.get_xaxis()[hsp_idx];

        return xticks, labels;


    def get_num_electrons(self):
        if self.nelec is None:
            root= self.tree.getroot()
            bs  = root.find("output/band_structure");
            self.nelec = int(float(next(bs.iter("nelec")).text));           
        return self.nelec

    def set_kpoints(self):
        root = self.tree.getroot()
        bs   = root.find("output/band_structure");
        kpoints= [(x.find("k_point").text).split(" ") for x in bs.iter("ks_energies")]
        self.kpoints = np.array(kpoints, dtype=float);
        return self;

    def get_kpoints(self):
        return self.kpoints;

    def get_fermi_energy(self):
        root = self.tree.getroot()
        bs   = root.find("output/band_structure");
        EF   = float( root.find("output/band_structure/fermi_energy").text )*c.Ry2eV;
        return EF;

    def get_eigenvalues(self):
        root = self.tree.getroot()
        bs   = root.find("output/band_structure");
        eigvals= [f.remove_empty((eigv.find("eigenvalues").text.replace("\n","")).split(" ")) for eigv in bs.iter("ks_energies")]
        eigvals= np.array(eigvals, dtype=float).T*c.Ry2eV;
        return eigvals;

    #Compute the X-axis of a band structure by summing the norm of the k steps
    def get_xaxis(self):
        kpath = self.kpoints
        return np.cumsum( np.linalg.norm(np.diff(kpath,axis=0,prepend=[kpath[0]]),axis=1) );

    def bands(self, shift_Efermi=False):
        eigvals = self.get_eigenvalues();
        xaxis   = self.get_xaxis();
        if shift_Efermi:
            eigvals = eigvals - self.get_fermi_energy();
            self.ylabel = self.ylabel.replace("(eV)","- EF (eV)")

        return np.array([ (xaxis, list(v)) for v in eigvals ])

    def get_ylabel(self):
        return self.ylabel;
