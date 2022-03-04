
import math as m

import qe_suite.builder as builder;


#import qe_suite.qe_suite as qes
#import qe_suite.io as qe_io


system = "graphene"

from qe_suite.namelists.control import Control

#Define Graphene structure
a= 2.467  # lattice constant in nm which is the default unit. 
cell = [ [a,0,0], [-a/2, a*m.sqrt(3)/2,0], [0,0,5*a] ];
fractional_positions = [ [0,0,0.0], [2/3,1/3,0.0] ]
atomic_symbols = ['C', 'C']
structure = builder.Structure(cell, fractional_positions, atomic_symbols);
#Evaluate the symmetries of the structure and symmetrize it
print("The space group is: ", structure.hm_symbol() )


#Create a configuration file based on graphene's structure
qe_input = builder.QEInput(name=system, structure= structure);
qe_input.namelists.control.calculation = "SCF"
qe_input.write("out");



#xyz  = qe_io.load_xyz(syst+".xyz");
#cell = qe_io.load_cell(syst+".uc");

#qes_handler =  qes.generate_from_xyz(xyz=xyz, cell=cell, two_dimensional=True);

#qes_handler.set_calculation("scf");
#qes_handler.use_symmetries();
#qes_handler.use_SSSP(functional="PBEsol", target="precision", path="../SSSP");

#qes_handler.write_input_file("QEsuite.scf");



#from qe_suite import format as f
#todict = "ibrav= None; celldm= None; A= None; B=None; C=None; cosAB=None; cosAC=None; cosBC = None; nat= None; ntyp=None; nbnd=None; ecutwfc=None; ecutrho= None; ecutfock = None; tot_charge= None; starting_charge= None; tot_magnetization= None; starting_magnetization= None; nr1 = None; nr2= None; nr3 = None;; nr1s= None; nr2s= None; nr3s = None;; nosym =None; nosym_evc = None; noinv =None; no_t_rev =None; force_symmorphic = None; use_all_frac = None; occupations =None; one_atom_occupations = None; starting_spin_angle = None; degauss = None; smearing = None; nspins = None; noncolin=None; ecfixed = None; qcutz= None; q2sigma=None; input_dft = None; ace = None; exx_fraction =None; screening_parameter = None; exxdiv_treatment =None; x_gamma_extrapolation = None; ecutvcut =None; nqx1= None; nqx2= None; nqx3 = None; localization_thr = None; lda_plus_u = None; lda_plus_u_kind = None; Hubbard_U = None; Hubbard_J0= None; Hubbard_V = None; Hubbard_alpha = None; Hubbard_beta= None; Hubbard_J = None; starting_ns_eigenvalue = None; U_projection_type = None; Hubbard_parameters = None; dmft = None; dmft_prefix= None; ensemble_energies = None; edir =None; emaxpos = None; eopreg=None; eamp = None; zgate = None; relaxz = None; block=None; block_1=None; block_2 = None; block_height=None; angle1= None; angle2 = None; lforcet = None; constrained_magnetization = None; fixed_magnetization = None; report = None; Lambda = None; lspinorb = None; assume_isolated = None; esm_bc = None; esm_w = None; esm_efield=None; esm_nfit = None; lgcscf	= None; gcscf_mu = None; gcscf_conv_thr = None; gcscf_beta = None; vdw_corr = None; london_s6= None; london_c6 = None; london_rvdw = None; london_rcut = None; dftd3_version = None; dftd3_threebody = None; ts_vdw_econv_thr = None; ts_vdw_isolated = None; xdm_a1 = None; xdm_a2 = None; space_group = None; uniqueb = None; origin_choice = None; rhombohedral = None;"
#for kv_pair in todict.replace("\t"," ").split(";"):
#    kv_pair = kv_pair.replace(" ", "");
#    if len(kv_pair) != 0:
#        k,v = kv_pair.split("=");
#        print("\""+k+"\":", "None," )
