import qe_suite.qe_suite as qes
import qe_suite.io as qe_io


import numpy as np
from  numpy.linalg import norm, det
from numpy import dot
import sympy as sp


from copy import copy

#Reference: https://www.johndcook.com/blog/2010/10/20/best-rational-approximation/


a_syst = "C2"
a_xyz  = qe_io.load_xyz(a_syst+".xyz");
a_cell = qe_io.load_cell(a_syst+".uc");

b_syst = "Pd2Se4"
b_xyz  = qe_io.load_xyz(b_syst+".xyz");
b_cell = qe_io.load_cell(b_syst+".uc");

print( get_vdw_cell( a_cell, b_cell, max_strain=0.2, strain_cell="a", max_size=10 ) )

