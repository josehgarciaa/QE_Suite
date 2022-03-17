# Package directory
import os
from re import M
import sys
sys.path.insert(0, os.path.abspath('../../'))
sys.path.insert(0, os.path.abspath('.'))

from math import sqrt
import qe_suite as qes
from qe_suite.builder import Structure, System, PWInput
from qe_suite.builder.calculation import SCF, Relaxation

#First we define the structure of the system
hexagonal= Structure( cell = [(2.467, 0.000, 0.0), (-1.234, 2.137, 0.0), (0.000, 0.000, 15.0)],
                      fractional_positions = [ (0.6667, 0.3333, 0.5), (0.3333, 0.6667, 0.5) ],
                      atomic_symbols = ["C", "C"], 
                      symmetrize=True);
hexagonal.set_as_2D();

#That structure is used to define the system to be study
graphene = System();
graphene.set_atomic_species( ["C"], use_SSSP= "../SSSP/SSSP_1.1.2_PBEsol_precision.json" );
graphene.set_structure(hexagonal, use_symmetries=True);

#Then we create a self consistent calculation to determine the density matrix
relax = Relaxation().set_pseudopot_dir("../SSSP").set_cell_do_free('ibrav');
relax.set_k_points("automatic", hexagonal.get_kpoints(type="automatic"))


pw_input = PWInput(calculation = relax, system=graphene  );
inpfile=pw_input.write("qe_suite.relax.inp")
#qes.run_pw(inpfile=pw_input.write("qe_suite.relax.inp"),shell=False );



#ibrav      structure                   celldm(2)-celldm(6)
# 12          Monoclinic P, unique axis c     celldm(2)=b/a
#                                             celldm(3)=c/a,
#                                             celldm(4)=cos(ab)
#      v1=(a,0,0), v2=(b*cos(gamma),b*sin(gamma),0),  v3 = (0,0,c)
#      where gamma is the angle between axis a and b.
#-12          Monoclinic P, unique axis b     celldm(2)=b/a
#                                             celldm(3)=c/a,
#                                             celldm(5)=cos(ac)
#      v1 = (a,0,0), v2 = (0,b,0), v3 = (c*cos(beta),0,c*sin(beta))
#      where beta is the angle between axis a and c

# 13          Monoclinic base-centered        celldm(2)=b/a
#             (unique axis c)                 celldm(3)=c/a,
#                                             celldm(4)=cos(gamma)
#      v1 = (  a/2,         0,          -c/2),
#      v2 = (b*cos(gamma), b*sin(gamma), 0  ),
#      v3 = (  a/2,         0,           c/2),
#      where gamma=angle between axis a and b projected on xy plane

#-13          Monoclinic base-centered        celldm(2)=b/a
#             (unique axis b)                 celldm(3)=c/a,
#                                             celldm(5)=cos(beta)
#      v1 = (  a/2,       b/2,             0),
#      v2 = ( -a/2,       b/2,             0),
#      v3 = (c*cos(beta),   0,   c*sin(beta)),
#      where beta=angle between axis a and c projected on xz plane
# IMPORTANT NOTICE: until QE v.6.4.1, axis for ibrav=-13 had a
# different definition: v1(old) =-v2(now), v2(old) = v1(now)

 #14          Triclinic                       celldm(2)= b/a,
 #                                            celldm(3)= c/a,
 #                                            celldm(4)= cos(bc),
 #                                            celldm(5)= cos(ac),
 #                                            celldm(6)= cos(ab)
 #     v1 = (a, 0, 0),
 #     v2 = (b*cos(gamma), b*sin(gamma), 0)
 #     v3 = (c*cos(beta),  c*(cos(alpha)-cos(beta)cos(gamma))/sin(gamma),
 #          c*sqrt( 1 + 2*cos(alpha)cos(beta)cos(gamma)
 #                    - cos(alpha)^2-cos(beta)^2-cos(gamma)^2 )/sin(gamma) )
 #     where alpha is the angle between axis b and c
 #            beta is the angle between axis a and c
 #           gamma is the angle between axis a and b

