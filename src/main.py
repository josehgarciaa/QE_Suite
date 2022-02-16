import qesuite as qes

import qe_io 


syst = "Pd2Se4"
xyz  = qe_io.load_xyz(syst+".xyz");
cell = qe_io.load_cell(syst+".uc");

qes_handler =  qes.generate_from_xyz(xyz=xyz, cell=cell);

qes_handler.set_calculation("scf");
qes_handler.use_symmetries();
qes_handler.use_SSSP(functional="PBE", target="precision", path="./SSSP");

qes_handler.write_input_file("input.cfs");



