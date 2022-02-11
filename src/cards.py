cards = ["ATOMIC_SPECIES","ATOMIC_POSITIONS","K_POINTS","CELL_PARAMETERS",
         "OCCUPATIONS","CONSTRAINTS","ATOMIC_VELOCITIES","ATOMIC_FORCES"];

class specie:
    symbol   = "X";    
    mass     = "Mass_X";
    pseudopot= "PseudoPot_X";
    positions=  [ [], [] , []]
    #ATOMIC_POSITIONS { alat | bohr | angstrom | crystal | crystal_sg }
    #   X 0.0  0.0  0.0  {if_pos(1) if_pos(2) if_pos(3)}

species = [ specie("X","Mass_X","PseudoPot_X"), 
            specie("Y","Mass_Y","PseudoPot_Z") ]; #list of species. 

class kpoints:
    """
    K_POINTS { tpiba | automatic | crystal | gamma | tpiba_b | crystal_b | tpiba_c | crystal_c }
    if (gamma)
    nothing to read
    if (automatic)
    nk1, nk2, nk3, k1, k2, k3
    if (not automatic)
    nks
    xk_x, xk_y, xk_z,  wk
    if (tpipa_b or crystal_b in a 'bands' calculation) see Doc/brillouin_zones.pdf
    """

class cell:
    """
    [ CELL_PARAMETERS { alat | bohr | angstrom }
    v1(1) v1(2) v1(3)
    v2(1) v2(2) v2(3)
    v3(1) v3(2) v3(3) ]
    """
