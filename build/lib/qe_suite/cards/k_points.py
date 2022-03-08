from .card import Card

class KPoints(Card):
    """The definition of the atomic species used in thecsimulation. 

      Attributes:
      ------------
        This class implement attributes that follow the naming convenction and with
        the same description as given in the Quantum Espresso Card: `\&ATOMIC_POSITIONS <https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm1311`_
        Unless otherwise specifed, the attributes are to `None` and won't be included in the input_file which 
        will force Quantum Espresso to use its default values.

        options can be gamma
        
        When options is automatic, only a tuple of six parameter is needed (nk1,nk2,nk3,sk1,sk2 ,sk3 )


        When used within the python API, the atomic position should be provided as a list of tuples:
         [ (kx_0,ky_0,kz_0,wk_0), (kx_1,ky_1,kz_1,wk_1), ...,(kx_n,ky_n,kz_n,wk_n)] where X_i the label that identifies the position of
         the ith atom that defines the crystalline structure and x_i,y_i,z_i its cartesian position. n here should be equal to number of atoms
         nat 

      Parameters:
      ------------
      option:
        A string defining the units of the atomic position: (tpiba , crystal , tpiba_b , crystal_b , tpiba_c , crystal_c)
      structure:
        An instance of a :py:class:`~qe_suite.builder.Structure`

      Example
      ------------
      >>> # qe_input.cards.atomic_positions.option	= "bohr";
    """

    def __init__(self):
        self.set_option("automatic");
        self.set_value( value=None, type=list );
        self.set_name( "K_POINTS" );

    def set_kpoints(self, kpoints, kptype="automatic", pbc=(True, True, True)):
        s = self
        s.kptype = kptype
        if kptype == "automatic":
            s.kpoints, s.shifts = kpoints
            s.kpoints = [k if per else 1 for per, k in zip(pbc, s.kpoints)]
            s.options["K_POINTS"] = s.kptype + \
                "\n{} {} {} {} {} {} ".format(*s.kpoints, *s.shifts)
        if kptype == "crystal_b":
            s.kpoints = []
            for k, v in kpoints.items():
                s.kpoints += list(v)
            # Save kpoints as text
            s.options["K_POINTS"] = s.kptype+"\n"+str(len(s.kpoints))
            for kp in s.kpoints:
                s.options["K_POINTS"] += "\n {} {} {} 1 ".format(*list(kp))

    def get_kpoints(self):
        s = self
        if s.kptype == "automatic":
            return (s.kpoints, s.shifts)


    def __str__(self):
        out= self.header();
        if self.option == "gamma":
            return out+"\n";

        if self.option == "automatic":
            return out+" {} {} {} {} {} {}\n".format(*self.get_value());

        #else tpiba  | crystal | tpiba_b | crystal_b | tpiba_c | crystal_c
        nks = len(self.get_value());
        out+= str(nks);
        for v in self.get_value():
            out+= " {} {} {} {}\n".format(*v);
        return out