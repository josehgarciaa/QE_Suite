class handler:
    """
    The species class storages and handle everything relate to different atomic species in the simulation.

    Attributes:

    """
    def __init__(self):
        self.kptype  = "automatic";
        self.kpoints = [1,1,1];
        self.shifts  = [1,1,1];
        self.options = dict();
        self.set_kpoints( ((10,10,10),(1,1,1)), kptype = "automatic")

    def set_kpoints(self, kpoints, kptype = "automatic", pbc=(True,True,True)):
        s = self;
        s.kptype = kptype;
        if kptype == "automatic":
            s.kpoints, s.shifts = kpoints
            s.kpoints =[ k if per else 1 for per,k in zip(pbc,s.kpoints) ];
            s.options["K_POINTS"]=s.kptype+"\n{} {} {} {} {} {} ".format(*s.kpoints,*s.shifts);
        if kptype == "crystal_b":
            s.kpoints = [];
            for k,v in kpoints.items():
                s.kpoints+= list(v)
            #Save kpoints as text
            s.options["K_POINTS"]=s.kptype+"\n"+str(len(s.kpoints));
            for kp in s.kpoints:
                s.options["K_POINTS"]+="\n {} {} {} 1 ".format(*list(kp));

    def get_kpoints(self):
        s = self;
        if s.kptype == "automatic":
            return (s.kpoints, s.shifts);

    def text(self):
        key = "\nK_POINTS ";
        return key + self.options["K_POINTS"];

    def print(self):
        print(self.text());




    
