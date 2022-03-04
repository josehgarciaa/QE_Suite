import qe_suite.io as qe_io

class Namelist():

    
    def __init__(self):
        self.name = ""

    def set_namelist_name(self,name):
        self.name = name;

    def __str__(self):
        out = self.name+"\n";
        for k, v in self.__dict__.items():
            if v is not None and k != "name":
                out += k+"="+qe_io.format(v)+"\n"
        out += "/\n"
        return out
