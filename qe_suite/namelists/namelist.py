import qe_suite.io as qe_io

class Namelist():

    
    def __init__(self):
        self.name = ""
        self.parameters = {};

    def set(self, **kwargs ):
        for k, v in kwargs.items():
            if (k not in self.parameters) or  (v is  None):
                raise ValueError("key:",k," does not exists or ",  v, "is None")
            self.parameters.update({k:v})
        self.__dict__.update(self.parameters);

    def set_name(self,name):
        self.name = name;

    def set_parameters(self, parameters):
        self.parameters = parameters;

    def update_parameters(self, parameters):
        for key in self.parameters.keys():
            if key in parameters:
                self.parameters.update({key:parameters[key]})

    def get_parameters(self):
        return self.parameters;

    def __str__(self):
        #Update parameters from the initialize parameters in the class 
        self.update_parameters(self.__dict__);
        #Check is there is at leas a  None, which means the namelist is not set
        if not all(v is None for v in self.get_parameters().values()):
            out = self.name+"\n";
            for k, v in self.get_parameters().items():
                if v is not None and k != "name":
                    out += qe_io.key_format(k)+"="+qe_io.format(v)+"\n"
            out += "/\n"
            return out
        #If all parameters are none it means the namelist is not present
        return ""
