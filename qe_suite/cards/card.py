


class Card:

    def __init__(self):
        self.option= "";
        self.value = None;
        self.name  = None;
        self.type  = None;

    def set_name(self,name):
        self.name = name;

    def set(self, option, value):
        self.option = option;
        self.value  = value;

    def set_value(self, value=None, type=None ):
        if value is not None:
            self.value = value;

        if type is not None:
            self.type = type;

    def set_option(self,option):
        self.option = option;

    def get_name(self):
        return self.name;

    def get_option(self):
        return self.option;

    def get_value(self):
        return self.value;

    def header(self):
        return self.get_name()+" " + self.get_option()+"\n";
