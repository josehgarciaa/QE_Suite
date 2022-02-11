namelists = ["&CONTROL","&SYSTEM","&ELECTRONS","&IONS","&CELL", "&FCP"]


def get_namelist_block(nl, data):
    if nl in data:
        beg = data.find(nl) + len(nl);
        data= data[beg:];
        end = data.find("\n/");
        data= data[:end];
        return data
    return None;
