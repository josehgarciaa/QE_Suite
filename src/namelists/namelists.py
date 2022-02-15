import _namelists.control as control
import _namelists.system  as systems
import _namelists.electrons as electrons
import _namelists.ions as ions
import _namelists.cell as cell
import _namelists.fcp as fcp


namelists = ["&CONTROL","&SYSTEM","&ELECTRONS","&IONS","&CELL", "&FCP"]


def get_namelist_block(nl, data):
    if nl in data:
        beg = data.find(nl) + len(nl);
        data= data[beg:];
        end = data.find("\n/");
        data= data[:end];
        return data
    return None;
