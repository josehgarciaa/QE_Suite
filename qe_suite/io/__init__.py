import qe_suite.format as f
import .wannier90 as wannier90

def key_format(x):
    
    if not isinstance(x, str):
        raise ValueError("Error value in key_format")

    if "celldm" in x:
        print(x)
        celldm, i = x.split("_");
        return celldm +"("+i+")"; 

    return x;

def format(x):
    
    if isinstance(x, str):
        return "\'"+x+"\'";

    if isinstance(x, bool):
        return ".TRUE." if x else ".FALSE.";

    try: 
        x = str(x)
    except ValueError:
        print(x," is not a valid string") 

    if "celldm" in x:
        print(x)
        celldm, i = x.split();
        return celldm +"("+i+")"; 

    return str(x);

def wannnier90():
    def key_format(x):
        
        if not isinstance(x, str):
            raise ValueError("Error value in key_format")

        return x.upper();

    def format(x):
        
        if isinstance(x, str):
            return "\'"+x+"\'";

        if isinstance(x, bool):
            return "true" if x else "false";

        try: 
            x = str(x)
        except ValueError:
            print(x," is not a valid string") 

        return str(x);

