import qe_suite.format as f

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

