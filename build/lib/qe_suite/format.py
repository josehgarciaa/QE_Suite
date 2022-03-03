def remove_double(c, string):
    c=str(c);
    while c+c in string:
        string = string.replace('  ', ' ')
    return string
def remove_empty(x):
    try:
        xiter = iter(x)
    except TypeError:
        print("Error")
    else:
        return [x for x in xiter if len(x)!=0 ];
