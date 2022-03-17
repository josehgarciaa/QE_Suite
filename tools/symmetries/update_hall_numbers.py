
import json
spghnum =None;#obtained form http://cci.lbl.gov/sginfo/itvb_2001_table_a1427_hall_symbols.html

#spglib convention in https://spglib.github.io/spglib/definition.html
#qe_convention in https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm200
with open('spg_hall_numbers.json') as fp:
    spghnum= json.loads(fp.read());

for k,v in spghnum.items():
    spgnum = int(v["Number"].split(":")[0]);
    spghnum[k]["ibrav"]= None
    spghnum[k]["qe_comp"]= False;
    spghnum[k]["spglib_comp"]= False;

    is_cubic = 195 <=spgnum <=230;
    is_hexagonal = 168 <=spgnum <=194;
    is_trigonal= 143 <=spgnum <=167;
    is_tetragonal= 75 <=spgnum <= 142;
    is_orthorhombic= 16 <=spgnum <= 74
    is_monoclinic= 3 <=spgnum <= 15;
    is_triclinic = 1 <=spgnum <= 2;    
    if is_cubic:
        if "P" in v["Hermann-Mauguin"]:
           spghnum[k]["ibrav"]= 1
           spghnum[k]["qe_comp"]= True;
           spghnum[k]["spglib_comp"]= True;

        if "F" in v["Hermann-Mauguin"]:
            spghnum[k]["ibrav"]= 2
            spghnum[k]["qe_comp"]= False;
            spghnum[k]["spglib_comp"]= False;
        if "I" in v["Hermann-Mauguin"]:
            spghnum[k]["ibrav"]=-3 #We are using the "most symmetric" axis
            spghnum[k]["qe_comp"]= False;
            spghnum[k]["spglib_comp"]= False;

    if is_hexagonal:
        if "P" in v["Hermann-Mauguin"]:
            spghnum[k]["ibrav"]= 4
            spghnum[k]["qe_comp"]= True;
            spghnum[k]["spglib_comp"]= True;

    if is_trigonal:
        if "P" in v["Hermann-Mauguin"]:
            spghnum[k]["ibrav"]= 4
            spghnum[k]["qe_comp"]= False;
            spghnum[k]["spglib_comp"]= False;
        if "R" in v["Hermann-Mauguin"]:
            spghnum[k]["ibrav"]=5
            spghnum[k]["qe_comp"]= False;
            spghnum[k]["spglib_comp"]= False;

    if is_tetragonal :
        if "P" in v["Hermann-Mauguin"]:
            spghnum[k]["ibrav"]= 6
            spghnum[k]["qe_comp"]= True;
            spghnum[k]["spglib_comp"]= True;
        if "I" in v["Hermann-Mauguin"]:
            spghnum[k]["ibrav"]= 7
            spghnum[k]["qe_comp"]= False;
            spghnum[k]["spglib_comp"]= False;

    if is_orthorhombic :
        if "P" in v["Hermann-Mauguin"]:
            spghnum[k]["ibrav"]= 8
            spghnum[k]["qe_comp"]= True;
            spghnum[k]["spglib_comp"]= True;
        if "C" in v["Hermann-Mauguin"]:
            spghnum[k]["ibrav"]= 9
            spghnum[k]["qe_comp"]= False;
            spghnum[k]["spglib_comp"]= False;

        if "A" in v["Hermann-Mauguin"]:
            spghnum[k]["ibrav"]= 91
            spghnum[k]["qe_comp"]= False;
            spghnum[k]["spglib_comp"]= False;
        if "F" in v["Hermann-Mauguin"]:
            spghnum[k]["ibrav"]= 10
            spghnum[k]["qe_comp"]= False;
            spghnum[k]["spglib_comp"]= False;
        if "I" in v["Hermann-Mauguin"]:
            spghnum[k]["ibrav"]= 11
            spghnum[k]["qe_comp"]= False;
            spghnum[k]["spglib_comp"]= False;

    if is_monoclinic:
        if "P" in v["Hermann-Mauguin"]:
            spghnum[k]["ibrav"]=-12
            spghnum[k]["qe_comp"]= True;
            spghnum[k]["spglib_comp"]= True;
        if "C" in v["Hermann-Mauguin"]:
            spghnum[k]["ibrav"]= 13
            spghnum[k]["qe_comp"]= False;
            spghnum[k]["spglib_comp"]= False;
        if "B" in v["Hermann-Mauguin"]:
            spghnum[k]["ibrav"]=-13
            spghnum[k]["qe_comp"]= False;
            spghnum[k]["spglib_comp"]= False;

    if is_triclinic:
        spghnum[k]["ibrav"]=14
        spghnum[k]["qe_comp"]= False;
        spghnum[k]["spglib_comp"]= False;


with open('spg_hall_numbers.json',"w") as fp:
    spghnum= json.dump(spghnum,fp);

