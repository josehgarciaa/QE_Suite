import qesuite as qes

import qe_io 

syst = "Pd2Se4"
xyz  = qe_io.load_xyz(syst+".xyz");
cell = qe_io.load_cell(syst+".uc");

qes_handler =  qes.generate_from_xyz(xyz=xyz, cell=cell);




from requests_html import HTMLSession

def get_sssp_pseudos():
    session = HTMLSession()

    url = "https://doi.org/10.24435/materialscloud:rz-77"
    r = session.get(url);

    exts = ('.tar.gz','json');
    unique_links= [];
    for ext in exts:
        for elem in r.html.find('a', containing=ext):
            for link in elem.links:
                if ext in link and link not in unique_links:
                    unique_links.append(link);

    base_url = "https://archive.materialscloud.org"
    print(base_url+unique_links[0])
