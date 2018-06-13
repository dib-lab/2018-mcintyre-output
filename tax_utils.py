import csv 
from sourmash.lca import lca_utils, command_index

# Makes acc_to_lineage a global variable 
acc_to_lineage = None 

def load_sourmash_csv(filename):
    with open(filename, 'rt') as fp:
        r = csv.DictReader(fp)
        rows = list(r)
    return rows

notfound = set()
def make_gather_lineages(filename):
    rows = load_sourmash_csv(filename)
    rows2 = []
    for d in rows:
        name = d['name']
        lineage = get_lineage_by_acc(name)
        if lineage is None:
            print('ZZZ found no lineage for {}'.format(name))
            notfound.add(name)
            continue
        lineage = [x for x in lineage if x.rank != "strain"]
        lineage = tuple(lineage)
        d2 = dict(d)
        d2['lineage'] = lineage
        rows2.append(d2)
        
    return rows2

# Gets lineage for an NCBI accession number 
def get_lineage_by_acc(acc):
    global acc_to_lineage
    acc = acc.split(' ')[0].split('.')[0]
    return acc_to_lineage.get(acc, None)

def load_lineages(filename): 
    global acc_to_lineage
    acc_to_lineage, num_rows = command_index.load_taxonomy_assignments(filename, start_column=3)


