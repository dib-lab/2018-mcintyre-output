from __future__ import print_function
import csv 
from sourmash.lca import lca_utils, command_index
import sys
sys.path.insert(0, '../2018-ncbi-lineages/')
import ncbi_taxdump_utils
import os
import matplotlib.pyplot as plt
import seaborn as sns
from pprint import pprint
from sourmash.lca import lca_utils, command_index
from clustergrammer_widget import *
from ipywidgets import interact, interactive, fixed, interact_manual
import matplotlib as mpl
import pandas as pd
import glob
import qgrid
import numpy as np

# Makes acc_to_lineage a global variable 

def get_ncbi_lineages():
    taxfoo = ncbi_taxdump_utils.NCBI_TaxonomyFoo()
    taxfoo.load_nodes_dmp('../2018-ncbi-lineages/genbank/nodes.dmp')
    taxfoo.load_names_dmp('../2018-ncbi-lineages/genbank/names.dmp')
    want_taxonomy = ['superkingdom', 'phylum', 'order', 'class', 'family', 'genus', 'species']
    return taxfoo, want_taxonomy

def format_lineage(lineage_tup):
    return ";".join(lca_utils.zip_lineage(lineage_tup))

def load_gather_lineages(filename):
    acc_to_lineage, num_rows = command_index.load_taxonomy_assignments(filename, start_column=3)
    return acc_to_lineage, num_rows
    
def load_sourmash_csv(filename):
    with open(filename, 'rt') as fp:
        r = csv.DictReader(fp)
        rows = list(r)
    return rows

def get_lineage_by_acc(acc):
    acc_to_lineage, num_rows = load_gather_lineages('gather-lineages.csv')
    acc = acc.split(' ')[0].split('.')[0]
    return acc_to_lineage.get(acc, None)

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

def load_truth_file(filename):
    with open(filename, 'rt') as fp:
        lines = fp.readlines()
        
    lines = [ x.strip() for x in lines ]
    lines = [ x.split('\t') for x in lines ]
    
    rows = []
    for x in lines:
        taxid, a, b, rank, name = x
        taxid = int(taxid)
        rows.append((taxid, a, b, rank, name))
    return rows

def make_lineage_from_taxid(taxid):
    taxfoo = ncbi_taxdump_utils.NCBI_TaxonomyFoo()
    taxfoo.load_nodes_dmp('../2018-ncbi-lineages/genbank/nodes.dmp')
    taxfoo.load_names_dmp('../2018-ncbi-lineages/genbank/names.dmp')
    want_taxonomy = ['superkingdom', 'phylum', 'order', 'class', 'family', 'genus', 'species']
    lineage_d = taxfoo.get_lineage_as_dict(taxid, want_taxonomy)
    
    lineage = []
    for rank in lca_utils.taxlist():
        if rank == 'strain': continue 
        name = lineage_d.get(rank, 'unassigned')
        lineage_pair = lca_utils.LineagePair(rank, name)
        lineage.append(lineage_pair)
    return tuple(lineage)

def make_truth_lineages(filename):
    rows = load_truth_file(filename)
    
    rows2 = []
    for (taxid, a, b, rank, name) in rows:
        lineage = make_lineage_from_taxid(taxid)
        rows2.append((taxid, a, b, rank, name, lineage))
        for lintup in lineage:
            if lintup.rank == rank:
                if lintup.name != name:
                    print('DISAGREE: ncbi={}, truthfile={}'.format(lintup.name, name))
    return rows2

truth_lineages = make_truth_lineages('truth_sets/species/Huttenhower_HC1_TRUTH.txt')

def make_lca_gather_lineages(filename):
    rows = load_sourmash_csv(filename)
    rows2 = []
    for d in rows:
        d2 = dict(d)
        
        lineage = []
        for rank in lca_utils.taxlist():
            if rank in d2:
                name = d2.get(rank)
                del d2[rank]
                lineage.append((rank, name))
                
        lineage = [ lca_utils.LineagePair(r, n) for (r, n) in lineage ]
        d2['lineage'] = tuple(lineage)
        rows2.append(d2)
        
    return rows2

def compare_gather_to_truth(truth_file, gather_csv):
    truth = make_truth_lineages(truth_file)
    gather = make_gather_lineages(gather_csv)
    
    truth_lineages = set([ t[5] for t in truth ])
    gather_lineages = set([ row['lineage'] for row in gather ])

    # Print intersection and union     
    print(len(truth_lineages.intersection(gather_lineages)))
    print(len(truth_lineages.union(gather_lineages)))
    
    print(len(truth_lineages))
    
    # Compare gather to truth and print 'in gather but not truth' and 'in truth but not gather'
    print('** in gather but not truth:')
    for diff in gather_lineages - truth_lineages:
        print('\t', format_lineage(diff))
    
    print('\n** in truth but not gather:')
    for diff in truth_lineages - gather_lineages:
        print('\t', format_lineage(diff))
        
def compare_lca_gather_to_truth(truth_file, gather_csv):
    truth = make_truth_lineages(truth_file)
    gather = make_lca_gather_lineages(gather_csv)
    
    truth_lineages = set([ t[5] for t in truth ])
    gather_lineages = set([ row['lineage'] for row in gather ])

    # Print intersection and union     
    print(len(truth_lineages.intersection(gather_lineages)))
    print(len(truth_lineages.union(gather_lineages)))
    
    print(len(truth_lineages))
    
    # Compare gather to truth and print 'in gather but not truth' and 'in truth but not gather'
    print('** in gather but not truth:')
    for diff in gather_lineages - truth_lineages:
        print('\t', format_lineage(diff))
    
    print('\n** in truth but not gather:')
    for diff in truth_lineages - gather_lineages:
        print('\t', format_lineage(diff))
        
def compare_lca_gather_to_gather(lca_gather_csv, gather_csv):
    reg_gather = make_gather_lineages(gather_csv)
    lca_gather = make_lca_gather_lineages(lca_gather_csv)
    
    reg_gather_lineages = set([ row['lineage'] for row in reg_gather ])
    lca_gather_lineages = set([ row['lineage'] for row in lca_gather ])
    
    print(len(lca_gather_lineages.intersection(reg_gather_lineages)))
    print(len(lca_gather_lineages.union(reg_gather_lineages)))
        
    print('** in lca_gather but not reg gather:')
    for diff in lca_gather_lineages - reg_gather_lineages:
        print('\t', format_lineage(diff))
    
    print('\n** in gather but not lca gather:')
    for diff in reg_gather_lineages - lca_gather_lineages:
        print('\t', format_lineage(diff))

def compare_lca_gather_reads_to_truth_fp(truth_file, gather_csv):
    truth = make_truth_lineages(truth_file)
    gather = make_lca_gather_lineages(gather_csv)
    
    truth_lineages = set([ t[5] for t in truth ])
    gather_lineages = set([ row['lineage'] for row in gather ])

    return len(gather_lineages - truth_lineages)

def compare_gather_reads_to_truth_fp(truth_file, gather_csv):
    truth = make_truth_lineages(truth_file)
    gather = make_lca_gather_lineages(gather_csv)
    
    truth_lineages = set([ t[5] for t in truth ])
    gather_lineages = set([ row['lineage'] for row in gather ])
    
    return len(gather_lineages - truth_lineages)

def compare_lca_gather_contigs_to_truth_fp(truth_file, gather_csv):
    truth = make_truth_lineages(truth_file)
    gather = make_lca_gather_lineages(gather_csv)
    
    truth_lineages = set([ t[5] for t in truth ])
    gather_lineages = set([ row['lineage'] for row in gather ])
    
    return len(gather_lineages - truth_lineages)

def compare_lca_gather_reads_to_truth_tp(truth_file, gather_csv):
    truth = make_truth_lineages(truth_file)
    gather = make_lca_gather_lineages(gather_csv)
    
    truth_lineages = set([ t[5] for t in truth ])
    gather_lineages = set([ row['lineage'] for row in gather ])
    
    return len(truth_lineages.intersection(gather_lineages))

def compare_gather_reads_to_truth_tp(truth_file, gather_csv):
    truth = make_truth_lineages(truth_file)
    gather = make_gather_lineages(gather_csv)
    
    truth_lineages = set([ t[5] for t in truth ])
    gather_lineages = set([ row['lineage'] for row in gather ])
    

    return len(truth_lineages.intersection(gather_lineages))

def compare_lca_gather_contigs_to_truth_tp(truth_file, gather_csv):
    truth = make_truth_lineages(truth_file)
    gather = make_lca_gather_lineages(gather_csv)
    
    truth_lineages = set([ t[5] for t in truth ])
    gather_lineages = set([ row['lineage'] for row in gather ])
    

    return len(truth_lineages.intersection(gather_lineages))

def compare_lca_gather_reads_to_truth_fn(truth_file, gather_csv):
    truth = make_truth_lineages(truth_file)
    gather = make_lca_gather_lineages(gather_csv)
    
    truth_lineages = set([ t[5] for t in truth ])
    gather_lineages = set([ row['lineage'] for row in gather ])
    
    return len((truth_lineages) - (truth_lineages.intersection(gather_lineages)))

def compare_lca_gather_contigs_to_truth_fn(truth_file, gather_csv):
    truth = make_truth_lineages(truth_file)
    gather = make_lca_gather_lineages(gather_csv)
    
    
    truth_lineages = set([ t[5] for t in truth ])
    gather_lineages = set([ row['lineage'] for row in gather ])
        
    return len((truth_lineages) - (truth_lineages.intersection(gather_lineages)))

def compare_gather_reads_to_truth_fn(truth_file, gather_csv):
    truth = make_truth_lineages(truth_file)
    gather = make_gather_lineages(gather_csv)

    truth_lineages = set([ t[5] for t in truth ])
    gather_lineages = set([ row['lineage'] for row in gather ])
    
    return len(truth_lineages - gather_lineages)
