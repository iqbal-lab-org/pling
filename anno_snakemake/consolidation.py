#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 16:45:19 2022
@author: daria
"""

import pandas as pd
from pafpy import PafFile
import os
from pathlib import Path
import math

def gene_locations(pafs,genomes):
    output = pd.DataFrame() #dataframe of output
    for genome in genomes:
        with PafFile(f"{pafs}/{genome}.paf") as paf:
            for record in paf:
                gene = str(record.strand) + str(record.qname)
                location = record.tstart
                info = pd.DataFrame.from_dict({"Genome": [genome], "Gene": [gene], "Location": [location]})
                output = pd.concat([output, info], ignore_index=True) #store info in output
    return output

def sort(data,genome_list): #input: data - output from gene_locations; output: dictionary of sorted genomes, keys are genome names
    sorted_genomes={}
    for genome in genome_list:
        sorting = data[data['Genome']==genome]
        sorting = sorting.sort_values(by='Location')
        sequence = sorting['Gene'].tolist()
        sorted_genomes[genome]=sequence
    return sorted_genomes

def to_int(sorting,data):     #read in dictionary of genomes, DataFrame of genes and locations, return dictionary of integer genomes, table from name to integer
    gene_set= list(data['Gene'])
    gene_set = list(dict.fromkeys(gene_set))
    for i in range(len(gene_set)):
        gene_set[i]=gene_set[i][1:]
    gene_set = set(gene_set)
    gene_list = [f"+{gene}" for gene in gene_set]
    rev_gene_list = [f"-{gene}" for gene in gene_set]
    table = {gene:i+1 for i,gene in enumerate(gene_list)}
    table2 = {gene:-i-1 for i,gene in enumerate(rev_gene_list)}
    table.update(table2)
    int_sequences = {}
    for genome in sorting.keys():
        int_sequences[genome] = [table[el] for el in sorting[genome]]
    return int_sequences, table

def make_unimog(unimogpath, int_sequences): #create unimog file of integer sequences
    with open(unimogpath, 'w') as unimog:
        for genome in int_sequences.keys():
            unimog.write(">"+genome+"\n")
            for gene in int_sequences[genome]:
                unimog.write(str(gene)+" ")
            unimog.write(")\n")

pafs=snakemake.params.pafs
genomes=snakemake.params.genomes

data = gene_locations(pafs,genomes)
sorting = sort(data,genomes)
ints, table = to_int(sorting,data)

unimogpath = snakemake.output.unimog
output_map = snakemake.output.map

make_unimog(unimogpath, ints)
pd.Series(table).to_csv(output_map, sep="\t", header=False)
