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
        with PafFile(pafs+'/' +genome+'.paf') as paf:
            for record in paf:
                gene = str(record.strand) + str(record.qname)
                location = record.tstart
                info = pd.DataFrame.from_dict({"Genome": [genome], "Gene": [gene], "Location": [location]})
                output = pd.concat([output, info], ignore_index=True) #store info in output
    return output


def sort_and_lineup_single(data,genome_id,lineup_id): #input: data - output from gene_locations, genome_id - name of genome as string, lineup_id: gene_id at which to align sorted genome; output: genome - a sorted list of genes
    #sort the genes of one genome, and then align the sorted list at one preselected gene
    sorting = data[data['Genome']==genome_id]
    sorting = sorting.sort_values(by='Location')
    genome = sorting['Gene'].tolist()
    lineup = 0
    if lineup_id is not None:
        try:
            lineup = genome.index('+'+lineup_id)
        except:
            lineup = genome.index('-'+lineup_id)
    genome1 = genome[:lineup]
    genome2 = genome[lineup:]
    genome = genome2 + genome1
    return genome

def sort_and_lineup(data,genome_list,lineup_id): #input: data - output from gene_locations, cluster - name of cluster; output: DataFrame of sorted genomes, row index is genome name
    sorted_genomes=pd.DataFrame()
    for el in genome_list:
        genome=pd.DataFrame([sort_and_lineup_single(data,el,lineup_id)],index=[el])
        sorted_genomes=pd.concat([sorted_genomes,genome])
    return sorted_genomes

def name_to_int(table, x):
        if type(x) == float:
            return math.nan
        else:
            return table[x]

def to_int(sorting,data):     #read in dataframe of genomes, return integer matrix of genomes, table from name to integer
    gene_set= list(data['Gene'])
    gene_set = list(dict.fromkeys(gene_set))
    for i in range(len(gene_set)):
        gene_set[i]=gene_set[i][1:]
    gene_list = ['+'+gene for gene in gene_set]
    rev_gene_list = ['-'+gene for gene in gene_set]
    table = {gene:i+1 for i,gene in enumerate(gene_list)}
    table2 = {gene:-i-1 for i,gene in enumerate(rev_gene_list)}
    table.update(table2)
    namemap = lambda x: name_to_int(table,x)
    sorting = sorting.applymap(namemap)
    return sorting, table

def lineup_id(gene_matrix, topology):
    blub = 1
    if topology == 'lin':
        print('Linear plasmid present in lineage.')
        return None
    for gene in gene_matrix.index:
        for genome in gene_matrix.columns:
            if gene_matrix.loc[gene, genome] == 0:
                blub = 0
        if blub==1:
            return gene
    print('No core genes!')
    return None


pafs=snakemake.params.pafs
cluster=snakemake.params.cluster
lineage=snakemake.params.lineage
genomes=snakemake.params.genomes
murraypath=snakemake.params.murraypath
topology=snakemake.params.topology

gene_matrix = pd.read_csv(murraypath+"/pangenome/PANAROO-per_lineage/results/c"+cluster+"/"+lineage+"/gene_presence_absence.Rtab", sep='\t', index_col=0)
data = gene_locations(pafs,genomes)
lineup_id=lineup_id(gene_matrix, topology)
sorting = sort_and_lineup(data,genomes,lineup_id)
ints, table = to_int(sorting,data)

output0 = snakemake.output[0]
output1 = snakemake.output[1]
output2 = snakemake.output[2]

#output0 = cluster+'_names.txt'
#output1 = cluster+'_ints.txt'
#output2 = cluster+'_map.txt'

#os.mkdir(cluster+'/genomes')
sorting.to_csv(output0, sep="\t", header=False)
ints.to_csv(output1, sep="\t", header=False)
pd.Series(table).to_csv(output2, sep="\t")
