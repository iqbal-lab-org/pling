#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 15:14:40 2023

@author: daria
"""
import pymummer
from pyfastaq import intervals
from pafpy import PafFile
import pandas as pd
import networkx as nx
import numpy as np

def read_in_unimog(inputpath_unimog):
    sequences = {}
    with open(inputpath_unimog) as unimog:
        genome = ""
        for line in unimog:
            if line[0]==">":
                genome = line.strip(">\n")
            else:
                sequence=line.strip(")\n").split()
                sequences[genome]=[int(el) for el in sequence]
    return sequences

def read_in_nucmer(inputpath_nucmer, genomes):
    block_matches = pd.DataFrame(data=[[[] for genome in genomes] for genome in genomes], index=genomes, columns=genomes)
    for match in pymummer.coords_file.reader(inputpath_nucmer):
        block_matches.loc[match.ref_name, match.qry_name].append(match)
    return block_matches

def read_in_pafs(inputpath_paf, genomes):
    output = {}
    for genome in genomes:
        genes = pd.DataFrame()
        with PafFile(f"{inputpath_paf}/{genome}.paf") as paf:
            for record in paf:
                gene = str(record.strand) + str(record.qname)
                start = record.tstart
                stop = record.tend
                info = pd.DataFrame.from_dict({"Gene": [gene], "Start": [start], "Stop":[stop]})
                genes = pd.concat([genes, info], ignore_index=True)
        genes = genes.sort_values(by='Start', ignore_index=True)
        output[genome] = genes
    return output

def construct_gene_blocks(block_match, gene_matches_of_genome, genome, ref_bool, name_to_int):
    gene_block = []
    for el in gene_matches_of_genome.index:
        gene_name = gene_matches_of_genome.loc[el,'Gene']
        strand = np.sign(name_to_int[gene_name])
        gene = abs(name_to_int[gene_name])
        gene_match = intervals.Interval(gene_matches_of_genome.loc[el,'Start'], gene_matches_of_genome.loc[el,'Stop'])
        if ref_bool:
            if (block_match.ref_coords().contains(gene_match)):
                gene_block.append((genome, gene, el, strand))
        else:
            if (block_match.qry_coords().contains(gene_match)):
                gene_block.append((genome, gene, el, strand))
    return gene_block

def match_dups(ref, query, name_to_int, block_matches, gene_matches, duplicates, G):
    matches = block_matches.loc[ref, query]
    genes_ref = gene_matches[ref]
    genes_query = gene_matches[query]
    for match in matches:
        gene_block_ref = construct_gene_blocks(match, genes_ref, ref, True, name_to_int)
        gene_block_query = construct_gene_blocks(match, genes_query, query, False, name_to_int)
        last_pos_query = -1
        for i in range(len(gene_block_ref)):
            if gene_block_ref[i][1] in duplicates:
                j = last_pos_query + 1
                matching_bool = True
                while j<len(gene_block_query) and matching_bool:
                    if gene_block_ref[i][1] == gene_block_query[j][1]:
                        G.add_edge(gene_block_ref[i], gene_block_query[j])
                        last_pos_query = j
                        matching_bool = False
                    j = j+1

def build_multipartite(name_to_int, block_matches, gene_matches, duplicates, genomes):
    G=nx.Graph()
    for i in range(len(genomes)):
        j=0
        while j<i:
            match_dups(genomes[i], genomes[j], name_to_int, block_matches, gene_matches, duplicates, G)
            j=j+1
    return G

def check_clique(subgraph):
    clique_bool = True
    nodes = subgraph.nodes
    n = len(nodes)
    if n>1:
        for node in nodes:
            if subgraph.degree(node)!=n-1:
                clique_bool = False
                return clique_bool
    else:
        clique_bool = False
    return clique_bool

def relabel_graph(int_to_name, G, label):
    S = [G.subgraph(c).copy() for c in nx.connected_components(G)]
    relabelled_G = G.copy()
    new_label = label+1
    counts = {}
    for subgraph in S:
        clique_bool = check_clique(subgraph)
        if clique_bool==True:
            mapping = {}
            for node in subgraph.nodes:
                mapping[node]=(node[0], new_label, node[2], node[3])
            name = int_to_name[list(subgraph.nodes)[0][1]]
            if name in counts.keys():
                counts[name] = counts[name] + 1
            else:
                counts[name] = 1
            dedup_name = f"{name}_dedup_{counts[name]}"
            int_to_name[new_label] = dedup_name
            int_to_name[-1*new_label] = dedup_name.replace('+','-')
            relabelled_G = nx.relabel_nodes(relabelled_G, mapping, copy=False)
            new_label = new_label + 1
    return relabelled_G

def relabel_sequences(relabelled_G, sequences):
    relabelled_sequences = {}
    for node in relabelled_G.nodes:
        genome = node[0]
        gene = node[1]
        position = node[2]
        strand = node[3]
        if not genome in relabelled_sequences.keys():
            relabelled_sequences[genome]=[gene for gene in sequences[genome]]
        relabelled_sequences[genome][position]=strand*gene
    return relabelled_sequences

def generate_dedup_unimog(outputpath, map_outputpath, relabelled_sequences, int_to_name):
    with open(outputpath, 'w') as new_unimog:
        for genome in relabelled_sequences.keys():
            new_unimog.write(">"+genome+"\n")
            for gene in relabelled_sequences[genome]:
                new_unimog.write(str(gene)+" ")
            new_unimog.write(")\n")
    with open(map_outputpath, 'w') as dedup_map:
        for integer in int_to_name.keys():
            dedup_map.write(int_to_name[integer]+"\t"+str(integer)+"\n")

def find_duplicates(sequences):
    duplicates = set()
    for genome in sequences.keys():
        absolute = [abs(el) for el in sequences[genome]]
        for el in absolute:
            if absolute.count(el)>1:
                duplicates.add(el)
    return duplicates

def deduplication(inputpath_unimog, inputpath_nucmer, inputpath_paf, inputpath_map, genomes, outputpath, map_outputpath):
    name_and_int = pd.read_csv(inputpath_map,index_col=None, sep='\t', header=None)
    name_and_int = name_and_int.values.tolist()
    name_to_int = {el[0]:el[1] for el in name_and_int}
    int_to_name = {el[1]:el[0] for el in name_and_int}
    sequences = read_in_unimog(unimog)
    block_matches = read_in_nucmer(nucmer, genomes)
    gene_matches = read_in_pafs(pafs, genomes)
    duplicates = find_duplicates(sequences)
    label = max([max([abs(el) for el in sequences[genome]]) for genome in genomes])
    if duplicates:
        multipartite = build_multipartite(name_to_int, block_matches, gene_matches, duplicates, genomes)
        relabelled_multipartite = relabel_graph(int_to_name, multipartite, label)
        relabelled_sequences = relabel_sequences(relabelled_multipartite,sequences)
        generate_dedup_unimog(outputpath, map_outputpath, relabelled_sequences, int_to_name)
    else:
        generate_dedup_unimog(outputpath, map_outputpath, sequences, int_to_name)

testing = False

if testing == True:
    OUTPUTPATH =  "/home/daria/Documents/projects/pling/tests/test1/annotation"
    fasta = f"{OUTPUTPATH}/tmp_files/community_fastas/0.fna"
    nucmer = f"{OUTPUTPATH}/tmp_files/nucmer/0.nucmer"
    nucmer_threshold = 98.5
    unimog = f"{OUTPUTPATH}/unimogs/0_anno.unimog"
    inputpath_map = f"{OUTPUTPATH}/0_map.txt"
    pafs = f"{OUTPUTPATH}/tmp_files/minimap/0/output"
    genomes = [el[0] for el in pd.read_csv("/home/daria/Documents/projects/pling/tests/test1/fastas/list.txt", header=None).values]
    outputpath = f"{OUTPUTPATH}/unimogs/relabelled/dedup/0_dedup.unimog"
    map_outputpath = f"{OUTPUTPATH}/unimogs/relabelled/dedup/0_map_dedup.txt"
else:
    fasta = snakemake.input.fasta
    nucmer = snakemake.output.nucmer
    nucmer_threshold = snakemake.params.nucmer_threshold
    unimog = snakemake.input.unimog
    inputpath_map = snakemake.input.map
    pafs = snakemake.params.pafs
    genomes = snakemake.params.genomes
    outputpath = snakemake.output.relabelled_unimog
    map_outputpath = snakemake.output.dedup_map

runner = pymummer.nucmer.Runner(
                    fasta,
                    fasta,
                    nucmer,
                    min_id=nucmer_threshold,     # minimum percent identity of match
                    promer=False,
                    maxmatch=True
                )
runner.run()

deduplication(unimog, nucmer, pafs, inputpath_map, genomes, outputpath, map_outputpath)
