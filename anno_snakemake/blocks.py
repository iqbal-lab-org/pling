#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import networkx as nx
import numpy as np
import os

def read_in_unimog(unimog_filepath):
    sequences = {}
    with open(unimog_filepath) as unimog:
        genome = ""
        for line in unimog:
            if line[0]==">":
                genome = line.strip(">\n")
            else:
                sequence=line.strip(")\n").split()
                sequences[genome]=[int(el) for el in sequence]
    return sequences

def read_in_map(map_filepath):
    int_to_name = {}
    with open(map_filepath, "r") as int_map:
        for line in int_map:
            name, integer = line.strip().split("\t")
            integer = int(integer)
            int_to_name[integer] = name
    return int_to_name

def construct_graph(seqs, community):
    graph = nx.DiGraph()
    for genome in community:
        for j in range(len(seqs[genome])-1):
            graph.add_edge(seqs[genome][j], seqs[genome][j+1])
        graph.add_edge(seqs[genome][-1], seqs[genome][0])
    return graph

def naming(block, name, block_maps, strand):
    if strand == '+':
        block_maps[name] = [block[0],block]
    elif strand == '-':
        block_maps[name] = [block[-1],block]

def reverse_path(path):
    rev_path = []
    for i in range(len(path)-1,-1,-1):
        rev_path.append(-path[i])
    return rev_path

def forward_block(graph, gene, nodes):
    neighbour = [n for n in graph[gene]]
    block = []
    block.append(gene)
    nodes.remove(gene)
    while graph.out_degree(gene)==1 and graph.in_degree(neighbour[0])==1 and neighbour[0] in nodes:
        gene = neighbour[0]
        neighbour = [n for n in graph[gene]]
        block.append(gene)
        nodes.remove(gene)
    return nodes, block

def reverse_block(graph, gene, nodes):
    block = []
    neighbour = [n for n in graph.predecessors(gene)]
    if graph.in_degree(gene)==1:
        while graph.in_degree(gene)==1 and graph.out_degree(neighbour[0])==1 and neighbour[0] in nodes:
            gene = neighbour[0]
            neighbour = [n for n in graph.predecessors(gene)]
            block.insert(0,gene)
            nodes.remove(gene)
    return nodes, block

def double_forward_block(graph, gene, nodes):
    neighbour = [n for n in graph[gene]]
    block = [gene]
    nodes.remove(gene)
    if graph.has_node(-gene)==False:
        return nodes, block
    neg_neighbour = [n for n in graph.predecessors(-gene)]
    nodes.remove(-gene)
    while graph.out_degree(gene)==1 and graph.in_degree(neighbour[0])==1 and graph.in_degree(-gene)==1 and graph.out_degree(neg_neighbour[0])==1 and neighbour[0]==-neg_neighbour[0]:
        gene = neighbour[0]
        neighbour = [n for n in graph[gene]]
        block.append(gene)
        nodes.remove(gene)
        if graph.has_node(-gene)==False:
            return nodes, block
        neg_neighbour = [n for n in graph.predecessors(-gene)]
        nodes.remove(-gene)
    return nodes, block

def double_reverse_block(graph, gene, nodes):
    block = []
    neighbour = [n for n in graph.predecessors(gene)]
    if graph.has_node(-gene)==False:
        return nodes, block
    neg_neighbour = [n for n in graph[-gene]]
    if graph.in_degree(gene)==1:
        while graph.in_degree(gene)==1 and graph.out_degree(neighbour[0])==1 and graph.out_degree(-gene)==1 and graph.in_degree(neg_neighbour[0])==1 and neighbour[0]==-neg_neighbour[0]:
            gene = neighbour[0]
            neighbour = [n for n in graph.predecessors(gene)]
            block.insert(0,gene)
            nodes.remove(gene)
            if graph.has_node(-gene)==False:
                return nodes, block
            neg_neighbour = [n for n in graph[-gene]]
            nodes.remove(-gene)
    return nodes, block

def blocks(seqs, community, int_to_name):
    graph = construct_graph(seqs, community)
    nodes = list(nx.nodes(graph))
    nodes_2 = []
    blocks=[]
    block_maps={}
    i=0
    while nodes != []:
        start = nodes[0]
        nodes, for_block = double_forward_block(graph, start, nodes)
        nodes, rev_block = double_reverse_block(graph, start, nodes)
        block = rev_block+for_block
        #renaming
        if len(block)==1:
            if graph.has_node(-block[0]):
                name=int_to_name[block[0]]
                naming(block, name, block_maps, '+')
                name=int_to_name[-block[0]]
                neg_block = [-block[0]]
                naming(neg_block, name, block_maps, '-')
                blocks.append(block)
            else:
                nodes_2.append(block[0])
        else:
            name = "+block_"+str(i)
            naming(block, name, block_maps, '+')
            name = "-block_"+str(i)
            neg_block = reverse_path(block)
            naming(neg_block, name, block_maps, '-')
            i=i+1
            blocks.append(block)
    #second pass over genes only present on one strand
    while nodes_2 != []:
        start = nodes_2[0]
        nodes_2, for_block = forward_block(graph, start, nodes_2)
        nodes_2, rev_block = reverse_block(graph, start, nodes_2)
        block = rev_block+for_block
        blocks.append(block)
        if len(block)==1: #renaming
            name=int_to_name[block[0]]
            naming(block, name, block_maps, '+')
        else:
            name = "+block_"+str(i)
            naming(block, name, block_maps, '+')
            i=i+1
    return blocks, block_maps

def output_files(seqs, community, blocks, block_maps, int_to_name, output_map, relabelled_unimog, directory):
    new_seqs={}
    for genome in community:
        new_seqs[genome]=[]
        j=0
        while j < len(seqs[genome]):
            n=1
            for block in blocks:
                if block[0]==seqs[genome][j]:
                    new_seqs[genome].append(block[0])
                    n = len(block)
                elif -block[-1]==seqs[genome][j]:
                    new_seqs[genome].append(-block[-1])
                    n=len(block)
            j = j+n
    maps = pd.DataFrame(data=block_maps)
    maps.T.to_csv(output_map, sep="\t", header=False)
    if not os.path.exists(directory):
        os.makedirs(directory)
    with open(relabelled_unimog, 'w') as f:
        for genome in community:
            unimog = " ".join(str(el) for el in new_seqs[genome])
            unimog = unimog + " )"
            f.write(f">{genome}\n{unimog}\n")

testing = False

if testing == True:
    OUTPUTPATH = "/home/daria/Documents/projects/pling/local_tests"
    map_filepath = f"{OUTPUTPATH}/0_map.txt"
    unimog_filepath = f"{OUTPUTPATH}/unimogs/0_anno.unimog"
    relabelled_unimog = f"{OUTPUTPATH}/unimogs/relabelled/blocks/0_blocks.unimog"
    output_map = f"{OUTPUTPATH}/unimogs/relabelled/blocks/0_map_blocks.txt"
    genomes = [el[0] for el in pd.read_csv("/home/daria/Documents/projects/Murray_Family/clusters_adrian/lists/cluster_4_ids.csv", header=None).values]
else:
    map_filepath = snakemake.input.map
    unimog_filepath = snakemake.input.unimog
    relabelled_unimog = snakemake.output.relabelled_unimog
    relabelled_dir = snakemake.params.relabelled_dir
    output_map = snakemake.output.blocks_map
    genomes = snakemake.params.genomes

int_to_name = read_in_map(map_filepath)
seqs = read_in_unimog(unimog_filepath)
blocks, block_maps = blocks(seqs, genomes, int_to_name)
output_files(seqs, genomes, blocks, block_maps, int_to_name, output_map, relabelled_unimog, relabelled_dir)
