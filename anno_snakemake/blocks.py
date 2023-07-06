#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import networkx as nx
import numpy as np
import math
import matplotlib.pyplot as plt

def read_in_unimog(unimog_filepath):
    seqs = {}
    with open(unimog_filepath, "r") as unimog:
        lines = unimog.read().split(")\n")
        for line in lines:
            content = line.split("\n")
            genome_id = content[0][1:]
            sequence = [int(el) for el in content[1].split()]
            seqs[genome_id] = sequence
    return seqs

def read_in_map(map_filepath):
    int_to_name = {}
    with open(map_filepath, "r") as map:
        for line in map:
            name, int = line.strip().split("\t")
            int = int(int)
            int_to_name[int] = name
    return int_to_name

def get_communities(communities_filepath):
    communities = {}
    with open(communities_filepath) as communities_fh:
        for community_index, line in enumerate(communities_fh):
            plasmids = line.strip().split()
            communities[community_index] = plasmids

def contruct_graph(seqs, community):
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

def output_files(seqs, community, blocks, block_maps, int_to_name, output_map, relabelled_unimog):
    new_seqs={}
    for genome in community:
        new_seqs[genome]=[]
        for el in seqs[genome]:
            new_seqs[genome].append(block_maps[int_to_name[el]][0])
    maps = pd.DataFrame(data=block_maps)
    maps.T.to_csv(output_map, sep="\t", header=False)
    with open(relabelled_unimog) as f:
        for genome in community:
            unimog = " ".join(str(el) for el in new_seqs[genome])
            unimog = unimog + " )"
            f.write(f">{genome}\n{unimog}\n")

map_filepath = snakemake.input.map
communties_filepath = snakemake.input.communties
unimog_filepath = snakemake.input.unimog
output_dir = snakemake.output.relabelled_dir

int_to_name = read_in_map(map_filepath)
communties = get_communities(communties_filepath)
seqs = read_in_unimog(unimog_filepath)
for community in communities:
    blocks, block_maps = blocks(seqs, communties[community], int_to_name)
    output_files(seqs, community, blocks, block_maps, int_to_name, f"{output_dir}/{community}_blocks_map.txt", f"{output_dir}/{community}_blocks.unimog")
