#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 17:14:44 2022

@author: daria
"""
import pandas as pd
import networkx as nx
import numpy as np
import math
import matplotlib.pyplot as plt


def naming(block, name, name_and_block, block_name_and_block_int, block_name_and_int, strand):
    if strand == '+':
        block_name_and_block_int.append([name,block[0]])
    elif strand == '-':
        block_name_and_block_int.append([name,block[-1]])
    name_and_block.append([name,block])
    for el in block:
        block_name_and_int.append([name,el])
        
def reverse_path(path):
    rev_path = []
    for i in range(len(path)-1,-1,-1):
        rev_path.append(-path[i])       
    return rev_path

class GeneSeq:
    
    def __init__(self, genomes, genome_ids, name_and_int = [], graph = False): #genomes - numpy matrix of genomes, genome_ids - ordered list of genomes (ordered the same as genomes matrix), name_and_int - list of tuples/lists (name,int)
        self.genomes = genomes
        self.genome_ids = genome_ids
        self.name_to_int = {}
        self.int_to_name = {}
        self.graph = nx.DiGraph()
        for el in name_and_int:
            self.name_to_int[el[0]] = el[1]
            self.int_to_name[el[1]] = el[0]
        if graph == True:
            self.construct_graph()            
            
    def __eq__(self, other):
        genomes_bool = np.all((self.genomes == other.genomes) | np.isnan(self.genomes) | np.isnan(other.genomes))
        ids_bool = (self.genome_ids, other.genome_ids)
        graph_bool = True
        if nx.is_empty(self.graph) == False and nx.is_empty(other.graph) == False:
            graph_bool = (self.graph.edges==other.graph.edges) and (self.graph.nodes==other.graph.nodes)
        return genomes_bool and ids_bool and graph_bool

    @classmethod
    def from_file(cls, path):
        pass
    
    @classmethod
    def from_cluster(cls, path_to_cluster, cluster, graph = False):
        path = path_to_cluster
        genomes = pd.read_csv(path+'/'+cluster+'_ints.txt', header = None, index_col=0, sep='\t')
        genome_ids = genomes.index
        name_to_int = pd.read_csv(path+'/'+cluster+'_map.txt',index_col=None, sep='\t')
        genomes=genomes.to_numpy()
        name_and_int = name_to_int.values.tolist()
        return cls(genomes, genome_ids, name_and_int, graph)
    
    def construct_graph(self, strand = True):
        if not nx.is_empty(self.graph):
            self.graph=nx.empty_graph()
        if strand == False:
            con_genomes = np.absolute(self.genomes)
        else:
            con_genomes = self.genomes
        for i in range(len(con_genomes)):
            for j in range(len(con_genomes[0])):
                if j==len(con_genomes[0])-1:
                    if not math.isnan(con_genomes[i][j]):
                        self.graph.add_edge(con_genomes[i][j],con_genomes[i][0])                        
                elif math.isnan(con_genomes[i][j+1]) and not math.isnan(con_genomes[i][j]):
                    self.graph.add_edge(con_genomes[i][j],con_genomes[i][0])
                elif not math.isnan(con_genomes[i][j+1]) and not math.isnan(con_genomes[i][j]):
                    self.graph.add_edge(con_genomes[i][j],con_genomes[i][j+1])
                    
    def graph_vis(self, strand = True, size = (50,50), options = None):
        if nx.is_empty(self.graph) == True or strand == False:
            self.construct_graph(strand)
        plt.figure(figsize=size)
        plt.plot()
        if options == None:
            options = {
            
                'node_color': 'red',
            
                'node_size': 100,
            
                'width': 2,
            
            }
        nx.draw_kamada_kawai(self.graph,**options)
    
    def to_cluster(self):
        pass
    
    def has_reverse_path(self, path):
        for i in range(len(path)-1,0,-1):
            tmp_bool = self.graph.has_edge(-path[i],-path[i-1])
            if tmp_bool == False:
                return False
        return True
    
    def forward_block(self, gene, nodes):
        neighbour = [n for n in self.graph[gene]]
        block = []
        block.append(gene)
        nodes.remove(gene)
        while self.graph.out_degree(gene)==1 and self.graph.in_degree(neighbour[0])==1 and neighbour[0] in nodes:
            gene = neighbour[0]
            neighbour = [n for n in self.graph[gene]]
            block.append(gene)
            nodes.remove(gene)
        return nodes, block

    def reverse_block(self, gene, nodes):
        block = []
        neighbour = [n for n in self.graph.predecessors(gene)]
        if self.graph.in_degree(gene)==1:
            while self.graph.in_degree(gene)==1 and self.graph.out_degree(neighbour[0])==1 and neighbour[0] in nodes:
                gene = neighbour[0]
                neighbour = [n for n in self.graph.predecessors(gene)]
                block.insert(0,gene)
                nodes.remove(gene)
        return nodes, block

    def double_forward_block(self, gene, nodes):
        neighbour = [n for n in self.graph[gene]]
        block = [gene]
        nodes.remove(gene)
        if self.graph.has_node(-gene)==False:
            return nodes, block
        neg_neighbour = [n for n in self.graph.predecessors(-gene)]
        nodes.remove(-gene)
        while self.graph.out_degree(gene)==1 and self.graph.in_degree(neighbour[0])==1 and self.graph.in_degree(-gene)==1 and self.graph.out_degree(neg_neighbour[0])==1 and neighbour[0]==-neg_neighbour[0]:
            gene = neighbour[0]
            neighbour = [n for n in self.graph[gene]]
            block.append(gene)
            nodes.remove(gene)
            if self.graph.has_node(-gene)==False:
                return nodes, block
            neg_neighbour = [n for n in self.graph.predecessors(-gene)]
            nodes.remove(-gene)
        return nodes, block

    def double_reverse_block(self, gene, nodes):
        block = []
        neighbour = [n for n in self.graph.predecessors(gene)]
        if self.graph.has_node(-gene)==False:
            return nodes, block
        neg_neighbour = [n for n in self.graph[-gene]]
        if self.graph.in_degree(gene)==1:
            while self.graph.in_degree(gene)==1 and self.graph.out_degree(neighbour[0])==1 and self.graph.out_degree(-gene)==1 and self.graph.in_degree(neg_neighbour[0])==1 and neighbour[0]==-neg_neighbour[0]:
                gene = neighbour[0]
                neighbour = [n for n in self.graph.predecessors(gene)]
                block.insert(0,gene)
                nodes.remove(gene)
                if self.graph.has_node(-gene)==False:
                    return nodes, block
                neg_neighbour = [n for n in self.graph[-gene]]
                nodes.remove(-gene)
        return nodes, block
    
    def blocks(self):
        if nx.is_empty(self.graph) == True:
            self.construct_graph()
        nodes = list(nx.nodes(self.graph))
        nodes_2 = []
        blocks=[]
        name_and_block=[]
        block_name_and_block_int=[]
        block_name_and_int=[]
        i=0
        while nodes != []:
            start = nodes[0]
            nodes, for_block = self.double_forward_block(start, nodes)
            nodes, rev_block = self.double_reverse_block(start, nodes)
            block = rev_block+for_block
            if len(block)==1:
                if self.graph.has_node(-block[0]):
                    name=self.int_to_name[block[0]]
                    naming(block, name, name_and_block, block_name_and_block_int, block_name_and_int, '+')
                    name=self.int_to_name[-block[0]]
                    neg_block = [-block[0]]
                    naming(neg_block, name, name_and_block, block_name_and_block_int, block_name_and_int, '-')
                    blocks.append(block)
                else:
                    nodes_2.append(block[0])
            else:    
                name = "+block_"+str(i)
                naming(block, name, name_and_block, block_name_and_block_int, block_name_and_int, '+')
                name = "-block_"+str(i)
                neg_block = reverse_path(block)
                naming(neg_block, name, name_and_block, block_name_and_block_int, block_name_and_int, '-')             
                i=i+1
                blocks.append(block)
        while nodes_2 != []:
            start = nodes_2[0]
            nodes, for_block = self.forward_block(start, nodes_2)
            nodes, rev_block = self.reverse_block(start, nodes_2)
            block = rev_block+for_block
            blocks.append(block)
            if len(block)==1:
                name=self.int_to_name[block[0]]
                naming(block, name, name_and_block, block_name_and_block_int, block_name_and_int, '+')
            else:    
                name = "+block_"+str(i)
                naming(block, name, name_and_block, block_name_and_block_int, block_name_and_int, '+')         
                i=i+1
        return blocks, name_and_block, block_name_and_block_int, block_name_and_int

    
class BlockSeq(GeneSeq):
    
    def __init__(self, genomes, genome_ids, name_and_int = [], graph = False, name_and_block = [], geneseq=None):
        super().__init__(genomes, genome_ids, name_and_int, graph)
        self.geneseq = geneseq
        self.name_to_block = {}
        self.int_to_block = {}
        for el in name_and_block:
            self.name_to_block[el[0]]=el[1]
            self.int_to_block[self.name_to_int[el[0]]]=el[1]
            
    
    @classmethod
    def from_GeneSeq(cls, gs, graph=False):
        if nx.is_empty(gs.graph) == True:
            gs.construct_graph()
        blocks, name_and_block, block_name_and_block_int, block_name_and_int = gs.blocks()
        block_genomes = []
        int_to_block_name = {}
        for el in block_name_and_int:
            int_to_block_name[el[1]] = el[0]
        block_name_to_block_int = {}
        for el in block_name_and_block_int:
            block_name_to_block_int[el[0]]=el[1]
        for i in range(len(gs.genomes)):
            block_genomes.append([])
            j=0
            while j <len(gs.genomes[0]):
                n=1
                for block in blocks:
                    if block[0]==gs.genomes[i][j]:
                        block_genomes[i].append(block_name_to_block_int[int_to_block_name[block[0]]])
                        n = len(block)                    
                    elif -block[-1]==gs.genomes[i][j]:
                        block_genomes[i].append(block_name_to_block_int[int_to_block_name[-block[-1]]])
                        n = len(block)
                j = j+n
        length = max([len(block_genomes[i]) for i in range(len(block_genomes))])
        for i in range(len(block_genomes)):
            n = len(block_genomes[i])
            for j in range(n,length):
                block_genomes[i].append(np.NAN)
        return cls(np.array(block_genomes), gs.genome_ids, block_name_and_block_int, graph, name_and_block, gs)
    
    def expand(self):
        pass
    
    def to_cluster(self, path, cluster):
        genomes = pd.DataFrame(data=self.genomes, index = self.genome_ids)
        name_genomes = genomes.applymap(lambda x: self.int_to_name[x] if not math.isnan(x) else pd.NA)
        genomes.to_csv(path+'/'+cluster+'_blocks_ints.txt', sep="\t", header=False)
        name_genomes.to_csv(path+'/'+cluster+'_blocks_names.txt', sep="\t", header=False)
        maps = pd.DataFrame(data=[self.name_to_int,self.name_to_block])
        maps.T.to_csv(path+'/'+cluster+'_blocks_maps.txt', sep="\t", header=False)
    
if __name__ == '__main__':
    cluster = 'cluster_4'
    path = cluster + '/genomes'
    cluster_4 = GeneSeq.from_cluster( path, cluster)
    #block_cluster_4 = BlockSeq.from_GeneSeq(cluster_4)
    cluster_4.graph_vis()
    #block_cluster_4.to_cluster(path,cluster)
