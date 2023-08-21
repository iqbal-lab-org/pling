#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 14:33:22 2023

@author: daria
"""

import pandas as pd
import gffutils
import os

gff_dir = "/home/daria/Documents/projects/test_ggcaller/test1/ggCaller_output/GFF"
csv = pd.read_csv("/home/daria/Documents/projects/test_ggcaller/test1/ggCaller_output/gene_presence_absence.csv")
FASTAFILES = [os.path.basename(el[0]) for el in pd.read_csv("/home/daria/Documents/projects/pling/tests/test1/fastas/list.txt", header=None).values]
FASTAEXT = {os.path.splitext(el)[0]:os.path.splitext(el)[1] for el in FASTAFILES}
genomes = list(FASTAEXT.keys())
db = {}
for genome in genomes:
    db[genome] = gffutils.create_db(f"{gff_dir}/{genome}.gff",  ":memory:")

gene_locations = {}
sorted_genomes = {}
#go through gene_presence_absence.csv columnwise (i.e. genomewise) and for each entry get gene name from "Gene" column, construct gff ID from ggcaller ID, and get location from gff file
for genome in genomes:
    presence_absence = csv[genome]
    gene_locations[genome] = pd.DataFrame()
    for i in range(len(presence_absence)):
        name = csv.loc[i, "Gene"]
        start = float('inf')
        for id in presence_absence[i].split(";"):
            number = id.split("_")[2]
            gff_id = f"{genome}_{number}"
            strand = db[genome][gff_id].strand
            start = min(start, db[genome][gff_id].start)
            gene = strand+name
        info = pd.DataFrame.from_dict({"Gene": [gene], "Start": [start]})
        gene_locations[genome] = pd.concat([gene_locations[genome], info], ignore_index=True)
    gene_locations[genome] = gene_locations[genome].sort_values(by='Start')
    sequence = gene_locations[genome]['Gene'].tolist()
    sorted_genomes[genome]=sequence
    print(sorted_genomes[genome])
