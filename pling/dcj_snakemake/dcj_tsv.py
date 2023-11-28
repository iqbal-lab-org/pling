import pandas as pd
import os

dist=0

genome1s = []
genome2s = []
distances = []

for genome1 in snakemake.params.genomes:
    for genome2 in snakemake.params.genomes:
        if genome1 == genome2:
            dist = 0
        else:
            file1 =  f"{snakemake.params.outputpath}/tmp_files/dists_pairwise/{genome1}~{genome2}.dist"
            file2 = f"{snakemake.params.outputpath}/tmp_files/dists_pairwise/{genome2}~{genome1}.dist"
            if os.path.exists(file1)==True and os.path.exists(file2)==False:
                file = file1
                f = open(file, 'r')
                line = f.readline()
                dist = int(line.strip().split(" ")[2])
                f.close()
            elif os.path.exists(file2)==True and os.path.exists(file1)==False:
                file = file2
                f = open(file, 'r')
                line = f.readline()
                dist = int(line.strip().split(" ")[2])
                f.close()
            else:
                dist = pd.NA

        genome1s.append(genome1)
        genome2s.append(genome2)
        distances.append(dist)


df = pd.DataFrame({"plasmid_1": genome1s, "plasmid_2": genome2s, "distance": distances})
df.to_csv(snakemake.output.tsv, sep="\t")
