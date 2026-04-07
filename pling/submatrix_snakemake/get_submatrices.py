import pandas as pd
import logging
from pathlib import Path

SUBCOMMUNITIESPATH = snakemake.input.subcom
outdir_1 = snakemake.output.outdir_1
outdir_2 = snakemake.output.outdir_2
Path(outdir_1).mkdir(parents=True, exist_ok=True)
Path(outdir_2).mkdir(parents=True, exist_ok=True)
subcommunities_df = pd.read_csv(SUBCOMMUNITIESPATH, sep='\t')
subcommunities_list = list(subcommunities_df['type'])
count = {subcommunity:subcommunities_list.count(subcommunity) for subcommunity in set(subcommunities_list)}
dcj_tsv = snakemake.input.tsv
dcj_dists = {}
with open(dcj_tsv) as f:
    next(f)
    for line in f:
        plasmid1, plasmid2, dist = line.strip().split("\t")
        dcj_dists[(plasmid1,plasmid2)] = int(dist)
        dcj_dists[(plasmid2,plasmid1)] = int(dist)
submatrices = {}

for subcommunity in count.keys():
    if count[subcommunity]>2:
        plasmids = list(subcommunities_df[subcommunities_df['type']==subcommunity]['plasmid'])
        distances = pd.DataFrame(index = [len(plasmids)] + plasmids, columns = plasmids)
        complete = True
        for plasmid1 in plasmids:
            for plasmid2 in plasmids:
                if (plasmid1,plasmid2) in dcj_dists.keys():
                    distances.loc[plasmid1,plasmid2]=dcj_dists[(plasmid1,plasmid2)]
                elif plasmid1 == plasmid2:
                    distances.loc[plasmid1,plasmid2]=0
                else:
                    distances.loc[plasmid1,plasmid2]=pd.NA
                    logging.warning(f"Plasmid pair {plasmid1}, {plasmid2} does not meet containment threshold")
                    complete = False
        if complete:
            distances.to_csv(f"{outdir_2}/{subcommunity}.dist", sep="\t", header=False)
        else:
            distances.to_csv(f"{outdir_1}/{subcommunity}_incomplete.dist", sep="\t", header=False)
