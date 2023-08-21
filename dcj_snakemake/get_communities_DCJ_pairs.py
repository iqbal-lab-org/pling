from pathlib import Path
from collections import defaultdict

testing = False

# for testing
if testing:
    gene_jaccard_filepath = "../test_data/small_gene_jaccard.tsv"
    communities_filepath = "../test_data/output/gene_jaccard_communities.txt"
    dcj_input_dir = "../test_data/output/dcj_input_dir"
    gene_jaccard_threshold = 0.4
else:
    gene_jaccard_filepath = snakemake.input.jaccard
    communities_filepath = snakemake.input.communities
    gene_jaccard_threshold = snakemake.params.jaccard_threshold
    dcj_input_dir = snakemake.output.DCJ_input_dir


def get_plasmid_to_community(communities_filepath):
    plasmid_to_community = {}
    with open(communities_filepath) as communities_fh:
        for community_index, line in enumerate(communities_fh):
            plasmids = line.strip().split()
            for plasmid in plasmids:
                plasmid_to_community[plasmid] = community_index
    return plasmid_to_community


plasmid_to_community = get_plasmid_to_community(communities_filepath)
community_id_to_pairs = defaultdict(list)

with open(gene_jaccard_filepath) as gene_jaccard_fh:
    for line in gene_jaccard_fh:
        plasmid_1, plasmid_2, gene_jaccard = line.strip().split("\t")
        gene_jaccard = float(gene_jaccard)
        if gene_jaccard >= gene_jaccard_threshold and plasmid_1 in plasmid_to_community and plasmid_2 in plasmid_to_community:
            community_1 = plasmid_to_community[plasmid_1]
            community_2 = plasmid_to_community[plasmid_2]
            if community_1 == community_2:
                community_id_to_pairs[community_1].append((plasmid_1, plasmid_2, str(gene_jaccard)))

dcj_input_dir = Path(dcj_input_dir)
dcj_input_dir.mkdir(parents=True)
for community_id, pairs in community_id_to_pairs.items():
    with open(dcj_input_dir / f"DCJ_pairs_community_{community_id}.txt", "w") as dcj_pairs_fh:
        print("plasmid1 plasmid2 jaccard", file=dcj_pairs_fh)
        for pair in pairs:
            print(" ".join(pair), file=dcj_pairs_fh)
