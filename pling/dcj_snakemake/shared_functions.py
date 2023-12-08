def get_plasmid_to_community(communitypath):
    plasmid_to_community = {}
    with open(communitypath) as communities_fh:
        for community_index, line in enumerate(communities_fh):
            plasmids = line.strip().split()
            for plasmid in plasmids:
                plasmid_to_community[plasmid] = community_index
    return plasmid_to_community

def get_jaccard_distances_for_batch(jaccard_tsv):
    jaccards = {}
    with open(jaccard_tsv, "r") as f:
        for line in f:
            plasmid_1, plasmid_2, jaccard = line.strip().split("\t")
            jaccard = float(jaccard)
            jaccards[(plasmid_1,plasmid_2)] = jaccard
    return jaccards

def get_unimog(outputpath, integerisation, plasmid_to_community, batch, genome1, genome2):
    unimog = ""
    if integerisation == "align":
        unimog = f"{outputpath}/unimogs/batch_{batch}_align.unimog"
    elif integerisation == "anno":
        community = plasmid_to_community[genome1]
        unimog = f"{outputpath}/unimogs/relabelled/blocks/{community}_blocks.unimog"
    return unimog

def get_entries(integerisation, genome_1, genome_2):
    if integerisation == "align":
        return f"{genome_1}~{genome_2}:{genome_1}", f"{genome_1}~{genome_2}:{genome_2}"
    elif integerisation == "anno":
        return genome_1, genome_2
