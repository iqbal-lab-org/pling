import networkx as nx
from networkx.algorithms.community import *
from networkx.algorithms.components import connected_components

testing = False

# for testing
if testing:
    gene_jaccard_filepath = "../test_data/small_gene_jaccard.tsv"
    communities_filepath = "../test_data/output/communities.txt"
    communities_sizes_filepath = "../test_data/output/communities_sizes.txt"
    gene_jaccard_threshold = 0.4
else:
    gene_jaccard_filepath = snakemake.input.gene_jaccard
    communities_filepath = snakemake.output.communities
    communities_sizes_filepath = snakemake.output.communities_sizes
    gene_jaccard_threshold = snakemake.params.gene_jaccard_threshold



def get_graph(gene_jaccard_filepath, gene_jaccard_threshold):
    graph = nx.Graph()
    with open(gene_jaccard_filepath) as gene_jaccard_fh:
        for line in gene_jaccard_fh:
            from_plasmid, to_plasmid, gene_jaccard = line.strip().split("\t")
            gene_jaccard = float(gene_jaccard)

            # TODO: is there a bug here? what if a plasmid is totally isolated (i.e. its gene jaccard to all others is 0, it won't be added to the graph!)
            if gene_jaccard >= gene_jaccard_threshold:
                graph.add_node(from_plasmid)
                graph.add_node(to_plasmid)
                graph.add_edge(from_plasmid, to_plasmid, weight=gene_jaccard)
    return graph



def main():
    print(f"Building graph...")
    graph = get_graph(gene_jaccard_filepath, gene_jaccard_threshold)
    print(f"Building graph - done!")

    print(f"Finding communities...")
    # communities = list(asyn_lpa_communities(G=graph, weight='weight', seed=42))
    communities = list(connected_components(graph))
    print(f"Finding communities - done!")

    # outputs communities
    with open(communities_filepath, "w") as communities_fh, open(communities_sizes_filepath, "w") as communities_sizes_fh:
        for community in communities:
            print(" ".join(community), file=communities_fh)
            print(f"{len(community)}", file=communities_sizes_fh)

main()
