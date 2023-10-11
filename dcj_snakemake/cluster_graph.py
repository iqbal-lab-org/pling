from commons import produce_full_visualization
from networkx.algorithms.community import *
import itertools
import shutil
import sys
import networkx as nx
import pandas as pd
from pathlib import Path
from glob import glob
import tarfile
import pickle
import logging
logging.basicConfig(stream=sys.stderr,
                    level=logging.DEBUG,
                    format='[%(asctime)s] (%(levelname)s) %(message)s')
logging.getLogger('matplotlib.font_manager').disabled = True


def get_dcj_dists(dcj_dist_matrix): #replace this by reading in dcj dist matrix
    dcj_dists_df = pd.read_csv(dcj_dist_matrix, sep="\t", header=None, index_col=0, skiprows=[0])
    dcj_dists_df = dcj_dists_df.rename(columns={i+1:genome for i, genome in enumerate(list(dcj_dists_df.index))})
    dcj_bool = dcj_dists_df.notna()
    dcj_dists = {}
    for plasmid_1 in dcj_dists_df.index:
        for plasmid_2 in dcj_dists_df.columns:
            if dcj_bool.loc[plasmid_1, plasmid_2]==True:
                    dcj_dists[(plasmid_1, plasmid_2)] = int(dcj_dists_df.loc[plasmid_1, plasmid_2])
    return dcj_dists


def get_communities(gene_jaccard_communities_file):
    logging.info("Getting communities...")
    communities = []
    plasmid_to_community = {}
    with open(gene_jaccard_communities_file) as gene_jaccard_communities_fh:
        for line in gene_jaccard_communities_fh:
            plasmids = line.strip().split()
            for plasmid in plasmids:
                plasmid_to_community[plasmid] = len(communities)
            communities.append(plasmids)
    return communities, plasmid_to_community


def get_node_to_subcommunity(subcommunities):
    node_to_subcommunity = {}
    for subcommunity_index, subcommunity in enumerate(subcommunities):
        for node in subcommunity:
            node_to_subcommunity[node] = subcommunity_index
    return node_to_subcommunity


def fix_small_subcommunities(graph, subcommunities, small_subcommunity_size_threshold):
    # sort subcommunities by size so that we can safely move smaller subcommunities into larger ones
    subcommunities = sorted(subcommunities, key=lambda subcommunity: len(subcommunity))

    for subcommunity_idx, subcommunity in enumerate(subcommunities):
        subcommunity_is_too_small = len(subcommunity)<=small_subcommunity_size_threshold
        if subcommunity_is_too_small:
            subcommunity_neighbours = set(neighbor for node in subcommunity for neighbor in graph.neighbors(node)) - set(subcommunity)
            subcommunity_has_no_neighbors = len(subcommunity_neighbours)==0
            if subcommunity_has_no_neighbors:
                continue

            node_to_subcommunity = get_node_to_subcommunity(subcommunities)
            candidate_subcommunities = set(node_to_subcommunity[neighbor] for neighbor in subcommunity_neighbours)
            largest_subcommunity_size = max(len(subcommunities[subcommunity]) for subcommunity in candidate_subcommunities)
            small_will_be_integrated_into_large = largest_subcommunity_size >= len(subcommunity)
            if small_will_be_integrated_into_large:
                largest_subcommunity_idx = next(filter(lambda idx: len(subcommunities[idx])==largest_subcommunity_size, candidate_subcommunities))
                subcommunities[largest_subcommunity_idx].update(subcommunity)
                subcommunities[subcommunity_idx]=set()
    return subcommunities

def get_subcommunities(graphs, small_subcommunity_size_threshold):
    logging.info(f"Getting subcommunities...")
    subcommunities = []
    plasmid_to_subcommunity = {}
    for graph_index, graph in enumerate(graphs):
        local_subcommunities = list(asyn_lpa_communities(G=graph, weight='weight', seed=42))
        local_subcommunities = fix_small_subcommunities(graph, local_subcommunities,
                                                  small_subcommunity_size_threshold=small_subcommunity_size_threshold)

        node_to_subcommunity = get_node_to_subcommunity(local_subcommunities)
        subcommunities.append(local_subcommunities)
        plasmid_to_subcommunity.update(node_to_subcommunity)
    return subcommunities, plasmid_to_subcommunity



def get_graphs(communities, dcj_dists, dcj_dist_threshold):
    logging.info(f"Getting graphs from communities...")
    graphs = []
    for community_index, community in enumerate(communities):
        graph = nx.Graph()
        graph.add_nodes_from(community)
        for plasmid_1, plasmid_2 in itertools.combinations(community, 2):
            if (plasmid_1, plasmid_2) in dcj_dists:
                dist = dcj_dists[(plasmid_1, plasmid_2)]
                if dist<=dcj_dist_threshold:
                    graph.add_edge(plasmid_1, plasmid_2, dist=dcj_dists[(plasmid_1, plasmid_2)], color="black")
        graphs.append(graph)
    return graphs


def produce_plasmid_communities_tsv(plasmid_communities_filepath, plasmid_to_community, plasmid_to_subcommunity):
    plasmids = []
    communities = []
    blackhole_index = 0
    for plasmid in plasmid_to_community:
        plasmids.append(plasmid)
        if plasmid in plasmid_to_subcommunity:
            subcommunity = plasmid_to_subcommunity[plasmid]
        else:
            subcommunity = f"blackhole_{blackhole_index}"
            blackhole_index += 1
        description = f"community_{plasmid_to_community[plasmid]}_subcommunity_{subcommunity}"
        communities.append(description)

    plasmid_communities_df = pd.DataFrame(data={
        "plasmid": plasmids,
        "community": communities,
    })
    plasmid_communities_df.to_csv(plasmid_communities_filepath, sep="\t", index=False)


def get_blackhole_plasmids(graphs, blackhole_connectivity_threshold, edge_density):
    logging.info("Getting blackhole plasmids...")
    blackhole_plasmids = []
    for graph in graphs:
        blackhole_plasmids_in_graph = []
        for node in graph.nodes:
            if graph.degree(node)>=blackhole_connectivity_threshold:
                neighbors = list(graph.neighbors(node))
                subgraph = nx.induced_subgraph(graph, neighbors)
                nb_of_edges_between_neighbours = subgraph.number_of_edges()
                max_nb_of_edges_between_neighbours = (len(neighbors) * (len(neighbors)-1))//2
                edge_rate = nb_of_edges_between_neighbours/max_nb_of_edges_between_neighbours
                if edge_rate <= edge_density:
                    blackhole_plasmids_in_graph.append(node)
                    logging.debug(f"{node} is a blackhole plasmid, REMOVED")
                else:
                    logging.debug(f"{node} is highly connected but does not connect unrelated plasmids, not removed")
        blackhole_plasmids.append(blackhole_plasmids_in_graph)
    return blackhole_plasmids


def remove_plasmids(graphs, plasmids_for_each_graph):
    for graph, plasmids in zip(graphs, plasmids_for_each_graph):
        graph.remove_nodes_from(plasmids)


def save_data_to_disk(filepath, data):
    with open(filepath, "wb") as fout:
        pickle.dump(data, fout)


def create_outdirs(root_outdir):
    outdir = Path(root_outdir)

    if outdir.exists():
        logging.info("Removing previous visualisations...")
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)

    visualisation_outdir = outdir / "visualisation"
    visualisation_outdir.mkdir()
    misc_outdir = outdir / "misc"
    misc_outdir.mkdir()
    return visualisation_outdir, misc_outdir



jaccard = snakemake.input.communities
dcj_dists_matrix = snakemake.input.matrix
outdir = snakemake.output.dcj_graph_outdir
dcj_threshold =snakemake.params.dcj_dist_threshold
bh_connectivity = snakemake.params.bh_connectivity
bh_neighbours_edge_density = snakemake.params.bh_neighbours_edge_density
small_subcommunity_size_threshold = snakemake.params.small_subcommunity_size_threshold

visualisation_outdir, misc_outdir = create_outdirs(outdir)

dcj_dists = get_dcj_dists(dcj_dists_matrix)

communities, plasmid_to_community = get_communities(jaccard)

graphs = get_graphs(communities, dcj_dists, int(dcj_threshold))
graphs_backup = [graph.copy() for graph in graphs]

blackhole_plasmids = get_blackhole_plasmids(graphs, blackhole_connectivity_threshold=bh_connectivity, edge_density=bh_neighbours_edge_density)
remove_plasmids(graphs, blackhole_plasmids)

subcommunities, plasmid_to_subcommunity = get_subcommunities(graphs, small_subcommunity_size_threshold)

produce_full_visualization(graphs_backup, plasmid_to_subcommunity, visualisation_outdir, blackhole_plasmids,
                           show_blackholes_filter=True)

produce_plasmid_communities_tsv(visualisation_outdir / "plasmid_communities.tsv", plasmid_to_community, plasmid_to_subcommunity)

save_data_to_disk(misc_outdir/"state.pickle", [graphs_backup, plasmid_to_subcommunity, blackhole_plasmids])

logging.info("All done!")
