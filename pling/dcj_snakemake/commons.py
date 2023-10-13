import collections
import shutil
import json
import networkx as nx
import logging
from collections import defaultdict
import random
from pathlib import Path
RNG = random.Random(42)


colours = ["darkred", "indianred", "fuchsia", "hotpink", "darkgreen", "chartreuse", "darkblue", "blue",
           "aqua", "beige", "darkviolet", "mediumpurple", "dimgrey", "darkorange", "maroon", "yellow", "orange",
           "teal"]
basic_colours = ["red", "magenta", "green", "blue", "black", "orange", "maroon", "cyan"]
RNG.shuffle(colours)


def get_script_parent_dir():
    return Path(__file__).parent


def get_templates_dir():
    return get_script_parent_dir() / "templates"


def get_libs_dir():
    return get_script_parent_dir() / "libs"


def fix_node_to_subcommunity_attributes(graph, node_to_subcommunity, blackhole_plasmids):
    for node, attrs in graph.nodes.items():
        if node in node_to_subcommunity:
            attrs["color"] = colours[node_to_subcommunity[node]%len(colours)]
        else:
            attrs["color"] = "black"

        if node in blackhole_plasmids:
            attrs["shape"] = "star"
            attrs["is_blackhole"] = True
        else:
            attrs["is_blackhole"] = False


def read_template(template_filepath):
    with open(template_filepath) as template_fh:
        template_src = template_fh.readlines()
    return map(lambda line: line.strip("\n"), template_src)


def get_plasmid_induced_components(graph, plasmids):
    subgraph = graph.subgraph(plasmids)
    return nx.connected_components(subgraph)


def get_simulation_time(graph):
    if graph.number_of_nodes()<=5 or graph.number_of_edges() <= 10:
        return 1000
    else:
        return 10000

def produce_visualisation(graph, outdir, label, visualisation_count, nb_of_black_holes, show_blackholes_filter=False, show_samples_filter=False,
                          sample_to_plasmids=None):
    outdir = outdir / str(visualisation_count//4000)
    outdir.mkdir(parents=True, exist_ok=True)

    visualisation_src = read_template(get_templates_dir() / "visualisation_template.html")

    graph_as_cy_dict = nx.cytoscape_data(graph)
    elements_as_cy_json = json.dumps(graph_as_cy_dict["elements"])

    if sample_to_plasmids:
        samples_selectors = []
        for sample, plasmids in sample_to_plasmids.items():
            induced_components = get_plasmid_induced_components(graph, plasmids)
            for component_index, component in enumerate(induced_components):
                node_selector = [f"node#{node}" for node in component]
                node_selector = ", ".join(node_selector)
                if component_index==0:
                    samples_selectors.append(f"sample_selector_nodes['{sample}'] = [];")
                samples_selectors.append(f"sample_selector_nodes['{sample}'].push(cy.elements('{node_selector}'));")
        samples_selectors_str = "\n".join(samples_selectors)

        sample_hits_checkboxes = []
        for sample_index, sample in enumerate(sample_to_plasmids):
            colour = basic_colours[sample_index % len(basic_colours)]
            sample_hits_checkboxes.append(
                f'<input type="checkbox" id="{sample}" name="{sample}" onclick="show_sample_hits(\'{sample}\', \'{colour}\')">'
                f'<label for="{sample}">{sample} ({len(sample_to_plasmids[sample])} hits) <span style="color:{colour}">&#9632;</span></label><br/>'
            )
    else:
        samples_selectors_str = ""
        sample_hits_checkboxes = []

    filters = []
    if show_blackholes_filter:
        filters.append(f'<label for="hide_blackholes">Hide blackhole plasmids ({nb_of_black_holes} present)</label>'
                       f'<input type="checkbox" id="hide_blackholes" name="hide_blackholes"><br/>')
    if show_samples_filter:
        filters.append("Show hits for samples:<br/>")
        filters.extend(sample_hits_checkboxes)

    custom_buttons = []
    if show_blackholes_filter:
        custom_buttons.append('<div><input type="submit" value="Redraw" onclick="redraw()"></div>')

    with open(outdir / f"{label}.html", "w") as visualisation_fh:
        for line in visualisation_src:
            line = line.replace("<samples_selectors>", samples_selectors_str)
            line = line.replace("<elements_tag>", elements_as_cy_json)
            line = line.replace("<movementThreshold>", str(graph.number_of_edges()))
            line = line.replace("<maxSimulationTime>", str(get_simulation_time(graph)))
            line = line.replace("<filters_tag>", "\n".join(filters))
            line = line.replace("<custom_buttons_tag>", "\n".join(custom_buttons))
            print(line, file=visualisation_fh)


def produce_index_file(outdir, graphs, graphs_descriptions, objects_description, blackhole_plasmids, graph_to_sample_to_plasmids=None):
    nb_of_elems_to_graph_indexes = collections.defaultdict(list)
    if graph_to_sample_to_plasmids is None:
        # sort by edges
        for graph_index, graph in enumerate(graphs):
            nb_of_elems_to_graph_indexes[graph.number_of_edges()].append(graph_index)
    else:
        # sort by number of sample hits
        for graph_index, sample_to_plasmids in graph_to_sample_to_plasmids.items():
            nb_of_sample_hits = len(sample_to_plasmids)
            nb_of_elems_to_graph_indexes[nb_of_sample_hits].append(graph_index)

    index_src = read_template(get_templates_dir() / "index_template.html")
    visualisation_links=[]
    for nb_of_elems in sorted(nb_of_elems_to_graph_indexes.keys(), reverse=True):
        for graph_index in nb_of_elems_to_graph_indexes[nb_of_elems]:
            graph = graphs[graph_index]
            blackhole_plasmids_for_graph = blackhole_plasmids[graph_index]
            if len(blackhole_plasmids_for_graph)>0:
                warning=" - WARNING: BLACKHOLE SPOTTED!"
            else:
                warning=""
            description = f'View {objects_description} {graphs_descriptions[graph_index]} '
            if graph_to_sample_to_plasmids is not None:
                description += f"- {len(graph_to_sample_to_plasmids[graph_index])} samples hit "
            description += f'({graph.number_of_edges()} edges, {graph.number_of_nodes()} nodes){warning}</a><br/>'
            visualisation_links.append(f'<a href="graphs/{graph_index//4000}/{graphs_descriptions[graph_index]}.html" target="_blank">{description}')

    with open(outdir / f"index.html", "w") as index_fh:
        for line in index_src:
            if graph_to_sample_to_plasmids is None:
                line = line.replace("<header_message>", f"<h3>Largest {objects_description} are shown first</h3>")
            else:
                line = line.replace("<header_message>", f"<h3>{objects_description} with more sample hits are shown first</h3>")
            line = line.replace("<communities_links_tag>", "\n".join(visualisation_links))
            line = line.replace("<objects_description>", objects_description)
            print(line, file=index_fh)

    shutil.copytree(get_libs_dir(), outdir / "libs")


def get_subgraphs(graph, plasmid_to_subcommunity, use_subgraphs):
    if use_subgraphs:
        subgraphs = {comp_index: graph.subgraph(component).copy() for comp_index, component in enumerate(nx.connected_components(graph))}
    else:
        components = defaultdict(list)
        for plasmid in graph:
            subcommunity = plasmid_to_subcommunity[plasmid]
            components[subcommunity].append(plasmid)
        subgraphs = {comp_index: graph.subgraph(component).copy() for comp_index, component in components.items()}

    return subgraphs


def produce_full_visualization(graphs, plasmid_to_subcommunity, visualisation_outdir, blackhole_plasmids,
                               use_subgraphs=False,
                               use_subcommunities=False,
                               show_blackholes_filter=False, show_samples_filter=False, graph_to_sample_to_plasmids=None):
    all_subgraphs = []
    subgraphs_blackhole_plasmids = []
    subgraph_to_sample_to_plasmids = defaultdict(lambda: defaultdict(set))
    descriptions = []

    visualisation_count=0
    for graph_index, graph in enumerate(graphs):
        logging.info(f"Producing graph {graph_index}/{len(graphs)}")
        blackhole_plasmids_for_graph = set(blackhole_plasmids[graph_index] if blackhole_plasmids is not None else [])
        fix_node_to_subcommunity_attributes(graph, plasmid_to_subcommunity, blackhole_plasmids_for_graph)
        sample_to_plasmids = graph_to_sample_to_plasmids.get(graph_index) if graph_to_sample_to_plasmids else None

        if use_subgraphs or use_subcommunities:

            # TODO: I am not it is a good idea to remove blackhole plasmids here, for now let's keep this...
            graph.remove_nodes_from(blackhole_plasmids_for_graph)

            subgraphs = get_subgraphs(graph, plasmid_to_subcommunity, use_subgraphs)

            for subgraph_index, subgraph in subgraphs.items():
                # populates all subgraphs variables
                all_subgraphs.append(subgraph)

                blackhole_plasmids_for_subgraph = blackhole_plasmids_for_graph.intersection(list(subgraph))
                subgraphs_blackhole_plasmids.append(blackhole_plasmids_for_subgraph)

                subgraph_sample_to_plasmids = None
                all_subgraphs_index = len(all_subgraphs)-1
                if sample_to_plasmids is not None:
                    subgraph_sample_to_plasmids = defaultdict(set)
                    for sample, plasmids in sample_to_plasmids.items():
                        for plasmid in plasmids:
                            if plasmid in subgraph:
                                subgraph_sample_to_plasmids[sample].add(plasmid)
                    subgraph_to_sample_to_plasmids[all_subgraphs_index] = subgraph_sample_to_plasmids

                description = f"graph_{graph_index}_subgraph_{subgraph_index}"
                descriptions.append(description)

                produce_visualisation(subgraph, visualisation_outdir / "graphs", description, visualisation_count,
                                      nb_of_black_holes=len(blackhole_plasmids_for_subgraph),
                                      show_blackholes_filter=show_blackholes_filter,
                                      show_samples_filter=show_samples_filter,
                                      sample_to_plasmids=subgraph_sample_to_plasmids)
                visualisation_count+=1
        else:  # normal graphs (no subgraph/subcommunities)
            description = f"graph_{graph_index}"
            descriptions.append(description)
            produce_visualisation(graph, visualisation_outdir / "graphs", description, visualisation_count,
                              nb_of_black_holes=len(blackhole_plasmids_for_graph),
                              show_blackholes_filter=show_blackholes_filter, show_samples_filter=show_samples_filter,
                              sample_to_plasmids=sample_to_plasmids)
            visualisation_count += 1

    logging.info("Producing index file...")
    if use_subgraphs or use_subcommunities:
        object_description = "subgraphs" if use_subgraphs else "subcommunity"
        produce_index_file(visualisation_outdir, all_subgraphs, descriptions, object_description, subgraphs_blackhole_plasmids, subgraph_to_sample_to_plasmids)
    else:
        produce_index_file(visualisation_outdir, graphs, descriptions, "community", blackhole_plasmids, graph_to_sample_to_plasmids)
    logging.info("Producing index file - done!")
