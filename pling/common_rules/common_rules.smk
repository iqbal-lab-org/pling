# common rules and functions shared by align_snakemake and jac_network_snakemake workflows
# TODO: refactor as this is a bit brittle because it requires the predefinition of the following global variables
# TODO: prior to including this:
# OUTPUTPATH
# GENOMES
from pling.utils import get_number_of_batches

rule create_genomes_tsv:
    output:
        genomes_tsv = f"{OUTPUTPATH}/tmp_files/genomes.tsv"
    params:
        genomes = GENOMES
    run:
        with open(output.genomes_tsv, "w") as f:
            f.write("plasmid\n")
            for genome in params.genomes:
                f.write(f"{genome}\n")
localrules: create_genomes_tsv

rule cat_containment:
    input:
        containments = expand(f"{OUTPUTPATH}/tmp_files/containment_batchwise/batch_{{batch}}_containment.tsv", batch=[str(i) for i in range(get_number_of_batches(OUTPUTPATH))])
    output:
        all_containment_distances = f"{OUTPUTPATH}/containment/all_pairs_containment_distance.tsv"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000*attempt
    run:
        with open(output.all_containment_distances, "w") as contain_out:
            contain_out.write("plasmid_1\tplasmid_2\tdistance\n")
            for file in input.containments:
                with open(file, "r") as f:
                    to_cat = f.read()
                contain_out.write(to_cat)

localrules: cat_containment

def get_metadata(metadata):
    if metadata == "None":
        return ""
    else:
        return f"--plasmid-metadata {metadata}"

rule get_communities:
    input:
        containment = f"{OUTPUTPATH}/containment/all_pairs_containment_distance.tsv",
        genomes = rules.create_genomes_tsv.output.genomes_tsv
    output:
        communities = directory(f"{OUTPUTPATH}/containment/containment_communities"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: config["get_communities_mem"]*attempt
    conda: "../envs/plasnet.yaml"
    params:
        containment_distance=config["seq_containment_distance"],
        bh_connectivity=config["bh_connectivity"],
        bh_neighbours_edge_density=config["bh_neighbours_edge_density"],
        metadata = get_metadata(config["metadata"])
    shell: """
        plasnet split \
            --distance-threshold {params.containment_distance} \
            --bh-connectivity {params.bh_connectivity} \
            --bh-neighbours-edge-density {params.bh_neighbours_edge_density} \
            --output-plasmid-graph \
            {params.metadata} \
            {input.genomes} \
            {input.containment} \
            {output.communities}
    """
