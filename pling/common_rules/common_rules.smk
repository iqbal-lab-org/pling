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

def get_not_pairs_jaccard_file():
    if config.get("sourmash", False):
        return f"{OUTPUTPATH}/tmp_files/jaccard_pairwise/not_pairs_jaccard_distance.tsv"
    else:
        return ""

rule cat_jaccard:
    input:
        jaccards = expand(f"{OUTPUTPATH}/tmp_files/jaccard_pairwise/batch_{{batch}}_jaccard.tsv", batch=[str(i) for i in range(get_number_of_batches(OUTPUTPATH))]),
        not_pairs = get_not_pairs_jaccard_file()
    output:
        all_jaccard_distances = f"{OUTPUTPATH}/jaccard/all_pairs_jaccard_distance.tsv"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000*attempt
    shell:
        """
        cat <(echo -e "plasmid_1\tplasmid_2\tdistance") {input.jaccards} {input.not_pairs}> {output.all_jaccard_distances}
        """

localrules: cat_jaccard
