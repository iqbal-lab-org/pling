# common rules and functions shared by align_snakemake and jac_network_snakemake workflows
# TODO: refactor as this is a bit brittle because it requires the predefinition of the following global variables
# TODO: prior to including this:
# OUTPUTPATH
# GENOMES

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


def get_pairs(GENOMES):
    genome_pairs=[]
    n = len(GENOMES)
    for i in range(n):
        j=0
        while j<i:
            genome_pairs.append([GENOMES[i], GENOMES[j]])
            j=j+1
    return genome_pairs

def get_files(type,OUTPUTPATH, PAIRS):
    files = []
    if type == "jaccard":
        dir = f"{OUTPUTPATH}/tmp_files/jaccard_pairwise"
        end = "jaccard_distance.tsv"
    elif type == "unimog":
        dir = f"{OUTPUTPATH}/unimogs"
        end = "align.unimog"
    for el in PAIRS:
        if len(el)>1:
            files.append(f"{dir}/{el[0]}~{el[1]}_{end}")
    return files


rule cat_jaccard:
    input:
        jaccards = get_files("jaccard", OUTPUTPATH, get_pairs(GENOMES))
    output:
        all_jaccard_distances = f"{OUTPUTPATH}/jaccard/all_pairs_jaccard_distance.tsv"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000*attempt
    shell:
        """
        cat <(echo -e "plasmid_1\tplasmid_2\tdistance") {input.jaccards} > {output.all_jaccard_distances}
        """
localrules: cat_jaccard
