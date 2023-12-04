# common rules and functions shared by align_snakemake and jac_network_snakemake workflows
# TODO: refactor as this is a bit brittle because it requires the predefinition of the following global variables
# TODO: prior to including this:
# OUTPUTPATH
# GENOMES
# number_of_batches

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
        jaccards = expand(f"{OUTPUTPATH}/tmp_files/jaccard_pairwise/batch_{{batch}}_jaccard.tsv", batch=[str(i) for i in range(number_of_batches)])
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

rule get_batches:
    output:
        [f"{OUTPUTPATH}/tmp_files/batches/batch_{i}.txt" for i in range(number_of_batches)]
    params:
        outputpath = OUTPUTPATH,
        genomes_list = config["genomes_list"],
        batch_size = batch_size,
        number_of_batches = number_of_batches,
        pling_root_dir = get_pling_root_dir()
    shell:
        """
        PATH="$CONDA_PREFIX"/bin:$PATH
        python {params.pling_root_dir}/pling/align_snakemake/get_batches.py \
            --genomes_list {params.genomes_list} \
            --batch_size {params.batch_size} \
            --number_of_batches {params.number_of_batches} \
            --outputpath {params.outputpath}
        """
