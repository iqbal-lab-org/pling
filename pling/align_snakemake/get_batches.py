import argparse
import pandas as pd
import os

def write_batch_file(output_dir, batches, batch):
    with open(f"{output_dir}/batch_{batch}.txt","w") as batch_list:
        batch_list.write("\n".join([str(el) for el in batches[str(batch)]]))

def get_pairs(GENOMES):
    genome_pairs=[]
    n = len(GENOMES)
    for i in range(n):
        j=0
        while j<i:
            genome_pairs.append([GENOMES[i], GENOMES[j]])
            j=j+1
    return genome_pairs

def main():
    parser = argparse.ArgumentParser(description="Generate lists of genome pairs per batch")

    parser.add_argument("--genomes_list")
    parser.add_argument("--batch_size", type=int)
    parser.add_argument("--number_of_batches", type=int)
    parser.add_argument("--outputpath")

    args = parser.parse_args()

    FASTAFILES = [el[0] for el in pd.read_csv(args.genomes_list, header=None).values]
    fastaext = {os.path.splitext(os.path.basename(el))[0]:os.path.splitext(os.path.basename(el))[1] for el in FASTAFILES}
    GENOMES = list(fastaext.keys())

    PAIRS = get_pairs(GENOMES)
    number_of_batches = args.number_of_batches
    batches = {}
    for i in range(number_of_batches-1):
        batches[str(i)] = [PAIRS[j] for j in range(i*args.batch_size, (i+1)*args.batch_size)]
        write_batch_file(f"{args.outputpath}/tmp_files/batches", batches, str(i))
    batches[str(number_of_batches-1)] = [PAIRS[j] for j in range((number_of_batches-1)*args.batch_size, len(PAIRS))]
    write_batch_file(f"{args.outputpath}/tmp_files/batches", batches, str(number_of_batches-1))

if __name__ == "__main__":
    main()
