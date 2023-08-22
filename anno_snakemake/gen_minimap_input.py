from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path
import os

seq_out = dict()
for file in os.listdir(snakemake.input.align_dir):
    file_name = str(Path(file))
    gene = str(Path(file).with_suffix(""))
    gene = gene.replace(".aln", "")
    seq_in = SeqIO.parse(f"{snakemake.input.align_dir}/{file_name}", "fasta")
    for record in seq_in:
        genome, prokka_gene = record.id.split(';')
        blah = (str(record.seq)).replace('-','')
        seq = Seq(blah)
        if genome not in seq_out:
            seq_out[genome] = [SeqRecord(seq, gene, "", "")]
        else:
            (seq_out[genome]).append(SeqRecord(seq, gene, "", ""))
for genome in seq_out:
    SeqIO.write(seq_out[genome], f"{snakemake.params.out_dir}/{genome}.fna", "fasta")
