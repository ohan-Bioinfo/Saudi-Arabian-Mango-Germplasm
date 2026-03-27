import allel
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# === Load VCF ===
vcf_path = "mango_LD_final.vcf"
callset = allel.read_vcf(vcf_path, fields=["samples", "calldata/GT"], alt_number=1)

samples = callset["samples"]
genotypes = allel.GenotypeArray(callset["calldata/GT"])

# === Build haploid consensus per sample (use first allele for simplicity) ===
fasta_records = []

for i, sample in enumerate(samples):
    # Extract genotype: shape (num_variants, 2), pick 1st allele (0 or 1), handle missing
    haploid = genotypes[:, i, 0]  # use first allele
    seq = "".join(['N' if g < 0 else str(g) for g in haploid])  # convert to "0", "1", or "N"
    # Translate 0/1 to reference/alt: simple for SNP trees
    seq = seq.replace('0', 'A').replace('1', 'G')  # arbitrarily use A/G for binary encoding
    record = SeqRecord(Seq(seq), id=sample, description="")
    fasta_records.append(record)

# === Write FASTA ===
SeqIO.write(fasta_records, "mango_snps_alignment.fasta", "fasta")
print("✅ FASTA written: mango_snps_alignment.fasta")
