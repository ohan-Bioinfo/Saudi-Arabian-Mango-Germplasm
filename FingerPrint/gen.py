import pandas as pd

# === Load data ===
genotypes = pd.read_csv("top150_mango_genotypes.raw", sep='\s+')
bim = pd.read_csv("mango_QC_v5_LD_final.bim", sep='\t', header=None,
                  names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
top_snps = pd.read_csv("top150_fingerprint_snps.csv")

# === Prepare genotype data ===
sample_info = genotypes[["FID", "IID"]]
geno_only = genotypes.drop(columns=["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"])

# Normalize SNP IDs from .raw (e.g., 'chr1:12345_A' -> 'chr1:12345')
normalized_ids = geno_only.columns.str.extract(r"(chr[\w\d]+:\d+)")[0].values
snp_col_map = dict(zip(geno_only.columns, normalized_ids))

# Extract allele info for top SNPs
bim_subset = bim[bim["SNP"].isin(top_snps["SNP"])].set_index("SNP")[["A1", "A2"]]

# Map alleles back to genotype columns
def convert_genotype(col):
    snp = snp_col_map[col.name]
    if snp not in bim_subset.index:
        return col.map(lambda x: "NA")
    a1, a2 = bim_subset.loc[snp]
    return col.map({0: f"{a1}/{a1}", 1: f"{a1}/{a2}", 2: f"{a2}/{a2}"}).fillna("NA")

# Apply conversion
genotype_alleles = geno_only.apply(convert_genotype, axis=0)

# Rename columns to SNP1, SNP2, ...
genotype_alleles.columns = [f"SNP{i+1}" for i in range(genotype_alleles.shape[1])]

# Add sample info
genotype_alleles.insert(0, "IID", sample_info["IID"])
genotype_alleles.insert(0, "FID", sample_info["FID"])

# Save to CSV
genotype_alleles.to_csv("top150_mango_genotypes_alleles.csv", index=False)
print("✅ Final allele file saved as 'top150_mango_genotypes_alleles.csv'")
