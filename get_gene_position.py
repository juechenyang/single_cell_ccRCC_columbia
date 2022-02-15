from gtfparse import read_gtf

df = read_gtf("gencode.v37.annotation.gtf")
df_genes = df[(df["feature"] == "gene") & (df["gene_type"] == "protein_coding")]
yes = df_genes[["gene_name", "seqname", "start", "end"]]
yes.to_csv("human_gene_pos.csv", index=False)