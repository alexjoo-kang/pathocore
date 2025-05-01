import pandas as pd
import os

# Example mapping from sequence ID to scientific species name
taxonomy_mapping = {
    "Insulin_pre": "Homo sapiens",
    "MBP_tagged": "Escherichia coli",
    "14_3_3_protein": "Mus musculus",
    "Control_repeat": "Arabidopsis thaliana",
    "Calmodulin": "Danio rerio",

    "CHD8_variant": "Homo sapiens",
    "SHANK3_fragment": "Macaca mulatta",
    "SCN2A_short": "Rattus norvegicus",
    "FOXP1_fragment": "Bos taurus",
    "CNTNAP2_repeat": "Gallus gallus",

    "TP53_Nterm": "Homo sapiens",
    "KRAS_tagged": "Mus musculus",
    "PIK3CA_fragment": "Rattus norvegicus",
    "BRAF_fragment": "Danio rerio",
    "EGFR_like": "Canis lupus familiaris",

    "Keratin_repeat": "Sus scrofa",
    "Control_fibrous": "Oryza sativa",
    "Histone_like": "Saccharomyces cerevisiae",
    "Hemoglobin_like": "Pan troglodytes",
    "HSP90_motif": "Drosophila melanogaster"
}

# Convert dictionary to a pandas DataFrame
df_taxonomy = pd.DataFrame(list(taxonomy_mapping.items()), columns=["sequence_id", "taxonomic_label"])

# Ensure the 'data' directory exists
os.makedirs("data", exist_ok=True)

# Save the taxonomy mapping to a TSV file
output_path = "data/kraken_labels.tsv"
df_taxonomy.to_csv(output_path, sep="\t", index=False)

print(f"✅ Taxonomic label file saved to: {output_path}")
