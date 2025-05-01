import os
from collections import Counter, defaultdict
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Define directory where cluster FASTA files are located
cluster_dir = "output/clusters_fasta"
output_path = "output/cluster_composition_summary.tsv"

# Define standard amino acids
standard_aa = list("ACDEFGHIKLMNPQRSTVWY")

# Store composition per cluster
composition_data = defaultdict(dict)

# Analyze each cluster file
for file in sorted(os.listdir(cluster_dir)):
    if file.endswith(".fasta"):
        cluster_id = file.replace("cluster_", "").replace(".fasta", "")
        total_counts = Counter()
        total_length = 0

        for record in SeqIO.parse(os.path.join(cluster_dir, file), "fasta"):
            aa_counts = Counter(str(record.seq))
            total_counts.update(aa_counts)
            total_length += len(record.seq)

        # Normalize to frequency
        for aa in standard_aa:
            freq = total_counts[aa] / total_length if total_length > 0 else 0
            composition_data[cluster_id][aa] = round(freq, 4)

# Convert to DataFrame
composition_df = pd.DataFrame.from_dict(composition_data, orient="index")
composition_df.index.name = "cluster_id"
composition_df = composition_df[standard_aa]  # Order columns
composition_df.to_csv(output_path, sep="\t")

print(f"\n✅ Amino acid composition per cluster saved to: {output_path}")

# Plot heatmap
plt.figure(figsize=(12, 6))
sns.heatmap(composition_df, annot=True, cmap="coolwarm", fmt=".2f", cbar_kws={"label": "Frequency"})
plt.title("Amino Acid Composition per Cluster")
plt.ylabel("Cluster ID")
plt.xlabel("Amino Acid")
plt.tight_layout()
plt.savefig("output/cluster_composition_heatmap.png")
plt.show()
