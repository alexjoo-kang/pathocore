import os
import pandas as pd
from collections import Counter
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns

# Define basic amino acid property groups
AA_PROPERTIES = {
    "Hydrophobic": set("AILMFWV"),
    "Polar": set("NQSTYC"),
    "Positive": set("KRH"),
    "Negative": set("DE"),
    "Special": set("GP")
}

# Prepare a structure to hold counts
cluster_properties = {}

# Path to cluster FASTA files
cluster_fasta_dir = "output/clusters_fasta"
fasta_files = [f for f in os.listdir(cluster_fasta_dir) if f.endswith(".fasta")]

for fasta_file in fasta_files:
    cluster_id = fasta_file.replace("cluster_", "").replace(".fasta", "")
    records = list(SeqIO.parse(os.path.join(cluster_fasta_dir, fasta_file), "fasta"))
    
    # Initialize counters
    property_counter = Counter()
    total_residues = 0
    
    for record in records:
        seq = str(record.seq)
        total_residues += len(seq)
        for aa in seq:
            for prop, aa_set in AA_PROPERTIES.items():
                if aa in aa_set:
                    property_counter[prop] += 1
                    break

    # Normalize by total residues
    property_freq = {prop: count / total_residues for prop, count in property_counter.items()}
    cluster_properties[cluster_id] = property_freq

# Convert to DataFrame
df_props = pd.DataFrame.from_dict(cluster_properties, orient="index").fillna(0)
df_props.index.name = "Cluster"

# Plot heatmap
plt.figure(figsize=(8, 6))
sns.heatmap(df_props, annot=True, cmap="coolwarm", fmt=".2f", cbar_kws={"label": "Proportion"})
plt.title("Average Amino Acid Property Proportions per Cluster")
plt.ylabel("Cluster")
plt.xlabel("Amino Acid Property")
plt.tight_layout()

# Save plot
plot_path = os.path.join("output", "cluster_aa_property_heatmap.png")
plt.savefig(plot_path)
plt.show()
