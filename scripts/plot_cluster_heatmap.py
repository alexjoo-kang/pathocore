import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Load all TSV cluster results
output_dir = "output"
tsv_files = [f for f in os.listdir(output_dir) if f.startswith("results_") and f.endswith(".tsv")]

# Combine all into one DataFrame
all_data = []
for file in tsv_files:
    group = file.replace("results_", "").replace(".tsv", "")
    df = pd.read_csv(os.path.join(output_dir, file), sep="\t")
    df["group"] = group
    all_data.append(df)

df_all = pd.concat(all_data)

# Create a cross-tab of group vs cluster_id
cluster_counts = pd.crosstab(df_all["group"], df_all["cluster_id"])

# Normalize row-wise to show proportion
cluster_props = cluster_counts.div(cluster_counts.sum(axis=1), axis=0)

# Plot as heatmap
plt.figure(figsize=(10, 6))
sns.heatmap(cluster_props, annot=True, cmap="viridis", fmt=".2f", cbar_kws={"label": "Proportion"})
plt.title("Proportional Cluster Distribution per Group")
plt.ylabel("Sample Group")
plt.xlabel("Cluster ID")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "cluster_proportion_heatmap.png"))
plt.show()
