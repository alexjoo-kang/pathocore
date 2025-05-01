import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load functional summary file
input_path = "output/functional_summary.tsv"
df = pd.read_csv(input_path, sep="\t")

# Pivot the data: rows = function, columns = cluster_id, values = counts
pivot_df = df.pivot(index="function", columns="cluster_id", values="count").fillna(0)

# Plot heatmap
plt.figure(figsize=(10, 6))
sns.heatmap(pivot_df, annot=True, fmt=".0f", cmap="YlGnBu", cbar_kws={"label": "Count"})
plt.title("Function Distribution Across Clusters")
plt.xlabel("Cluster ID")
plt.ylabel("Function")
plt.tight_layout()

# Save and show
plot_path = "output/functional_distribution_heatmap.png"
plt.savefig(plot_path)
plt.show()

print(f"✅ Functional distribution heatmap saved to: {plot_path}")
