import os
import pandas as pd
import matplotlib.pyplot as plt

# Load all cluster results from the output directory
output_dir = "output"
tsv_files = [f for f in os.listdir(output_dir) if f.endswith(".tsv")]

# Initialize a dataframe to collect group and cluster data
all_data = []

for file in tsv_files:
    group_name = file.replace("results_", "").replace(".tsv", "")
    df = pd.read_csv(os.path.join(output_dir, file), sep="\t")
    df["group"] = group_name
    all_data.append(df)

# Combine all into one dataframe
df_all = pd.concat(all_data)

# Group by and count cluster assignments per group
summary = df_all.groupby(["group", "cluster_id"]).size().unstack(fill_value=0)

# Plot
summary.plot(kind="bar", stacked=True, colormap="tab10", figsize=(10, 6))
plt.title("Cluster Distribution per Group")
plt.xlabel("Sample Group")
plt.ylabel("Number of Sequences")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("output/cluster_distribution_per_group.png")
plt.show()
