import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

# Paths
data_dir = "data"
function_file = os.path.join(data_dir, "function_labels.tsv")
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

# Load function label mapping
df_func = pd.read_csv(function_file, sep="\t")
func_map = dict(zip(df_func["sequence_id"], df_func["function"]))

# Collect function counts per group
function_summary = []

for fasta_file in os.listdir(data_dir):
    if fasta_file.endswith(".fasta"):
        group = fasta_file.replace(".fasta", "")
        records = SeqIO.parse(os.path.join(data_dir, fasta_file), "fasta")
        count_dict = {}

        for record in records:
            func = func_map.get(record.id, "Unknown")
            count_dict[func] = count_dict.get(func, 0) + 1

        for func, count in count_dict.items():
            function_summary.append({
                "group": group,
                "function": func,
                "count": count
            })

# Convert to DataFrame
df_summary = pd.DataFrame(function_summary)

# Pivot and normalize
pivot_df = df_summary.pivot(index="group", columns="function", values="count").fillna(0)
pivot_prop = pivot_df.div(pivot_df.sum(axis=1), axis=0)

# Plot
plt.figure(figsize=(10, 8))  # Increase figure height
sns.heatmap(pivot_prop, annot=True, cmap="magma", fmt=".2f", cbar_kws={"label": "Proportion"})

plt.title("Proportional Function Distribution by Sample Group")
plt.ylabel("Sample Group", fontsize=12)
plt.xlabel("Function", fontsize=12)
plt.yticks(rotation=0, fontsize=10)  # Horizontal labels
plt.xticks(rotation=45, ha="right", fontsize=10)
plt.tight_layout()

# Save
plot_path = os.path.join(output_dir, "function_distribution_by_group.png")
plt.savefig(plot_path)
plt.show()
print(f"✅ Saved heatmap to: {plot_path}")
