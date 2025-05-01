import os
import pandas as pd
from collections import defaultdict
from Bio import SeqIO

# Paths
data_dir = "data"
output_dir = "output"
cluster_dir = os.path.join(output_dir, "clusters_fasta")
function_file = os.path.join(data_dir, "function_labels.tsv")
out_summary_path = os.path.join(output_dir, "functional_summary.tsv")

# Load function annotation table
df_func = pd.read_csv(function_file, sep="\t")
func_map = dict(zip(df_func["sequence_id"], df_func["function"]))

# Analyze function distribution per cluster
function_summary = []

for fname in os.listdir(cluster_dir):
    if fname.endswith(".fasta"):
        cluster_id = fname.replace("cluster_", "").replace(".fasta", "")
        records = list(SeqIO.parse(os.path.join(cluster_dir, fname), "fasta"))
        func_count = defaultdict(int)

        for record in records:
            func = func_map.get(record.id, "Unknown")
            func_count[func] += 1

        for func, count in func_count.items():
            function_summary.append({
                "cluster_id": int(cluster_id),
                "function": func,
                "count": count
            })

# Save summary
df_summary = pd.DataFrame(function_summary)
df_summary.sort_values(by=["cluster_id", "count"], ascending=[True, False], inplace=True)
df_summary.to_csv(out_summary_path, sep="\t", index=False)

print(f"\u2705 Functional summary saved to: {out_summary_path}")
