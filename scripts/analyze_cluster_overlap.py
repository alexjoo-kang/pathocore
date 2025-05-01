import os
import pandas as pd
from collections import defaultdict

output_dir = "output"
tsv_files = [f for f in os.listdir(output_dir) if f.endswith(".tsv") and not f.startswith("cluster_overlap_summary")]

cluster_map = defaultdict(set)

for file in tsv_files:
    group = file.replace("results_", "").replace(".tsv", "")
    df = pd.read_csv(os.path.join(output_dir, file), sep="\t")
    for cid in df["cluster_id"].unique():
        cluster_map[int(cid)].add(group)  

summary_data = [
    {"cluster_id": cid, "num_groups": len(groups), "groups": ", ".join(sorted(groups))}
    for cid, groups in sorted(cluster_map.items())
]

summary_df = pd.DataFrame(summary_data).sort_values(by="num_groups", ascending=False)
summary_df.to_csv("output/cluster_overlap_summary.tsv", sep="\t", index=False)

print("✅ Cluster overlap analysis complete!")
print("📄 Saved to: output/cluster_overlap_summary.tsv")
print(summary_df)
