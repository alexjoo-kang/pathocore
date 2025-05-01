import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

# === Paths ===
data_dir = "data"
output_dir = "output"
cluster_fasta_dir = os.path.join(output_dir, "clusters_fasta")
os.makedirs(cluster_fasta_dir, exist_ok=True)

# === Step 1: Load all amino acid sequences from the data/ directory ===
sequence_lookup = {}
for fasta_file in os.listdir(data_dir):
    if fasta_file.endswith(".fasta"):
        fasta_path = os.path.join(data_dir, fasta_file)
        for record in SeqIO.parse(fasta_path, "fasta"):
            sequence_lookup[record.id] = record

# === Step 2: Load all cluster assignments from output/ ===
cluster_to_records = {}
for result_file in os.listdir(output_dir):
    if result_file.startswith("results_") and result_file.endswith(".tsv"):
        df = pd.read_csv(os.path.join(output_dir, result_file), sep="\t")
        for _, row in df.iterrows():
            seq_id = row["sequence_id"]
            cluster_id = row["cluster_id"]
            if cluster_id not in cluster_to_records:
                cluster_to_records[cluster_id] = []
            if seq_id in sequence_lookup:
                cluster_to_records[cluster_id].append(sequence_lookup[seq_id])

# === Step 3: Write sequences into cluster-specific FASTA files ===
for cluster_id, records in cluster_to_records.items():
    out_path = os.path.join(cluster_fasta_dir, f"cluster_{cluster_id}.fasta")
    with open(out_path, "w") as handle:
        SeqIO.write(records, handle, "fasta")

print(f"✅ Cluster-separated FASTA files saved to '{cluster_fasta_dir}/'")
