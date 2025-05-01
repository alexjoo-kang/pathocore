from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import pandas as pd

# Paths
data_dir = "data"
output_dir = "output"
cluster_fasta_dir = os.path.join(output_dir, "clusters_fasta")
os.makedirs(cluster_fasta_dir, exist_ok=True)

# Load all sequence data and build a lookup
sequence_lookup = {}
for fasta_file in os.listdir(data_dir):
    if fasta_file.endswith(".fasta"):
        records = list(SeqIO.parse(os.path.join(data_dir, fasta_file), "fasta"))
        for record in records:
            sequence_lookup[record.id] = record

# Load all cluster assignments and write grouped FASTA files
cluster_to_records = {}

for tsv_file in os.listdir(output_dir):
    if tsv_file.startswith("results_") and tsv_file.endswith(".tsv"):
        df = pd.read_csv(os.path.join(output_dir, tsv_file), sep="\t")
        for _, row in df.iterrows():
            seq_id = row["sequence_id"]
            cluster_id = row["cluster_id"]
            if cluster_id not in cluster_to_records:
                cluster_to_records[cluster_id] = []
            if seq_id in sequence_lookup:
                cluster_to_records[cluster_id].append(sequence_lookup[seq_id])

# Write one FASTA file per cluster
for cluster_id, records in cluster_to_records.items():
    out_path = os.path.join(cluster_fasta_dir, f"cluster_{cluster_id}.fasta")
    with open(out_path, "w") as handle:
        SeqIO.write(records, handle, "fasta")

import os
cluster_files = os.listdir(cluster_fasta_dir)
cluster_files
