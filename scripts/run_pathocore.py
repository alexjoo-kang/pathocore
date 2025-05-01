import os
from core.io_utils import load_fasta
from core.embedding import embed_proteins
from core.clustering import cluster_embeddings
import numpy as np
import pandas as pd

# Constants
data_dir = "data"
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

def process_fasta_file(filepath, n_clusters=3):
    filename = os.path.basename(filepath)
    group_name = os.path.splitext(filename)[0]

    print(f"\n📂 Processing {group_name}...")
    
    # Step 1: Load sequences
    sequences = load_fasta(filepath)
    sequence_ids = list(sequences.keys())
    print(f"✅ Loaded {len(sequences)} sequences.")

    # Step 2: Embed sequences
    embeddings = embed_proteins(sequences)
    print(f"✅ Embedding shape: {embeddings.shape}")

    # Step 3: Cluster embeddings
    labels = cluster_embeddings(embeddings, n_clusters=n_clusters)
    print(f"✅ Assigned {len(set(labels))} clusters.")

    # Step 4: Save results
    results_df = pd.DataFrame({
        "sequence_id": sequence_ids,
        "cluster_id": labels
    })
    out_path = os.path.join(output_dir, f"results_{group_name}.tsv")
    results_df.to_csv(out_path, sep="\t", index=False)
    print(f"📄 Saved results to: {out_path}")

def main():
    fasta_files = [f for f in os.listdir(data_dir) if f.endswith(".fasta")]
    if not fasta_files:
        print("❌ No FASTA files found in the 'data/' directory.")
        return

    for fasta in fasta_files:
        process_fasta_file(os.path.join(data_dir, fasta), n_clusters=4)

if __name__ == "__main__":
    main()
