import os
import pandas as pd
from collections import defaultdict

def load_taxonomy_mapping(taxonomy_file="data/kraken_labels.tsv"):
    """
    Load taxonomy mapping from TSV file. Assumes columns: sequence_id, taxonomic_label
    """
    if not os.path.exists(taxonomy_file):
        raise FileNotFoundError(f"❌ Taxonomy file not found: {taxonomy_file}")
    
    df = pd.read_csv(taxonomy_file, sep="\t")
    return dict(zip(df["sequence_id"], df["taxonomic_label"]))

def summarize_taxonomy_by_cluster(cluster_fasta_dir="output/clusters_fasta", taxonomy_file="data/kraken_labels.tsv"):
    """
    Summarize the distribution of taxonomic labels across clusters.
    """
    from Bio import SeqIO

    taxonomy_map = load_taxonomy_mapping(taxonomy_file)
    cluster_taxonomy_summary = defaultdict(list)

    for fasta_file in os.listdir(cluster_fasta_dir):
        if fasta_file.endswith(".fasta"):
            cluster_id = fasta_file.replace("cluster_", "").replace(".fasta", "")
            fasta_path = os.path.join(cluster_fasta_dir, fasta_file)
            for record in SeqIO.parse(fasta_path, "fasta"):
                tax_label = taxonomy_map.get(record.id, "Unknown")
                cluster_taxonomy_summary[cluster_id].append(tax_label)

    # Summarize into DataFrame
    summary_rows = []
    for cluster_id, labels in cluster_taxonomy_summary.items():
        label_counts = pd.Series(labels).value_counts().to_dict()
        for label, count in label_counts.items():
            summary_rows.append({
                "cluster_id": cluster_id,
                "taxonomic_label": label,
                "count": count
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_path = "output/taxonomy_by_cluster.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    print("✅ Taxonomy summary saved to:", summary_path)
    return summary_df
