import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from core.io_utils import load_fasta
from core.embedding import embed_proteins

# Set paths
data_dir = "data"
output_dir = "output"

# Initialize containers
all_embeddings = []
all_labels = []
all_groups = []

# Load and process each FASTA file and matching clustering results
for fasta_file in sorted(os.listdir(data_dir)):
    if not fasta_file.endswith(".fasta"):
        continue

    group_name = fasta_file.replace(".fasta", "")
    fasta_path = os.path.join(data_dir, fasta_file)
    tsv_path = os.path.join(output_dir, f"results_{group_name}.tsv")

    sequences = load_fasta(fasta_path)
    embeddings = embed_proteins(sequences)
    df_clusters = pd.read_csv(tsv_path, sep="\t")

    all_embeddings.append(embeddings)
    all_labels.extend(df_clusters["cluster_id"].tolist())
    all_groups.extend([group_name] * len(embeddings))

# Combine all data
X = np.vstack(all_embeddings)
y = np.array(all_labels)
groups = np.array(all_groups)

# Standardize
X_scaled = StandardScaler().fit_transform(X)

# t-SNE
tsne = TSNE(n_components=2, random_state=42, perplexity=5)
X_tsne = tsne.fit_transform(X_scaled)

# Plot
plt.figure(figsize=(10, 7))
scatter = plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=y, cmap="tab10", alpha=0.8)

# Optional: annotate by group
for i, group in enumerate(groups):
    plt.annotate(group, (X_tsne[i, 0] + 0.5, X_tsne[i, 1]), fontsize=7, alpha=0.6)

plt.title("t-SNE Visualization of AAC Embeddings by Cluster")
plt.xlabel("t-SNE 1")
plt.ylabel("t-SNE 2")
plt.colorbar(scatter, label="Cluster ID")
plt.tight_layout()

# Save plot
plot_path = os.path.join(output_dir, "tsne_cluster_visualization.png")
plt.savefig(plot_path)
plt.show()

print(f"✅ t-SNE plot saved to: {plot_path}")
