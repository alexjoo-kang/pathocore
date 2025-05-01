from sklearn.cluster import KMeans
import numpy as np

def cluster_embeddings(embeddings, n_clusters=3, random_state=42):
    """
    Cluster protein embeddings using KMeans.

    Args:
        embeddings (np.ndarray): Array of shape (n_samples, n_features)
        n_clusters (int): Number of clusters to create
        random_state (int): Seed for reproducibility

    Returns:
        List[int]: Cluster labels for each embedding
    """
    kmeans = KMeans(n_clusters=n_clusters, random_state=random_state)
    labels = kmeans.fit_predict(embeddings)
    return labels
