import numpy as np
from tqdm import tqdm

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

def aac_vector(sequence):
    """
    Converts a protein sequence to an amino acid composition vector (20-dim).
    Each element is the frequency of an amino acid in the sequence.
    """
    sequence = sequence.upper()
    length = len(sequence)
    counts = [sequence.count(aa) / length for aa in AMINO_ACIDS]
    return np.array(counts)

def embed_proteins(sequences, batch_size=8):
    """
    Embed proteins using amino acid composition (AAC) vectors.

    Args:
        sequences (List[str]): Protein sequences
        batch_size (int): Ignored, kept for API consistency

    Returns:
        np.ndarray: shape (n_sequences, 20)
    """
    embeddings = []
    for seq in tqdm(sequences, desc="Embedding with AAC"):
        embeddings.append(aac_vector(seq))
    return np.array(embeddings)
