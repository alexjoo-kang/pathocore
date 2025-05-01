from Bio import SeqIO

def load_fasta(filepath):
    """
    Load sequences from a FASTA file and return them as a dictionary.

    Args:
        filepath (str): Path to the FASTA file.

    Returns:
        dict: A dictionary with sequence IDs as keys and sequences (str) as values.
    """
    sequences = {}
    for record in SeqIO.parse(filepath, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences
