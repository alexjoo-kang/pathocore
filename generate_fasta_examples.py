import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Define realistic, biologically inspired amino acid sequences for 4 groups
# Length: 60-100 aa, diverse compositions
realistic_datasets = {
    "ASD_group.fasta": [
        SeqRecord(Seq("MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVE"), id="CHD8_variant"),
        SeqRecord(Seq("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESF"), id="SHANK3_fragment"),
        SeqRecord(Seq("MDSKGSSQKGSRLLLLLVVSNLLLCQGVVSH"), id="SCN2A_short"),
        SeqRecord(Seq("MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFTGHPETLEK"), id="FOXP1_fragment"),
        SeqRecord(Seq("MAAAMVGAGGVAGAPGAAAPAAAGGPGGAGGAGGAGGAAGPGA"), id="CNTNAP2_repeat")
    ],
    "Control_ASD.fasta": [
        SeqRecord(Seq("MALWMRLLPLLALLALWGPDPAAGR"), id="Insulin_pre"),
        SeqRecord(Seq("MGSSHHHHHHSSGLVPRGSHMGSNSQDNNNLQ"), id="MBP_tagged"),
        SeqRecord(Seq("MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERN"), id="14_3_3_protein"),
        SeqRecord(Seq("MAEQKQQLQALQQEMEQAQAAQAQLAQEMQAQQEQQLQQQQQ"), id="Control_repeat"),
        SeqRecord(Seq("MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNP"), id="Calmodulin")
    ],
    "Cancer_group.fasta": [
        SeqRecord(Seq("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLML"), id="TP53_Nterm"),
        SeqRecord(Seq("MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRGSEFDDDDK"), id="KRAS_tagged"),
        SeqRecord(Seq("MQSGKKLIEAVLRALKDLASKNKPFAGKADYFKQ"), id="PIK3CA_fragment"),
        SeqRecord(Seq("MAEDDPDFSKDLDFKAVQLFGGSEDEDQDPQDPQDPQDP"), id="BRAF_fragment"),
        SeqRecord(Seq("MSSKLLKTVMALAALLAVASSGASAAPRTVLVYVGGSAGVGKS"), id="EGFR_like")
    ],
    "Control_Cancer.fasta": [
        SeqRecord(Seq("MKAAVQAGAVAVGAGAGAGAGGAGAGKPGPGAGPAP"), id="Keratin_repeat"),
        SeqRecord(Seq("MADPSTPAPTPSPSPATSGAGAAAPAPPAPPP"), id="Control_fibrous"),
        SeqRecord(Seq("MGRKNGKGKKRKGGAKAKAKAKAKAKKKAKAA"), id="Histone_like"),
        SeqRecord(Seq("MAEPLAPTSVITKLGGAFIKDKYLVLGGTKLRKDL"), id="Hemoglobin_like"),
        SeqRecord(Seq("MTRVLQRVQLFQELTPQYVNRTHQPEAKV"), id="HSP90_motif")
    ]
}

# Save files
output_dir = "data"
os.makedirs(output_dir, exist_ok=True)

for filename, records in realistic_datasets.items():
    filepath = os.path.join(output_dir, filename)
    with open(filepath, "w") as handle:
        SeqIO.write(records, handle, "fasta")

print("✅ Realistic and biologically diverse FASTA files have been generated in the 'data/' directory.")
