import pandas as pd
import os

# Simulated function annotations for each known sequence ID
function_mapping = {
    "Insulin_pre": "Hormone",
    "MBP_tagged": "Transport",
    "14_3_3_protein": "Signal transduction",
    "Control_repeat": "Structural protein",
    "Calmodulin": "Calcium binding",

    "CHD8_variant": "Transcription regulation",
    "SHANK3_fragment": "Synaptic scaffolding",
    "SCN2A_short": "Ion channel",
    "FOXP1_fragment": "Transcription factor",
    "CNTNAP2_repeat": "Cell adhesion",

    "TP53_Nterm": "Tumor suppressor",
    "KRAS_tagged": "Signal transduction",
    "PIK3CA_fragment": "Kinase activity",
    "BRAF_fragment": "Kinase activity",
    "EGFR_like": "Receptor tyrosine kinase",

    "Keratin_repeat": "Cytoskeleton",
    "Control_fibrous": "Structural protein",
    "Histone_like": "Chromatin organization",
    "Hemoglobin_like": "Oxygen transport",
    "HSP90_motif": "Protein folding"
}

# Convert to DataFrame
df_function = pd.DataFrame(list(function_mapping.items()), columns=["sequence_id", "function"])

# Ensure the data directory exists
os.makedirs("data", exist_ok=True)

# Save as TSV
df_function.to_csv("data/function_labels.tsv", sep="\t", index=False)

print("✅ Mock functional annotation file saved to: data/function_labels.tsv")
