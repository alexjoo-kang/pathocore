import subprocess

steps = [
    "scripts.export_cluster_fastas",
    "scripts.analyze_cluster_properties",
    "scripts.analyze_cluster_aa_properties",
    "scripts.analyze_functional_summary",
    "scripts.analyze_taxonomy",
    "scripts.plot_tsne",
    "scripts.plot_cluster_distribution",
    "scripts.plot_cluster_heatmap",
    "scripts.plot_function_heatmap",
    "scripts.plot_function_by_group"
]

for step in steps:
    print(f"\n🔁 Running {step} ...")
    subprocess.run(["python3", "-m", step], check=True)

print("\n✅ All analysis and plots completed!")
