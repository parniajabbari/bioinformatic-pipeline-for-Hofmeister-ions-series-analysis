import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

# ---------------------------
# --- ARGUMENT PARSER ------
# ---------------------------
parser = argparse.ArgumentParser(description="Plot residue counts for a specific ion")
parser.add_argument("ion", help="Ion symbol, e.g. LCP, GAI")
parser.add_argument("--output", default=os.getcwd(), help="Output directory containing CSV and to save plots")
args = parser.parse_args()

ION_SYMBOL = args.ion
OUTPUT_DIR = args.output

# ---------------------------
# --- CSV PATH -------------
# ---------------------------
csv_path = os.path.join(OUTPUT_DIR, "interaction_statistics_with_residue_counts.csv")

if not os.path.exists(csv_path):
    raise FileNotFoundError(f"❌ File not found: {csv_path}")

# ---------------------------
# --- LOAD CSV -------------
# ---------------------------
df = pd.read_csv(csv_path)
print(f"✅ Loaded data for ion '{ION_SYMBOL}' from:\n{csv_path}\n")
print(df.head())

# ---------------------------
# --- PROCESSING -----------
# ---------------------------
if "Residue Name" not in df.columns or "Unique Residues in Shell" not in df.columns:
    raise KeyError("❌ Expected columns not found! Ensure CSV has 'Residue Name' and 'Unique Residues in Shell'.")

df_grouped = (
    df.groupby("Residue Name", as_index=False)["Unique Residues in Shell"]
      .sum()
      .rename(columns={"Unique Residues in Shell": "Residue Count"})
)

standard_aas_order = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR",
    "TRP", "TYR", "VAL"
]
df_filtered = df_grouped[df_grouped["Residue Name"].isin(standard_aas_order)]
df_filtered["Residue Name"] = pd.Categorical(df_filtered["Residue Name"], categories=standard_aas_order, ordered=True)
df_filtered = df_filtered.sort_values("Residue Name")

# ---------------------------
# --- PLOT -----------------
# ---------------------------
plt.figure(figsize=(14, 7))
plt.plot(df_filtered["Residue Name"], df_filtered["Residue Count"], marker='o', label='Residue Count')
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.xlabel("Residues", fontsize=12, fontweight='bold')
plt.ylabel("Count unique residues in the shell", fontsize=12, fontweight='bold')
plt.title(f"Residue Frequency around {ION_SYMBOL} Ion", fontsize=14, fontweight='bold')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()

# ---------------------------
# --- SAVE & SHOW ----------
# ---------------------------
output_path = os.path.join(OUTPUT_DIR, f"{ION_SYMBOL}_residue_count_plot_summed.png")
plt.savefig(output_path, dpi=300)
plt.show()

print(f"\n📊 Summed and ordered residue plot saved successfully at:\n{output_path}")

