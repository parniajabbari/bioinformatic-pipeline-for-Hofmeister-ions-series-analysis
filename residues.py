import pandas as pd
import matplotlib.pyplot as plt
import os

# ======================================================
# USER SETTINGS
# ======================================================
ion_symbol = "GAI"  # ← change this only
base_path = "/Users/respina/desktop"
# ======================================================

# Construct CSV path
csv_path = os.path.join(base_path, ion_symbol, "interaction_statistics_with_residue_counts.csv")

if not os.path.exists(csv_path):
    raise FileNotFoundError(f"❌ File not found: {csv_path}")

# Read CSV
df = pd.read_csv(csv_path)

print(f"✅ Loaded data for ion '{ion_symbol}' from:\n{csv_path}\n")
print(df.head())

# ======================================================
# PROCESSING
# ======================================================
# Ensure required columns exist
if "Residue Name" not in df.columns or "Unique Residues in Shell" not in df.columns:
    raise KeyError("❌ Expected columns not found! Ensure CSV has 'Residue Name' and 'Unique Residues in Shell'.")

# Group by residue name and sum counts
df_grouped = (
    df.groupby("Residue Name", as_index=False)["Unique Residues in Shell"]
      .sum()
      .rename(columns={"Unique Residues in Shell": "Residue Count"})
)

# Filter for standard amino acids (optional, remove this block if you want all)
standard_aas_order = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR",
    "TRP", "TYR", "VAL"
]
df_filtered = df_grouped[df_grouped["Residue Name"].isin(standard_aas_order)]

# Sort residues by predefined amino acid order
df_filtered["Residue Name"] = pd.Categorical(df_filtered["Residue Name"], categories=standard_aas_order, ordered=True)
df_filtered = df_filtered.sort_values("Residue Name")

# ======================================================
# PLOT
# ======================================================
plt.figure(figsize=(14, 7))
plt.plot(df_filtered["Residue Name"], df_filtered["Residue Count"], marker='o', label='Residue Count')
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.xlabel("Residues", fontsize=12, fontweight='bold')
plt.ylabel("Count unique residues in the shell", fontsize=12, fontweight='bold')
plt.title(f"Residue Frequency around {ion_symbol} Ion", fontsize=14, fontweight='bold')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()

# ======================================================
# SAVE & SHOW
# ======================================================
output_path = os.path.join(base_path, ion_symbol, f"{ion_symbol}_residue_count_plot_summed.png")
plt.savefig(output_path, dpi=300)
plt.show()

print(f"\n📊 Summed and ordered residue plot saved successfully at:\n{output_path}")

