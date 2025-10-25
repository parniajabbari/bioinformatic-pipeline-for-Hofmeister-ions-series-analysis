import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse

# -----------------------------
# 1. Command-line arguments
# -----------------------------
parser = argparse.ArgumentParser(
    description="Generate pie chart of interaction types around a given ion."
)
parser.add_argument(
    "ion",
    type=str,
    help="Ion symbol (e.g., F, CL, SCN, SO4, etc.)"
)
parser.add_argument(
    "--output",
    type=str,
    default=None,
    help="Optional output directory (default: ~/Desktop/{ION}/)"
)
args = parser.parse_args()
ION_SYMBOL = args.ion.upper()

# -----------------------------
# 2. Define directories
# -----------------------------
BASE_DESKTOP = os.path.expanduser("~/Desktop")
BASE_DIR = args.output or os.path.join(BASE_DESKTOP, ION_SYMBOL)

input_file = os.path.join(BASE_DIR, "interaction_statistics_with_residue_counts.csv")
output_file = os.path.join(BASE_DIR, f"{ION_SYMBOL}_interaction_type_distribution.png")

# -----------------------------
# 3. Interaction color mapping
# -----------------------------
interaction_colors = {
    'H-bond with Main Chain': '#90EE90',  # Light Green
    'H-bond with Side Chain': '#87CEFA',  # Light Blue
    'Aliphatic': '#FFD700',               # Light Orange/Gold
    'Aromatic': '#FFFFE0',                # Light Yellow
    'Metal': '#FFA500',                    # Orange
    'Salt Bridge': '#D2B48C',             # Tan
    'Cofactor/Ligand': '#FFB6B9',         # Soft Pink
    'H-bond with Water': '#D8BFD8',       # Thistle (Light Purple)
    'Ion': '#B0E0E6'                       # Powder Blue
}

desired_order = [
    'Aromatic', 'Aliphatic', 'Metal', 'H-bond with Main Chain', 
    'H-bond with Side Chain', 'Cofactor/Ligand', 'H-bond with Water', 
    'Ion', 'Salt Bridge'
]

# -----------------------------
# 4. Main processing
# -----------------------------
if not os.path.exists(input_file):
    print(f"Error: Input file not found: {input_file}")
else:
    try:
        data = pd.read_csv(input_file)

        if 'Interaction Type' not in data.columns or 'Approx. Percentage in Shell' not in data.columns:
            raise KeyError("'Interaction Type' or 'Approx. Percentage in Shell' column missing.")

        data['Interaction Type'] = data['Interaction Type'].str.strip()

        interaction_percentage = data.groupby('Interaction Type')['Approx. Percentage in Shell'].sum()
        interaction_percentage = interaction_percentage[interaction_percentage > 0]

        interaction_percentage = interaction_percentage.reindex(
            [interaction for interaction in desired_order if interaction in interaction_percentage.index]
        )

        colors = [interaction_colors.get(interaction, '#808080') for interaction in interaction_percentage.index]

        plt.figure(figsize=(10, 10))
        plt.pie(
            interaction_percentage,
            labels=interaction_percentage.index,
            autopct=lambda p: f'{p:.2f}%' if p > 0 else '',
            startangle=140,
            colors=colors
        )
        plt.title(f"Distribution of Interaction Types around {ION_SYMBOL}", fontsize=16)
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.show()

        print(f"✅ Pie chart saved to: {output_file}")

    except Exception as e:
        print(f"An error occurred: {e}")

