import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse

# ------------------------
# CLI Arguments
# ------------------------
parser = argparse.ArgumentParser(description="Plot secondary structure distribution per ion")
parser.add_argument("ion", help="Ion symbol (e.g., F, SCN, AL, LCP)")
parser.add_argument("--output", help="Output directory (default: ~/Desktop/<ION>)")
args = parser.parse_args()

ION_SYMBOL = args.ion.upper()
BASE_DIR = args.output or os.path.join(os.path.expanduser("~"), "Desktop", ION_SYMBOL)
os.makedirs(BASE_DIR, exist_ok=True)

# ------------------------
# Input and output files
# ------------------------
INPUT_FILE = os.path.join(BASE_DIR, "result.csv")
OUTPUT_FILE = os.path.join(BASE_DIR, "secondary_structure_distribution.png")  # fixed name

if not os.path.exists(INPUT_FILE):
    raise FileNotFoundError(f"Input CSV not found: {INPUT_FILE}")

# ------------------------
# Load Data
# ------------------------
data = pd.read_csv(INPUT_FILE)

# Filter by distance <= 5 Å and drop NaNs
data_filtered = data[data['Distance (Å)'] <= 5].dropna(subset=['Secondary Structure'])

# Normalize case for matching keys
data_filtered['Secondary Structure'] = data_filtered['Secondary Structure'].str.lower()

# ------------------------
# Structure Order and Colors
# ------------------------
structure_order = ['coil', 'beta bridge', '3-10 helix', 'turn', 'pi helix',
                   'bend', 'extended strand', 'alpha helix', 'polyproline ii helix']

structure_colors = {
    'coil': '#D3D3D3',  # lighter grey
    'beta bridge': '#ADD8E6',  # lighter light blue
    '3-10 helix': '#98FB98',  # lighter green
    'turn': '#FFCC99',  # lighter orange
    'pi helix': '#D8BFD8',  # lighter thistle
    'bend': '#FFB6C1',  # lighter pink
    'extended strand': '#E0FFFF',  # lighter cyan
    'alpha helix': '#FF7F50',  # lighter coral
    'polyproline ii helix': '#FFFACD'  # lighter lemon chiffon
}

structure_labels = {
    'coil': 'coil',
    'beta bridge': 'β-bridge',
    '3-10 helix': '3-10 helix',
    'turn': 'turn',
    'pi helix': 'π-helix',
    'bend': 'bend',
    'extended strand': 'extended strand',
    'alpha helix': 'α-helix',
    'polyproline ii helix': 'PPII-helix'
}

# ------------------------
# Filter and reorder
# ------------------------
available_structures = [s for s in structure_order if s in data_filtered['Secondary Structure'].unique()]
structure_counts = data_filtered['Secondary Structure'].value_counts(normalize=True) * 100
structure_counts = structure_counts[available_structures]
colors = [structure_colors[s] for s in available_structures]
labels = [structure_labels[s] for s in available_structures]

# ------------------------
# Plot Pie Chart
# ------------------------
plt.figure(figsize=(12, 12))
wedges, texts, autotexts = plt.pie(
    structure_counts,
    labels=labels,
    autopct=lambda p: f'{p:.1f}%' if p > 0 else '',
    startangle=90,
    colors=colors,
    textprops={'fontsize': 20, 'fontweight': 'bold'}
)

for autotext in autotexts:
    autotext.set_fontsize(20)
    autotext.set_fontweight('bold')

plt.axis('equal')  # Ensure circular pie

# Save to fixed output file
plt.savefig(OUTPUT_FILE, dpi=600, bbox_inches='tight')
plt.show()

print(f"✅ Secondary structure pie chart saved to: {OUTPUT_FILE}")
