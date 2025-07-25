import pandas as pd
import matplotlib.pyplot as plt

# File path
file_path = '/home/parnia/parnia/F/result.csv'

# Load the data
data = pd.read_csv(file_path)

# Filter by distance <= 5 Å and drop NaNs in secondary structure
data_filtered = data[data['Distance (Å)'] <= 5].dropna(subset=['Secondary Structure'])

# Normalize case for matching keys
data_filtered['Secondary Structure'] = data_filtered['Secondary Structure'].str.lower()

# Define the order and lighter colors (one degree lighter)
structure_order = ['coil', 'beta bridge', '3-10 helix', 'turn', 'pi helix', 'bend', 'extended strand', 'alpha helix', 'polyproline ii helix']
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

# Label mapping for plot
structure_labels = {
    'coil': 'coil',
    'beta bridge': 'β-bridge',
    '3-10 helix': '3-10 helix',
    'turn': 'turn',
    'pi helix': 'π-helix ',
    'bend': 'bend',
    'extended strand': 'extended strand',
    'alpha helix': 'α-helix',
    'polyproline ii helix': 'PPII-helix'
}

# Filter and reorder structure counts
available_structures = [s for s in structure_order if s in data_filtered['Secondary Structure'].unique()]
structure_counts = data_filtered['Secondary Structure'].value_counts(normalize=True) * 100
structure_counts = structure_counts[available_structures]
colors = [structure_colors[s] for s in available_structures]
labels = [structure_labels[s] for s in available_structures]

# Plot the pie chart and capture label and autopct objects
plt.figure(figsize=(12, 12))
wedges, texts, autotexts = plt.pie(
    structure_counts,
    labels=labels,
    autopct=lambda p: f'{p:.1f}%' if p > 0 else '',
    startangle=90,
    colors=colors,
    textprops={'fontsize': 20, 'fontweight': 'bold'}
)

# Set font size and weight for percentage text explicitly
for autotext in autotexts:
    autotext.set_fontsize(20)
    autotext.set_fontweight('bold')


plt.axis('equal')  # Equal aspect ratio ensures the pie chart is circular

# Save the plot
output_path = '/home/parnia/parnia/F/secondary_structure_distribution.png'
plt.savefig(output_path, dpi=600, bbox_inches='tight')

# Show the plot
plt.show()


