import pandas as pd
import matplotlib.pyplot as plt
import os

# Define the path to the directory and the file name
directory = "/Users/respina/desktop/GAI"
file_name = "interaction_statistics_with_residue_counts.csv"
file_path = os.path.join(directory, file_name)

# Define a lighter color mapping for interaction types
interaction_colors = {
    'H-bond with Main Chain': '#90EE90',  # Light Green
    'H-bond with Side Chain': '#87CEFA',  # Light Blue
    'Aliphatic': '#FFD700',  # Light Orange/Gold
    'Aromatic': '#FFFFE0',  # Light Yellow
    'Metal': '#FFA500',     # Orange
    'Salt Bridge': '#D2B48C',  # Tan
    'Cofactor/Ligand': '#FFB6B9',  # Soft Pink
    'H-bond with Water': '#D8BFD8',  # Thistle (Light Purple)
    'Ion': '#B0E0E6'  # Powder Blue
}

# Desired order of interactions
desired_order = [
    'Aromatic', 'Aliphatic', 'Metal', 'H-bond with Main Chain', 
    'H-bond with Side Chain', 'Cofactor/Ligand', 'H-bond with Water', 
    'Ion', 'Salt Bridge'
]

# Check if the file exists
if not os.path.exists(file_path):
    print(f"Error: The file '{file_name}' was not found in the directory '{directory}'.")
else:
    try:
        # Read the CSV file into a DataFrame
        data = pd.read_csv(file_path)

        # Check if required columns are present
        if 'Interaction Type' in data.columns and 'Approx. Percentage in Shell' in data.columns:

            # Remove any whitespace or unexpected characters in column values
            data['Interaction Type'] = data['Interaction Type'].str.strip()

            # Group by 'Interaction Type' and sum the percentages for each type
            interaction_percentage = data.groupby('Interaction Type')['Approx. Percentage in Shell'].sum()

            # Filter out interaction types with 0% to avoid showing them on the chart
            interaction_percentage = interaction_percentage[interaction_percentage > 0]

            # Ensure that the remaining interaction types appear in the desired order
            interaction_percentage = interaction_percentage.reindex(
                [interaction for interaction in desired_order if interaction in interaction_percentage.index]
            )

            # Ensure that no interaction type is missing from the color mapping
            colors = [
                interaction_colors.get(interaction, '#808080')  # Default: Gray
                for interaction in interaction_percentage.index
            ]

            # Plot the pie chart with more precise percentages
            plt.figure(figsize=(10, 10))
            plt.pie(
                interaction_percentage,
                labels=interaction_percentage.index,
                autopct=lambda p: f'{p:.2f}%' if p > 0 else '',
                startangle=140,
                colors=colors
            )

            # Add a title to the chart
            plt.title("Distribution of Interaction Types around ION", fontsize=16)

            # Save the figure to the same directory
            output_path = os.path.join(directory, "interaction_type_distribution.png")
            plt.savefig(output_path, dpi=300, bbox_inches='tight')

            # Display the plot
            plt.show()

        else:
            print("Error: 'Interaction Type' or 'Approx. Percentage in Shell' column not found in the file.")

    except Exception as e:
        print(f"An error occurred: {e}")
