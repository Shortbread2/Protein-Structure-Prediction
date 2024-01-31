import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

HP_AMINO_ARRAY = "hpph"

def generate_protein_structure(coordinates):
    global HP_AMINO_ARRAY
    # Define amino acid colors and markers based on the HP amino acids
    aa_settings = {'h': {'color': 'red', 'marker': '^', 'label': 'Hydrophobic'},
                   'p': {'color': 'blue', 'marker': 'o', 'label': 'Polar'}}

    # Plot the 3D structure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i, aa_type in enumerate(HP_AMINO_ARRAY):
        x, y, z = coordinates[i]
        settings = aa_settings[aa_type]
        ax.scatter(x, y, z, s=100, c=settings['color'], marker=settings['marker'], label=settings['label'])

    # Draw lines between consecutive protein coordinates
    for i in range(len(coordinates) - 1):
        ax.plot([coordinates[i][0], coordinates[i + 1][0]],
                [coordinates[i][1], coordinates[i + 1][1]],
                [coordinates[i][2], coordinates[i + 1][2]], c='blue', linestyle='dashed')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Protein Structure')
    ax.legend(['Polar', 'Hydrophobic'])
    plt.show()

protein_coordinates = np.array([(0, 0, 0), (-1, 0, 1), (-1,0,0), (-2,1,0)])

generate_protein_structure(protein_coordinates)
