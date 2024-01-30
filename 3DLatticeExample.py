import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def generate_fcc_lattice(a, n):
    coordinates = []

    for i in range(n):
        for j in range(n):
            for k in range(n):
                coordinates.append([i * a, j * a, k * a])
                coordinates.append([(i + 0.5) * a, j * a, (k + 0.5) * a])
                coordinates.append([i * a, (j + 0.5) * a, (k + 0.5) * a])
                coordinates.append([(i + 0.5) * a, (j + 0.5) * a, k * a])

    return np.array(coordinates)

def generate_protein_structure():
    protein_coordinates = np.array([(0, 0, 0), (-1, 0, 1), (-1,0,0), (-2,1,0)])
    return protein_coordinates

a = 1.0  # lattice constant
n = 5    # number of unit cells along each axis

# Generate FCC lattice coordinates
fcc_coordinates = generate_fcc_lattice(a, n)

# Generate protein structure coordinates
protein_coordinates = generate_protein_structure()

# Plot the 3D FCC lattice with the protein structure and lines
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(protein_coordinates[:, 0], protein_coordinates[:, 1], protein_coordinates[:, 2], s=100, c='blue', marker='^', label='Protein')

# Draw lines between consecutive protein coordinates
for i in range(len(protein_coordinates) - 1):
    ax.plot([protein_coordinates[i, 0], protein_coordinates[i + 1, 0]],
            [protein_coordinates[i, 1], protein_coordinates[i + 1, 1]],
            [protein_coordinates[i, 2], protein_coordinates[i + 1, 2]], c='blue')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('3D FCC Lattice with Protein and Lines')
ax.legend()
plt.show()
