import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Generate protein structure coordinates
protein_coordinates = np.array([(0, 0, 0), (-1, 0, 1), (-1,0,0), (-2,1,0)])

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
