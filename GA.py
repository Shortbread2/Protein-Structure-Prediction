# initial format is heavily inspired from geaks for geaks - https://www.geeksforgeeks.org/genetic-algorithms/
# based off of the paper 'An Enhanced Genetic Algorithm for Ab Initio Protein Structure Prediction' - http://www.cis.umassd.edu/~fkhatib/Papers/IEEE_TEC.pdf

"""
notes:
abcdefghijkl - each letter will represent a direction on a 3d plane eg (-1,0,1)
the HP_AMINO_ACID will be the actual sequance of amino acids e.g (hpphhphpph)
TODO:
need to edit the create_a_gene (in other words the chromosome creation) to have "Self-avoiding walk"
need to change fitness, fitness function, mutation and crossover
"""

import sys
import random
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import distance
  
# Number of individuals in each generation 
POPULATION_SIZE = 1
  
# the valid letters that makes up a gene (when a gene is created it will check this varaible for what letters its allowed to have) 
GENES = 'abcdefghijkl'
  
# the array of amino acids, only in Hydrophobic or polar
HP_AMINO_ACID = "hpphphhphhp"

# for debuging
DEBUG_MODE = False
  
class Individual(object): 
    ''' 
    Class representing individual in population 
    '''
    def __init__(self, chromosome,coordinates): 
        self.chromosome = chromosome
        self.coordinates = coordinates 
        self.fitness = self.cal_fitness() 
  
    @classmethod
    def create_a_gene(self): 
        ''' 
        this actually builds the gene/ individual 
        '''
        global GENES 
        gene = random.choice(GENES) 
        return gene 
  
    @classmethod
    def create_gnome(self): 
        ''' 
        create chromosome or string of genes 
        '''
        global HP_AMINO_ACID 
        gnome_len = len(HP_AMINO_ACID) 
        return [self.create_a_gene() for _ in range(gnome_len)] 
  
    def mate(self, par2): 
        print("Mate chromosome")
        ''' 
        Perform mating and produce new offspring 
        '''
  
        # chromosome for offspring 
        child_chromosome = [] 
        for gp1, gp2 in zip(self.chromosome, par2.chromosome):     
  
            # random probability   
            prob = random.random() 
  
            # if prob is less than 0.45, insert gene 
            # from parent 1  
            if prob < 0.45: 
                child_chromosome.append(gp1) 
  
            # if prob is between 0.45 and 0.90, insert 
            # gene from parent 2 
            elif prob < 0.90: 
                child_chromosome.append(gp2) 
  
            # otherwise insert random gene(mutate),  
            # for maintaining diversity 
            else: 
                child_chromosome.append(self.create_a_gene()) 
  
        # create new Individual(offspring) using  
        # generated chromosome for offspring 
        return Individual(child_chromosome) 
  
    def cal_fitness(self):
        '''
        psudo code from the paper 'An Enhanced Genetic Algorithm for Ab Initio Protein Structure Prediction'
        
        
        AA: Amino acid array
        N: No. of amino acids in the sequence
        
        AA ←getAminoAcid(X1);
        for (i ← 0 to N - 1) do
            for (j ← i + 2 to N - 1) do
                if (AcidT ype[i] = AcidT ype[j] = 'H') then
                    pointI ← AA[i];
                    pointJ ← AA[j];
                    sqrD ← getSqrDist(pointI, pointJ);
                    if (sqrD = 2) then
                        fitness ← fitness - 1;
        update(X1, fitness);
        return X1;
        '''
        global HP_AMINO_ACID
        gnome_len = len(HP_AMINO_ACID)
        fitness = 0
        
        
        #code to offset lack of self avoiding walk for now, gives a penalty to coordinates that collide (code is temporary)
        for i in range(len(self.coordinates)):
            for j in range(i + 1, len(self.coordinates)):
                if  self.coordinates[i] ==  self.coordinates[j]:
                    if DEBUG_MODE: print("overlap with :", i ," and ", j)
                    fitness -= 5
        
        # checks if any h amino acids are too close not including the chemically bonded h's
        for i in range(gnome_len):
            for j in range(i + 2, gnome_len):
                if (HP_AMINO_ACID[i] == 'h') & (HP_AMINO_ACID[i] == HP_AMINO_ACID[j]):
                    sqrD = distance.euclidean(self.coordinates[i], self.coordinates[j])
                    if sqrD <= 2:
                        if DEBUG_MODE: print("index1: ", i, "index2: ", j)
                        fitness -= 1
        
        return fitness
    
def main():
    
    if "--help" in sys.argv[1:]:
        print("-----------------")
        print("Commands: --help, input:hpphphphphp (for example), input:random (this will prompt a user input for a number which can only be higher than 2)")
        print("population: (after this command enter in a number > 0 without spaces e.g. population:20), -debug (enables debug messages and stuff)")
        print("\nwill add feature to add option to use amino acid names for input instead of just having the option to enter in the denotions h or p")
        print("-----------------")
        sys.exit(0)
        
    global HP_AMINO_ACID
    global POPULATION_SIZE
    global DEBUG_MODE
    
    for command in sys.argv[1:]:
        if DEBUG_MODE: print(command)
        validCommand = False
        if "input:" in command:
            validCommand = True
            value = command.split(":")[1].lower()
            if value == "random":
                while True:
                    userNumInput = input("\nplease enter in how many amino acids you want to generate\n")
                    if userNumInput.isdigit():
                        userNumInput = int(userNumInput)
                        if userNumInput <= 2:
                            print("invalid user input, please enter a number greater than 2")
                        else:
                            HP_AMINO_ACID = ''.join(random.choice(['p', 'h']) for _ in range(userNumInput))
                            break
                    else:
                        print("invalid user input, only numbers are allowed")
            else:
                print(value)
                if all(char in {'p', 'h'} for char in value):
                    HP_AMINO_ACID = value
                else:
                    print("Invalid input. Please enter 'p' or 'h'.")
                    sys.exit(0)
        if "population:" in command:
            validCommand = True
            value = command.split(":")[1].lower()
            if value.isdigit():
                value = int(value)
                if value < 1:
                    print("invalid population input, please enter a number greater than 0")
                    sys.exit(0)
                else:
                    POPULATION_SIZE = value
                    continue
            else:
                print("invalid population input, only numbers are allowed")
                sys.exit(0)
        if "-debug" in command:
            validCommand = True
            DEBUG_MODE = True
            print(DEBUG_MODE)
        if not validCommand:
            print("invalid command:", command)
            sys.exit(0)

    coordinates = {
        'A': (1, 1, 0),
        'B': (-1, -1, 0),
        'C': (-1, 1, 0),
        'D': (1, -1, 0),
        'E': (0, 1, 1),
        'F': (0, -1, -1),
        'G': (0, 1, -1),
        'H': (0, -1, 1),
        'I': (-1, 0, -1),
        'J': (1, 0, 1),
        'K': (-1, 0, 1),
        'L': (1, 0, -1),
    }
  
    generation = 1
  
    found = False
    population = [] 
    tempCoordArray = []
  
    # create initial population 
    for _ in range(POPULATION_SIZE): 
                gnome = Individual.create_gnome()
                if DEBUG_MODE: print(gnome)
                for char in gnome:
                    #print(coordinates[char.upper()])
                    if tempCoordArray == []:
                        tempCoordArray.append(coordinates[char.upper()])
                    else:
                        tempCoordArray.append(tuple(map(lambda i, j: i + j, tempCoordArray[-1], coordinates[char.upper()])))
                if DEBUG_MODE: print(tempCoordArray)
                population.append(Individual(gnome,tempCoordArray))
                tempCoordArray = []
                
    population = sorted(population, key = lambda x:-x.fitness)
    if DEBUG_MODE == True:
        print("------------------")
        print(population[0].chromosome)
        print(population[0].coordinates)
        print(population[0].fitness)
    
    generate_protein_structure(population[0].coordinates)
    

def generate_protein_structure(coordinates):
    global HP_AMINO_ACID
    # Define amino acid colors and markers based on the HP amino acids
    aa_settings = {'h': {'color': 'red', 'marker': '^', 'label': 'Hydrophobic'},
                   'p': {'color': 'blue', 'marker': 'o', 'label': 'Polar'}}

    # Plot the 3D structure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i, aa_type in enumerate(HP_AMINO_ACID):
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


if __name__ == '__main__': 
    main() 