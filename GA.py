# initial format is heavily inspired from geaks for geaks - https://www.geeksforgeeks.org/genetic-algorithms/
# based off of the paper 'An Enhanced Genetic Algorithm for Ab Initio Protein Structure Prediction' - http://www.cis.umassd.edu/~fkhatib/Papers/IEEE_TEC.pdf

"""
notes:
abcdefghijkl - each letter will represent a direction on a 3d plane eg (-1,0,1)
the target will be the actual sequance of amino acids e.g (hpphhphpph)
TODO:
need to edit the mutated_genes (in other words the chromosome creation) to have "Self-avoiding walk"
need to change fitness, fitness function, mutation and crossover
"""

import random
  
# Number of individuals in each generation 
POPULATION_SIZE = 10
  
# Valid genes 
GENES = 'abcdefghijkl'
  
# Target string to be generated 
TARGET = "abfghijkl"
  
class Individual(object): 
    ''' 
    Class representing individual in population 
    '''
    def __init__(self, chromosome,coordinates): 
        self.chromosome = chromosome
        self.coordinates = coordinates 
        self.fitness = self.cal_fitness() 
  
    @classmethod
    def mutated_genes(self): 
        ''' 
        create random genes for mutation 
        '''
        global GENES 
        gene = random.choice(GENES) 
        return gene 
  
    @classmethod
    def create_gnome(self): 
        ''' 
        create chromosome or string of genes 
        '''
        global TARGET 
        gnome_len = len(TARGET) 
        return [self.mutated_genes() for _ in range(gnome_len)] 
  
    def mate(self, par2): 
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
                child_chromosome.append(self.mutated_genes()) 
  
        # create new Individual(offspring) using  
        # generated chromosome for offspring 
        return Individual(child_chromosome) 
  
    def cal_fitness(self):
        print("TODO-fitnessArea") 
        ''' 
        Calculate fitness score, it is the number of 
        characters in string which differ from target 
        string. 
        '''
        global TARGET 
        fitness = 0
        for gs, gt in zip(self.chromosome, TARGET): 
            if gs != gt: fitness+= 1
        return fitness 
  
# Driver code 
def main(): 
    global POPULATION_SIZE 
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
  
    #current generation 
    generation = 1
  
    found = False
    population = [] 
    tempCoordArray = []
  
    # create initial population 
    for _ in range(POPULATION_SIZE): 
                gnome = Individual.create_gnome()
                print(gnome)
                for char in gnome:
                    #print(coordinates[char.upper()])
                    if tempCoordArray == []:
                        tempCoordArray.append(coordinates[char.upper()])
                    else:
                        tempCoordArray.append(tuple(map(lambda i, j: i + j, tempCoordArray[-1], coordinates[char.upper()])))
                print(tempCoordArray)
                population.append(Individual(gnome,tempCoordArray))
                tempCoordArray = []
                
    population = sorted(population, key = lambda x:x.fitness)
    print("------------------")
    print(population[0].chromosome)
    print(population[0].coordinates)
    print(population[0].fitness)
  
if __name__ == '__main__': 
    main() 