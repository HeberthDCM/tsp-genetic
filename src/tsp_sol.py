# Python3 program to implement traveling salesman
# problem using naive approach.
from sys import maxsize
from itertools import permutations
import random
V = 4

# implementation of traveling Salesman Problem
def travellingSalesmanProblem(graph, s):

	# store all vertex apart from source vertex
	vertex = []
	for i in range(V):
		if i != s:
			vertex.append(i)

	# store minimum weight Hamiltonian Cycle
	min_path = maxsize
	next_permutation=permutations(vertex)
	for i in next_permutation:

		# store current Path weight(cost)
		current_pathweight = 0

		# compute current path weight
		k = s
		for j in i:
			current_pathweight += graph[k][j]
			k = j
		current_pathweight += graph[k][s]

		# update minimum
		min_path = min(min_path, current_pathweight)
	return min_path


# Driver Code
if __name__ == "__main__":

        # constant values
        cities = {0: 'Delhi', 1: 'Mumbai', 2: 'Pune', 3: 'Kolkata', 
                4: 'Gujrat', 5: 'Banglore', 6: 'J&K', 7: 'Nagaland'}

        # Distance between each pair of cities
        graph = [[100, 500, 450, 480, 200, 560, 180, 400],
        [500, 100, 80, 300, 250, 500, 590, 600],
        [450, 80, 100, 400, 226, 400, 500, 500],
        [480, 300, 400, 100, 500, 400, 500, 300],
        [200, 250, 226, 500, 100, 500, 400, 600],
        [560, 500, 400, 400, 500, 100, 600, 600],
        [180, 590, 500, 500, 400, 600, 100, 600],
        [400, 600, 500, 500, 600, 600, 600, 100]]
	
        # matrix representation of graph
        s = 0
        print(travellingSalesmanProblem(graph, s))

##### crossover
def single_point(gene_one, gene_two):
    p = random.randint(0, min(len(gene_one), len(gene_two)))
    gene_one[p: ], gene_two[p: ] = gene_two[p: ], gene_one[p: ]
    return gene_one, gene_two

def two_point(gene_one, gene_two):
    size = min(len(gene_one), len(gene_two))
    p1, p2 = random.randint(0, size), random.randint(0, size -1)
    if p1 > p2:
        p1, p2 = p2, p1
    else:
        p1, p2 = p1, p2 + 1
    gene_one[p1: p2], gene_two[p1: p2] = gene_two[p1: p2], gene_one[p1: p2]
    return gene_one, gene_two

def uniform(gene_one, gene_two, prob=0.4):
    for i in range(min(len(gene_one), len(gene_two))):
        if random.random() < prob:
            gene_one[i], gene_two[i] = gene_two[i], gene_one[i]
    return gene_one, gene_two

def partially_mapped(gene_one, gene_two):
    p1, p2 = random.randint(0, len(gene_one) - 1), random.randint(0, len(gene_one) - 1)
    if p1 >= p2:
        p1, p2 = p2, p1 + 1
    pos_one = {value: idx for idx, value in enumerate(gene_one)}
    pos_two = {value: idx for idx, value in enumerate(gene_two)}
    for i in range(p1, p2):
        value1, value2 = gene_one[i], gene_two[i]
        pos1, pos2 = pos_one[value2], pos_two[value1]
        gene_one[i], gene_one[pos1] = gene_one[pos1], gene_one[i]
        gene_two[i], gene_two[pos2] = gene_two[pos2], gene_two[i]
        pos_one[value1], pos_one[value2] = pos1, i
        pos_two[value1], pos_two[value2] = i, pos2

    return gene_one, gene_two

def cyclic(gene_one, gene_two):
    lookup = {v:i for i,v in enumerate(gene_one)}
    cycles = [-1] * len(gene_one)
    cyclestart = (i for i,v in enumerate(cycles) if v < 0)
    for cycle_no, pos in enumerate(cyclestart, 1):
        while cycles[pos] < 0:
            cycles[pos] = cycle_no
            pos = lookup[gene_two[pos]]

    for cycle in cycles:
        gene_one[cycle] = gene_two[gene_one[cycle]]
    return gene_one, gene_two
######evolution
class GeneticSolver:
    def __init__(self, fitness_func, crossover_func, mutation_func, pop_size, mut_prob, tourn_size, reverse=True):
        self._fitness_func, self._crossover_func, self._mutation_func = fitness_func, crossover_func, mutation_func
        self._pop_size, self._mut_prob, self._tourn_size, self._rev = pop_size, mut_prob, tourn_size, reverse

    def _tournament(self, pop, min_max):
        tourn = [pop[random.randint(0, len(pop) - 1)] for _ in range(self._tourn_size)]
        fitness = [self._fitness_func(gene) for gene in tourn]
        return pop[fitness.index(min_max(fitness))]

    def evolve(self, pop, gen_size):
        min_max = max if self._rev else min
        best_sol, best_fit = [], sys.maxsize

        for g in range(gen_size + 1):
            pop = sorted(pop, key= lambda gene: self._fitness_func(gene), reverse=self._rev)
            if (x := self._fitness_func(pop[0])) < best_fit:
                best_fit, best_sol = x, pop[0]

            print(f"Generation {g} :: {pop[0]}")
            print(f"Fitness :: {self._fitness_func(pop[0])} \n")

            new_pop = pop[0: 5]
            for i in range(int(len(pop) / 2) - 1):
                parent_one, parent_two = self._tournament(pop, min_max), self._tournament(pop, min_max)
                child_one, child_two = self._crossover_func(parent_one, parent_two)
                child_one, child_two = self._mutation_func(child_one, self._mut_prob), self._mutation_func(child_two, self._mut_prob)
                new_pop += [child_one, child_two]
            pop = new_pop

        print(f"Best so far :: {best_sol} {best_fit}")
        return sorted(pop, key= lambda gene: self._fitness_func(gene), reverse=self._rev)[0]

###########
def swap(gene, prob):
    if random.random() < prob:
        p1, p2 = random.randint(0, len(gene) - 2), random.randint(0, len(gene) - 1)
        if p1 >= p2:
            p1, p2 = p2, p1 + 1
        gene[p1], gene[p2] = gene[p2], gene[p1]
    return gene

def reverse(gene, prob):
    if random.random() < prob:
        p1, p2 = random.randint(0, len(gene) - 2), random.randint(0, len(gene) - 1)
        if p1 >= p2:
            p1, p2 = p2, p1 + 1
        gene[p1: p2] = gene[p1: p2][:: -1]
    return gene

def random_(gene, prob):
    if random.random() < prob:
        p1 = random.randint(0, len(gene) - 1)
        gene[p1] = random.choice(gene)
    return gene