import random


class Gene():
    def __init__(self, top, center, bottom):
        self.top = top
        self.center = center
        self.bottom = bottom
        self.points = None

    def shake(self, max_move, max_width):
        new_top = self.top + random.randint(-1*max_move, max_move)
        if new_top > max_width: new_top = max_width
        if new_top < 1: new_top = 1
        self.top = new_top

        new_center = self.center + random.randint(-1*max_move, max_move)
        if new_center > max_width: new_center = max_width
        if new_center < 1: new_center = 1
        self.center = new_center

        new_bottom = self.bottom + random.randint(-1*max_move, max_move)
        if new_bottom > max_width: new_bottom = max_width
        if new_bottom < 1: new_bottom = 1
        self.bottom = new_bottom

        self.points = None

    def get_points(self, max_height):
        if self.points is not None:
            return self.points
        
        self.points = []

        horizon = (max_height + 1) // 2
        above_horizon = range(1, horizon)
        under_horizon = range(horizon + 1, max_height + 1)

        for y in above_horizon:
            x = int(self.center - (self.center - self.top) * (horizon - y) / (horizon - 1))
            self.points.append((x, y))

        self.points.append(((self.center), horizon))

        for y in under_horizon:
            x = int(self.bottom - (self.bottom - self.center) * (max_height - y) / (max_height - horizon))
            self.points.append((x, y))
        
        return self.points

    def min_points(self):
        return min(self.top, self.center, self.bottom)

    def max_points(self):
        return max(self.top, self.center, self.bottom)


class Chromosome():
    def __init__(self, genes=None):
        self.genes = genes if genes is not None else []

    def add_gene(self, gene):
        self.genes.append(gene)

    def get_cost(self, image):
        max_height, _ = image.shape

        cost = 0
        for gene in self.genes:
            total_intensity = 0
            count = 0
            
            for point in gene.get_points(max_height):
                x, y = point
                intensity = image[y-1, x-1]
                if intensity == 0:
                    continue
                
                total_intensity += intensity
                count += 1

            if count > 0:
                cost += total_intensity / count

        return cost

    def length(self):
        return len(self.genes)

    def mutate(self, max_width):
        random_gene = random.choice(self.genes)
        random_gene.shake(10, max_width)

    def init_random(length, max_width):
        chromosome = Chromosome()
        range_widths = length * [0]
        range_index = 0
        remaining_pixels = max_width
        while remaining_pixels > 0:
            remaining_pixels -= 1
            range_widths[range_index] += 1
            range_index = (range_index + 1) % length

        vertical_offset = 1
        for range_width in range_widths:
            points_in_range = range(vertical_offset, vertical_offset + range_width)
            
            random_top = random.choice(points_in_range)
            random_center = random.choice(points_in_range)
            random_bottom = random.choice(points_in_range)
            
            chromosome.add_gene(Gene(random_top, random_center, random_bottom))
            
            vertical_offset += range_width

        return chromosome


class Individual():
    def __init__(self, chromosome, parent1=None, parent2=None):
        self.chromosome = chromosome
        self.parent1 = parent1
        self.parent2 = parent2

    def mutate(self, max_width):
        self.chromosome.mutate(max_width)
    
    def get_cost(self, image):
        return self.chromosome.get_cost(image)

    def init_random(chromosome_length, max_width):
        return Individual(Chromosome.init_random(chromosome_length, max_width))

    def init_crossover_onepoint(parent1, parent2):
        if parent1.chromosome.length() != parent2.chromosome.length():
            raise Exception('GA:EXCEPTION >>> parents must have chromosomes with the same length')
        
        chromosome_length = parent1.chromosome.length()

        cross_point = random.randint(1, chromosome_length)

        child1 = Individual(Chromosome(parent1.chromosome.genes[:cross_point] + parent2.chromosome.genes[cross_point:]), parent1, parent2)
        child2 = Individual(Chromosome(parent2.chromosome.genes[:cross_point] + parent1.chromosome.genes[cross_point:]), parent1, parent2)

        return child1, child2

    def init_crossover_onepoint_best_parts(parent1, parent2, image):
        if parent1.chromosome.length() != parent2.chromosome.length():
            raise Exception('GA:EXCEPTION >>> parents must have chromosomes with the same length')
        
        chromosome_length = parent1.chromosome.length()

        cross_point = random.randint(1, chromosome_length)

        left_genes_parent_1 = parent1.chromosome.genes[:cross_point]
        left_genes_parent_2 = parent2.chromosome.genes[:cross_point]
        right_genes_parent_1 = parent1.chromosome.genes[cross_point:]
        right_genes_parent_2 = parent2.chromosome.genes[cross_point:]

        left_1_cost = Chromosome(left_genes_parent_1).get_cost(image)
        left_2_cost = Chromosome(left_genes_parent_2).get_cost(image)
        right_1_cost = Chromosome(right_genes_parent_1).get_cost(image)
        right_2_cost = Chromosome(right_genes_parent_2).get_cost(image)

        left_genes = left_genes_parent_1 if left_1_cost < left_2_cost else left_genes_parent_2
        right_genes = right_genes_parent_1 if right_1_cost < right_2_cost else right_genes_parent_2

        return Individual(Chromosome(left_genes + right_genes), parent1, parent2)


class Population():
    def __init__(self, gen):
        self.gen = gen
        self.individuals = []

    def add_individual(self, individual):
        self.individuals.append(individual)

    def add_individuals(self, individuals):
        self.individuals.extend(individuals)

    def validate_size(self, required_size):
        if len(self.individuals) < required_size:
            raise Exception('GA:EXCEPTION >>> lack of individuals int the population')
        
        self.individuals = self.individuals[:required_size]

    def get_random_individual(self):
        return random.choice(self.individuals)

    def get_half_top(self, image):
        costs_individuals = []
        for individual in self.individuals:
            costs_individuals.append((individual.get_cost(image), individual))
        costs_individuals.sort(key= lambda t: t[0])
        
        return [t[1] for t in costs_individuals[:(len(costs_individuals) + 1) // 2]]

    def get_best_individual(self, image):
        best_individual = None
        best_cost = None

        for individual in self.individuals:
            cost = individual.get_cost(image)
            if best_individual is None or cost < best_cost:
                best_individual = individual
                best_cost = cost

        return best_individual

    def init_random(pop_size, chromosome_length, max_width):
        pop = Population(0)
        for _ in range(pop_size):
            pop.add_individual(Individual.init_random(chromosome_length, max_width))

        return pop


class GeneticCore():
    def __init__(self, max_gen, pop_size, chromosome_length, mutation_probability, image):
        self.max_gen = max_gen
        self.pop_size = pop_size
        self.chromosome_length = chromosome_length
        self.mutation_probability = mutation_probability
        self.image = image
        self.populations = []

    def run(self):
        # initialize first generation with random population
        current_pop = Population.init_random(self.pop_size, self.chromosome_length, self.image.shape[1])
        self.populations.append(current_pop)
       
        best_individual = current_pop.get_best_individual(self.image)
        best_cost = current_pop.get_best_individual(self.image).get_cost(self.image)
        best_gen = current_pop.gen

        # log
        print(f'GA:INFO >>> gen:{current_pop.gen} cost:{best_cost:.2f} best_cost:{best_cost:.2f} best_gen:{best_gen}')

        # simulate evolution
        while current_pop.gen < self.max_gen:
            current_half_top = current_pop.get_half_top(self.image)
            
            # next generation population
            next_pop = Population(current_pop.gen + 1)
            
            # generating half of the next population
            for _ in range(self.pop_size // 2 + 1):
                # select parents
                parent1 = random.choice(current_half_top)
                parent2 = random.choice(current_half_top)
                
                # crossover parents
                ### child1, child2 = Individual.init_crossover_onepoint(parent1, parent2)
                child = Individual.init_crossover_onepoint_best_parts(parent1, parent2, self.image)

                # mutate newly generated individual
                ### if random.random() < self.mutation_probability: child1.mutate(self.image.shape[1])
                ### if random.random() < self.mutation_probability: child2.mutate(self.image.shape[1])
                if random.random() < self.mutation_probability: child.mutate(self.image.shape[1])

                # add to new population
                ### next_pop.add_individual(child1)
                ### next_pop.add_individual(child2)
                next_pop.add_individual(child)


            # forwarding half top of the current generation to the next population
            next_pop.add_individuals(current_half_top)

            # discard extra individuals
            next_pop.validate_size(self.pop_size)

            # proceed to next generation
            self.populations.append(next_pop)
            current_pop = next_pop

            # remember best answer
            cost = current_pop.get_best_individual(self.image).get_cost(self.image)
            if cost < best_cost:
                best_individual = current_pop.get_best_individual(self.image)
                best_cost = cost
                best_gen = current_pop.gen

            # log
            print(f'GA:INFO >>> gen:{current_pop.gen} cost:{cost:.2f} best_cost:{best_cost:.2f} best_gen:{best_gen}')


        return best_individual.chromosome