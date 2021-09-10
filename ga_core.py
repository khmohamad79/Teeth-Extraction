class Chromosome():
    def __init__(self, length):
        self.length = length


class Individual():
    def __init__(self, chromosome):
        self.chromosome = chromosome


class Population():
    def __init__(self, gen):
        self.gen = gen


class GeneticCore():
    def __init__(self, max_gen):
        self.max_gen = max_gen

