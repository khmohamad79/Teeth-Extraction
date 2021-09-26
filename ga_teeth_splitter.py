import cv2
import time
import copy
import numpy as np
import os
from ga_core import GeneticCore
from preprocessing import CLAHE

generations = 50
population = 100
chromosome_length = 30 # at least 16 # best number is 30
mutation_probability = 0.02



def gene_avg_intensity(gene, image, max_height):
    total_intensity = 0
    count = 0
    for point in gene.get_points(max_height):
        x, y = point
        intensity = image[y-1, x-1]
        if intensity == 0:
            continue
        
        total_intensity += intensity
        count += 1

    return total_intensity // count

def mean_avg_intensity(genes, image, max_height):
    costs = []
    for gene in genes:
        cost = gene_avg_intensity(gene, image, max_height)
        costs.append(cost)

    return np.mean(costs)



for i in range(1, 4):
    print('')
    print('========================= number %d ==========================' % i)
    
    print('loading upper jaw number %d' % i)
    img_address = './upper_jaws/%d_upper_clahe_sauvola.bmp' % i
    img = cv2.imread(img_address, 0)
    final_image = copy.deepcopy(img)
    enhanced_image = CLAHE(final_image)

    print('applying genetic algorithm in order to find best tooth separator lines')
    t0 = time.time()
    ga = GeneticCore(generations, population, chromosome_length, mutation_probability, enhanced_image)
    chromosome = ga.run()
    t1 = time.time()
    print('elapsed time for GA: %.2f secs' % (t1 - t0))

    print('drawing lines')
    final_image = cv2.cvtColor(final_image, cv2.COLOR_GRAY2RGB)

    genes_avg_intensities = []
    for gene in chromosome.genes:
        genes_avg_intensities.append(gene_avg_intensity(gene, enhanced_image, enhanced_image.shape[0]))

    regions = []

    for j in range(len(chromosome.genes)):
        gene = chromosome.genes[j]
        
        color = [255, 0, 0]
        if j > 0 and j < (len(chromosome.genes) - 1):
            if genes_avg_intensities[j] > genes_avg_intensities[j-1] and genes_avg_intensities[j] > genes_avg_intensities[j+1]:
                color = [0, 0, 127]
            else:
                regions.append((previous_gene.min_points(), gene.max_points()))
                previous_gene = gene
                previous_gene_j = j
        elif j == 0:
            regions.append((1, gene.max_points()))
            previous_gene = gene
            previous_gene_j = j
        elif j == (len(chromosome.genes) - 1):
            regions.append((gene.min_points(), final_image.shape[1]))

        for point in gene.get_points(final_image.shape[0]):
            x, y = point
            final_image[y-1, x-1] = color

    # Export teeth
    if not os.path.exists('./upper_jaws/%d' % i):
        os.mkdir('./upper_jaws/%d' % i)
    y1, y2 = 0, final_image.shape[0]-1
    number = 1
    for region in regions:
        x1, x2 = region[0]-1, region[1]-1
        cropped = img[y1:y2+1, x1:x2+1]
        cv2.imwrite(f'./upper_jaws/{i}/{number}.bmp', cropped)
        number += 1

    print('saving result image')
    cv2.imwrite('./upper_jaws/%d_upper_clahe_sauvola_result.bmp' % i, final_image)

    #####################

    print('loading lower jaw number %d' % i)
    img_address = './lower_jaws/%d_lower_clahe_sauvola.bmp' % i
    img = cv2.imread(img_address, 0)
    final_image = copy.deepcopy(img)
    enhanced_image = CLAHE(final_image)

    print('applying genetic algorithm in order to find best tooth separator lines')
    t0 = time.time()
    ga = GeneticCore(generations, population, chromosome_length, mutation_probability, enhanced_image)
    chromosome = ga.run()
    t1 = time.time()
    print('elapsed time for GA: %.2f secs' % (t1 - t0))

    mean_cost = mean_avg_intensity(chromosome.genes, enhanced_image, enhanced_image.shape[0])

    print('drawing lines')
    final_image = cv2.cvtColor(final_image, cv2.COLOR_GRAY2RGB)
    
    genes_avg_intensities = []
    for gene in chromosome.genes:
        genes_avg_intensities.append(gene_avg_intensity(gene, enhanced_image, enhanced_image.shape[0]))

    regions = []

    for j in range(len(chromosome.genes)):
        gene = chromosome.genes[j]
        
        color = [255, 0, 0]
        if j > 0 and j < (len(chromosome.genes) - 1):
            if genes_avg_intensities[j] > genes_avg_intensities[j-1] and genes_avg_intensities[j] > genes_avg_intensities[j+1]:
                color = [0, 0, 127]
            else:
                regions.append((previous_gene.min_points(), gene.max_points()))
                previous_gene = gene
                previous_gene_j = j
        elif j == 0:
            regions.append((1, gene.max_points()))
            previous_gene = gene
            previous_gene_j = j
        elif j == (len(chromosome.genes) - 1):
            regions.append((gene.min_points(), final_image.shape[1]))

        for point in gene.get_points(final_image.shape[0]):
            x, y = point
            final_image[y-1, x-1] = color

    # Export teeth
    if not os.path.exists('./lower_jaws/%d' % i):
        os.mkdir('./lower_jaws/%d' % i)
    y1, y2 = 0, final_image.shape[0]-1
    number = 1
    for region in regions:
        x1, x2 = region[0]-1, region[1]-1
        cropped = img[y1:y2+1, x1:x2+1]
        cv2.imwrite(f'./lower_jaws/{i}/{number}.bmp', cropped)
        number += 1

    print('saving result image')
    cv2.imwrite('./lower_jaws/%d_lower_clahe_sauvola_result.bmp' % i, final_image)


print('process finished')