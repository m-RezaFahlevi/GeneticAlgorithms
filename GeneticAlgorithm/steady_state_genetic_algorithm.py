"""
Author      :   Muhammad Reza Fahlevi
Dated       :   June 7th, 2021
Affiliation :   Departemen Ilmu Komputer,
                Fakultas Ilmu Komputer - Teknologi Informasi,
                Universitas Sumatera Utara, Indonesia
References  :   Benjamin Baka. Python Data Structures and Algorithms. 2017. Packt Publishing.
                García-Martínez C., Rodriguez F.J., Lozano M. (2018) Genetic Algorithms. In: Martí R., Pardalos P., Resende M. (eds) Handbook of Heuristics. Springer, Cham. https://doi.org/10.1007/978-3-319-07124-4_28
                https://docs.python.org/3.8/tutorial/index.html
                J M García et al 2020 J. Phys.: Conf. Ser. 1448 012020
                Stuart Russell and Peter Norvig. 2020. Artificial Intelligence: A Modern Approach, 4th Edition. University of California at Berkeley. Pearson.
"""

# Prepare the tools
import math
import pandas as pd
from numpy import random

# phi = lambda x,y: math.exp(-0.5 * ((x ** 2) + (y ** 2)))

def phi(x, y):
    first_denominator = 1 + ((x - 1) ** 2) + ((y - 1) ** 2)
    second_denominator = 1 + (1/4) * (((x + 0.5) ** 2)) + (1/36) * ((y - 1) ** 2)
    third_denominator = 1 + ((x - 2) ** 2) + ((y - 2) ** 2)
    fourth_denominator = 1 + ((x - 1) ** 2) + ((y + 1) ** 2)
    exp_term = (3 / 2) * math.exp(1 / first_denominator) + (-5 / 2) * math.exp(1 / second_denominator) + 2 * math.exp(1 / third_denominator) + 2 * math.exp(1 / fourth_denominator)
    return exp_term

xvect, yvect = random.random(10), random.random(10)
print(xvect, yvect)
print([phi(x, y) for x, y in zip(xvect, yvect)])

def competition(df_pop):
    st_candidate = df_pop.iloc[random.randint(9), 0:3]
    nd_candidate = df_pop.iloc[random.randint(9), 0:3]
    winner = st_candidate if st_candidate['fitness'] > nd_candidate['fitness'] else nd_candidate
    return winner

def pbxalpha(domain_interval, p1, p2, alpha):
    constant = abs(p1['xgenotype'] - p2['xgenotype'])
    p1_lower_bound = max(domain_interval[0], p1['xgenotype'] - (alpha * constant))
    p1_upper_bound = min(domain_interval[1], p1['xgenotype'] + (alpha* constant))
    
    constant = abs(p1['ygenotype'] - p2['ygenotype'])
    p2_lower_bound = max(domain_interval[0], p2['ygenotype'] - (alpha * constant))
    p2_upper_bound = min(domain_interval[1], p2['ygenotype'] + (alpha* constant))
    
    p1_offspring = [p1_lower_bound, p1_upper_bound]
    p2_offspring = [p2_lower_bound, p2_upper_bound]
    
    return p1_offspring, p2_offspring

gamma = lambda t, tmax, beta: (1 - (t / tmax)) ** beta

def delta(domain_interval, genotype, computed_gamma):
    get_probability = random.random()
    probability_criteria = 1 / 2
    normvar = random.random()
    if get_probability == probability_criteria:
        delta_value = (domain_interval[1] - genotype) * ((1 - normvar) ** computed_gamma)
    else:
        delta_value = (domain_interval[0] - genotype) * ((1 - normvar) ** computed_gamma)
    return delta_value

# Define the constant
NP = 10
Pm = 0.3
MAX_GENERATION = 50
problem_interval = [-5, 5]

def config_df(x_vect, y_vect):
    dframe = pd.DataFrame(
        {
            "xgenotype": x_vect,
            "ygenotype": y_vect,
            "fitness": [phi(x, y) for x, y in zip(x_vect, y_vect)]
        }, index = [(idx + 1) for idx in range(len(x_vect))]
    )
    data_frame = dframe.sort_values('fitness')
    data_frame = data_frame.iloc[1:10]
    return data_frame    

def steady_state_genetic_algorithm(population_size, mutation_probability, max_generation):
    generations = []
    x_vector = random.random(population_size)
    y_vector = random.random(population_size)
    
    for generation in range(max_generation):
        df_population = config_df(x_vector, y_vector)
        generations.append(df_population)
        x_vector = [] # start a new population
        y_vector = []
        x_vector.append(df_population.iloc[population_size - 2, 0]) # add the fittest to the next generation
        y_vector.append(df_population.iloc[population_size - 2, 1])
        vect_size = len(x_vector)
        while vect_size != population_size:
            p1, p2 = competition(df_population), competition(df_population) # Selection
            offsprings = pbxalpha(problem_interval, p1, p2, alpha=0.5) # Crossover, return tuple of list
            current_mutation_rate = random.random()
            if current_mutation_rate < mutation_probability:
                p11 = offsprings[0][0] + delta(problem_interval, offsprings[0][0], gamma(generation, max_generation, beta=9))
                p12 = offsprings[0][1] + delta(problem_interval, offsprings[0][1], gamma(generation, max_generation, beta=9))
                p21 = offsprings[1][0] + delta(problem_interval, offsprings[1][0], gamma(generation, max_generation, beta=9))
                p22 = offsprings[1][1] + delta(problem_interval, offsprings[1][1], gamma(generation, max_generation, beta=9))
                offsprings = tuple([[p11, p12], [p21, p22]]) # mutation occur, return tuple of list
            else:
                offsprings = offsprings # no mutation occurs
            # Evaluation
            first_value = phi(offsprings[0][0], offsprings[0][1])
            second_value = phi(offsprings[1][0], offsprings[1][1])
            survivor = [offsprings[0][0], offsprings[0][1]] if first_value > second_value else [offsprings[1][0], offsprings[1][1]]
            x_vector.append(survivor[0])
            y_vector.append(survivor[1])
            vect_size += 1
    return generations

for i in range(25):
    print(steady_state_genetic_algorithm(NP, Pm, MAX_GENERATION)[MAX_GENERATION - 1])
