/* Real Code Genetic Algorithm for
 * Solving Real Parameter Optimization Problem
 *
 * Author	: Muhammad Reza Fahlevi
 * Email	: muhammadrezafahlevi666@gmail.com
 * GitHub	: https://github.com/m-RezaFahlevi
 * Dated	: 10th April 2023
 *
 * Problem:
 * 	arg min spherefun(x, y)
 * 	s.t. x = (l1, u1); y = (l2, y2)
 *
 * References:
 * [1] García-Martínez, C., Rodriguez, F.J., Lozano, M. (2018). 
 * Genetic Algorithms. In: Martí, R., Pardalos, P., Resende, M. (eds) Handbook of Heuristics. 
 * Springer, Cham. https://doi.org/10.1007/978-3-319-07124-4_28
 *
 * [2] Artificial Intelligence: A Modern Approach, Third Edition, ISBN 9780136042594, 
 * by Stuart J. Russell and Peter Norvig published by Pearson Education © 2010.
 *
 * [3] Sean Luke, 2013, Essentials of Metaheuristics, Lulu, second edition,
 * available at http://cs.gmu.edu/∼sean/book/metaheuristics/ 
 *
 * [4] Körner, T.W. 2013. Vectors, pure and applied : a general introduction to linear algebra. 
 * Cambridge. Cambridge University Press.
 *
 * [5] Peter Gottschling. 2022. Discovering Modern C++: An Intensive Course for 
 * Scientists, Engineers, and Programmers (second edition). Pearson Education, Inc.
 *
 * [6] Bjarne Stroustrup. 2018. A Tour of C++ (second edition). Pearson Education, Inc.
 *
 * [7] https://www.cplusplus.com/
 *
 */
#include <bits/stdc++.h>
#include <chrono>
#define NPOPULATION 15
#define NTOURNAMENT 2 // at least 2
#define ALPHA 0.85
#define MUTATION_PROBABILITY 0.78
#define BETA 0.67
#define DISP_EVOL false // boolean parameter
#define MAX_GENERATION 500
#define N_DIMENSION 2
using namespace std;

/* Set seed that will be used as a generator in
 * pseudo random number generator (PRNG) function.
 */
const int seed = chrono::system_clock::now().time_since_epoch().count();
mt19937 generator(seed);

/* runif stand for (discrete) random uniform, it is
 * PRNG that return an integer number range from
 * at_least to at_most.
 */
int runif(int at_least, int at_most) {
	uniform_int_distribution<int> randint(at_least, at_most);
	return randint(generator);
}

/* overloading runif, runif stand for (continuous) random
 * uniform, it is PRNG that return a real number range
 * from at_least to at_most.
 */
double runif(double at_least, double at_most) {
	uniform_real_distribution<double> randu(at_least, at_most);
	return randu(generator);
}

const double boundary[] = {-5.12, 5.12};

/* spherefun stand for sphere function that map from
 * R^2 to R s.t. f(x,y) = x^2 + y^2
 */
double spherefun(double x, double y) {
	return pow(x, 2.0) + pow(y, 2.0);
}

class Individual {
	vector<double> real_code;
	public:
		Individual(); // default constructor
		Individual(vector<double>);
		void print_gene();
		double fitness();
		vector<double> export_gene();
};

// default constructor
Individual::Individual() {
	for (int rcode = 0; rcode < N_DIMENSION; ++rcode) {
		double gen = runif(boundary[0], boundary[1]);
		real_code.push_back(gen);
	}
}

Individual::Individual(vector<double> genes) {
	real_code = genes;
}

void Individual::print_gene() {
	cout << "(";
	for (double gene: real_code)
		cout << gene << ", ";
	cout << "\b\b \b)";
}

double Individual::fitness() {
	double fitness_val = spherefun(real_code.at(0), real_code.at(1));
	return fitness_val;
}

vector<double> Individual::export_gene() {
	return real_code;
}

/* Define global variables track_best and track_worst
 * as an index to track the worst and the best individual in
 * a population for each generation t = 1, 2, ... , MAX_GENERATION.
 */
int track_best = 0;
int track_worst = 0;

/* Population class have 3 private members :
 * individuals are member of population.
 * worst_indiv (worst individual), in arg max fun problem the worst in-
 * dividual is individual for which has the lowest fitness value. In
 * arg min fun problem, the worst individual is individual for which has
 * the highest fitness value.
 * best_indiv (best individual), the opposite of worst individual.
 */
class Population {
	vector<Individual> individuals;
	Individual worst_indiv, best_indiv;
	public:
		Population();
		Individual tournament(int t); // t >= 2 is size of tournament
		void update_population(Individual offspring);
		Individual export_worst();
		Individual export_best();
		void print_population();
		void print_worst();
		void print_best();
};

Population::Population() {
	double curr_worst = -numeric_limits<double>::infinity();
	double curr_best = numeric_limits<double>::infinity();
	int search_worst = 0;
	int search_best = 0;
	for (int indiv = 0; indiv < NPOPULATION; ++indiv) {
		Individual member_pop; // default constructor;
		if (member_pop.fitness() < curr_best) {
			curr_best = member_pop.fitness();
			best_indiv = member_pop;
			track_best = search_best;
		} if (member_pop.fitness() > curr_worst) {
			curr_worst = member_pop.fitness();
			worst_indiv = member_pop;
			track_worst = search_worst;
		}
		individuals.push_back(member_pop);
		search_worst++;
		search_best++;
	}
}

// int t >= 2 is total tournament
Individual Population::tournament(int t) {
	int get_num = runif(0, NPOPULATION - 1);
	Individual curr_suvivor = individuals.at(get_num);
	for (int nth_tournament = 2; nth_tournament <= t; ++nth_tournament) {
		get_num = runif(0, NPOPULATION - 1);
		Individual gladiator = individuals.at(get_num);
		if (gladiator.fitness() < curr_suvivor.fitness())
			curr_suvivor = gladiator;
	}
	return curr_suvivor;
}

void Population::update_population(Individual offspring) {
	individuals.at(track_worst) = offspring;
	worst_indiv = offspring;
	// track best and worst individual (linear search), O(n)
	int find_worst = 0;
	int find_best = 0;
	for (int i = 0; i < NPOPULATION; ++i) {
		Individual scanned_indiv = individuals.at(i);
		if (scanned_indiv.fitness() < best_indiv.fitness()) {
			best_indiv = scanned_indiv;
			track_best = find_best;
		} if (scanned_indiv.fitness() > worst_indiv.fitness()) {
			worst_indiv = scanned_indiv;
			track_worst = find_worst;
		}
		find_worst++;
		find_best++;
	}
}

Individual Population::export_best() {
	return best_indiv;
}

Individual Population::export_worst() {
	return worst_indiv;
}

void Population::print_population() {
	int counter = 0;
	for (Individual creature: individuals) {
		cout << counter + 1 << ". f"; creature.print_gene();
		cout << " = " << creature.fitness() << "\n";
		counter++;
	}
	cout << "Worst individual : " << track_worst + 1 << ". f"; worst_indiv.print_gene();
	cout << " = " << worst_indiv.fitness() << endl;
	cout << "Best individual : " << track_best + 1 << ". f"; best_indiv.print_gene();
	cout << " = " << best_indiv.fitness() << endl;
}

void Population::print_worst() {
	Individual the_worst = individuals.at(track_worst);
	cout << "Worst individual: "; the_worst.print_gene();
	cout << "\nFitness value : " << the_worst.fitness() << endl;
}

void Population::print_best() {
	Individual the_best = individuals.at(track_best);
	cout << "Best individual: "; the_best.print_gene();
	cout << "\nFitness value: " << the_best.fitness() << endl;
}

/* blx_alpha recombination have 3 parameters:
 * gene1 is gene from parent 1.
 * gene2 is gene from parent 2.
 * conf_beta stand for configured beta, this 
 * parameter need to be tuned, its value range from 0.0 to 1.0.
 */
double blx_alpha(double gene1, double gene2, double conf_alpha) {
	double child_gene = gene1 - gene2;
	child_gene = abs(child_gene);
	double l_bound = max(boundary[0], min(gene1, gene2) - (conf_alpha * child_gene));
	double u_bound = min(boundary[1], max(gene1, gene2) + conf_alpha * child_gene);
	child_gene = runif(l_bound, u_bound);
	return child_gene;
}

/* nonunif_mutation stand for nonuniform mutation, this function
 * take 3 parameters:
 * the_gene is gene that will be mutated.
 * curr_generation stand for current generation.
 * conf_beta stand for configured beta, this parameter need to
 * be tuned, its value range from 0.0 to 1.0.
 */
double nonunif_mutation(double the_gene, int curr_generation, double conf_beta) {
	double comp_gamma = pow(1 - curr_generation / MAX_GENERATION, conf_beta);
	double p_mm = runif(0.0, 1.0);
	double k = p_mm <= 0.5 ? (boundary[1] - the_gene) : (boundary[0] - the_gene);
	double comp_delta = pow(k * (1 - runif(0.0, 1.0)), comp_gamma);
	the_gene += comp_delta;
	return the_gene;
}

/* steady_state_rcga contain 2 parameters that need
 * to be tuned:
 * p_mutation is the probability that mutation will occur.
 * disp_evol is boolean, if disp_evol is true, then print the population
 * for each generation. If disp_evol is false, then only print population of
 * the the first and the last generation.
 */
Individual steady_state_rcga(double p_mutation, bool disp_evol) {
	Population population; // initial population
	cout << "Initial population\n";
	population.print_population();
	// evolution
	for (int nth_generation = 0; nth_generation < MAX_GENERATION; ++nth_generation) {
		if (disp_evol) {
			cout << nth_generation << "th generation: \n";
			population.print_population();
			cout << endl;
		}
		for (int epoch = 0; epoch < NPOPULATION; ++epoch) {
			// tournament selection
			vector<double> parent1 = population.tournament(NTOURNAMENT).export_gene();
			vector<double> parent2 = population.tournament(NTOURNAMENT).export_gene();
			vector<double> offspring;
			for (int i = 0; i < N_DIMENSION; ++i) {
				// BLX-alpha recombination
				double offspring_gene = blx_alpha(parent1.at(i), parent2.at(i), ALPHA);
				// nonuniform mutation, 
				double p_mut = runif(0.0, 1.0);
				if (p_mut <= p_mutation) {
					offspring_gene = nonunif_mutation(offspring_gene, nth_generation, BETA);
				}
				offspring.push_back(offspring_gene);
			}
			Individual an_offspring(offspring);
			if (an_offspring.fitness() < population.export_worst().fitness()) {
				population.update_population(an_offspring);
			}
		}
	}
	cout << endl;
	cout << MAX_GENERATION << "th generation:\n";
	population.print_population();
	return population.export_best();
}

void prelude() {
	cout << "Real code genetic algorithm to solve arg min spherefun(x,y) = x^2 + y^2\n";
	cout << "Population\t: there are " << NPOPULATION << " individual\n";
	cout << "Selection\t: tournament\n";
	cout << "\t\t  tuned parameter NTOURNAMENT = " << NTOURNAMENT << endl;
	cout << "Recombination\t: BLX-alpha\n";
	cout << "\t\t  tuned parameter ALPHA = " << ALPHA << endl;
	cout << "Mutation\t: non-uniform mutation\n";
	cout << "\t\t  tuned parameter BETA = " << BETA << endl;
	cout << "Evolution model\t: steady state\n";
	cout << "Stop criterion\t: fixed number of iteration\n";
	cout << "\t\t  tuned parameter MAX_GENERATION = " << MAX_GENERATION << "\n\n";
}

int main(void) {
	prelude();
	Individual init_sol;
	cout << "Initial solution : "; init_sol.print_gene();
	cout << "\nFitness value : " << init_sol.fitness() << endl;

	Individual obt_sol = steady_state_rcga(MUTATION_PROBABILITY, DISP_EVOL);
	cout << "\nObtained solution : "; obt_sol.print_gene();
	cout << "\nFitness value : " << obt_sol.fitness() << endl;
}
