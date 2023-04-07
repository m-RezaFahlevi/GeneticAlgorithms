/*
 * Steady State Genetic Algorithm
 * Author: Muhammad Reza Fahlevi
 * Problem:
 * 	arg min spherefun(x, y)
 * 	s.t. x = (l1, u1); y = (l2, y2)
 *
 */
#include <bits/stdc++.h>
#include <chrono>
#define MU 0.0
#define STDEV 1.0
#define NPOPULATION 15
#define ALPHA 0.85
#define BETA 0.67
#define MAX_GENERATION 50
#define N_DIMENSION 2
using namespace std;

const int seed = chrono::system_clock::now().time_since_epoch().count();
mt19937 generator(seed);

normal_distribution<double> rnorm(MU, STDEV);

int runif(int at_least, int at_most) {
	uniform_int_distribution<int> randint(at_least, at_most);
	return randint(generator);
}

double runif(double at_least, double at_most) {
	uniform_real_distribution<double> randu(at_least, at_most);
	return randu(generator);
}

const double boundary[] = {-5.12, 5.12};

double spherefun(double x, double y) {
	return pow(x, 2.0) + pow(y, 2.0);
}

class Individual {
	vector<double> real_code;
	public:
		Individual();
		void print_gene();
		double fitness();
};

Individual::Individual() {
	for (int rcode = 0; rcode < N_DIMENSION; ++rcode) {
		double gen = runif(boundary[0], boundary[1]);
		real_code.push_back(gen);
	}
}

void Individual::print_gene() {
	for (double gene: real_code)
		cout << gene << " ";
}

double Individual::fitness() {
	double fitness_val = spherefun(real_code.at(0), real_code.at(1));
	return fitness_val;
}

/*
 * Population class have 3 private members :
 * individuals are member of population
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
		void print_population();
};

Population::Population() {
	double curr_worst = -99999.0;
	double curr_best = 99999.0;
	for (int indiv = 0; indiv < NPOPULATION; ++indiv) {
		Individual member_pop;
		if (member_pop.fitness() < curr_best) {
			curr_best = member_pop.fitness();
			best_indiv = member_pop;
		} if (member_pop.fitness() > curr_worst) {
			curr_worst = member_pop.fitness();
			worst_indiv = member_pop;
		}
		individuals.push_back(member_pop);
	}
}

void Population::print_population() {
	cout << "Individual \t\t finess\n";
	for (Individual creature: individuals) {
		creature.print_gene();
		cout << "\t " << creature.fitness() << "\n";
	}
	cout << "Worst individual : "; worst_indiv.print_gene();
	cout << "\nFitness value : " << worst_indiv.fitness() << endl;
	cout << "Best individual : "; best_indiv.print_gene();
	cout << "\nFitness value : " << best_indiv.fitness() << endl;
}

int main(void) {
	Individual individual;
	Population population;
	
	cout << rnorm(generator) << endl;
	cout << boundary[0] << endl;
	individual.print_gene();
	cout << individual.fitness() << "\n\n";

	population.print_population();
}
