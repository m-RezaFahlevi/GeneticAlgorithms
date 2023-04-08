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
int track_best = 0;
int track_worst = 0;

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
		Individual tournament(int t); // t >= 1 is size of tournament
		void update_population(Individual offspring);
		void print_population();
		void print_worst();
		void print_best();
};

Population::Population() {
	double curr_worst = -99999.0;
	double curr_best = 99999.0;
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

// int t >= 2 is total of tournament
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

void Population::print_population() {
	cout << " Individual \t\t finess\n";
	int counter = 0;
	for (Individual creature: individuals) {
		cout << counter << ". "; creature.print_gene();
		cout << "\t " << creature.fitness() << "\n";
		counter++;
	}
	cout << "Worst individual : "; worst_indiv.print_gene();
	cout << "\nFitness value : " << worst_indiv.fitness() << endl;
	cout << "Best individual : "; best_indiv.print_gene();
	cout << "\nFitness value : " << best_indiv.fitness() << endl;
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

int main(void) {
	Individual individual;
	Population population;
	Individual parent1 = population.tournament(2);
	Individual parent2 = population.tournament(2);

	cout << "Parent1 : "; parent1.print_gene();
	cout << "\nParent2 : "; parent2.print_gene();
	cout << endl;
	
	cout << rnorm(generator) << endl;
	cout << boundary[0] << endl;
	individual.print_gene();
	cout << individual.fitness() << "\n\n";

	population.print_population();
	population.print_worst();
	cout << "track best : " << track_best << endl;
	cout << "track worst : " << track_worst << endl;

	Individual achild;
	population.update_population(achild);

	cout << "\noffspring : "; achild.print_gene();
	cout << "\n";
	population.print_population();
	population.print_worst();
	cout << "track best : " << track_best << endl;
	cout << "track worst : " << track_worst << endl;
}
