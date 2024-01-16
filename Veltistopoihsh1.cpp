#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>

using namespace std;

double objectiveFunction(const vector<double>& x);
vector<int> tournamentSelection(const vector<double>& fitness, int tournamentSize);
void mutation(vector<double>& offspring, int iter);
void gradientDescent(vector<double>& point, int iter);
void printPopulation(const char* m,vector<vector<double>> p);

bool debug = false;
int N = 10;//number of chromosomes
int cSize=3;//size of chromosomes
int m=4;//number of offsprings to be created(m<=N)
double mProb = 0.4;//mutation probability(initialize from 0 to 1)
int ITERMAX = 200;//max iterations
double gradientStep=0.4;
const int tournamentSize=5;

double upperBound=1.0;
double lowerBound=0.0;


int main() {
	
	srand (time(nullptr));

// Step 1: Initialization
	vector<vector<double>> S (N, vector<double>(cSize));
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < cSize; ++j) {
			S[i][j] = lowerBound + static_cast<double>(rand()) / RAND_MAX * (upperBound - lowerBound);
		}
	}
	if(debug)
		printPopulation("Original population:",S);
	
	int iter = 0;
//the loop
	while(iter<ITERMAX) {
// Step 2: Evaluation
		vector<double> fitness(N);
		for (int i = 0; i < N; ++i) {
			fitness[i] = objectiveFunction(S[i]);
		}
// Step 3: Termination check
		if ((fitness[N-1] - fitness[0]) <= 1e-6) {
			break;
		}
// Step 4: Genetic operations
		vector<vector<double>> offspring(N, vector<double>(cSize));
		for (int i = 0; i < m; i+=2) {
		// Selection
			vector<int> selectedIndices = tournamentSelection(fitness, tournamentSize);
			int p1 = selectedIndices[0];
			int p2 = selectedIndices[1];
		// Crossover
			double a = -0.5 + static_cast<double>(rand()) / RAND_MAX * 2.0;
			for(int j=0;j<cSize;j++) {
				offspring[i][j]=a*S[p1][j]+(1-a)*S[p2][j];
				offspring[i+1][j]=a*S[p2][j]+(1-a)*S[p1][j];
			}
		}
		// Mutation
		for (int i = 0; i < N; i++) {
			double r= 0 + static_cast<double>(rand()) / RAND_MAX * (1 - 0);
			if (r < mProb) {
				mutation(offspring[i], iter);
			}
		}
// Step 5: Replacement
		for (int i = 0; i < m; ++i) {
			for(int j=0;j<cSize;j++) {
				S[i][j] = offspring[i][j];
			}
		}
//Step 6 (local technique)(Gradient descent):
		for (int i = 0; i < N; ++i) {
			// Local search using gradient descent
			gradientDescent(S[i], iter);
			// Evaluate the objective function for the updated solution
			double updatedFitness = objectiveFunction(S[i]);
			// Replace the original solution if the updated one is better
			if (updatedFitness < fitness[i]) {
				fitness[i] = updatedFitness;
			}
		}
		
		iter++;
	}
	
	if(debug)
		printPopulation("\n\nFinal population:",S);
	
	cout << "Best solution found: " << S[0][0] << " " << S[0][1] << endl <<  "Objective function value: " << objectiveFunction(S[0]) << endl;

	return 0;
}


double objectiveFunction(const vector<double>& x) {
	//return 1.0 / (1.0 + x[0]*x[0] + x[1]*x[1]);
	//Exponential function
	double sum = 0.0;
	for (double xi : x) {
		sum += xi * xi;
	}
	return exp(-sum);
}


vector<int> tournamentSelection(const vector<double>& fitness, int tournamentSize) {
	vector<int> selectedIndices;

	for (int i = 0; i < 2; ++i) {
		int bestIndex = rand() % fitness.size();
		for (int j = 1; j < tournamentSize; ++j) {
			int candidateIndex = rand() % fitness.size();
			if (fitness[candidateIndex] < fitness[bestIndex]) {
				bestIndex = candidateIndex;
			}
		}
		selectedIndices.push_back(bestIndex);
	}

	return selectedIndices;
}


double calcD(int iter,double y) {
	int r=rand() % 2;
	return y*(1.0-pow(r, 1.0 - static_cast<double>(iter) / ITERMAX));
}

void mutation(vector<double>& offspring, int iter) {
	//using
	//x'i = {xi + D(iter, ri-xi), t=0,
	//		{xi - D(iter, xi-li), t=1,
	//where t = 0 or 1, r = random number in [0,1]
	//and D(i,y)=y*(1-r^(1-i/iMax))
	for (int i = 0; i < cSize; ++i) {
		int t = rand() % 2;
		double r= 0 + static_cast<double>(rand()) / RAND_MAX * (1 - 0);
		if(t==0) {
			offspring[i] = offspring[i] + calcD(iter,r - offspring[i]);
		}else {
			offspring[i] = offspring[i] - calcD(iter,offspring[i] - offspring[0]);
		}
	}
}


void gradientDescent(vector<double>& point, int iter) {
	vector<double> gradient(point.size(), 0.0);

	for (int i = 0; i < point.size(); ++i) {
		// Compute the derivative of the objective function with respect to the ith variable
		double perturbation = 1e-6;
		double originalValue = objectiveFunction(point);
		point[i] += perturbation;
		double perturbedValue = objectiveFunction(point);
		point[i] -= perturbation;

		// Update the gradient
		gradient[i] = (perturbedValue - originalValue) / perturbation;
	}

	// Update the solution using gradient descent
	for (int i = 0; i < point.size(); ++i) {
		point[i] -= gradientStep * gradient[i];
	}
}


void printPopulation(const char* m,vector<vector<double>> p) {
	cout << m << endl;
	int k=-0;
	for (int i = 0; i < N; ++i) {
		cout << k++ << ": ";
		for(int j=0;j<cSize;j++) {
			cout << p[i][j] << "  ";
		}
		cout << " Fitness: " << objectiveFunction(p[i]) << endl;
	}
	cout << endl;
}
