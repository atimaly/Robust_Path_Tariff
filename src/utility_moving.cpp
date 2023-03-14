#include <bits/stdc++.h>
#include <cstdio>
#include <ilconcert/iloenv.h>
#include <ilconcert/iloexpression.h>
#include <limits>
#include <random>
#include <lemon/list_graph.h>
#include <lemon/time_measure.h>
#include <lemon/dijkstra.h>
#include <lemon/adaptors.h>
#include <lemon/concepts/path.h>
#include <lemon/concepts/graph.h>

#include <ilcplex/ilocplex.h>


#include "robust_energy_cplex.h"

using namespace lemon;
using namespace std;


void Paths::MoveInUtilitySpace(const vector<vector<double>> &x_flow, const vector<double> &alpha_plus, const vector<double> &alpha_negative, 
	const vector<map<int,double>> &beta, int index_of_utility, const vector<vector<double>> &utility) {
	//Moving in the utility space while keeping the x_flow the optimal answer
	IloArray<IloNumArray> new_utility(env, people_n_); //vector<double>(edge_number_, 0)
		

	set_of_utilities_.push_back(new_utility);
}

bool Paths::SubstantiallyDifferentyUtility(double delta , int index_of_utility, const vector<vector<double>> &utility) {
	//Check wether the utility at the of index_of_utility in set_of_utilities_ is substantially different than utility
	FOR(i,people_n_) {
		FOR(e,edge_number_) {
			if(edge_number_*delta <= std::abs(utility[i][e] - set_of_utilities_[index_of_utility][i][e])) {
				return true;
			}
		}
	}	
	
	return false;
}