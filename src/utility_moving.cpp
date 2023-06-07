#include "robust_energy_cplex.h"

#include <algorithm>
#include <cstdio>
//#include <ilconcert/iloenv.h>
//#include <ilconcert/iloexpression.h>
#include <limits>
#include <ostream>
#include <random>
#include <lemon/list_graph.h>
#include <lemon/time_measure.h>
#include <lemon/dijkstra.h>
#include <lemon/adaptors.h>
#include <lemon/concepts/path.h>
#include <lemon/concepts/graph.h>


using namespace lemon;

double Paths::MoveInUtilitySpace(const vector<vector<double>> &x_flow, int index_of_utility, vector<vector<double>> &u_sol, std::ostream &os) {
	//Moving in the utility space while keeping the x_flow the optimal answer
	//IloArray<IloNumArray> new_utility(env, people_n_); //vector<double>(edge_number_, 0)
	
	IloModel model(env);
	char varname[20];

	//add Q poly
		//model.add(polyhedra_q_);
	
	//Q poly for the u to from which we move
		NumVarMatrix u_sub(env);
		u_sub.setSize(people_n_); 
		FOR(i,people_n_) {
			u_sub[i] = IloNumVarArray(env, edge_number_);
			FOR(j,edge_number_) {
				sprintf(varname, "u_sub_%d%d", i, j);
				u_sub[i][j] = IloNumVar(env, 0., +IloInfinity, ILOFLOAT, varname);
			}
		}
		int lines = static_cast<int>(defining_polyhedra_u_.size());
		FOR(i, lines) {
			int variab = defining_polyhedra_u_[i][0];
			IloExpr expr(env);
			int person = defining_polyhedra_u_[i][1];
			int indi_vs = 2;
			FOR(_j, variab) {
				#if _DEBUG_EXTRA
				cerr << "i: " << i << " _j: " << _j << endl;
				#endif
				int u = defining_polyhedra_u_[i][indi_vs]; ++indi_vs;
				int v = defining_polyhedra_u_[i][indi_vs]; ++indi_vs;
				double coeff = defining_polyhedra_u_[i][indi_vs]; ++indi_vs;
				expr += coeff*u_sub[person][pair_to_arc[make_pair(u, v)]];
			}
			double maxi = defining_polyhedra_u_[i].back();
			model.add(expr <= maxi);
			expr.end();
		}
	
	//add U poly
		model.add(polyhedra_u_);
	
	//delta
		sprintf(varname, "delta");
		IloNumVar Delta(env, -IloInfinity, +IloInfinity, ILOFLOAT, varname);
		
	//defining \nu
		NumMatrix nu(env, people_n_);//set_of_utilities_[static_cast<int>(set_of_utilities_.size())-1];
		FOR(j,people_n_) {
			nu[j] = IloNumArray(env, edge_number_);
			FOR(e,edge_number_) {
				nu[j][e] = - static_cast<float>(set_of_utilities_[static_cast<int>(set_of_utilities_.size())-1][j][e]) + static_cast<float>(set_of_utilities_[index_of_utility][j][e]);
			}
		}

	//u_0 + \nu * delta \in U
		FOR(j,people_n_) {
			FOR(e,edge_number_) {
				model.add(u_sub[j][e] == Delta*nu[j][e] + static_cast<float>(set_of_utilities_[static_cast<int>(set_of_utilities_.size())-1][j][e])); // set_of_utilities_[index_of_utility][j][e]);
				//model.add(u_[j][e] == Delta*nu[j][e] + set_of_utilities_[index_of_utility][j][e]);
				// u_sub[j][e]);
			}
		}
	
	//Do not need new x_flow for the new utility
	IloExpr val_sol(env);
	FOR(i,people_n_) {
		FOR(e,edge_number_) {
			val_sol += (q_[e] - arc_buy_p_[g.arcFromId(e)])*x_flow[i][e];
		}
	}
	model.add(val_sol <= leader_max_earn_);

	IloNumVarArray alpha_plus(env, people_n_);
	IloNumVarArray alpha_negative(env, people_n_);
	FOR(i,people_n_) {
		sprintf(varname, "ap_%d", i);
		alpha_plus[i] = IloNumVar(env, 0., +INFINITY, ILOFLOAT, varname);

		//a_j^plus is zero or not
		double alp{0};
		for (ListDigraph::InArcIt a(g, g.nodeFromId(paths_[i].second)); a != INVALID; ++a) {
			alp += x_flow[i][g.id(a)];
		}
		if(alp - 1 > 0) {
			model.add(alpha_plus[i] == 0);
		}

		sprintf(varname, "an_%d", i);
		alpha_negative[i] = IloNumVar(env, 0., +INFINITY, ILOFLOAT, varname);
		//a_j^negative is zero or not
		double aln{0};
		for (ListDigraph::OutArcIt a(g, g.nodeFromId(paths_[i].first)); a != INVALID; ++a) {
			aln += x_flow[i][g.id(a)];
		}
		if(1 - aln > 0) {
			model.add(alpha_negative[i] == 0);
		}
	}

	vector<map<int,IloNumVar>> beta(people_n_);
	//vector<map<int,IloNumVar>> beta_helpers(people_n_);
	FOR(j,people_n_) {
		FOR(v,n_) {
			ListDigraph::Node curr_v = g.nodeFromId(v);
			if(paths_[j].first != v && paths_[j].second != v) {
				sprintf(varname, "b_%d_%d", j, v);
				beta[j][v] = IloNumVar(env, 0., +IloInfinity, ILOFLOAT, varname);
				double val_diff{0};
				int is_there{0};
				for(ListDigraph::InArcIt a(g, curr_v); a != INVALID; ++a) {
					val_diff+= x_flow[j][g.id(a)];
					++is_there;
				}
				for (ListDigraph::OutArcIt a(g, curr_v); a != INVALID; ++a) {
					val_diff -= x_flow[j][g.id(a)];
					++is_there;
				}
				if(is_there && val_diff > 0) {
					model.add(beta[j][v] == 0);
				}
			}
		}
	}

	FOR(j,people_n_) {
		FOR(e,edge_number_) {
			if(x_flow[j][e] == 0) {
				IloExpr expr(env);
				int source = g.id(g.source(g.arcFromId(e))); int target = g.id(g.target(g.arcFromId(e)));
				if(target == paths_[j].second) expr += alpha_plus[j];
				if(source == paths_[j].first) expr -= alpha_negative[j];
				
				if(paths_[j].first != target && paths_[j].second != target) {
					expr += beta[j][target];
				}
				if(paths_[j].first != source && paths_[j].second != source) {
					expr -= beta[j][source];
				}
				expr += q_[e] + u_sub[j][e];
				model.add(expr <= 0);
				
			}
		}
	}

	
	IloObjective obj(env, Delta, IloObjective::Maximize);
	model.add(obj);

	IloCplex cplex(model);

	cplex.solve();
	double obj_value{0};
	auto stat_solv = cplex.getStatus();
	switch (stat_solv)
	{

		case IloAlgorithm::Optimal:
			obj_value = cplex.getObjValue();
			#if _DEBUG_EXTRA
			os << "THE OPT SOL is " << obj_value << endl;
			#endif
			FOR(i,people_n_) {
				u_sol[i] = vector<double>(edge_number_, 0);
				FOR(e,edge_number_) {
					try{
						u_sol[i][e] = static_cast<double>(cplex.getValue(u_sub[i][e]));
					}
					catch(IloAlgorithm::NotExtractedException) {
						//cerr << "Didn' t get the vertex " << i << endl;
						u_sol[i][e] = 0;
					}
				}
				//cplex.getValues(u_sol[i], u_sub[i]);
			}
			break;
		case IloAlgorithm::Unbounded:
			os << "THE PROBLEM IS UNBOUNDED:\n";

			break;
		case IloAlgorithm::Infeasible:
			os << "THE PROBLEM IS Infeasible:\n";

			break;
		
		case IloAlgorithm::Error:
			os << "THE PROBLEM HAS Error:\n";
			break;

		case IloAlgorithm::Feasible:
			os << "THE PROBLEM is FEASIBLE:\n";
			break;
		
		case IloAlgorithm::InfeasibleOrUnbounded:
			os << "THE PROBLEM is InFEASIBLE or Unbounded:\n";
			break;
		
		case IloAlgorithm::Unknown:
			os << "The Problem is unknown\n";
			break;
		default:
			os << "Something has happened.\n";
			break;
	}
	

	FOR(j, people_n_) {
		for(auto m : beta[j])
			m.second.end();
	}
	beta.end();
	
	FOR(i, people_n_) {
		alpha_plus[i].end();
		alpha_negative[i].end();
	}
	alpha_plus.end();
	alpha_negative.end();
	Delta.end();
	FOR(j,people_n_) {
		FOR(e,edge_number_) {
			u_sub[j][e].end();
		}
	}
	u_sub.end();
	model.end();
	cplex.end();
	if(stat_solv != IloAlgorithm::Optimal) return std::numeric_limits<double>::min();
	return obj_value;
}


double Paths::MoveInUtilitySpaceRandomDirection(const vector<vector<double>> &x_flow, RandomUnitVecGen &uni_gen, vector<vector<double>> &u_sol, int which_dir, std::ostream &os) {
	//Moving in the utility space while keeping the x_flow the optimal answer
	
	IloModel model(env);
	char varname[20];

		
	//Q poly for the u to from which we move
		NumVarMatrix u_sub(env);
		u_sub.setSize(people_n_); 
		FOR(i,people_n_) {
			u_sub[i] = IloNumVarArray(env, edge_number_);
			FOR(j,edge_number_) {
				sprintf(varname, "u_sub_%d%d", i, j);
				u_sub[i][j] = IloNumVar(env, 0., +IloInfinity, ILOFLOAT, varname);
			}
		}
		int lines = static_cast<int>(defining_polyhedra_u_.size());
		FOR(i, lines) {
			int variab = defining_polyhedra_u_[i][0];
			IloExpr expr(env);
			int person = defining_polyhedra_u_[i][1];
			int indi_vs = 2;
			FOR(_j, variab) {
				#if _DEBUG_EXTRA
				cerr << "i: " << i << " _j: " << _j << endl;
				#endif
				int u = defining_polyhedra_u_[i][indi_vs]; ++indi_vs;
				int v = defining_polyhedra_u_[i][indi_vs]; ++indi_vs;
				double coeff = defining_polyhedra_u_[i][indi_vs]; ++indi_vs;
				expr += coeff*u_sub[person][pair_to_arc[make_pair(u, v)]];
			}
			double maxi = defining_polyhedra_u_[i].back();
			model.add(expr <= maxi);
			expr.end();
		}
	
	//add U poly
		model.add(polyhedra_u_);
	
	//delta
		sprintf(varname, "delta");
		IloNumVar Delta(env, -IloInfinity, +IloInfinity, ILOFLOAT, varname);
		
	//defining \nu
		NumMatrix nu(env, people_n_);//set_of_utilities_[static_cast<int>(set_of_utilities_.size())-1];
				/*
		if(which_dir == 0) random_dir = uni_gen.GenerateDUnitVec(edge_number_);
		if(which_dir == 1) random_dir = uni_gen.GenerateDUnitVecPositive(edge_number_);
		if(which_dir == 2) random_dir = uni_gen.GenerateDUnitVecNegative(edge_number_);*/

		FOR(j,people_n_) {
			nu[j] = IloNumArray(env, edge_number_);
			vector<double> random_dir; // =  uni_gen.GenerateDUnitVec(edge_number_);
			if(which_dir == 0) random_dir = uni_gen.GenerateDUnitVec(edge_number_);
			if(which_dir == 1) random_dir = uni_gen.GenerateDUnitVecPositive(edge_number_);
			if(which_dir == 2) random_dir = uni_gen.GenerateDUnitVecNegative(edge_number_);
			FOR(e,edge_number_) {
				nu[j][e] = static_cast<float>(random_dir[e]);
			}
		}

	//u_0 + \nu * delta \in U
		FOR(j,people_n_) {
			FOR(e,edge_number_) {
				model.add(u_sub[j][e] == Delta*nu[j][e] + static_cast<float>(set_of_utilities_[static_cast<int>(set_of_utilities_.size())-1][j][e]));
				//model.add(u_[j][e] == Delta*nu[j][e] + set_of_utilities_[index_of_utility][j][e]);
				// u_sub[j][e]);
			}
		}
	
	//Do not need new x_flow for the new utility
	IloExpr val_sol(env);
	FOR(i,people_n_) {
		FOR(e,edge_number_) {
			val_sol += (q_[e] - arc_buy_p_[g.arcFromId(e)])*x_flow[i][e];
		}
	}
	model.add(val_sol <= leader_max_earn_);

	IloNumVarArray alpha_plus(env, people_n_);
	IloNumVarArray alpha_negative(env, people_n_);
	FOR(i,people_n_) {
		sprintf(varname, "ap_%d", i);
		alpha_plus[i] = IloNumVar(env, 0., +INFINITY, ILOFLOAT, varname);

		//a_j^plus is zero or not
		double alp{0};
		for (ListDigraph::InArcIt a(g, g.nodeFromId(paths_[i].second)); a != INVALID; ++a) {
			alp += x_flow[i][g.id(a)];
		}
		if(alp - 1 > 0) {
			model.add(alpha_plus[i] == 0);
		}

		sprintf(varname, "an_%d", i);
		alpha_negative[i] = IloNumVar(env, 0., +INFINITY, ILOFLOAT, varname);
		//a_j^negative is zero or not
		double aln{0};
		for (ListDigraph::OutArcIt a(g, g.nodeFromId(paths_[i].first)); a != INVALID; ++a) {
			aln += x_flow[i][g.id(a)];
		}
		if(1 - aln > 0) {
			model.add(alpha_negative[i] == 0);
		}
	}

	vector<map<int,IloNumVar>> beta(people_n_);
	FOR(j,people_n_) {
		FOR(v,n_) {
			ListDigraph::Node curr_v = g.nodeFromId(v);
			if(paths_[j].first != v && paths_[j].second != v) {
				sprintf(varname, "b_%d_%d", j, v);
				beta[j][v] = IloNumVar(env, 0., +IloInfinity, ILOFLOAT, varname);
				double val_diff{0};
				int is_there{0};
				for(ListDigraph::InArcIt a(g, curr_v); a != INVALID; ++a) {
					val_diff+= x_flow[j][g.id(a)];
					++is_there;
				}
				for (ListDigraph::OutArcIt a(g, curr_v); a != INVALID; ++a) {
					val_diff -= x_flow[j][g.id(a)];
					++is_there;
				}
				if(is_there && val_diff > 0) {
					model.add(beta[j][v] == 0);
				}
			}
		}
	}

	FOR(j,people_n_) {
		FOR(e,edge_number_) {
			if(x_flow[j][e] == 0) {
				IloExpr expr(env);
				int source = g.id(g.source(g.arcFromId(e))); int target = g.id(g.target(g.arcFromId(e)));
				if(target == paths_[j].second) expr += alpha_plus[j];
				if(source == paths_[j].first) expr -= alpha_negative[j];
				
				if(paths_[j].first != target && paths_[j].second != target) {
					expr += beta[j][target];
				}
				if(paths_[j].first != source && paths_[j].second != source) {
					expr -= beta[j][source];
				}
				expr += q_[e] + u_sub[j][e];
				model.add(expr <= 0);
				
			}
		}
	}

	
	IloObjective obj(env, Delta, IloObjective::Maximize);
	model.add(obj);

	IloCplex cplex(model);

	cplex.solve();
	double obj_value{0};
	auto stat_solv = cplex.getStatus();
	switch (stat_solv)
	{

		case IloAlgorithm::Optimal:
			obj_value = cplex.getObjValue();
			#if _DEBUG_EXTRA
			os << "THE OPT SOL is " << obj_value << endl;
			#endif
			FOR(i,people_n_) {
				u_sol[i] = vector<double>(edge_number_, 0);
				FOR(e,edge_number_) {
					try{
						u_sol[i][e] = static_cast<double>(cplex.getValue(u_sub[i][e]));
					}
					catch(IloAlgorithm::NotExtractedException) {
						//cerr << "Didn' t get the vertex " << i << endl;
						u_sol[i][e] = 0;
					}
				}
				//cplex.getValues(u_sol[i], u_sub[i]);
			}
			break;
		case IloAlgorithm::Unbounded:
			os << "THE PROBLEM IS UNBOUNDED:\n";
			break;
		case IloAlgorithm::Infeasible:
			os << "THE PROBLEM IS Infeasible:\n";
			break;
		
		case IloAlgorithm::Error:
			os << "THE PROBLEM HAS Error:\n";
			break;

		case IloAlgorithm::Feasible:
			os << "THE PROBLEM is FEASIBLE:\n";
			break;
		
		case IloAlgorithm::InfeasibleOrUnbounded:
			os << "THE PROBLEM is InFEASIBLE or Unbounded:\n";
			break;
		
		case IloAlgorithm::Unknown:
			os << "The Problem is unknown\n";
			break;
		default:
			os << "Something has happened.\n";
			break;
	}
	

	FOR(j, people_n_) {
		for(auto m : beta[j])
			m.second.end();
	}
	beta.end();
	
	FOR(i, people_n_) {
		alpha_plus[i].end();
		alpha_negative[i].end();
	}
	alpha_plus.end();
	alpha_negative.end();
	Delta.end();
	FOR(j,people_n_) {
		FOR(e,edge_number_) {
			u_sub[j][e].end();
		}
	}
	u_sub.end();
	model.end();
	cplex.end();
	if(stat_solv != IloAlgorithm::Optimal) return std::numeric_limits<double>::min();
	
	return obj_value;
}

bool Paths::SubstantiallyDifferentyUtility(double delta , int index_of_utility) {
	//Check wether the utility at the of index_of_utility in set_of_utilities_ is substantially different than utility
	//The utility we check is the last one in the set_of_utilities_  
	int last_util = static_cast<int>(set_of_utilities_.size())-1;
	double all_diff{0};
	FOR(i,people_n_) {
		FOR(e,edge_number_) {
			all_diff += std::abs(set_of_utilities_[last_util][i][e] - set_of_utilities_[index_of_utility][i][e]);
			/*
			if(delta <= std::abs(set_of_utilities_[last_util][i][e] - set_of_utilities_[index_of_utility][i][e])) {
				return true;
			}*/
		}
	}	
	
	if(all_diff <= delta) {
		return false;
	}

	return true;
}

bool Paths::UtilityMovingIfDifferent(vector<vector<double>> &x_flow, int random_dir_tries, std::ostream &os) {
	
	int current_utility_number = static_cast<int>(set_of_utilities_.size());
	//random_dir_tries *= 3; //Because we try to move in positive, negative normal direction
	vector<double> movement_delta(current_utility_number+random_dir_tries); //Store the optimal deltas
	vector<vector<vector<double>>> u_vals(current_utility_number+random_dir_tries-1, vector<vector<double>>(people_n_)); //Store the created new utilities

	//Try Moving in given direction
	FOR(si,current_utility_number-1) {
		//u_vals[si] = vector<vector<double>>(people_n_); //NumMatrix(env, people_n_); //Save the utility value for later use
		#if _DEBUG_EXTRA
		os << "Utility move calc: " << si << endl;
		#endif
		movement_delta[si] = MoveInUtilitySpace(x_flow, si, u_vals[si]);
		#if _DEBUG_EXTRA
		os << "The current delta for: " << si << " is " << movement_delta[si] << endl;
		os << "The new utility's value.";
		Print_MatrixClear(u_vals[si]);
		#endif
	}

	//Try Moving in Random Unit vector Direction 
	//Try moving in Random positive unit vector direction
	//Try moving in random negative unit vector direction
	RandomUnitVecGen uni_gen;
	current_utility_number += random_dir_tries;
	#if _DEBUG_EXTRA
	os << "current_utility_number: " << current_utility_number << endl;
	#endif
	std::discrete_distribution<> d({1, 1, 1});
	for(int si = current_utility_number-random_dir_tries-1; si < current_utility_number-1; ++si) {
		movement_delta[si] = MoveInUtilitySpaceRandomDirection(x_flow, uni_gen, u_vals[si], static_cast<int>(d(gen))); //, static_cast<int>(d(gen)));
		#if _DEBUG_EXTRA
		os << "The current delta for: " << si << " is " << movement_delta[si] << endl;
		os << "The new utility's value.";
		Print_MatrixClear(u_vals[si]);
		#endif

	}
	

	//Can we move in any direction?
	auto maxi_delta = std::max_element(movement_delta.begin(), movement_delta.end());
	int indi_delta = std::distance(movement_delta.begin(), maxi_delta);
	if(*maxi_delta == std::numeric_limits<double>::min()) {
		//set_of_utilities_.pop_back();
		return false;
	}
	#if _DEBUG_EXTRA
	os << "Euclidean distance calculating." << endl;
	#endif

	//Euclid dist between the new utilities
	vector<double> utility_euclid_distances(current_utility_number-1, 0.);
	vector<double> utility_euclid_min_distances(current_utility_number-1, std::numeric_limits<double>::max());
	FOR(i,current_utility_number-1) {
		if(movement_delta[i] != std::numeric_limits<double>::min()) //It can be solved
		FOR(t,current_utility_number-random_dir_tries-1) {
			#if _DEBUG_EXTRA
			os << "Want to reach " << i << ", " << t << endl;
			#endif
			double temp_dist = DifferenceBetweenUtilityies(u_vals[i], t);
			utility_euclid_distances[i] += temp_dist;
			utility_euclid_min_distances[i] = std::min(utility_euclid_min_distances[i], temp_dist);
		}
	}

	#if _DEBUG_EXTRA
	os << "For the created utilities the current euclidean distances from set_of_utilities_ utilities\n";
	Print_vector(utility_euclid_distances);
	#endif
	
	auto most_different_util = std::max_element(utility_euclid_distances.begin(), utility_euclid_distances.end());
	int indi_most_diff_util = std::distance(utility_euclid_distances.begin(), most_different_util);
	//auto min_different_util = std::min_element(utility_euclid_distances.begin(), utility_euclid_distances.end());
	
	if(utility_euclid_min_distances[indi_most_diff_util] <= 0.01) return false; //If the new utility is too close to a previous one then find a new one
	set_of_utilities_.pop_back();
	NumMatrix best_new_utility(env, people_n_);
	FOR(j,people_n_) {
		best_new_utility[j] = IloNumArray(env, edge_number_);
		FOR(e,edge_number_) {
			best_new_utility[j][e] = u_vals[indi_most_diff_util][j][e]; //This is where we choose which utility we use
		}
	}
	set_of_utilities_.push_back(best_new_utility);
	return true;

}

double Paths::DifferenceBetweenUtilityies(vector<vector<double>> &u_vals, int index_of_utility, std::ostream &os) {
	//This function calculates the euclidian distance between the u_vals utility and index_of_utility-th utility in set_of_utilities_
	#if _DEBUG_EXTRA
	os << "DifferenceBetweenUtilityies\n";
	#endif
	double euclid_diff{0};
	FOR(j,people_n_) {
		FOR(e,edge_number_) {
			#if _DEBUG_EXTRA
			os << "Distance at point " << j << ", " << e << endl;
			os << "Value of u_vals[j][e] " << u_vals[j][e] << endl;
			#endif

			euclid_diff += std::abs(u_vals[j][e] - set_of_utilities_[index_of_utility][j][e]);
		}
	}

	return euclid_diff;
}


bool HalvingUtilityMoving(vector<int> &index_to_be_vertices) {
	//This program checks if the index_of_utility-th set_of_utilities_ and the last of set_of_utilities_ half is in U
	
	return 0;
}

