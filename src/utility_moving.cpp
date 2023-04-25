#include <bits/stdc++.h>
#include <cstdio>
#include <ilconcert/iloenv.h>
#include <ilconcert/iloexpression.h>
#include <limits>
#include <ostream>
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

double Paths::MoveInUtilitySpace(const vector<vector<double>> &x_flow, int index_of_utility, NumMatrix &u_sol, std::ostream &os) {
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
				nu[j][e] = static_cast<float>(set_of_utilities_[static_cast<int>(set_of_utilities_.size())-1][j][e]) - static_cast<float>(set_of_utilities_[index_of_utility][j][e]);
			}
		}

	//u_0 + \nu * delta \in U
		FOR(j,people_n_) {
			FOR(e,edge_number_) {
				model.add(u_sub[j][e] == Delta*nu[j][e] + set_of_utilities_[index_of_utility][j][e]);
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
	switch (cplex.getStatus())
	{

		case IloAlgorithm::Optimal:
			obj_value = cplex.getObjValue();
			os << "THE OPT SOL is " << obj_value << endl;
			FOR(i,people_n_) {
				u_sol[i] = IloNumArray(env, edge_number_);
				FOR(e,edge_number_) {
					try{
						u_sol[i][e] = cplex.getValue(u_sub[i][e]);
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
			os << "THE PROBLEM IS UNBOUNDED:\n";return std::numeric_limits<double>::min();

			break;
		case IloAlgorithm::Infeasible:
			os << "THE PROBLEM IS Infeasible:\n";return std::numeric_limits<double>::min();

			break;
		
		case IloAlgorithm::Error:
			os << "THE PROBLEM HAS Error:\n";
			break;

		case IloAlgorithm::Feasible:
			os << "THE PROBLEM is FEASIBLE:\n";return std::numeric_limits<double>::min();

			break;
		
		case IloAlgorithm::InfeasibleOrUnbounded:
			os << "THE PROBLEM is InFEASIBLE or Unbounded:\n";
			return std::numeric_limits<double>::min();
			break;
		
		case IloAlgorithm::Unknown:
			os << "The Problem is unknown\n";
			break;
		default:
			os << "Something has happened.\n";
			break;
	}

	beta.end();
	alpha_negative.end();
	alpha_plus.end();
	u_sub.end();
	Delta.end();
	model.end();
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
			
			if(delta <= std::abs(set_of_utilities_[last_util][i][e] - set_of_utilities_[index_of_utility][i][e])) {
				return true;
			}
		}
	}	
	/*
	if(all_diff <= delta) {
		return false;
	}*/

	return false;
}

void Paths::UtilityMovingIfDifferent(vector<vector<double>> &x_flow, std::ostream &os) {
	
	vector<double> movement_delta(static_cast<int>(set_of_utilities_.size()));
	NumMatrix3D u_vals(static_cast<int>(set_of_utilities_.size())-1);
	FOR(si,static_cast<int>(set_of_utilities_.size())-1) {
		u_vals[si] = NumMatrix(env, people_n_); //Save the utility value for later use
		movement_delta[si] = MoveInUtilitySpace(x_flow, si, u_vals[si]);
		os << "The current delta for: " << si << " is " << movement_delta[si] << endl;
	}
	auto maxi_delta = std::max_element(movement_delta.begin(), movement_delta.end());
	int indi_delta = std::distance(movement_delta.begin(), maxi_delta);
	if(*maxi_delta == std::numeric_limits<double>::min()) {
		set_of_utilities_.pop_back();
		return;
	}
	/*
	NumMatrix nu_best(env, people_n_);
	FOR(j,people_n_) {
		nu_best[j] = IloNumArray(env, edge_number_);
		FOR(e,edge_number_) {
			nu_best[j][e] = + set_of_utilities_[static_cast<int>(set_of_utilities_.size())-1][j][e] - set_of_utilities_[indi_delta][j][e];
		}
	}

	NumMatrix new_util(env, people_n_);
	FOR(j,people_n_) {
		new_util[j] = IloNumArray(env, edge_number_);
		FOR(e,edge_number_) {
			new_util[j][e] = set_of_utilities_[static_cast<int>(set_of_utilities_.size())-1][j][e] + (*maxi_delta)*nu_best[j][e];
		}
	}*/

	set_of_utilities_.pop_back();
	set_of_utilities_.push_back(u_vals[indi_delta]);

}

bool HalvingUtilityMoving(vector<int> &index_to_be_vertices) {
	//This program checks if the index_of_utility-th set_of_utilities_ and the last of set_of_utilities_ half is in U
	

}

