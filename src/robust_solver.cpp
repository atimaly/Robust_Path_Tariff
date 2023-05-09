#include "robust_energy_cplex.h"
#include <algorithm>
#include <cstdio>
#include <limits>
#include <random>
#include <set>
#include <lemon/list_graph.h>
#include <lemon/time_measure.h>
#include <lemon/dijkstra.h>
#include <lemon/adaptors.h>
#include <lemon/concepts/path.h>
#include <lemon/concepts/graph.h>
 //#include <lemon/johnson.h>

#include "random_unit_vec.cpp"


using namespace lemon;
using namespace std;

void Paths::PerturbationOfq(vector<double> &q_tariff, const double delta, std::ostream &os) const{
	#if _DEBUG
	os << "\n\n-------PerturbationOfq BEGIN-------" << endl;
	#endif

	std::vector<pair<double, int>> q_indexed;
	int ind{0};
	std::transform(all(q_tariff), std::back_inserter(q_indexed),
		[&ind, this](double &a)
		{auto p = std::make_pair(a-arc_buy_p_[g.arcFromId(ind)], ind);
		 ++ind;
		 return p;
		});

	#if _DEBUG
	os << "\n q_tariff-p indexed\n";
	Print_vector_pairs(q_indexed);
	#endif

	std::sort(all(q_indexed));

	#if _DEBUG
	os << "\n Sorted q_indexed-p\n";
	Print_vector_pairs(q_indexed);
	#endif

	int first_non_negative = std::distance(q_indexed.begin(), std::find_if(all(q_indexed), [](auto &a){return a.first > 0;}));
	//assert( != end());

	#if _DEBUG
	os << "\nFirst index T start: " << first_non_negative << endl;
	#endif

	FOR(i,static_cast<int>(q_indexed.size())) {
		q_indexed[i].first = q_indexed[i].first + arc_buy_p_[g.arcFromId(i)] + (i-first_non_negative)*delta;
	}
	/*
	std::transform(all(q_indexed), q_indexed.begin(), 
		[this, first_non_negative, delta](auto &a)
		{return std::make_pair(a.first + arc_buy_p_[g.arcFromId(a.second)] + (a.second-first_non_negative)*delta, a.second); });*/
	
	#if _DEBUG
	os << "\n Applied perturbation\n";
	Print_vector_pairs(q_indexed);
	#endif

	std::sort(all(q_indexed), [](auto &a, auto &b){return a.second < b.second; });

	#if _DEBUG
	os << "\nSort according to original order\n";
	Print_vector_pairs(q_indexed);
	#endif

	std::transform(all(q_indexed), q_tariff.begin(), [this](auto &a){return a.first;});

	#if _DEBUG
    	os << "\n\n-------PerturbationOfq END-------" << endl;
    	#endif
}

void Paths::InitialQValue(vector<double> &q_tariff, std::ostream &os) const{
	//Getting a q tariff that satisfies the constraints

	#if _DEBUG
    cerr << "\n\n-------InitialQValue BEGIN-------" << endl;
    #endif

	IloModel poly_q(env); poly_q.add(polyhedra_q_); 

	IloExpr expr(env); FOR(i,q_.getSize()) expr += q_[i];
	IloObjective obj(env, expr, IloObjective::Maximize);
	poly_q.add(obj);
	
	IloCplex cplex_q(poly_q);

	#if _DEBUG
    //cplex_q.exportModel("PROBLEM_InitialQValue.lp");
    #endif

	cplex_q.solve();
	if(cplex_q.getStatus() == IloAlgorithm::Optimal) {
        IloNumArray qr(env); cplex_q.getValues(qr, q_);
		#if _DEBUG
		os << "The solution is for a sample q value:\n" << qr << endl;
		#endif
		FOR(i,qr.getSize()) {
			//arc_cost_q_[g.arcFromId(i)] = qr[i];
			q_tariff.push_back(qr[i]);
		}
		qr.end();
    }
	else{os << "CAN NOT FIND SOLUTION FOR SettingQValue():" << endl;}

	expr.end();
	obj.end();
	poly_q.end();
	cplex_q.end();

	#if _DEBUG
    cerr << "\n\n-------InitialQValue END-------" << endl;
    #endif
}

double Paths::MinimizeLeadersEarning(const vector<double> &q_tariff, const int big_M, vector<double> &alpha_pl_return, vector<double> &alpha_neg_return,
		vector<map<int,double>> &beta_return, vector<vector<double>> &x_flow_return, std::ostream &os) {
	//Given the tariff's on the roads, it gives the worst case for the leader.
	//Rerturn's the leader's minimal earning
	#if _DEBUG
    cerr << "\n\n-------MinimizeLeadersEarning BEGIN-------" << endl;
    #endif
	char varname[20]; //Nameing the variables of cplex
	int bound_flow_binary{2}; //for inequalities with binary constraints bounding the flow type expressions

	IloModel model(env);
	// u \in U
		model.add(polyhedra_u_); //Hindering utilities
	#if _DEBUG_EXTRA
    cerr << "-------IMPORTED POLYHEDRA U-------" << endl;
    #endif

	NumVarMatrix x(env, people_n_); //Flow
	FOR(i,people_n_) {
		x[i] = IloNumVarArray(env, edge_number_);
		FOR(j,edge_number_) {
			sprintf(varname, "x_%d_%d", i, j);
			x[i][j] = IloNumVar(env, 0., 1., ILOFLOAT, varname);
		}
	}

	#if _DEBUG_EXTRA
    cerr << "-------DEFINING X-------" << endl;
    #endif

	IloNumVarArray alpha_plus(env, people_n_);
	IloNumVarArray alpha_negative(env, people_n_);
	IloNumVarArray helper_alpha(env, people_n_);
	IloNumVarArray helper2_alpha(env, people_n_);
	FOR(i,people_n_) {
		sprintf(varname, "ap_%d", i);
		alpha_plus[i] = IloNumVar(env, 0., +INFINITY, ILOFLOAT, varname);
		sprintf(varname, "an_%d", i);
		alpha_negative[i] = IloNumVar(env, 0., +INFINITY, ILOFLOAT, varname);

		//a_plus >= 0 and x_j(\rho(t_j)) >= 1
			IloExpr expr(env);
			//x_j(\rho(t_j)) >= 1
			for (ListDigraph::InArcIt a(g, g.nodeFromId(paths_[i].second)); a != INVALID; ++a) {
				expr += x[i][g.id(a)];
			}
			model.add(expr >= 1);

			sprintf(varname, "hap_1_%d", i);
			helper_alpha[i] =IloNumVar(env, 0, 1, ILOINT, varname);
			model.add(alpha_plus[i] <= helper_alpha[i] *big_M);
			model.add(expr-1 <= (1- helper_alpha[i])*bound_flow_binary);
			expr.end();
			



		//a_minus >= 0 and -x_j(\delta(t_j))+1 >= 0
			IloExpr expr_2(env);
			//-x_j(\delta(t_j))+1 >= 0
			for (ListDigraph::OutArcIt a(g, g.nodeFromId(paths_[i].first)); a != INVALID; ++a) {
				expr_2 -= x[i][g.id(a)];
			}
			model.add(expr_2+1 >= 0);
			
			sprintf(varname, "han_2_%d", i);
			helper2_alpha[i] = IloNumVar(env, 0, 1, ILOINT, varname);
			model.add(alpha_negative[i] <= helper2_alpha[i] *big_M);
			model.add(expr_2+1 <= (1- helper2_alpha[i])*bound_flow_binary);
			expr_2.end();
			

	}

	#if _DEBUG_EXTRA
    cerr << "-------ALPHAS DEFINED-------" << endl;
    #endif
	
	vector<map<int,IloNumVar>> beta(people_n_);
	vector<map<int,IloNumVar>> beta_helpers(people_n_);

	FOR(j,people_n_) {
		FOR(v,n_) {
			ListDigraph::Node curr_v = g.nodeFromId(v);
			if(paths_[j].first != v && paths_[j].second != v) {
				sprintf(varname, "b_%d_%d", j, v);
				#if _DEBUG
				//cerr << "Created variable for beta with: " << j << ", " << v << endl;
				#endif
				beta[j][v] = IloNumVar(env, 0., +IloInfinity, ILOFLOAT, varname);
				//beta[i].insert({j,IloNumVar(env, 0., +IloInfinity, ILOFLOAT)});
				IloExpr expr(env);
				int is_there{0};
				for(ListDigraph::InArcIt a(g, curr_v); a != INVALID; ++a) {
					expr += x[j][g.id(a)];
					++is_there;
				}
				for (ListDigraph::OutArcIt a(g, curr_v); a != INVALID; ++a) {
					expr -= x[j][g.id(a)];
					++is_there;
				}
				if(is_there) {
					model.add(expr >= 0);

					sprintf(varname, "hb_%d%d", j, v);
					beta_helpers[j][v] =  IloNumVar(env, 0, 1, ILOINT, varname);
					model.add(beta[j][v] <= beta_helpers[j][v]*big_M);
					model.add(expr <= (1-beta_helpers[j][v])*bound_flow_binary);
					//helper.end();
				}
				
				expr.end();
			}
		}
	}
	

	#if _DEBUG_EXTRA
    cerr << "-------BETA DEFINED-------" << endl;
    #endif
	
    	vector<vector<IloNumVar>> duality_helpers(people_n_);
	FOR(j,people_n_) {
		duality_helpers[j].resize(edge_number_);
		FOR(e,edge_number_) {
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
			expr -= q_tariff[e] + u_[j][e];
			model.add(expr <= 0);

			sprintf(varname, "hdu_%d_%d", j, e);
			duality_helpers[j][e] = IloNumVar(env, 0, 1, ILOINT, varname);
			model.add(x[j][e] <= duality_helpers[j][e]*bound_flow_binary);
			model.add(-expr <= (1-duality_helpers[j][e])*big_M);
			expr.end();
		}
	}

	#if _DEBUG_EXTRA
    cerr << "-------DUAL OPTIMALITY DEFINED-------" << endl;
    #endif

	
	IloExpr expr_obj(env);
	FOR(i,people_n_) {
		FOR(e,edge_number_) {
			expr_obj += (q_tariff[e]-arc_buy_p_[g.arcFromId(e)])*x[i][e];
		}
	}
	IloObjective obj(env, expr_obj, IloObjective::Minimize);
	model.add(obj);

	#if _DEBUG_EXTRA
    cerr << "-------OBJECTIVE DEFINED-------" << endl;
    #endif


	IloCplex cplex(model);
	#if _DEBUG
    cplex.exportModel("PROBLEM_MinimizeLeadersEarning.lp");
    #endif

	cplex.solve();
	double obj_value{0};
	switch (cplex.getStatus())
	{

		case IloAlgorithm::Optimal:
			try{
				NumMatrix ur(env); ur.setSize(people_n_);
				FOR(i,people_n_) {
					ur[i] = IloNumArray(env, edge_number_);
					cplex.getValues(ur[i], u_[i]);
				}
				if(true) {//TODO meddle with the value of u
					set_of_utilities_.push_back(ur);
				}

				#if _DEBUG
				os << "The solution for u is:\n";
				FOR(i,people_n_) {
					os << "\tFor the " << i << "th people\n\t\t" << ur[i] << endl;
				}
				#endif

				NumMatrix xr(env); xr.setSize(people_n_);
				FOR(i,people_n_) {
					xr[i] = IloNumArray(env, edge_number_);
					cplex.getValues(xr[i], x[i]);
					FOR(e,edge_number_) {
						x_flow_return[i][e] = xr[i][e];
					}
				}

				#if _DEBUG
				os << "The solution for x is:\n";
				FOR(i,people_n_) {
					os << "\tFor the " << i << "th people\n\t\t" << xr[i] << endl;
				}
				#endif

				#if _DEBUG
				IloNumArray alpha_p(env, n_); cplex.getValues(alpha_p, alpha_plus);
				os << "Alpha Plus values are :\n" << alpha_p << "\n";
				IloNumArray alpha_n(env, n_); cplex.getValues(alpha_n, alpha_negative);
				os << "Alpha Negativa values are :\n" << alpha_n << "\n";

				FOR(i,people_n_) alpha_pl_return[i] = alpha_p[i];
				FOR(i,people_n_) alpha_neg_return[i] = alpha_n[i];
				
				
				vector<IloNumArray> betar_all(people_n_);
				FOR(j,people_n_) {
					//os << "BETA: " << j << endl;
					FOR(i,n_) {
						betar_all[j] = IloNumArray(env, n_);
						if(beta[j].find(i) != beta[j].end()) {
							//cerr << "Beta asking: " << i << " vertex for 0. person" << endl;
							IloNumVar refb = beta[j][i];
							try{
								betar_all[j][i] = cplex.getValue(refb);
							}
							catch(IloAlgorithm::NotExtractedException) {
								//cerr << "Didn' t get the vertex " << i << endl;
								betar_all[j][i] = 0;
							}
						}
						else betar_all[j][i] = 0;
					}
				}
				os << "Beta values are :\n";
				FOR(j,people_n_)
					os << "For the person: " <<  j << "\n\t"  << betar_all[j] << "\n";
				
				/*
				FOR(p,people_n_) {
					#if _DEBUG
					os << "p: "<< endl;
					#endif
					FOR(i,n_) {
						if(beta[p].find(i) != beta[0].end()) {
							IloNumVar refb = beta[p][i];
							beta_return[p][i] = cplex.getValue(refb);
							#if _DEBUG
							os << "i : " << i << " value: "<< beta_return[p][i] << endl;
							#endif
						}
					}
				}
				*/
				#endif

				obj_value = cplex.getObjValue();
				#if _DEBUG
				os << "The objective value is : " << obj_value << "\n";
				#endif

			}
			catch(IloAlgorithm::NotExtractedException) {
				cerr << "MinimizeLeadersEarning VARIABLES ARE NOT RELATED TO THE OBJECTIVE THEY HAVE BEEN DELETED BECAUSE OF THE REDUCTION\n"; 
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

	expr_obj.end();
	cplex.end();

	FOR(i, people_n_)
	{
		FOR(j, edge_number_) {
			x[i][j].end();
		}
		x[i].end();
	}
	x.end();

	FOR(i, people_n_) {
		alpha_plus[i].end();
		alpha_negative[i].end();
	}
	alpha_plus.end();
	alpha_negative.end();
	FOR(j, people_n_) {
		for(auto m : beta[j])
			m.second.end();
	}
	beta.end();
	FOR(j, people_n_) {
		for(auto m: beta_helpers[j])
			m.second.end();
	}
	beta_helpers.end();

	FOR(j, people_n_) {
		FOR(e, edge_number_) {
			duality_helpers[j][e].end();
		}
	}

	FOR(i, people_n_) {
		helper_alpha[i].end();
		helper2_alpha[i].end();
	}
	helper_alpha.end();
	helper2_alpha.end();
	model.end();

	#if _DEBUG
    cerr << "\n\n-------MinimizeLeadersEarning END-------" << endl;
    #endif

	return obj_value;
}

double Paths::FindingTariffWithFiniteUtilities(vector<double> &q_tariff, const int big_M, std::ostream &os) {
	//Given a discrete set of utility vectors (set_of_utilities_) it determines the leader's best response

	IloModel model(env);
	char varname[20]; //Nameing the variables of cplex
	int bound_flow_binary{2}; //for inequalities with binary constraints bounding the flow type expressions
	int utility_quantity = static_cast<int>(set_of_utilities_.size());

	// q \in Q
		model.add(polyhedra_q_);

	sprintf(varname, "z");
	IloNumVar Z(env, -IloInfinity, +IloInfinity, ILOFLOAT, varname);

	#if _DEBUG
    cerr << "\n\n-------FindingTariffWithFiniteUtilities BEGIN-------" << endl;
    #endif

	#if _DEBUG_EXTRA
    cerr << "-------IMPORTED POLYHEDRA U-------" << endl;
    #endif

	NumVar3D x(env, utility_quantity); //The flow
	FOR(k,utility_quantity) {
		x[k] = NumVarMatrix(env, people_n_);
		FOR(i,people_n_) {
			x[k][i] = IloNumVarArray(env, edge_number_);
			FOR(j,edge_number_) {
				sprintf(varname, "x_%d_%d", i, j);
				x[k][i][j] = IloNumVar(env, 0., 1., ILOFLOAT, varname);
			}
		}
	}

	#if _DEBUG_EXTRA
    cerr << "-------DEFINING X-------" << endl;
    #endif

	NumVarMatrix alpha_plus(env, utility_quantity);
	NumVarMatrix alpha_negative(env, utility_quantity);
	FOR(k,utility_quantity) {
		alpha_plus[k] = IloNumVarArray(env, people_n_);
		alpha_negative[k] = IloNumVarArray(env, people_n_);

		FOR(i,people_n_) {
			sprintf(varname, "ap_%d_%d", k, i);
			alpha_plus[k][i] = IloNumVar(env, 0., +INFINITY, ILOFLOAT, varname);
			sprintf(varname, "an_%d_%d", k, i);
			alpha_negative[k][i] = IloNumVar(env, 0., +INFINITY, ILOFLOAT, varname);

			//a_plus >= 0 and x_j(\rho(t_j)) >= 1
				IloExpr expr(env);
				//x_j(\rho(t_j)) >= 1
				for (ListDigraph::InArcIt a(g, g.nodeFromId(paths_[i].second)); a != INVALID; ++a) {
					expr += x[k][i][g.id(a)];
				}
				model.add(expr >= 1);

				sprintf(varname, "hap_1_%d_%d", k, i);
				IloNumVar helper(env, 0, 1, ILOINT, varname);
				model.add(alpha_plus[k][i] <= helper*big_M);
				model.add(expr-1 <= (1-helper)*bound_flow_binary);
				//helper.end();
				expr.end();
				



			//a_minus >= 0 and -x_j(\delta(t_j))+1 >= 0
				IloExpr expr_2(env);
				//-x_j(\delta(t_j))+1 >= 0
				for (ListDigraph::OutArcIt a(g, g.nodeFromId(paths_[i].first)); a != INVALID; ++a) {
					expr_2 -= x[k][i][g.id(a)];
				}
				model.add(expr_2+1 >= 0);
				
				sprintf(varname, "han_2_%d_%d", k, i);
				IloNumVar helper_2(env, 0, 1, ILOINT, varname);
				model.add(alpha_negative[k][i] <= helper_2*big_M);
				model.add(expr_2+1 <= (1-helper_2)*bound_flow_binary);
				expr_2.end();
				

		}
	}

	#if _DEBUG_EXTRA
    cerr << "-------ALPHAS DEFINED-------" << endl;
    #endif
	
	vector<vector<map<int,IloNumVar>>> beta(utility_quantity, vector<map<int,IloNumVar>>(people_n_));
	FOR(k,utility_quantity) {
		FOR(j,people_n_) {
			FOR(v,n_) {
				ListDigraph::Node curr_v = g.nodeFromId(v);
				if(paths_[j].first != v && paths_[j].second != v) {
					sprintf(varname, "b_%d_%d_%d", k, j, v);
					beta[k][j][v] = IloNumVar(env, 0., +IloInfinity, ILOFLOAT, varname);
					
					IloExpr expr(env);
					int is_there{0};
					for(ListDigraph::InArcIt a(g, curr_v); a != INVALID; ++a) {
						expr += x[k][j][g.id(a)];
						++is_there;
					}
					for (ListDigraph::OutArcIt a(g, curr_v); a != INVALID; ++a) {
						expr -= x[k][j][g.id(a)];
						++is_there;
					}
					if(is_there) {
						model.add(expr >= 0);

						sprintf(varname, "hb_%d_%d_%d", k, j, v);
						IloNumVar helper(env, 0, 1, ILOINT, varname);
						model.add(beta[k][j][v] <= helper*big_M);
						model.add(expr <= (1-helper)*bound_flow_binary);
						//helper.end();
					}
					
					expr.end();
				}
			}
		}
	}	

	#if _DEBUG_EXTRA
    cerr << "-------BETA DEFINED-------" << endl;
    #endif
	FOR(k,utility_quantity) {
		FOR(j,people_n_) {
			FOR(e,edge_number_) {
				IloExpr expr(env);
				int source = g.id(g.source(g.arcFromId(e))); int target = g.id(g.target(g.arcFromId(e)));
				if(target == paths_[j].second) expr += alpha_plus[k][j];
				if(source == paths_[j].first) expr -= alpha_negative[k][j];

				if(paths_[j].first != target && paths_[j].second != target) {
					expr += beta[k][j][target];
				}
				if(paths_[j].first != source && paths_[j].second != source) {
					expr -= beta[k][j][source];
				}
				#if _DEBUG_EXTRA
				cerr << "Trying to reach the " << k << "th utility" << endl;
				#endif
				expr -= q_[e] + set_of_utilities_[k][j][e];
				model.add(expr <= 0);

				sprintf(varname, "hdu_%d_%d_%d", k, j, e);
				IloNumVar helper(env, 0, 1, ILOINT, varname);
				model.add(x[k][j][e] <= helper*bound_flow_binary);
				model.add(-expr <= (1-helper)*big_M);
				expr.end();
			}
		}
	}

	#if _DEBUG_EXTRA
    cerr << "-------DUAL OPTIMALITY DEFINED-------" << endl;
    #endif

	FOR(k,utility_quantity) {
		IloExpr expr(env);
		FOR(j,people_n_) {
			FOR(e,edge_number_) {
				//IloExpr temp_expr(env);
				//int source = g.id(g.source(g.arcFromId(e))); int target = g.id(g.target(g.arcFromId(e)));
				//if(target == paths_[j].second) temp_expr += alpha_plus[j];
				//if(source == paths_[j].first) temp_expr -= alpha_negative[j];

				/*if(paths_[j].first != target && paths_[j].second != target) {
					temp_expr += beta[j][target];
				}
				if(paths_[j].first != source && paths_[j].second != source) {
					temp_expr -= beta[j][source];
				}*/
				
				//temp_expr -= set_of_utilities_[i][j][e]*x[j][e];
				//temp_expr -= q_[e] + set_of_utilities_[i][j][e];
				//expr += temp_expr;

				//temp_expr.end();
				expr -= set_of_utilities_[k][j][e]*x[k][j][e];
				expr -= arc_buy_p_[g.arcFromId(e)]*x[k][j][e];
			}
			expr += alpha_plus[k][j]; expr -= alpha_negative[k][j];
		}
		/*
		FOR(e,edge_number_) {
			expr -= arc_buy_p_[g.arcFromId(e)]*x[j][e];
		}
		*/
		model.add(Z <= expr);
		expr.end();
	}

	#if _DEBUG_EXTRA
    cerr << "-------CONSTRAINTS ON Z DEFINED-------" << endl;
    #endif

	IloObjective obj(env, Z, IloObjective::Maximize);
	model.add(obj);

	#if _DEBUG_EXTRA
	cerr << "-------OBJECTIVE DEFINED-------" << endl;
	#endif

	IloCplex cplex(model);

	#if _DEBUG
    cplex.exportModel("PROBLEM_FindingTariffWithFiniteUtilities.lp");
    #endif

	cplex.solve();
	double obj_value{0};
	switch (cplex.getStatus())
	{

		case IloAlgorithm::Optimal:
			try{
			NumMatrix3D xr(utility_quantity); //xr.setSize(people_n_);
			FOR(k,utility_quantity) {
				xr[k] = NumMatrix(env, people_n_);
				FOR(i,people_n_) {
					xr[k][i] = IloNumArray(env, edge_number_);
					cplex.getValues(xr[k][i], x[k][i]);
				}
			}

			#if _DEBUG
			os << "The solution for x is:\n";
			FOR(k,utility_quantity) {
				os << "For the " << k << "th utility\n";
				FOR(i,people_n_) {
					os << "\tFor the " << i << "th people\n";
					os << "\t\t" << xr[k][i] << endl;
				}
			}
			#endif

			IloNumArray qr(env, edge_number_);
			cplex.getValues(qr, q_);
			#if _DEBUG
			os << "The solution for q is:\n";
			os << qr << endl;
			#endif

			FOR(i,qr.getSize()) {
				//arc_cost_q_[g.arcFromId(i)] = qr[i];
				q_tariff[i] = qr[i];
			}
			qr.end();

			#if _DEBUG
			NumMatrix alpha_p(env, utility_quantity); 
			FOR(k,utility_quantity) {
				alpha_p[k] = IloNumArray(env, n_);
				cplex.getValues(alpha_p[k], alpha_plus[k]);
			}
			FOR(k,utility_quantity) os << "Alpha Plus values are :\n" << alpha_p[k] << "\n";

			//IloNumArray alpha_n(env, n_); cplex.getValues(alpha_n, alpha_negative);
			NumMatrix alpha_n(env, utility_quantity); 
			FOR(k,utility_quantity) {
				alpha_n[k] = IloNumArray(env, n_);
				cplex.getValues(alpha_n[k], alpha_negative[k]);
			}
			FOR(k,utility_quantity) os << "Alpha Negativa values are :\n" << alpha_n[k] << "\n";


			/*
			vector<vector<double>> betar(utility_quantity, vector<double>(n_));
			FOR(k,utility_quantity) {
				FOR(i,n_) {
					if(beta[k][0].find(i) != beta[k][0].end()) {
						IloNumVar refb = beta[k][0][i];
						betar[k][i] = cplex.getValue(refb);
					}
					else betar[k][i] = 0;
				}
			}
			FOR(k,utility_quantity) FOR(i,n_) os << "Beta values are :\n" << betar[k][i] << "\n";*/


			#endif
			
			obj_value = cplex.getObjValue();
			#if _DEBUG
			os << "Maximum objective value is: " << cplex.getObjValue();
			#endif
			}
			catch(IloAlgorithm::NotExtractedException) {
				cerr << "FindingTariffWithFiniteUtilities VARIABLES ARE NOT RELATED TO THE OBJECTIVE\n"; 
			}
			break;
	
		case IloAlgorithm::Unbounded:
			os << "THE PROBLEM IS UNBOUNDED:\n";
			throw INFEASIBLE{};
			break;
		case IloAlgorithm::Infeasible:
			os << "THE PROBLEM IS Infeasible:\n";
			break;
		
		case IloAlgorithm::Error:
			os << "THE PROBLEM HAS Error:\n";
			break;
		case IloAlgorithm::InfeasibleOrUnbounded:
			cplex.exportModel("Infeasible_Fixed_Utility_Prob.lp");
			os << "THE PROBLEM is InFEASIBLE or Unbounded:\n";
			break;
		
		case IloAlgorithm::Unknown:
			os << "The Problem is unknown\n";
			break;
		default:
			os << "Something has happened.\n";
			break;

	}


	Z.end();
	//expr_obj.end();
	cplex.end();
	FOR(k,utility_quantity) {
		FOR(i,people_n_) {
			FOR(j,edge_number_) {
				x[k][i][j].end();
			}
			x[k][i].end();
		}
		x[k].end();
	}
	x.end();
	FOR(k,utility_quantity) {
		FOR(i,people_n_) {
			alpha_plus[k][i].end();
			alpha_negative[k][i].end();
		}
		alpha_plus[k].end();
		alpha_negative[k].end();
	}
	alpha_plus.end();
	alpha_negative.end();
	FOR(k,utility_quantity) {
		FOR(j,people_n_) {
			for(auto m: beta[k][j]) {
				m.second.end();
			}
			beta[k][j].end();
		}
		beta[k].end();
	}
	beta.end();
	model.end();

	#if _DEBUG
    cerr << "\n\n-------FindingTariffWithFiniteUtilities END-------" << endl;
    #endif

	return obj_value;
}

void Paths::FindingOptimalCost(std::ostream &os) {
	//Getting an initial q value saving it in arc_cost_q_
		vector<double> q_tariff;
		InitialQValue(q_tariff);
			
	//double leader_max_earn_{-1};
	int iteration{50};
	set<double> all_of_leaders_earnings;
	set<vector<double>> all_of_q_tariffs;
	FOR(i,iteration) {
		#if _DEBUG
		os << "\n\n------------------------------------------------------------------------------------\n";
		os << "---------------------ITERATION: " << i << " --------------------------------------------\n";
		os << "------------------------------------------------------------------------------------\n\n";
		#endif
		//Petrubate q_tariff
			PerturbationOfq(q_tariff, 1e-5);

		all_of_q_tariffs.insert(q_tariff); //To check if it changes

		//Section 3.4. Solver Finding worst case for leader
			int big_M = 300;
			vector<double> alpha_pl_return(people_n_); vector<double> alpha_neg_return(people_n_);
			vector<map<int,double>> beta_return(people_n_); vector<vector<double>> x_flow_return(people_n_, vector<double>(edge_number_)); 
			double current_leader_earning = MinimizeLeadersEarning(q_tariff, big_M, alpha_pl_return, alpha_neg_return, beta_return, x_flow_return); //hyperparam TODO to-tune
		
		//Section 3.5.
			try{
				double leader_earn = FindingTariffWithFiniteUtilities(q_tariff, big_M);
				leader_max_earn_ = std::min(leader_earn, leader_max_earn_);
				all_of_leaders_earnings.insert(floor(leader_earn * 1000.0) / 1000.0);
			}
			catch(INFEASIBLE) {
				os << "The program has become Infeasible ending it here.\n"; 
				break;
			}
			
			//Is the current q_tariff an optimal solution?
			if(leader_max_earn_ <= current_leader_earning + 0.001) { //TO-TUNE
				os << "\n\n-----------------------Found the optimal solution in the iteration: " << i << " and the optimal value is: " << leader_max_earn_  << "    -------------------------" << endl << endl;
				optimal_q_val = q_tariff;
				break;
			}

			int how_different{0};
			FOR(si,static_cast<int>(set_of_utilities_.size())-1) {
				if(SubstantiallyDifferentyUtility(0.1, si)) {
					++how_different;
				}
			}
			#if _DEBUG
			os << "The new utility is " << how_different << " different and there are " << static_cast<int>(set_of_utilities_.size())-1 << " amount of utilities." << endl;
			#endif

			os << "\n\nCurrent set of utilities: " << endl;
			Print_vector(set_of_utilities_, os);
			
			//The new utility is not that different
			if(how_different != static_cast<int>(set_of_utilities_.size())-1) {
				UtilityMovingIfDifferent(x_flow_return, 100);
			}

			os << "\n\nCurrent set of utilities: " << endl;
			Print_vector(set_of_utilities_, os);

			os << "Leader's current maximum profit: " << leader_max_earn_ << endl << endl;
			os << "Every leader's profit we have encountered" << endl;
			Print_vector(all_of_leaders_earnings, os);
			
			/*
			if(static_cast<int>(all_of_leaders_earnings.size()) >= 3) {
				break;
			}*/
			//PrintData();
	}
	os << "-------------------------------------------------ALL OF THE ITERATIONS HAS ENDED-----------------------------------------" << endl << endl << endl;
	//Print out the encountered q_tariffs
	os << "The encountered q_tariffs." << endl;
	for(auto &v : all_of_q_tariffs) {
		Print_vector(v, os);
	}

	//PrintData();
}
