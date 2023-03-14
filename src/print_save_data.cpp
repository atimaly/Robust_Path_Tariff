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


void Paths::PrintData(std::ostream &os) const {
	os << "\n\n-----------------PRINTDATA----------------" << endl;
    os << "The number of vertices: " << n_ << endl;
	os << "The number of arcs: " << edge_number_ << endl;
	os << "The Arcs: " << endl;
	for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
		os << "	(" << g.id(g.source(e)) << ", " << g.id(g.target(e)) << ")" << " the vendor is buying it originally at price: " << arc_buy_p_[e] << endl;
	}
	os << "Number of people travelling: " << people_n_ << endl;
	os << "Their starting and end points: " << endl;
	Print_vector_pairs(paths_, os);
	
	os << "The defining polyhedra for q:" << endl;
	os << polyhedra_q_ << endl;

	os << "The defining polyhedra for u:" << endl;
	os << polyhedra_u_ << endl;

	os << "Current set of utilities: " << endl;
	Print_vector(set_of_utilities_, os);

	os << "Leader's current maximum profit: " << leader_max_earn_ << endl;

	os << "-----------------PRINTDATA----------------" << endl;
}

void Paths::PrintDataRaw(std::ostream &os) const {
	os << n_ << endl;
	FOR(i,n_) {
		int out = 0;
		for (ListDigraph::OutArcIt a(g, g.nodeFromId(i)); a != INVALID; ++a) {
			++out;
		}
		os << out << " ";
		for (ListDigraph::OutArcIt a(g, g.nodeFromId(i)); a != INVALID; ++a) {
			os << g.id(g.target(a)) << " ";
		}
		os << endl;
	}
	for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
		os << arc_buy_p_[e] << " ";
	}
	os << people_n_ << endl;
	Print_vector_pairs_raw(paths_, os);
	/*Print_Matrix(defining_polyhedra_q_, os);
	Print_Matrix(defining_polyhedra_u_, os);*/
}
