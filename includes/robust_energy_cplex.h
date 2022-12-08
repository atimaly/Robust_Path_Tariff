#include <bits/stdc++.h>
#include <ilconcert/iloenv.h>
#include <random>
#include <lemon/list_graph.h>
#include <lemon/time_measure.h>
#include <lemon/dijkstra.h>
#include <lemon/adaptors.h>
#include <lemon/concepts/path.h>
#include <lemon/concepts/graph.h>

#include <ilcplex/ilocplex.h>

using namespace lemon;
using namespace std;


#ifndef ROBUST_PATH
#define ROBUST_PATH

struct INFEASIBLE{};

using NumVarMatrix = IloArray<IloNumVarArray>;
using NumVar3D = IloArray<NumVarMatrix>;
using NumMatrix = IloArray<IloNumArray>;
using NumMatrix3D = vector<NumMatrix>;

class Graph{
	protected:
		int n_;
		int edge_number_;
		double erdos_edge_possible_;
		ListDigraph g;
		
		Graph(int n, double erdos_p); //You read the input more specially D graph

		Graph(std::istream &is);
};


class Paths : public Graph {
    private:
        IloEnv env;
		int people_n_;
		vector<pair<int,int>> paths_; //Starting and Ending vertices of the people travelling
		ListDigraph::ArcMap<double> arc_cost_q_{g};
		ListDigraph::ArcMap<double> arc_buy_p_{g};

        IloModel polyhedra_q_{env}; IloNumVarArray q_{env};
        IloModel polyhedra_u_{env}; NumVarMatrix u_{env};
		NumMatrix3D set_of_utilities_;

		//vector<vector<int>> defining_polyhedra_q_;
		//vector<vector<int>> defining_polyhedra_u_;
		double leader_max_earn_;
		
		pair<int,int> RandomPath();
		
		void RandomPaths();

		void PerturbationOfq(vector<double> &q_tariff, const double delta, std::ostream &os = std::cerr) const;

		void InitialQValue(vector<double> &q_tariff, std::ostream &os = std::cerr) const;

		double MinimizeLeadersEarning(const vector<double> &q_tariff, const int big_M, std::ostream &os = std::cerr);

		double FindingTariffWithFiniteUtilities(vector<double> &q_tariff, const int big_M, std::ostream &os = std::cerr);

	public:

		Paths(const int people, const int n, const double erdos_p);

		Paths(std::istream &is); //You read the input, more specifically D graph, paths of the people, arc_cost_buy_p, Q polyhedra, U polyhedra

        ~Paths();

		void GenerateProblem(const int seed);

		void FindingOptimalCost(std::ostream &os = std::cerr);

		void PrintData(std::ostream &os = std::cerr) const;
        
		void PrintDataRaw(std::ostream &os) const;

		void SaveGenerated(std::ostream &os);
    
    
};

#endif //ROBUST_PATH


#ifndef UTILITY_TOOLS
#define UTILITY_TOOLS
using ll = long long int;

const long long int INF = std::numeric_limits<long long int>::max();
const bool DEBUG = true;

#define all(x) begin(x), end(x)
#define FOR(i,n) for(int i = 0; i < (n); ++i)
#define FORO(i,n) for(int i = 1; i < (n); ++i)

template <typename C>
void Print_vector(const C &Original, std::ostream &os = std::cerr) {
	for(const auto &v : Original) {
	    os << v << " ";
	}
	os << endl;
}

template <typename C>
void Print_Matrix(const C &M, std::ostream &os = std::cerr) {
	for(auto &v : M) {
		os << v.size() << " ";
		for(auto &u : v) {
			os << u; os << " ";            
			}
	    os << endl;
	}
}

template<class T, class C>
void Print_pair(const pair<T,C> &M, std::ostream &os = std::cerr) {
    os << "(" << M.first << " , " << M.second << " ) ";
}

template <class C>
void Print_vector_pairs(const C &Original, std::ostream &os = std::cerr) {
	for(const auto &v : Original) {
	    Print_pair(v, os);
	}
	os << endl;
}

template <class C>
void Print_vector_pairs_raw(const C &Original, std::ostream &os = std::cerr) {
	for(const auto &v : Original) {
	    os << v.first << " " << v.second << " ";
	}
	os << endl;
}
#endif //UTILITY_TOOLS
