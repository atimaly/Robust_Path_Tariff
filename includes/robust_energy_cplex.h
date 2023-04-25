#include <bits/stdc++.h>
#include <ilconcert/iloenv.h>
#include <lemon/adaptors.h>
#include <lemon/concepts/graph.h>
#include <lemon/concepts/path.h>
#include <lemon/dijkstra.h>
#include <lemon/list_graph.h>
#include <lemon/time_measure.h>
#include <limits>
#include <random>

#include <ilcplex/ilocplex.h>
#include <vector>

using namespace lemon;
using namespace std;

#ifndef ROBUST_PATH
#define ROBUST_PATH

struct INFEASIBLE {};

using NumVarMatrix = IloArray<IloNumVarArray>;
using NumVar3D = IloArray<NumVarMatrix>;
using NumMatrix = IloArray<IloNumArray>;
using NumMatrix3D = vector<NumMatrix>;

class Graph {
       protected:
        //std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen{std::random_device{}()};  //Standard mersenne_twister_engine seeded with rd()
        int n_;
        int edge_number_;
        double erdos_edge_possible_;
        ListDigraph g;
        map<pair<int, int>, int> pair_to_arc;  //Given a pair in our original indexing it gives back the ListDigraph's matching arc's id
        //map<int, pair<int,int>> arc_to_pair; //Given the arc's id
	vector<vector<int>> adjacency_list_generated;

        Graph(int n, double erdos_p);  //You read the input more specially D graph

        Graph(std::istream& is);
};

struct PolyCreator {
        //The meta data to create a Subset Polyhedra (every constraint is just a subset of variables summed and bounded by an upper value)
        //Every variable is non negative
        int row_numb_;
        int col_numb_;
        double prob_in_subset_;
        double max_upper_bound_var_;     //The upper bound for the variables it is a uniform distribution between 0 and max_upper_bound_var_
        double max_upper_bound_subset_;  //The upper bound is a normal distribution around max_upper_bound_ with standard deviation 1
};

class Paths : public Graph {
       private:
        IloEnv env;
        int people_n_;
        vector<pair<int, int>> paths_;  //Starting and Ending vertices of the people travelling
        ListDigraph::ArcMap<double> arc_cost_q_{g};
        ListDigraph::ArcMap<double> arc_buy_p_{g};

        IloModel polyhedra_q_{env};
        IloNumVarArray q_{env};
        IloModel polyhedra_u_{env};
        NumVarMatrix u_{env};
        NumMatrix3D set_of_utilities_;

        vector<vector<double>> defining_polyhedra_q_; //After random problem generation it is filled with the necessary information to write out
        vector<vector<double>> defining_polyhedra_u_; //After random problem generation it is filled with the necessary information to write out

        double leader_max_earn_;

        void SubsetPolyhedra();

        pair<int, int> RandomPath(); //Create Random Start, Destination for a person 

        void RandomPaths();

        void PerturbationOfq(vector<double>& q_tariff, const double delta, std::ostream& os = std::cerr) const;

        void InitialQValue(vector<double>& q_tariff, std::ostream& os = std::cerr) const;

        double MinimizeLeadersEarning(const vector<double>& q_tariff, const int big_M, vector<double>& alpha_pl_return,
                                      vector<double>& alpha_neg_return, vector<map<int, double>>& beta_return, vector<vector<double>>& x_flow_return,
                                      std::ostream& os = std::cerr);

        double FindingTariffWithFiniteUtilities(vector<double>& q_tariff, const int big_M, std::ostream& os = std::cerr);

        bool SubstantiallyDifferentyUtility(double delta, int index_of_utility);

        double MoveInUtilitySpace(const vector<vector<double>>& x_flow, int index_of_utility, NumMatrix &u_sol, std::ostream& os = std::cerr);

	void UtilityMovingIfDifferent(vector<vector<double>> &x_flow, std::ostream &os = std::cerr);

        //For Constructors
        void CreateRandomPeoplePaths(int people, int n);


        void PolyhedronPrices(PolyCreator metad);  //Creates a random polyhedron for polyhedra_q_ (prices)
        void PolyhedronUtility(PolyCreator metad);

       public:
        Paths(const int people, const int n, const double erdos_p, PolyCreator prices_metad, PolyCreator utility_metad);

        Paths(std::istream& is);  //You read the input, more specifically D graph, paths of the people, arc_cost_buy_p, Q polyhedra, U polyhedra

        ~Paths();

        void GenerateProblem(const int seed);

        void FindingOptimalCost(std::ostream& os = std::cerr);

        void PrintData(std::ostream& os = std::cerr) const;

        void PrintDataRaw(std::ostream& os) const;

        void SaveGenerated(std::ostream& os) const;
};

#endif  //ROBUST_PATH

#ifndef UTILITY_TOOLS
#define UTILITY_TOOLS
using ll = long long int;

const long long int INF = std::numeric_limits<long long int>::max();
const bool DEBUG = true;

#define all(x) begin(x), end(x)
#define FOR(i, n) for (int i = 0; i < (n); ++i)
#define FORO(i, n) for (int i = 1; i < (n); ++i)

template <typename C>
void Print_vector(const C& Original, std::ostream& os = std::cerr) {
        for (const auto& v : Original) {
                os << v << " ";
        }
        os << endl;
}

template <typename C>
void Print_Matrix(const C& M, std::ostream& os = std::cerr) {
        for (auto& v : M) {
                os << v.size() << " ";
                for (auto& u : v) {
                        os << u;
                        os << " ";
                }
                os << endl;
        }
}

template <typename C>
void Print_MatrixClear(const C& M, std::ostream& os = std::cerr) {
        for (auto& v : M) {
                for (auto& u : v) {
                        os << u << " ";
                }
                os << endl;
        }
}

template <class T, class C>
void Print_pair(const pair<T, C>& M, std::ostream& os = std::cerr) {
        os << "(" << M.first << " , " << M.second << " ) ";
}

template <class C>
void Print_vector_pairs(const C& Original, std::ostream& os = std::cerr) {
        for (const auto& v : Original) {
                Print_pair(v, os);
        }
        os << endl;
}

template <class C>
void Print_vector_pairs_raw(const C& Original, std::ostream& os = std::cerr) {
        for (const auto& v : Original) {
                os << v.first << " " << v.second << " ";
        }
        os << endl;
}
#endif  //UTILITY_TOOLS
