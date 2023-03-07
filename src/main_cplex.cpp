#include <bits/stdc++.h>
#include <random>
#include <lemon/list_graph.h>
#include <lemon/time_measure.h>
#include <lemon/dijkstra.h>
#include <lemon/adaptors.h>
#include <lemon/concepts/path.h>
#include <lemon/concepts/graph.h>

#include "robust_energy_cplex.h"


using namespace std;
using namespace lemon;


int main(int argc, char** argv) {
	#if _GENERATE
		PolyCreator prices_metad{4, 0, 0.1, 3, 5};
		
		PolyCreator utility_metad = prices_metad;
		Paths Test(2, 10, 0.1, prices_metad, utility_metad);
		//ofstream fout("output_robust_path.txt");
		//Test.SaveGenerated(fout);
		//fout.close();
		Test.PrintData();
		Test.FindingOptimalCost();

	#endif

	#if ! defined(_GENERATE)
	assert(argc >= 2);
	ifstream fin(argv[1]);
	Paths Test(fin);
	fin.close();
	Test.FindingOptimalCost();
	#endif
}
