#include "robust_energy_cplex.h"
#include <random>
#include <lemon/list_graph.h>
#include <lemon/time_measure.h>
#include <lemon/dijkstra.h>
#include <lemon/adaptors.h>
#include <lemon/concepts/path.h>
#include <lemon/concepts/graph.h>



using namespace lemon;


int main(int argc, char** argv) {
	
	if(argc < 2) {
		PolyCreator prices_metad{4, 0, 0.1, 3, 5};
		
		PolyCreator utility_metad = prices_metad;
		Paths Test(3, 10, 0.1, prices_metad, utility_metad);
		//ofstream fout("output_robust_path.txt");
		//Test.SaveGenerated(fout);
		//fout.close();
		Test.PrintData();
		Test.FindingOptimalCost();
		Test.PrintData();
		cerr << "\n\nSaveGenerated: " << endl;

		ofstream fout("generated_test.txt");
		Test.SaveGenerated(fout);
		fout.close();
	}
	else {
		assert(argc >= 2);
		ifstream fin(argv[1]);
		Paths Test(fin);
		fin.close();
		Test.FindingOptimalCost();
		Test.PrintData();
	}
}
