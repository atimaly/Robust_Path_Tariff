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
	assert(argc >= 2);
	#if _GENERATE
		Paths Test(1, 10, 0.1);
		ofstream fout("output_robust_path.txt");
		Test.SaveGenerated(fout);
		fout.close();
	#endif
	ifstream fin(argv[1]);
	Paths Test(fin);
	fin.close();
	Test.FindingOptimalCost();
	Test.PrintData();
	//Test.PrintData();
	
	/*Paths Test(1, 10, 0.1);
	Test.PrintData();
	Test.FindingOptimalCost(4, 0.3, 5, 4, 0.2, 10, std::cerr);*/
	/*
	ofstream fout("out_robust_data.txt");
	Test.PrintData(fout);
	fout.close();
	*/
}
