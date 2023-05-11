#include <algorithm>
#include <numeric>
#include <random>

using namespace std;

#ifndef UNITVEC_GEN
#define UNITVEC_GEN

class RandomUnitVecGen {
	
	std::random_device rd{};
    	std::mt19937 gen{rd()};
	std::normal_distribution<> d{0, 1};
	
	public:

	vector<double> GenerateDUnitVec(int dim) {
		vector<double> unit_v(dim);
		std::generate(unit_v.begin(), unit_v.end(), [&](){return d(gen);}); //Filling it up
		double normi = sqrt(std::inner_product(unit_v.begin(), unit_v.end(), unit_v.begin(), 0.0L)); //The norm of the vector
		std::for_each(unit_v.begin(), unit_v.end(), [normi](double &c){c /= normi;}); //Norming to 1
		return unit_v;
	} 

	vector<double> GenerateDUnitVecPositive(int dim) {
		vector<double> unit_v(dim);
		std::generate(unit_v.begin(), unit_v.end(), [&](){return std::abs(d(gen));}); //Filling it up
		double normi = sqrt(std::inner_product(unit_v.begin(), unit_v.end(), unit_v.begin(), 0.0L)); //The norm of the vector
		std::for_each(unit_v.begin(), unit_v.end(), [normi](double &c){c /= normi;}); //Norming to 1
		return unit_v;
	}

	vector<double> GenerateDUnitVecNegative(int dim) {
		vector<double> unit_v(dim);
		std::generate(unit_v.begin(), unit_v.end(), [&](){return -std::abs(d(gen));}); //Filling it up
		double normi = sqrt(std::inner_product(unit_v.begin(), unit_v.end(), unit_v.begin(), 0.0L)); //The norm of the vector
		std::for_each(unit_v.begin(), unit_v.end(), [normi](double &c){c /= normi;}); //Norming to 1
		return unit_v;
	}
	
};
#endif //RANDOM_UNIT_VEC_GEN
