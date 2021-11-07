#pragma once
#include "pch.h"
#include "Random.h"
#include "Identity.h"
//class Identity;
class Species;

class System
{
public:
	std::vector<Species*> Sp;
	int sr;
	bool collect_data;
	unsigned int mutants;
	long long int t;
	unsigned int size;
	float nu;
	float gama0;
	float delta;
	Identity Focal;
	Identity Env;
	System(unsigned int N, float GAMMA, float DELTA, float NU);
	~System();
	Species* MRCA;

	void addToCdf(unsigned int n, double m) { 
		if (fit_cdf.size())
		{
			double m_cdf = fit_cdf.back() + m;
			fit_cdf.push_back(m_cdf);
		//	int n_cdf = grid_cdf.back() + n;
		//	grid_cdf.push_back(n_cdf);
		}
		else
		{
			fit_cdf.push_back(m);
		//	grid_cdf.push_back(n);
		}
	};
	void update_n(unsigned int i, unsigned int new_n);
	void update_fit(unsigned int i, unsigned int new_fit);
	void update_cdf(int n); 
	void env_change();
	void death_process();
	void birth_process(); 
	void mutation_process(Identity *id);
	void ForcedMutation(Identity* identity);
	void Extinction();
	void Fixation();
	double inline gen_fit() { return fit_dist(); };
	void update();
	void log_Sp();
	void print_system();
private:
	void InitFileNames();
	Random rand_gen;
	double change_num;
	unsigned int missing_deaths;
	unsigned int deaths_num;
	std::vector<Species*> death_list;
	void inline mutation_process(unsigned int idx);
	std::function<double()> fit_dist;
	std::function<double()> rand;
	std::vector<double> fit_cdf;
	//std::vector<double> grid_cdf;
	std::vector<int> new_gen;
	std::vector<Identity *> new_mutations;
	//	double* create_cdf();
	std::vector<long double>  SAD;
	std::vector<long double>  total_wins;
	std::vector<long double>  total_loses;
	std::vector<long double>  wins;
	std::vector<long double>  loses;
	std::vector<long double>  t_wins;
	std::vector<long double>  t_loses;
	std::vector<long double>  visits;
	std::vector<long double>  fit_visits;
	std::vector<long long int>  invasions;
	Identity* current_focal;
	std::vector<long double> mm1;
	std::vector<long double> mm2;

	std::vector<long double> p;
	std::vector<long double> mu;
	std::vector<long double> sig2;
	std::vector<long double> SR;
	std::vector<long double> x1;
	std::vector<long double> x2;
	std::vector<long double> g11;
	std::vector<long double> g21;
	std::vector<long double> g12;
	std::vector<long double> g22;
	std::vector<long double> SR_env;
	std::vector<long double> x1_env;
	std::vector<long double> x2_env;
	std::vector<long double> g11_env;
	std::vector<long double> g21_env;
	std::vector<long double> g12_env;
	std::vector<long double> g22_env;
	//	std::vector<long double> g31;
//	std::vector<long double> g13;
	std::vector<long long int> xvisits;

	std::string DIR = ".\\";
	std::string file_names[32];

};

template <typename T>
void print_vec(std::string file_name, std::vector<T> vec, int size)
{
	std::ofstream File(file_name);
	auto itr = vec.begin();
	for (int x = 0; x < size; x++)
	{
		File << *itr << "\n";
		++itr;
	}
	File.close();
}
