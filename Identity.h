#pragma once

class System;

class Identity
{

public:
	Identity(System *sys);
	long long int last_focal_n;
	unsigned long long int n;
	unsigned long long int prev_n;

	std::vector<long long int> times;
	std::vector<long long int> visits;
	//std::vector<long long int> x1;
	//std::vector<long long int> x2;
	//std::vector<long long int> g11;
	//std::vector<long long int> g21;
	//std::vector<long long int> g12;
	//std::vector<long long int> g22;
	//std::vector<long long int> g31;
	//std::vector<long long int> g13;
	//std::vector<long long int> xvisits;

	inline void set_n(int new_n) { 
		n = new_n; 
	};
	inline void set_m(int new_m) { 
		m = new_m; 
	};
	inline void change_n(int dn) { 
		n += dn; 
	};
	inline void change_m(double dm) { 
		m += dm; 
	};
	inline unsigned int get_n() const { return n; };
	inline double get_m() const { return m; };
	void save();
private:
	System* sys;

	double m;

	
	std::unordered_map<unsigned int , unsigned int > m_Index;
};