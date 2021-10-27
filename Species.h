#pragma once
#include "pch.h"
#include "Identity.h"

class System;


class Species
{

public:
	System* sys;
	//Species* father;
	unsigned int n;
	double m;
	//std::vector<Species* > kids;
	//unsigned int kid_id;
	unsigned int Sp_id;
	double fit;
	long long unsigned int birth;
	Identity* lineage;
	inline void set_n(unsigned int new_n) { n = new_n; m = fit * n; };
	inline void change_n(int dn) { 
		n += dn; m = fit * n; 
		lineage->change_n(dn); lineage->change_m(fit * dn);
	};

	inline unsigned int get_n() const { return n; };
	inline double get_fit() const { return fit; };
	inline double get_m() const { return m; };

	Species(System* sys, unsigned int n0 = 0);
	Species(System* sys, Species *s, unsigned int n0 = 0);
	void Revive(System* sys, Species* s, unsigned int n0 = 0);
	~Species();

	//inline int kid_num() const { return kids.size(); };

	inline void set_fit(double new_fit)  { fit = new_fit; m = fit * n; };
	inline void set_lineage(Identity* branch) { lineage = branch; };

	void extinction();


};

