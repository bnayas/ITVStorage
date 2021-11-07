#include "Species.h"
#include "System.h"
#include "Identity.h"

Species::Species(System* sys, unsigned int n0 ) : 
	n(n0), birth(sys->t)
{
	this->sys = sys;
	fit = sys->gen_fit();
	m = fit;
	sys->sr++;
	sys->Sp.push_back(this);
	sys->addToCdf(n, fit);
//	father = nullptr;
}
Species::Species(System* sys, Species* s, unsigned int n0 ) : 
	n(n0), birth(sys->t)
{
//	father = s;
	lineage = s->lineage;
	lineage->change_n(1);
	this->sys = sys;
	fit = sys->gen_fit();
	m = fit;
	lineage->change_m(fit);
	Sp_id = sys->Sp.size();
	//kid_id = s->kids.size();
	sys->sr++;
	sys->Sp.push_back(this);
	//s->kids.push_back(this);
	sys->addToCdf(n, fit);
}
Species::Species(System* sys, Identity* id, unsigned int n0) :
	n(n0), birth(sys->t)
{
	//	father = s;
	lineage = id;
	lineage->change_n(1);
	this->sys = sys;
	fit = sys->gen_fit();
	m = fit;
	lineage->change_m(fit);
	Sp_id = sys->Sp.size();
	//kid_id = s->kids.size();
	sys->sr++;
	sys->Sp.push_back(this);
	//s->kids.push_back(this);
	sys->addToCdf(n, fit);
}
void Species::Revive(System* sys, Species* s, unsigned int n0) 
{
	n = n0;
	birth = sys->t;
	//father = s;
	lineage = s->lineage;
	lineage->change_n(1);
	this->sys = sys;
	fit = sys->gen_fit();
	m = fit;
	lineage->change_m(fit);
	Sp_id = sys->Sp.size();
	//kid_id = s->kids.size();
	sys->sr++;
	sys->Sp.push_back(this);
	//s->kids.push_back(this);
	sys->addToCdf(n, fit);
}

void Species::Revive(System* sys, Identity* id, unsigned int n0)
{
	n = n0;
	birth = sys->t;
	//father = s;
	lineage = id;
	lineage->change_n(1);
	this->sys = sys;
	fit = sys->gen_fit();
	m = fit;
	lineage->change_m(fit);
	Sp_id = sys->Sp.size();
	//kid_id = s->kids.size();
	sys->sr++;
	sys->Sp.push_back(this);
	sys->addToCdf(n, fit);
}


Species::~Species()
{
	/*if (father)
	{
		if (kid_id == father->kids.size() - 1)
		{
			father->kids.pop_back();
		}
		else
		{
			auto it = father->kids.begin() + kid_id;
			*it = std::move(father->kids.back());
			father->kids.pop_back();
			father->kids[kid_id]->kid_id = kid_id;
		}
	}*/
}


void Species::extinction()
{
	if(!sys->collect_data)
	{
		sys->collect_data = true;
		auto s = sys->Sp.begin();
		while (sys->collect_data && s != sys->Sp.end())
		{
			if ((*s)->birth == 0)
				sys->collect_data = false;
			s++;
		}
	}
//	if(father)
//		if (father->kid_num() == 0 && father->n == 0)
//			father->extinction();
//	//if (lineage->get_n() == 0)
//	//	sys->extinct_lineage();
////	delete this;

}


