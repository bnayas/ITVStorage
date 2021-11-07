#include "pch.h"
#include "System.h"
#include "Species.h"
#include "Identity.h"
#include "Utility.h"

System::System(unsigned int N, float GAMMA , float DELTA, float NU) :
  size(N), gama0(GAMMA), delta(DELTA), nu(NU), change_num(delta * N), t(0), 
	Env(this), Focal(this), collect_data(false)
{
	InitFileNames();
	fit_dist = (rand_gen.real_unif_dist(-gama0, gama0));
	rand = (rand_gen.real_unif_dist(0, 1));
	double m = 0;
	for (int i = 0; i < N; i++)
	{
		Species* s = new Species(this, 1);
		s->set_n(1);
		s->set_lineage(&Env);
		new_gen.push_back(0);
		m += s->get_m();
	}
	Env.set_m(m);
	Env.set_n(N);
	Focal.set_n(0);
	SAD.resize(size, 0);
	total_wins.resize(size, 0);
	total_loses.resize(size, 0);
	wins.resize(size, 0);
	loses.resize(size, 0);
	t_wins.resize(size, 0);
	t_loses.resize(size, 0);
	fit_visits.resize(size, 0);
	visits.resize(size, 0);
	mm1.resize(size, 0);
	mm2.resize(size, 0);
	x1.resize(size, 0);
	x2.resize(size, 0);
	g11.resize(size, 0);
	g21.resize(size, 0);
	g12.resize(size, 0);
	g22.resize(size, 0);
	//g31.resize(size, 0);
	//g13.resize(size, 0);
	SR.resize(size, 0);
	p.resize(size, 0);
	mu.resize(size, 0);
	sig2.resize(size, 0);
	xvisits.resize(size, 0);
	invasions.resize(size, 0);
	current_focal = &Focal;
	invasions[0]++;
	current_focal->last_focal_n = 1;

	x2_env.resize(size, 0);
	x1_env.resize(size, 0);
	g11_env.resize(size, 0);
	g21_env.resize(size, 0);
	g12_env.resize(size, 0);
	g22_env.resize(size, 0);
	SR_env.resize(size, 0);

}

System::~System()
{
	while (Sp.size())
	{
		delete Sp.back();
		Sp.pop_back();
	}
}

void System::update_n(unsigned int i,unsigned int new_n)
{
	Sp[i]->set_n(new_n);
}

void System::update_fit(unsigned int i, unsigned int new_fit)
{
	Sp[i]->set_fit(new_fit);
}

void System::env_change()
{
	auto s = Sp.begin();
	auto cdf_itr = fit_cdf.begin();
	double cdf = 0;
	//(*s)->set_fit(fit_dist());
	//*(x + 1) = *x + ((*s)->get_fit() * (*s)->get_n());
	//x++; s++;
	Focal.set_m(0);
	Env.set_m(0);
	long double G1 = 0, G2 = 0, G3 = 0, X1 = 0, X2 = 0, sr = 0;
	long double G1_o = 0, G2_o = 0, G3_o = 0, X1_o = 0, X2_o = 0, sr_o = 0;
	//double m1 = 0, m2 = 0;
	for (; s != Sp.end(); ++s)
	{
		double gamma_i = fit_dist();
		(*s)->set_fit(exp(gamma_i));
		if (collect_data)
		{
			if ((*s)->lineage == &Focal)
			{
				double xi = (*s)->get_n();
				X1 += xi;
				X2 += pow(xi, 2);
				//	m1 += xi * exp(gamma_i);
				G1 += gamma_i * xi;
				G2 += pow(gamma_i, 2) * xi;
				//		G3 += pow(gamma_i, 3) * xi;
				sr++;
			}
			else
			{
				double xi = (*s)->get_n();
				X1_o += xi;
				X2_o += pow(xi, 2);
				//	m1 += xi * exp(gamma_i);
				G1_o += gamma_i * xi;
				G2_o += pow(gamma_i, 2) * xi;
				//		G3 += pow(gamma_i, 3) * xi;
				sr_o++;
			}
		}
		(*s)->lineage->change_m((*s)->get_m());
		cdf += (*s)->get_m() / size;
		*cdf_itr = cdf;
		++cdf_itr;
		/*if (abs(m1 - Focal.get_m()) > 0.000001)
			std::cout << "m1 Error\n";
		if (abs(m2 - Env.get_m()) > 0.000001)
			std::cout << "m2 Error\n";*/

//		*x = *(x-1) + ((*s)->get_fit() * (*s)->get_n());
//		x++;
	}
	if (collect_data)
	{
		int x = Focal.get_n();
		if (x > 0)
		{
			g11[x - 1] += (G1 / x);
			g12[x - 1] += pow(G1 / x, 2);
			//		g13[x - 1] += pow(G1 / x, 3);
			g21[x - 1] += ((double)G2 / (double)x);
			//		g31[x - 1] += (G3 / x);
			x1[x - 1] += (X1 / x);
			x2[x - 1] += X2 / pow(x, 2);
			SR[x - 1] += sr;
			double m1 = Focal.get_m(), m2 = Env.get_m();
			double N1 = (double)x, N2 = (double)size - (double)x;
			double p0 = ((m1 / N1) - (m2 / N2)) / (Focal.get_m() + Env.get_m());
			p[x - 1] += p0;
			fit_visits[x - 1] ++;
			/*if (abs((size * x1[x - 1] / xvisits[x - 1]) - (xvisits[x - 1]*x/SR[x-1])) > 0.001)
				std::cout << "Error";*/
		}
		int z = size - x;
		if (z > 0)
		{
			g11_env[z - 1] += (G1_o / z);
			g12_env[z - 1] += pow(G1_o / z, 2);
			g21_env[z - 1] += ((double)G2_o / (double)z);
			x1_env[z - 1] += (X1_o / z);
			x2_env[z - 1] += X2_o / pow(z, 2);
			SR_env[z - 1] += sr_o;
			/*if (abs(_env(size * x1[x - 1] / xvisits[x - 1]) - (xvisits[x - 1]*x/SR[x-1])) > 0.001)
				std::c_envout << "Error";*/
		}
	}
	while (fit_cdf.size() > Sp.size())
		fit_cdf.pop_back();

}

//void System::update_cdf(int n)
//{
//	auto fit_cdf_itr = fit_cdf.begin() + n; auto cdf_itr = grid_cdf.begin() + n;
//	for (auto s = Sp.begin() + n ; s != Sp.end(); ++s)
//	{
//		*(fit_cdf_itr + 1) = *fit_cdf_itr + ((*s)->get_fit() * (*s)->get_n());
//		*(cdf_itr + 1) = *cdf_itr + (*s)->get_n();
//		fit_cdf_itr++;
//		cdf_itr++;
//	}
//}

void System::death_process()
{
//	std::cout << "before death Sp: " << Sp.size() << " , new_gen: " << new_gen.size() << " , death_list: " << death_list.size() << " , mutations: " << new_mutations.size() << "\n";
	//int sum_n = 0;
	//for (auto s = Sp.begin(); s != Sp.end(); s++)
	//	sum_n += (*s)->get_n();
	//if (sum_n != size)
	//	std::cout << "Error";

	deaths_num = 0;
	auto  s = Sp.begin();
	for (unsigned int idx = 0; s != Sp.end(); ++idx)
	{
		for (int i = 0; i < (*s)->get_n(); i++)
		{
			if (rand_gen.rand() < delta)
				//if (rand_gen.rand() < delta)
				new_gen[idx] --;
		}
		deaths_num -= new_gen[idx];
	//	idx++;
		++s;
	}
	
	//std::cout << "Death count:\n";
	s = Sp.begin();
	auto dn = new_gen.begin();
	for (; dn != new_gen.end(); ++dn)
	{
		int n = (*s)->get_n();
		//std::cout << n << " , " << *dn << "\n";
		if (n  + *dn < 0)
			std::cout << " Strange";
		s++;
	}
	//std::cout << "\n";
	int sum_n = 0;
	/*dn = new_gen.begin();
	*/
	for (s = Sp.begin(); s != Sp.end(); s++)
	{
		sum_n += (*s)->get_n() ;
//		dn++;
	}
	if (sum_n != size)
		std::cout << "Error";

	//std::cout << "after death Sp: " << Sp.size() << " , new_gen: " << new_gen.size() << " , death_list: " << death_list.size() << "\n";
}

void System::birth_process()
{
	//std::cout << "before birth Sp: " << Sp.size() << " , new_gen: " << new_gen.size() << " , death_list: " << death_list.size() << " , mutations: " << new_mutations.size() << "\n";

	auto max = fit_cdf.back();
	unsigned int idx = 0;
	mutants = 0;
	int sum_n = 0;
	for (auto s = Sp.begin(); s != Sp.end(); s++)
		sum_n += (*s)->get_n();
	if (sum_n != size)
		std::cout << "Error";

	//double s = 0 , times = 100000000;
	//for (int i = 0; i < times; i++)
	////	if (rand_gen.rand() < nu)
	//	if (rand() < nu)
	//		s++;
	//std::cout << s / times << " , " << nu << "\n";
	/*for (auto f : fit_cdf)
		std::cout << f << " ";*/
	//std::cout << "\n";
	for (int i = 0;   i < deaths_num; i++)
	{
//		idx = find_in_sorted(rand_gen.rand(), fit_cdf, Sp.size()-mutants);
		idx = find_in_sorted(fit_cdf.back()*rand_gen.rand(), fit_cdf, Sp.size());
		if (rand() < nu)
			//if (rand_gen.rand() < nu)
		{
			new_mutations.push_back(Sp[idx]->lineage);
//			mutation_process(idx);
			mutants++;
		}
		else
			new_gen[idx] ++;

	}
	//std::cout << "Death count (suppose to be  " << change_num - missing_deaths << " ):\n";
	auto s = Sp.begin();
	sum_n = 0;
	auto dn = new_gen.begin();
	for (; dn != new_gen.end(); ++dn)
	{
		sum_n += (*s)->get_n() + *(dn);
		int n = (*s)->get_n();
	//	std::cout << n << " , " << *dn <<"\n";
		if (n + *dn < 0)
			std::cout << "Strange";
		s++;
	}
	if (sum_n + mutants != size)
		std::cout << "Error";
//	std::cout << "\n";
//	std::cout << "after birth Sp: " << Sp.size() << " , new_gen: " << new_gen.size() << " , death_list: " << death_list.size() << " , mutations: " << new_mutations.size() << "\n";

}

void System::mutation_process(unsigned int idx)
{
	Species* s;
	if(death_list.size() < 1)
		s = new Species(this, Sp[idx], 1);
	else
	{
		s = death_list.back();
		s->Revive(this, Sp[idx], 1);
		death_list.pop_back();
	//	s->set_n(1);
	//	s->set_fit(exp(gen_fit()));
	//	s->Sp_id = Sp.size();
	//	Sp.push_back(s);
	//	s->kid_id = Sp[idx]->kids.size();
	//	sr++;
	////	Sp[idx]->kids.push_back(s);
	//	addToCdf(1, s->get_fit());
	//	s->birth = t;
	//	s->set_lineage(Sp[idx]->lineage);
	//	s->lineage->change_n(1);
	}
	new_gen.push_back(0);
//	std::cout << t <<": Mutation procees: " << idx << " " << s->get_n() << "\n";
}

void System::mutation_process(Identity *id)
{
	Species* s;
	if (death_list.size() < 1)
	{
		s = new Species(this, id, 1);
	//	std::cout << "New species \n";
	}
	else
	{
		s = death_list.back();
		s->Revive(this, id, 1);
		death_list.pop_back();
		//	s->set_n(1);
		//	s->set_fit(exp(gen_fit()));
		//	s->Sp_id = Sp.size();
		//	Sp.push_back(s);
		//	s->kid_id = Sp[idx]->kids.size();
		//	sr++;
		////	Sp[idx]->kids.push_back(s);
		//	addToCdf(1, s->get_fit());
		//	s->birth = t;
		//	s->set_lineage(Sp[idx]->lineage);
		//	s->lineage->change_n(1);
	}
	new_gen.push_back(0);
	//	std::cout << t <<": Mutation procees: " << idx << " " << s->get_n() << "\n";
}

void System::update()
{
//	std::cout << "before update Sp: " << Sp.size() << " , new_gen: " << new_gen.size() << " , death_list: " << death_list.size()  << " , mutations: " << new_mutations.size() << "\n";

	auto dn = new_gen.begin();
	unsigned int n = 0;
	int i = 0;
	//	int j = 0;
	auto s = Sp.begin();
	//std::vector<Species*> deads;
	for ( ; s !=  Sp.end() ; )
	{	
		(*s)->change_n(*dn);
		*dn = 0;

		n = (*s)->get_n();

		if (n < 1)
		{
	//		deads.push_back(*s);
		/*	std::cout << t << ": Death process of " << i << " \n ";
			for (auto verf_itr = Sp.begin(); verf_itr != Sp.end(); ++verf_itr)
				std::cout << (*verf_itr)->get_n() << " ";
			std::cout << "\n";*/
		/*	if ((*s)->kid_num() == 0)
				(*s)->extinction();*/
			(*s)->extinction();
			if (s != Sp.end() - 1)
			{

				// std::cout << "Not last. \n ";
				//*m_itr = std::move(fit_cdf.back());
				//*n_itr = std::move(grid_cdf.back());
				int sp_id = (*s)->Sp_id;
				death_list.push_back(*s);
				*dn = std::move(new_gen.back());
				*s = std::move(Sp.back());
				(*s)->Sp_id = sp_id;
			//	grid_cdf.pop_back();
			//	fit_cdf.pop_back();
				Sp.pop_back();
				new_gen.pop_back();
				sr--;
			//	mutants--;

			}
			else
			{
			//	grid_cdf.pop_back();
			//	fit_cdf.pop_back();
				death_list.push_back(Sp.back());
				Sp.pop_back();
				new_gen.pop_back();
				sr--;

				break;
			}
			//if(s != Sp.begin())
			//	s--;
			
			
			/*for (auto verf_itr = Sp.begin(); verf_itr != Sp.end(); ++verf_itr)
				if ((*verf_itr)->get_n() > size)
					std::cout << "Deleted species \n";*/
		}
		else
		{
		//	cum_n += n;
		//	sum_m += (*s)->get_m();
		//	*n_itr = cum_n;
		//	*m_itr = cum_m;
			dn++;  
			//m_itr++;  n_itr++;  
			s++;
			i++;
		}
//		j++;
	}
	//for (auto d : deads)
	//{
	//	bool found = false;
	//	for (auto sd = death_list.begin(); sd != death_list.end(); sd++)
	//	{
	//		if (*sd == d)
	//			found = true;
	//	}
	//	if (!found)
	//		std::cout << "Error";
	//}
	//std::cout << "before mutations birth Sp: " << Sp.size() << " , new_gen: " << new_gen.size() << " , death_list: " << death_list.size() << " , mutations: " << new_mutations.size() << "\n";
	while ( new_mutations.size() > 0)
	{
		Identity* sp = new_mutations.back();
		mutation_process(sp);
		new_mutations.pop_back();
		mutants--;
	}
//	std::cout << "after mutations birth Sp: " << Sp.size() << " , new_gen: " << new_gen.size() << " , death_list: " << death_list.size() << " , mutations: " << new_mutations.size() << "\n";

	s = Sp.begin();
	double cum_m = 0;
	//m_itr = fit_cdf.begin();
	//for (auto m = fit_cdf.begin(); m != fit_cdf.end(); ++m)
	//{

	//	cum_m += (*s)->get_m() / sum_m;
	//	*m = cum_m;
	//	++s;
	//}
//	if (cum_n != size || grid_cdf.back() != size  )
	//	std::cout << "Error";
	
	int n1 = Focal.get_n();
	int n2 = Env.get_n();
	if (collect_data)
	{
		
		mu[Focal.prev_n] += (double)n1 - Focal.prev_n;
		sig2[Focal.prev_n] += pow((double)n1 - Focal.prev_n, 2);
		xvisits[Focal.prev_n] ++;
		Focal.prev_n = n1;
		mu[Env.prev_n] += (double)n2 - Env.prev_n;
		sig2[Env.prev_n] += pow((double)n2 - Env.prev_n, 2);
		xvisits[Env.prev_n] ++;
		Env.prev_n = n2;
	}
	int sum_n1 = 0;
	for (s = Sp.begin(); s != Sp.end(); s++)
	{
		if ((*s)->get_n() == 0)
			std::cout << "Error: Dead\n";
		if ((*s)->lineage == &Focal)
			sum_n1 += ((*s)->get_n());
	}
	if (n1 != sum_n1)
		std::cout << "Error: x problem\n";
	int sum_l = n1 + n2;
	if (current_focal->get_n() > current_focal->last_focal_n)
	{
		for (int inv_idx = current_focal->last_focal_n; inv_idx < current_focal->get_n(); inv_idx++)
			invasions[inv_idx]++;
		current_focal->last_focal_n = current_focal->get_n();
	}
	if (n1 == size)
		Fixation();
	else if (n1 > 0 && collect_data)
	{
		Focal.save();
		double df = ((Focal.get_m()) / n1) - ((Env.get_m() / n2));
		mm1[n1 - 1] += df;
		mm2[n1 - 1] += pow(df,2);
		visits[n1 - 1] ++;
		if (n1 < size)
			SAD[n1 - 1] ++;
	}
	if (sum_l != size || n1 > size || n2 > size )
		std::cout << "Error: sum";
	/*if(Sp.size() > size)
		std::cout << "Error";*/
	

}


void System::log_Sp()
{
	std::cout << t << ":\n";
	std::cout << Focal.get_n() << " , " << Env.get_n() << " :  ";
	int n = 0;
	for (auto verf_itr = Sp.begin(); verf_itr != Sp.end(); ++verf_itr)
	{
		n = (*verf_itr)->get_n();
		std::cout << n << " ";
	}
	std::cout << "\n";
}


void System::ForcedMutation(Identity *identity)
{
	Species* s_max = Sp.back();
	int idx = 0, max = 0, i = 0;;
	for (auto s : Sp)
	{
		if (s->get_n() > max)
		{
			max = s->get_n();
			s_max = s;
		}
		i++;
	}
	s_max->change_n(-1);
	Species* s;
	if(death_list.size() == 0)
		 s = new Species(this, 1);
	else
	{
		s = death_list.back();
		s->Revive(this, identity, 1);
		death_list.pop_back();
	}
	identity->set_n(1);
	identity->set_m( s->get_m());
	s->set_lineage(identity);
	current_focal = identity;
	invasions[0]++;
	current_focal->last_focal_n = 1;
	Focal.prev_n = Focal.n;
	Env.prev_n = Env.n;
	new_gen.push_back(0);
}


void System::Extinction()
{
	if (collect_data)
	{
		auto lose_itr = loses.begin();
		auto tlose_itr = t_loses.begin();
		auto total_lose_itr = total_loses.begin();

		auto s_times_itr = Focal.times.begin();
		auto s_visits_itr = Focal.visits.begin();

		for (int i = 0; i < size; i++)
		{
			if (*s_times_itr > 0)
			{
				//(*lose_itr)++;
				//(*tlose_itr) += t - (*s_times_itr);
				(*lose_itr) += (*s_visits_itr);
				(*tlose_itr) += t * (*s_visits_itr) - (*s_times_itr);
				*total_lose_itr += (*s_visits_itr);
				(*s_times_itr) = 0;
				(*s_visits_itr) = 0;
			}
			++lose_itr;
			++tlose_itr;
			++total_lose_itr;
			++s_times_itr;
			++s_visits_itr;
		}
	}
	else
	{
		collect_data = true;
		auto s = Sp.begin();
		while (collect_data && s != Sp.end() )
		{
			if ((*s)->birth == 0)
				collect_data = false;
			s++;
		}
	}
	ForcedMutation(&Focal);
}
void System::Fixation()
{
	if (collect_data)
	{
		auto win_itr = wins.begin();
		auto twin_itr = t_wins.begin();
		auto total_win_itr = total_wins.begin();

		auto s_times_itr = Focal.times.begin();
		auto s_visits_itr = Focal.visits.begin();

		for (int i = 0; i < size; i++)
		{
			if (*s_times_itr > 0)
			{
				//(*twin_itr) += t - (*s_times_itr);
				//(*win_itr)++;
				(*twin_itr) += (*s_visits_itr) * t - (*s_times_itr);
				(*win_itr) += (*s_visits_itr);

				*total_win_itr += (*s_visits_itr);
				(*s_times_itr) = 0;
				(*s_visits_itr) = 0;
			}
			++win_itr;
			++twin_itr;
			++total_win_itr;
			++s_times_itr;
			++s_visits_itr;
		}
	}
	ForcedMutation(&Env);
}


void System::print_system() {
//	std::string NAME[22] = { "mm1","mm2","visits",
//"Wins","Loses","TWins","TLoses","Pi","SAD","SR","x1","x2","g11","g12","g13","g21","g22","g31","p","mu","sig2","invasions","xvisits","fvisits" };
	if (collect_data)
	{
		if (death_list.size() > 10)
		{
			for (int i = 0; i < floor(death_list.size() / 2); i++)
			{
				delete death_list.back();
				death_list.pop_back();
			}
		}
		print_vec(file_names[0], mm1, size);
		print_vec(file_names[1], mm2, size);
		print_vec(file_names[2], visits, size);
		print_vec(file_names[3], wins, size);
		print_vec(file_names[4], loses, size);
		print_vec(file_names[5], t_wins, size);
		print_vec(file_names[6], t_loses, size);
		//	print_vec(file_names[7], mm1, size);
		print_vec(file_names[8], SAD, size);
		print_vec(file_names[9], SR, size);
		print_vec(file_names[10], x1, size);
		print_vec(file_names[11], x2, size);
		print_vec(file_names[12], g11, size);
		print_vec(file_names[13], g12, size);
		//print_vec(file_names[14], g13, size);
		print_vec(file_names[15], g21, size);
		print_vec(file_names[16], g22, size);
		//print_vec(file_names[17], g31, size);
		print_vec(file_names[18], p, size);
		print_vec(file_names[19], mu, size);
		print_vec(file_names[20], sig2, size);
		print_vec(file_names[21], invasions, size);
		print_vec(file_names[22], xvisits, size);
		print_vec(file_names[23], fit_visits, size);
		print_vec(file_names[24], SR_env, size);
		print_vec(file_names[25], x1_env, size);
		print_vec(file_names[26], x2_env, size);
		print_vec(file_names[27], g11_env, size);
		print_vec(file_names[28], g12_env, size);
		print_vec(file_names[29], g21_env, size);
		print_vec(file_names[30], g22_env, size);

		/*std::ofstream m1Fil(file_names[0]);
		std::ofstream m2File(file_names[1]);
		std::ofstream vFile(file_names[2]);
		std::ofstream winsFile(file_names[3]);
		std::ofstream losesFile(file_names[4]);
		std::ofstream twinsFile(file_names[5]);
		std::ofstream tlosesFile(file_names[6]);
		std::ofstream piFile(file_names[7]);
		std::ofstream SADFile(file_names[8]);
		std::ofstream SRFile(file_names[9]);
		std::ofstream x1File(file_names[10]);
		std::ofstream x2File(file_names[11]);
		std::ofstream g11File(file_names[12]);
		std::ofstream g21File(file_names[13]);
		std::ofstream g12File(file_names[14]);
		std::ofstream g22File(file_names[15]);
		std::ofstream g31File(file_names[16]);
		std::ofstream g13File(file_names[17]);
		std::ofstream pFile(file_names[18]);
		std::ofstream muFile(file_names[19]);
		std::ofstream sig2File(file_names[20]);
		std::ofstream invasionsFile(file_names[21]);
		auto m1_itr = mm1.begin();
		auto m2_itr = mm2.begin();
		auto visits_itr = visits.begin();
		auto win_itr = wins.begin();
		auto lose_itr = loses.begin();
		auto twin_itr = t_wins.begin();
		auto tlose_itr = t_loses.begin();
		auto totalwin_itr = total_wins.begin();
		auto totallose_itr = total_loses.begin();
		auto SAD_itr = SAD.begin();
		auto p_itr = p.begin();
		auto SR_itr = SR.begin();
		auto x1_itr = x1.begin();
		auto x2_itr = x2.begin();
		auto mu_itr = mu.begin();
		auto sig2_itr = sig2.begin();
		auto invasions_itr = invasions.begin();

		auto g11_itr = g11.begin();
		auto g21_itr = g21.begin();
		auto g12_itr = g12.begin();
		auto g22_itr = g22.begin();
		auto g31_itr = g31.begin();
		auto g13_itr = g13.begin();
		auto xvisits_itr = xvisits.begin();
		auto fit_visits_itr = fit_visits.begin();
		for (int x = 0; x < size; x++)
		{
			double v_fit = *fit_visits_itr;
			double x_vis = *xvisits_itr;
			m1File << *m1_itr << "\n";
			m2File << *m2_itr << "\n";
			vFile << *visits_itr << "\n";
			winsFile << *win_itr << "\n";
			losesFile << *lose_itr << "\n";
			twinsFile << *twin_itr << "\n";
			tlosesFile << *tlose_itr << "\n";
			SADFile << *SAD_itr << "\n";
			piFile << ((double) *totalwin_itr) / (*totalwin_itr + *totallose_itr) << "\n";
			SRFile << *SR_itr / v_fit << "\n";
			x1File << *x1_itr / v_fit << "\n";
			x2File << *x2_itr / v_fit << "\n";
			invasionsFile << *invasions_itr << "\n";
			muFile <<   *mu_itr  / x_vis << "\n";
			sig2File << *sig2_itr / x_vis << "\n";
			g11File << *g11_itr / v_fit << "\n";
			g21File << *g21_itr / v_fit << "\n";
			g12File << *g12_itr / v_fit << "\n";
			g22File << *g22_itr / v_fit << "\n";
			g31File << *g31_itr / v_fit << "\n";
			g13File << *g13_itr / v_fit << "\n";
			pFile << *p_itr / v_fit << "\n";

			++mu_itr;
			++sig2_itr;
			++invasions_itr;
			++SAD_itr;
			++SR_itr ;
			++x1_itr ;
			++x2_itr ;
			++g11_itr;
			++g21_itr;
			++g12_itr;
			++g22_itr;
			++g31_itr;
			++g13_itr;
			++xvisits_itr;
			++m1_itr;
			++m2_itr;

			++visits_itr;
			++win_itr;
			++lose_itr;
			++twin_itr;
			++tlose_itr;
			++totalwin_itr;
			++totallose_itr;
			++SAD_itr;
			++p_itr;

		}
		m1File.close();
		m2File.close();
		vFile.close();
		winsFile.close();
		losesFile.close();
		twinsFile.close();
		tlosesFile.close();
		piFile.close();
		SADFile.close();
		SRFile.close();
		x1File.close();
		x2File.close();
		g11File.close();
		g21File.close();
		g12File.close();
		g22File.close();
		g31File.close();
		g13File.close();
		pFile.close();
		invasionsFile.close();
		muFile.close();
		sig2File.close();*/
	}
};

void System::InitFileNames()
{
	std::string NAME[31] = { "mm1","mm2","visits",
"Wins","Loses","TWins","TLoses","Pi","SAD","SR","x1","x2","g11","g12","g13","g21","g22","g31","p","mu","sig2","invasions","xvisits","fvisits","SR_env","x1_env","x2_env","g11_env","g12_env","g21_env","g22_env" };

	std::string suffix = "_N" + std::to_string(size) +
		"gama" + round(gama0, 4) + "delta" + round(delta, 4) +
		"nu" + round(nu, 5) + ".dat";
	file_names[0] = DIR + NAME[0] + suffix;

	int x = 0;

	if (fileExists(file_names[0].c_str()))
	{
		while (!(fileExists(file_names[0].c_str())))
		{
			file_names[0] = DIR + NAME[0] + suffix;
			x++;
		}
		for (int i = 1; i < 31; i++)
			file_names[i] = DIR + NAME[i] + suffix;
	}
	else
		for (int i = 0; i < 31; i++)
			file_names[i] = DIR + NAME[i] + suffix;

}