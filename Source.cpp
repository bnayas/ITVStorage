#include "pch.h"
#include "Species.h"
#include "System.h"
#include "Random.h"
#include "Utility.h"
#define _CRT_SECURE_NO_WARNINGS

int create_cdf(int* nums, int* cdf, int size)
{
	int N = 0;
	for (int i = 0; i < size; i++)
	{
		*cdf = N + *nums;
		N += *nums;
		cdf++;
		nums++;
	}
	return N;
}


int main(int argc, char* argv[])
{
	param parameters(argc, argv);
	long long int print_next = parameters.print_each;
	Random rand_gen;
	System Sys(parameters.N, parameters.gama0, parameters.delta, parameters.nu);
	
	for (unsigned int j = 0; j < parameters.runs; j++)
	{
		for (int i = 0; i < Sys.size; i++)
		{
			Sys.t++;
		//	if (rand_gen.rand() < Sys.delta)
			Sys.env_change();
			if (Sys.t > print_next)
			{
				Sys.print_system();
				print_next = Sys.t + parameters.print_each;
				Sys.log_Sp();
			}
			if (Sys.Focal.get_n() == 0)
				Sys.Extinction();
			if (Sys.Focal.get_n() == Sys.size)
				Sys.Fixation();
			Sys.death_process();
			Sys.birth_process();
			Sys.update();
		
		}
	}
}