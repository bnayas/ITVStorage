#include "pch.h"
#include "Identity.h"
#include "System.h"

Identity::Identity(System* s) :
	sys(s),n(0), m(0)
{
	visits.resize(sys->size, 0);
	times.resize(sys->size, 0);
	//x1.resize(sys->size, 0);
	//x2.resize(sys->size, 0);
	//g11.resize(sys->size, 0);
	//g21.resize(sys->size, 0);
	//g12.resize(sys->size, 0);
	//g22.resize(sys->size, 0);
	//g31.resize(sys->size, 0);
	//g13.resize(sys->size, 0);
	//xvisits.resize(sys->size, 0);
}
void Identity::save()
{
	//if(times[n] == 0)
	//	times[n] = sys->t;
	times[n] += sys->t;
	visits[n] ++;
}
