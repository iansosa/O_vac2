#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <string>
#include <cassert>
#include <vector>
#include <iostream>
#include <armadillo>
#include <utility>
#include <omp.h>
#include <boost/numeric/odeint.hpp>
#include <boost/random.hpp>


int main()
{
	using namespace boost::numeric::odeint;

    double amp;

	FILE *r= fopen("stability.txt", "w");
	FILE *s1= fopen("stability_s1.txt", "w");
	FILE *s0= fopen("stability_s0.txt", "w");

	double aux1,aux2,aux3,aux4;
	int aux5;

	while(fscanf(r, "%lf %lf %lf %lf %d", &aux1,&aux2,&aux3,&aux4,&aux5) != EOF)
	{
		if(aux5==0)
		{
			fprintf(s0, "%lf   %lf   %lf  \n",aux1,aux2,aux3);
		}
		if(aux5==1)
		{
			fprintf(s1, "%lf   %lf   %lf  \n",aux1,aux2,aux3);
		}
	}


	return 0;
}