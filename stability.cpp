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

double analize(int i)
{
	char buf[100];
	sprintf(buf, "Amp_r%.2f.txt", i/100.0);
	FILE *f= fopen(buf, "r");

	double aux1,aux2,aux3,aux4,aux5;
	double amp;
	while(fscanf(f, "%lf %lf %lf %lf %lf", &aux1,&aux2,&aux3,&aux4,&aux5) != EOF)
	{
		amp=aux5;
	}

	fclose(f);
	return(amp);
}


int main()
{
	using namespace boost::numeric::odeint;

    double w;
    printf("w: ");
    std::cin >>w; 

	char buf[100];
	sprintf(buf, "stability_w%.2f.txt", w);
    double amp;


	FILE *f= fopen(buf, "w");

	for (int i = 0; i < 100; ++i)
	{
		amp=analize(i);
		if(amp>2)
		{
		fprintf(f, "%lf   %lf   %lf   %d   \n",i/100.0,amp,w,1);
		}
		if(amp<2)
		{
		fprintf(f, "%lf   %lf   %lf   %d   \n",i/100.0,amp,w,0);
		}


	}
	
	return 0;
}