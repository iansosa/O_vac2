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

typedef std::vector< double > state_type;

class harm_osc 
{

    std::vector< double >& m_G;
    std::vector< double >& m_B;
    int m_N;
    std::vector< double >& m_M;
    arma::Mat<double> &m_K;
    std::vector< double >& m_L;
    std::vector< double >& m_F;
    std::vector< double >& m_W;
    boost::mt19937 &m_rng;


	public:
    harm_osc( std::vector< double > &G,std::vector< double > &B, int N, std::vector< double > &M,arma::Mat<double> &K,std::vector< double > &L,std::vector< double > &F,std::vector< double > &W ,boost::mt19937 &rng) : m_G(G) , m_B(B) , m_N(N), m_M(M) , m_K(K) , m_L(L), m_F(F), m_W(W), m_rng(rng) { }

    void operator() ( const state_type &x , state_type &dxdt , const double t  )
    {
    	boost::uniform_real<> unif( -1, 1);//la distribucion de probabilidad uniforme entre cero y 2pi
    	boost::variate_generator< boost::mt19937&, boost::uniform_real<> > gen( m_rng , unif );//gen es una funcion que toma el engine y la distribucion y devuelve el numero random
    	double sum=0;
        for (int i = 0; i < m_N; ++i)
        {
        	sum=0;
        	for (int j = 0; j < m_N; ++j)
    		{
    			sum=sum+m_K(i,j)*(x[2*i]-x[2*j]);
    		}
    		sum=sum/(m_N*m_M[i]);


    	
    		if(t<=2000000)
    		{
        		dxdt[2*i]=x[2*i+1];
        		dxdt[2*i+1]= -(m_G[i]/m_M[i])*x[2*i+1]-(m_K(i,i)/m_M[i])*(x[2*i])-(m_B[i]/m_M[i])*pow(x[2*i],3)+m_F[i]*cos(m_W[i]*t)-sum+gen()*0.001;
    		}
    		else
    		{
        		dxdt[2*i]=x[2*i+1];
        		dxdt[2*i+1]= -(m_G[i]/m_M[i])*x[2*i+1]-(m_K(i,i)/m_M[i])*(x[2*i])-(m_B[i]/m_M[i])*pow(x[2*i],3)-sum;
    		}

        }
    }
};

struct push_back_state_and_time
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times ) : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void inicialcond(state_type &x,int N,boost::mt19937 &rng,int caso,double r)
{
    boost::uniform_real<> unif( 0, 5);//la distribucion de probabilidad uniforme entre cero y 2pi
    boost::variate_generator< boost::mt19937&, boost::uniform_real<> > gen( rng , unif );//gen es una funcion que toma el engine y la distribucion y devuelve el numero random
	char buf[100];
	sprintf(buf, "Xi_r%.2f.txt", r);
    if(caso==0)
    {
    	FILE *w= fopen(buf, "w");
    	for (int i = 0; i < N; ++i)
		{
			if(i < (int) N*r -0.1)
			{
				fprintf(w, "%f  ",5.0);
				fprintf(w, "%f\n",0.0);
			}
			else
			{
				fprintf(w, "%f  ",0.0);
				fprintf(w, "%f\n",0.0);
			}
		}
		fclose(w);
		FILE *re= fopen(buf, "r");
		for (int i = 0; i < N; ++i)
		{
			fscanf(re, "%lf", &x[2*i]); // posicion inicial i
			fscanf(re, "%lf", &x[2*i+1]); // momento inicial i
		}
		fclose(re);
    }
    if(caso==1)
    {
    	FILE *re= fopen(buf, "r");
		for (int i = 0; i < N; ++i)
		{
			fscanf(re, "%lf", &x[2*i]); // posicion inicial i
			fscanf(re, "%lf", &x[2*i+1]); // momento inicial i
		}
		fclose(re);
    }
}

void fillK(arma::Mat<double> &K,int N,boost::mt19937 &rng,double k)
{
    	for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j <= i; ++j)
			{
				if(i==j)
				{
					K(i,j)=1.0;
				}
				else
				{
					K(i,j)=k;
				}

			}
		}
    	for (int i = 0; i < N; ++i)
		{
			for (int j = N-1; j > i; --j)
			{
				K(i,j)=K(j,i);
			}
		}
}

void fillM(std::vector<double> &M,int N,boost::mt19937 &rng)
{
		for (int i = 0; i < N; ++i)
		{
			M[i]=1.0;
		}
}

void fillG(std::vector<double> &G,int N,boost::mt19937 &rng)
{

		for (int i = 0; i < N; ++i)
		{
			G[i]=0.1;
		}
}

void fillB(std::vector<double> &B,int N,boost::mt19937 &rng)
{
		for (int i = 0; i < N; ++i)
		{
			B[i]=0.1*(4.0/3.0);
		}
}

void fillL(std::vector<double> &L,int N,boost::mt19937 &rng)
{
		for (int i = 0; i < N; ++i)
		{
			L[i]=0.0;
		}
}

void fillF(std::vector<double> &F,int N,boost::mt19937 &rng)
{
		for (int i = 0; i < N; ++i)
		{
			F[i]=1.0;
		}
}

void fillFw(std::vector<double> &Fw,int N,boost::mt19937 &rng)
{
		for (int i = 0; i < N; ++i)
		{
			Fw[i]=1.5;
		}
}

double calcampi(size_t steps,std::vector< state_type > &x_vec,int j)
{
	double max=0;
	double min=0;
	for (int i = (int)steps*9.0/10.0; i < (int)steps; ++i)
	{
		if(max<=x_vec[i][2*j])
		{
			max=x_vec[i][2*j];
		}
		if(min>=x_vec[i][2*j])
		{
			min=x_vec[i][2*j];
		}
	}
	return((max-min)/2.0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double calcamp(size_t steps,std::vector< state_type > &x_vec, int N, double r, double &Amax, double &Amin)
{
	std::vector<double> Amp(N);

	for (int i = 0; i < N; ++i)
	{
		Amp[i]=calcampi(steps,x_vec,i);
	}
	double maxabs=0;
	double max=0;
	for (int i = 0; i < N; ++i)
	{
		if(maxabs<=fabs(Amp[0]-Amp[i]))
		{
			maxabs=fabs(Amp[0]-Amp[i]);
		}
		if(max<=Amp[i])
		{
			max=Amp[i];
		}
	}
	Amax=max;
	Amin=max-maxabs;
	double truer=0;
	for (int i = 0; i < N; ++i)
	{
		if(Amp[i]>=max-0.2)
		{
			truer=truer+1;
		}
	}
	truer=(double)truer/(N);
	return(truer);
}


void printci(std::vector< state_type > &x_vec,int time, int N,double r)
{
		char buf[100];
		sprintf(buf, "Xi_r%.2f.txt", r);
    	FILE *w= fopen(buf, "w");
    	for (int i = 0; i < N; ++i)
		{
			fprintf(w, "%f  ",x_vec[time][2*i] );
			fprintf(w, "%f\n",x_vec[time][2*i+1] );
		}
		fclose(w);
}

void printAmp(double k, double r, double truer, double Amax, double Amin)
{
		char buf[100];
		sprintf(buf, "Amp_r%.2f.txt", r);
    	FILE *w= fopen(buf, "a");
		fprintf(w, "%lf   %lf   %lf   %lf   %lf\n",k,r,truer,Amax,Amin );
		fclose(w);
}

double iterateK(int N, double k,boost::mt19937 &rng,int caso,double r)
{
    using namespace std;
    using namespace boost::numeric::odeint;
    arma::Mat<double> K(N,N);
    std::vector<double> M(N);
    std::vector<double> G(N);
    std::vector<double> B(N);
    std::vector<double> L(N);
    std::vector<double> F(N);
    std::vector<double> W(N);
	state_type x(2*N); //condiciones iniciales
    std::vector<state_type> x_vec;
    std::vector<double> times;
    double Amax,Amin,truer;
    ////////////////////////////////////////////////////////////////////
	fillK(K,N,rng,k);
	fillM(M,N,rng);
	fillG(G,N,rng);
	fillB(B,N,rng);
	fillL(L,N,rng);
	fillF(F,N,rng);
	fillFw(W,N,rng);
	if(caso==0)
	{
		inicialcond(x,N,rng,0,r);
	}
	else
	{
		inicialcond(x,N,rng,1,r);
	}
    ////////////////////////////////////////////////////////////////////
    harm_osc ho(G,B,N,M,K,L,F,W,rng);
    runge_kutta4 < state_type > stepper;
	size_t steps = integrate_adaptive(stepper, ho, x , 0.0 , 50.0*2.0*M_PI/1.5 , 0.01, push_back_state_and_time( x_vec , times )); //1 funcion. 2 condiciones iniciales. 3 tiempo inicial. 4 tiempo final. 5 dt inicial. 6 vector de posicion y tiempo
	printci(x_vec,steps,N,r);
	truer=calcamp(steps,x_vec,N,r,Amax,Amin);
	printAmp(k,r,truer,Amax,Amin);
	return(truer);
}

void printkcrit(double k, double r)
{
    	FILE *w= fopen("Kcrit.txt", "a");
		fprintf(w, "%lf   %lf   \n",r,k);
		fclose(w);
}


int main()
{
    using namespace std;
    using namespace boost::numeric::odeint;

    boost::mt19937 rng(static_cast<unsigned int>(std::time(0))); 

///////////////////////////////////////////////////////////////////////
    int N;
    printf("N: ");
    std::cin >>N; 

    double dk;
    printf("K dk: ");
    std::cin >>dk; 

    int R_step;
    printf("R steps: ");
    std::cin >>R_step; 

    double dr;
    printf("dr: ");
    std::cin >>dr; 

    double R_start;
    printf("R start: ");
    std::cin >>R_start; 
////////////////////////////////////////////////////////////////////////////
    double k=0;
	char buf[100];
	double truer;
	FILE *v= fopen("Kcrit.txt", "w");
    fclose(v);

	#pragma omp parallel for schedule(static) private(buf,truer,k,rng)
    for (int i = 0; i < R_step; ++i)
    {
		sprintf(buf, "Amp_r%.2f.txt", R_start+i*dr);
		FILE *w= fopen(buf, "w");
        fclose(w);
        int j=0;
        int itera=0;
    	while(1==1)
    	{
    		k=dk*j;
			truer=iterateK(N,k,rng,j,R_start+i*dr);
			if(R_start+i*dr < truer-0.00001 || R_start+i*dr > truer+0.00001 || R_start+i*dr > 1.0-0.00001 ||  R_start+i*dr < 0.00001)
			{
				if(itera==0)
				{
					printkcrit(k,R_start+i*dr);
					printf("finished r=%lf\n",R_start+i*dr);
				}
				if(itera==2)
				{
					break;
				}
				itera++;
			}
			j++;    	
		}
    }
    	return 0;

}