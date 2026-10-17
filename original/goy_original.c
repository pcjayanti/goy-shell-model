/*
 * GOY shell model for 2D turbulence -- ORIGINAL THESIS CODE (2016)
 *
 * Transcribed verbatim from Appendix A of:
 *   "Study of fluid turbulence using the GOY shell model and semi-analytical
 *    modelling of nonlinear interactions", IIT Madras, May 2016.
 *
 * DO NOT FIX ANYTHING IN THIS FILE. It is the reference baseline used to
 * reproduce the 2016 results. All known defects are preserved deliberately
 * and marked [BUG-n] below. Fixes belong in the C++ rewrite.
 *
 * Deviations from the printed listing (unavoidable, mechanical):
 *   - line-wrap artifacts from the PDF rejoined (e.g. "v[i-\n1]" -> "v[i-1]")
 *   - page numbers / form feeds removed
 *   - the `iter[]` declaration is uncommented: it appears as "//int\niter[]={...}"
 *     in the thesis, which cannot compile, yet main() uses iter[] and
 *     sizeof(iter). The printed values are used as-is.
 *
 * Build:  gcc -O3 -o goy_original goy_original.c -lm
 * Run:    ./goy_original            (appends to ./energyspectrum)
 */

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>

double nu=0.0000000001,nu1=0.009,f1=0.002,f2=0.002,k0=0.0009765625,c0=0.0009765625,dt=0.0001;

int N=29,q=2;
int iter[]={5000000,10000000,15000000,20000000,25000000,30000000,40000000,50000000,60000000,80000000,100000000,120000000,240000000,480000000};

double f2D(int j,int i,double v[N],double k[N],double c[N][3],int delta1[N],int delta2[N],int theta2D[N])
{
	double ans;
	double a=0;
	if(i==1)
	{
		ans=-nu*pow(k[i],2)*v[i]-nu1*theta2D[i]*v[i]/k[i]+c[i][0]*v[i+1]*v[i+2]+c[i][1]*v[i-1]*v[i+1]+f1*delta1[i]+f2*delta2[i];
	}

	else if(i==0)
	{
		ans=-nu*pow(k[i],2)*v[i]-nu1*theta2D[i]*v[i]/k[i]+c[i][0]*v[i+1]*v[i+2]+f1*delta1[i]+f2*delta2[i];
	}

	else if(i==(N-2))
	{
		ans=-nu*pow(k[i],2)*v[i]-nu1*theta2D[i]*v[i]/k[i]+c[i][1]*v[i-1]*v[i+1]+c[i][2]*v[i-1]*v[i-2]+f1*delta1[i]+f2*delta2[i];
	}

	else if(i==(N-1))
	{
		ans=-nu*pow(k[i],2)*v[i]-nu1*theta2D[i]*v[i]/k[i]+c[i][2]*v[i-1]*v[i-2]+f1*delta1[i]+f2*delta2[i];
	}

	else
	{
		ans=-nu*pow(k[i],2)*v[i]-nu1*theta2D[i]*v[i]/k[i]+c[i][0]*v[i+1]*v[i+2]+c[i][1]*v[i-1]*v[i+1]+c[i][2]*v[i-1]*v[i-2]+f1*delta1[i]+f2*delta2[i];
	}

	return ans;
};

double randomnumber()
{
	int i,n=15,alpha=3;
	//time_t t;
	double rno=0.0;

	/* Initializes random number generator */
	/* [BUG-1] srand() is called on EVERY invocation. The 29 initialization
	   calls all land inside the same wall-clock second, so time(NULL) is
	   identical, rand() restarts the same sequence, and every shell gets the
	   SAME magnitude -- differing only by the pow(-1,i) sign flip applied by
	   the caller. The initial condition is therefore not random across shells;
	   it is a deterministic alternating-sign constant vector. */
	srand((unsigned) time(NULL));

	for(i=0;i<n;i++)
	{
		rno=rno+(rand()%10)*pow(10,-alpha-i);
	}
	return rno;
};

int main()
{
	time_t t;

	double k1,k2,k3,k4;
	int p,i,j,l1,t1,m;

	//Heaviside theta for Ekman dissipation
	int theta2D[N];
	for(i=0;i<N;i++)
	{
		if (i<=5)
			theta2D[i]=1;

		else
			theta2D[i]=0;
	}

	//Kronecker delta for forcing terms
	int delta1[N],delta2[N];
	for(i=0;i<N;i++)
	{
		delta1[i]=0;
		delta2[i]=0;
	}
	delta1[19]=1;
	delta2[20]=1;

	//wavenumbers
	double k[N],knorm[N];
	for(i=0;i<N;i++)
	{
		k[i]=k0*pow(q,i);
	}

	//coupling coefficients
	double c[N][3];
	for(i=0;i<N;i++)
	{
		c[i][0]=c0*pow(q,i);
	}
	c[0][1]=0;
	c[0][2]=0;
	c[1][2]=0;
	c[N-2][0]=0;
	c[N-1][0]=0;
	c[N-1][1]=0;

	for(i=1;i<N-1;i++)
	{
		c[i][1]=-c[i-1][0]*(pow(k[i-1],2)-pow(k[i+1],2))/(pow(k[i],2)-pow(k[i+1],2));
	}

	for(i=2;i<N;i++)
	{
		c[i][2]=c[i-2][0]*(pow(k[i-2],2)-pow(k[i-1],2))/(pow(k[i-1],2)-pow(k[i],2));
	}
	c[0][1]=0;
	c[0][2]=0;
	c[1][2]=0;
	c[N-2][0]=0;
	c[N-1][0]=0;
	c[N-1][1]=0;

	//u is initialized as random numbers of the order of 1e-3
	double u[N][2];
	for(i=0;i<N;i++)
	{
		u[i][0]=pow(-1,i)*randomnumber();
	}


	//Energy accumulators
	/* [BUG-2] Enorm is sized by `i`, a loop variable, not by N. It is an
	   accidental VLA that happens to have length 29 only because the preceding
	   loop exited with i==N. Any reordering of the blocks above silently
	   corrupts memory. */
	double E[N],Enorm[i],E1[N];
	for(i=0;i<N;i++)
	{
		E[i]=0.0;
		E1[i]=0.0;
	}


	double K[N][4],v[N],a[N];
	for (i=0; i<N; i++)
	{
		for (j=0;j<4;j++)
		{
			K[i][j]=0;
		}
	}

	int r,itno;
	int itmax=sizeof(iter)/sizeof(int);

	for (itno=0;itno<itmax;itno++)
	{
		if(itno==0)
		{
			m=iter[itno];
		}
		else
		{
			m=iter[itno]-iter[itno-1];
		}

		//time-marching starts
		for(t1=0;t1<m-1;t1++)
		{
			//K1
			for (i=0;i<N;i++)
			{
				v[i]=0;
				a[i]=0;
			}
			for(i=0;i<N;i++)
			{
				v[i]=u[i][0]+a[i];
			}
			for(r=0;r<N;r++)
			{
				K[r][0]=dt*f2D(0,r,v,k,c,delta1,delta2,theta2D);
			}

			//K2
			for (i=0;i<N;i++)
			{
				a[i]=K[i][0]/2;
			}
			for(i=0;i<N;i++)
			{
				v[i]=u[i][0]+a[i];
			}
			for(r=0;r<N;r++)
			{
				K[r][1]=dt*f2D(0,r,v,k,c,delta1,delta2,theta2D);
			}

			//K3
			for (i=0; i<N; i++)
			{
				a[i]=K[i][1]/2;
			}
			for(i=0;i<N;i++)
			{
				v[i]=u[i][0]+a[i];
			}
			for(r=0;r<N;r++)
			{
				K[r][2]=dt*f2D(0,r,v,k,c,delta1,delta2,theta2D);
			}

			//K4
			for (i=0; i<N; i++)
			{
				a[i]=K[i][2];
			}
			for(i=0;i<N;i++)
			{
				v[i]=u[i][0]+a[i];
			}
			for(r=0;r<N;r++)
			{
				K[r][3]=dt*f2D(0,r,v,k,c,delta1,delta2,theta2D);
			}

			for(i=0;i<N;i++)
			{
				u[i][1]=u[i][0]+(K[i][0]+2*(K[i][1]+K[i][2])+K[i][3])/6;
			}

			/* [BUG-4, performance/precision] naive accumulation of a small
			   increment into a large running total over up to 4.8e8 steps. */
			for(i=0;i<N;i++)
			{
				E1[i]=E1[i]+pow(u[i][1],2);
			}

			for(i=0;i<N;i++)
			{
				u[i][0]=u[i][1];
			}
		} /*time marching end*/

		//average E
		for(i=0;i<N;i++)
		{
			E[i]=E1[i]/(2*(iter[itno]-1)*k[i]);
		}

		//Kolmogorov scale and other derived parameters
		double eta=0.0,kd,eps=0.0,kf;
		for(i=0;i<N;i++)
		{
			eta=eta+2*nu*E[i]*pow(k[i],4);
		}
		kd=pow(eta,0.166666)*pow(nu,-0.5);

		//normalized quantities
		for(i=0;i<N;i++)
		{
			knorm[i]=k[i]/kd;
			/* [BUG-3] -1/6 is INTEGER division -> 0, so pow(eta,0)==1 and the
			   eta^(-1/6) factor of Eq. 2.11 silently vanishes. Slopes (BKL,
			   Kolmogorov) are unaffected; the vertical placement of each
			   spectrum is not comparable across times. Note 0.166666 was
			   written correctly two lines above in kd. */
			Enorm[i]=E[i]*pow(nu,-1.5)*pow(eta,-1/6);
		}

		for(i=0;i<N;i++)
		{
			printf("\n%lf\t%lf",log10(knorm[i]),log10(Enorm[i]));
			//printf("\n%lf",Enorm[i]);
		}
		printf("\n\n\n");
		printf("\n\nDissipation wave number is: %lf\n\n",kd);

		//writing to a file
		FILE *f = fopen("energyspectrum", "a");
		fprintf(f,"\nTime=%lf seconds\n",iter[itno]*dt);
		for(i=0;i<N;i++)
		{
			fprintf(f,"\n%lf\t%lf",log10(knorm[i]),log10(Enorm[i]));
		}
		fprintf(f,"\n");
		fprintf(f,"\n%lf\t%lf",kd,eta);
		fprintf(f,"\n\n\n");
		fclose(f);
	}
	return 0;
}
