#include <stdio.h>
#include <math.h>
#include "MatToolBox.h"
#include <omp.h>

bool updateRow(int &x,int &y);

//**************************************************
//  Global Variables
//**************************************************

// the step size in each direction and default Temp
const double dw = 100, // X
	   		 dh = 100, // Y
	    	 dl = 100, // Z
			 defltTemp = 0;

// Physical Dimensions of the box
const double w = 2, // X
	   		 h = 2, // Y
	    	 l = 2; // Z

// Number of nodes in each direction
const int i = (int)(w*dw),
	  	  j = (int)(h*dh),
		  k = (int)(l*dl);

CMatToolBox<double> M;

CMatrix<double> dA(k-2,k-2);


double qo = 100.0,
	   qf = 0.0,
	   minChg = 1E-2;

CVector<double> db(k-2);
CVector<double> dx(k-2);

// Array of all the nodes
double nodes[k][j][i] = { defltTemp };

// Thermal defusivity and time step
double alpha=1, dt=0.00003;


int cnt = 0;

//**************************************************
//  Main Program
//**************************************************
int main(int argc, char** argv)
{
	// Open file to which data will be written
//	FILE *Ofile;
//	Ofile = fopen("DataC.txt","w");

	// coordinate of the middle of the structure
	int midI = (int)floor(double(i)/2.0);
	int midJ = (int)floor(double(j)/2.0);

	double	lastTemp = defltTemp,
			diff = 0,
			minDiff = 2E-5,
			curNodeTemp = 0;

	// Flags
	bool changed = false,
		 done	 = false,
		 endLoop = false;

	// Set initial temperature
	#pragma omp parallel for
	for(int z = 0; z < k; z+=(k-1))
	{
		for(int y = 0; y < j; y++)
		{
			for(int x = 0; x < i; x++)
			{
				if(z == 0)
				{
					nodes[z][y][x] = qo;
				}
				else
				{
					nodes[z][y][x] = qf;
				}
			}
		}
	}

	int temp = 0;
	while(!done)
	{
		cnt++;

		if(!updateRow(temp, temp))
		{
			printf("instability in Ax = b");
			return 1;
		}

		if( !changed )
		{
			if(nodes[k-2][0][0] > defltTemp + minChg  || nodes[k-2][0][0] < defltTemp - minChg)
			{
				changed = true;
//				printf("Changed with node value of %f\n",nodes[k-2][0][0]);
			}
		}
		//  ASUMPTION that the loop in which it initially changes is not the same 
		//  loop that it reaches the steady-state otherwise we will get one more
		//  itteration
		else
		{
			diff = fabs(lastTemp - nodes[k-2][0][0]);
			lastTemp = nodes[k-2][0][0];
			if(diff < minDiff)
			{
				done = true;
//				printf("diff is %f", diff);
			}
		}
	}

	#pragma omp parallel for
	for(int z = 1; z < k-1; z++)
	{
		for(int y = 0; y < j; y++)
		{
			for(int x = 0; x < i; x++)
			{
				nodes[z][y][x] = dx(z + 1);
//				if(x == 0 && y == 0)
//					fprintf(Ofile,"%i, %f\n",cnt, dx(z));
			}
		}
	}

//	fclose(Ofile);

	return 0;
}

bool updateRow(int &x, int &y)
{
	double zP1  = 0,
		   zM1  = 0,
		   ztM1 = 0,
		   ztCur = 0,
		   ztP1 = 0,
		   zCur = nodes[0][y][x],
		   r = 0;

	r = (dt*alpha*dl*dl*0.5);

	#pragma omp parallel for
	for(int q = 1; q < k-1; q++)
	{
		zP1  = nodes[q+1][y][x];
		zCur = nodes[q][y][x];
		zM1  = nodes[q-1][y][x];

		if(q == 1)
		{
			ztM1 = zM1;
			db(q) = (1-2*r)*zCur+r*(zM1+zP1+ztM1);
			dA(q,q)   = (1+2*r);
			dA(q+1,q) = (-1*r);
		}
		else if(q == k-2)
		{
			ztP1 = zP1;
			db(q) = (1-2*r)*zCur+r*(zM1+zP1+ztP1);
			dA((q-1),q) = (-1*r);
			dA((q),q)   = (1+2*r);
		}
		else
		{
			dA(q-1,q) = -r;
			dA(q,q)   = (1+2*r);
			dA(q+1,q) = -r;
			db(q) = (1-2*r)*zCur + r*(zM1+zP1);
		}
	}

//	M.Display( "Matrix dA: ", dA); 

	if(!M.LDLTFactorization (dA, 1e-9))//M.AxEqb(dA, dx, db, 1e-15))
	{
		printf("Failed Factorizing\n");
		return false;
	}

	if(!M.LDLTSolve(dA, dx, db))//M.AxEqb(dA, dx, db, 1e-15))
	{
		printf("Failed Solving Ax = b\n");
		return false;
	}

	#pragma omp parallel for
	for(int q = 1; q < k-1; q++)
	{
		if(dx(q) < defltTemp || dx(q) > qo)
			printf("Became Unstable with dx(%i) = %f at cnt %i\n", q, dx(q), cnt);
		nodes[q][y][x] = dx(q);
	}

	return true;
}

