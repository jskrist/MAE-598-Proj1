#include <stdio.h>
#include <math.h>

double update(int &x,int &y,int &z);

//**************************************************
//  Global Variables
//**************************************************

// Physical Dimensions of the box and default Temp
const double w = 1,
	   		 h = 1,
	    	 l = 1,
			 defltTemp = 0;

// Number of nodes in each direction
const int i = (int)(w*50),
	  	  j = (int)(h*50),
		  k = (int)(l*50);

// Array of all the nodes
double nodes[k][j][i] = { defltTemp };

// Thermal defusivity and time step
double alpha=1, dt=0.00001;

//**************************************************
//  Main Program
//**************************************************
int main(int argc, char** argv)
{
	// Open file to which data will be written
	FILE *Ofile;
	Ofile = fopen("tempData.txt","w");

	// coordinate of the middle of the structure
	int midI = (int)floor(double(i)/2.0);
	int midJ = (int)floor(double(j)/2.0);
	int midK = (int)floor(double(k)/2.0);

//	printf("midI = %i, midJ = %i, midK = %i", midI, midJ, midK);

	int cnt = 0;
	double qo = 1000.0;
	double qf = 0.0;

	double lastTemp = defltTemp;
	double diff = 0,
		   minDiff = 1E-6,
		   minChg = 1E-2;

	// Flags
	bool changed = false,
		 done	 = false;

	// Set initial temperature
	for(int z = 0; z < k; z+=(k-1))
	{
		for(int y = 0; y < j; y++)
		{
			for(int x = 0; x < i; x++)
			{
				if(z == 0)
					nodes[z][y][x] = qo;
				else if (z == k-1)
					nodes[z][y][x] = qf;
			}
		}
	}

	while(!done)
	{
		cnt++;
		for(int z = 1; z < k-1; z++)
		{
			for(int y = 0; y < j; y++)
			{
				for(int x = 0; x < i; x++)
				{
					nodes[z][y][x] = update(x,y,z);

					if((nodes[k-2][midJ][midI] > defltTemp + minChg  ||
						nodes[k-2][midJ][midI] < defltTemp - minChg) && !changed)
					{
						printf("midNode = %f\n", nodes[midK][midJ][midI]);
						changed = true;
					}
					if(x == midI && y == midJ)
						fprintf(Ofile,"%d,%f\n", cnt, nodes[z][y][x]);
				}
			}
		}
		if (changed)
		{
			diff = fabs(lastTemp - nodes[k-2][midJ][midI]);
			lastTemp = nodes[k-2][midJ][midI];
			if(diff < minDiff)
			{
				printf("diff = %f, minDiff = %f\n", diff, minDiff);
				done = true;
			}
		}
	}
	fclose(Ofile);

	return 0;
}

double update(int &x, int &y, int &z)
{
	double result = 0;

	double xP1=0,
		   xM1=0,
		   yP1=0,
		   yM1=0,
		   zP1=0,
		   zM1=0,
		   nCur=nodes[z][y][x];

	if(x+1 < i)
		xP1=nodes[z][y][x+1];
	else
		xP1=0;

	if(x-1 >= 0)
		xM1=nodes[z][y][x-1];
	else
		xM1=0;

	if(y+1 < j)
		yP1=nodes[z][y+1][x];
	else
		yP1=0;

	if(y-1 >= 0)
		yM1=nodes[z][y-1][x];
	else
		yM1=0;

	zP1=nodes[z+1][y][x];
	zM1=nodes[z-1][y][x];

	result = alpha*dt*(pow((double)i,2)*(xM1-2*nCur+xP1)
					 + pow((double)j,2)*(yM1-2*nCur+yP1)
					 + pow((double)k,2)*(zM1-2*nCur+zP1)) + nCur;

	return result;
}
