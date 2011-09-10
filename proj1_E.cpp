#include <stdio.h>
#include <math.h>

double update(int &x,int &y,int &z);

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

// Array of all the nodes
double nodes[k][j][i] = { defltTemp };

// Thermal defusivity and time step
double alpha=1, dt=0.0001;

//**************************************************
//  Main Program
//**************************************************
int main(int argc, char** argv)
{
	// Open file to which data will be written
	FILE *Ofile;
	Ofile = fopen("Data.txt","w");

	// coordinate of the middle of the structure
	int midI = (int)floor(double(i)/2.0);
	int midJ = (int)floor(double(j)/2.0);
	int midK = (int)floor(double(k)/2.0);

//	printf("midI = %i, midJ = %i, midK = %i", midI, midJ, midK);

	int cnt = 0;
	double qo = 100.0;
	double qf = 0.0;

	double	lastTemp = defltTemp,
			diff = 0,
			minDiff = 1E-4,
			minChg = 1E-2,
			curNodeTemp = 0;

	// Flags
	bool changed = false,
		 done	 = false,
		 endLoop = false;

	// Set initial temperature
	for(int z = 0; z < k; z+=(k-1))
	{
		for(int y = 0; y < j; y++)
		{
			for(int x = 0; x < i; x++)
			{
				if(z == 0)
					nodes[z][y][x] = qo;
				else
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
					curNodeTemp = update(x,y,z);
					nodes[z][y][x] = curNodeTemp;

					if(curNodeTemp > (qo + qf))
					{
						printf("Became unstable");

						return 1;
					}
					else if(curNodeTemp <= (defltTemp + minChg*minChg) )
					{
						endLoop = true;
						break;
					}

					if( x == midI && y == midJ )
					{
						if( !changed && z == (k-2) )
						{
							if(curNodeTemp > defltTemp + minChg  || curNodeTemp < defltTemp - minChg)
								changed = true;
						}

						if(cnt %100 == 0)
							fprintf(Ofile,"%d,%f\n", cnt, nodes[z][y][x]);
					}
				}
				if(endLoop)
					break;
			}
			if(endLoop)
				break;
		}
		if(endLoop)
			endLoop = false;

		if (changed)
		{
			diff = fabs(lastTemp - nodes[k-2][midJ][midI]);
			lastTemp = nodes[k-2][midJ][midI];
			if(diff < minDiff)
			{
//				printf("diff = %f, minDiff = %f\n", diff, minDiff);
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
		xP1=nCur;

	if(x-1 >= 0)
		xM1=nodes[z][y][x-1];
	else
		xM1=nCur;

	if(y+1 < j)
		yP1=nodes[z][y+1][x];
	else
		yP1=nCur;

	if(y-1 >= 0)
		yM1=nodes[z][y-1][x];
	else
		yM1=nCur;

	zP1=nodes[z+1][y][x];
	zM1=nodes[z-1][y][x];

	result = alpha*dt*(pow((double)dw,2)*(xM1-2*nCur+xP1)
					 + pow((double)dh,2)*(yM1-2*nCur+yP1)
					 + pow((double)dl,2)*(zM1-2*nCur+zP1)) + nCur;

	return result;
}