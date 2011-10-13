#include <iostream>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string>

double update(int &x,int &y,int &z);
void   setDims(int &nSize);

//**************************************************
//  Global Variables
//**************************************************

// the step size in each direction and default Temp
const double dw = 200, // X
	   		 dh = 200, // Y
	    	 dl = 200, // Z
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
double alpha=1, dt=0.00003;

// Cartesian communicator variables
MPI_Comm oldComm, cartComm;
int ndims = 3,
	dims[3],
	periods[3] = {0, 0, 0},
	reorder = 0,
	coords[3];

//**************************************************
//  Main Program
//**************************************************
int main(int argc, char** argv)
{
	int rank, size;
	int nLoopX = 0, nLoopY = 0, nLoopZ = 0;

	MPI_Init(&argc, &argv);
	oldComm = MPI_COMM_WORLD;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	setDims(size);

	int dX = i / dims[0],
		dY = j / dims[1],
		dZ = k / dims[2];
	
	// Creates a cartesian communicator
	MPI_Cart_create(oldComm, ndims, dims, 
					periods, reorder, &cartComm);

	// Get the coordinates of this processor
	MPI_Cart_coords(cartComm, rank, 3, coords);

	// Set the start points for the loops
	nLoopX = coords[0] * dX;
	nLoopY = coords[1] * dY;
	nLoopZ = coords[2] * dZ;

	// Open file to which data will be written
	FILE *Ofile;
	// Only have the first processor control the file
	if(rank == size - 1)
	{
		Ofile = fopen("Data.txt","w");
	}

	// coordinate of the middle of the structure
	int midI = (int)floor(double(i) * 0.5);
	int midJ = (int)floor(double(j) * 0.5);

	int cnt = 0;
	double qo = 100.0;
	double qf = 0.0;

	double	lastTemp = defltTemp,
			diff = 0,
			minDiff = 2E-5,
			minChg = 1E-2,
			curNodeTemp = 0;

	// Flags
	bool endLoop = false;
	int  done 	 = 0,
		 changed = 0;

	// Set initial temperature
	if(nLoopZ == 0)
	{
		for(int y = nLoopY; y < nLoopY + dY; y++)
		{
			for(int x = nLoopX; x < nLoopX + dX; x++)
			{
				nodes[0][y][x] = qo;
			}
		}
	}
	else if (nLoopZ + dZ == k)
	{
		for(int y = nLoopY; y < nLoopY + dY; y++)
		{
			for(int x = nLoopX; x < nLoopX + dX; x++)
			{
				nodes[k-1][y][x] = qf;
			}
		}
	}

	// Make sure all of the initalization is done
//	MPI_Barrier(cartComm);
	MPI_Allreduce( MPI_IN_PLACE, &nodes, i*j*k, MPI_DOUBLE, MPI_MAX, cartComm );

	while(done == 0)
	{
//std::cout << rank << " is still working\n";
		cnt++;
		for(int z = nLoopZ; z < nLoopZ + dZ; z++)
		{
			if(z != 0 && z != k-1)
			{
				for(int y = nLoopY; y < nLoopY + dY; y++)
				{
					for(int x = nLoopX; x < nLoopX + dX; x++)
					{
						curNodeTemp = update(x,y,z);
						nodes[z][y][x] = curNodeTemp;

						if(curNodeTemp > (qo + qf))
						{
							printf("Became unstable");
							MPI_Finalize();

							return 1;
						}
						else if(curNodeTemp < (defltTemp + minDiff) )
						{
							endLoop = true;
							break;
						}
					}
					if(endLoop)
						break;
				}
				if(endLoop)
					break;
			}
		}
		endLoop = false;

		if( changed == 0 )
		{
			if(nodes[k-2][nLoopY][nLoopX] > defltTemp + minChg || nodes[k-2][nLoopY][nLoopX] < defltTemp - minChg)
			{
				changed = 1;
			}
		}
		//  ASUMPTION that the loop in which it initially changes is not the same 
		//  loop that it reaches the steady-state otherwise we will get one more
		//  itteration
		else
		{
			diff = fabs(lastTemp - nodes[k-2][nLoopY][nLoopX]);
			lastTemp = nodes[k-2][nLoopY][nLoopX];
			if(diff < minDiff)
			{
				done = 1;
			}
		}

		if(cnt%100 == 0 || done > 0)
		{
			// If one processor meets the ending requirements then send that to
			// all processors so they can quit working
			MPI_Allreduce( MPI_IN_PLACE, &done, 1, MPI_INT, MPI_MAX, cartComm);

			// Make sure all data is up to date at the end of every 100th time step
			MPI_Allreduce( MPI_IN_PLACE, &nodes, i*j*k, MPI_DOUBLE, MPI_MAX, cartComm);
			MPI_Allreduce( MPI_IN_PLACE, &changed, 1, MPI_INT, MPI_MAX, cartComm);	
		}
	}
	// Note when a processor exits the loop
	// std::cout << "Node " << rank << " is out\n";
	
	// When the last node exits have it write out to the data file
	if(rank == size - 1)
	{
		for(int z = 0; z < k; z++)
		{
			fprintf(Ofile,"%i,%f,%i\n", cnt, nodes[z][midJ][midI], rank);
		}
		fclose(Ofile);
	}
	// Make sure all the nodes are ready before we finalize and exit
	//MPI_Barrier(oldComm);

	MPI_Finalize();
	return 0;
}

void   setDims(int &nSize)
{
	int tmpDim = (int)floor(pow((double)nSize,0.33333333333333333333));
	int tmpSize = nSize / tmpDim;

	dims[0] = tmpDim;

	tmpDim = (int)floor(pow((double)tmpSize,0.5));

	dims[1] = tmpDim;

	tmpDim = tmpSize / tmpDim;

	dims[2] = tmpDim;	
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
