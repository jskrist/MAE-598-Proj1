#include <iostream>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

double update(int &x,int &y,int &z, double &nodes);
double updateFace(int &x, int &y, int &z, double &nodes, Face &fFace);
void   setDims(int &nSize);
void   setNeighbors();
//**************************************************
//  Global Variables
//**************************************************

// 1 over the step size in each direction 
// and default Temp
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
//double nodes[k][j][i] = { defltTemp };

// Thermal defusivity and time step
double alpha=1, dt=1E-5;

// Cartesian communicator variables
MPI_Comm oldComm, cartComm;
int ndims = 3,
	dims[3],
	periods[3] = {0, 0, 0},
	reorder = 0,
	coords[3];

// processor rank and group size
int rank, size;

// flags to check neighbors existence
bool left	= false,
	 right	= false,
	 top	= false,
	 bottom	= false,
	 front	= false,
	 back	= false;

// rank of neighbors will be stored here
int  nLeft		= -1,
	 nRight		= -1,
	 nTop		= -1,
	 nBottom	= -1,
	 nFront		= -1,
	 nBack		= -1;

// Enumerated type to send to updateFace to know what face to update
enum Face { LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK };

//**************************************************
//  Main Program
//**************************************************
int main(int argc, char** argv)
{
	int nLoopX = 0, nLoopY = 0, nLoopZ = 0;

	MPI_Init(&argc, &argv);
	oldComm = MPI_COMM_WORLD;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	setDims(size);

	int dX = i / dims[0],
		dY = j / dims[1],
		dZ = k / dims[2];
	
	// Dynamically create a 3d array to store values
	double ***core; // = (double***)malloc(( dX-2 )*sizeof(double**));
	core = new double**[dZ-2];
	for (int r = 0; r < dZ-2; ++r)
	{
		core[r] = new double*[dY-2];

		for (int q = 0; q < dY-2; ++q)
 			core[r][q] = new double[dX-2];
	}
	
/*	for (int r = 0; r < (dZ-2); r++)
	{
    	core[r] = (double**)malloc((dZ-2)*sizeof(double*));
    	for (int q = 0; q < (dY-2); q++)
		{
      		core[r][q] = (double*)malloc((dY-2)*sizeof(double));
    	}
  	}
*/
	// initialize core array
	for(int r = 0; r < dZ-2; r++)
	{
		for(int q = 0; q < dY-2; q++)
		{
			for(int p = 0; p < dX-2; p++)
			{
				core[r][q][p] = 0;
			}
		}
	}

	// Creates a cartesian communicator
	MPI_Cart_create(oldComm, ndims, dims, 
					periods, reorder, &cartComm);

	// Get the coordinates of this processor
	MPI_Cart_coords(cartComm, rank, 3, coords);

	// Set the start points for the loops
	nLoopX = coords[0] * dX;
	nLoopY = coords[1] * dY;
	nLoopZ = coords[2] * dZ;

	// Determine what neighbors are around
	setNeihbors();

	// Open file to which data will be written
	FILE *Ofile;

	// Only have the first processor control the file
	if(rank == size - 1)
	{
		Ofile = fopen("Data.txt","w");
	}

	// coordinate of the middle of the structure
	//	int midI = (int)floor(double(i) * 0.5);
	//	int midJ = (int)floor(double(j) * 0.5);

	int cnt = 0;		// loop counter
	double qo = 100.0;	// initial front temp
	double qf = 0.0;	// initial back temp

	double	lastTemp = defltTemp,
			diff = 0,
			minDiff = 2E-5,
			minChg = 1E-2,
			curNodeTemp = 0;

	// Flags
	bool endLoop = false;
	int  done 	 = 0,
		 changed = 0;

	// Dynamically create send data arrays
	double** sendLeft;
	double** sendRight;
	double** sendTop;
	double** sendBottom;
	double** sendFront;
	double** sendBack;

	// Dynamically create recv data array
	double** recvLeft;
	double** recvRight;
	double** recvTop;
	double** recvBottom;
	double** recvFront;
	double** recvBack;

	sendLeft 	= new double*[dZ];
	sendRight	= new double*[dZ];
	for (int r = 0; r < dZ; ++r)
	{
    	sendLeft[r]  = new double[dY];
    	sendRight[r] = new double[dY];
	}

	sendTop		= new double*[dZ];
	sendBottom	= new double*[dZ];
	for (int r = 0; r < dZ; ++r)
	{
    	sendTop[r]    = new double[dX];
    	sendBottom[r] = new double[dX];
	}

	sendFront 	= new double*[dY];
	sendBack	= new double*[dY];
	for (int q = 0; q < dY; ++q)
	{
    	sendFront[q] = new double[dX];
    	sendBack[q]	 = new double[dX];
	}

/*	// Dynamically create send data arrays
	double* sendLeft	= (double*) malloc(dY * dZ * sizeof(double));
	double* sendRight	= (double*) malloc(dY * dZ * sizeof(double));
	double* sendTop		= (double*) malloc(dX * dZ * sizeof(double));
	double* sendBottom	= (double*) malloc(dX * dZ * sizeof(double));
	double* sendFront	= (double*) malloc(dY * dX * sizeof(double));
	double* sendBack	= (double*) malloc(dY * dX * sizeof(double));

	// Dynamically create recv data array
	double* recvLeft	= (double*) malloc(dY * dZ * sizeof(double));
	double* recvRight	= (double*) malloc(dY * dZ * sizeof(double));
	double* recvTop		= (double*) malloc(dX * dZ * sizeof(double));
	double* recvBottom	= (double*) malloc(dX * dZ * sizeof(double));
	double* recvFront	= (double*) malloc(dY * dX * sizeof(double));
	double* recvBack	= (double*) malloc(dY * dX * sizeof(double));
*/

	for(int z = 0; z < dZ; z++)
	{
		for(int y = 0; y < dY; y++)
		{
			sendLeft[z][y]  = 0;
			sendRight[z][y] = 0;
			recvLeft[z][y]  = 0;
			recvRight[z][y] = 0;
		}
	}
	for(int z = 0; z < dZ; z++)
	{
		for(int x = 0; x < dX; x++)
		{
			sendTop[z][x]	 = 0;
			sendBottom[z][x] = 0;
			recvTop[z][x]	 = 0;
			recvBottom[z][x] = 0;
		}
	}
	for(int y = 0; y < dY; y++)
	{
		for(int x = 0; x < dX; x++)
		{
			sendFront[y][x]	= 0;
			sendBack[y][x]	= 0;
			recvFront[y][x]	= 0;
			recvBack[y][x]	= 0;
		}
	}

	// Set initial temperature
	if(nLoopZ == 0)
	{
		for(int y = 0; y < dY; y++)
		{
			for(int x = 0; x < dX; x++)
			{
				sendFront[y][x] = qo;
			}
		}
	}
	else if (nLoopZ + dZ == k)
	{
		for(int y = 0; y < dY; y++)
		{
			for(int x = 0; x < dX; x++)
			{
				sendBack[y][x] = qf;
			}
		}
	}

	// Make sure all of the initalization is done
	// MPI_Allreduce( MPI_IN_PLACE, &nodes, i*j*k, MPI_DOUBLE, MPI_MAX, cartComm );

	while(done == 0)
	{
//std::cout << rank << " is still working\n";
		cnt++;
		for(int z = 0; z < dZ-2; z++)
		{
			for(int y = 0; y < dY-2; y++)
			{
				for(int x = 0; x < dX-2; x++)
				{
					// Update core 3D array
					curNodeTemp = updateCore(x, y, z, core);
					// Update each face
					updateFace(x, y, z, sendLeft, LEFT);
					updateFace(x, y, z, sendRight, RIGHT);
					updateFace(x, y, z, sendTop, TOP);
					updateFace(x, y, z, sendBottom, BOTTOM);
					// Only update front and back if they are not
					// the very front or very back
					if(nLoopZ != 0)
						updateFace(x, y, z, sendFront, FRONT);
					if(nLoopZ + dZ != k)
						updateFace(x, y, z, sendBack, BACK);

					core[z][y][x] = curNodeTemp;

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
		endLoop = false;

		// Update each face
		for(int z = 0; z < dZ; z++)
		{
			for(int y = 0; y < dY; y++)
			{
				updateFace(0, y, z, sendLeft, LEFT);
				updateFace(dX, y, z, sendRight, RIGHT);
			}
		}
		for(int z = 0; z < dZ; z++)
		{
			for(int x = 0; x < dX; x++)
			{
				updateFace(x, dY, z, sendTop, TOP);
				updateFace(x, 0, z, sendBottom, BOTTOM);
			}
		}
		for(int y = 0; y < dY; y++)
		{
			for(int x = 0; x < dX; x++)
			{
				// Only update front and back if they are not
				// the very front or very back
				if(nLoopZ != 0)
					updateFace(x, y, 0, sendFront, FRONT);
				if(nLoopZ + dZ != k)
					updateFace(x, y, dZ, sendBack, BACK);
			}
		}

		if(nLoopZ + dZ == k)
		{
			if( changed == 0 )
			{
				if(sendBack[0][0] > defltTemp + minChg || sendBack[0][0] < defltTemp - minChg)
				{
					changed = 1;
				}
			}
			//  ASUMPTION that the loop in which it initially changes is not the same 
			//  loop that it reaches the steady-state otherwise we will get one more
			//  itteration
			else
			{
				diff = fabs(lastTemp - sendBack[0][0]);
				lastTemp = sendBack[0][0];
				if(diff < minDiff)
				{
					done = 1;
				}
			}
		}

		if(cnt%100 == 0 || done > 0)
		{
			// If one processor meets the ending requirements then send that to
			// all processors so they can quit working
			MPI_Allreduce( MPI_IN_PLACE, &done, 1, MPI_INT, MPI_MAX, cartComm);
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

/*	for(int r = 0; r < dZ-2; r++)
	{
		for(int q = 0; q < dY-2; q++)
		{
			for(int p = 0; p < dX-2; p++)
			{
				free(core[r][q][p]);
			}
		}
	}
	free(core);
*/
	// Release the allocated Memory
	for (int q = 0; q < dY; ++q)
	{
		delete [] sendLeft[r][q];
		delete [] sendRight[r][q];
	}

	for (int p = 0; p < dX; ++p)
	{
		delete [] sendTop[r][p];
		delete [] sendBottom[r][p];
	}

	for (int p = 0; p < dX; ++p)
	{
		delete [] sendFront[q][p];
		delete [] sendBack[q][p];
	}

	for (int q = 0; q < dY; ++q)
	{
		for (int p = 0; p < dX; ++p)
		{
			delete [] core[q][p];
		}

		delete [] core[q];
	}

	delete [] core;
	delete [] sendLeft;
	delete [] sendRight;
	delete [] sendTop;
	delete [] sendBottom;
	delete [] sendFront;
	delete [] sendBottom;

	MPI_Finalize();
	return 0;
}

void setNeighbors()
{
	// left
	MPI_Cart_shift(cartComm, 0, -1, rank, nLeft);
	if(nLeft >=0 && nLeft <=size)
		left = true;

	// right
	MPI_Cart_shift(cartComm, 0, 1, rank, nRight);
	if(nRight >=0 && nRight <=size)
		right = true;

	// top
	MPI_Cart_shift(cartComm, 1, 1, rank, nTop);
	if(nTop >=0 && nTop <=size)
		top = true;

	// bottom
	MPI_Cart_shift(cartComm, 1, -1, rank, nBottom);
	if(nBottom >=0 && nBottom <=size)
		bottom = true;

	// front
	MPI_Cart_shift(cartComm, 2, -1, rank, nFront);
	if(nFront >=0 && nFront <=size)
		front = true;

	// back
	MPI_Cart_shift(cartComm, 2, 1, rank, nBack);
	if(nBack >=0 && nBack <=size)
		back = true;
}

void setDims(int &nSize)
{
	int tmpDim = (int)floor(pow((double)nSize,0.33333333333333333333));
	int tmpSize = nSize / tmpDim;

	dims[0] = tmpDim;

	tmpDim = (int)floor(pow((double)tmpSize,0.5));

	dims[1] = tmpDim;

	tmpDim = tmpSize / tmpDim;

	dims[2] = tmpDim;	
}

double update(int &x, int &y, int &z, double &nodes)
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

double updateFace(int &x, int &y, int &z, double &nodes, Face &fFace);
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
