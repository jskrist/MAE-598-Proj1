#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "MatToolBox.h"

bool updateRow(int &x,int &y);
void setDims(int &nSize);

//**************************************************
//  Global Variables
//**************************************************

// the step size in each direction and default Temp
const double dw = 100, // X
	   		 dh = 50, // Y
	    	 dl = 50, // Z
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

CVector<double> db(k-2);
CVector<double> dx(k-2);

// Array of all the nodes
double nodes[k][j][i] = { defltTemp };

// Thermal defusivity and time step
double alpha=1, dt=0.00003;

double qo = 100.0,
	   qf = 0.0,
	   minChg = 1E-2;

// Cartesian communicator variables
MPI_Comm oldComm, cartComm;
int ndims = 3,
	dims[3],
	periods[3] = {0, 0, 0},
	reorder = 0,
	coords[3];

int cnt = 0;

//**************************************************
//  Main Program
//**************************************************
int main(int argc, char** argv)
{
	int rank, size;
	int nLoopX = 0, nLoopY = 0, nLoopZ = 0;
	int bufInt[4] = {0,0,0,0};
	double bufDouble[1] = {0.0};
	int fileSize = 0;

//	MPI_Init(&argc, &argv);
	MPI::Init();

	oldComm = MPI_COMM_WORLD;
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
	rank = MPI::COMM_WORLD.Get_rank();
    size = MPI::COMM_WORLD.Get_size();

	setDims(size);

	int dX = i / dims[0],
		dY = j / dims[1],
		dZ = k;
	
	// Creates a cartesian communicator
	MPI_Cart_create(oldComm, ndims, dims, 
					periods, reorder, &cartComm);

	// Get the coordinates of this processor
	MPI_Cart_coords(cartComm, rank, 3, coords);

	// Set the start points for the loops
	nLoopX = coords[0] * dX;
	nLoopY = coords[1] * dY;
	nLoopZ = coords[2] * dZ;
	
	MPI::File dataFile;

	dataFile = MPI::File::Open( oldComm,
								"/scratch/jskristo/LUx200_2",
								MPI::MODE_CREATE | MPI::MODE_RDWR,
								MPI::INFO_NULL);

	fileSize = dataFile.Get_size();

	// coordinate of the middle of the structure
	int midI = (int)floor(double(dX)/2.0);
	int midJ = (int)floor(double(dY)/2.0);

	double	lastTemp = defltTemp,
			diff = 0,
			minDiff = 2E-5,
			curNodeTemp = 0;

	// Flags
	bool changed = false,
		 endLoop = false;
	int	 done	 = 0;
	
	// Set initial temperature
	if(fileSize == 0)
	{
		for(int y = nLoopY; y < nLoopY + dY; y++)
		{
			for(int x = nLoopX; x < nLoopX + dX; x++)
			{
				nodes[0][y][x]   = qo;
				nodes[k-1][y][x] = qf;
			}
		}
	}
	else
	{
		dataFile.Seek(rank*dX*dY*dZ*(4*sizeof(int)+sizeof(double)), MPI_SEEK_SET);
		for(int z = nLoopZ; z < nLoopZ + dZ; z++)
		{
			for(int y = nLoopY; y < nLoopY + dY; y++)
			{
				for(int x = nLoopX; x < nLoopX + dX; x++)
				{
					dataFile.Read(bufInt, 4, MPI_INT);
					dataFile.Read(bufDouble, 1, MPI_DOUBLE);
					nodes[bufInt[3]][bufInt[2]][bufInt[1]] = bufDouble[0];
				}
			}
		}
	}

	// Make sure all of the initalization is done
	MPI_Allreduce( MPI_IN_PLACE, &nodes, i*j*k, MPI_DOUBLE, MPI_MAX, cartComm );

	while(done == 0)
	{
		cnt++;

		if(!updateRow(nLoopX, nLoopY))
		{
			printf("%i has instability in Ax = b\n", rank);
			return 1;
		}
		if( !changed )
		{
			if(nodes[k-2][nLoopY][nLoopX] > defltTemp + minChg  || nodes[k-2][nLoopY][nLoopX] < defltTemp - minChg)
			{
				changed = true;
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
		// If one processor meets the ending requirements then send that to
		// all processors so they can quit working
		MPI_Allreduce( MPI_IN_PLACE, &done, 1, MPI_INT, MPI_MAX, cartComm);
		if(cnt%500 == 0 || done > 0)
		{
			dataFile.Seek(rank*dX*dY*dZ*(4*sizeof(int)+sizeof(double)), MPI_SEEK_SET);
			for(int z = nLoopZ; z < nLoopZ + dZ; z++)
			{
				for(int y = nLoopY; y < nLoopY + dY; y++)
				{
					for(int x = nLoopX; x < nLoopX + dX; x++)
					{
						bufInt[0] = cnt;
						bufInt[1] = x;
						bufInt[2] = y;
						bufInt[3] = z;
						bufDouble[0] = nodes[z][y][x];
						dataFile.Write(bufInt, 4, MPI_INT);
						dataFile.Write(bufDouble, 1, MPI_DOUBLE);
					}
				}
			}
		}
	}

//printf("%i has exited the While loop\n", rank);

	for(int z = nLoopZ; z < nLoopZ + dZ - 1; z++)
	{
		for(int y = nLoopY; y < nLoopY + dY; y++)
		{
			for(int x = nLoopX; x < nLoopX + dX; x++)
			{
				nodes[z][y][x] = dx(z);
			}
		}
	}

	// Make sure all data is up to date at the end
	MPI_Allreduce( MPI_IN_PLACE, &nodes, i*j*k, MPI_DOUBLE, MPI_MAX, cartComm);

	dataFile.Seek(rank*dX*dY*dZ*(4*sizeof(int)+sizeof(double)), MPI_SEEK_SET);
	for(int z = nLoopZ; z < nLoopZ + dZ; z++)
	{
		for(int y = nLoopY; y < nLoopY + dY; y++)
		{
			for(int x = nLoopX; x < nLoopX + dX; x++)
			{
				bufInt[0] = cnt;
				bufInt[1] = x;
				bufInt[2] = y;
				bufInt[3] = z;
				bufDouble[0] = nodes[z][y][x];
				dataFile.Write(bufInt, 4, MPI_INT);
				dataFile.Write(bufDouble, 1, MPI_DOUBLE);
			}
		}
	}
	dataFile.Close();

	MPI::Finalize();

	return 0;
}

void   setDims(int &nSize)
{
	int tmpDim = (int)floor(pow((double)nSize,0.5));
	int tmpSize = nSize / tmpDim;

	dims[0] = tmpDim;

	dims[1] = tmpSize;

	dims[2] = 1;	
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

	for(int q = 1; q < k-1; q++)
	{
		nodes[q][y][x] = dx(q);
	}
	return true;
}

