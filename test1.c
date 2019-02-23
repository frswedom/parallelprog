#include <mpi.h> 
#include <stdio.h> 
int main(int argc, char *argv[])
{ 
	int ProcNum, ProcRank, tmp; 
	MPI_Status status; 
	MPI_Init(&argc, &argv); 
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum); 
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank); 
	printf("Hello world from process %i \n", ProcRank); 
	MPI_Finalize(); 
	return 0; 
} 