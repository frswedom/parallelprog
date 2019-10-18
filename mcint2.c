#include "stdlib.h"
#include "stdio.h"
#include "mpi.h"
#include "time.h"
#include "math.h"

double f(double a)
{
    return a*a;
}

int main(int argc, char const *argv[])
{
    int id, procnum;

    double a_min = -1., a_max = 1.;
    double a_range = a_max - a_min;
    double eps = 1e-4;
    double err;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &procnum);

    double total_time = MPI_Wtime();
    //-------------------------------------------------------------------------
    // naive monte carlo approach
    //double mc_time_start = MPI_Wtime();

    srand(time(NULL) + id);
    long int mc_total_num = 200000, runs = 80;
    double mc_int[runs];

    for(int j = 0; j < runs; j++)
    {
        double mc_int_val = 0;
        mc_int[j] = 0;

        for(long int i = 0; i < mc_total_num; i++)
        {
            double a = ((double)rand() / (double)RAND_MAX) * a_range + a_min;
            mc_int_val += f(a);
        }
        MPI_Reduce(&mc_int_val, &mc_int[j], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    
    total_time  = MPI_Wtime() - total_time;
    MPI_Finalize();
    if (!id)
    {
       double mc_int_mean = 0, mc_int_meansq = 0;
       for(int j= 0; j< runs; j++)
       {    
            mc_int[j] *= a_range / (procnum * mc_total_num); 
            mc_int_mean += mc_int[j];
            mc_int_meansq += mc_int[j] * mc_int[j];       
}
       double stdev = sqrt((mc_int_meansq - mc_int_mean * mc_int_mean / runs) / (runs - 1));
       printf("MC int mean: %.10f\n", mc_int_mean / runs);
       printf("MC stdev   : %.2e\n", stdev);
    }


    if (!id)
    {
        printf("Error: %.2e\n rect_size: %d\n", eps, rect_size);
        double rect_extr = total_rect_int_val + (total_rect_int_val - old_int_val) / 3;
        printf("Rectangle : %.8f Error: %.2e\n", rect_extr, err);
        printf("Time: %.8lf secs %d procs\n", total_time, procnum);
    }
    return 0;
}
