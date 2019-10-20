#include "stdlib.h"
#include "stdio.h"
#include "mpi.h"
#include "time.h"
#include "math.h"

double f(double a)
{
    return a*a;
}

void zrhqr(float a[], int m, float rtr[], float rti[])
{
    void balanc(float **a, int n);
    void hqr(float **a, int n, float wr[], float wi[]);
    int j,k;
    float **hess,xr,xi;

    hess=matrix(1,MAXM,1,MAXM);
    if (m > MAXM || a[m] == 0.0) nrerror("bad args in zrhqr");
    for (k=1;k<=m;k++) {
        hess[1][k] = -a[m-k]/a[m];
        for (j=2;j<=m;j++) hess[j][k]=0.0;
        if (k != m) hess[k+1][k]=1.0;
    }
    balanc(hess,m);
    hqr(hess,m,rtr,rti);
    for (j=2;j<=m;j++) {
        xr=rtr[j];
        xi=rti[j];
        for (k=j-1;k>=1;k--) {
            if (rtr[k] <= xr) break;
            rtr[k+1]=rtr[k];
            rti[k+1]=rti[k];
        }
        rtr[k+1]=xr;
        rti[k+1]=xi;
    }
    free_matrix(hess,1,MAXM,1,MAXM);
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
    

    //-------------------------------------------------------------------------
    // simple rectangle. only outer loop will be paralleled
    //double rect_time_start = MPI_Wtime();
    int rect_size = 2;
    double total_rect_int_val = 0, old_int_val = 0;
    double a_start = a_min + (id * a_range) / procnum;
    double a_step = a_range / (rect_size * procnum);
    do{
        MPI_Barrier(MPI_COMM_WORLD);
        old_int_val = total_rect_int_val;
        double rect_int_val = 0.;
        double a = a_start;
        for(int l = 0; l < rect_size; l++, a += a_step)
        {
            rect_int_val += f(a);
        }
        total_rect_int_val = 0;
        MPI_Allreduce(&rect_int_val, &total_rect_int_val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        total_rect_int_val *= a_range / (rect_size * procnum);
        rect_size *= 2;
        a_step /= 2;
        err = fabs(total_rect_int_val - old_int_val)/3;
    } while (err > eps);

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
