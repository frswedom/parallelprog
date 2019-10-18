#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <stdio.h>

double func(double x)
{
    return sin(x);
}

double integral(double start, int n, double step)
{
    double x = start, value = 0.0;
    value = func(x) / 2.;
    
    for (int j = 0; j < n - 1; j++)
    {
        x += step;
        value += func(x);
    }
    value += func(x + step) / 2.;
    return value * step;
}

double integral_SR(double start, int n, double step) //Simpson integration rule
{
    double x = start, value = 0.0;
    value = func(x) / 2.;
    int N = n / 3;

    for (int j = 0; j < n - 1; j++)
    {
        x += step;
        value += func(x);
    }
    value += func(x + step) / 2.;
    return value * step;
}

void RungeKuttaMethod( double (*f)(double, double), double minValue, double maxValue, double initialValue)
{
    const int steps = 20;
    
    double y[steps + 1];
    double x[steps + 1];
    
    double h = (maxValue - minValue) / steps;
    x[0] = minValue;
    x[steps] = maxValue;
    y[0] = initialValue;
    for (int i = 0; i < steps; i++)
    {
        x[i] = x[0] + i*h;
        double k1 = h*f(x[i], y[i]);
        double k2 = h*f(x[i] + h / 2, y[i] + k1 / 2);
        double k3 = h*f(x[i] + h / 2, y[i] + k2 / 2);
        double k4 = h*f(x[i] + h, y[i] + k3);

        y[i + 1] = y[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }
    cout << "RungeKuttaMethod\n";
    for (int i = 0; i < steps; i++)
    {
        cout << x[i] << ", ";
    }
    cout << x[steps];
    cout << '\n';
    for (int i = 0; i < steps; i++)
    {
        cout << y[i] << ", ";
    }
    cout << y[steps];
    cout << '\n';
    return;
}

int main(int argc, char **argv)
{   
	int process_num, my_id;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &process_num);

    int total_iterations, my_num;
    double a, b, step, result, my_a, range, starttime;    
    
    total_iterations = 1 << 11; 
    a = 0.;
    b = M_PI;
    step = (b - a) / total_iteratoins;
    result = 0.;
    
    if(my_id == 0)
        starttime = MPI_Wtime();
    
    my_num = total_iterations / process_num;
    range = (b - a) / process_num;
    
    my_a = a + my_id * range;
    
    double time_thread = MPI_Wtime();
    double my_result = integral(my_a, my_num, step);
    time_thread = (MPI_Wtime() - time_thread) * 1000;

    printf("Process %d. Time used %f ms. My result = %f\n", my_id, time_thread, my_result);

    MPI_Reduce(&my_result, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (my_id == 0)
    {
        printf("Integral = %f\n", result);
        time_total = ( MPI_Wtime() - starttime ) * 1000.;
        printf("Total time is%f ms \n",time_total);

    }
    MPI_Finalize();
    return 0;
}
