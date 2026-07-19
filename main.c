#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
double function_y_star (double h, double x, double y);
double f_star(double x, double y_IEM, double h);
double f_normal(double x, double y_IEM);
double runge_kutta_method (double x, double y_rk4, double h);
void graphic_functions (void);
void graphic_errors (void);
int main (void)
{
    double h[5], x = 0.0, argument, y_predict, f1, f2, error_EM = 0.0, error_IEM = 0.0, error_RK4 = 0.0;
    double max_error_EM, max_error_IEM, max_error_RK4;
    double *y_EM, *y_IEM, *y_real, *y_rk4;
    double range [2];
    int i, j, steps [5];
    FILE *doc_graphic, *doc_error;
    doc_graphic = fopen ("results/result_graphic.csv", "w");
    if (doc_graphic == NULL)
    {
        printf ("Error");
        return 1;
    }
    doc_error = fopen ("results/result_errors.csv", "w");
    if (doc_error == NULL)
    {
        printf ("Error");
        return 1;
    }
    printf ("Solving: dy/dx = xy\nInitial condition: y(0)=1\n");
    for (i=0;i<5;i++)
    {
        printf ("Give the value of step %d: ", i+1);
        scanf ("%lf", &h[i]);
    }
    for (i=0;i<2;i++)
    {
        printf ("Give the x%d to determine range: ", i+1);
        scanf ("%lf", range+i);
    }
    for (i=0;i<5;i++)
    {
    steps [i] = (range[1] - range[0])/h[i];
    }
    fprintf (doc_graphic,"0;1;1;1;1\n");
    for (j=0;j<5;j++)
    {
        max_error_EM = max_error_IEM = max_error_RK4 = 0.0;
        y_EM = malloc((steps[j]+1)*sizeof(double));
        if (y_EM == NULL)
        {
            printf ("Error");
            return 1;
        }   
        y_IEM = malloc((steps[j]+1)*sizeof(double));
        if (y_IEM == NULL)
        {
            printf ("Error");
            return 1;
        }   
        y_real = malloc((steps[j]+1)*sizeof(double));
        if (y_real == NULL)
        {
            printf ("Error");
            return 1;
        }
        y_rk4 = malloc((steps[j]+1)*sizeof(double));
        if (y_rk4 == NULL)
        {
            printf ("Error");
            return 1;
        }
        y_real [0] = y_IEM[0] = y_EM[0] = y_rk4 [0] = 1;
        x = range [0];
        for (i=0;i<steps[j];i++)
        {
            y_EM[i+1] = function_y_star (h[j], x, y_EM[i]);
            y_predict = function_y_star(h[j], x, y_IEM[i]);
            f1 = f_star(x, y_predict, h[j]);
            f2 = f_normal(x, y_IEM[i]);
            y_IEM[i+1] = 0.5*h[j]*(f1 + f2) + y_IEM[i];
            y_rk4[i+1] = runge_kutta_method (x, y_rk4[i], h[j]);
            argument = ((x+h[j])*(x+h[j]))/2;
            y_real[i+1] = exp(argument);
            error_EM = fabs(y_real[i+1] - y_EM[i+1]);
            error_IEM = fabs(y_real[i+1] - y_IEM[i+1]);
            error_RK4 = fabs(y_real[i+1] - y_rk4[i+1]);
            if(error_EM > max_error_EM) 
            {
                max_error_EM = error_EM;
            }
            if(error_IEM > max_error_IEM)
            {
                max_error_IEM = error_IEM;
            }
            if(error_RK4 > max_error_RK4)
            {
                max_error_RK4 = error_RK4;
            }
            x = range[0] + (i+1)*h[j];
            if (j==0)
            {
                fprintf (doc_graphic, "%lf;%lf;%lf;%lf;%lf\n", x, y_rk4[i+1], y_IEM[i+1], y_EM[i+1], y_real [i+1]);
            }
        }
        if (j==0)
        {
            fclose (doc_graphic);
            graphic_functions();
            system("gnuplot results/plot_data_functions.gnuplot");
            doc_graphic = fopen ("results/result_graphic.csv", "a");
        }
        fprintf (doc_error, "%lf;%lf;%lf;%lf\n", h[j], max_error_RK4, max_error_IEM, max_error_EM);
        free(y_real);
        free(y_EM);
        free(y_IEM);
        free(y_rk4);
    }
    fclose (doc_graphic);
    fclose (doc_error);
    graphic_errors();
    system("gnuplot results/plot_data_errors.gnuplot");    
    return 0;
}

double function_y_star (double h, double x, double y)
{
    double product;
    product = x*y;
    double result;
    result = h*product + y;
    return result;
}

double f_star(double x, double y_IEM, double h)
{
    return (x+h)*y_IEM;
}

double f_normal(double x, double y_IEM) 
{
    return x*y_IEM;
}

double runge_kutta_method (double x, double y_rk4, double h)
{
    double k1, k2, k3, k4;
    k1 = f_normal (x, y_rk4);
    k2 = f_normal (x+0.5*h, y_rk4+(0.5*h*k1));
    k3 = f_normal (x+0.5*h, y_rk4+(0.5*h*k2));
    k4 = f_normal (x+h, y_rk4+(h*k3));
    return y_rk4 + (1.0/6.0)*h*(k1 + 2*k2 + 2*k3 + k4);
}

void graphic_functions (void)
{
    int success;
    FILE *fp;
    fp = fopen ("results/plot_data_functions.gnuplot","w");
    if (fp == NULL)
    {
        printf("Error");
        exit (1);
    }
    fprintf(fp, "set datafile separator \";\"\n");
    fprintf(fp, "set title \"Graphic\"\n");
    fprintf(fp, "set xlabel \"X\"\n");
    fprintf(fp, "set ylabel \"Y\"\n");
    fprintf(fp, "set grid\n");
    fprintf(fp, "set xrange [0:*]\n");
    fprintf(fp, "set yrange [0:*]\n");
    fprintf(fp, "set terminal png size 1280, 720\n");
    fprintf(fp, "set output \"results/differential eq. solving.png\"\n");
    fprintf(fp, "plot \"results/result_graphic.csv\" using 1:2 title \"RK4\" with lines, \\\n");
    fprintf(fp, "\"results/result_graphic.csv\" using 1:3 title \"Improved Euler\" with lines, \\\n");
    fprintf(fp, "\"results/result_graphic.csv\" using 1:4 title \"Euler\" with lines, \\\n");
    fprintf(fp, "\"results/result_graphic.csv\" using 1:5 title \"real value\" with lines, \n");
    success = fclose (fp);
    if (success != 0)
    {
        printf ("Error");
        exit (1);
    }
    else 
    {
        printf ("Done 1\n");
    }
}

void graphic_errors (void)
{
    int success;
    FILE *fp;
    fp = fopen ("results/plot_data_errors.gnuplot","w");
    if (fp == NULL)
    {
        printf("Error");
        exit (1);
    }
    fprintf(fp, "set datafile separator \";\"\n");
    fprintf(fp, "set title \"Errors\"\n");
    fprintf(fp, "set xlabel \"h\"\n");
    fprintf(fp, "set ylabel \"Err\"\n");
    fprintf(fp, "set logscale xy\n");
    fprintf(fp, "set grid\n");
    fprintf(fp, "set terminal png size 1280, 720\n");
    fprintf(fp, "set output \"results/Errors vs h.png\"\n");
    fprintf(fp, "plot \"results/result_errors.csv\" using 1:2 title \"errors RK4\" with linespoints, \\\n");
    fprintf(fp, "\"results/result_errors.csv\" using 1:3 title \"errors IΕΜ\" with linespoints, \\\n");
    fprintf(fp, "\"results/result_errors.csv\" using 1:4 title \"errors ΕΜ\" with linespoints, \n");
    success = fclose (fp);
    if (success != 0)
    {
        printf ("Error");
        exit (1);
    }
    else 
    {
        printf ("Done 2");
    }
}