#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <float.h>
#include <algorithm>

#ifdef WIN32
__int64 timerStart = 0;
#else
struct timeval timerStart;
#endif

#include "timer.h"
#include "utils.h"
#include "tecplot.h"

// area is [AxB]x[CxD]
static double A;
static double B;
static double C;
static double D;
static int OX_LEN;
static int OY_LEN;
static int OX_LEN_1;
static int OY_LEN_1;
static int XY_LEN;
static double TAU;
static int TIME_STEP_CNT;
static int JAK_ITER_CNT;
static double HX;
static double HY;
static double R_SQ; // radius of circle in second power
static double INN_DENSITY; // density inside circle
static double OUT_DENSITY; // density out of circle boundary

#define EPS 10e-7

inline static double analytical_solution_circle(double x, double y)
{
    double x0 = A + OX_LEN*HX/2.;
    double y0 = C + OY_LEN*HY/2.;
    double value = (x-x0)*(x-x0) + (y-y0)*(y-y0);
    if (value <= R_SQ) return INN_DENSITY;
    return OUT_DENSITY;
}

inline static double func_u(double time_value, double x, double y) { return 0; }

inline static double func_v(double time_value, double x, double y) { return 0; }

static double get_phi(int ii, int jj, double* prev_density, double time_value)
{
    double x1 = A + ii*HX - HX/2.;
    double y1 = C + jj*HY - HY/2.;
    double x2 = A + ii*HX + HX/2.;
    double y2 = C + jj*HY - HY/2.;
    double x3 = A + ii*HX + HX/2.;
    double y3 = C + jj*HY + HY/2.;
    double x4 = A + ii*HX - HX/2.;
    double y4 = C + jj*HY + HY/2.;

    double u = func_u(time_value, x1, y1);
    double v = func_v(time_value, x1, y1);
    x1 = x1 - TAU * u;
    y1 = y1 - TAU * v;
    u = func_u(time_value, x2, y2);
    v = func_v(time_value, x2, y2);
    x2 = x2 - TAU * u;
    y2 = y2 - TAU * v;
    u = func_u(time_value, x3, y3);
    v = func_v(time_value, x3, y3);
    x3 = x3 - TAU * u;
    y3 = y3 - TAU * v;
    u = func_u(time_value, x4, y4);
    v = func_v(time_value, x4, y4);
    x4 = x4 - TAU * u;
    y4 = y4 - TAU * v;

    int nx = 10;
    int ny = 10;
    double x_step = 1./nx;
    double y_step = 1./ny;

    // get right part for jakoby
    double phi = 0.;
    double mes = x_step * y_step;
    for(int i = 0; i < nx; ++i)
    {
        for(int j = 0; j < ny; ++j)
        {
            double ideal_x = i*x_step + x_step/2.;
            double ideal_y = j*y_step + y_step/2.;

            double real_x = x1 + (x2 - x1)*ideal_x + (x4 - x1)*ideal_y
                            + (x1 + x3 - x2 - x4)*ideal_x*ideal_y;
            double real_y = y1 + (y2 - y1)*ideal_x + (y4 - y1)*ideal_y
                            + (y1 + y3 - y2 - y4)*ideal_x*ideal_y;

            // find out in which square real point was placed
            int sq_i = (int)((real_x - A) / HX);
            int sq_j = (int)((real_y - C) / HY);
            double x = A + sq_i * HX;
            double y = C + sq_j * HY;

            // formula 4
            double dens = prev_density[sq_i * OY_LEN_1 + sq_j]*(1 - (real_x - x)/HX)*(1 - (real_y - y)/HY)
                          + prev_density[(sq_i + 1) * OY_LEN_1 + sq_j]*((real_x - x)/HX)*(1 - (real_y - y)/HY)
                          + prev_density[(sq_i + 1) * OY_LEN_1 + sq_j + 1]*((real_x - x)/HX)*((real_y - y)/HY)
                          + prev_density[sq_i * OY_LEN_1 + sq_j + 1]*(1 - (real_x - x)/HX)*((real_y - y)/HY);

            double a11 = (x2 - x1) + (x1 + x3 - x2 - x4) * ideal_y;
            double a12 = (x4 - x1) + (x1 + x3 - x2 - x4) * ideal_x;
            double a21 = (y2 - y1) + (y1 + y3 - y2 - y4) * ideal_y;
            double a22 = (y4 - y1) + (y1 + y3 - y2 - y4) * ideal_x;
            double jakob = a11*a22 - a21*a12;

            phi += mes * dens * jakob;
        }
    }

    phi *= 16./9.;
    return phi;
}

static double* solve(double &tme)
{
    StartTimer();

    fflush(stdout);

    double* phi = new double[XY_LEN];
    double* prev_density = new double[XY_LEN];
    double* density = new double[XY_LEN];
    double maxErr = FLT_MIN;

    for (int i = 0; i < OX_LEN_1; i++)
        for (int j = 0; j < OY_LEN_1; j++)
            prev_density[i * OY_LEN_1 + j] = analytical_solution_circle(HX*i, HY*j);

    memcpy(density, prev_density, XY_LEN * sizeof(double));

    for (int tl = 1; tl <= TIME_STEP_CNT; tl++)
    {
        // with usage of prev_density we calculate phi function values
        for (int i = 1; i < OX_LEN; ++i)
            for (int j = 1; j < OY_LEN; ++j)
                phi[i * OY_LEN_1 + j] = get_phi(i, j, prev_density, TAU * tl);

        if (tl == TIME_STEP_CNT)
            print_surface_as_v("phi", OX_LEN, OY_LEN, HX, HY, tl, A, C, phi);

        // fill inner matrix of prev_density by zero
        // because we will use it in Jakoby method
        for (int i = 1; i < OX_LEN; ++i)
            for (int j = 1; j < OY_LEN; ++j)
                prev_density[i * OY_LEN_1 + j] = 0.;

        int ic = 0;
        maxErr = FLT_MAX;
        while(maxErr > EPS && ic < JAK_ITER_CNT)
        {
            for (int i = 1; i < OX_LEN; ++i)
            {
                for (int j = 1; j < OY_LEN; ++j)
                {
                    density[i * OY_LEN_1 + j] = -1./9.*(
                            1.5*(
                                    prev_density[OY_LEN_1 * i + j - 1] + // left
                                    prev_density[OY_LEN_1 * (i - 1) + j] + // upper
                                    prev_density[OY_LEN_1 * i + j + 1] + // right
                                    prev_density[OY_LEN_1 * (i + 1) + j] // bottom
                            ) + 0.25*(
                                    prev_density[OY_LEN_1 * (i + 1) + j + 1] + // bottom right
                                    prev_density[OY_LEN_1 * (i + 1) + j - 1] + // bottom left
                                    prev_density[OY_LEN_1 * (i - 1) + j - 1] + // upper right
                                    prev_density[OY_LEN_1 * (i - 1) + j + 1] // upper left
                            )) + phi[OY_LEN_1 * i + j]/(HX*HY);
                }
            }
            ++ic;
            // print_matrix(density, OX_LEN_1, OY_LEN_1, 4);
            // printf("%s\n", "\n");
            // print_matrix(prev_density, OX_LEN_1, OY_LEN_1, 4);
            maxErr = FLT_MIN;
            for (int i = 1; i < OX_LEN; ++i)
            {
                for (int j = 1; j < OY_LEN; ++j)
                {
                    double val = fabs(density[i * OY_LEN_1 + j] - prev_density[i * OY_LEN_1 + j]);
                    if (val > maxErr) {
                        maxErr = val;
                        //printf("%le\n", val);
                    }
                }
            }
            memcpy(prev_density, density, XY_LEN * sizeof(double));
        }

        // if (tl == 1) {
        // 	printf("EPS: %le max jak = %d\n", EPS, JAK_ITER_CNT);
        // 	printf("Iterations: %d max err.: %le\n", ic, maxErr);
        // 	if (ic >= JAK_ITER_CNT)
        // 	{
        // 		printf("%s\n", "Break on max Jakoby iteration exceed");
        // 	}
        // 	else if (maxErr < EPS)
        // 	{
        // 		printf("%s\n", "Break on l_inf exceed");
        // 	}
        // }
    }

    delete[] prev_density;
    delete[] phi;
    tme = GetTimer()/1000;
    return density;
}

double* calc_error(double hx, double hy, double* solution)
{
    double* res = new double[XY_LEN];
    for (int i = 0; i < OX_LEN_1; i++)
        for (int j = 0; j < OY_LEN_1; j++)
            res[i*OY_LEN_1 + j] = fabs(solution[i * OY_LEN_1 + j]
                                       - analytical_solution_circle(A + hx*i, C + hy*j));
    return res;
}

int main()
{
    double tme = 0.;

    double d = 0;

    for (int i = 0; i < 5; ++i)
    {
        switch(i)
        {
            case 0:
                d = 50.;
                break;
            case 1:
                d = 100.;
                break;
            case 2:
                d = 200.;
                break;
            case 3:
                d = 400.;
                break;
            case 4:
                d = 800.;
                break;
            case 5:
                d = 1600.;
                break;
        }

        A = 0.;
        B = 1.;
        C = 0.;
        D = 1.;
        R_SQ = 0.25*0.25;
        INN_DENSITY = 1.;
        OUT_DENSITY = 0.;
        JAK_ITER_CNT = 3000;

        TAU = 1./d;
        OX_LEN = (int)d;
        OY_LEN = (int)d;
        OX_LEN_1 = OX_LEN + 1;
        OY_LEN_1 = OY_LEN + 1;
        HX = (B - A)/OX_LEN;
        HY = (D - C)/OY_LEN;
        TIME_STEP_CNT = (int)d;
        XY_LEN = OX_LEN_1 * OY_LEN_1;
        printf("OX_LEN = %d OY_LEN = %d\n", OX_LEN, OY_LEN);
        double* density = solve(tme);
        double* err = calc_error(HX, HY, density);
        print_surface_as_v("rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, density);
        print_surface_as_v("err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, err);
        double l1 = get_l1_norm(HX, HY, OX_LEN_1, OY_LEN_1, err);
        double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
        printf("l1 %le \n", l1);
        printf("l_inf %le\n", l_inf);
        delete[] density;
        delete[] err;
    }

    return 0;
}