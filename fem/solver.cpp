#include <math.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include "consts.h"
#include "solver.h"
#include "tecplot.h"
#include <algorithm>

#ifdef WIN32
__int64 timerStart = 0;
#else
struct timeval timerStart;
#endif

#include "timer.h"


double A = 0.;
double B = 0.;
double C = 0.;
double D = 0.;
int OX_LEN = 0;
int OY_LEN = 0;
int OX_LEN_1 = 0;
int OY_LEN_1 = 0;
int XY_LEN = 0;
double TAU = 0.;
int TIME_STEP_CNT = 0;
int JAK_ITER_CNT = 0;
double HX = 0.;
double HY = 0.;
double R_SQ = 0.;
double INN_DENSITY = 0.;
double OUT_DENSITY = 0.;

inline static double analytical_solution_circle(double x, double y) {
    double x0 = A + OX_LEN * HX / 2.;
    double y0 = C + OY_LEN * HY / 2.;
    double value = (x - x0) * (x - x0) + (y - y0) * (y - y0);
    if (value <= R_SQ) return INN_DENSITY;
    return OUT_DENSITY;
}

inline static double func_u(double time_value, double x, double y) { return 0; }

inline static double func_v(double time_value, double x, double y) { return 0; }

static double get_phi(int ii, int jj, double *prev_density, double time_value) {
    double x1 = A + ii * HX - HX / 2.;
    double y1 = C + jj * HY - HY / 2.;
    double x2 = A + ii * HX + HX / 2.;
    double y2 = C + jj * HY - HY / 2.;
    double x3 = A + ii * HX + HX / 2.;
    double y3 = C + jj * HY + HY / 2.;
    double x4 = A + ii * HX - HX / 2.;
    double y4 = C + jj * HY + HY / 2.;

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
    double x_step = 1. / nx;
    double y_step = 1. / ny;

    // point 0, 0

    // point 1, 0

    // point 0, 1

    // point 1, 1

    // G1 left boundary

    // G2 bottom boundary

    // G3 right boundary

    // G4 top boundary

    // get right part for jakoby
    double phi = 0.;
    double mes = x_step * y_step;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {

            double ideal_x = i * x_step + x_step / 2.;
            double ideal_y = j * y_step + y_step / 2.;
            double real_x = x1 + (x2 - x1) * ideal_x + (x4 - x1) * ideal_y
                            + (x1 + x3 - x2 - x4) * ideal_x * ideal_y;
            double real_y = y1 + (y2 - y1) * ideal_x + (y4 - y1) * ideal_y
                            + (y1 + y3 - y2 - y4) * ideal_x * ideal_y;

            if (ii > 0 && ii < OX_LEN && jj > 0 && jj < OY_LEN) {

            }
            else if (ii == 0 && jj == 0) { // point (0,0)

            }
            else if (ii == 0 && jj == OY_LEN) { // point (0,1)

            }
            else if (ii == OX_LEN && jj == 0) { // point (1,0)

            }
            else if (ii == OX_LEN && jj == OY_LEN) { // point (1,1)

            }
            else if (ii == 0 && jj > 0 && jj < OY_LEN) {  // G1 left boundary

            }
            else if (jj == 0 && ii > 0 && ii < OX_LEN) { // G2 bottom boundary

            }
            else if (ii == OX_LEN && jj > 0 && jj < OY_LEN) { // G3 right boundary

            }
            else if (jj == OY_LEN && ii > 0 && ii < OX_LEN) { // G4 top boundary

            }

            // find out in which square real point was placed
            int sq_i = (int) ((real_x - A) / HX);
            int sq_j = (int) ((real_y - C) / HY);
            double x = A + sq_i * HX;
            double y = C + sq_j * HY;

            // formula 4
            double dens = prev_density[sq_i * OY_LEN_1 + sq_j] * (1 - (real_x - x) / HX) * (1 - (real_y - y) / HY)
                          + prev_density[(sq_i + 1) * OY_LEN_1 + sq_j] * ((real_x - x) / HX) * (1 - (real_y - y) / HY)
                          + prev_density[(sq_i + 1) * OY_LEN_1 + sq_j + 1] * ((real_x - x) / HX) * ((real_y - y) / HY)
                          + prev_density[sq_i * OY_LEN_1 + sq_j + 1] * (1 - (real_x - x) / HX) * ((real_y - y) / HY);

            double a11 = (x2 - x1) + (x1 + x3 - x2 - x4) * ideal_y;
            double a12 = (x4 - x1) + (x1 + x3 - x2 - x4) * ideal_x;
            double a21 = (y2 - y1) + (y1 + y3 - y2 - y4) * ideal_y;
            double a22 = (y4 - y1) + (y1 + y3 - y2 - y4) * ideal_x;
            double jakob = a11 * a22 - a21 * a12;

            phi += mes * dens * jakob;
        }
    }

    phi *= 16. / 9.;
    return phi;
}

double *solve(double &tme) {
    StartTimer();

    fflush(stdout);

    double *phi = new double[XY_LEN];
    double *prev_density = new double[XY_LEN];
    double *density = new double[XY_LEN];

    for (int i = 0; i < OX_LEN_1; ++i)
        for (int j = 0; j < OY_LEN_1; ++j)
            prev_density[OY_LEN_1 * i + j] = analytical_solution_circle(HX * i, HY * j);

    memcpy(density, prev_density, XY_LEN * sizeof(double));

    for (int tl = 1; tl <= TIME_STEP_CNT; tl++) {

        // with usage of prev_density we calculate phi function values
        for (int i = 0; i < OX_LEN_1; ++i)
            for (int j = 0; j < OY_LEN_1; ++j)
                phi[i * OY_LEN_1 + j] = get_phi(i, j, prev_density, TAU * tl);

//        if (tl == TIME_STEP_CNT)
//            print_surface_as_v("phi", OX_LEN, OY_LEN, HX, HY, tl, A, C, phi);

        if (tl == 1) {
            // fill inner matrix of prev_density by zero
            // because we will use it in Jakoby method
            for (int i = 0; i < OX_LEN_1; ++i)
                for (int j = 0; j < OY_LEN_1; ++j)
                    prev_density[i * OY_LEN_1 + j] = 0.;
        }
        else {
            // none
        }

        int ic = 0;
        double maxErr = FLT_MAX;
        while (maxErr > EPS && ic < JAK_ITER_CNT) {
            double rpCoef = 64. / (9. * HX * HY);

            // point 0,0
            density[OY_LEN_1 * 0 + 0] = -1. / 3. * (prev_density[OY_LEN_1 * 0 + 1] + prev_density[OY_LEN_1 * 1 + 0]) -
                                        1. / 9. * prev_density[OY_LEN_1 * 1 + 1]
                                        + rpCoef * phi[OY_LEN_1 * 0 + 0];

            // point 1,0 (i = OX_LEN, j = 0)
            density[OY_LEN_1 * 0 + OX_LEN] = -1. / 3. * (prev_density[OY_LEN_1 * 0 + (OX_LEN - 1)]
                                                         + prev_density[OY_LEN_1 * 1 + OX_LEN])
                                             - 1. / 9. * prev_density[OY_LEN_1 * 1 + (OX_LEN - 1)]
                                             + rpCoef * phi[OY_LEN_1 * 0 + OX_LEN];

            // point 1,1
            density[OY_LEN_1 * OX_LEN + OX_LEN] = -1. / 3. * (prev_density[OY_LEN_1 * (OX_LEN - 1) + OX_LEN]
                                                              + prev_density[OY_LEN_1 * OX_LEN + OX_LEN - 1])
                                                  - 1. / 9. * prev_density[OY_LEN_1 * (OX_LEN - 1) + (OX_LEN - 1)]
                                                  + rpCoef * phi[OY_LEN_1 * OX_LEN + OX_LEN];

            // point 0,1 i = 0, j = OY_LEN
            density[OY_LEN_1 * OX_LEN + 0] = -1. / 3. * (prev_density[OY_LEN_1 * OX_LEN + 1]
                                                         + prev_density[OY_LEN_1 * (OX_LEN - 1) + 0])
                                             - 1. / 9. * prev_density[OY_LEN_1 * (OX_LEN - 1) + 1]
                                             + rpCoef * phi[OY_LEN_1 * OX_LEN + 0];

            // print_matrix_to_file(OX_LEN_1, OY_LEN_1, prev_density, "prev_density_test.dat");

            double bdCoef = 32. / (9. * HX * HY);
            // G1 left boundary
            for (int j = 1; j < OY_LEN; ++j) {
                density[OY_LEN_1 * j + 0] = -3. * prev_density[OY_LEN_1 * j + 1]
                                            - 1.5
                                              * (prev_density[OY_LEN_1 * (j + 1) + 0]
                                                 + prev_density[OY_LEN_1 * (j - 1) + 0]
                                              )
                                            - 0.5 * (prev_density[OY_LEN_1 * (j + 1) +
                                                                  1] +
                                                     prev_density[OY_LEN_1 * (j - 1) +
                                                                  1]) + bdCoef * phi[OY_LEN_1 * j + 0];
            }

            // G2 bottom boundary
            for (int i = 1; i < OX_LEN; ++i) {
                density[OY_LEN_1 * 0 + i] = -3. * prev_density[OY_LEN_1 * 1 + i]
                                            - 1.5 *
                                              (prev_density[OY_LEN_1 * 0 + i + 1] + prev_density[OY_LEN_1 * 0 + i - 1])
                                            - 0.5 *
                                              (prev_density[OY_LEN_1 * 1 + i + 1] + prev_density[OY_LEN_1 * 1 + i - 1])
                                            + bdCoef * phi[OY_LEN_1 * 0 + i];
            }

            // G3 right boundary
            for (int j = 1; j < OY_LEN; ++j) {
                density[OY_LEN_1 * j + OX_LEN] = -3. * prev_density[OY_LEN_1 * j + OX_LEN - 1]
                                                 - 1.5 * (prev_density[OY_LEN_1 * (j + 1) + OX_LEN] +
                                                          prev_density[OY_LEN_1 * (j - 1) + OX_LEN])
                                                 - 0.5 * (prev_density[OY_LEN_1 * (j + 1) + OX_LEN - 1] +
                                                          prev_density[OY_LEN_1 * (j - 1) + OX_LEN - 1])
                                                 + bdCoef * phi[OY_LEN_1 * j + OX_LEN];
            }

            // G4 top boundary
            for (int i = 1; i < OX_LEN; ++i) {
                density[OY_LEN_1 * OX_LEN + i] = -3. * prev_density[OY_LEN_1 * (OX_LEN - 1) + i]
                                                 - 1.5 * (prev_density[OY_LEN_1 * OX_LEN + i - 1] +
                                                          prev_density[OY_LEN_1 * OX_LEN + i + 1])
                                                 - 0.5 * (prev_density[OY_LEN_1 * (OX_LEN - 1) + i - 1] +
                                                          prev_density[OY_LEN_1 * (OX_LEN - 1) + i + 1])
                                                 + bdCoef * phi[OY_LEN_1 * OX_LEN + i];
            }


            for (int i = 1; i < OX_LEN; ++i) {
                for (int j = 1; j < OY_LEN; ++j) {
                    density[OY_LEN_1 * i + j] = -1. / 9. * (
                            1.5 * (
                                    prev_density[OY_LEN_1 * i + j - 1] + // left
                                    prev_density[OY_LEN_1 * (i - 1) + j] + // upper
                                    prev_density[OY_LEN_1 * i + j + 1] + // right
                                    prev_density[OY_LEN_1 * (i + 1) + j] // bottom
                            ) + 0.25 * (
                                    prev_density[OY_LEN_1 * (i + 1) + j + 1] + // bottom right
                                    prev_density[OY_LEN_1 * (i + 1) + j - 1] + // bottom left
                                    prev_density[OY_LEN_1 * (i - 1) + j - 1] + // upper right
                                    prev_density[OY_LEN_1 * (i - 1) + j + 1] // upper left
                            )) + phi[OY_LEN_1 * i + j] / (HX * HY);
                }
            }
            ++ic;

            maxErr = FLT_MIN;
            for (int i = 1; i < OX_LEN; ++i) {
                for (int j = 1; j < OY_LEN; ++j) {
                    double val = fabs(density[i * OY_LEN_1 + j] - prev_density[i * OY_LEN_1 + j]);
                    if (val > maxErr) { maxErr = val; }
                }
            }

            if (tl == 1 && ic == 1) {
                print_surface_as_v("prev_density", OX_LEN, OY_LEN, HX, HY, tl, A, C, prev_density);
                print_surface_as_v("density", OX_LEN, OY_LEN, HX, HY, tl, A, C, density);
            }

            memcpy(prev_density, density, XY_LEN * sizeof(double));

        }
    }

    delete[] prev_density;
    delete[] phi;
    tme = GetTimer() / 1000;
    return density;
}

double *calc_error(double hx, double hy, double *solution) {
    double *res = new double[XY_LEN];
    for (int i = 0; i < OX_LEN_1; i++)
        for (int j = 0; j < OY_LEN_1; j++)
            res[i * OY_LEN_1 + j] = fabs(solution[i * OY_LEN_1 + j]
                                         - analytical_solution_circle(A + hx * i, C + hy * j));
    return res;
}

