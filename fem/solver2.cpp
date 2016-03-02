#include <math.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include "consts.h"
#include "solver2.h"
#include "timer.h"
#include "utils.h"

double get_center_x_2() { return A + 0.3; }

double get_center_y_2() { return C + 0.3; }

inline static double func_u(double t, double x, double y) { return U_VELOCITY; }

inline static double func_v(double t, double x, double y) { return V_VELOCITY; }

inline static double analytical_solution_circle(double t, double x, double y) {
    double x0 = get_center_x_2() + t * func_u(t, x, y);
    double y0 = get_center_y_2() + t * func_v(t, x, y);
    double value = (x - x0) * (x - x0) + (y - y0) * (y - y0);
    //if (value <= R_SQ) return INN_DENSITY;
    if (value < R_SQ) return INN_DENSITY;
    return OUT_DENSITY;
}

static double get_phi(int ii, int jj, double *density, double time_value) {
    double x1 = 0.;
    double y1 = 0.;
    double x2 = 0.;
    double y2 = 0.;
    double x3 = 0.;
    double y3 = 0.;
    double x4 = 0.;
    double y4 = 0.;

    if (ii > 0 && ii < OX_LEN && jj > 0 && jj < OY_LEN) {
        // p1
        x1 = A + ii * HX - HX / 2.;
        y1 = C + jj * HY - HY / 2.;
        // p2
        x2 = A + ii * HX - HX / 2.;
        y2 = C + jj * HY + HY / 2.;
        // p3
        x3 = A + ii * HX + HX / 2.;
        y3 = C + jj * HY + HY / 2.;
        // p4
        x4 = A + ii * HX + HX / 2.;
        y4 = C + jj * HY - HY / 2.;
    }
    else if (ii == OX_LEN && jj == OY_LEN) { // point (1,1)  omega_{i-1,j-1}
        // p1
        x1 = B - HX / 2.;
        y1 = D - HY / 2.;
        // p2
        x2 = A - HX / 2.;
        y2 = D;
        // p3
        x3 = B;
        y3 = D;
        // p4
        x4 = B;
        y4 = D - HY / 2.;
    }
    else if (jj == OY_LEN && ii > 0 && ii < OX_LEN) { // G3 right boundary
        // p1
        x1 = A + ii * HX - HX / 2.;
        y1 = D - HY / 2.;
        // p2
        x2 = A + ii * HX - HX / 2.;
        y2 = D;
        // p3
        x3 = A + ii * HX + HX / 2.;
        y3 = D;
        // p4
        x4 = A + ii * HX - HX / 2.;
        y4 = D - HY / 2.;
    }
    else if (ii == 0 && jj > 0 && jj < OY_LEN) { // G4 top boundary
        // p1
        x1 = B - HX / 2.;
        y1 = C + jj * HY - HY / 2.;
        // p2
        x2 = B - HX / 2.;
        y2 = C + jj * HY + HY / 2.;
        // p3
        x3 = B;
        y3 = C + jj * HY + HY / 2.;
        // p4
        x4 = B;
        y4 = C + jj * HY - HY / 2.;
    }

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

    int nx = 64;
    int ny = 64;
    int nx_1 = nx + 1;
    int ny_1 = ny + 1;

    double x_step = 1. / nx;
    double y_step = 1. / ny;

    // get right part for jakoby
    double phi = 0.;
    double mes = x_step * y_step;
    for (int i = 0; i < nx_1; ++i) {
        for (int j = 0; j < ny_1; ++j) {

            double ideal_x = i * x_step + x_step / 2.;
            double ideal_y = j * y_step + y_step / 2.;
            double real_x = x1 + (x2 - x1) * ideal_x + (x4 - x1) * ideal_y
                            + (x1 + x3 - x2 - x4) * ideal_x * ideal_y;
            double real_y = y1 + (y2 - y1) * ideal_x + (y4 - y1) * ideal_y
                            + (y1 + y3 - y2 - y4) * ideal_x * ideal_y;

            // find out in which square real point was placed
            int sq_i = (int) ((real_x - A) / HX);
            int sq_j = (int) ((real_y - C) / HY);
            double x = A + sq_i * HX;
            double y = C + sq_j * HY;

            // formula 4
            double dens = density[sq_i * OY_LEN_1 + sq_j] * (1 - (real_x - x) / HX) * (1 - (real_y - y) / HY)
                          + density[(sq_i + 1) * OY_LEN_1 + sq_j] * ((real_x - x) / HX) * (1 - (real_y - y) / HY)
                          + density[(sq_i + 1) * OY_LEN_1 + sq_j + 1] * ((real_x - x) / HX) * ((real_y - y) / HY)
                          + density[sq_i * OY_LEN_1 + sq_j + 1] * (1 - (real_x - x) / HX) * ((real_y - y) / HY);

            double a11 = (x2 - x1) + (x1 + x3 - x2 - x4) * ideal_y;
            double a12 = (x4 - x1) + (x1 + x3 - x2 - x4) * ideal_x;
            double a21 = (y2 - y1) + (y1 + y3 - y2 - y4) * ideal_y;
            double a22 = (y4 - y1) + (y1 + y3 - y2 - y4) * ideal_x;
            double jakob = a11 * a22 - a21 * a12;

            phi += mes * dens * jakob;
        }
    }

    return phi;
}

double *solve_2(double &tme) {
    StartTimer();

    fflush(stdout);

    int ic = 0;
    double *phi = new double[XY_LEN];
    double *prev_density = new double[XY_LEN];
    double *density = new double[XY_LEN];

    for (int i = 0; i < OX_LEN_1; ++i) {
        for (int j = 0; j < OY_LEN_1; ++j) {
            density[OY_LEN_1 * i + j] = 0.;
            prev_density[OY_LEN_1 * i + j] = 0.;
        }
    }

    // G1
    for (int i = 0; i < OX_LEN_1; ++i)
        prev_density[OY_LEN_1 * i] = 0.;

    // G2
    for (int j = 0; j < OY_LEN_1; ++j)
        prev_density[j] = 0.;

    // G3
    for (int i = 1; i < OX_LEN; ++i)
        prev_density[OY_LEN_1 * i + OY_LEN] = analytical_solution_circle(0., A + HX * i, C + HY * OY_LEN);

    // G4
    for (int j = 1; j < OY_LEN; ++j)
        prev_density[OY_LEN_1 * OX_LEN + j] = analytical_solution_circle(0., HX * OX_LEN, HY * j);

    // (1,1)
    prev_density[OY_LEN_1 * OY_LEN + OX_LEN] = analytical_solution_circle(0., HX * OX_LEN, HY * OY_LEN);

    // inner points
    for (int i = 1; i < OX_LEN; ++i)
        for (int j = 1; j < OY_LEN; ++j)
            prev_density[OY_LEN_1 * i + j] = analytical_solution_circle(0., A + HX * i, C + HY * j);

    for (int tl = 1; tl <= TIME_STEP_CNT; tl++) {
        // with usage of prev_density we calculate phi function values

        // G3
        for (int i = 1; i < OX_LEN; ++i)
            phi[OY_LEN_1 * i + OY_LEN] = get_phi(i, OY_LEN, prev_density, TAU * tl);

        // G4
        for (int j = 1; j < OY_LEN; ++j)
            phi[j] = get_phi(OX_LEN, j, prev_density, TAU * tl);

        // point (1,1)
        phi[OY_LEN_1 * OX_LEN + OY_LEN] = get_phi(OX_LEN, OY_LEN, prev_density, TAU * tl);

        // inner points
        for (int i = 1; i < OX_LEN; ++i)
            for (int j = 1; j < OY_LEN; ++j)
                phi[OY_LEN_1 * i + j] = get_phi(i, j, prev_density, TAU * tl);

        ic = 0;
        double maxErr = FLT_MAX;
        while (maxErr > EPS) {
            double rpCoef = 64. / (9. * HX * HY);

            // point 1,1
            density[OY_LEN_1 * OX_LEN + OY_LEN] = -1. / 3. * (prev_density[OY_LEN_1 * OX_LEN + OY_LEN - 1]
                                                              + prev_density[OY_LEN_1 * (OX_LEN - 1) + OY_LEN])
                                                  - 1. / 9. * prev_density[OY_LEN_1 * (OX_LEN - 1) + OY_LEN - 1]
                                                  + rpCoef * phi[OY_LEN_1 * OX_LEN + OY_LEN];

            // G3 right boundary
            for (int i = 1; i < OX_LEN; ++i) {
                density[OY_LEN_1 * i + OY_LEN] = -2. / 9. * prev_density[OY_LEN_1 * i + OY_LEN - 1]
                                                 - 1. / 6. * (prev_density[OY_LEN_1 * (i + 1) + OY_LEN] +
                                                              prev_density[OY_LEN_1 * (i - 1) + OY_LEN])
                                                 - 1. / 18. * (prev_density[OY_LEN_1 * (i + 1) + OY_LEN - 1] +
                                                               prev_density[OY_LEN_1 * (i - 1) + OY_LEN - 1])
                                                 + (32. / (9. * HX * HY)) * phi[OY_LEN_1 * i + OY_LEN];
            }

            // G4 top boundary
            for (int j = 1; j < OY_LEN; ++j) {
                density[OY_LEN_1 * OX_LEN + j] = -2. / 9. * prev_density[OY_LEN_1 * (OX_LEN - 1) + j]
                                                 - 1. / 6. * (prev_density[OY_LEN_1 * OX_LEN + j - 1] +
                                                              prev_density[OY_LEN_1 * OX_LEN + j + 1])
                                                 - 1. / 18. * (prev_density[OY_LEN_1 * (OX_LEN - 1) - 1] +
                                                               prev_density[OY_LEN_1 * (OX_LEN - 1) + 1])
                                                 + (32. / (9. * HX * HY)) * phi[OY_LEN_1 * OX_LEN + j];
            }

            for (int i = 1; i < OX_LEN; ++i) {
                for (int j = 1; j < OY_LEN; ++j) {
                    density[OY_LEN_1 * i + j] = -1. / 6. * (
                            prev_density[OY_LEN_1 * i + j - 1] + // left
                            prev_density[OY_LEN_1 * (i - 1) + j] + // upper
                            prev_density[OY_LEN_1 * i + j + 1] + // right
                            prev_density[OY_LEN_1 * (i + 1) + j] // bottom
                    ) - 1. / 36. * (
                            prev_density[OY_LEN_1 * (i + 1) + j + 1] + // bottom right
                            prev_density[OY_LEN_1 * (i + 1) + j - 1] + // bottom left
                            prev_density[OY_LEN_1 * (i - 1) + j + 1] + // upper right
                            prev_density[OY_LEN_1 * (i - 1) + j - 1]  // upper left
                    ) + (16. * phi[OY_LEN_1 * i + j]) / (9. * HX * HY);
                }
            }
            ++ic;

            maxErr = FLT_MIN;
            for (int i = 0; i < OX_LEN_1; ++i)
                for (int j = 0; j < OY_LEN_1; ++j) {
                    double val = fabs(density[i * OY_LEN_1 + j] - prev_density[i * OY_LEN_1 + j]);
                    if (val > maxErr) maxErr = val;
                }

            memcpy(prev_density, density, XY_LEN * sizeof(double));
        }
        printf("Current tl = %d Iteration count = %d\n", tl, ic);
        fflush(stdout);
    }
    double *err = calc_error_2(HX, HY, density);
    double l1_err = get_l1_norm(HX, HY, OX_LEN_1, OY_LEN_1, err);
    append_statistics(OX_LEN_1, OY_LEN_1, TAU, ic, l1_err, TIME_STEP_CNT);
    delete[] prev_density;
    delete[] phi;
    delete[] err;
    tme = GetTimer() / 1000;
    return density;
}

double *calc_error_2(double hx, double hy, double *solution) {
    double *res = new double[XY_LEN];
    for (int i = 0; i < OX_LEN_1; i++)
        for (int j = 0; j < OY_LEN_1; j++)
            res[i * OY_LEN_1 + j] = fabs(solution[i * OY_LEN_1 + j]
                                         - analytical_solution_circle(TAU * TIME_STEP_CNT, A + hx * i, C + hy * j));
    return res;
}

double* get_exact_solution_2(double hx, double hy, double t)
{
    double *res = new double[XY_LEN];
    for (int i = 0; i < OX_LEN_1; i++)
        for (int j = 0; j < OY_LEN_1; j++)
            res[i * OY_LEN_1 + j] = fabs(analytical_solution_circle(t, A + hx * i, C + hy * j));
    return res;
}