#include <math.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include "consts.h"
#include "timer.h"
#include "utils.h"
#include "common.h"
#include <map>
#include <iostream>

using namespace std;

inline void print_data_to_files(double *phi, double *density, double *residual, int tl) {
//    print_surface("phi", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
//                  U_VELOCITY, V_VELOCITY, phi);
    print_surface("rho", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                  U_VELOCITY, V_VELOCITY, density);
//    print_surface("res", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
//                  U_VELOCITY, V_VELOCITY, residual);
    double *err_lock = calc_error_3(HX, HY, tl * TAU, density);
    print_surface("err-l", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(),
                  TAU, U_VELOCITY, V_VELOCITY, err_lock);
    delete[] err_lock;
}

inline static double func_u(double t, double x, double y) { return (-y + CENTER_OFFSET_Y) * U_VELOCITY; }

inline static double func_v(double t, double x, double y) { return (x - CENTER_OFFSET_X) * V_VELOCITY; }

inline static double analytical_solution_circle(double t, double x, double y) {
    double x0 = get_center_x() + t * func_u(t, x, y);
    double y0 = get_center_y() + t * func_v(t, x, y);
    double value = (x - x0) * (x - x0) + (y - y0) * (y - y0);
    if (value <= R_SQ) return INN_DENSITY;
    return OUT_DENSITY;
}

static const int dot_c00 = 1;
static const int dot_c10 = 2;
static const int dot_c01 = 3;
static const int dot_c11 = 4;
static const int dot_upper = 5;
static const int dot_left = 6;
static const int dot_right = 7;
static const int dot_bott = 8;
static const int dot_inner = 9;

// fill for inner pt
static void fill_coef(int ii, int jj, double *density, std::map<int, double> *phiMap,
                      std::map<int, double> &coef, int type) {

    int sq_i = (int) ((ii * HX - A) / HX);
    int sq_j = (int) ((jj * HY - C) / HY);
    int key = OY_LEN_1 * sq_i + sq_j;

    double sum = 0.;
    for (int i = 0; i < XY_LEN; ++i)
        if (phiMap[i].count(key) > 0) sum += phiMap[i][key];

    double real_integral_value = 0;
    double part = 0.25;
    double integ = 0;
    switch (type) {
        case dot_c00: // CP00
            //break;
        case dot_c10:  // CP10
            //break;
        case dot_c01: // CP01
            //    break;
        case dot_c11: // CP11
//            break;
        case dot_upper: // upper
//            break;
        case dot_left: // left
//            break;
        case dot_right: // right
//            break;
        case dot_bott: // bott
//            break;
        case 9: // inner
            integ = (density[sq_i * OY_LEN_1 + sq_j]
                     + density[(sq_i + 1) * OY_LEN_1 + sq_j]
                     + density[(sq_i + 1) * OY_LEN_1 + sq_j + 1]
                     + density[sq_i * OY_LEN_1 + sq_j + 1]);
            break;
    }

    real_integral_value = HX * HY * part * integ;

    double k = sum != 0 ? real_integral_value / sum : 0.;

    if (coef.count(key) > 0 && coef[key] != k)
        printf("\nAlarm! Key %d exists in coef map!\nValue in map -> %f new value -> %f\n", key, coef[key], k);

    coef[key] = k;
}

static void fill_phi_map(int ii, int jj, double *density, double time_value, std::map<int, double> *phiMap) {
    double x1 = 0.;
    double y1 = 0.;
    double x2 = 0.;
    double y2 = 0.;
    double x3 = 0.;
    double y3 = 0.;
    double x4 = 0.;
    double y4 = 0.;

    get_coordinates_on_curr(ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);

//    printf("POINT: %d   %d :  x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le **
//                   x4=%.8le * y4=%.8le\n", ii,jj, x1,y1, x2,y2, x3,y3, x4,y4);

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
    /*
     if (x1 <= A || x1 >= B || x2 <= A || x2 >= B || x3 <= A || x3 >= B || x4 <= A || x4 >= B
        || y1 <= C || y1 >= D || y2 <= C || y2 >= D || y3 <= C || y3 >= D || y4 <= C || y4 >= D)
        printf("PREV Time level %.8le! ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", time_value, ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
*/
    int nx = IDEAL_SQ_SIZE_X;
    int ny = IDEAL_SQ_SIZE_Y;

    double x_step = 1. / nx;
    double y_step = 1. / ny;

    // get right part for jakoby
    double phi = 0.;
    double mes = x_step * y_step;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {

            double ideal_x = i * x_step + x_step / 2.;
            double ideal_y = j * y_step + y_step / 2.;

            double a11 = (x2 - x1) + (x1 + x3 - x2 - x4) * ideal_y;
            double a12 = (x4 - x1) + (x1 + x3 - x2 - x4) * ideal_x;
            double a21 = (y2 - y1) + (y1 + y3 - y2 - y4) * ideal_y;
            double a22 = (y4 - y1) + (y1 + y3 - y2 - y4) * ideal_x;
            double jakob = a11 * a22 - a21 * a12;

            double real_x = x1 + (x2 - x1) * ideal_x + (x4 - x1) * ideal_y
                            + (x1 + x3 - x2 - x4) * ideal_x * ideal_y;
            double real_y = y1 + (y2 - y1) * ideal_x + (y4 - y1) * ideal_y
                            + (y1 + y3 - y2 - y4) * ideal_x * ideal_y;

            // find out in which square real point was placed
            int sq_i = (int) ((real_x - A) / HX);
            int sq_j = (int) ((real_y - C) / HY);

            int mapIndex = OY_LEN_1 * sq_i + sq_j;

            double x = A + sq_i * HX;
            double y = C + sq_j * HY;

            // formula 4
            double dens = density[sq_i * OY_LEN_1 + sq_j] * (1 - (real_x - x) / HX) * (1 - (real_y - y) / HY)
                          + density[(sq_i + 1) * OY_LEN_1 + sq_j] * ((real_x - x) / HX) * (1 - (real_y - y) / HY)
                          + density[(sq_i + 1) * OY_LEN_1 + sq_j + 1] * ((real_x - x) / HX) * ((real_y - y) / HY)
                          + density[sq_i * OY_LEN_1 + sq_j + 1] * (1 - (real_x - x) / HX) * ((real_y - y) / HY);

            phiMap[OY_LEN_1 * ii + jj][mapIndex] = mes * dens * jakob;
        }
    }
}

static double get_phi_integ_midpoint(int ii, int jj, double *density, double time_value) {
    double x1 = 0.;
    double y1 = 0.;
    double x2 = 0.;
    double y2 = 0.;
    double x3 = 0.;
    double y3 = 0.;
    double x4 = 0.;
    double y4 = 0.;

    get_coordinates_on_curr(ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);

//    printf("POINT: %d   %d :  x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le **
//                   x4=%.8le * y4=%.8le\n", ii,jj, x1,y1, x2,y2, x3,y3, x4,y4);

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
    /*
     if (x1 <= A || x1 >= B || x2 <= A || x2 >= B || x3 <= A || x3 >= B || x4 <= A || x4 >= B
        || y1 <= C || y1 >= D || y2 <= C || y2 >= D || y3 <= C || y3 >= D || y4 <= C || y4 >= D)
        printf("PREV Time level %.8le! ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", time_value, ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
*/
    int nx = IDEAL_SQ_SIZE_X;
    int ny = IDEAL_SQ_SIZE_Y;

    double x_step = 1. / nx;
    double y_step = 1. / ny;

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

            // find out in which square real point was placed
            int sq_i = (int) ((real_x - A) / HX);
            int sq_j = (int) ((real_y - C) / HY);
            double x = A + sq_i * HX;
            double y = C + sq_j * HY;

            double a11 = (x2 - x1) + (x1 + x3 - x2 - x4) * ideal_y;
            double a12 = (x4 - x1) + (x1 + x3 - x2 - x4) * ideal_x;
            double a21 = (y2 - y1) + (y1 + y3 - y2 - y4) * ideal_y;
            double a22 = (y4 - y1) + (y1 + y3 - y2 - y4) * ideal_x;
            double jakob = a11 * a22 - a21 * a12;

            // formula 4
            double dens = density[sq_i * OY_LEN_1 + sq_j] * (1 - (real_x - x) / HX) * (1 - (real_y - y) / HY)
                          + density[(sq_i + 1) * OY_LEN_1 + sq_j] * ((real_x - x) / HX) * (1 - (real_y - y) / HY)
                          + density[(sq_i + 1) * OY_LEN_1 + sq_j + 1] * ((real_x - x) / HX) * ((real_y - y) / HY)
                          + density[sq_i * OY_LEN_1 + sq_j + 1] * (1 - (real_x - x) / HX) * ((real_y - y) / HY);

            phi += mes * dens * jakob;
        }
    }

    if (fabs(phi) < fabs(DBL_MIN_TRIM)) phi = 0;
    return phi;
}

static double get_phi_integ_midpoint(int ii, int jj, double *density, double time_value, std::map<int, double> coef) {
    double x1 = 0.;
    double y1 = 0.;
    double x2 = 0.;
    double y2 = 0.;
    double x3 = 0.;
    double y3 = 0.;
    double x4 = 0.;
    double y4 = 0.;

    get_coordinates_on_curr(ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);

//    printf("POINT: %d   %d :  x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le **
//                   x4=%.8le * y4=%.8le\n", ii,jj, x1,y1, x2,y2, x3,y3, x4,y4);

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
    /*
     if (x1 <= A || x1 >= B || x2 <= A || x2 >= B || x3 <= A || x3 >= B || x4 <= A || x4 >= B
        || y1 <= C || y1 >= D || y2 <= C || y2 >= D || y3 <= C || y3 >= D || y4 <= C || y4 >= D)
        printf("PREV Time level %.8le! ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", time_value, ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
*/
    int nx = IDEAL_SQ_SIZE_X;
    int ny = IDEAL_SQ_SIZE_Y;

    double x_step = 1. / nx;
    double y_step = 1. / ny;

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

            // find out in which square real point was placed
            int sq_i = (int) ((real_x - A) / HX);
            int sq_j = (int) ((real_y - C) / HY);
            double x = A + sq_i * HX;
            double y = C + sq_j * HY;

            double a11 = (x2 - x1) + (x1 + x3 - x2 - x4) * ideal_y;
            double a12 = (x4 - x1) + (x1 + x3 - x2 - x4) * ideal_x;
            double a21 = (y2 - y1) + (y1 + y3 - y2 - y4) * ideal_y;
            double a22 = (y4 - y1) + (y1 + y3 - y2 - y4) * ideal_x;
            double jakob = a11 * a22 - a21 * a12;

            // formula 4
            double dens = density[sq_i * OY_LEN_1 + sq_j] * (1 - (real_x - x) / HX) * (1 - (real_y - y) / HY)
                          + density[(sq_i + 1) * OY_LEN_1 + sq_j] * ((real_x - x) / HX) * (1 - (real_y - y) / HY)
                          + density[(sq_i + 1) * OY_LEN_1 + sq_j + 1] * ((real_x - x) / HX) * ((real_y - y) / HY)
                          + density[sq_i * OY_LEN_1 + sq_j + 1] * (1 - (real_x - x) / HX) * ((real_y - y) / HY);

            int key = OY_LEN_1 * sq_i + sq_j;
            double k = coef[key];

            phi += mes * dens * jakob * k;
        }
    }

    if (fabs(phi) < fabs(DBL_MIN_TRIM)) phi = 0;
    return phi;
}

double *solve_5(double &tme) {
    StartTimer();

    fflush(stdout);

    int ic = 0;
    double *phi = new double[XY_LEN];
    double *prev_density = new double[XY_LEN];
    double *density = new double[XY_LEN];
    double *residual = new double[XY_LEN];

    std::map<int, double> *phiMap = new std::map<int, double>[XY_LEN];
    for (int k = 0; k < XY_LEN; ++k) phiMap[k] = std::map<int, double>();

    std::map<int, double> *coef = new std::map<int, double>();

    //<editor-fold desc="Fill initial data">

    for (int i = 0; i < OX_LEN_1; ++i) {
        for (int j = 0; j < OY_LEN_1; ++j) {
            density[OY_LEN_1 * i + j] = 0.;
            prev_density[OY_LEN_1 * i + j] = 0.;
            residual[OY_LEN_1 * i + j] = 0.;
            phi[OY_LEN_1 * i + j] = 0.;
        }
    }

    // G1 -- (x_i, 0=C) -- bottom boundary
    for (int i = 0; i < OX_LEN_1; ++i) {
        prev_density[OY_LEN_1 * i] = analytical_solution_circle(0., A + HX * i, C);
        if (fabs(prev_density[OY_LEN_1 * i]) < fabs(DBL_MIN_TRIM)) prev_density[OY_LEN_1 * i] = 0;

    }

    // G2 -- (OX_LEN=B, y_j) -- right boundary
    for (int j = 1; j < OY_LEN; ++j) {
        prev_density[OY_LEN_1 * OX_LEN + j] = analytical_solution_circle(0., A + HX * OX_LEN, C + HY * j);
        if (fabs(prev_density[OY_LEN_1 * OX_LEN + j]) < fabs(DBL_MIN_TRIM)) prev_density[OY_LEN_1 * OX_LEN + j] = 0;

    }

    // G3 -- (x_i, OY_LEN=D) -- top boundary
    for (int i = 0; i < OX_LEN_1; ++i) {
        prev_density[OY_LEN_1 * i + OY_LEN] = analytical_solution_circle(0., A + HX * i, C + HY * OY_LEN);
        if (fabs(prev_density[OY_LEN_1 * i + OY_LEN]) < fabs(DBL_MIN_TRIM)) prev_density[OY_LEN_1 * i + OY_LEN] = 0;

    }

    // G4 -- (0=A, y_j) -- left boundary
    for (int j = 1; j < OY_LEN; ++j) {
        prev_density[j] = analytical_solution_circle(0., A, C + HY * j);
        if (fabs(prev_density[j]) < fabs(DBL_MIN_TRIM)) prev_density[j] = 0;
    }

    memcpy(density, prev_density, XY_LEN * sizeof(double));

    // inner points
    for (int i = 1; i < OX_LEN; ++i) {
        for (int j = 1; j < OY_LEN; ++j) {
            prev_density[OY_LEN_1 * i + j] = analytical_solution_circle(0., A + HX * i, C + HY * j);
            //if (fabs(prev_density[OY_LEN_1 * i + j]) < fabs(DBL_MIN_TRIM)) prev_density[OY_LEN_1 * i + j] = 0;
        }
    }

    //</editor-fold>

    printf("SUM RHO INIT = %le\n", calc_array_sum(prev_density, OX_LEN_1, OY_LEN_1, 0));
    printf("SUM ABS RHO INIT= %le\n", calc_array_sum(prev_density, OX_LEN_1, OY_LEN_1, 1));
    fflush(stdout);

    double maxRes = FLT_MAX;
    double *extrems;

    for (int tl = 1; tl <= TIME_STEP_CNT; tl++) {

        for (int k = 0; k < XY_LEN; ++k) {
            phiMap[k].clear();
        }
        coef->clear();

        //<editor-fold desc="Calculate phi">

        // with usage of prev_density we calculate phi function values

        // G1 -- (x_i, 0=C) -- bottom boundary
        for (int i = 1; i < OX_LEN; ++i) {
            if (G1[i] == 1) {
                fill_phi_map(i, 0, prev_density, TAU * tl, phiMap);
            }
        }
        for (int i = 1; i < OX_LEN; ++i) {
            if (G1[i] == 1) {
                fill_coef(i, 0, prev_density, phiMap, *coef, dot_bott);
            }
        }
        for (int i = 1; i < OX_LEN; ++i) {
            if (G1[i] == 1) {
                phi[OY_LEN_1 * i] = get_phi_integ_midpoint(i, 0, prev_density, TAU * tl, *coef);
            }
        }
        // end G1

        // G2 -- (OX_LEN=B, y_j) -- right boundary
        for (int j = 1; j < OY_LEN; ++j) {
            if (G2[j] == 1) {
                fill_phi_map(OX_LEN, j, prev_density, TAU * tl, phiMap);
            }
        }
        for (int j = 1; j < OY_LEN; ++j) {
            if (G2[j] == 1) {
                fill_coef(OX_LEN, j, prev_density, phiMap, *coef, dot_right);
            }
        }
        for (int j = 1; j < OY_LEN; ++j) {
            if (G2[j] == 1) {
                phi[OY_LEN_1 * OX_LEN + j] = get_phi_integ_midpoint(OX_LEN, j, prev_density, TAU * tl, *coef);
            }
        }
        // end G2

        // G3 -- (x_i, OY_LEN=D) -- upper boundary
        for (int i = 1; i < OX_LEN; ++i) {
            if (G3[i] == 1) {
                fill_phi_map(i, OY_LEN, prev_density, TAU * tl, phiMap);
            }
        }
        for (int i = 1; i < OX_LEN; ++i) {
            if (G3[i] == 1) {
                fill_coef(i, OY_LEN, prev_density, phiMap, *coef, dot_upper);
            }
        }
        for (int i = 1; i < OX_LEN; ++i) {
            if (G3[i] == 1) {
                phi[OY_LEN_1 * i + OY_LEN] = get_phi_integ_midpoint(i, OY_LEN, prev_density, TAU * tl, *coef);
            }
        }
        // end G3

        // G4 -- (0=A, y_j) -- left boundary
        for (int j = 1; j < OY_LEN; ++j) {
            if (G4[j] == 1) {
                fill_phi_map(0, j, prev_density, TAU * tl, phiMap);
            }
        }
        for (int j = 1; j < OY_LEN; ++j) {
            if (G4[j] == 1) {
                fill_coef(0, j, prev_density, phiMap, *coef, dot_bott);
            }
        }
        for (int j = 1; j < OY_LEN; ++j) {
            if (G4[j] == 1) {
                phi[j] = get_phi_integ_midpoint(0, j, prev_density, TAU * tl, *coef);
            }
        }
        // G4

        // point (0.0)
        if (CP00 == 1) {
            fill_phi_map(0, 0, prev_density, TAU * tl, phiMap);
            fill_coef(0, 0, prev_density, phiMap, *coef, dot_c00);
            phi[0] = get_phi_integ_midpoint(0, 0, prev_density, TAU * tl, *coef);
        }

        // point (1.0)
        if (CP10 == 1) {
            fill_phi_map(OX_LEN, 0, prev_density, TAU * tl, phiMap);
            fill_coef(OX_LEN, 0, prev_density, phiMap, *coef, dot_c10);
            phi[OY_LEN_1 * OX_LEN] = get_phi_integ_midpoint(OX_LEN, 0, prev_density, TAU * tl, *coef);
        }

        // point (0.1)
        if (CP01 == 1) {
            fill_phi_map(0, OY_LEN, prev_density, TAU * tl, phiMap);
            fill_coef(0, OY_LEN, prev_density, phiMap, *coef, dot_c01);
            phi[OY_LEN] = get_phi_integ_midpoint(0, OY_LEN, prev_density, TAU * tl, *coef);
        }

        // point (1,1)
        if (CP11 == 1) {
            fill_phi_map(OX_LEN, OY_LEN, prev_density, TAU * tl, phiMap);
            fill_coef(OX_LEN, OY_LEN, prev_density, phiMap, *coef, dot_c11);
            phi[OY_LEN_1 * OX_LEN + OY_LEN] = get_phi_integ_midpoint(OX_LEN, OY_LEN, prev_density, TAU * tl, *coef);
        }

        // inner points

        for (int i = 1; i < OX_LEN; ++i) {
            for (int j = 1; j < OY_LEN; ++j) {
                fill_phi_map(i, j, prev_density, TAU * tl, phiMap);
            }
        }

//        for (int i = 1; i < OX_LEN; ++i) {
//            for (int j = 1; j < OY_LEN; ++j) {
//                for (map<int, double>::const_iterator it = phiMap[OY_LEN_1 * i + j].begin();
//                        it != phiMap[OY_LEN_1 * i + j].end(); ++it) {
//                    std::cout << " key -> " << it->first << " value -> " << it->second << std::endl;
//                }
//            }
//        }

        for (int i = 1; i < OX_LEN; ++i) {
            for (int j = 1; j < OY_LEN; ++j) {
                fill_coef(i, j, prev_density, phiMap, *coef, dot_inner);
            }
        }

//        std::cout << " COEF " << std::endl;
//
//        for (map<int, double>::const_iterator it = coef->begin();
//             it != coef->end(); ++it) {
//            std::cout << " COEF key -> " << it->first << " value -> " << it->second << std::endl;
//        }

        for (int i = 1; i < OX_LEN; ++i) {
            for (int j = 1; j < OY_LEN; ++j) {
                phi[OY_LEN_1 * i + j] = get_phi_integ_midpoint(i, j, prev_density, TAU * tl, *coef);
            }
        }

        print_surface("phi1", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                      U_VELOCITY, V_VELOCITY, phi);

        //</editor-fold>

            // G1 -- (x_i, 0=C) -- bottom boundary
            double rpCoef = 2. / (HX * HY);
            for (int i = 1; i < OX_LEN; ++i) {
                if (G1[i] == 1) {
                    density[OY_LEN_1 * i] = rpCoef * phi[OY_LEN_1 * i];
                    if (fabs(density[OY_LEN_1 * i]) < fabs(DBL_MIN_TRIM)) density[OY_LEN_1 * i] = 0;
                }
            }

            // G2 -- (OX_LEN=B, y_j) -- right boundary
            for (int j = 1; j < OY_LEN; ++j) {
                if (G2[j] == 1) {
                    density[OY_LEN_1 * OX_LEN + j] = rpCoef * phi[OY_LEN_1 * OX_LEN + j];
                    if (fabs(density[OY_LEN_1 * OX_LEN + j]) < fabs(DBL_MIN_TRIM)) density[OY_LEN_1 * OX_LEN + j] = 0;
                }
            }

            // G3 -- (x_i, OY_LEN=D) -- top boundary
            for (int i = 1; i < OX_LEN; ++i) {
                if (G3[i] == 1) {
                    density[OY_LEN_1 * i + OY_LEN] = rpCoef * phi[OY_LEN_1 * i + OY_LEN];
                    if (fabs(density[OY_LEN_1 * i + OY_LEN]) < fabs(DBL_MIN_TRIM)) density[OY_LEN_1 * i + OY_LEN] = 0;
                }
            }

            // G4 -- (0=A, y_j) -- left boundary
            for (int j = 1; j < OY_LEN; ++j) {
                if (G4[j] == 1) { // проверить коэф-ты
                    density[j] = rpCoef * phi[j];
                    if (fabs(density[j]) < fabs(DBL_MIN_TRIM)) density[j] = 0;
                }
            }

            rpCoef = 4. / (HX * HY);

            // point (0,0)
            if (CP00 == 1) {
                density[0] = phi[0];
                if (fabs(density[0]) < fabs(DBL_MIN_TRIM)) density[0] = 0;
            }

            // point (1,0)
            if (CP10 == 1) {
                density[OY_LEN_1 * OX_LEN] = rpCoef * phi[OY_LEN_1 * OX_LEN];
                if (fabs(density[OY_LEN_1 * OX_LEN]) < fabs(DBL_MIN_TRIM))
                    density[OY_LEN_1 * OX_LEN] = 0;
            }

            // point (0,1)
            if (CP01 == 1) {
                density[OY_LEN] = rpCoef * phi[OY_LEN];
                if (fabs(density[OY_LEN]) < fabs(DBL_MIN_TRIM))
                    density[OY_LEN] = 0;
            }

            // point (1,1)
            if (CP11 == 1) {
                density[OY_LEN_1 * OX_LEN + OY_LEN] = rpCoef * phi[OY_LEN_1 * OX_LEN + OY_LEN];
                if (fabs(density[OY_LEN_1 * OX_LEN + OY_LEN]) < fabs(DBL_MIN_TRIM))
                    density[OY_LEN_1 * OX_LEN + OY_LEN] = 0;
            }

            rpCoef = 1. / (HX * HY);
            for (int i = 1; i < OX_LEN; ++i) {
                for (int j = 1; j < OY_LEN; ++j) {
                    density[OY_LEN_1 * i + j] = rpCoef * phi[OY_LEN_1 * i + j];
                    if (fabs(density[OY_LEN_1 * i + j]) < fabs(DBL_MIN_TRIM)) density[OY_LEN_1 * i + j] = 0;
                }
            }

        memcpy(prev_density, density, XY_LEN * sizeof(double));

        printf("tl = %d Sum(Rho) = %le Sum(absRho) = %le\n",
               tl, calc_array_sum(density, OX_LEN_1, OY_LEN_1, 0), calc_array_sum(density, OX_LEN_1, OY_LEN_1, 1));
        fflush(stdout);

        if (tl % 5 == 0) {
            print_data_to_files(phi, density, residual, tl);
            /*int fixed_x = (int) (get_center_x() / HX);
            int fixed_y = (int) (get_center_y() / HY);
            print_line_along_x("rho", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                               U_VELOCITY, V_VELOCITY, density, fixed_y);
            print_line_along_y("rho", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                               U_VELOCITY, V_VELOCITY, density, fixed_x);*/
        }
    }

    double *err = calc_error_5(HX, HY, TAU * TIME_STEP_CNT, density);
    double l1_err_vec = get_l1_norm_vec(OX_LEN_1, OY_LEN_1, err);
    double l1_err_tr = get_l1_norm_int_trapezoidal(HX, HY, OX_LEN, OY_LEN, err); // note! a loop boundary
//    append_statistics(OX_LEN_1, OY_LEN_1, TAU, ic, l1_err_vec, l1_err_tr, maxRes, TIME_STEP_CNT);
    extrems = calc_array_extrems(density, OX_LEN_1, OY_LEN_1);
    append_statistics(OX_LEN_1, OY_LEN_1, TAU, ic, l1_err_vec, l1_err_tr, maxRes, extrems,
                      extrems, TIME_STEP_CNT); // !!!!!!!! tmp stab


    delete[] prev_density;
    delete[] phi;
    delete[] err;
    delete[] residual;
    delete[] extrems;
    delete[] phiMap;
    delete coef;
    tme = GetTimer() / 1000;
    return density;
}

double *calc_error_5(double hx, double hy, double tt, double *solution) {
    double *res = new double[XY_LEN];
    for (int i = 0; i < OX_LEN_1; i++)
        for (int j = 0; j < OY_LEN_1; j++)
            res[i * OY_LEN_1 + j] = fabs(solution[i * OY_LEN_1 + j]
                                         - analytical_solution_circle(tt, A + hx * i, C + hy * j));
    return res;
}

double *get_exact_solution_5(double hx, double hy, double t) {
    double *res = new double[XY_LEN];
    for (int i = 0; i < OX_LEN_1; i++)
        for (int j = 0; j < OY_LEN_1; j++)
            res[i * OY_LEN_1 + j] = fabs(analytical_solution_circle(t, A + hx * i, C + hy * j));
    return res;
}