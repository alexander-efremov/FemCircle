#include <math.h>#include <stdio.h>#include <string.h>#include <float.h>#include "consts.h"#include "timer.h"#include "utils.h"#include "common.h"#include "algorithm"using namespace std;inline void print_data_to_files(double *phi, double *density, double *residual, int tl) {    print_surface("phi", NX, NY, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,                  U, V, phi);    print_surface("rho", NX, NY, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,                  U, V, density);//    print_surface("res", NX, NY, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,//                  U, V, residual);    double *err_lock = calc_error_6(HX, HY, tl * TAU, density);    print_surface("err-l", NX, NY, HX, HY, tl, A, C, get_center_x(), get_center_y(),                  TAU, U, V, err_lock);    delete[] err_lock;}inline static double func_u(double t, double x, double y) { return U; }inline static double func_v(double t, double x, double y) { return V; }inline static double analytical_solution_circle(double t, double x, double y) {    double x0 = get_center_x() + t * func_u(t, x, y);    double y0 = get_center_y() + t * func_v(t, x, y);    double value = (x - x0) * (x - x0) + (y - y0) * (y - y0);    if (value <= R_SQ) return INN_DENSITY;    return OUT_DENSITY;}static inline double get_distance(double real_x, double real_y, double x, double y) {    double sqX = (real_x - x) * (real_x - x);    double sqY = (real_y - y) * (real_y - y);    double r = sqrt(sqX + sqY);    return r;}static double calc_phi(int ii, int jj, double *density, double time_value, double hx, double hy) {    double x1 = 0.;    double y1 = 0.;    double x2 = 0.;    double y2 = 0.;    double x3 = 0.;    double y3 = 0.;    double x4 = 0.;    double y4 = 0.;    get_coordinates_on_curr(ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);//    printf("POINT: %d   %d :  x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le **//                   x4=%.8le * y4=%.8le\n", ii,jj, x1,y1, x2,y2, x3,y3, x4,y4);    double u = func_u(time_value, x1, y1);    double v = func_v(time_value, x1, y1);    x1 = x1 - TAU * u;    y1 = y1 - TAU * v;    u = func_u(time_value, x2, y2);    v = func_v(time_value, x2, y2);    x2 = x2 - TAU * u;    y2 = y2 - TAU * v;    u = func_u(time_value, x3, y3);    v = func_v(time_value, x3, y3);    x3 = x3 - TAU * u;    y3 = y3 - TAU * v;    u = func_u(time_value, x4, y4);    v = func_v(time_value, x4, y4);    x4 = x4 - TAU * u;    y4 = y4 - TAU * v;    /*     if (x1 <= A || x1 >= B || x2 <= A || x2 >= B || x3 <= A || x3 >= B || x4 <= A || x4 >= B        || y1 <= C || y1 >= D || y2 <= C || y2 >= D || y3 <= C || y3 >= D || y4 <= C || y4 >= D)        printf("PREV Time level %.8le! ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "                       "x4=%.8le * y4=%.8le\n ", time_value, ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);*/    int nx = IDEAL_SQ_SIZE_X;    int ny = IDEAL_SQ_SIZE_Y;    double x_step = 1. / nx;    double y_step = 1. / ny;    // get right part for jakoby    double phi = 0.;    double mes = x_step * y_step;    for (int i = 0; i < nx; ++i) {        for (int j = 0; j < ny; ++j) {            double ideal_x = i * x_step + x_step / 2.;            double ideal_y = j * y_step + y_step / 2.;            double real_x = x1 + (x2 - x1) * ideal_x + (x4 - x1) * ideal_y                            + (x1 + x3 - x2 - x4) * ideal_x * ideal_y;            double real_y = y1 + (y2 - y1) * ideal_x + (y4 - y1) * ideal_y                            + (y1 + y3 - y2 - y4) * ideal_x * ideal_y;            // find out in which square real point was placed            int sq_i = (int) (((real_x - A) + FLT_MIN) / hx);            int sq_j = (int) (((real_y - C) + FLT_MIN) / hy);            double x = A + sq_i * HX;            double y = C + sq_j * HY;            double a11 = (x2 - x1) + (x1 + x3 - x2 - x4) * ideal_y;            double a12 = (x4 - x1) + (x1 + x3 - x2 - x4) * ideal_x;            double a21 = (y2 - y1) + (y1 + y3 - y2 - y4) * ideal_y;            double a22 = (y4 - y1) + (y1 + y3 - y2 - y4) * ideal_x;            double jakob = a11 * a22 - a21 * a12;            // left bottom            x1 = x;            y1 = y;            // right bottom            x2 = x + hx;            y2 = y;            // right top            x3 = x + hx;            y3 = y + hy;            // left top            x4 = x;            y4 = y + hy;            double d1 = get_distance(real_x, real_y, x1, y1);            double d2 = get_distance(real_x, real_y, x2, y2);            double d3 = get_distance(real_x, real_y, x3, y3);            double d4 = get_distance(real_x, real_y, x4, y4);            double minD = std::min(d1, std::min(d2, std::min(d3, d4)));            double dens = 0.;            if (minD == d1) // left bottom            {                dens = density[sq_i * NY_1 + sq_j];            }            else if (minD == d2) // right bottom            {                dens = density[(sq_i + 1) * NY_1 + sq_j];            }            else if (minD == d3) // right top            {                dens = density[(sq_i + 1) * NY_1 + sq_j + 1];            }            else if (minD == d4) // left top            {                dens = density[sq_i * (NY_1) + sq_j + 1];            }            phi += mes * dens * jakob;        }    }    if (fabs(phi) < fabs(DBL_MIN_TRIM)) phi = 0;    return phi;}static double numSolAtPoint(double *density, int *gridPr, double x, double y) {    if (x < A || x > B || y < C || y > D) {        return 0.;    }    int lev = 0;    int p = 0;    int q = 0;    do {        double currHx = (B - A) / (NX * std::pow(3., lev));        double currHy = (B - A) / (NX * std::pow(3., lev));        int coef = std::pow(3., R - lev);        p = (int) (x + (currHx / 2.)) / currHx;        q = (int) (y + (currHy / 2.)) / currHy;        int d = gridPr[p * coef * NY3_1 + q * coef];        if (d == lev) { return density[d]; }        lev++;    } while (lev < R);    return 0.;}static void gridVertChecking(int lx, int rx, int by, int uy, int lev, int *grid, double *density);static void gridHorChecking(int lx, int rx, int by, int uy, int lev, int *grid, double *density);static void gridHorRefine(int ii, int jj, int lev, int *grid, double *density) {    int coef = std::pow(3., R_LVL - lev);    int coefNew = R_LVL - (lev + 1);    int mult = std::pow(3., coefNew);    // a) если обе окрестности еще ни разу не мельчились    if (grid[ii * coef * NY3_1 + jj * coef] == lev && grid[(ii + 1) * coef * NY3_1 + jj * coef] == lev) {        // изменим обе окрестности        for (int p = -1; p < 2; ++p) { // p = -1, 0, 1            for (int q = -1; q < 5; ++q) { // q = -1, 0, 1, 2, 3, 4                int i = (3 * ii + p) * NY3_1 * mult + (3 * jj + q) * mult;                grid[i] = lev + 1; // изменим сетку                // вычислим значения                // density[i] !!!!            }        }        // осуществим проверку для вновь созданных ячеек (окрестностей)        gridVertChecking((3 * ii - 1) * mult, (3 * ii + 4) * mult, (3 * jj - 1) * mult, (3 * jj + 1) * mult,                         lev + 1, grid, density);        gridHorChecking((3 * ii - 1) * mult, (3 * ii + 4) * mult, (3 * jj - 1) * mult, (3 * jj + 1) * mult,                        lev + 1, grid, density);        return;    } // если обе окрестности еще ни разу не мельчились    // б) сетки уже адаптированные. Просто переходим на проверку на более низком уровне    if (grid[ii * coef * NY3_1 + jj * coef] != lev && grid[(ii + 1) * coef * NY3_1 + jj * coef] != lev) {        gridVertChecking((3 * ii - 1) * mult, (3 * ii + 4) * mult, (3 * jj - 1) * mult, (3 * jj + 1) * mult,                         lev + 1, grid, density);        gridHorChecking((3 * ii - 1) * mult, (3 * ii + 4) * mult, (3 * jj - 1) * mult, (3 * jj + 1) * mult,                        lev + 1, grid, density);        return;    } // сетки уже адаптированные. Просто переходим на проверку на более низком уровне    // в) одна из сеток не была адаптирована    if (grid[ii * coef * NY3_1 + jj * coef] == lev || grid[(ii + 1) * coef * NY3_1 + jj * coef] == lev) {        int ii_ad = 0;        if (grid[ii * coef * NY3_1 + jj * coef] == lev) {            ii_ad = ii + 1;        }        else {            ii_ad = ii + 1;        }        // теперь мельчим нужную окрестность        for (int p = -1; p < 2; ++p) { // p = -1, 0, 1            for (int q = -1; q < 2; ++q) { // q = -1, 0, 1                grid[(3 * ii_ad + p) * NY3_1 * mult + (3 * jj + q) * mult] = lev + 1;                // вычислим значения                // density[i] !!!!            }        }        // осуществим проверку для вновь созданных ячеек (окрестностей)        gridVertChecking((3 * ii - 1) * mult, (3 * ii + 1) * mult, (3 * jj - 1) * mult, (3 * jj + 4) * mult,                         lev + 1, grid, density);        gridHorChecking((3 * ii - 1) * mult, (3 * ii + 1) * mult, (3 * jj - 1) * mult, (3 * jj + 4) * mult,                        lev + 1, grid, density);        return;    } // одна из сеток не была адаптирована}static void gridVertRefine(int ii, int jj, int lev, int *grid, double *density) {    int coef = std::pow(3., R_LVL - lev);    int coefNew = R_LVL - (lev + 1);    int mult = std::pow(3., coefNew);    // a) если обе окрестности еще ни разу не мельчились    if (grid[ii * coef * NY3_1 + jj * coef] == lev && grid[ii * coef * NY3_1 + (jj + 1) * coef] == lev) {        // изменим обе окрестности        for (int p = -1; p < 2; ++p) { // p = -1, 0, 1            for (int q = -1; q < 5; ++q) { // q = -1, 0, 1, 2, 3, 4                int i = (3 * ii + p) * NY3_1 * mult + (3 * jj + q) * mult;                grid[i] = lev + 1; // изменим сетку                // вычислим значения                // density[i] !!!!            }        }        // осуществим проверку для вновь созданных ячеек (окрестностей)        gridVertChecking((3 * ii - 1) * mult, (3 * ii + 1) * mult, (3 * jj - 1) * mult, (3 * jj + 4) * mult,                         lev + 1, grid, density);        gridHorChecking((3 * ii - 1) * mult, (3 * ii + 1) * mult, (3 * jj - 1) * mult, (3 * jj + 4) * mult,                        lev + 1, grid, density);        return;    } // если обе окрестности еще ни разу не мельчились    // б) сетки уже адаптированные. Просто переходим на проверку на более низком уровне    if (grid[ii * coef * NY3_1 + jj * coef] != lev && grid[ii * coef * NY3_1 + (jj + 1) * coef] != lev) {        gridVertChecking((3 * ii - 1) * mult, (3 * ii + 1) * mult, (3 * jj - 1) * mult, (3 * jj + 4) * mult,                         lev + 1, grid, density);        gridHorChecking((3 * ii - 1) * mult, (3 * ii + 1) * mult, (3 * jj - 1) * mult, (3 * jj + 4) * mult,                        lev + 1, grid, density);        return;    } // сетки уже адаптированные. Просто переходим на проверку на более низком уровне    // в) одна из сеток не была адаптирована    if (grid[ii * coef * NY3_1 + jj * coef] == lev || grid[ii * coef * NY3_1 + (jj + 1) * coef] == lev) {        int jj_ad = 0;        if (grid[ii * coef * NY3_1 + jj * coef] == lev) {            jj_ad = jj + 1;        }        else {            jj_ad = jj + 1;        }        // теперь мельчим нужную окрестность        for (int p = -1; p < 2; ++p) { // p = -1, 0, 1            for (int q = -1; q < 2; ++q) { // q = -1, 0, 1                grid[(ii + p) * NY3_1 * mult + (jj_ad + q) * mult] = lev + 1;                // вычислим значения                // density[i] !!!!            }        }        // осуществим проверку для вновь созданных ячеек (окрестностей)        gridVertChecking((3 * ii - 1) * mult, (3 * ii + 1) * mult, (3 * jj - 1) * mult, (3 * jj + 4) * mult,                         lev + 1, grid, density);        gridHorChecking((3 * ii - 1) * mult, (3 * ii + 1) * mult, (3 * jj - 1) * mult, (3 * jj + 4) * mult,                        lev + 1, grid, density);        return;    } // одна из сеток не была адаптирована} // конец gridVertRefinestatic void gridVertChecking(int lx, int rx, int by, int uy, int lev, int *grid, double *density) {    if (lev < R_LVL) {        for (int i = lx; i < rx + 1; ++i) {            for (int j = by; j < uy; ++j) {                int coef = std::pow(3., R_LVL - lev);                if (std::abs(density[i * coef * NY3_1 + coef] - density[i * coef * NY3_1 + (j + 1) * coef]) >                    RES_EPS) //                {                    gridVertRefine(i, j, lev, grid, density);                }            }        }    }}static void gridHorChecking(int lx, int rx, int by, int uy, int lev, int *grid, double *density) {    if (lev < R_LVL) {        for (int j = by; j < uy + 1; ++j) {            for (int i = lx; i < rx; ++i) {                int coef = std::pow(3., R_LVL - lev);                if (std::abs(density[i * coef * NY3_1 + coef] - density[(i + 1) * coef * NY3_1 + j * coef]) >                    RES_EPS) //                {                    gridHorRefine(i, j, lev, grid, density);                }            }        }    }}double *solve_9(double &tme, int *grid, int *gridPr) {    StartTimer();    fflush(stdout);    double *phi = new double[XY];    double *prev_density = new double[XY];    double *density = new double[XY];    double *residual = new double[XY];    //<editor-fold desc="Init data">    for (int i = 0; i < NX_1; ++i) {        for (int j = 0; j < NY_1; ++j) {            density[NY_1 * i * R + j * R] = -1.;            prev_density[NY_1 * i * R + j * R] = -1.;            residual[NY_1 * i * R + j * R] = 0.;            phi[NY_1 * i * R + j * R] = 0.;        }    }    // G1 -- (x_i, 0=C) -- bottom boundary    for (int i = 0; i < NX_1; ++i) {        prev_density[NY_1 * i] = analytical_solution_circle(0., A + HX * i, C);        if (fabs(prev_density[NY_1 * i]) < fabs(DBL_MIN_TRIM)) prev_density[NY_1 * i] = 0;    }    // G2 -- (NX=B, y_j) -- right boundary    for (int j = 1; j < NY; ++j) {        prev_density[NY_1 * NX + j] = analytical_solution_circle(0., A + HX * NX, C + HY * j);        if (fabs(prev_density[NY_1 * NX + j]) < fabs(DBL_MIN_TRIM)) prev_density[NY_1 * NX + j] = 0;    }    // G3 -- (x_i, NY=D) -- top boundary    for (int i = 0; i < NX_1; ++i) {        prev_density[NY_1 * i + NY] = analytical_solution_circle(0., A + HX * i, C + HY * NY);        if (fabs(prev_density[NY_1 * i + NY]) < fabs(DBL_MIN_TRIM)) prev_density[NY_1 * i + NY] = 0;    }    // G4 -- (0=A, y_j) -- left boundary    for (int j = 1; j < NY; ++j) {        prev_density[j] = analytical_solution_circle(0., A, C + HY * j);        if (fabs(prev_density[j]) < fabs(DBL_MIN_TRIM)) prev_density[j] = 0;    }    memcpy(density, prev_density, XY * sizeof(double));    // inner points    for (int i = 1; i < NX; ++i) {        for (int j = 1; j < NY; ++j) {            //prev_density[NY_1 * i + j] = analytical_solution_circle(0., A + HX * i, C + HY * j);            //if (fabs(prev_density[NY_1 * i + j]) < fabs(DBL_MIN_TRIM)) prev_density[NY_1 * i + j] = 0;        }    }    //</editor-fold>    printf("SUM RHO INIT = %le\n", calc_array_sum(prev_density, NX_1, NY_1, 0));    printf("SUM ABS RHO INIT= %le\n", calc_array_sum(prev_density, NX_1, NY_1, 1));    fflush(stdout);    double *extrems;    double *extrems_err;    double *err;    double l1_err_vec;    double l1_err_tr;    for (int tl = 1; tl <= TIME_STEP_CNT; tl++) {        memcpy(gridPr, grid, XY * sizeof(double));        //<editor-fold desc="calculate phi">        for (int i = 1; i < NX; ++i) // G1 -- (x_i, 0=C) -- bottom boundary            if (G1[i] == 1) phi[NY_1 * i] = calc_phi(i, 0, prev_density, TAU * tl, HX, HY);        for (int j = 1; j < NY; ++j) // G2 -- (NX=B, y_j) -- right boundary            if (G2[j] == 1) phi[NY_1 * NX + j] = calc_phi(NX, j, prev_density, TAU * tl, HX, HY);        for (int i = 1; i < NX; ++i) // G3 -- (x_i, NY=D) -- top boundary            if (G3[i] == 1) phi[NY_1 * i + NY] = calc_phi(i, NY, prev_density, TAU * tl, HX, HY);        for (int j = 1; j < NY; ++j) // G4 -- (0=A, y_j) -- left boundary            if (G4[j] == 1) phi[j] = calc_phi(0, j, prev_density, TAU * tl, HX, HY);        if (CP00 == 1) phi[0] = calc_phi(0, 0, prev_density, TAU * tl, HX, HY); // point (0.0)        if (CP10 == 1) phi[NY_1 * NX] = calc_phi(NX, 0, prev_density, TAU * tl, HX, HY); // point (1.0)        if (CP01 == 1) phi[NY] = calc_phi(0, NY, prev_density, TAU * tl, HX, HY); // point (0.1)        if (CP11 == 1) phi[NY_1 * NX + NY] = calc_phi(NX, NY, prev_density, TAU * tl, HX, HY); // point (1,1)        for (int i = 1; i < NX; ++i) // inner points            for (int j = 1; j < NY; ++j) phi[NY_1 * i + j] = calc_phi(i, j, prev_density, TAU * tl, HX, HY);        //</editor-fold>        double rpCoef = 2. / (HX * HY);        for (int i = 1; i < NX; ++i) // G1 -- (x_i, 0=C) -- bottom boundary            if (G1[i] == 1) {                density[NY_1 * i] = rpCoef * phi[NY_1 * i];                if (fabs(density[NY_1 * i]) < fabs(DBL_MIN_TRIM)) density[NY_1 * i] = 0.;            }        for (int j = 1; j < NY; ++j) // G2 -- (NX=B, y_j) -- right boundary            if (G2[j] == 1) {                density[NY_1 * NX + j] = rpCoef * phi[NY_1 * NX + j];                if (fabs(density[NY_1 * NX + j]) < fabs(DBL_MIN_TRIM)) density[NY_1 * NX + j] = 0.;            }        for (int i = 1; i < NX; ++i) // G3 -- (x_i, NY=D) -- top boundary            if (G3[i] == 1) {                density[NY_1 * i + NY] = rpCoef * phi[NY_1 * i + NY];                if (fabs(density[NY_1 * i + NY]) < fabs(DBL_MIN_TRIM)) density[NY_1 * i + NY] = 0.;            }        for (int j = 1; j < NY; ++j) // G4 -- (0=A, y_j) -- left boundary            if (G4[j] == 1) {                density[j] = rpCoef * phi[j];                if (fabs(density[j]) < fabs(DBL_MIN_TRIM)) density[j] = 0.;            }        rpCoef = 4. / (HX * HY);        if (CP00 == 1) { // point (0,0)            density[0] = phi[0];            if (fabs(density[0]) < fabs(DBL_MIN_TRIM)) density[0] = 0.;        }        if (CP10 == 1) { // point (1,0)            density[NY_1 * NX] = rpCoef * phi[NY_1 * NX];            if (fabs(density[NY_1 * NX]) < fabs(DBL_MIN_TRIM)) density[NY_1 * NX] = 0.;        }        if (CP01 == 1) { // point (0,1)            density[NY] = rpCoef * phi[NY];            if (fabs(density[NY]) < fabs(DBL_MIN_TRIM)) density[NY] = 0.;        }        if (CP11 == 1) { // point (1,1)            density[NY_1 * NX + NY] = rpCoef * phi[NY_1 * NX + NY];            if (fabs(density[NY_1 * NX + NY]) < fabs(DBL_MIN_TRIM)) density[NY_1 * NX + NY] = 0.;        }        rpCoef = 1. / (HX * HY);        for (int i = 1; i < NX; ++i)            for (int j = 1; j < NY; ++j) {                density[NY_1 * i + j] = rpCoef * phi[NY_1 * i + j];                if (fabs(density[NY_1 * i + j]) < fabs(DBL_MIN_TRIM)) density[NY_1 * i + j] = 0.;            }        memcpy(prev_density, density, XY * sizeof(double));        gridVertChecking(0, NX, 0, NY, 0, grid, density);        gridHorChecking(0, NX, 0, NY, 0, grid, density);        if (tl % 5 == 0) {            err = calc_error_6(HX, HY, TAU * tl, density);            l1_err_vec = get_l1_norm_vec(NX_1, NY_1, err);            l1_err_tr = get_l1_norm_int_trapezoidal(HX, HY, NX, NY, err); // note! a loop boundary            extrems = calc_array_extrems(density, NX_1, NY_1);            extrems_err = calc_array_extrems(err, NX_1, NY_1);            printf("tl = %d Sum(Rho)= %le  ERR_VEC= %le  ERR_TR= %le  MAX_RHO= %le"                           "  MAX_ERR= %le\n", tl, calc_array_sum(density, NX_1, NY_1, 0),                   l1_err_vec, l1_err_tr, extrems[1], extrems_err[1]);            fflush(stdout);            print_data_to_files(phi, density, residual, tl);            int fixed_x = (int) (get_center_x() / HX);            int fixed_y = (int) (get_center_y() / HY);            //print_line_along_x("rho", NX, NY, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,            //                  U, V, density, fixed_y);            //print_line_along_y("rho", NX, NY, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,            //                   U, V, density, fixed_x);        }    } // time loop    err = calc_error_6(HX, HY, TAU * TIME_STEP_CNT, density);    l1_err_vec = get_l1_norm_vec(NX_1, NY_1, err);    l1_err_tr = get_l1_norm_int_trapezoidal(HX, HY, NX, NY, err); // note! a loop boundary    extrems = calc_array_extrems(density, NX_1, NY_1);    extrems_err = calc_array_extrems(err, NX_1, NY_1);    append_statistics_expl(NX_1, NY_1, TAU, l1_err_vec, l1_err_tr, extrems,                           extrems_err, TIME_STEP_CNT);    delete[] prev_density;    delete[] phi;    delete[] err;    delete[] residual;    delete[] extrems;    delete[] extrems_err;    tme = GetTimer() / 1000;    return density;}double *calc_error_9(double hx, double hy, double tt, double *solution) {    double *res = new double[XY];    for (int i = 0; i < NX_1; i++)        for (int j = 0; j < NY_1; j++)            res[i * NY_1 + j] = fabs(solution[i * NY_1 + j]                                     - analytical_solution_circle(tt, A + hx * i, C + hy * j));    return res;}double *get_exact_solution_9(double hx, double hy, double t) {    double *res = new double[XY];    for (int i = 0; i < NX_1; i++)        for (int j = 0; j < NY_1; j++)            res[i * NY_1 + j] = fabs(analytical_solution_circle(t, A + hx * i, C + hy * j));    return res;}