#include <math.h>#include <stdio.h>#include <string.h>#include <float.h>#include "consts.h"#include "timer.h"#include "utils.h"#include "common.h"#include "algorithm"using namespace std;inline void print_data_to_files(double *phi, double *density, int tl, int nx, int ny, double hx, double hy, double u,                                double v, double a, double c, double tau) {    print_surface("phi", nx, ny, hx, hy, tl, a, c, get_center_x(), get_center_y(), tau,                  u, v, phi);    print_surface("rho", nx, ny, hx, hy, tl, a, c, get_center_x(), get_center_y(), tau,                  u, v, density);}inline void print_surface_0(const char *filename, int ox_len, int oy_len, int ny3_1, int r,                          double hx, double hy, int t, double a, double c, double x0, double y0,                          double tau, double u, double v, double *data) {    char name[650];    sprintf(name, "%s_nx=%d_ny=%d_hx=%f_hy=%f_t=%d_x0=%f_y0=%f_tau=%f_u=%f_v=%f_a=%f_c=%f.dat",            filename, ox_len + 1, oy_len + 1, hx, hy, t, x0, y0, tau, u, v, a, c);    FILE *file = fopen(name, "w");    fprintf(file, "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = 'x' 'y' %s\nZONE T='SubZone'",filename);    fprintf(file, "\nI=%d J=%d K=%d ZONETYPE=Ordered", oy_len + 1, ox_len + 1, 1);    fprintf(file, "\nDATAPACKING=POINT\nDT=(SINGLE SINGLE SINGLE)");    for (int i = 0; i < ox_len + 1; i++)        for (int j = 0; j < oy_len + 1; j++)            fprintf(file, "\n%-30.20g  %-30.20g %-30.20g", i * hx, j * hy,                    data[ny3_1 * i * r + j * r]);    fclose(file);}inline static double func_u(double t, double x, double y) { return U; }inline static double func_v(double t, double x, double y) { return V; }inline static double analytical_solution_circle(double t, double x, double y) {    double x0 = get_center_x() + t * func_u(t, x, y);    double y0 = get_center_y() + t * func_v(t, x, y);    double value = (x - x0) * (x - x0) + (y - y0) * (y - y0);    if (value <= R_SQ) return INN_DENSITY;    return OUT_DENSITY;}static inline double get_distance(double real_x, double real_y, double x, double y) {    double sqX = (real_x - x) * (real_x - x);    double sqY = (real_y - y) * (real_y - y);    double r = sqrt(sqX + sqY);    return r;}static double calc_phi(int ii, int jj, double *density, double time_value, double hx, double hy) {    double x1 = 0.;    double y1 = 0.;    double x2 = 0.;    double y2 = 0.;    double x3 = 0.;    double y3 = 0.;    double x4 = 0.;    double y4 = 0.;    get_coordinates_on_curr(ii, jj, x1, y1, x2, y2, x3, y3, x4, y4, hx, hy, NX, NY);//    printf("POINT: %d   %d :  x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le **//                   x4=%.8le * y4=%.8le\n", ii,jj, x1,y1, x2,y2, x3,y3, x4,y4);    double u = func_u(time_value, x1, y1);    double v = func_v(time_value, x1, y1);    x1 = x1 - TAU * u;    y1 = y1 - TAU * v;    u = func_u(time_value, x2, y2);    v = func_v(time_value, x2, y2);    x2 = x2 - TAU * u;    y2 = y2 - TAU * v;    u = func_u(time_value, x3, y3);    v = func_v(time_value, x3, y3);    x3 = x3 - TAU * u;    y3 = y3 - TAU * v;    u = func_u(time_value, x4, y4);    v = func_v(time_value, x4, y4);    x4 = x4 - TAU * u;    y4 = y4 - TAU * v;    /*     if (x1 <= A || x1 >= B || x2 <= A || x2 >= B || x3 <= A || x3 >= B || x4 <= A || x4 >= B        || y1 <= C || y1 >= D || y2 <= C || y2 >= D || y3 <= C || y3 >= D || y4 <= C || y4 >= D)        printf("PREV Time level %.8le! ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "                       "x4=%.8le * y4=%.8le\n ", time_value, ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);*/    int nx = IDEAL_SQ_SIZE_X;    int ny = IDEAL_SQ_SIZE_Y;    double x_step = 1. / nx;    double y_step = 1. / ny;    double phi = 0.;    double mes = x_step * y_step;    for (int i = 0; i < nx; ++i) {        for (int j = 0; j < ny; ++j) {            double ideal_x = i * x_step + x_step / 2.;            double ideal_y = j * y_step + y_step / 2.;            double real_x = x1 + (x2 - x1) * ideal_x + (x4 - x1) * ideal_y                            + (x1 + x3 - x2 - x4) * ideal_x * ideal_y;            double real_y = y1 + (y2 - y1) * ideal_x + (y4 - y1) * ideal_y                            + (y1 + y3 - y2 - y4) * ideal_x * ideal_y;            // find out in which square real point was placed            int sq_i = (int) (((real_x - A + hx/2.)) / hx);            int sq_j = (int) (((real_y - C + hy/2.)) / hy);            double a11 = (x2 - x1) + (x1 + x3 - x2 - x4) * ideal_y;            double a12 = (x4 - x1) + (x1 + x3 - x2 - x4) * ideal_x;            double a21 = (y2 - y1) + (y1 + y3 - y2 - y4) * ideal_y;            double a22 = (y4 - y1) + (y1 + y3 - y2 - y4) * ideal_x;            double jakob = a11 * a22 - a21 * a12;            double dens = density[R * sq_i * NY3_1 + R * sq_j];            phi += mes * dens * jakob;        }    }    if (fabs(phi) < fabs(DBL_MIN_TRIM)) phi = 0.;    return phi;}static double num_sol_at_point(double *prev_rho, int *grid_prev, double x, double y) {    if (x < A || x > B || y < C || y > D) {        return 0.;    }    int lev = 0;    int p = 0;    int q = 0;    do {        double hx = (B - A) / (NX * std::pow(3., lev));        double hy = (D - C) / (NX * std::pow(3., lev));        int coef = std::pow(3., R_LVL - lev);        p = (int) (x + (hx / 2.)) / hx;        q = (int) (y + (hy / 2.)) / hy;        int d = grid_prev[p * coef * NY3_1 + q * coef];        if (d == lev) { return prev_rho[d]; }        lev++;    } while (lev < R_LVL);    return 0.;}static double get_numeric_solution(double x, double y, double hx, double hy) {    return 0.; // !!!!!!!!!!!!!!!!!!!!!!!}static void grid_vert_checking(int lx, int rx, int by, int uy, int lev, int *grid, double *density);static void grid_hor_checking(int lx, int rx, int by, int uy, int lev, int *grid, double *density);static void grid_hor_refine(int ii, int jj, int lev, int *grid, double *rho) {    int coef = std::pow(3., R_LVL - lev);    int coefNew = std::pow(3., R_LVL - lev + 1);    // a) если обе окрестности еще ни разу не мельчились    if (grid[ii * coef * NY3_1 + jj * coef] == lev && grid[(ii + 1) * coef * NY3_1 + jj * coef] == lev) {        // изменим обе окрестности        for (int p = -1; p < 2; ++p) { // p = -1, 0, 1            for (int q = -1; q < 5; ++q) { // q = -1, 0, 1, 2, 3, 4                int i = (3 * ii + p) * NY3_1 * coefNew + (3 * jj + q) * coefNew;                grid[i] = lev + 1; // изменим сетку                double currHx = (B - A) / (NX * std::pow(3., lev + 1));                double currHy = (D - C) / (NX * std::pow(3., lev + 1));                double x = (3 * ii + p) * coefNew * currHx;                double y = (3 * jj + q) * coefNew * currHy;                rho[i] = get_numeric_solution(x, y, currHx, currHy); // вычислим значение            }        }        // осуществим проверку для вновь созданных ячеек (окрестностей)        grid_vert_checking((3 * ii - 1) * coefNew, (3 * ii + 4) * coefNew, (3 * jj - 1) * coefNew, (3 * jj + 1) *                                                                                                   coefNew, lev + 1,                           grid, rho);        grid_hor_checking((3 * ii - 1) * coefNew, (3 * ii + 4) * coefNew, (3 * jj - 1) * coefNew, (3 * jj + 1) *                                                                                                  coefNew, lev + 1,                          grid, rho);        return;    } // если обе окрестности еще ни разу не мельчились    // б) сетки уже адаптированные. Просто переходим на проверку на более низком уровне    if (grid[ii * coef * NY3_1 + jj * coef] != lev && grid[(ii + 1) * coef * NY3_1 + jj * coef] != lev) {        grid_vert_checking((3 * ii - 1) * coefNew, (3 * ii + 4) * coefNew, (3 * jj - 1) * coefNew, (3 * jj + 1) *                                                                                                   coefNew, lev + 1,                           grid, rho);        grid_hor_checking((3 * ii - 1) * coefNew, (3 * ii + 4) * coefNew, (3 * jj - 1) * coefNew, (3 * jj + 1) *                                                                                                  coefNew, lev + 1,                          grid, rho);        return;    } // сетки уже адаптированные. Просто переходим на проверку на более низком уровне    // в) одна из сеток не была адаптирована    if (grid[ii * coef * NY3_1 + jj * coef] == lev || grid[(ii + 1) * coef * NY3_1 + jj * coef] == lev) {        int ii_ad = 0;        if (grid[ii * coef * NY3_1 + jj * coef] == lev) {            ii_ad = ii;        }        else {            ii_ad = ii + 1;        }        // теперь мельчим нужную окрестность        for (int p = -1; p < 2; ++p) { // p = -1, 0, 1            for (int q = -1; q < 2; ++q) { // q = -1, 0, 1                int i = (3 * ii_ad + p) * NY3_1 * coefNew + (3 * jj + q) * coefNew;                grid[i] = lev + 1;                double currHx = (B - A) / (NX * std::pow(3., lev + 1));                double currHy = (D - C) / (NX * std::pow(3., lev + 1));                double x = (3 * ii_ad + p) * coefNew * currHx;                double y = (3 * jj + q) * coefNew * currHy;                rho[i] = get_numeric_solution(x, y, currHx, currHy); // вычислим значение            }        }        // осуществим проверку для вновь созданных ячеек (окрестностей)        grid_vert_checking((3 * ii - 1) * coefNew, (3 * ii + 1) * coefNew, (3 * jj - 1) * coefNew, (3 * jj + 4) *                                                                                                   coefNew,                           lev + 1, grid, rho);        grid_hor_checking((3 * ii - 1) * coefNew, (3 * ii + 1) * coefNew, (3 * jj - 1) * coefNew, (3 * jj + 4) *                                                                                                  coefNew,                          lev + 1, grid, rho);        return;    } // одна из сеток не была адаптирована}static void grid_vert_refine(int ii, int jj, int lev, int *grid, double *density) {    int coef = std::pow(3., R_LVL - lev);    int coefNew = std::pow(3., R_LVL - lev + 1);    // a) если обе окрестности еще ни разу не мельчились    if (grid[ii * coef * NY3_1 + jj * coef] == lev && grid[ii * coef * NY3_1 + (jj + 1) * coef] == lev) {        // изменим обе окрестности        for (int p = -1; p < 2; ++p) { // p = -1, 0, 1            for (int q = -1; q < 5; ++q) { // q = -1, 0, 1, 2, 3, 4                int i = (3 * ii + p) * NY3_1 * coefNew + (3 * jj + q) * coefNew;                grid[i] = lev + 1; // изменим сетку                double currHx = (B - A) / (NX * std::pow(3., lev + 1));                double currHy = (D - C) / (NX * std::pow(3., lev + 1));                double x = (3 * ii + p) * coefNew * currHx;                double y = (3 * jj + q) * coefNew * currHy;                density[i] = get_numeric_solution(x, y, currHx, currHy); // вычислим значение            }        }        // осуществим проверку для вновь созданных ячеек (окрестностей)        grid_vert_checking((3 * ii - 1) * coefNew, (3 * ii + 1) * coefNew, (3 * jj - 1) * coefNew, (3 * jj + 4) *                                                                                                   coefNew,                           lev + 1, grid, density);        grid_hor_checking((3 * ii - 1) * coefNew, (3 * ii + 1) * coefNew, (3 * jj - 1) * coefNew, (3 * jj + 4) *                                                                                                  coefNew,                          lev + 1, grid, density);        return;    } // если обе окрестности еще ни разу не мельчились    // б) сетки уже адаптированные. Просто переходим на проверку на более низком уровне    if (grid[ii * coef * NY3_1 + jj * coef] != lev && grid[ii * coef * NY3_1 + (jj + 1) * coef] != lev) {        grid_vert_checking((3 * ii - 1) * coefNew, (3 * ii + 1) * coefNew, (3 * jj - 1) * coefNew, (3 * jj + 4) *                                                                                                   coefNew,                           lev + 1, grid, density);        grid_hor_checking((3 * ii - 1) * coefNew, (3 * ii + 1) * coefNew, (3 * jj - 1) * coefNew, (3 * jj + 4) *                                                                                                  coefNew,                          lev + 1, grid, density);        return;    } // сетки уже адаптированные. Просто переходим на проверку на более низком уровне    // в) одна из сеток не была адаптирована    if (grid[ii * coef * NY3_1 + jj * coef] == lev || grid[ii * coef * NY3_1 + (jj + 1) * coef] == lev) {        int jj_ad = 0;        if (grid[ii * coef * NY3_1 + jj * coef] == lev) {            jj_ad = jj;        }        else {            jj_ad = jj + 1;        }        // теперь мельчим нужную окрестность        for (int p = -1; p < 2; ++p) { // p = -1, 0, 1            for (int q = -1; q < 2; ++q) { // q = -1, 0, 1                int i = (ii + p) * NY3_1 * coefNew + (jj_ad + q) * coefNew;                grid[i] = lev + 1;                double currHx = (B - A) / (NX * std::pow(3., lev + 1));                double currHy = (D - C) / (NX * std::pow(3., lev + 1));                double x = (ii + p) * coefNew * currHx;                double y = (jj_ad + q) * coefNew * currHy;                density[i] = get_numeric_solution(x, y, currHx, currHy); // вычислим значение            }        }        // осуществим проверку для вновь созданных ячеек (окрестностей)        grid_vert_checking((3 * ii - 1) * coefNew, (3 * ii + 1) * coefNew, (3 * jj - 1) * coefNew, (3 * jj + 4) *                                                                                                   coefNew,                           lev + 1, grid, density);        grid_hor_checking((3 * ii - 1) * coefNew, (3 * ii + 1) * coefNew, (3 * jj - 1) * coefNew, (3 * jj + 4) *                                                                                                  coefNew,                          lev + 1, grid, density);        return;    } // одна из сеток не была адаптирована} // конец grid_vert_refinestatic void grid_vert_checking(int lx, int rx, int by, int uy, int lev, int *grid, double *density) {    if (lev >= R_LVL) return;    int coef = std::pow(3., R_LVL - lev);    for (int i = lx; i < rx + 1; ++i) {        for (int j = by; j < uy; ++j) {            if (std::abs(density[i * coef * NY3_1 + coef] - density[i * coef * NY3_1 + (j + 1) * coef]) >                RES_EPS) {                grid_vert_refine(i, j, lev, grid, density);            }        }    }}static void grid_hor_checking(int lx, int rx, int by, int uy, int lev, int *grid, double *density) {    if (lev >= R_LVL) return;    int coef = std::pow(3., R_LVL - lev);    for (int j = by; j < uy + 1; ++j) {        for (int i = lx; i < rx; ++i) {            if (std::abs(density[i * coef * NY3_1 + coef] - density[(i + 1) * coef * NY3_1 + j * coef]) >                RES_EPS) {                grid_hor_refine(i, j, lev, grid, density);            }        }    }}double *solve_9(double &tme, int *grid, int *gridPr) {    StartTimer();    fflush(stdout);    double *phi = new double[XY];    double *prev_rho = new double[XY];    double *rho = new double[XY];    double *residual = new double[XY];    //<editor-fold desc="Init data">    for (int i = 0; i < XY; ++i) {        grid[i] = -1;        gridPr[i] = -1;        rho[i] = -1.;        prev_rho[i] = -1.;        residual[i] = 0.;        phi[i] = 0.;    }    for (int i = 0; i < NX_1; ++i) {        for (int j = 0; j < NY_1; ++j) {            grid[NY3_1 * i * R + j * R] = 0;            gridPr[NY3_1 * i * R + j * R] = 0;        }    }    for (int i = 0; i < NX_1; ++i) // G1 -- (x_i, 0=C) -- bottom boundary        prev_rho[NY3_1 * i * R] = analytical_solution_circle(0., A + HX * i, C);    for (int j = 1; j < NY; ++j) // G2 -- (NX=B, y_j) -- right boundary        prev_rho[NY3_1 * NX3 + j * R] = analytical_solution_circle(0., A + HX * NX, C + HY * j);    for (int i = 0; i < NX_1; ++i) // G3 -- (x_i, NY=D) -- top boundary        prev_rho[NY3_1 * i * R + NY3] = analytical_solution_circle(0., A + HX * i, C + HY * NY);    for (int j = 1; j < NY; ++j) // G4 -- (0=A, y_j) -- left boundary        prev_rho[j * R] = analytical_solution_circle(0., A, C + HY * j);    memcpy(rho, prev_rho, XY * sizeof(double));    for (int i = 1; i < NX; ++i) // inner points        for (int j = 1; j < NY; ++j)            prev_rho[NY3_1 * i * R + j * R] = analytical_solution_circle(0., A + HX * i, C + HY * j);    print_surface_0("rho", NX, NY, NY3_1, R, HX, HY, 0, A, C, get_center_x(), get_center_y(), TAU, U, V, rho);    print_surface_0("prev_rho", NX, NY, NY3_1, R, HX, HY, 0, A, C, get_center_x(), get_center_y(), TAU, U, V, prev_rho);    //</editor-fold>    printf("SUM RHO INIT = %le\n", calc_array_sum(prev_rho, NX_1, NY_1, 0));    printf("SUM ABS RHO INIT= %le\n", calc_array_sum(prev_rho, NX_1, NY_1, 1));    fflush(stdout);//    double *extrems;//    double *extrems_err;    double *err;    double l1_err_vec;    double l1_err_tr;    for (int tl = 1; tl <= TIME_STEP_CNT; tl++) {        //memcpy(gridPr, grid, XY * sizeof(int));        //<editor-fold desc="calculate phi">        for (int i = 1; i < NX; ++i) // G1 -- (x_i, 0=C) -- bottom boundary            if (G1[i] == 1) phi[NY3_1 * i * R] = calc_phi(i, 0, prev_rho, TAU * tl, HX, HY);        for (int j = 1; j < NY; ++j) // G2 -- (NX=B, y_j) -- right boundary            if (G2[j] == 1) phi[NY3_1 * NX3 + j * R] = calc_phi(NX, j, prev_rho, TAU * tl, HX, HY);        for (int i = 1; i < NX; ++i) // G3 -- (x_i, NY=D) -- top boundary            if (G3[i] == 1) phi[NY3_1 * i * R + NY3] = calc_phi(i, NY, prev_rho, TAU * tl, HX, HY);        for (int j = 1; j < NY; ++j) // G4 -- (0=A, y_j) -- left boundary            if (G4[j] == 1) phi[j * R] = calc_phi(0, j, prev_rho, TAU * tl, HX, HY);        if (CP00 == 1) phi[0] = calc_phi(0, 0, prev_rho, TAU * tl, HX, HY); // point (0.0)        if (CP10 == 1) phi[NY3_1 * NX3] = calc_phi(NX, 0, prev_rho, TAU * tl, HX, HY); // point (1.0)        if (CP01 == 1) phi[NY3] = calc_phi(0, NY, prev_rho, TAU * tl, HX, HY); // point (0.1)        if (CP11 == 1) phi[NY3_1 * NX3 + NY3] = calc_phi(NX, NY, prev_rho, TAU * tl, HX, HY); // point (1,1)        for (int i = 1; i < NX; ++i) // inner points            for (int j = 1; j < NY; ++j)                phi[NY3_1 * i * R + j * R] = calc_phi(i, j, prev_rho, TAU * tl, HX, HY);        //</editor-fold>        double rpCoef = 2. / (HX * HY);        for (int i = 1; i < NX; ++i) // G1 -- (x_i, 0=C) -- bottom boundary            if (G1[i] == 1) {                rho[NY3_1 * i * R] = rpCoef * phi[NY3_1 * i * R];                if (fabs(rho[NY3_1 * i * R]) < fabs(DBL_MIN_TRIM)) rho[NY3_1 * i * R] = 0.;            }        for (int j = 1; j < NY; ++j) // G2 -- (NX=B, y_j) -- right boundary            if (G2[j] == 1) {                rho[NY3_1 * NX3 + j * R] = rpCoef * phi[NY3_1 * NX3 + j * R];                if (fabs(rho[NY3_1 * NX3 + j * R]) < fabs(DBL_MIN_TRIM)) rho[NY3_1 * NX3 + j * R] = 0.;            }        for (int i = 1; i < NX; ++i) // G3 -- (x_i, NY=D) -- top boundary            if (G3[i] == 1) {                rho[NY3_1 * i * R + NY3] = rpCoef * phi[NY3_1 * i * R + NY3];                if (fabs(rho[NY3_1 * i * R + NY3]) < fabs(DBL_MIN_TRIM)) rho[NY3_1 * i * R + NY3] = 0.;            }        for (int j = 1; j < NY; ++j) // G4 -- (0=A, y_j) -- left boundary            if (G4[j] == 1) {                rho[j * R] = rpCoef * phi[j * R];                if (fabs(rho[j * R]) < fabs(DBL_MIN_TRIM)) rho[j * R] = 0.;            }        rpCoef = 4. / (HX * HY);        if (CP00 == 1) { // point (0,0)            rho[0] = rpCoef * phi[0];            if (fabs(rho[0]) < fabs(DBL_MIN_TRIM)) rho[0] = 0.;        }        if (CP10 == 1) { // point (1,0)            rho[NY3_1 * NX3] = rpCoef * phi[NY3_1 * NX3];            if (fabs(rho[NY3_1 * NX3]) < fabs(DBL_MIN_TRIM)) rho[NY3_1 * NX3] = 0.;        }        if (CP01 == 1) { // point (0,1)            rho[NY3] = rpCoef * phi[NY3];            if (fabs(rho[NY3]) < fabs(DBL_MIN_TRIM)) rho[NY3] = 0.;        }        if (CP11 == 1) { // point (1,1)            rho[NY3_1 * NX3 + NY3] = rpCoef * phi[NY3_1 * NX3 + NY3];            if (fabs(rho[NY3_1 * NX3 + NY3]) < fabs(DBL_MIN_TRIM)) rho[NY3_1 * NX3 + NY3] = 0.;        }        rpCoef = 1. / (HX * HY);        for (int i = 1; i < NX; ++i)            for (int j = 1; j < NY; ++j) {                rho[NY3_1 * i * R + j * R] = rpCoef * phi[NY3_1 * i * R + j * R];                if (fabs(rho[NY3_1 * i * R + j * R]) < fabs(DBL_MIN_TRIM)) rho[NY3_1 * i * R + j * R] = 0.;            }        memcpy(prev_rho, rho, XY * sizeof(double));//        grid_vert_checking(0, NX, 0, NY, 0, grid, rho);//        grid_hor_checking(0, NX, 0, NY, 0, grid, rho);        if (tl % 1 == 0) {//            err = calc_error_9(HX, HY, TAU * tl, rho, NX3_1, NY3_1);//            l1_err_vec = get_l1_norm_vec(NX3_1, NY3_1, err);//            l1_err_tr = get_l1_norm_int_trapezoidal(HX, HY, NX3, NY3, err); // note! a loop boundary//            extrems = calc_array_extrems(rho, NX3_1, NY3_1);//            extrems_err = calc_array_extrems(err, NX3_1, NY3_1);//            printf("tl = %d Sum(Rho)= %le  ERR_VEC= %le  ERR_TR= %le  MAX_RHO= %le MAX_ERR= %le\n",//                   tl, calc_array_sum(rho, NX3_1, NY3_1, 0), l1_err_vec, l1_err_tr, extrems[1], extrems_err[1]);            fflush(stdout);            //double *phi, double *density, int tl, int nx, int ny, double hx, double hy, double u, double v, double a, double c, double tau           // print_data_to_files(phi, rho, tl, NX3, NY3, (B - A) / NX3, (D - C) / NY3, U, V, A, C, TAU);            print_surface_0("phi", NX, NY, NY3_1, R, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU, U, V, phi);        }    } // time loop//    err = calc_error_9(HX, HY, TAU * TIME_STEP_CNT, rho);//    l1_err_vec = get_l1_norm_vec(NX3_1, NY3_1, err);//    l1_err_tr = get_l1_norm_int_trapezoidal(HX, HY, NX3, NY3, err); // note! a loop boundary//    extrems = calc_array_extrems(rho, NX3_1, NY3_1);//    extrems_err = calc_array_extrems(err, NX3_1, NY3_1);////    append_statistics_expl(NX_1, NY_1, TAU, l1_err_vec, l1_err_tr, extrems,//                           extrems_err, TIME_STEP_CNT);    delete[] prev_rho;    delete[] phi;   // delete[] err;    delete[] residual;//    delete[] extrems;//    delete[] extrems_err;    tme = GetTimer() / 1000;    return rho;}double *calc_error_9(double hx, double hy, double tt, double *solution, int nx, int ny) {    double *res = new double[XY];    for (int i = 0; i < nx; i++)        for (int j = 0; j < ny; j++)            res[i * ny + j] = fabs(solution[i * ny + j]                                   - analytical_solution_circle(tt, A + hx * i, C + hy * j));    return res;}double *get_exact_solution_9(double hx, double hy, double t, int nx, int ny) {    double *res = new double[nx*ny];    for (int i = 0; i < nx; i++)        for (int j = 0; j < ny; j++)            res[i * ny + j] = fabs(analytical_solution_circle(t, A + hx * i, C + hy * j));    return res;}