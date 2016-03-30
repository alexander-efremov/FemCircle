#ifndef FEM_CIRCLE_COMMON_H
#define FEM_CIRCLE_COMMON_H

#include "consts.h"

inline double get_center_x()
{
    return A + CENTER_OFFSET_X;
}

inline double get_center_y()
{
    return C + CENTER_OFFSET_Y;
}

double *solve_1(double &tme);
double *calc_error_1(double hx, double hy, double *solution);

double *solve_2(double &tme);
double *calc_error_2(double hx, double hy, double tt, double *solution);
double *get_exact_solution_2(double hx, double hy, double t);

double *solve_3(double &tme);
double *calc_error_3(double hx, double hy, double tt, double *solution);
double *get_exact_solution_3(double hx, double hy, double t);

#endif //FEM_CIRCLE_COMMON_H
