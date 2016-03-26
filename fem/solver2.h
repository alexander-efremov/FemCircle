#ifndef FEM_CIRCLE_SOLVER_2_H
#define FEM_CIRCLE_SOLVER_2_H

double *solve_2(double &tme);

double *calc_error_2(double hx, double hy, double tt, double *solution);

double *get_exact_solution_2(double hx, double hy, double t);

#endif //FEM_CIRCLE_SOLVER_2_H