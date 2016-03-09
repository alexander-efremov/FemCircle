#ifndef FEM_CIRCLE_SOLVER_3_H
#define FEM_CIRCLE_SOLVER_3_H

double *solve_3(double &tme);
double *calc_error_3(double hx, double hy, double *solution);
double get_center_x_3();
double get_center_y_3();
double* get_exact_solution_3(double hx, double hy, double t);

#endif //FEM_CIRCLE_SOLVER_3_H
