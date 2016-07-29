#ifndef FEM_CIRCLE_COMMON_H
#define FEM_CIRCLE_COMMON_H

#include "consts.h"

inline double get_center_x() {
    return A + CENTER_OFFSET_X;
}

inline double get_center_y() {
    return C + CENTER_OFFSET_Y;
}

inline void get_left_area(int ii, int jj, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3, double &x4,
                   double &y4, double hy, double hx_add, double hy_add, double a, double b, double c,
                   double d) {
    // p1 (A, y_{j-1/2})
    x1 = a;
    y1 = c + jj * hy - hy_add / 2.;
    // p2 (x_{1/2}, y_{j-1/2})
    x2 = a + hx_add / 2.;
    y2 = c + jj * hy - hy_add / 2.;
    //p3 (x_{1/2}, y_{j+1/2})
    x3 = a + hx_add / 2.;
    y3 = c + jj * hy + hy_add / 2.;
    //p4 (A, y_{j+1/2})
    x4 = a;
    y4 = c + jj * hy + hy_add / 2.;
    if (x1 < a || x1 >= b || x2 <= a || x2 >= b || x3 <= a || x3 >= b || x4 < a || x4 >= b
        || y1 <= c || y1 >= d || y2 <= c || y2 >= d || y3 <= c || y3 >= d || y4 <= c || y4 >= d) {
        printf("9. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
}

inline void get_top_area(int ii, int jj, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3, double &x4,
                  double &y4, double hx, double hx_add, double hy_add, double a, double b, double c,
                  double d) {
    // p1 (x_{i-1/2}, y_{NY-1/2})
    x1 = a + ii * hx - hx_add / 2.;
    y1 = d - hy_add / 2.;
    // p2 (x_{i+1/2}, y_{NY-1/2})
    x2 = a + ii * hx + hx_add / 2.;
    y2 = d - hy_add / 2.;
    //p3 (x_{i+1/2}, D)
    x3 = a + ii * hx + hx_add / 2.;
    y3 = d;
    //p4 (x_{i-1/2}, D)
    x4 = a + ii * hx - hx_add / 2.;
    y4 = d;
    if (x1 <= a || x1 >= b || x2 <= a || x2 >= b || x3 <= a || x3 >= b || x4 <= a || x4 >= b
        || y1 <= c || y1 >= d || y2 <= c || y2 >= d || y3 <= c || y3 > d || y4 <= c || y4 > d) {
        printf("8. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
}

inline void get_right_area(int ii, int jj, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3,
                    double &x4, double &y4, double hy, double hx_add, double hy_add, double a, double b, double c,
                    double d) {
    // p1 (x_{NX-1/2}, y_{j-1/2})
    x1 = b - hx_add / 2.;
    y1 = c + jj * hy - hy_add / 2.;
    // p2 (B, y_{j-1/2})
    x2 = b;
    y2 = c + jj * hy - hy_add / 2.;
    // p3 (B, y_{j+1/2})
    x3 = b;
    y3 = c + jj * hy + hy_add / 2.;
    // p4 (x_{NX-1/2}, y_{j+1/2})
    x4 = b - hx_add / 2.;
    y4 = c + jj * hy + hy_add / 2.;
    if (x1 <= a || x1 >= b || x2 <= a || x2 > b || x3 <= a || x3 > b || x4 <= a || x4 >= b
        || y1 <= c || y1 >= d || y2 <= c || y2 >= d || y3 <= c || y3 >= d || y4 <= c || y4 >= d) {
        printf("7. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
}

inline void get_bottom_area(int ii, int jj, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3, double &x4,
                     double &y4, double hx, double hx_add, double hy_add, double a, double b, double c,
                     double d) {
    // p1 (x_{i-1/2}, C)
    x1 = a + ii * hx - hx_add / 2.;
    y1 = c;
    // p2 (x_{i+1/2}, C)
    x2 = a + ii * hx + hx_add / 2.;
    y2 = c;
    // p3 (x_{i+1/2}, y_{1/2})
    x3 = a + ii * hx + hx_add / 2.;
    y3 = c + hy_add / 2.;
    // p4 (x_{i-1/2}, y_{1/2})
    x4 = a + ii * hx - hx_add / 2.;
    y4 = c + hy_add / 2.;
    if (x1 <= a || x1 >= b || x2 <= a || x2 >= b || x3 <= a || x3 >= b || x4 <= a || x4 >= b
        || y1 < c || y1 >= d || y2 < c || y2 >= d || y3 <= c || y3 >= d || y4 <= c || y4 >= d) {
        printf("6. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
}

inline void get_point_11_area(int ii, int jj, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3,
                       double &x4, double &y4, double hx_add, double hy_add, double a, double b, double c,
                       double d) {
    // p1 (x_{NX-1/2}, y_{NY-1/2})
    x1 = b - hx_add / 2.;
    y1 = d - hy_add / 2.;
    // p2 (B, y_{NY-1/2})
    x2 = b;
    y2 = d - hy_add / 2.;
    // p3 (B, D)
    x3 = b;
    y3 = d;
    // p4 (x_{NX-1/2}, D)
    x4 = b - hx_add / 2.;
    y4 = d;
    if (x1 <= a || x1 >= b || x2 <= a || x2 > b || x3 <= a || x3 > b || x4 <= a || x4 >= b
        || y1 <= c || y1 >= d || y2 <= c || y2 >= d || y3 <= c || y3 > d || y4 <= c || y4 > d) {
        printf("5. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
}

inline void get_point_00_area(int ii, int jj, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3,
                       double &x4, double &y4, double hx_add, double hy_add, double a, double b, double c,
                       double d) {
    // p1 (A, C)
    x1 = a;
    y1 = c;
    // p2 (x_{1/2}, C)
    x2 = a + hx_add / 2.;
    y2 = c;
    // p3 (x_{1/2}, y_{1/2})
    x3 = a + hx_add / 2.;
    y3 = c + hy_add / 2.;
    // p4 (A, y_{1/2})
    x4 = a;
    y4 = c + hy_add / 2.;
    if (x1 < a || x1 >= b || x2 <= a || x2 >= b || x3 <= a || x3 >= b || x4 < a || x4 >= b
        || y1 < c || y1 >= d || y2 < c || y2 >= d || y3 <= c || y3 >= d || y4 <= c || y4 >= d) {
        printf("4. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
}

inline void get_point_01_area(int ii, int jj, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3,
                       double &x4, double &y4, double hx_add, double hy_add, double a, double b, double c,
                       double d) {
    // p1 (A, y_{NY-1/2})
    x1 = a;
    y1 = d - hy_add / 2.;
    // p2 (x_{1/2}, y_{NY-1/2})
    x2 = a + hx_add / 2.;
    y2 = d - hy_add / 2.;
    // p3 (x_{1/2}, D)
    x3 = a + hx_add / 2.;
    y3 = d;
    // p4 (A, D)
    x4 = a;
    y4 = d;
    if (x1 < a || x1 >= b || x2 <= a || x2 >= b || x3 <= a || x3 >= b || x4 < a || x4 >= b
        || y1 <= c || y1 >= d || y2 <= c || y2 >= d || y3 <= c || y3 > d || y4 <= c || y4 > d) {
        printf("3. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
}

inline void get_point_10_area(int ii, int jj, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3,
                       double &x4, double &y4, double hx_add, double hy_add, double a, double b, double c,
                       double d) {
    // p1 (x_{NX-1/2}, C)
    x1 = b - hx_add / 2.;
    y1 = c;
    // p2 (B, C)
    x2 = b;
    y2 = c;
    // p3 (B, y_{1/2})
    x3 = b;
    y3 = c + hy_add / 2.;
    // p4 (x_{NX-1/2}, y_{1/2})
    x4 = b - hx_add / 2.;
    y4 = c + hy_add / 2.;
    if (x1 <= a || x1 >= b || x2 <= a || x2 > b || x3 <= a || x3 > b || x4 <= a || x4 >= b
        || y1 < c || y1 >= d || y2 < c || y2 >= d || y3 <= c || y3 >= d || y4 <= c || y4 >= d)
        printf("2. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
}

inline void get_inner_area(int ii, int jj, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3, double &x4,
                    double &y4, double hx, double hy, double hx_add, double hy_add, double a, double b, double c,
                    double d) {
    // p1 (x_{i-1/2}, y_{j-1/2})
    x1 = a + ii * hx - hx_add / 2.; // a + ii * hx_small - hx_lev / 2.
    y1 = c + jj * hy - hy_add / 2.;
    // p2 (x_{i+1/2}, y_{j-1/2})
    x2 = a + ii * hx + hx_add / 2.;
    y2 = c + jj * hy - hy_add / 2.;
    // p3 (x_{i+1/2}, y_{j+1/2})
    x3 = a + ii * hx + hx_add / 2.;
    y3 = c + jj * hy + hy_add / 2.;
    // p4 (x_{i-1/2}, y_{j+1/2})
    x4 = a + ii * hx - hx_add / 2.;
    y4 = c + jj * hy + hy_add / 2.;
    if (x1 <= a || x1 >= b || x2 <= a || x2 >= b || x3 <= a || x3 >= b || x4 <= a || x4 >= b
        || y1 <= c || y1 >= d || y2 <= c || y2 >= d || y3 <= c || y3 >= d || y4 <= c || y4 >= d) {
        printf("1. Inner point, ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
}

inline void get_inner_area(double x, double y, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3, double &x4,
                    double &y4, double hx_add, double hy_add, double a, double b, double c,
                    double d) {
    // p1 (x_{i-1/2}, y_{j-1/2})
    x1 = x - hx_add / 2.; // a + ii * hx_small - hx_lev / 2.
    y1 = y - hy_add / 2.;
    // p2 (x_{i+1/2}, y_{j-1/2})
    x2 = x + hx_add / 2.;
    y2 = y - hy_add / 2.;
    // p3 (x_{i+1/2}, y_{j+1/2})
    x3 = x + hx_add / 2.;
    y3 = y + hy_add / 2.;
    // p4 (x_{i-1/2}, y_{j+1/2})
    x4 = x - hx_add / 2.;
    y4 = y + hy_add / 2.;
    if (x1 <= a || x1 >= b || x2 <= a || x2 >= b || x3 <= a || x3 >= b || x4 <= a || x4 >= b
        || y1 <= c || y1 >= d || y2 <= c || y2 >= d || y3 <= c || y3 >= d || y4 <= c || y4 >= d) {
        printf("1. Inner point, ERROR INDEX x=%.8le y=%.8le : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", x, y, x1, y1, x2, y2, x3, y3, x4, y4);
    }
}

inline void get_coordinates_on_curr(
        int ii, int jj,
        double &x1, double &y1,
        double &x2, double &y2,
        double &x3, double &y3,
        double &x4, double &y4, double hx, double hy, double hx_add, double hy_add, int nx, int ny, double a, double b,
        double c, double d
);

inline void get_coordinates_on_curr(
        int ii, int jj,
        double &x1, double &y1,
        double &x2, double &y2,
        double &x3, double &y3,
        double &x4, double &y4
) {
    get_coordinates_on_curr(ii, jj, x1, y1, x2, y2, x3, y3, x4, y4, HX, HY, HX, HY, NX, NY, A, B, C, D);
}

inline void get_coordinates_on_curr(
        int ii, int jj,
        double &x1, double &y1,
        double &x2, double &y2,
        double &x3, double &y3,
        double &x4, double &y4, double hx, double hy, double hx_add, double hy_add, int nx, int ny, double a, double b,
        double c, double d
) {
    if (ii > 0 && ii < nx && jj > 0 && jj < ny) {
        get_inner_area(ii, jj, x1, y1, x2, y2, x3, y3, x4, y4, hx, hy, hx_add, hy_add, a, b, c, d);
    }
    else if (ii == nx && jj == 0) { // point (1,0)  omega_{i-1,j}
        get_point_10_area(ii, jj, x1, y1, x2, y2, x3, y3, x4, y4, hx_add, hy_add, a, b, c, d);
    }
    else if (ii == 0 && jj == ny) { // point (0,1)  omega_{i,j-1}
        get_point_01_area(ii, jj, x1, y1, x2, y2, x3, y3, x4, y4, hx_add, hy_add, a, b, c, d);
    }
    else if (ii == 0 && jj == 0) { // point (0,0)  omega_{i,j}
        get_point_00_area(ii, jj, x1, y1, x2, y2, x3, y3, x4, y4, hx_add, hy_add, a, b, c, d);
    }
    else if (ii == nx && jj == ny) { // point (1,1)  omega_{i-1,j-1}
        get_point_11_area(ii, jj, x1, y1, x2, y2, x3, y3, x4, y4, hx_add, hy_add, a, b, c, d);
    }
    else if (ii > 0 && ii < nx && jj == 0) { // G1 -- bottom boundary
        get_bottom_area(ii, jj, x1, y1, x2, y2, x3, y3, x4, y4, hx, hx_add, hy_add, a, b, c, d);
    }
    else if (ii == nx && jj > 0 && jj < ny) { // G2 -- right boundary
        get_right_area(ii, jj, x1, y1, x2, y2, x3, y3, x4, y4, hy, hx_add, hy_add, a, b, c, d);
    }
    else if (jj == ny && ii > 0 && ii < nx) { // G3 -- top boundary
        get_top_area(ii, jj, x1, y1, x2, y2, x3, y3, x4, y4, hx, hx_add, hy_add, a, b, c, d);
    }
    else if (ii == 0 && jj > 0 && jj < ny) { // G4 -- left boundary
        get_left_area(ii, jj, x1, y1, x2, y2, x3, y3, x4, y4, hy, hx_add, hy_add, a, b, c, d);
    }
    else {
        printf("ERROR! INDEX i=%d j=%d ", ii, jj);
    }
}

double *solve_1(double &tme);

double *calc_error_1(double hx, double hy, double *solution);

double *solve_2(double &tme);

double *calc_error_2(double hx, double hy, double tt, double *solution);

double *get_exact_solution_2(double hx, double hy, double t);

double *solve_3(double &tme);

double *calc_error_3(double hx, double hy, double tt, double *solution);

double *get_exact_solution_3(double hx, double hy, double t);

double *solve_4(double &tme);

double *calc_error_4(double hx, double hy, double tt, double *solution);

double *get_exact_solution_4(double hx, double hy, double t);

double *solve_5(double &tme);

double *calc_error_5(double hx, double hy, double tt, double *solution);

double *get_exact_solution_5(double hx, double hy, double t);

double *solve_6(double &tme);

double *calc_error_6(double hx, double hy, double tt, double *solution);

double *get_exact_solution_6(double hx, double hy, double t);

double *solve_7(double &tme);

double *calc_error_7(double hx, double hy, double tt, double *solution);

double *get_exact_solution_7(double hx, double hy, double t);

double *solve_8(double &tme);

double *calc_error_8(double hx, double hy, double tt, double *solution);

double *get_exact_solution_8(double hx, double hy, double t);

double *solve_9(double &tme, int *grid, int *gridPr);

double *calc_error_9(double hx, double hy, double tt, double *solution, int nx, int ny);

double *get_exact_solution_9(double hx, double hy, double t, int nx, int ny);

#endif //FEM_CIRCLE_COMMON_H