#ifndef FEM_CIRCLE_TECPLOT_H
#define FEM_CIRCLE_TECPLOT_H

void print_surface(const char *filename, int ox_len, int oy_len, double hx, double hy,
                   int t, double a, double c, double x0, double y0, double tau, double u, double v, double *data);

void print_line_along_x(const char *filename, int ox_len, int oy_len, double hx, double hy,
                        int t, double a, double c, double x0, double y0, double tau, double u, double v, double *data,
                        int fixed_y);

void print_line_along_y(const char *filename, int ox_len, int oy_len, double hx, double hy,
                        int t, double a, double c, double x0, double y0, double tau, double u, double v, double *data,
                        int fixed_x);

#endif //FEM_CIRCLE_TECPLOT_H
