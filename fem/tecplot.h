#ifndef FEM_CIRCLE_TECPLOT_H
#define FEM_CIRCLE_TECPLOT_H

bool print_surface(const char *filename, int ox_len, int oy_len, double hx, double hy,
                   int t, double a, double c, double x0, double y0, double tau, double u, double v, double *data);

bool print_line(const char *filename, int ox_len, int oy_len, double hx, double hy,
                   int t, double a, double c, double x0, double y0, double tau, double u, double v, double *data);

#endif //FEM_CIRCLE_TECPLOT_H
