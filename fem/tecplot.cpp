#include <iostream>
#include "tecplot.h"

void print_surface_internal(const char *filename, int ox_len, int oy_len,
                            double hx, double hy, double *data) {
    FILE *file = fopen(filename, "w");
    fprintf(file, "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = 'X' 'Y' 'E'\nZONE T='SubZone'");
    fprintf(file, "\nI=%d J=%d K=%d ZONETYPE=Ordered", oy_len + 1, ox_len + 1, 1);
    fprintf(file, "\nDATAPACKING=POINT\nDT=(SINGLE SINGLE SINGLE)");
    for (int i = 0; i < ox_len + 1; i++)
        for (int j = 0; j < oy_len + 1; j++)
            fprintf(file, "\n%-30.20g  %-30.20g %-30.20g", i * hx, j * hy,
                    data[(oy_len + 1) * i + j]);

    fclose(file);
}

void print_line_by_x_internal(const char *filename, int ox, int oy, int fixed_y,
                              double hx, double hy, double *data) {
    FILE *file = fopen(filename, "w");
    fprintf(file, "TITLE = \"XY LINE\"\nVARIABLES = \"X\", \"F\"\nZONE T=\"Only Zone\",");
    fprintf(file, " I=%d, F=POINT", ox + 1);
    for (int i = 0; i < ox + 1; i++)
        fprintf(file, "\n%-30.20g  %-30.20g", i * hx, data[(oy + 1) * i + fixed_y]);

    fclose(file);
}

void print_line_by_y_internal(const char *filename, int oy, int fixed_x,
                              double hx, double hy, double *data) {
    FILE *file = fopen(filename, "w");
    fprintf(file, "TITLE = \"XY LINE\"\nVARIABLES = \"X\", \"F\"\nZONE T=\"Only Zone\",");
    fprintf(file, " I=%d, F=POINT", oy + 1);
    for (int j = 0; j < oy + 1; j++)
        fprintf(file, "\n%-30.20g  %-30.20g", j * hy, data[(oy + 1) * fixed_x + j]);

    fclose(file);
}

void print_surface(const char *filename, int ox_len, int oy_len,
                   double hx, double hy, int t, double a, double c, double x0, double y0,
                   double tau, double u, double v, double *data) {
    char name[650];
    sprintf(name, "%s_nx=%d_ny=%d_hx=%f_hy=%f_t=%d_x0=%f_y0=%f_tau=%f_u=%f_v=%f_a=%f_c=%f.dat",
            filename, ox_len + 1, oy_len + 1, hx, hy, t, x0, y0, tau, u, v, a, c);
    print_surface_internal(name, ox_len, oy_len, hx, hy, data);
}

void print_line_along_x(const char *filename, int ox_len, int oy_len,
                        double hx, double hy, int t, double a, double c, double x0, double y0,
                        double tau, double u, double v, double *data, int fixed_x) {
    char name[650];
    sprintf(name, "line_by_x_%s_nx=%d_ny=%d_hx=%f_hy=%f_t=%d_x0=%f_y0=%f_tau=%f_u=%f_v=%f_a=%f_c=%f.dat",
            filename, ox_len + 1, oy_len + 1, hx, hy, t, x0, y0, tau, u, v, a, c);
    print_line_by_x_internal(name, ox_len, oy_len, fixed_x, hx, hy, data);
}

void print_line_along_y(const char *filename, int ox_len, int oy_len,
                        double hx, double hy, int t, double a, double c, double x0, double y0,
                        double tau, double u, double v, double *data, int fixed_y) {
    char name[650];
    sprintf(name, "line_by_y_%s_nx=%d_ny=%d_hx=%f_hy=%f_t=%d_x0=%f_y0=%f_tau=%f_u=%f_v=%f_a=%f_c=%f.dat",
            filename, ox_len + 1, oy_len + 1, hx, hy, t, x0, y0, tau, u, v, a, c);
    print_line_by_y_internal(name, oy_len, fixed_y, hx, hy, data);
}

