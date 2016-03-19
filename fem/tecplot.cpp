#include <iostream>
#include "tecplot.h"

void print_surface_internal(const char *filename, int ox_len, int oy_len,
                            double hx, double hy, double *data) {
    FILE *pfile = fopen(filename, "w");
    fprintf(pfile, "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = 'X' 'Y' 'E'\nZONE T='SubZone'");
    fprintf(pfile, "\nI=%d J=%d K=%d ZONETYPE=Ordered", oy_len + 1, ox_len + 1, 1);
    fprintf(pfile, "\nDATAPACKING=POINT\nDT=(SINGLE SINGLE SINGLE)");
    for (int i = 0; i < ox_len + 1; i++)
        for (int j = 0; j < oy_len + 1; j++)
            fprintf(pfile, "\n%-30.20g  %-30.20g %-30.20g", i * hx, j * hy,
                    data[(oy_len + 1) * i + j]);

    fclose(pfile);
}

void print_line_internal(const char *filename, int ox_len, int oy_len, int y,
                            double hx, double hy, double *data) {
    FILE *pfile = fopen(filename, "w");
    fprintf(pfile, "TITLE = \"XY LINE\"\nVARIABLES = \"X\", \"Y\"\nZONE T=\"Only Zone\",");
    fprintf(pfile, " I=%d, F=POINT", ox_len + 1);
    for (int i = 0; i < ox_len + 1; i++)
       // for (int j = 0; j < oy_len + 1; j++)
            fprintf(pfile, "\n%-30.20g  %-30.20g %-30.20g", i * hx, y,
                    data[(oy_len + 1) * i + y]);

    fclose(pfile);
}

bool print_surface(const char *filename, int ox_len, int oy_len,
                   double hx, double hy, int t, double a, double c, double x0, double y0,
                   double tau, double u, double v, double *data) {
    char name[550];
    sprintf(name, "%s_nx=%d_ny=%d_hx=%f_hy=%f_t=%d_x0=%f_y0=%f_tau=%f_u=%f_v=%f.dat",
            filename, ox_len + 1, oy_len + 1, hx, hy, t, x0, y0, tau, u, v);
    print_surface_internal(name, ox_len, oy_len, hx, hy, data);
    return true;
}

bool print_line(const char *filename, int ox_len, int oy_len,
                   double hx, double hy, int t, double a, double c, double x0, double y0,
                   double tau, double u, double v, double *data) {
    char name[550];
    sprintf(name, "line_%s_nx=%d_ny=%d_hx=%f_hy=%f_t=%d_x0=%f_y0=%f_tau=%f_u=%f_v=%f.dat",
            filename, ox_len + 1, oy_len + 1, hx, hy, t, x0, y0, tau, u, v);
    print_line_internal(name, ox_len, oy_len, 1, hx, hy, data);
    return true;
}

