#include <iostream>
#include "tecplot.h"

void print_surface_as_v1(const char *filename, int ox_len, int oy_len,
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

bool print_surface_as_v(const char *filename, int ox_len, int oy_len,
                        double hx, double hy, int t, double a, double c, double x0, double y0, double *data) {
    char name[80];
    sprintf(name, "%s_nx=%d_ny=%d_hx=%f_hy=%f_t=%d_x0=%f_y0=%f.dat",
            filename, ox_len + 1, oy_len + 1, hx, hy, t, x0, y0);
    print_surface_as_v1(name, ox_len, oy_len, hx, hy, data);
    return true;
}

