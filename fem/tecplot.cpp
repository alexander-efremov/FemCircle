#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tecplot.h"

using namespace std;

void print_vector_by_3d(const char *filename, double *mas_0x, int len_x, double *mas_0y, int len_y, double *mas_0z) {
    FILE *pfile;
    pfile = fopen(filename, "w");
    if (pfile == NULL) {
        cout << "\nError can not creat-open file ", filename;
        exit(1);
    }
    fprintf(pfile, "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'");
    fprintf(pfile, "\nVARIABLES = 'X' 'Y' 'E' ");
    fprintf(pfile, "\nZONE T='SubZone'");
    fprintf(pfile, "\nI=%d J=%d K=%d ZONETYPE=Ordered", len_y, len_x, 1);
    fprintf(pfile, "\nDATAPACKING=POINT");
    fprintf(pfile, "\nDT=(SINGLE SINGLE SINGLE )");

    for (int i = 0; i < len_x; i++)
        for (int j = 0; j < len_y; j++)
            fprintf(pfile, "\n%-30.20g  %-30.20g %-30.20g", mas_0x[i], mas_0y[j], mas_0z[j * len_x + i]);

    fclose(pfile);
}

bool print_vector_by_3d(const char *fileName, int n, int numOfCurrTimeIter, int numOfSolOrd, double numOfTimeSteps,
                        double *masOx, double *VectorOfData, int dimOfVect) {
    char strOfName[50];
    char str1[3] = "N=", str2[7];
    char str3[3] = "k=", str4[7];
    char str5[5] = ".dat", str6[2];

    strcpy(strOfName, fileName);
    //itoa( numOfSolOrd, str6, 10);
    strcat(strOfName, str6);
    strcat(strOfName, str1);
    //itoa(n,str2, 10);
    strcat(strOfName, str2);
    strcat(strOfName, str3);
    //itoa( numOfCurrTimeIter, str4, 10);
    strcat(strOfName, str4);
    strcat(strOfName, str5);
    double *mas_Empty = new double[1];
    mas_Empty[0] = ((int) (numOfCurrTimeIter * 1000. / numOfTimeSteps)) / 10;
    print_vector_by_3d(strOfName, masOx, dimOfVect, mas_Empty, 1, VectorOfData);
    delete mas_Empty;
    return true;
}

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
                        double hx, double hy, int t, double a, double c, double *data) {
    char name[80];
    double x0 = a + ox_len * hx / 2.;
    double y0 = c + oy_len * hy / 2.;
    sprintf(name, "%s_nx=%d_ny=%d_hx=%f_hy=%f_t=%d_x0=%f_y0=%f.dat",
            filename, ox_len + 1, oy_len + 1, hx, hy, t, x0, y0);
    print_surface_as_v1(name, ox_len, oy_len, hx, hy, data);
    return true;
}

