#ifndef FEM_CIRCLE_TECPLOT_H
#define FEM_CIRCLE_TECPLOT_H

void print_vector_by_3d(const char *filename, double *mas_0x, int len_x, double *mas_0y, int len_y,
                        double *mas_0z);

bool print_vector_by_3d(const char *fileName, int n, int numOfCurrTimeIter, int numOfSolOrd, double numOfTimeSteps,
                        double *masOx, double *VectorOfData, int dimOfVect);

bool print_surface_as_v(const char *filename, int ox_len, int oy_len, double hx, double hy,
                        int t, double a, double b, double *data);

#endif //FEM_CIRCLE_TECPLOT_H
