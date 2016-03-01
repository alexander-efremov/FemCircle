#ifndef FEM_CIRCLE_UTILS_H
#define FEM_CIRCLE_UTILS_H

#include <float.h>
#include <stdio.h>
#include <math.h>

inline void print_matrix_to_file(int n, int m, double *data, const char *filename) {
    FILE *f = fopen(filename, "w");
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            int k = i * n + j;
            fprintf(f, "%20.14le ", data[k]);
        }
        fprintf(f, "\n ");
    }
    fclose(f);
}

inline void print_int(const char *str, int i) {
    printf("%s %d\n", str, i);
}

inline void print_int_double(const char *str, int i, double d) {
    printf("%s %d %f\n", str, i, d);
}

inline void print_double(const char *str, double d) {
    printf("%s %f\n", str, d);
}

inline void print_double(const char *str, double d1, double d2) {
    printf("%s %f %f\n", str, d1, d2);
}

inline void print_float(const char *str, float d) {
    printf("%s %f\n", str, d);
}

inline void print_double_exp(const char *str, double d) {
    printf("%s %e\n", str, d);
}

inline void print_float_exp(const char *str, float d) {
    printf("%s %e\n", str, d);
}

inline void print_matrix(double *a, int n, int m, int precision = 8) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            int k = i * n + j;
            switch (precision) {
                case 1:
                    printf("%.1le ", a[k]);
                    break;
                case 2:
                    printf("%.2le ", a[k]);
                    break;
                case 3:
                    printf("%.3le ", a[k]);
                    break;
                case 4:
                    printf("%.4le ", a[k]);
                    break;
                case 5:
                    printf("%.5le ", a[k]);
                    break;
                case 6:
                    printf("%.6le ", a[k]);
                    break;
                case 7:
                    printf("%.7le ", a[k]);
                    break;
                case 8:
                    printf("%.8le ", a[k]);
                    break;
            }
        }
        printf("\n");
    }
}

inline void print_matrix(double *a, int n, int m, const char *text, int precision = 8) {
    printf("\n%s\n", text);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            int k = i * n + j;
            switch (precision) {
                case 1:
                    printf("%.1le ", a[k]);
                    break;
                case 2:
                    printf("%.2le ", a[k]);
                    break;
                case 3:
                    printf("%.3le ", a[k]);
                    break;
                case 4:
                    printf("%.4le ", a[k]);
                    break;
                case 5:
                    printf("%.5le ", a[k]);
                    break;
                case 6:
                    printf("%.6le ", a[k]);
                    break;
                case 7:
                    printf("%.7le ", a[k]);
                    break;
                case 8:
                    printf("%.8le ", a[k]);
                    break;
            }
        }
        printf("\n");
    }
}

inline void print_vector(double *a, int n, int precision = 8) {
    for (int k = 0; k < n; ++k) {
        switch (precision) {
            case 1:
                printf("%.1f ", a[k]);
                break;
            case 2:
                printf("%.2f ", a[k]);
                break;
            case 3:
                printf("%.3f ", a[k]);
                break;
            case 4:
                printf("%.4f ", a[k]);
                break;
            case 5:
                printf("%.5f ", a[k]);
                break;
            case 6:
                printf("%.6f ", a[k]);
                break;
            case 7:
                printf("%.7f ", a[k]);
                break;
            case 8:
                printf("%.8f ", a[k]);
                break;
        }
    }
}

inline double get_l1_norm(double hx, double hy, int x_len, int y_len, double *data) {
    double r = 0.;
    for (int i = 0; i < x_len; ++i)
        for (int j = 0; j < y_len; ++j)
            r += fabs(data[x_len * i + j]);
    return r * hx * hy;
}

inline double get_l_inf_norm(int x_len, int y_len, double *data) {
    double max = FLT_MIN;
    for (int i = 0; i < x_len; ++i)
        for (int j = 0; j < y_len; ++j)
            if (fabs(data[x_len * i + j]) > max)
                max = fabs(data[x_len * i + j]);
    return max;
}

inline bool is_empty_file(FILE *f) {
    long savedOffset = ftell(f);
    fseek(f, 0, SEEK_END);
    if (ftell(f) == 0) {
        return true;
    }
    fseek(f, savedOffset, SEEK_SET);
    return false;
}

inline void append_statistics(int ox_len, int oy_len, double tau, int iterCount, double err_l1) {
    FILE *file;
    const char* filename = "/home/jane/ClionProjects/fem_circle/statistics.dat";
    file = fopen(filename, "a");
    if (file == NULL) {
        perror("Error opening file.");
        return;
    }
    if (is_empty_file(file)) {
        fprintf(file, "%s\t%s\t%s\t%s\t%s\n", "OX", "OY", "TAU", "ITERCOUNT", "L1ERR");
    }

    fprintf(file, "%d\t%d\t%le\t%d\t%le\n", ox_len, oy_len, tau, iterCount, err_l1);

    fclose(file);
}

#endif //FEM_CIRCLE_UTILS_H
