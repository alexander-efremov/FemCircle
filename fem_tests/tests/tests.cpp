#include <utils.h>
#include <cmath>
#include <common.h>
#include <dirent.h>
#include <fstream>
#include "gtest/gtest.h"
#include <algorithm>
#include <locale>

class FemFixture : public ::testing::Test {
protected:
    virtual void TearDown() {
        if (G1 != NULL)
            delete[] G1;
        if (G2 != NULL)
            delete[] G2;
        if (G3 != NULL)
            delete[] G3;
        if (G4 != NULL)
            delete[] G4;
    }

    virtual void SetUp() {
    }

public:
    FemFixture() : Test() {
    }
};

void print_params() {
    printf("\nNXxNY = %dx%d\n", NX, NY);
    printf("(U, V) = (%le, %le)\n", U, V);
    printf("(HX, HY) = (%le, %le)\n", HX, HY);
    printf("TAU = %le\n", TAU);
    printf("TIME_STEP_CNT = %d\n", TIME_STEP_CNT);
    printf("INTEGR_TYPE = %d\n", INTEGR_TYPE);
    printf("IDEAL_SQ_SIZE = %dx%d\n", IDEAL_SQ_SIZE_X, IDEAL_SQ_SIZE_Y);
    printf("CENTER_OFFSET = %le, %le\n", CENTER_OFFSET_X, CENTER_OFFSET_Y);
}

void init_boundary_arrays_and_cp() {
    G1 = new int[NX_1];
    G2 = new int[NY_1];
    G3 = new int[NX_1];
    G4 = new int[NY_1];
    for (int i = 0; i < NX_1; ++i) {
        G1[i] = 0;
        G3[i] = 0;
    }
    for (int j = 0; j < NY_1; ++j) {
        G2[j] = 0;
        G4[j] = 0;
    }
    CP00 = 0;
    CP10 = 0;
    CP01 = 0;
    CP11 = 0;
}

/*function... might want it in some class?*/
int getdir(std::string dir, std::vector<std::string> &files) {
    DIR *dp;
    struct dirent *dirp;
    if ((dp = opendir(dir.c_str())) == NULL) {
        std::cout << "Error opening " << dir << std::endl;
        return -1;
    }

    while ((dirp = readdir(dp)) != NULL) {
        files.push_back(std::string(dirp->d_name));
    }
    closedir(dp);
    return 0;
}

bool starts_with(const std::string &s1, const std::string &s2) {
    return s2.size() <= s1.size() && s1.compare(0, s2.size(), s2) == 0;
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

// trim from start
static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}

double *getArr(const int skipLines, std::string file, int size) {
    int idx = 0;
    double *arr = new double[size];
    std::ifstream input(file.c_str());
    std::string line;

    for (int j = 0; j < skipLines; ++j) {
        std::getline(input, line);
    }

    while (std::getline(input, line)) {
        std::vector<std::string> splt = split(line, ' ');
        splt[2] = trim(splt[2]);
        double d = atof(splt[2].c_str());
        arr[idx] = d;
    }
    return arr;
}


int get_tl(std::string &basic_string) {
    std::vector<std::string> splt = split(basic_string, '_');
    std::vector<std::string> t = split(splt[5], '=');
    return atoi(t[1].c_str());
}

TEST_F(FemFixture, tecplot2dTo1d) {
    const int skipLines = 6;
    A = 0.;
    B = 1.;
    C = 0.;
    D = 1.;
    NX = 400;
    NY = 400;
    NX_1 = NX + 1;
    NY_1 = NY + 1;
    XY = NX_1 * NY_1;
    R_SQ = 0.099 * 0.099;
    INN_DENSITY = 1.;
    OUT_DENSITY = 0.;

    HX = (B - A) / NX;
    HY = (D - C) / NY;
    IDEAL_SQ_SIZE_X = 128 * (3 + 1);
    IDEAL_SQ_SIZE_Y = 128 * (3 + 1);

    CENTER_OFFSET_X = 0.5;
    CENTER_OFFSET_Y = 0.5;

    INTEGR_TYPE = 1;

    U = 1.;
    V = 1.;
    OMEGA = 1.;
    TAU = 2.5e-3;

    TIME_STEP_CNT = 2700;

    double fixed_y = 0;
    double fixed_x = 0;

    std::string dir = std::string("/home/jane/ClionProjects/fem_circle/fem_tests/test3_1/");
    std::vector<std::string> files = std::vector<std::string>();

    getdir(dir, files);

    std::vector<std::string> errFiles = std::vector<std::string>();
    std::vector<std::string> rhoFiles = std::vector<std::string>();
    for (unsigned int i = 0; i < files.size(); i++) {
        if (starts_with(files[i], "err")) {
            errFiles.push_back(files[i]);
        }
    }

    for (unsigned int i = 0; i < files.size(); i++) {
        if (starts_with(files[i], "rho")) {
            rhoFiles.push_back(files[i]);
        }
    }
    std::cout << "Starting rho files convertsion..." << std::endl;
    for (unsigned int i = 0; i < rhoFiles.size(); i++) {
        std::cout << rhoFiles[i] << std::endl;
        int tl = get_tl(rhoFiles[i]);
        double *arr = getArr(skipLines, rhoFiles[i], XY);
        print_line_along_x("rho_x_1d", NX, NY, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                           U, V, arr, fixed_y);
        print_line_along_y("rho_y_1d", NX, NY, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                           U, V, arr, fixed_x);
        delete[] arr;
    }
    std::cout << "Starting err files convertsion..." << std::endl;
    for (unsigned int i = 0; i < errFiles.size(); i++) {
        std::cout << errFiles[i] << std::endl;
        int tl = get_tl(errFiles[i]);
        double *arr = getArr(skipLines, errFiles[i], XY);
        print_line_along_x("err_x_1d", NX, NY, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                           U, V, arr, fixed_y);
        print_line_along_y("err_x_1d", NX, NY, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                           U, V, arr, fixed_x);
        delete[] arr;
    }
}

TEST_F(FemFixture, test1) {
    double tme = 0.;

    double d = 0;

    for (int i = 0; i < 2; ++i) {
        switch (i) {
            case 0:
                d = 50.;
                break;
            case 1:
                d = 100.;
                break;
            case 2:
                d = 200.;
                break;
            case 3:
                d = 400.;
                break;
            case 4:
                d = 800.;
                break;
            case 5:
                d = 1600.;
                break;
            default:
                break;
        }

        A = 0.;
        B = 1.;
        C = 0.;
        D = 1.;
        R_SQ = 0.25 * 0.25;
        INN_DENSITY = 1.;
        OUT_DENSITY = 0.;
        JAK_ITER_CNT = 3000;

        TAU = 1. / d;
        NX = (int) d;
        NY = (int) d;
        NX_1 = NX + 1;
        NY_1 = NY + 1;
        HX = (B - A) / NX;
        HY = (D - C) / NY;
        TIME_STEP_CNT = (int) d;
        XY = NX_1 * NY_1;
        U = 0.;
        V = 0.;
        CENTER_OFFSET_X = NX * HX / 2.;
        CENTER_OFFSET_Y = NY * HY / 2.;

        print_params();

        double *density = solve_1(tme);
        double *err = calc_error_1(HX, HY, density);
        double x0 = get_center_x();
        double y0 = get_center_y();
        print_surface("test1_rho", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                      V, density);
        print_surface("test1_err", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                      V, err);
        double l1 = get_l1_norm(HX, HY, NX_1, NY_1, err);
        double l_inf = get_l_inf_norm(NX_1, NY_1, err);
        printf("l1 %le \n", l1);
        printf("l_inf %le\n", l_inf);
        delete[] density;
        delete[] err;
    }
}

TEST_F(FemFixture, test2) {
    double tme = 0.;

    double d = 0;

    for (int i = 2; i < 4; ++i) {
        switch (i) {
            case 0:
                d = 50.;
                break;
            case 1:
                d = 100.;
                break;
            case 2:
                d = 200.;
                break;
            case 3:
                d = 400.;
                break;
            case 4:
                d = 800.;
                break;
            case 5:
                d = 1600.;
                break;
            default:
                break;
        }

        A = 0.;
        B = 1.;
        C = 0.;
        D = 1.;
        R_SQ = 0.1 * 0.1;
        INN_DENSITY = 1.;
        OUT_DENSITY = 0.;

        NX = (int) d;
        NY = (int) d;
        NX_1 = NX + 1;
        NY_1 = NY + 1;
        HX = (B - A) / NX;
        HY = (D - C) / NY;

        // (u, v) = (a, a)
        double a = 1;
        U = a;
        V = a;
        TAU = std::min(HX / (3. * a), HY / (3. * a));

        print_params();

        TIME_STEP_CNT = (int) 1;
        XY = NX_1 * NY_1;

        printf("NX = %d NY = %d\n", NX, NY);
        double *density = solve_2(tme);
        double *err = calc_error_2(HX, HY, TAU * TIME_STEP_CNT, density);
        double x0 = get_center_x();
        double y0 = get_center_y();
        print_surface("test2_rho", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                      V, density);
        print_surface("test2_err", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                      V, err);
        double l1 = get_l1_norm(HX, HY, NX_1, NY_1, err);
        double l_inf = get_l_inf_norm(NX_1, NY_1, err);
        printf("l1 %le \n", l1);
        printf("l_inf %le\n", l_inf);
        delete[] density;
        delete[] err;
    }
}

TEST_F(FemFixture, test2_1) {
    double tme = 0.;
    A = 0.;
    B = 1.;
    C = 0.;
    D = 1.;
    R_SQ = 0.099 * 0.099;
    INN_DENSITY = 1.;
    OUT_DENSITY = 0.;

    NX = (int) 400;
    NY = (int) 400;
    NX_1 = NX + 1;
    NY_1 = NY + 1;
    HX = (B - A) / NX;
    HY = (D - C) / NY;
    IDEAL_SQ_SIZE_X = 64;
    IDEAL_SQ_SIZE_Y = 64;
    CENTER_OFFSET_X = 0.3;
    CENTER_OFFSET_Y = 0.3;

    INTEGR_TYPE = 1;

    U = 1.;
    V = 1.;
    TAU = 1.e-3;
    //TIME_STEP_CNT = (int) ((1 - get_center_x_2() - get_center_y_2()) / TAU);
    TIME_STEP_CNT = 2;
    XY = NX_1 * NY_1;

    print_params();

    double *density = solve_2(tme);
    double *err = calc_error_2(HX, HY, TAU * TIME_STEP_CNT, density);
    double *exact0 = get_exact_solution_2(HX, HY, 0);
    double *exactT = get_exact_solution_2(HX, HY, TAU * TIME_STEP_CNT);

    double y0 = get_center_y();
    double x0 = get_center_x();
    print_surface("rho", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                  V, density);
    print_surface("err", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                  V, err);
    print_surface("exact", NX, NY, HX, HY, 0, A, C, x0, y0, TAU, U,
                  V, exact0);
    print_surface("exact", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                  V, exactT);
    double l1 = get_l1_norm(HX, HY, NX_1, NY_1, err);
    double l_inf = get_l_inf_norm(NX_1, NY_1, err);
    printf("l1 %le\n", l1);
    printf("l_inf %le\n", l_inf);
    delete[] density;
    delete[] exact0;
    delete[] exactT;
    delete[] err;
}

TEST_F(FemFixture, test2_2) {
    double tme = 0.;
    for (int iter = 0; iter < 1; ++iter) {

        double d = 0;
        for (int i = 3; i < 4; ++i) {
            switch (i) {
                case 0:
                    d = 50.;
                    break;
                case 1:
                    d = 100.;
                    break;
                case 2:
                    d = 200.;
                    break;
                case 3:
                    d = 400.;
                    break;
                case 4:
                    d = 800.;
                    break;
                case 5:
                    d = 1600.;
                    break;
                default:
                    return;
            }

            A = 0.;
            B = 1.;
            C = 0.;
            D = 1.;
            R_SQ = 0.099 * 0.099;
            INN_DENSITY = 1.;
            OUT_DENSITY = 0.;

            NX = (int) d;
            NY = (int) d;
            NX_1 = NX + 1;
            NY_1 = NY + 1;
            HX = (B - A) / NX;
            HY = (D - C) / NY;
            IDEAL_SQ_SIZE_X = 4;
            IDEAL_SQ_SIZE_Y = 4;
            CENTER_OFFSET_X = 0.3;
            CENTER_OFFSET_Y = 0.3;

            INTEGR_TYPE = 1;

            U = 1.;
            V = 1.;
            //TAU = 16. / pow(2., (i + 1));
            TAU = 1.2e-3;
            TIME_STEP_CNT = (int) pow(2., i);
            TIME_STEP_CNT = 20;
            XY = NX_1 * NY_1;

            print_params();

            double *density = solve_2(tme);
            double *err = calc_error_2(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_2(HX, HY, 0);
            double *exactT = get_exact_solution_2(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, density);
            print_surface("err", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, err);
            print_surface("exact", NX, NY, HX, HY, 0, A, C, x0, y0, TAU, U,
                          V, exact0);
            print_surface("exact", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, exactT);

            double l1 = get_l1_norm(HX, HY, NX_1, NY_1, err);
            double l_inf = get_l_inf_norm(NX_1, NY_1, err);
            printf("l1 %le \n", l1);
            printf("l_inf %le\n", l_inf);
            delete[] density;
            delete[] exact0;
            delete[] exactT;
            delete[] err;
        }
    }
}

// тестируем вылет пятна за границу
TEST_F(FemFixture, test2_3) {
    double tme = 0.;
    for (int iter = 0; iter < 1; ++iter) {

        double d = 0;
        for (int i = 3; i < 4; ++i) {
            switch (i) {
                case 0:
                    d = 50.;
                    break;
                case 1:
                    d = 100.;
                    break;
                case 2:
                    d = 200.;
                    break;
                case 3:
                    d = 400.;
                    break;
                case 4:
                    d = 800.;
                    break;
                case 5:
                    d = 1600.;
                    break;
                default:
                    return;
            }

            A = 0.;
            B = 1.;
            C = 0.;
            D = 1.;
            R_SQ = 0.099 * 0.099;
            INN_DENSITY = 1.;
            OUT_DENSITY = 0.;

            NX = (int) d;
            NY = (int) d;
            NX_1 = NX + 1;
            NY_1 = NY + 1;
            HX = (B - A) / NX;
            HY = (D - C) / NY;
            IDEAL_SQ_SIZE_X = 128 * (iter + 1);
            IDEAL_SQ_SIZE_Y = 128 * (iter + 1);

            CENTER_OFFSET_X = 0.85;
            CENTER_OFFSET_Y = 0.5;

            INTEGR_TYPE = 1;

            U = 1.;
            V = 1.;
            TAU = 16. / pow(2., (i + 1));
            TAU *= 1.e-3;
            //TIME_STEP_CNT = (int) pow(2., i);
            TIME_STEP_CNT = 1;
            XY = NX_1 * NY_1;

            print_params();

            double *density = solve_2(tme);
            double *err = calc_error_2(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_2(HX, HY, 0);
            double *exactT = get_exact_solution_2(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, density);
            print_surface("err", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, err);
            print_surface("exact", NX, NY, HX, HY, 0, A, C, x0, y0, TAU, U,
                          V, exact0);
            print_surface("exact", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, exactT);

            double l1 = get_l1_norm(HX, HY, NX_1, NY_1, err);
            double l_inf = get_l_inf_norm(NX_1, NY_1, err);
            printf("l1 %le \n", l1);
            printf("l_inf %le\n", l_inf);
            delete[] density;
            delete[] exact0;
            delete[] exactT;
            delete[] err;
        }
    }
}

// тестируем проход по диагонали с переходом по ячейкам
TEST_F(FemFixture, test2_4) {
    double tme = 0.;
    IDEAL_SQ_SIZE_X = 128;
    IDEAL_SQ_SIZE_Y = 128;
    for (int iter = 0; iter < 1; ++iter) {

        double d = 0;
        for (int i = 3; i < 4; ++i) {
            switch (i) {
                case 0:
                    d = 50.;
                    break;
                case 1:
                    d = 100.;
                    break;
                case 2:
                    d = 200.;
                    break;
                case 3:
                    d = 400.;
                    break;
                case 4:
                    d = 800.;
                    break;
                case 5:
                    d = 1600.;
                    break;
                default:
                    return;
            }

            A = 0.;
            B = 1.;
            C = 0.;
            D = 1.;
//            R_SQ = 0.099 * 0.099;
            R_SQ = 0.0999 * 0.0999;
            INN_DENSITY = 1.;
            OUT_DENSITY = 0.;

            NX = (int) d;
            NY = (int) d;
            NX_1 = NX + 1;
            NY_1 = NY + 1;
            HX = (B - A) / NX;
            HY = (D - C) / NY;
            //IDEAL_SQ_SIZE_X = 128 * (iter + 1);
            //IDEAL_SQ_SIZE_Y = 128 * (iter + 1);


            CENTER_OFFSET_X = 0.3;
            CENTER_OFFSET_Y = 0.3;

            INTEGR_TYPE = 4;

            U = 1.;
            V = 1.;
//            TAU = HX;
            TAU = 0.001;

            //TIME_STEP_CNT = (int) pow(2., i);
            //TIME_STEP_CNT = 80;
            TIME_STEP_CNT = 1;
            //TIME_STEP_CNT = 1 + (0.7-CENTER_OFFSET_X)/TAU; // 4*2.5e-3 = 0.01 0.7-0.3=0.4/2.5e-3=160

            XY = NX_1 * NY_1;

            print_params();

            double *density = solve_2(tme);
            double *err = calc_error_2(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_2(HX, HY, 0);
            double *exactT = get_exact_solution_2(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, density);
            print_surface("err", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, err);
            print_surface("exact", NX, NY, HX, HY, 0, A, C, x0, y0, TAU, U,
                          V, exact0);
            print_surface("exact", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, exactT);

            double l1 = get_l1_norm_vec(NX_1, NY_1, err);
            printf("l1_err_vec %le \n", l1);
            l1 = get_l1_norm_int_trapezoidal(HX, HY, NX, NY, err); // note! a loop boundary
            double l_inf = get_l_inf_norm(NX_1, NY_1, err);
            printf("l1_err_tr %le \n", l1);
            printf("l_inf %le\n", l_inf);
            delete[] density;
            delete[] exact0;
            delete[] exactT;
            delete[] err;
            IDEAL_SQ_SIZE_X *= 2;
            IDEAL_SQ_SIZE_Y *= 2;
        }
    }
}

// тестируем проход по диагонали с переходом по ячейкам
// шаг по времени = hx/2
// интегрирование приближ
TEST_F(FemFixture, test2_5) {
    double tme = 0.;
    for (int iter = 0; iter < 1; ++iter) {

        double d = 0;
        for (int i = 3; i < 4; ++i) {
            switch (i) {
                case 0:
                    d = 50.;
                    break;
                case 1:
                    d = 100.;
                    break;
                case 2:
                    d = 200.;
                    break;
                case 3:
                    d = 400.;
                    break;
                case 4:
                    d = 800.;
                    break;
                case 5:
                    d = 1600.;
                    break;
                default:
                    return;
            }

            A = 0.;
            B = 1.;
            C = 0.;
            D = 1.;
            R_SQ = 0.099 * 0.099;
            INN_DENSITY = 1.;
            OUT_DENSITY = 0.;

            NX = (int) d;
            NY = (int) d;
            NX_1 = NX + 1;
            NY_1 = NY + 1;
            HX = (B - A) / NX;
            HY = (D - C) / NY;
            IDEAL_SQ_SIZE_X = 128 * (iter + 1);
            IDEAL_SQ_SIZE_Y = 128 * (iter + 1);

            CENTER_OFFSET_X = 0.3;
            CENTER_OFFSET_Y = 0.3;

            INTEGR_TYPE = 1;

            U = 1.;
            V = 1.;
            TAU = HX / 2.;
            //TIME_STEP_CNT = (int) pow(2., i);
            //TIME_STEP_CNT = 1;
            TIME_STEP_CNT = (0.7 - CENTER_OFFSET_X) / TAU; // 4*2.5e-3 = 0.01 0.7-0.3=0.4/2.5e-3=160

            XY = NX_1 * NY_1;

            print_params();

            double *density = solve_2(tme);
            double *err = calc_error_2(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_2(HX, HY, 0);
            double *exactT = get_exact_solution_2(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, density);
            print_surface("err", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, err);
            print_surface("exact", NX, NY, HX, HY, 0, A, C, x0, y0, TAU, U,
                          V, exact0);
            print_surface("exact", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, exactT);

            double l1 = get_l1_norm(HX, HY, NX_1, NY_1, err);
            double l_inf = get_l_inf_norm(NX_1, NY_1, err);
            printf("l1 %le \n", l1);
            printf("l_inf %le\n", l_inf);
            delete[] density;
            delete[] exact0;
            delete[] exactT;
            delete[] err;
        }
    }
}

// тестируем проход по диагонали с переходом по ячейкам
// шаг по времени = hx/2
// интегрирование точное
TEST_F(FemFixture, test2_6) {
    double tme = 0.;
    for (int iter = 0; iter < 1; ++iter) {

        double d = 0;
        for (int i = 3; i < 4; ++i) {
            switch (i) {
                case 0:
                    d = 50.;
                    break;
                case 1:
                    d = 100.;
                    break;
                case 2:
                    d = 200.;
                    break;
                case 3:
                    d = 400.;
                    break;
                case 4:
                    d = 800.;
                    break;
                case 5:
                    d = 1600.;
                    break;
                default:
                    return;
            }

            A = 0.;
            B = 1.;
            C = 0.;
            D = 1.;
            R_SQ = 0.099 * 0.099;
            INN_DENSITY = 1.;
            OUT_DENSITY = 0.;

            NX = (int) d;
            NY = (int) d;
            NX_1 = NX + 1;
            NY_1 = NY + 1;
            HX = (B - A) / NX;
            HY = (D - C) / NY;
            IDEAL_SQ_SIZE_X = 128 * (iter + 1);
            IDEAL_SQ_SIZE_Y = 128 * (iter + 1);

            CENTER_OFFSET_X = 0.3;
            CENTER_OFFSET_Y = 0.3;

            INTEGR_TYPE = 3;

            U = 1.;
            V = 1.;
            TAU = HX / 2.;
            //TIME_STEP_CNT = (int) pow(2., i);
            //TIME_STEP_CNT = 1;
            TIME_STEP_CNT = (0.7 - CENTER_OFFSET_X) / TAU; // 4*2.5e-3 = 0.01 0.7-0.3=0.4/2.5e-3=160

            XY = NX_1 * NY_1;

            print_params();

            double *density = solve_2(tme);
            double *err = calc_error_2(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_2(HX, HY, 0);
            double *exactT = get_exact_solution_2(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, density);
            print_surface("err", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, err);
            print_surface("exact", NX, NY, HX, HY, 0, A, C, x0, y0, TAU, U,
                          V, exact0);
            print_surface("exact", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, exactT);

            double l1 = get_l1_norm(HX, HY, NX_1, NY_1, err);
            double l_inf = get_l_inf_norm(NX_1, NY_1, err);
            printf("l1 %le \n", l1);
            printf("l_inf %le\n", l_inf);
            delete[] density;
            delete[] exact0;
            delete[] exactT;
            delete[] err;
        }
    }
}

// тестируем третий случай - движение по кругу
TEST_F(FemFixture, test3_1) {
    double tme = 0.;
    for (int iter = 0; iter < 1; ++iter) {

        double d = 0;
        for (int i = 1; i < 2; ++i) {
            switch (i) {
                case 0:
                    d = 50.;
                    break;
                case 1:
                    d = 100.;
                    break;
                case 2:
                    d = 200.;
                    break;
                case 3:
                    d = 400.;
                    break;
                case 4:
                    d = 800.;
                    break;
                case 5:
                    d = 1600.;
                    break;
                default:
                    return;
            }

            A = 0.;
            B = 1.;
            C = 0.;
            D = 1.;
            R_SQ = 0.099 * 0.099;
            INN_DENSITY = 1.;
            OUT_DENSITY = 0.;

            NX = (int) d;
            NY = (int) d;
            NX_1 = NX + 1;
            NY_1 = NY + 1;
            HX = (B - A) / NX;
            HY = (D - C) / NY;
            IDEAL_SQ_SIZE_X = 128 * (iter + 1);
            IDEAL_SQ_SIZE_Y = 128 * (iter + 1);

            CENTER_OFFSET_X = 0.5;
            CENTER_OFFSET_Y = 0.5;

            INTEGR_TYPE = 1;

            U = 1.;
            V = 1.;
            OMEGA = 1.;
            TAU = 2.5e-3;

            TIME_STEP_CNT = 1;
            XY = NX_1 * NY_1;

            init_boundary_arrays_and_cp();

            int midIndexX = NX_1 / 2;
            int midIndexY = NY_1 / 2;
            for (int i = 0; i < NX_1; ++i) {
                if (i * HX + A < .5 && i > 0)
                    G1[i] = 1;
            }
            for (int j = 0; j < NY_1; ++j) {
                if (j * HY + C < .5 && j > 0)
                    G2[j] = 1;
            }
            for (int i = 0; i < NX_1; ++i) {
                if (i * HX + A > .5 && i < NX_1 - 1)
                    G3[i] = 1;
            }
            for (int j = 0; j < NY_1; ++j) {
                if (j * HY + C > .5 && j > 0 && j < NX_1 - 1) {
                    G4[j] = 1;
                }
            }
/*
            CP00 = 1;
            CP10 = 1;
            CP01 = 1;
            CP11 = 1;

            */

            CP00 = 0;
            CP10 = 0;
            CP01 = 0;
            CP11 = 0;

            print_params();
            printf("rel = %le\n", HX / (-HY + 1.));
            printf("midIndexX = %d\n", midIndexX);
            printf("midIndexY = %d\n", midIndexY);
            printf("G1\n");
            print_vector(G1, NX_1);
            printf("G2\n");
            print_vector(G2, NY_1);
            printf("G3\n");
            print_vector(G3, NX_1);
            printf("G4\n");
            print_vector(G4, NY_1);

            double *density = solve_3(tme);
            double *err = calc_error_3(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_3(HX, HY, 0);
            double *exactT = get_exact_solution_3(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, density);
            print_surface("err", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, err);
            print_surface("exact", NX, NY, HX, HY, 0, A, C, x0, y0, TAU, U,
                          V, exact0);
            print_surface("exact", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, exactT);

            double l1 = get_l1_norm_vec(NX_1, NY_1, err);
            //           double l1 = get_l1_norm_int_middle(HX, HY, NX_1, NY_1, err);
            double l_inf = get_l_inf_norm(NX_1, NY_1, err);
            printf("l1 %le \n", l1);
            printf("l_inf %le\n", l_inf);
            delete[] density;
            delete[] exact0;
            delete[] exactT;
            delete[] err;
        }
    }
}

// тестируем третий случай - движение по кругу
TEST_F(FemFixture, test3_2) {
    double tme = 0.;
    for (int iter = 0; iter < 1; ++iter) {

        double d = 0;
        for (int i = 2; i < 3; ++i) {
            switch (i) {
                case 0:
                    d = 50.;
                    break;
                case 1:
                    d = 100.;
                    break;
                case 2:
                    d = 200.;
                    break;
                case 3:
                    d = 400.;
                    break;
                case 4:
                    d = 800.;
                    break;
                case 5:
                    d = 1600.;
                    break;
                default:
                    return;
            }

            A = 0.;
            B = 1.;
            C = 0.;
            D = 1.;
            R_SQ = 0.099 * 0.099;
            INN_DENSITY = 1.;
            OUT_DENSITY = 0.;

            NX = (int) d;
            NY = (int) d;
            NX_1 = NX + 1;
            NY_1 = NY + 1;
            HX = (B - A) / NX;
            HY = (D - C) / NY;
            IDEAL_SQ_SIZE_X = 128 * (iter + 1);
            IDEAL_SQ_SIZE_Y = 128 * (iter + 1);

            CENTER_OFFSET_X = 0.5;
            CENTER_OFFSET_Y = 0.5;

            INTEGR_TYPE = 1;

            U = 1.;
            V = 1.;
            OMEGA = 1.;
            TAU = 2.5e-3;

            TIME_STEP_CNT = 50;
            XY = NX_1 * NY_1;

            double x0 = get_center_x();
            double y0 = get_center_y();

            double *exact0 = get_exact_solution_3(HX, HY, 0.);
            double *exactT = get_exact_solution_3(HX, HY, TAU * TIME_STEP_CNT);
            for (int j = 0; j < TIME_STEP_CNT; ++j) {
                double *exact0 = get_exact_solution_3(HX, HY, j * TAU);
                print_surface("exact", NX, NY, HX, HY, j, A, C, x0, y0, TAU, U,
                              V, exact0);
            }

            print_surface("exact", NX, NY, HX, HY, 0, A, C, x0, y0, TAU, U,
                          V, exact0);
            print_surface("exact", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, exactT);


            delete[] exact0;
            delete[] exactT;
        }
    }
}

// тестируем третий случай - движение по кругу
// сходимость
TEST_F(FemFixture, test3_3) {
    double tme = 0.;
    for (int iter = 0; iter < 1; ++iter) {

        double d = 0;
        for (int i = 0; i < 5; ++i) {
            switch (i) {
                case 0:
                    d = 50.;
                    break;
                case 1:
                    d = 100.;
                    break;
                case 2:
                    d = 200.;
                    break;
                case 3:
                    d = 400.;
                    break;
                case 4:
                    d = 800.;
                    break;
                case 5:
                    d = 1600.;
                    break;
                default:
                    return;
            }

            A = 0.;
            B = 1.;
            C = 0.;
            D = 1.;
            R_SQ = 0.099 * 0.099;
            INN_DENSITY = 1.;
            OUT_DENSITY = 0.;

            NX = (int) d;
            NY = (int) d;
            NX_1 = NX + 1;
            NY_1 = NY + 1;
            HX = (B - A) / NX;
            HY = (D - C) / NY;
            IDEAL_SQ_SIZE_X = 128 * (iter + 1);
            IDEAL_SQ_SIZE_Y = 128 * (iter + 1);

            CENTER_OFFSET_X = 0.5;
            CENTER_OFFSET_Y = 0.7;

            INTEGR_TYPE = 1;

            U = 1.;
            V = 1.;
            OMEGA = 1.;

            TAU = 16. / pow(2., (i + 1));
            TAU *= 1.e-3;

            TIME_STEP_CNT = (int) pow(2., i);
            XY = NX_1 * NY_1;

            init_boundary_arrays_and_cp();

            int midIndexX = NX_1 / 2;
            int midIndexY = NY_1 / 2;
            for (int i = 0; i < NX_1; ++i) {
                if (i * HX + A < .5 && i > 0)
                    G1[i] = 1;
            }
            for (int j = 0; j < NY_1; ++j) {
                if (j * HY + C < .5 && j > 0)
                    G2[j] = 1;
            }
            for (int i = 0; i < NX_1; ++i) {
                if (i * HX + A > .5 && i < NX_1 - 1)
                    G3[i] = 1;
            }
            for (int j = 0; j < NY_1; ++j) {
                if (j * HY + C > .5 && j > 0 && j < NX_1 - 1) {
                    G4[j] = 1;
                }
            }

            CP00 = 1;
            CP10 = 1;
            CP01 = 1;
            CP11 = 1;

            /*

             CP00 = 0;
             CP10 = 0;
             CP01 = 0;
             CP11 = 0;

             */

            print_params();
            printf("rel = %le\n", HX / (-HY + 1.));
            printf("midIndexX = %d\n", midIndexX);
            printf("midIndexY = %d\n", midIndexY);
            printf("G1\n");
            print_vector(G1, NX_1);
            printf("G2\n");
            print_vector(G2, NY_1);
            printf("G3\n");
            print_vector(G3, NX_1);
            printf("G4\n");
            print_vector(G4, NY_1);

            double *density = solve_3(tme);
            double *err = calc_error_3(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_3(HX, HY, 0);
            double *exactT = get_exact_solution_3(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, density);
            print_surface("err", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, err);
            print_surface("exact", NX, NY, HX, HY, 0, A, C, x0, y0, TAU, U,
                          V, exact0);
            print_surface("exact", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, exactT);

            double l1 = get_l1_norm(HX, HY, NX_1, NY_1, err);
            double l_inf = get_l_inf_norm(NX_1, NY_1, err);
            printf("l1 %le \n", l1);
            printf("l_inf %le\n", l_inf);
            delete[] density;
            delete[] exact0;
            delete[] exactT;
            delete[] err;
        }
    }
}

TEST_F(FemFixture, test4_1) {
    double tme = 0.;
    for (int iter = 0; iter < 1; ++iter) {

        double d = 0;
        for (int i = 3; i < 4; ++i) {
            switch (i) {
                case 0:
                    d = 50.;
                    break;
                case 1:
                    d = 100.;
                    break;
                case 2:
                    d = 200.;
                    break;
                case 3:
                    d = 400.;
                    break;
                case 4:
                    d = 800.;
                    break;
                case 5:
                    d = 1600.;
                    break;
                default:
                    return;
            }

            A = 0.;
            B = 1.;
            C = 0.;
            D = 1.;
            R_SQ = 0.099 * 0.099;
            INN_DENSITY = 1.;
            OUT_DENSITY = 0.;

            NX = (int) d;
            NY = (int) d;
            NX_1 = NX + 1;
            NY_1 = NY + 1;
            HX = (B - A) / NX;
            HY = (D - C) / NY;
            IDEAL_SQ_SIZE_X = 32 * (iter + 1);
            IDEAL_SQ_SIZE_Y = 32 * (iter + 1);

            CENTER_OFFSET_X = 0.3;
            CENTER_OFFSET_Y = 0.3;

            INTEGR_TYPE = 1;

            U = 1.;
            V = 1.;

            //TAU = 16. / pow(2., (i + 1));
            TAU = 1.e-3;

            //TIME_STEP_CNT = (int) pow(2., i);
            TIME_STEP_CNT = 2;
            XY = NX_1 * NY_1;

            init_boundary_arrays_and_cp();

            int midIndexX = NX_1 / 2;
            int midIndexY = NY_1 / 2;
            for (int i = 1; i < NX; ++i) {
                G1[i] = 0;
                G3[i] = 1;
            }
            for (int j = 1; j < NY; ++j) {
                G2[j] = 1;
                G4[j] = 0;
            }

            CP00 = 0;
            CP10 = 0;
            CP01 = 0;
            CP11 = 1;

            print_params();
            printf("rel = %le\n", HX / (-HY + 1.));
            printf("midIndexX = %d\n", midIndexX);
            printf("midIndexY = %d\n", midIndexY);
            printf("G1\n");
            print_vector(G1, NX_1);
            printf("G2\n");
            print_vector(G2, NY_1);
            printf("G3\n");
            print_vector(G3, NX_1);
            printf("G4\n");
            print_vector(G4, NY_1);

            double *density = solve_4(tme);
            double *err = calc_error_4(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_4(HX, HY, 0);
            double *exactT = get_exact_solution_4(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, density);
            print_surface("err", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, err);
            print_surface("exact", NX, NY, HX, HY, 0, A, C, x0, y0, TAU, U,
                          V, exact0);
            print_surface("exact", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, exactT);

            double l1 = get_l1_norm(HX, HY, NX_1, NY_1, err);
            double l_inf = get_l_inf_norm(NX_1, NY_1, err);
            printf("l1 %le \n", l1);
            printf("l_inf %le\n", l_inf);
            delete[] density;
            delete[] exact0;
            delete[] exactT;
            delete[] err;
        }
    }
}

// тестируем солвер 5 движение из угла в угол
// подход с накоплением коэффициентов
// оказался медленным и размазывает решение
TEST_F(FemFixture, test5_1) {
    double tme = 0.;
    for (int iter = 0; iter < 1; ++iter) {

        double d = 0;
        for (int i = 3; i < 4; ++i) {
            switch (i) {
                case 0:
                    d = 50.;
                    break;
                case 1:
                    d = 100.;
                    break;
                case 2:
                    d = 200.;
                    break;
                case 3:
                    d = 400.;
                    break;
                case 4:
                    d = 800.;
                    break;
                case 5:
                    d = 1600.;
                    break;
                default:
                    return;
            }

            A = 0.;
            B = 1.;
            C = 0.;
            D = 1.;
            R_SQ = 0.099 * 0.099;
            INN_DENSITY = 1.;
            OUT_DENSITY = 0.;

            NX = (int) d;
            NY = (int) d;
            NX_1 = NX + 1;
            NY_1 = NY + 1;
            HX = (B - A) / NX;
            HY = (D - C) / NY;
            IDEAL_SQ_SIZE_X = 64;
            IDEAL_SQ_SIZE_Y = 64;

            CENTER_OFFSET_X = 0.3;
            CENTER_OFFSET_Y = 0.3;

            INTEGR_TYPE = 1;

            U = 1.;
            V = 1.;
            OMEGA = 1.;
            TAU = 1.2475e-3;

            TIME_STEP_CNT = 322;
            XY = NX_1 * NY_1;

            init_boundary_arrays_and_cp();

            int midIndexX = NX_1 / 2;
            int midIndexY = NY_1 / 2;

            for (int i = 0; i < NX_1; ++i) {
                if (i * HX + A < .5 && i > 0)
                    G1[i] = 0;
            }
            for (int j = 0; j < NY_1; ++j) {
                if (j * HY + C < .5 && j > 0)
                    G2[j] = 0;
            }
            for (int i = 0; i < NX_1; ++i) {
                if (i * HX + A > .5 && i < NX_1 - 1)
                    G3[i] = 0;
            }
            for (int j = 0; j < NY_1; ++j) {
                if (j * HY + C > .5 && j > 0 && j < NX_1 - 1) {
                    G4[j] = 0;
                }
            }

            CP00 = 0;
            CP10 = 0;
            CP01 = 0;
            CP11 = 0;

            print_params();
            printf("rel = %le\n", HX / (-HY + 1.));
            printf("midIndexX = %d\n", midIndexX);
            printf("midIndexY = %d\n", midIndexY);

            double *density = solve_5(tme);
            double *err = calc_error_5(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_5(HX, HY, 0);
            double *exactT = get_exact_solution_5(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, density);
            print_surface("err", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, err);
            print_surface("exact", NX, NY, HX, HY, 0, A, C, x0, y0, TAU, U,
                          V, exact0);
            print_surface("exact", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, exactT);

            double l1 = get_l1_norm_vec(NX_1, NY_1, err);
            double l_inf = get_l_inf_norm(NX_1, NY_1, err);
            printf("l1 %le \n", l1);
            printf("l_inf %le\n", l_inf);
            delete[] density;
            delete[] exact0;
            delete[] exactT;
            delete[] err;
        }
    }
}

// тестируем солвер 6 движение из угла в угол
// подход с выбором четверти
//
TEST_F(FemFixture, test6_1) {
    double tme = 0.;
    for (int iter = 0; iter < 1; ++iter) {

        double d = 0;
        for (int i = 3; i < 4; ++i) {
            switch (i) {
                case 0:
                    d = 50.;
                    break;
                case 1:
                    d = 100.;
                    break;
                case 2:
                    d = 200.;
                    break;
                case 3:
                    d = 400.;
                    break;
                case 4:
                    d = 800.;
                    break;
                case 5:
                    d = 1600.;
                    break;
                default:
                    return;
            }

            A = 0.;
            B = 1.;
            C = 0.;
            D = 1.;
            R_SQ = 0.099 * 0.099;
            INN_DENSITY = 1.;
            OUT_DENSITY = 0.;

            NX = (int) d;
            NY = (int) d;
            NX_1 = NX + 1;
            NY_1 = NY + 1;
            HX = (B - A) / NX;
            HY = (D - C) / NY;
            IDEAL_SQ_SIZE_X = 64;
            IDEAL_SQ_SIZE_Y = 64;

            CENTER_OFFSET_X = 0.3;
            CENTER_OFFSET_Y = 0.3;

            INTEGR_TYPE = 1;

            U = 1.;
            V = 1.;
            OMEGA = 1.;
            TAU = 1.2475e-3;

            TIME_STEP_CNT = 322;
            XY = NX_1 * NY_1;

            init_boundary_arrays_and_cp();

            int midIndexX = NX_1 / 2;
            int midIndexY = NY_1 / 2;

            for (int i = 0; i < NX_1; ++i) {
                if (i * HX + A < .5 && i > 0)
                    G1[i] = 0;
            }
            for (int j = 0; j < NY_1; ++j) {
                if (j * HY + C < .5 && j > 0)
                    G2[j] = 1;
            }
            for (int i = 0; i < NX_1; ++i) {
                if (i * HX + A > .5 && i < NX_1 - 1)
                    G3[i] = 1;
            }
            for (int j = 0; j < NY_1; ++j) {
                if (j * HY + C > .5 && j > 0 && j < NX_1 - 1) {
                    G4[j] = 0;
                }
            }

            CP00 = 0;
            CP10 = 0;
            CP01 = 0;
            CP11 = 1;

            print_params();
            printf("rel = %le\n", HX / (-HY + 1.));
            printf("midIndexX = %d\n", midIndexX);
            printf("midIndexY = %d\n", midIndexY);

            double *density = solve_6(tme);
            double *err = calc_error_6(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_6(HX, HY, 0);
            double *exactT = get_exact_solution_6(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, density);
            print_surface("err", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, err);
            print_surface("exact", NX, NY, HX, HY, 0, A, C, x0, y0, TAU, U,
                          V, exact0);
            print_surface("exact", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, exactT);

            double l1 = get_l1_norm_vec(NX_1, NY_1, err);
            double l_inf = get_l_inf_norm(NX_1, NY_1, err);
            printf("l1 %le \n", l1);
            printf("l_inf %le\n", l_inf);
            delete[] density;
            delete[] exact0;
            delete[] exactT;
            delete[] err;
        }
    }
}

// тестируем солвер 7 движение по кругу
// подход с выбором четверти
//
TEST_F(FemFixture, test7_1) {
    double tme = 0.;
    for (int iter = 0; iter < 1; ++iter) {

        double d = 0;
        for (int i = 4; i < 5; ++i) {
            switch (i) {
                case 0:
                    d = 10.;
                    break;
                case 1:
                    d = 100.;
                    break;
                case 2:
                    d = 200.;
                    break;
                case 3:
                    d = 400.;
                    break;
                case 4:
                    d = 800.;
                    break;
                case 5:
                    d = 1600.;
                    break;
                default:
                    return;
            }

            A = 0.;
            B = 1.;
            C = 0.;
            D = 1.;
            R_SQ = 0.099 * 0.099;
            INN_DENSITY = 1.;
            OUT_DENSITY = 0.;

            NX = (int) d;
            NY = (int) d;
            NX_1 = NX + 1;
            NY_1 = NY + 1;
            HX = (B - A) / NX;
            HY = (D - C) / NY;
            IDEAL_SQ_SIZE_X = 64;
            IDEAL_SQ_SIZE_Y = 64;

            CENTER_OFFSET_X = 0.5;
            CENTER_OFFSET_Y = 0.5;

            INTEGR_TYPE = 1;

            U = 1.;
            V = 1.;
            OMEGA = 1.;
            TAU = 5. * HX;

            TIME_STEP_CNT = 270;
            XY = NX_1 * NY_1;

            init_boundary_arrays_and_cp();

            int midIndexX = NX_1 / 2;
            int midIndexY = NY_1 / 2;

            for (int i = 0; i < NX_1; ++i) {
                if (i * HX + A < .5 && i > 0)
                    G1[i] = 1;
            }
            for (int j = 0; j < NY_1; ++j) {
                if (j * HY + C < .5 && j > 0)
                    G2[j] = 1;
            }
            for (int i = 0; i < NX_1; ++i) {
                if (i * HX + A > .5 && i < NX_1 - 1)
                    G3[i] = 1;
            }
            for (int j = 0; j < NY_1; ++j) {
                if (j * HY + C > .5 && j > 0 && j < NX_1 - 1) {
                    G4[j] = 1;
                }
            }

            CP00 = 0;
            CP10 = 0;
            CP01 = 0;
            CP11 = 0;

            print_params();
            printf("rel = %le\n", HX / (-HY + 1.));
            printf("midIndexX = %d\n", midIndexX);
            printf("midIndexY = %d\n", midIndexY);
            printf("G1\n");
            print_vector(G1, NX_1);
            printf("G2\n");
            print_vector(G2, NY_1);
            printf("G3\n");
            print_vector(G3, NX_1);
            printf("G4\n");
            print_vector(G4, NY_1);

            double *exact0 = get_exact_solution_7(HX, HY, 0);
            //print_matrix_to_file( NX_1, NY_1, exact0, "/home/jane/tmp.txt");
            double *exactT = get_exact_solution_7(HX, HY, TAU * TIME_STEP_CNT);

            double *density = solve_7(tme);
            double *err = calc_error_7(HX, HY, TAU * TIME_STEP_CNT, density);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, density);
            print_surface("err", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, err);
            print_surface("exact", NX, NY, HX, HY, 0, A, C, x0, y0, TAU, U,
                          V, exact0);
            print_surface("exact", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                          V, exactT);

            double l1 = get_l1_norm_vec(NX_1, NY_1, err);
            double l_inf = get_l_inf_norm(NX_1, NY_1, err);

            printf("l1 %le \n", l1);
            printf("l_inf %le\n", l_inf);

            delete[] density;
            delete[] exact0;
            delete[] exactT;
            delete[] err;
        }
    }
}

// тестируем солвер 8 гладкая задача
// подход с выбором четверти
//
TEST_F(FemFixture, test8_1) {
    double tme = 0.;

    double d = 0;
    for (int i = 2; i < 3; ++i) {
        switch (i) {
            case 0:
                d = 50.;
                break;
            case 1:
                d = 100.;
                break;
            case 2:
                d = 200.;
                break;
            case 3:
                d = 400.;
                break;
            case 4:
                d = 800.;
                break;
            case 5:
                d = 1600.;
                break;
            default:
                return;
        }

        A = 0.;
        B = 1.;
        C = 0.;
        D = 1.;

        NX = (int) d;
        NY = (int) d;
        NX_1 = NX + 1;
        NY_1 = NY + 1;
        XY = NX_1 * NY_1;

        IDEAL_SQ_SIZE_X = 64;
        IDEAL_SQ_SIZE_Y = 64;

        HX = (B - A) / NX;
        HY = (D - C) / NY;
        TAU = 7e-4;
        TIME_STEP_CNT = 1;

        init_boundary_arrays_and_cp();


        for (int i = 0; i < NX_1; ++i) {
            G1[i] = 1;
        }
        for (int j = 0; j < NY_1; ++j) {
            G2[j] = 1;
        }
        for (int i = 0; i < NX_1; ++i) {
            G3[i] = 1;
        }
        for (int j = 0; j < NY_1; ++j) {
            G4[j] = 0;
        }

        CP00 = 0;
        CP10 = 1;
        CP01 = 0;
        CP11 = 1;

        print_params();
        printf("rel = %le\n", HX / (-HY + 1.));
        printf("G1\n");
        print_vector(G1, NX_1);
        printf("G2\n");
        print_vector(G2, NY_1);
        printf("G3\n");
        print_vector(G3, NX_1);
        printf("G4\n");
        print_vector(G4, NY_1);

        double *density = solve_8(tme);
        double *err = calc_error_8(HX, HY, TAU * TIME_STEP_CNT, density);
        double *exact0 = get_exact_solution_8(HX, HY, 0);
        double *exactT = get_exact_solution_8(HX, HY, TAU * TIME_STEP_CNT);

        double x0 = get_center_x();
        double y0 = get_center_y();
        print_surface("rho", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                      V, density);
        print_surface("err", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                      V, err);
        print_surface("exact", NX, NY, HX, HY, 0, A, C, x0, y0, TAU, U,
                      V, exact0);
        print_surface("exact", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
                      V, exactT);

        double l1 = get_l1_norm_vec(NX_1, NY_1, err);
        double l_inf = get_l_inf_norm(NX_1, NY_1, err);

        printf("l1 %le \n", l1);
        printf("l_inf %le\n", l_inf);

        delete[] density;
        delete[] exact0;
        delete[] exactT;
        delete[] err;
    }
}

// тестируем солвер 9 движение из угла в угол
// подход с выбором четверти и адаптивной сеткой
//
TEST_F(FemFixture, test9_1) {
    double tme = 0.;
    for (int iter = 0; iter < 1; ++iter) {

        double d = 0;
        for (int i = 1; i < 2; ++i) {
            switch (i) {
                case 0:
                    d = 50.;
                    break;
                case 1:
                    d = 100.;
                    break;
                case 2:
                    d = 200.;
                    break;
                case 3:
                    d = 400.;
                    break;
                case 4:
                    d = 800.;
                    break;
                case 5:
                    d = 1600.;
                    break;
                default:
                    return;
            }

            A = 0.;
            B = 1.;
            C = 0.;
            D = 1.;
            R_SQ = 0.099 * 0.099;
            INN_DENSITY = 1.;
            OUT_DENSITY = 0.;
            R_LVL = 1;

           // d=5.;

            int sz = (int) d;
            sz = sz * ((int) std::pow(3., R_LVL));

            NX = (int) d;
            NY = (int) d;
            NX_1 = NX + 1;
            NY_1 = NY + 1;
            HX = (B - A) / NX;
            HY = (D - C) / NY;
            IDEAL_SQ_SIZE_X = 64;
            IDEAL_SQ_SIZE_Y = 64;
            EPS_GRID = 0.5;
            RES_EPS = 1.e-9;

            NX3 = sz;
            NY3 = sz;
            NX3_1 = NX3 + 1;
            NY3_1 = NY3 + 1;
            R = (int) std::pow(3., R_LVL);


            CENTER_OFFSET_X = 0.3;
            CENTER_OFFSET_Y = 0.3;

            U = 1.;
            V = 1.;
            OMEGA = 1.;
            TAU = 4.9e-3;

            TIME_STEP_CNT = 1;
            XY = NX3_1 * NY3_1;

            init_boundary_arrays_and_cp();
            print_params();
            printf("NX3 = %d\n", NX3);
            printf("NX3_1 = %d\n", NX3_1);
            printf("NY3 = %d\n", NY3);
            printf("NY3_1 = %d\n", NY3_1);
            printf("XY = %d\n", XY);
            printf("R_LVL = %d\n", R_LVL);
            printf("R = %d\n", R);
            printf("EPS_GRID = %e\n", EPS_GRID);
            printf("RES_EPS = %e\n", RES_EPS);

            int* grid = new int[XY];
            int* gridPr = new int[XY];

            double *density = solve_9(tme, grid, gridPr);
//            double *err = calc_error_9(HX, HY, TAU * TIME_STEP_CNT, density, NX3_1, NY3_1);
            double *exact0 = get_exact_solution_9(HX, HY, 0, NX_1, NY_1);
            double *exactT = get_exact_solution_9(HX, HY, TAU * TIME_STEP_CNT, NX_1, NY_1);

            double x0 = get_center_x();
            double y0 = get_center_y();
//            print_surface("rho", NX3, NY3, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
//                          V, density);
//            print_surface("err", NX3, NY3, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U,
//                          V, err);
            print_surface("exact", NX, NY, HX, HY, 0, A, C, x0, y0, TAU, U, V, exact0);
            print_surface("exact", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U, V, exactT);

//            double l1 = get_l1_norm_vec(NX3_1, NY3_1, err);
//            double l_inf = get_l_inf_norm(NX3_1, NY3_1, err);
//            printf("l1 %le \n", l1);
//            printf("l_inf %le\n", l_inf);
//            delete[] density;
            delete[] exact0;
            delete[] exactT;
//            delete[] err;
            delete[] grid;
            delete[] gridPr;
        }
    }
}