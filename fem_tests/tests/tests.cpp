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
    printf("\nOX_LENxOY_LEN = %dx%d\n", OX_LEN, OY_LEN);
    printf("(U, V) = (%le, %le)\n", U_VELOCITY, V_VELOCITY);
    printf("(HX, HY) = (%le, %le)\n", HX, HY);
    printf("TAU = %le\n", TAU);
    printf("TIME_STEP_CNT = %d\n", TIME_STEP_CNT);
    printf("INTEGR_TYPE = %d\n", INTEGR_TYPE);
    printf("IDEAL_SQ_SIZE = %dx%d\n", IDEAL_SQ_SIZE_X, IDEAL_SQ_SIZE_Y);
    printf("CENTER_OFFSET = %le, %le\n", CENTER_OFFSET_X, CENTER_OFFSET_Y);
}

void init_boundary_arrays_and_cp() {
    G1 = new int[OX_LEN_1];
    G2 = new int[OY_LEN_1];
    G3 = new int[OX_LEN_1];
    G4 = new int[OY_LEN_1];
    for (int i = 0; i < OX_LEN_1; ++i) {
        G1[i] = 0;
        G3[i] = 0;
    }
    for (int j = 0; j < OY_LEN_1; ++j) {
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
    OX_LEN = 400;
    OY_LEN = 400;
    OX_LEN_1 = OX_LEN + 1;
    OY_LEN_1 = OY_LEN + 1;
    XY_LEN = OX_LEN_1 * OY_LEN_1;
    R_SQ = 0.099 * 0.099;
    INN_DENSITY = 1.;
    OUT_DENSITY = 0.;

    HX = (B - A) / OX_LEN;
    HY = (D - C) / OY_LEN;
    IDEAL_SQ_SIZE_X = 128 * (3 + 1);
    IDEAL_SQ_SIZE_Y = 128 * (3 + 1);

    CENTER_OFFSET_X = 0.5;
    CENTER_OFFSET_Y = 0.5;

    INTEGR_TYPE = 1;

    U_VELOCITY = 1.;
    V_VELOCITY = 1.;
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
        double *arr = getArr(skipLines, rhoFiles[i], XY_LEN);
        print_line_along_x("rho_x_1d", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                           U_VELOCITY, V_VELOCITY, arr, fixed_y);
        print_line_along_y("rho_y_1d", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                           U_VELOCITY, V_VELOCITY, arr, fixed_x);
        delete[] arr;
    }
    std::cout << "Starting err files convertsion..." << std::endl;
    for (unsigned int i = 0; i < errFiles.size(); i++) {
        std::cout << errFiles[i] << std::endl;
        int tl = get_tl(errFiles[i]);
        double *arr = getArr(skipLines, errFiles[i], XY_LEN);
        print_line_along_x("err_x_1d", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                           U_VELOCITY, V_VELOCITY, arr, fixed_y);
        print_line_along_y("err_x_1d", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                           U_VELOCITY, V_VELOCITY, arr, fixed_x);
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
        OX_LEN = (int) d;
        OY_LEN = (int) d;
        OX_LEN_1 = OX_LEN + 1;
        OY_LEN_1 = OY_LEN + 1;
        HX = (B - A) / OX_LEN;
        HY = (D - C) / OY_LEN;
        TIME_STEP_CNT = (int) d;
        XY_LEN = OX_LEN_1 * OY_LEN_1;
        U_VELOCITY = 0.;
        V_VELOCITY = 0.;
        CENTER_OFFSET_X = OX_LEN * HX / 2.;
        CENTER_OFFSET_Y = OY_LEN * HY / 2.;

        print_params();

        double *density = solve_1(tme);
        double *err = calc_error_1(HX, HY, density);
        double x0 = get_center_x();
        double y0 = get_center_y();
        print_surface("test1_rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                      V_VELOCITY, density);
        print_surface("test1_err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                      V_VELOCITY, err);
        double l1 = get_l1_norm(HX, HY, OX_LEN_1, OY_LEN_1, err);
        double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
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

        OX_LEN = (int) d;
        OY_LEN = (int) d;
        OX_LEN_1 = OX_LEN + 1;
        OY_LEN_1 = OY_LEN + 1;
        HX = (B - A) / OX_LEN;
        HY = (D - C) / OY_LEN;

        // (u, v) = (a, a)
        double a = 1;
        U_VELOCITY = a;
        V_VELOCITY = a;
        TAU = std::min(HX / (3. * a), HY / (3. * a));

        print_params();

        TIME_STEP_CNT = (int) 1;
        XY_LEN = OX_LEN_1 * OY_LEN_1;

        printf("OX_LEN = %d OY_LEN = %d\n", OX_LEN, OY_LEN);
        double *density = solve_2(tme);
        double *err = calc_error_2(HX, HY, TAU * TIME_STEP_CNT, density);
        double x0 = get_center_x();
        double y0 = get_center_y();
        print_surface("test2_rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                      V_VELOCITY, density);
        print_surface("test2_err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                      V_VELOCITY, err);
        double l1 = get_l1_norm(HX, HY, OX_LEN_1, OY_LEN_1, err);
        double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
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

    OX_LEN = (int) 400;
    OY_LEN = (int) 400;
    OX_LEN_1 = OX_LEN + 1;
    OY_LEN_1 = OY_LEN + 1;
    HX = (B - A) / OX_LEN;
    HY = (D - C) / OY_LEN;
    IDEAL_SQ_SIZE_X = 64;
    IDEAL_SQ_SIZE_Y = 64;
    CENTER_OFFSET_X = 0.3;
    CENTER_OFFSET_Y = 0.3;

    INTEGR_TYPE = 1;

    U_VELOCITY = 1.;
    V_VELOCITY = 1.;
    TAU = 1.e-3;
    //TIME_STEP_CNT = (int) ((1 - get_center_x_2() - get_center_y_2()) / TAU);
    TIME_STEP_CNT = 2;
    XY_LEN = OX_LEN_1 * OY_LEN_1;

    print_params();

    double *density = solve_2(tme);
    double *err = calc_error_2(HX, HY, TAU * TIME_STEP_CNT, density);
    double *exact0 = get_exact_solution_2(HX, HY, 0);
    double *exactT = get_exact_solution_2(HX, HY, TAU * TIME_STEP_CNT);

    double y0 = get_center_y();
    double x0 = get_center_x();
    print_surface("rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                  V_VELOCITY, density);
    print_surface("err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                  V_VELOCITY, err);
    print_surface("exact", OX_LEN, OY_LEN, HX, HY, 0, A, C, x0, y0, TAU, U_VELOCITY,
                  V_VELOCITY, exact0);
    print_surface("exact", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                  V_VELOCITY, exactT);
    double l1 = get_l1_norm(HX, HY, OX_LEN_1, OY_LEN_1, err);
    double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
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

            OX_LEN = (int) d;
            OY_LEN = (int) d;
            OX_LEN_1 = OX_LEN + 1;
            OY_LEN_1 = OY_LEN + 1;
            HX = (B - A) / OX_LEN;
            HY = (D - C) / OY_LEN;
            IDEAL_SQ_SIZE_X = 4;
            IDEAL_SQ_SIZE_Y = 4;
            CENTER_OFFSET_X = 0.3;
            CENTER_OFFSET_Y = 0.3;

            INTEGR_TYPE = 1;

            U_VELOCITY = 1.;
            V_VELOCITY = 1.;
            //TAU = 16. / pow(2., (i + 1));
            TAU = 1.2e-3;
            TIME_STEP_CNT = (int) pow(2., i);
            TIME_STEP_CNT = 20;
            XY_LEN = OX_LEN_1 * OY_LEN_1;

            print_params();

            double *density = solve_2(tme);
            double *err = calc_error_2(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_2(HX, HY, 0);
            double *exactT = get_exact_solution_2(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, density);
            print_surface("err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, err);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, 0, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exact0);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exactT);

            double l1 = get_l1_norm(HX, HY, OX_LEN_1, OY_LEN_1, err);
            double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
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

            OX_LEN = (int) d;
            OY_LEN = (int) d;
            OX_LEN_1 = OX_LEN + 1;
            OY_LEN_1 = OY_LEN + 1;
            HX = (B - A) / OX_LEN;
            HY = (D - C) / OY_LEN;
            IDEAL_SQ_SIZE_X = 128 * (iter + 1);
            IDEAL_SQ_SIZE_Y = 128 * (iter + 1);

            CENTER_OFFSET_X = 0.85;
            CENTER_OFFSET_Y = 0.5;

            INTEGR_TYPE = 1;

            U_VELOCITY = 1.;
            V_VELOCITY = 1.;
            TAU = 16. / pow(2., (i + 1));
            TAU *= 1.e-3;
            //TIME_STEP_CNT = (int) pow(2., i);
            TIME_STEP_CNT = 1;
            XY_LEN = OX_LEN_1 * OY_LEN_1;

            print_params();

            double *density = solve_2(tme);
            double *err = calc_error_2(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_2(HX, HY, 0);
            double *exactT = get_exact_solution_2(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, density);
            print_surface("err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, err);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, 0, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exact0);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exactT);

            double l1 = get_l1_norm(HX, HY, OX_LEN_1, OY_LEN_1, err);
            double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
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

            OX_LEN = (int) d;
            OY_LEN = (int) d;
            OX_LEN_1 = OX_LEN + 1;
            OY_LEN_1 = OY_LEN + 1;
            HX = (B - A) / OX_LEN;
            HY = (D - C) / OY_LEN;
            //IDEAL_SQ_SIZE_X = 128 * (iter + 1);
            //IDEAL_SQ_SIZE_Y = 128 * (iter + 1);


            CENTER_OFFSET_X = 0.3;
            CENTER_OFFSET_Y = 0.3;

            INTEGR_TYPE = 4;

            U_VELOCITY = 1.;
            V_VELOCITY = 1.;
//            TAU = HX;
            TAU = 0.001;

            //TIME_STEP_CNT = (int) pow(2., i);
            //TIME_STEP_CNT = 80;
            TIME_STEP_CNT = 1;
            //TIME_STEP_CNT = 1 + (0.7-CENTER_OFFSET_X)/TAU; // 4*2.5e-3 = 0.01 0.7-0.3=0.4/2.5e-3=160

            XY_LEN = OX_LEN_1 * OY_LEN_1;

            print_params();

            double *density = solve_2(tme);
            double *err = calc_error_2(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_2(HX, HY, 0);
            double *exactT = get_exact_solution_2(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, density);
            print_surface("err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, err);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, 0, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exact0);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exactT);

            double l1 = get_l1_norm_vec(OX_LEN_1, OY_LEN_1, err);
            printf("l1_err_vec %le \n", l1);
            l1 = get_l1_norm_int_trapezoidal(HX, HY, OX_LEN, OY_LEN, err); // note! a loop boundary
            double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
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

            OX_LEN = (int) d;
            OY_LEN = (int) d;
            OX_LEN_1 = OX_LEN + 1;
            OY_LEN_1 = OY_LEN + 1;
            HX = (B - A) / OX_LEN;
            HY = (D - C) / OY_LEN;
            IDEAL_SQ_SIZE_X = 128 * (iter + 1);
            IDEAL_SQ_SIZE_Y = 128 * (iter + 1);

            CENTER_OFFSET_X = 0.3;
            CENTER_OFFSET_Y = 0.3;

            INTEGR_TYPE = 1;

            U_VELOCITY = 1.;
            V_VELOCITY = 1.;
            TAU = HX/2.;
            //TIME_STEP_CNT = (int) pow(2., i);
            //TIME_STEP_CNT = 1;
            TIME_STEP_CNT = (0.7-CENTER_OFFSET_X)/TAU; // 4*2.5e-3 = 0.01 0.7-0.3=0.4/2.5e-3=160

            XY_LEN = OX_LEN_1 * OY_LEN_1;

            print_params();

            double *density = solve_2(tme);
            double *err = calc_error_2(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_2(HX, HY, 0);
            double *exactT = get_exact_solution_2(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, density);
            print_surface("err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, err);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, 0, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exact0);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exactT);

            double l1 = get_l1_norm(HX, HY, OX_LEN_1, OY_LEN_1, err);
            double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
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

            OX_LEN = (int) d;
            OY_LEN = (int) d;
            OX_LEN_1 = OX_LEN + 1;
            OY_LEN_1 = OY_LEN + 1;
            HX = (B - A) / OX_LEN;
            HY = (D - C) / OY_LEN;
            IDEAL_SQ_SIZE_X = 128 * (iter + 1);
            IDEAL_SQ_SIZE_Y = 128 * (iter + 1);

            CENTER_OFFSET_X = 0.3;
            CENTER_OFFSET_Y = 0.3;

            INTEGR_TYPE = 3;

            U_VELOCITY = 1.;
            V_VELOCITY = 1.;
            TAU = HX/2.;
            //TIME_STEP_CNT = (int) pow(2., i);
            //TIME_STEP_CNT = 1;
            TIME_STEP_CNT = (0.7-CENTER_OFFSET_X)/TAU; // 4*2.5e-3 = 0.01 0.7-0.3=0.4/2.5e-3=160

            XY_LEN = OX_LEN_1 * OY_LEN_1;

            print_params();

            double *density = solve_2(tme);
            double *err = calc_error_2(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_2(HX, HY, 0);
            double *exactT = get_exact_solution_2(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, density);
            print_surface("err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, err);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, 0, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exact0);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exactT);

            double l1 = get_l1_norm(HX, HY, OX_LEN_1, OY_LEN_1, err);
            double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
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

            OX_LEN = (int) d;
            OY_LEN = (int) d;
            OX_LEN_1 = OX_LEN + 1;
            OY_LEN_1 = OY_LEN + 1;
            HX = (B - A) / OX_LEN;
            HY = (D - C) / OY_LEN;
            IDEAL_SQ_SIZE_X = 128 * (iter + 1);
            IDEAL_SQ_SIZE_Y = 128 * (iter + 1);

            CENTER_OFFSET_X = 0.5;
            CENTER_OFFSET_Y = 0.5;

            INTEGR_TYPE = 1;

            U_VELOCITY = 1.;
            V_VELOCITY = 1.;
            OMEGA = 1.;
            TAU = 2.5e-3;

            TIME_STEP_CNT = 2700;
            XY_LEN = OX_LEN_1 * OY_LEN_1;

            init_boundary_arrays_and_cp();

            int midIndexX = OX_LEN_1 / 2;
            int midIndexY = OY_LEN_1 / 2;
            for (int i = 0; i < OX_LEN_1; ++i) {
                if (i * HX + A < .5 && i > 0)
                    G1[i] = 1;
            }
            for (int j = 0; j < OY_LEN_1; ++j) {
                if (j * HY + C < .5 && j > 0)
                    G2[j] = 1;
            }
            for (int i = 0; i < OX_LEN_1; ++i) {
                if (i * HX + A > .5 && i < OX_LEN_1 - 1)
                    G3[i] = 1;
            }
            for (int j = 0; j < OY_LEN_1; ++j) {
                if (j * HY + C > .5 && j > 0 && j < OX_LEN_1 - 1) {
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
            print_vector(G1, OX_LEN_1);
            printf("G2\n");
            print_vector(G2, OY_LEN_1);
            printf("G3\n");
            print_vector(G3, OX_LEN_1);
            printf("G4\n");
            print_vector(G4, OY_LEN_1);

            double *density = solve_3(tme);
            double *err = calc_error_3(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_3(HX, HY, 0);
            double *exactT = get_exact_solution_3(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, density);
            print_surface("err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, err);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, 0, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exact0);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exactT);

            double l1 = get_l1_norm_vec(OX_LEN_1, OY_LEN_1, err);
            //           double l1 = get_l1_norm_int_middle(HX, HY, OX_LEN_1, OY_LEN_1, err);
            double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
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

            OX_LEN = (int) d;
            OY_LEN = (int) d;
            OX_LEN_1 = OX_LEN + 1;
            OY_LEN_1 = OY_LEN + 1;
            HX = (B - A) / OX_LEN;
            HY = (D - C) / OY_LEN;
            IDEAL_SQ_SIZE_X = 128 * (iter + 1);
            IDEAL_SQ_SIZE_Y = 128 * (iter + 1);

            CENTER_OFFSET_X = 0.5;
            CENTER_OFFSET_Y = 0.5;

            INTEGR_TYPE = 1;

            U_VELOCITY = 1.;
            V_VELOCITY = 1.;
            OMEGA = 1.;
            TAU = 2.5e-3;

            TIME_STEP_CNT = 50;
            XY_LEN = OX_LEN_1 * OY_LEN_1;

            double x0 = get_center_x();
            double y0 = get_center_y();

            double *exact0 = get_exact_solution_3(HX, HY, 0.);
            double *exactT = get_exact_solution_3(HX, HY, TAU * TIME_STEP_CNT);
            for (int j = 0; j < TIME_STEP_CNT; ++j) {
                double *exact0 = get_exact_solution_3(HX, HY, j * TAU);
                print_surface("exact", OX_LEN, OY_LEN, HX, HY, j, A, C, x0, y0, TAU, U_VELOCITY,
                              V_VELOCITY, exact0);
            }

            print_surface("exact", OX_LEN, OY_LEN, HX, HY, 0, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exact0);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exactT);


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

            OX_LEN = (int) d;
            OY_LEN = (int) d;
            OX_LEN_1 = OX_LEN + 1;
            OY_LEN_1 = OY_LEN + 1;
            HX = (B - A) / OX_LEN;
            HY = (D - C) / OY_LEN;
            IDEAL_SQ_SIZE_X = 128 * (iter + 1);
            IDEAL_SQ_SIZE_Y = 128 * (iter + 1);

            CENTER_OFFSET_X = 0.5;
            CENTER_OFFSET_Y = 0.7;

            INTEGR_TYPE = 1;

            U_VELOCITY = 1.;
            V_VELOCITY = 1.;
            OMEGA = 1.;

            TAU = 16. / pow(2., (i + 1));
            TAU *= 1.e-3;

            TIME_STEP_CNT = (int) pow(2., i);
            XY_LEN = OX_LEN_1 * OY_LEN_1;

            init_boundary_arrays_and_cp();

            int midIndexX = OX_LEN_1 / 2;
            int midIndexY = OY_LEN_1 / 2;
            for (int i = 0; i < OX_LEN_1; ++i) {
                if (i * HX + A < .5 && i > 0)
                    G1[i] = 1;
            }
            for (int j = 0; j < OY_LEN_1; ++j) {
                if (j * HY + C < .5 && j > 0)
                    G2[j] = 1;
            }
            for (int i = 0; i < OX_LEN_1; ++i) {
                if (i * HX + A > .5 && i < OX_LEN_1 - 1)
                    G3[i] = 1;
            }
            for (int j = 0; j < OY_LEN_1; ++j) {
                if (j * HY + C > .5 && j > 0 && j < OX_LEN_1 - 1) {
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
            print_vector(G1, OX_LEN_1);
            printf("G2\n");
            print_vector(G2, OY_LEN_1);
            printf("G3\n");
            print_vector(G3, OX_LEN_1);
            printf("G4\n");
            print_vector(G4, OY_LEN_1);

            double *density = solve_3(tme);
            double *err = calc_error_3(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_3(HX, HY, 0);
            double *exactT = get_exact_solution_3(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, density);
            print_surface("err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, err);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, 0, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exact0);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exactT);

            double l1 = get_l1_norm(HX, HY, OX_LEN_1, OY_LEN_1, err);
            double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
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

            OX_LEN = (int) d;
            OY_LEN = (int) d;
            OX_LEN_1 = OX_LEN + 1;
            OY_LEN_1 = OY_LEN + 1;
            HX = (B - A) / OX_LEN;
            HY = (D - C) / OY_LEN;
            IDEAL_SQ_SIZE_X = 32 * (iter + 1);
            IDEAL_SQ_SIZE_Y = 32 * (iter + 1);

            CENTER_OFFSET_X = 0.3;
            CENTER_OFFSET_Y = 0.3;

            INTEGR_TYPE = 1;

            U_VELOCITY = 1.;
            V_VELOCITY = 1.;

            //TAU = 16. / pow(2., (i + 1));
            TAU = 1.e-3;

            //TIME_STEP_CNT = (int) pow(2., i);
            TIME_STEP_CNT = 2;
            XY_LEN = OX_LEN_1 * OY_LEN_1;

            init_boundary_arrays_and_cp();

            int midIndexX = OX_LEN_1 / 2;
            int midIndexY = OY_LEN_1 / 2;
            for (int i = 1; i < OX_LEN; ++i) {
                G1[i] = 0;
                G3[i] = 1;
            }
            for (int j = 1; j < OY_LEN; ++j) {
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
            print_vector(G1, OX_LEN_1);
            printf("G2\n");
            print_vector(G2, OY_LEN_1);
            printf("G3\n");
            print_vector(G3, OX_LEN_1);
            printf("G4\n");
            print_vector(G4, OY_LEN_1);

            double *density = solve_4(tme);
            double *err = calc_error_4(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_4(HX, HY, 0);
            double *exactT = get_exact_solution_4(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, density);
            print_surface("err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, err);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, 0, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exact0);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exactT);

            double l1 = get_l1_norm(HX, HY, OX_LEN_1, OY_LEN_1, err);
            double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
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

            OX_LEN = (int) d;
            OY_LEN = (int) d;
            OX_LEN_1 = OX_LEN + 1;
            OY_LEN_1 = OY_LEN + 1;
            HX = (B - A) / OX_LEN;
            HY = (D - C) / OY_LEN;
            IDEAL_SQ_SIZE_X = 64;
            IDEAL_SQ_SIZE_Y = 64;

            CENTER_OFFSET_X = 0.3;
            CENTER_OFFSET_Y = 0.3;

            INTEGR_TYPE = 1;

            U_VELOCITY = 1.;
            V_VELOCITY = 1.;
            OMEGA = 1.;
            TAU = 1.2475e-3;

            TIME_STEP_CNT = 322;
            XY_LEN = OX_LEN_1 * OY_LEN_1;

            init_boundary_arrays_and_cp();

            int midIndexX = OX_LEN_1 / 2;
            int midIndexY = OY_LEN_1 / 2;

            for (int i = 0; i < OX_LEN_1; ++i) {
                if (i * HX + A < .5 && i > 0)
                    G1[i] = 0;
            }
            for (int j = 0; j < OY_LEN_1; ++j) {
                if (j * HY + C < .5 && j > 0)
                    G2[j] = 0;
            }
            for (int i = 0; i < OX_LEN_1; ++i) {
                if (i * HX + A > .5 && i < OX_LEN_1 - 1)
                    G3[i] = 0;
            }
            for (int j = 0; j < OY_LEN_1; ++j) {
                if (j * HY + C > .5 && j > 0 && j < OX_LEN_1 - 1) {
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
            print_surface("rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, density);
            print_surface("err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, err);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, 0, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exact0);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exactT);

            double l1 = get_l1_norm_vec(OX_LEN_1, OY_LEN_1, err);
            double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
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

            OX_LEN = (int) d;
            OY_LEN = (int) d;
            OX_LEN_1 = OX_LEN + 1;
            OY_LEN_1 = OY_LEN + 1;
            HX = (B - A) / OX_LEN;
            HY = (D - C) / OY_LEN;
            IDEAL_SQ_SIZE_X = 64;
            IDEAL_SQ_SIZE_Y = 64;

            CENTER_OFFSET_X = 0.3;
            CENTER_OFFSET_Y = 0.3;

            INTEGR_TYPE = 1;

            U_VELOCITY = 1.;
            V_VELOCITY = 1.;
            OMEGA = 1.;
            TAU = 1.2475e-3;

            TIME_STEP_CNT = 322;
            XY_LEN = OX_LEN_1 * OY_LEN_1;

            init_boundary_arrays_and_cp();

            int midIndexX = OX_LEN_1 / 2;
            int midIndexY = OY_LEN_1 / 2;

            for (int i = 0; i < OX_LEN_1; ++i) {
                if (i * HX + A < .5 && i > 0)
                    G1[i] = 0;
            }
            for (int j = 0; j < OY_LEN_1; ++j) {
                if (j * HY + C < .5 && j > 0)
                    G2[j] = 1;
            }
            for (int i = 0; i < OX_LEN_1; ++i) {
                if (i * HX + A > .5 && i < OX_LEN_1 - 1)
                    G3[i] = 1;
            }
            for (int j = 0; j < OY_LEN_1; ++j) {
                if (j * HY + C > .5 && j > 0 && j < OX_LEN_1 - 1) {
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
            print_surface("rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, density);
            print_surface("err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, err);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, 0, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exact0);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exactT);

            double l1 = get_l1_norm_vec(OX_LEN_1, OY_LEN_1, err);
            double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
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

            OX_LEN = (int) d;
            OY_LEN = (int) d;
            OX_LEN_1 = OX_LEN + 1;
            OY_LEN_1 = OY_LEN + 1;
            HX = (B - A) / OX_LEN;
            HY = (D - C) / OY_LEN;
            IDEAL_SQ_SIZE_X = 64;
            IDEAL_SQ_SIZE_Y = 64;

            CENTER_OFFSET_X = 0.3;
            CENTER_OFFSET_Y = 0.3;

            INTEGR_TYPE = 1;

            U_VELOCITY = 1.;
            V_VELOCITY = 1.;
            OMEGA = 1.;
            TAU = 1.2475e-3;

            TIME_STEP_CNT = 322;
            XY_LEN = OX_LEN_1 * OY_LEN_1;

            init_boundary_arrays_and_cp();

            int midIndexX = OX_LEN_1 / 2;
            int midIndexY = OY_LEN_1 / 2;

            for (int i = 0; i < OX_LEN_1; ++i) {
                if (i * HX + A < .5 && i > 0)
                    G1[i] = 0;
            }
            for (int j = 0; j < OY_LEN_1; ++j) {
                if (j * HY + C < .5 && j > 0)
                    G2[j] = 1;
            }
            for (int i = 0; i < OX_LEN_1; ++i) {
                if (i * HX + A > .5 && i < OX_LEN_1 - 1)
                    G3[i] = 1;
            }
            for (int j = 0; j < OY_LEN_1; ++j) {
                if (j * HY + C > .5 && j > 0 && j < OX_LEN_1 - 1) {
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

            double *density = solve_7(tme);
            double *err = calc_error_7(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_7(HX, HY, 0);
            double *exactT = get_exact_solution_7(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, density);
            print_surface("err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, err);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, 0, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exact0);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exactT);

            double l1 = get_l1_norm_vec(OX_LEN_1, OY_LEN_1, err);
            double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);

            printf("l1 %le \n", l1);
            printf("l_inf %le\n", l_inf);

            delete[] density;
            delete[] exact0;
            delete[] exactT;
            delete[] err;
        }
    }
}