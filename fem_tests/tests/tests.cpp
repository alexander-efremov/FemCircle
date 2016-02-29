#include <utils.h>
#include "gtest/gtest.h"
#include "consts.h"
#include "solver1.h"
#include "solver2.h"
#include "solver3.h"
#include "tecplot.h"

class FemFixture : public ::testing::Test {
protected:
    virtual void TearDown() {
    }

    virtual void SetUp() {
    }

public:
    FemFixture() : Test() {
    }
};

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
        printf("\nOX_LEN = %d OY_LEN = %d\n", OX_LEN, OY_LEN);
        double *density = solve_1(tme);
        double *err = calc_error_1(HX, HY, density);
        print_surface_as_v("test1_rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, density);
        print_surface_as_v("test1_err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, err);
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
        printf("\nOX_LEN = %d OY_LEN = %d\n", OX_LEN, OY_LEN);
        double *density = solve_2(tme);
        double *err = calc_error_2(HX, HY, density);
        print_surface_as_v("test2_rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, density);
        print_surface_as_v("test2_err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, err);
        double l1 = get_l1_norm(HX, HY, OX_LEN_1, OY_LEN_1, err);
        double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
        printf("l1 %le \n", l1);
        printf("l_inf %le\n", l_inf);
        delete[] density;
        delete[] err;
    }
}

TEST_F(FemFixture, test3) {
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
        printf("\nOX_LEN = %d OY_LEN = %d\n", OX_LEN, OY_LEN);
        double *density = solve_3(tme);
        double *err = calc_error_3(HX, HY, density);
        print_surface_as_v("test3_rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, density);
        print_surface_as_v("test3_err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, err);
        double l1 = get_l1_norm(HX, HY, OX_LEN_1, OY_LEN_1, err);
        double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
        printf("l1 %le \n", l1);
        printf("l_inf %le\n", l_inf);
        delete[] density;
        delete[] err;
    }
}
