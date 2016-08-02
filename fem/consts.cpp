#include "consts.h"
#include <algorithm>

#ifdef WIN32
__int64 timerStart = 0;
#else
struct timeval timerStart;
#endif

double A = 0.;
double B = 0.;
double C = 0.;
double D = 0.;
int NX = 0;
int NY = 0;
int NX_1 = 0;
int NY_1 = 0;
int NX3 = 0;
int NY3 = 0;
int NX3_1 = 0;
int NY3_1 = 0;
int XY = 0;
double TAU = 0.;
int TIME_STEP_CNT = 0;
int JAK_ITER_CNT = 0;
double HX = 0.;
double HY = 0.;
double R_SQ = 0.;
double INN_DENSITY = 0.;
double OUT_DENSITY = 0.;
double U = 0.;
double V = 0.;
double DBL_MIN_TRIM = 1.e-16;
double RES_EPS = 1.e-14;
double EPS = 1.e-8;
int INTEGR_TYPE = 1;
int IDEAL_SQ_SIZE_X = 0;
int IDEAL_SQ_SIZE_Y = 0;
double CENTER_OFFSET_X = 0.;
double CENTER_OFFSET_Y = 0.;
double OMEGA = 0.;
int* G1 = NULL;
int* G2 = NULL;
int* G3 = NULL;
int* G4 = NULL;
int CP00 = 0;
int CP10 = 0;
int CP11 = 0;
int CP01 = 0;
int R_LVL = 0;
int R = 0;
double EPS_GRID = 0.;
double HX_SMALLEST = 0.;
double HY_SMALLEST = 0.;
int APPROX_TYPE  = 1;