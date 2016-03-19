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
int OX_LEN = 0;
int OY_LEN = 0;
int OX_LEN_1 = 0;
int OY_LEN_1 = 0;
int XY_LEN = 0;
double TAU = 0.;
int TIME_STEP_CNT = 0;
int JAK_ITER_CNT = 0;
double HX = 0.;
double HY = 0.;
double R_SQ = 0.;
double INN_DENSITY = 0.;
double OUT_DENSITY = 0.;
double U_VELOCITY = 0.;
double V_VELOCITY = 0.;
double DBL_MIN_TRIM = 0.;
double RES_EPS = 0.;
int INTEGR_TYPE = 0;
int IDEAL_SQ_SIZE_X = 0;
int IDEAL_SQ_SIZE_Y = 0;

