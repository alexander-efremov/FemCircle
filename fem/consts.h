#ifndef FEM_CIRCLE_CONSTS_H
#define FEM_CIRCLE_CONSTS_H

// area is [AxB]x[CxD]
extern double A;
extern double B;
extern double C;
extern double D;
extern int OX_LEN;
extern int OY_LEN;
extern int OX_LEN_1;
extern int OY_LEN_1;
extern int XY_LEN;
extern double TAU;
extern int TIME_STEP_CNT;
extern int JAK_ITER_CNT;
extern double HX;
extern double HY;
extern double R_SQ; // radius of circle in second power
extern double INN_DENSITY; // density inside circle
extern double OUT_DENSITY; // density out of circle boundary

#define EPS 10e-8

#endif //FEM_CIRCLE_CONSTS_H
