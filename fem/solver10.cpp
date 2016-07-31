#include <math.h>#include <stdio.h>#include <string.h>#include <assert.h>#include "consts.h"#include "timer.h"#include "utils.h"#include "common.h"#include "algorithm"using namespace std;// только грубая сеткаinline void print_surface_0(const char *filename, int ox_len, int oy_len, int ny3_1, int r,                            double hx, double hy, int t, double a, double c, double x0, double y0,                            double tau, double u, double v, double *data) {    char name[650];    sprintf(name, "%s_nx=%d_ny=%d_hx=%f_hy=%f_t=%d_x0=%f_y0=%f_tau=%f_u=%f_v=%f_a=%f_c=%f.dat",            filename, ox_len + 1, oy_len + 1, hx, hy, t, x0, y0, tau, u, v, a, c);    FILE *file = fopen(name, "w");    fprintf(file, "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = 'x' 'y' %s\nZONE T='SubZone'",            filename);    fprintf(file, "\nI=%d J=%d K=%d ZONETYPE=Ordered", oy_len + 1, ox_len + 1, 1);    fprintf(file, "\nDATAPACKING=POINT\nDT=(SINGLE SINGLE SINGLE)");    for (int i = 0; i < ox_len + 1; i++)        for (int j = 0; j < oy_len + 1; j++)            fprintf(file, "\n%-30.20g  %-30.20g %-30.20g", i * hx, j * hy,                    data[ny3_1 * i * r + j * r]);    fclose(file);}// вывод ВСЕЙ сетки уровнями вложенностиinline void print_surface_0(const char *filename, int ox_len, int oy_len, int ny3_1,                            double hx, double hy, int t, double a, double c, double x0, double y0,                            double tau, double u, double v, int *data) {    char name[650];    sprintf(name, "%s_nx=%d_ny=%d_hx=%f_hy=%f_t=%d_x0=%f_y0=%f_tau=%f_u=%f_v=%f_a=%f_c=%f.dat",            filename, ox_len + 1, oy_len + 1, hx, hy, t, x0, y0, tau, u, v, a, c);    FILE *file = fopen(name, "w");    fprintf(file, "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = 'x' 'y' %s\nZONE T='SubZone'",            filename);    fprintf(file, "\nI=%d J=%d K=%d ZONETYPE=Ordered", oy_len + 1, ox_len + 1, 1);    fprintf(file, "\nDATAPACKING=POINT\nDT=(SINGLE SINGLE SINGLE)");    for (int i = 0; i < ox_len + 1; i++)        for (int j = 0; j < oy_len + 1; j++)            fprintf(file, "\n%-30.20g  %-30.20g %d", i * hx, j * hy,                    data[ny3_1 * i + j]);    fclose(file);}// вывод расчетной функции ВО ВСЕХ ТОЧКАХinline void print_surface_0(const char *filename, int ox_len, int oy_len, int ny3_1,                            double hx, double hy, int t, double a, double c, double x0, double y0,                            double tau, double u, double v, int *data, double *values) {    char name[650];    sprintf(name, "%s_nx=%d_ny=%d_hx=%f_hy=%f_t=%d_x0=%f_y0=%f_tau=%f_u=%f_v=%f_a=%f_c=%f.dat",            filename, ox_len + 1, oy_len + 1, hx, hy, t, x0, y0, tau, u, v, a, c);    FILE *file = fopen(name, "w");    fprintf(file, "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = 'x' 'y' %s\nZONE T='SubZone'",            filename);    fprintf(file, "\nI=%d J=%d K=%d ZONETYPE=Ordered", oy_len + 1, ox_len + 1, 1);    fprintf(file, "\nDATAPACKING=POINT\nDT=(SINGLE SINGLE SINGLE)");    for (int i = 0; i < ox_len + 1; i++) {        for (int j = 0; j < oy_len + 1; j++) {            if (data[ny3_1 * i + j] > -1) {                fprintf(file, "\n%-30.20g  %-30.20g %-30.20g", i * hx, j * hy,                        values[ny3_1 * i + j]);            }            else                fprintf(file, "\n%-30.20g  %-30.20g %-30.20g", i * hx, j * hy,                        values[ny3_1 * i + j]);        }    }    fclose(file);}// вывод расчетной функции в РАСЧЕТНЫХ ТОЧКАХ, во всех остальных выводится уровень вложенностиinline void print_surface_1(const char *filename, int ox_len, int oy_len, int ny3_1,                            double hx, double hy, int t, double a, double c, double x0, double y0,                            double tau, double u, double v, int *data, double *values) {    char name[650];    sprintf(name, "%s_nx=%d_ny=%d_hx=%f_hy=%f_t=%d_x0=%f_y0=%f_tau=%f_u=%f_v=%f_a=%f_c=%f.dat",            filename, ox_len + 1, oy_len + 1, hx, hy, t, x0, y0, tau, u, v, a, c);    FILE *file = fopen(name, "w");    fprintf(file, "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = 'x' 'y' %s\nZONE T='SubZone'",            filename);    fprintf(file, "\nI=%d J=%d K=%d ZONETYPE=Ordered", oy_len + 1, ox_len + 1, 1);    fprintf(file, "\nDATAPACKING=POINT\nDT=(SINGLE SINGLE SINGLE)");    for (int i = 0; i < ox_len + 1; i++) {        for (int j = 0; j < oy_len + 1; j++) {            if (data[ny3_1 * i + j] > -1) {                fprintf(file, "\n%-30.20g  %-30.20g %-30.20g", i * hx, j * hy,                        values[ny3_1 * i + j]);            }            else                fprintf(file, "\n%-30.20g  %-30.20g %d", i * hx, j * hy,                        data[ny3_1 * i + j]);        }    }}inline static double func_u(double t, double x, double y) { return U; }inline static double func_v(double t, double x, double y) { return V; }inline static double analytical_solution_circle(double t, double x, double y) {    double x0 = get_center_x() + t * func_u(t, x, y);    double y0 = get_center_y() + t * func_v(t, x, y);    double value = (x - x0) * (x - x0) + (y - y0) * (y - y0);    if (value <= R_SQ) return INN_DENSITY;    return OUT_DENSITY;}// определение количества точек по направлению в зависимости от уровняinline int get_n3_1(int n, int lvl) { // не обязательно max_lvl == MAX_LVL, но важно чтобы n == NX || n == NY    assert(n == NX || n == NY);    return (int) (n * std::pow(3., lvl) + 1);}// получаем самый мелкий шагinline void get_hx_hy_smallest(int nx, int ny, int max_lvl, double &hx,                               double &hy) {//!!ВАЖНО: nx,ny,max_lvl соответствовало NX,NY,MAX_LVL    assert(nx == NX);    assert(ny == NY);    assert(max_lvl == R_LVL);    hx = (B - A) / (nx * std::pow(3., max_lvl));    hy = (D - C) / (ny * std::pow(3., max_lvl));}// 1) знаем i-абсолютное, j-абсолютное -> найти lev, hinline void get_hx_hy(int i, int j, int *grid, int nx, int ny,                      int max_lvl, //!!ВАЖНО: nx,ny,max_lvl соответствовало NX,NY,MAX_LVL                      double &hx, double &hy, int &lev) {    assert(nx == NX);    assert(ny == NY);    assert(max_lvl == R_LVL);    int ny3_1 = get_n3_1(ny, max_lvl);    lev = grid[i * ny3_1 + j];    if (lev != -1) {        hx = (B - A) / (nx * std::pow(3., lev));        hy = (D - C) / (ny * std::pow(3., lev));    }    else {        hx = 0.;        hy = 0.;    }}// 1) знаем i-абсолютное, j-абсолютное -> найти lev, hinline void get_hx_hy(int i, int j, int *grid, int nx, int ny, int max_lvl, double &hx, double &hy) {    int lv = 0;    get_hx_hy(i, j, grid, nx, ny, max_lvl, hx, hy, lv);}inline static double get_h(int n, int lev, double left_bound, double right_bound) {    assert(n == NX || n == NY);    double h = (left_bound - right_bound) / (n * std::pow(3., lev));    return h;}inline static int get_p(double x, int lev, int nx, double a, double b) { //!!ВАЖНО: nx соответствовало NX    assert(nx == NX);    assert(a == A);    assert(b == B);    double hx = get_h(nx, lev, b, a);    int p = (int) ((x + (hx / 2.)) / hx);    return p;}inline static int get_q(double y, int lev, int ny, double c, double d) {    assert(ny == NY);    assert(c == C);    assert(d == D);    double hy = get_h(ny, lev, d, c);    int q = (int) ((y + (hy / 2.)) / hy);    return q;}// 2) знаем x,y -> найти i-абсолютное, j-абсолютное, lev, hstatic double num_sol_at_point(double *rhoPr, int *gridPr, double x, double y, int max_lvl,                               int nx, int ny, int &lev, // что такое здесь nx, ny, lev?                               int &p,                               int &q) {    if (x < A || x > B || y < C || y > D) {        printf("\nnum_sol_at_point ERROR out if bounds x=%e y=%e\n", x, y);        return -2.;//!!Временно для отладки, Саша имел ввиду, наверное, что если из-за границы. то прет 0, но я хочу прка видеть где выход за границу был...    }    assert(nx == NX);    assert(ny == NY);    assert(max_lvl == R_LVL);    assert(lev == 0);    assert(p == 0);    assert(q == 0);    //int ny3_1 = get_n3_1(ny, max_lvl); // потом вернуться к этому...    int ny3_1 = NY3_1;    int d;    do {        p = get_p(x, lev, nx, A, B);        q = get_q(y, lev, ny, C, D);        int coef = (int) std::pow(3, max_lvl - lev);        d = gridPr[p * coef * ny3_1 + q * coef];        if (d == lev) { return rhoPr[p * coef * ny3_1 + q * coef]; }        lev++;    } while (lev <= max_lvl);    printf("\nnum_sol_at_point ERROR RHO MINUS ONE %e %e\n", x, y);    return -1.;}static double calc_phi(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4,                       double *rhoPr, double time, int *gridPr, int nx, int ny, double tau, int nx_ideal, int ny_ideal,                       int max_lvl) {//    printf("POINT: %d   %d :  x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le **//                   x4=%.8le * y4=%.8le\n", ii,jj, x1,y1, x2,y2, x3,y3, x4,y4);    assert(nx == NX);    assert(ny == NY);    assert(max_lvl == R_LVL);    double u = func_u(time, x1, y1);    double v = func_v(time, x1, y1);    x1 = x1 - tau * u;    y1 = y1 - tau * v;    u = func_u(time, x2, y2);    v = func_v(time, x2, y2);    x2 = x2 - tau * u;    y2 = y2 - tau * v;    u = func_u(time, x3, y3);    v = func_v(time, x3, y3);    x3 = x3 - tau * u;    y3 = y3 - tau * v;    u = func_u(time, x4, y4);    v = func_v(time, x4, y4);    x4 = x4 - tau * u;    y4 = y4 - tau * v;    /*     if (x1 <= A || x1 >= B || x2 <= A || x2 >= B || x3 <= A || x3 >= B || x4 <= A || x4 >= B        || y1 <= C || y1 >= D || y2 <= C || y2 >= D || y3 <= C || y3 >= D || y4 <= C || y4 >= D)        printf("PREV Time level %.8le! ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "                       "x4=%.8le * y4=%.8le\n ", time, ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);*/    double x_step = 1. / nx_ideal;    double y_step = 1. / ny_ideal;    double mes = x_step * y_step;    double phi = 0.;    for (int i = 0; i < nx_ideal; ++i) {        for (int j = 0; j < ny_ideal; ++j) {            double ideal_x = i * x_step + x_step / 2.;            double ideal_y = j * y_step + y_step / 2.;            double real_x = x1 + (x2 - x1) * ideal_x + (x4 - x1) * ideal_y                            + (x1 + x3 - x2 - x4) * ideal_x * ideal_y;            double real_y = y1 + (y2 - y1) * ideal_x + (y4 - y1) * ideal_y                            + (y1 + y3 - y2 - y4) * ideal_x * ideal_y;            double a11 = (x2 - x1) + (x1 + x3 - x2 - x4) * ideal_y;            double a12 = (x4 - x1) + (x1 + x3 - x2 - x4) * ideal_x;            double a21 = (y2 - y1) + (y1 + y3 - y2 - y4) * ideal_y;            double a22 = (y4 - y1) + (y1 + y3 - y2 - y4) * ideal_x;            double jakob = a11 * a22 - a21 * a12;            int lev = 0, p = 0, q = 0;            double dens = num_sol_at_point(rhoPr, gridPr, real_x, real_y, max_lvl, nx, ny, lev, p, q);            phi += mes * dens * jakob;        }    }    return fabs(phi) < fabs(DBL_MIN_TRIM) ? 0. : phi;}static double get_numeric_solution(double x, double y, double hx_lev, double hy_lev, double *rhoPr, int *gridPr,                                   int tl) {    double x1 = 0., y1 = 0., x2 = 0., y2 = 0., x3 = 0., y3 = 0., x4 = 0., y4 = 0.;    get_inner_area(x, y, x1, y1, x2, y2, x3, y3, x4, y4, hx_lev, hy_lev, A, B, C, D);    double val = calc_phi(x1, y1, x2, y2, x3, y3, x4, y4, rhoPr, TAU * tl, gridPr, NX, NY,                          TAU, IDEAL_SQ_SIZE_X, IDEAL_SQ_SIZE_Y, R_LVL);    double coef = 1. / (hx_lev * hy_lev);    val = coef * val;    return val;}static void grid_vert_checking(int lx, int rx, int by, int uy, int lev, int *grid, double *rho, double *rhoPr,                               int *gridPr, int tl);static void grid_hor_checking(int lx, int rx, int by, int uy, int lev, int *grid, double *rho, double *rhoPr,                              int *gridPr, int tl);static void grid_hor_refine(int ii, int jj, int lev, int *grid, double *rho, double *rhoPr, int *gridPr, int tl) {    assert(lev + 1 <= R_LVL);    int coef = (int) std::pow(3., R_LVL - lev);    int coefNew = (int) std::pow(3., R_LVL - (lev + 1));    // a) если обе окрестности еще ни разу не мельчились Е.Д.: нет, это если обе окрестности на одном уровне вложенности lev - это важно для потом, пока пропустим и вернемся    if (grid[ii * coef * NY3_1 + jj * coef] == lev && grid[(ii + 1) * coef * NY3_1 + jj * coef] == lev) {        // изменим обе окрестности        for (int p = -1; p < 5; ++p) { // p = -1, 0, 1, 2, 3, 4  Е.Д.: это по x изменение изменения индексов            // ок, я оставлю комменты пока, в скайпе в две минуты объясню потом            for (int q = -1; q < 2; ++q) { // q = -1, 0, 1 Е.Д.: а это по y изменение изменения индексов,                // т.к. мельчим две клетки вправо и только текущую по y                //смысл в том что мельчим от ii jj две ячейки вправо по x поэтому индексов больше и только текущую по y                int i = (3 * ii + p) * coefNew * NY3_1 +                        (3 * jj + q) * coefNew;// вернула назад, но не уверена по-прежнему                grid[i] = lev + 1; // изменим сетку                double hx_lev = HX_SMALLEST * std::pow(3., R_LVL - (lev + 1));                double hy_lev = HX_SMALLEST * std::pow(3., R_LVL - (lev + 1));                double x = (3 * ii + p) * coefNew * HX_SMALLEST;                double y = (3 * jj + q) * coefNew * HY_SMALLEST;                // ладно, оставим пока, я потом подумаю, может я перемудрила                rho[i] = get_numeric_solution(x, y, hx_lev, hy_lev, rhoPr, gridPr, tl); // вычислим значение            }        }        // осуществим проверку для вновь созданных ячеек (окрестностей)        grid_hor_checking(                (3 * ii - 1), // остался один вопрос: когда начинаем передаем NX, NY не умножая на 3^R_LVL,                // а здесь умножаем... может надо умножить там?                (3 * ii + 4),                (3 * jj - 1),                (3 * jj + 1),                lev + 1,                grid, rho, rhoPr, gridPr, tl);        grid_vert_checking((3 * ii - 1),                           (3 * ii + 4),                           (3 * jj - 1),                           (3 * jj + 1),                           lev + 1,                           grid, rho, rhoPr, gridPr, tl);        return; // как-то не уютно, когда реккурсия ничего не возвращает...        // что по этому поводу думаете? проверять сложно.. . подумайте, потом как-нибудь, если не пойдет, что можно вернуть чтобы убеждаться что из цепочки вышли и как        // я видел такие рекурсии поэтому абсолютно спокойно к этому отношусь    } // если обе окрестности еще ни разу не мельчились    // б) сетки уже адаптированные. Просто переходим на проверку на более точной сетке    // Е.Д.: не совсем уловила смысл этого, как такое может быть? здесь вопрос к АВ, можно в скайпе    // Е.Д.: может ли быть точка с уровнем ниже заданного? надо ли это проверить? что, если да? ...    if (grid[ii * coef * NY3_1 + jj * coef] != lev && grid[(ii + 1) * coef * NY3_1 + jj * coef] != lev) {        grid_hor_checking((3 * ii - 1),                          (3 * ii + 4),                          (3 * jj - 1),                          (3 * jj + 1),                          lev + 1,                          grid, rho, rhoPr, gridPr, tl);        grid_vert_checking((3 * ii - 1),                           (3 * ii + 4),                           (3 * jj - 1),                           (3 * jj + 1),                           lev + 1,                           grid, rho, rhoPr, gridPr, tl);        return;    } // сетки уже адаптированные. Просто переходим на проверку на более низком уровне    // в) одна из сеток не была адаптирована    if (grid[ii * coef * NY3_1 + jj * coef] == lev || grid[(ii + 1) * coef * NY3_1 + jj * coef] == lev) {        /* ( grid[ii * coef * NY3_1 + jj * coef] == lev && grid[(ii + 1) * coef * NY3_1 + jj * coef] != lev)         * ||         * ( grid[(ii + 1) * coef * NY3_1 + jj * coef] == lev && grid[ii * coef * NY3_1 + jj * coef] != lev)         *         * grid[ii * coef * NY3_1 + jj * coef] == lev || grid[(ii + 1) * coef * NY3_1 + jj * coef] == lev         */        int ii_ad = 0;        if (grid[ii * coef * NY3_1 + jj * coef] == lev) { // точно не может быть здесь и "то" и "то"?!            ii_ad = ii;        }        else {            ii_ad = ii + 1;        }        // теперь мельчим нужную окрестность        for (int p = -1; p < 2; ++p) { // p = -1, 0, 1            for (int q = -1; q < 2; ++q) { // q = -1, 0, 1                int i = (3 * ii_ad + p) * NY3_1 * coefNew + (3 * jj + q) * coefNew;                grid[i] = lev + 1;                double hx_lev = HX_SMALLEST * std::pow(3., R_LVL - (lev + 1));                double hy_lev = HX_SMALLEST * std::pow(3., R_LVL - (lev + 1));                double x = (3 * ii_ad + p) * coefNew * HX_SMALLEST;                double y = (3 * jj + q) * coefNew * HY_SMALLEST;                rho[i] = get_numeric_solution(x, y, hx_lev, hy_lev, rhoPr, gridPr, tl); // вычислим значение            }        }        // осуществим проверку для вновь созданных ячеек (окрестностей)        grid_hor_checking((3 * ii - 1),                          (3 * ii + 1),                          (3 * jj - 1),                          (3 * jj + 4),                          lev + 1,                          grid, rho, rhoPr, gridPr, tl);        grid_vert_checking((3 * ii - 1),                           (3 * ii + 1),                           (3 * jj - 1),                           (3 * jj + 4),                           lev + 1,                           grid, rho, rhoPr, gridPr, tl);        return;    } // одна из сеток не была адаптирована}static void grid_vert_refine(int ii, int jj, int lev, int *grid, double *rho, double *rhoPr, int *gridPr, int tl) {    assert(lev + 1 <= R_LVL);    int coef = (int) std::pow(3., R_LVL - lev);    int coefNew = (int) std::pow(3., R_LVL - (lev + 1));    // a) если обе окрестности еще ни разу не мельчились    if (grid[ii * NY3_1 * coef + jj * coef] == lev && grid[ii * NY3_1 * coef + (jj + 1) * coef] == lev) {        // изменим обе окрестности        for (int p = -1; p < 2; ++p) {            // p = -1, 0, 1            for (int q = -1; q < 5; ++q) {                // q = -1, 0, 1, 2, 3, 4                int i = (3 * ii + p) * NY3_1 * coefNew + (3 * jj + q) * coefNew;                grid[i] = lev + 1; // изменим сетку                double hx_lev = HX_SMALLEST * std::pow(3., R_LVL - (lev + 1));                double hy_lev = HX_SMALLEST * std::pow(3., R_LVL - (lev + 1));                double x = (3 * ii + p) * coefNew * HX_SMALLEST;                double y = (3 * jj + q) * coefNew * HY_SMALLEST;                rho[i] = get_numeric_solution(x, y, hx_lev, hy_lev, rhoPr, gridPr, tl); // вычислим значение            }        }        // осуществим проверку для вновь созданных ячеек (окрестностей)        grid_hor_checking((3 * ii - 1),                          (3 * ii + 1),                          (3 * jj - 1),                          (3 * jj + 4),                          lev + 1,                          grid, rho, rhoPr, gridPr, tl);        grid_vert_checking((3 * ii - 1),                           (3 * ii + 1),                           (3 * jj - 1),                           (3 * jj + 4),                           lev + 1,                           grid, rho, rhoPr, gridPr, tl);        return;    } // если обе окрестности еще ни разу не мельчились    // б) сетки уже адаптированные. Просто переходим на проверку на более низком уровне    if (grid[ii * coef * NY3_1 + jj * coef] != lev && grid[ii * coef * NY3_1 + (jj + 1) * coef] != lev) {        grid_hor_checking((3 * ii - 1),                          (3 * ii + 1),                          (3 * jj - 1),                          (3 * jj + 4),                          lev + 1,                          grid, rho, rhoPr, gridPr, tl);        grid_vert_checking((3 * ii - 1),                           (3 * ii + 1),                           (3 * jj - 1),                           (3 * jj + 4),                           lev + 1,                           grid, rho, rhoPr, gridPr, tl);        return;    } // сетки уже адаптированные. Просто переходим на проверку на более низком уровне    // в) одна из сеток не была адаптирована    if (grid[ii * coef * NY3_1 + jj * coef] == lev || grid[ii * coef * NY3_1 + (jj + 1) * coef] == lev) {        int jj_ad = 0;        if (grid[ii * coef * NY3_1 + jj * coef] == lev) {            jj_ad = jj;        }        else {            jj_ad = jj + 1;        }        // теперь мельчим нужную окрестность        for (int p = -1; p < 2; ++p) { // p = -1, 0, 1            for (int q = -1; q < 2; ++q) { // q = -1, 0, 1                int i = (ii * 3 + p) * NY3_1 * coefNew + (jj_ad * 3 + q) * coefNew;                grid[i] = lev + 1;                double hx_lev = HX_SMALLEST * std::pow(3., R_LVL - (lev + 1));                double hy_lev = HX_SMALLEST * std::pow(3., R_LVL - (lev + 1));                double x = (ii * 3 + p) * coefNew * HX_SMALLEST;                double y = (jj_ad * 3 + q) * coefNew * HY_SMALLEST;                rho[i] = get_numeric_solution(x, y, hx_lev, hy_lev, rhoPr, gridPr, tl); // вычислим значение            }        }        // осуществим проверку для вновь созданных ячеек (окрестностей)        grid_hor_checking((3 * ii - 1),                          (3 * ii + 1),                          (3 * jj - 1),                          (3 * jj + 4),                          lev + 1,                          grid, rho, rhoPr, gridPr, tl);        grid_vert_checking((3 * ii - 1),                           (3 * ii + 1),                           (3 * jj - 1),                           (3 * jj + 4),                           lev + 1,                           grid, rho, rhoPr, gridPr, tl);        return;    } // одна из сеток не была адаптирована} // конец grid_vert_refine//проверяем область [lx, rx] X [by, uy] включая границы этой области, а не нашего прямоугольника!static void grid_vert_checking(int lx, int rx, int by, int uy, int lev, int *grid, double *rho, double *rhoPr,                               int *gridPr, int tl) {    if (lev > R_LVL) {        printf("\ngrid_vert_checking: %d >lev_max ", lev);        return;    }    if (lev == R_LVL)        return;    int coef = (int) std::pow(3., R_LVL - lev);    for (int i = lx; i < rx + 1; ++i) {        for (int j = by; j < uy; ++j) {            if (std::abs(rho[i * NY3_1 * coef + j * coef] - rho[i * NY3_1 * coef + (j + 1) * coef]) > EPS_GRID) {                grid_vert_refine(i, j, lev, grid, rho, rhoPr, gridPr, tl);            }        }    }}//проверяем область [lx, rx] X [by, uy] включая границы этой области, а не нашего прямоугольника!static void grid_hor_checking(int lx, int rx, int by, int uy, int lev,                              int *grid, double *rho, double *rhoPr, int *gridPr, int tl) {    if (lev > R_LVL) {        printf("\ngrid_hor_checking: %d >lev_max ", lev);        return;    }    if (lev == R_LVL)        return;    int coef = (int) std::pow(3., R_LVL - lev);    for (int i = lx; i < rx; ++i) { // проверяем от lx до rx - 1 (диапозон [lx,...,rx]) я тут поправила опять щас подумаю        for (int j = by; j < uy + 1; ++j) { // проверяем от by до uy (диапозон [by,...,uy])            if (std::abs(rho[i * coef * NY3_1 + j * coef] - rho[(i + 1) * coef * NY3_1 + j * coef]) > EPS_GRID) {                grid_hor_refine(i, j, lev, grid, rho, rhoPr, gridPr, tl);            }        }    }}inline double get_l1(int *grid, double *rho, int nx3_1, int ny3_1, double hx_smallest, double hy_smallest,                     int max_lvl) {    double res = 0.;    for (int i = 0; i < nx3_1; ++i) {        for (int j = 0; j < ny3_1; ++j) {            int lev = grid[i * ny3_1 + j];            if (lev >= 0) {                double hx_lev = hx_smallest * std::pow(3., max_lvl - lev);                double hy_lev = hy_smallest * std::pow(3., max_lvl - lev);                double r = rho[i * ny3_1 + j] * hx_lev * hy_lev;                res += r;            }        }    }    return res;}double *solve_10(double &tme, int *grid, int *gridPr) {    StartTimer();    fflush(stdout);    double *rhoPr = new double[XY];    double *rho = new double[XY];    //<editor-fold desc="Init data">    for (int i = 0; i < XY; ++i) {        grid[i] = -1;        gridPr[i] = -1;        rho[i] = -1.;        rhoPr[i] = -1.;    }    for (int i = 0; i < NX_1; ++i) {        for (int j = 0; j < NY_1; ++j) {            grid[NY3_1 * i * R + j * R] = 0;            gridPr[NY3_1 * i * R + j * R] = 0;        }    }    for (int i = 0; i < NX_1; ++i) // G1 -- (x_i, 0=C) -- bottom boundary        rhoPr[NY3_1 * i * R] = analytical_solution_circle(0., A + HX * i, C);    for (int j = 1; j < NY; ++j) // G2 -- (NX=B, y_j) -- right boundary        rhoPr[NY3_1 * NX3 + j * R] = analytical_solution_circle(0., A + HX * NX, C + HY * j);    for (int i = 0; i < NX_1; ++i) // G3 -- (x_i, NY=D) -- top boundary        rhoPr[NY3_1 * i * R + NY3] = analytical_solution_circle(0., A + HX * i, C + HY * NY);    for (int j = 1; j < NY; ++j) // G4 -- (0=A, y_j) -- left boundary        rhoPr[j * R] = analytical_solution_circle(0., A, C + HY * j);    memcpy(rho, rhoPr, XY * sizeof(double));    for (int i = 1; i < NX; ++i) // inner points        for (int j = 1; j < NY; ++j)            rhoPr[NY3_1 * i * R + j * R] = analytical_solution_circle(0., A + HX * i, C + HY * j);    //</editor-fold>    printf("SUM RHO INIT = %le\n", calc_array_sum(rhoPr, NX_1, NY_1, 0));    printf("SUM ABS RHO INIT = %le\n", calc_array_sum(rhoPr, NX_1, NY_1, 1));    printf("L1 WIEGHTED INIT = %le\n", get_l1(grid, rhoPr, NX3_1, NY3_1, HX_SMALLEST, HY_SMALLEST, R_LVL));    print_surface_0("grid-0", NX3_1, NY3_1, NY3_1, HX_SMALLEST, HY_SMALLEST, 0, A, C, get_center_x(), get_center_y(), TAU, U, V, grid);//    print_surface_0("grid-pr-0", NX3_1, NY3_1, NY3_1, HX_SMALLEST, HY_SMALLEST, 0, A, C, get_center_x(), get_center_y(), TAU, U, V, grid);    print_surface_0("rho_init-all", NX3, NY3, NY3_1, HX_SMALLEST, HY_SMALLEST, 0, A, C, get_center_x(),                    get_center_y(), TAU, U, V, gridPr, rhoPr);    fflush(stdout);    for (int tl = 1; tl <= TIME_STEP_CNT; tl++) {        for (int i = 1; i < NX3; ++i) // G1 -- (x_i, 0=C) -- bottom boundary        {            if (G1[i] == 1) {                double val;                if (grid[NY3_1 * i] == -1) {                    val = -1.;                }                else {                    double hy_lev, hx_lev, hx_smallest, hy_smallest;                    get_hx_hy(i, 0, grid, NX, NY, R_LVL, hx_lev, hy_lev);                    get_hx_hy_smallest(NX, NY, R_LVL, hx_smallest, hy_smallest);                    double x1 = 0., y1 = 0., x2 = 0., y2 = 0., x3 = 0., y3 = 0., x4 = 0., y4 = 0.;                    get_bottom_area(i, 0, x1, y1, x2, y2, x3, y3, x4, y4, hx_smallest, hx_lev, hy_lev, A, B, C, D);                    val = calc_phi(x1, y1, x2, y2, x3, y3, x4, y4, rhoPr, TAU * tl, gridPr, NX, NY, TAU,                                   IDEAL_SQ_SIZE_X, IDEAL_SQ_SIZE_Y, R_LVL);                    double coef = 2. / (hx_lev * hy_lev);                    val = coef * val;                }                rho[NY3_1 * i] = val;            }        }        for (int j = 1; j < NY3; ++j) // G2 -- (NX=B, y_j) -- right boundary        {            if (G2[j] == 1) {                double val;                if (grid[NY3_1 * NX3 + j] == -1) {                    val = -1.;                }                else {                    double hy_lev, hx_lev, hx_smallest, hy_smallest;                    get_hx_hy(NX, j, grid, NX, NY, R_LVL, hx_lev, hy_lev);                    get_hx_hy_smallest(NX, NY, R_LVL, hx_smallest, hy_smallest);                    double x1 = 0., y1 = 0., x2 = 0., y2 = 0., x3 = 0., y3 = 0., x4 = 0., y4 = 0.;                    get_right_area(NX, 0, x1, y1, x2, y2, x3, y3, x4, y4, hy_smallest, hx_lev, hy_lev, A, B, C, D);                    val = calc_phi(x1, y1, x2, y2, x3, y3, x4, y4, rhoPr, TAU * tl, gridPr, NX, NY, TAU,                                   IDEAL_SQ_SIZE_X, IDEAL_SQ_SIZE_Y, R_LVL);                    double coef = 2. / (hx_lev * hy_lev);                    val = coef * val;                }                rho[NY3_1 * NX3 + j] = val;            }        }        for (int i = 1; i < NX3; ++i) // G3 -- (x_i, NY=D) -- top boundary        {            if (G3[i] == 1) {                double val;                if (grid[NY3_1 * i + NY3] == -1) {                    val = -1.;                }                else {                    double hy_lev, hx_lev, hx_smallest, hy_smallest;                    get_hx_hy(i, NY, grid, NX, NY, R_LVL, hx_lev, hy_lev);                    get_hx_hy_smallest(NX, NY, R_LVL, hx_smallest, hy_smallest);                    double x1 = 0., y1 = 0., x2 = 0., y2 = 0., x3 = 0., y3 = 0., x4 = 0., y4 = 0.;                    get_top_area(i, NY, x1, y1, x2, y2, x3, y3, x4, y4, hx_smallest, hx_lev, hy_lev, A, B, C, D);                    val = calc_phi(x1, y1, x2, y2, x3, y3, x4, y4, rhoPr, TAU * tl, gridPr, NX, NY, TAU,                                   IDEAL_SQ_SIZE_X, IDEAL_SQ_SIZE_Y, R_LVL);                    double coef = 2. / (hx_lev * hy_lev);                    val = coef * val;                }                rho[NY3_1 * i + NY3] = val;            }        }        for (int j = 1; j < NY3; ++j) // G4 -- (0=A, y_j) -- left boundary        {            if (G4[j] == 1) {                double val;                if (grid[j] == -1) {                    val = -1.;                }                else {                    double hy_lev, hx_lev, hx_smallest, hy_smallest;                    get_hx_hy(0, j, grid, NX, NY, R_LVL, hx_lev, hy_lev);                    get_hx_hy_smallest(NX, NY, R_LVL, hx_smallest, hy_smallest);                    double x1 = 0., y1 = 0., x2 = 0., y2 = 0., x3 = 0., y3 = 0., x4 = 0., y4 = 0.;                    get_left_area(0, j, x1, y1, x2, y2, x3, y3, x4, y4, hy_smallest, hx_lev, hy_lev, A, B, C, D);                    val = calc_phi(x1, y1, x2, y2, x3, y3, x4, y4, rhoPr, TAU * tl, gridPr, NX, NY, TAU,                                   IDEAL_SQ_SIZE_X, IDEAL_SQ_SIZE_Y, R_LVL);                    double coef = 2. / (hx_lev * hy_lev);                    val = coef * val;                }                rho[j] = val;            }        }        if (CP00 == 1) {            double val;            if (grid[0] == -1) {                val = -1.;            }            else {                double hy_lev, hx_lev, hx_smallest, hy_smallest;                get_hx_hy(0, 0, grid, NX, NY, R_LVL, hx_lev, hy_lev);                get_hx_hy_smallest(NX, NY, R_LVL, hx_smallest, hy_smallest);                double x1 = 0., y1 = 0., x2 = 0., y2 = 0., x3 = 0., y3 = 0., x4 = 0., y4 = 0.;                get_point_00_area(0, 0, x1, y1, x2, y2, x3, y3, x4, y4, hx_lev, hy_lev, A, B, C, D);                val = calc_phi(x1, y1, x2, y2, x3, y3, x4, y4, rhoPr, TAU * tl, gridPr, NX, NY, TAU,                               IDEAL_SQ_SIZE_X, IDEAL_SQ_SIZE_Y, R_LVL);                double coef = 4. / (hx_lev * hy_lev);                val = coef * val;            }            rho[0] = val; // point (0.0)        }        if (CP10 == 1) {            double val;            if (grid[NY3_1 * NX3] == -1) {                val = -1.;            }            else {                double hy_lev, hx_lev, hx_smallest, hy_smallest;                get_hx_hy(NX, 0, grid, NX, NY, R_LVL, hx_lev, hy_lev);                get_hx_hy_smallest(NX, NY, R_LVL, hx_smallest, hy_smallest);                double x1 = 0., y1 = 0., x2 = 0., y2 = 0., x3 = 0., y3 = 0., x4 = 0., y4 = 0.;                get_point_10_area(NX, 0, x1, y1, x2, y2, x3, y3, x4, y4, hx_lev, hy_lev, A, B, C, D);                val = calc_phi(x1, y1, x2, y2, x3, y3, x4, y4, rhoPr, TAU * tl, gridPr, NX, NY, TAU,                               IDEAL_SQ_SIZE_X, IDEAL_SQ_SIZE_Y, R_LVL);                double coef = 4. / (hx_lev * hy_lev);                val = coef * val;            }            rho[NY3_1 * NX3] = val; // point (1.0)        }        if (CP01 == 1) {            double val;            if (grid[NY3] == -1) {                val = -1.;            }            else {                double hy_lev, hx_lev, hx_smallest, hy_smallest;                get_hx_hy(0, NY, grid, NX, NY, R_LVL, hx_lev, hy_lev);                get_hx_hy_smallest(NX, NY, R_LVL, hx_smallest, hy_smallest);                double x1 = 0., y1 = 0., x2 = 0., y2 = 0., x3 = 0., y3 = 0., x4 = 0., y4 = 0.;                get_point_01_area(0, NY, x1, y1, x2, y2, x3, y3, x4, y4, hx_lev, hy_lev, A, B, C, D);                val = calc_phi(x1, y1, x2, y2, x3, y3, x4, y4, rhoPr, TAU * tl, gridPr, NX, NY, TAU,                               IDEAL_SQ_SIZE_X, IDEAL_SQ_SIZE_Y, R_LVL);                double coef = 4. / (hx_lev * hy_lev);                val = coef * val;            }            rho[NY3] = val; // point (0.1)        }        if (CP11 == 1) {            double val;            if (grid[NY3_1 * NX3 + NY3] == -1) {                val = -1.;            }            else {                double hy_lev, hx_lev, hx_smallest, hy_smallest;                get_hx_hy(NX, NY, grid, NX, NY, R_LVL, hx_lev, hy_lev);                get_hx_hy_smallest(NX, NY, R_LVL, hx_smallest, hy_smallest);                double x1 = 0., y1 = 0., x2 = 0., y2 = 0., x3 = 0., y3 = 0., x4 = 0., y4 = 0.;                get_point_11_area(NX, NY, x1, y1, x2, y2, x3, y3, x4, y4, hx_lev, hy_lev, A, B, C, D);                val = calc_phi(x1, y1, x2, y2, x3, y3, x4, y4, rhoPr, TAU * tl, gridPr, NX, NY, TAU,                               IDEAL_SQ_SIZE_X, IDEAL_SQ_SIZE_Y, R_LVL);                double coef = 4. / (hx_lev * hy_lev);                val = coef * val;            }            rho[NY3_1 * NX3 + NY3] = val; // point (1,1)        }        for (int i = 1; i < NX3; ++i) // inner points        {            for (int j = 1; j < NY3; ++j) {                double val;                if (grid[NY3_1 * i + j] == -1) {                    val = -1.;                }                else {                    double hy_lev, hx_lev, hx_smallest, hy_smallest;                    get_hx_hy(i, j, grid, NX, NY, R_LVL, hx_lev, hy_lev);                    get_hx_hy_smallest(NX, NY, R_LVL, hx_smallest, hy_smallest);                    double x1 = 0., y1 = 0., x2 = 0., y2 = 0., x3 = 0., y3 = 0., x4 = 0., y4 = 0.;                    get_inner_area(i, j, x1, y1, x2, y2, x3, y3, x4, y4, hx_smallest, hy_smallest, hx_lev, hy_lev, A, B,                                   C, D);                    val = calc_phi(x1, y1, x2, y2, x3, y3, x4, y4, rhoPr, TAU * tl, gridPr, NX, NY, TAU,                                   IDEAL_SQ_SIZE_X, IDEAL_SQ_SIZE_Y, R_LVL);                    double coef = 1. / (hx_lev * hy_lev);                    val = coef * val;                }                rho[NY3_1 * i + j] = val;            }        }        //print_surface_0("rho-adapt_minus-all", NX3, NY3, NY3_1, HX_SMALLEST, HY_SMALLEST, tl, A, C, get_center_x(),        //                get_center_y(), TAU, U, V, grid, rho);        //print_surface_0("rho-adapt-minus", NX, NY, NY3_1, R, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU, U, V, rho);        grid_hor_checking(1, NX - 1, 1, NY - 1, 0, grid, rho, rhoPr, gridPr, tl);        grid_vert_checking(1, NX - 1, 1, NY - 1, 0, grid, rho, rhoPr, gridPr, tl);        memcpy(rhoPr, rho, XY * sizeof(double));        memcpy(gridPr, grid, XY * sizeof(int));        printf("L1 WIEGHTED = %le\n", get_l1(grid, rhoPr, NX3_1, NY3_1, HX_SMALLEST, HY_SMALLEST, R_LVL));        fflush(stdout);        if (tl % 5 == 0) {            print_surface_0("rho", NX, NY, NY3_1, R, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU, U, V, rho);            print_surface_0("grid-all", NX3, NY3, NY3_1, HX_SMALLEST, HY_SMALLEST,                            tl, A, C, get_center_x(), get_center_y(), TAU, U, V, grid);            print_surface_0("rho-all", NX3, NY3, NY3_1, HX_SMALLEST, HY_SMALLEST, tl, A, C, get_center_x(),                            get_center_y(), TAU, U, V, grid, rho);        }    } // time loop    delete[] rhoPr;    tme = GetTimer() / 1000;    return rho;}double *calc_error_10(double hx, double hy, double tt, double *solution, int nx, int ny) {    double *res = new double[XY];    for (int i = 0; i < nx; i++)        for (int j = 0; j < ny; j++)            res[i * ny + j] = fabs(solution[i * ny + j]                                   - analytical_solution_circle(tt, A + hx * i, C + hy * j));    return res;}double *get_exact_solution_10(double hx, double hy, double t, int nx, int ny) {    double *res = new double[nx * ny];    for (int i = 0; i < nx; i++)        for (int j = 0; j < ny; j++)            res[i * ny + j] = fabs(analytical_solution_circle(t, A + hx * i, C + hy * j));    return res;}