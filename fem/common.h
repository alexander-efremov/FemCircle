#ifndef FEM_CIRCLE_COMMON_H
#define FEM_CIRCLE_COMMON_H

#include "consts.h"

inline double get_center_x()
{
    return A + CENTER_OFFSET_X;
}

inline double get_center_y()
{
    return C + CENTER_OFFSET_Y;
}

#endif //FEM_CIRCLE_COMMON_H
