#ifndef BALANCE_H
#define BALANCE_H

#define g 9.80665

#include "GeneralFunc.h"
#include <aris.h>

void GetAngleFromAcc(double *acc, double *angle);
void GetDeltaAnglePID(double *fce, double *delta_angle);

#endif
