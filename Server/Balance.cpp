#include "Balance.h"

void GetAngleFromAcc(double *acc, double *angle)
{
    double tilt_angle=atan((g-acc[1])/sqrt(acc[0]*acc[0]+acc[2]*acc[2]));
    double tilt_axis[3];

    double y_axis[3] {0,1,0};
    double acc_tmp[3] {0};
    acc_tmp[0]=acc[0];
    acc_tmp[2]=acc[2];
    aris::dynamic::s_cro3(y_axis,acc_tmp,tilt_axis);

    double Pmb[4][4];
    GeneralFunc::GetPmOfU(tilt_angle,tilt_axis,*Pmb);
    aris::dynamic::s_pm2pe(*Pmb,angle,"213");
}

void GetDeltaAnglePID(double *fce, double *delta_angle)
{

}
