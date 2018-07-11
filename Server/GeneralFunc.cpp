#include "GeneralFunc.h"
#include <cmath>

namespace GeneralFunc
{
void GetPmOfU(double angle, double *axis, double *pm)
{
    double s=sin(angle);
    double c=cos(angle);
    double v=1-c;

    double u_len=sqrt(axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2]);
    double ux,uy,uz;
    if(u_len!=1)
    {
        ux=axis[0]/u_len;
        uy=axis[1]/u_len;
        uz=axis[2]/u_len;
    }

    pm[0]=ux*ux*v+c;
    pm[1]=ux*uy*v-uz*s;
    pm[2]=ux*uz*v+uy*s;
    pm[3]=0;

    pm[4]=ux*uy*v+uz*s;
    pm[5]=uy*uy*v+c;
    pm[6]=uy*uz*v-ux*s;
    pm[7]=0;

    pm[8]=ux*uz*v-uy*s;
    pm[9]=uy*uz*v+ux*s;
    pm[10]=uz*uz*v+c;
    pm[11]=0;

    pm[12]=0;
    pm[13]=0;
    pm[14]=0;
    pm[15]=1;
}

}
