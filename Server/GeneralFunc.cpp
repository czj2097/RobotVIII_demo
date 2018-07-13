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

    //Bezier form(Cited from Chen Xianbao): c0*(1-t)^5+c1*t(1-t)^4+c2*t^2*(1-t)^3+c3*t^3*(1-t)^2+c4*t^4*(1-t)+c5*t^5, t=(curT-startT)/totalT
    void Get5thPolynomial(double startP, double startV, double startA, double endP, double endV, double endA, double endT, double *c)
    {
        c[0]=startP;
        c[1]=(endT*startV+5*startP)/5;
        c[2]=(20*startP+endT*endT*startA+8*endT*endV)/20;
        c[3]=(20*endP+endT*endT*endA-8*endT*endV)/20;
        c[4]=(-endT*endV+5*endP)/5;
        c[5]=endP;
    }

    void inv3(double * mtrx,double * inv_mtrx)
    {
        double a1=mtrx[0];
        double b1=mtrx[1];
        double c1=mtrx[2];
        double a2=mtrx[3];
        double b2=mtrx[4];
        double c2=mtrx[5];
        double a3=mtrx[6];
        double b3=mtrx[7];
        double c3=mtrx[8];

        double value_mtrx=a1*(b2*c3-c2*b3)-a2*(b1*c3-c1*b3)+a3*(b1*c2-c1*b2);

        inv_mtrx[0]=(b2*c3-c2*b3)/value_mtrx;
        inv_mtrx[1]=(c1*b3-b1*c3)/value_mtrx;
        inv_mtrx[2]=(b1*c2-c1*b2)/value_mtrx;
        inv_mtrx[3]=(c2*a3-a2*c3)/value_mtrx;
        inv_mtrx[4]=(a1*c3-c1*a3)/value_mtrx;
        inv_mtrx[5]=(a2*c1-a1*c2)/value_mtrx;
        inv_mtrx[6]=(a2*b3-a3*b2)/value_mtrx;
        inv_mtrx[7]=(b1*a3-a1*b3)/value_mtrx;
        inv_mtrx[8]=(a1*b2-b1*a2)/value_mtrx;
    }

    double norm(double *vec_in)
    {
        return	sqrt(vec_in[0]*vec_in[0]+vec_in[1]*vec_in[1]+vec_in[2]*vec_in[2]);
    }
}
