#include <cstring>
#include "Calman2dot1.h"
#include "GeneralFunc.h"

void Calman2dot1::init(double *mtrx_F, double *noise_Q, double *mtrx_H, double noise_R, double *input_u, double input_z)
{
    memcpy(*predictMtrx_F,mtrx_F,4*sizeof(double));
    memcpy(*predictNoise_Q,noise_Q,4*sizeof(double));
    memcpy(*measureMtrx_H,mtrx_H,4*sizeof(double));
    measureNoise_R=noise_R;
    memcpy(predictInput_u,input_u,2*sizeof(double));
    measureInput_z=input_z;
}

void Calman2dot1::init(double *mtrx_F, double *noise_Q, double *mtrx_H, double noise_R)
{
    memcpy(*predictMtrx_F,mtrx_F,4*sizeof(double));
    memcpy(*predictNoise_Q,noise_Q,4*sizeof(double));
    memcpy(*measureMtrx_H,mtrx_H,4*sizeof(double));
    measureNoise_R=noise_R;
}

void Calman2dot1::SetCalmanParam(double *input_u, double input_z)
{
    memcpy(predictInput_u,input_u,2*sizeof(double));
    measureInput_z=input_z;
}

void Calman2dot1::doFilter(double *lst_x, double *lst_P)
{
    for(int i=0;i<2;i++)
    {
        predict_x[i]=predictInput_u[i];
        for(int j=0;j<2;j++)
        {
            predict_x[i]+=predictMtrx_F[i][j]*lst_x[j];

            predict_P[i][j]=predictNoise_Q[i][j];
            double sum {0};
            for(int m=0;m<2;m++)
            {
                for(int n=0;n<2;n++)
                {
                    sum+=predictMtrx_F[i][m]*lst_P[2*m+n]*predictMtrx_F[j][n];
                }
            }
            predict_P[i][j]+=sum;
        }
    }

    double gain_tmp {0};
    double inv_tmp {0};
    gain_tmp=measureNoise_R;
    double sum {0};
    for(int m=0;m<2;m++)
    {
        for(int n=0;n<2;n++)
        {
            sum+=measureMtrx_H[m]*predict_P[m][n]*measureMtrx_H[n];
        }
    }
    gain_tmp+=sum;
    inv_tmp=1.0/gain_tmp;

    for(int i=0;i<2;i++)
    {
        calmanGain_K[i]=0;
        double sum {0};
        for(int m=0;m<2;m++)
        {
            for(int n=0;n<2;n++)
            {
                sum+=predict_P[i][m]*measureMtrx_H[m]*inv_tmp;
            }
        }
        calmanGain_K[i]+=sum;
    }

    double z_tmp;
    z_tmp=measureInput_z-(measureMtrx_H[0]*predict_x[0]+measureMtrx_H[1]*predict_x[1]);

    for(int i=0;i<2;i++)
    {
        update_x[i]=predict_x[i]+calmanGain_K[i]*z_tmp;
    }

    for(int i=0;i<2;i++)
    {
        for(int j=0;j<2;j++)
        {
            update_P[i][j]=predict_P[i][j];
            double sum {0};
            for(int n=0;n<2;n++)
            {
                sum-=calmanGain_K[i]*measureMtrx_H[n]*predict_P[n][j];
            }
            update_P[i][j]+=sum;
        }
    }

    memcpy(lst_x,update_x,2*sizeof(double));
    memcpy(lst_P,*update_P,4*sizeof(double));
}
