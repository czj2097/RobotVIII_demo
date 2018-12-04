#include <cstring>
#include "Calman2d.h"
#include "GeneralFunc.h"

void Calman2d::init(double *mtrx_F, double *noise_Q, double *mtrx_H, double *noise_R, double *input_u, double *input_z)
{
    memcpy(*predictMtrx_F,mtrx_F,4*sizeof(double));
    memcpy(*predictNoise_Q,noise_Q,4*sizeof(double));
    memcpy(*measureMtrx_H,mtrx_H,4*sizeof(double));
    memcpy(*measureNoise_R,noise_R,4*sizeof(double));
    memcpy(predictInput_u,input_u,2*sizeof(double));
    memcpy(measureInput_z,input_z,2*sizeof(double));
}

void Calman2d::init(double *mtrx_F, double *noise_Q, double *mtrx_H, double *noise_R)
{
    memcpy(*predictMtrx_F,mtrx_F,4*sizeof(double));
    memcpy(*predictNoise_Q,noise_Q,4*sizeof(double));
    memcpy(*measureMtrx_H,mtrx_H,4*sizeof(double));
    memcpy(*measureNoise_R,noise_R,4*sizeof(double));
}

void Calman2d::SetCalmanParam(double *input_u, double *input_z)
{
    memcpy(predictInput_u,input_u,2*sizeof(double));
    memcpy(measureInput_z,input_z,2*sizeof(double));
}

void Calman2d::doFilter(double *lst_x, double *pre_x, double *lst_P)
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
                    sum+=predictMtrx_F[i][m]*lst_P[m][n]*predictMtrx_F[j][n];
                }
            }
            predict_P[i][j]+=sum;
        }
    }

    double gain_tmp[2][2] {0};
    double inv_tmp[2][2] {0};
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<2;j++)
        {
            gain_tmp[i][j]=measureNoise_R[i][j];
            double sum {0};
            for(int m=0;m<2;m++)
            {
                for(int n=0;n<2;n++)
                {
                    sum+=measureMtrx_H[i][m]*predict_P[m][n]*measureMtrx_H[j][n];
                }
            }
            gain_tmp[i][j]+=sum;
        }
    }
    GeneralFunc::inv2(*gain_tmp,*inv_tmp);

    for(int i=0;i<2;i++)
    {
        for(int j=0;j<2;j++)
        {
            calmanGain_K[i][j]=0;
            double sum {0};
            for(int m=0;m<2;m++)
            {
                for(int n=0;n<2;n++)
                {
                    sum+=predict_P[i][m]*measureMtrx_H[n][m]*inv_tmp[n][j];
                }
            }
            calmanGain_K[i][j]+=sum;
        }
    }

    double z_tmp[2];
    for(int i=0;i<2;i++)
    {
        z_tmp[i]=measureInput_z[i]-(measureMtrx_H[i][0]*predict_x[0]+measureMtrx_H[i][1]*predict_x[1]);
    }
    for(int i=0;i<2;i++)
    {
        update_x[i]=predict_x[i]+calmanGain_K[i][0]*z_tmp[0]+calmanGain_K[i][1]*z_tmp[1];
    }

    for(int i=0;i<2;i++)
    {
        for(int j=0;j<2;j++)
        {
            update_P[i][j]=predict_P[i][j];
            double sum {0};
            for(int m=0;m<2;m++)
            {
                for(int n=0;n<2;n++)
                {
                    sum-=calmanGain_K[i][m]*measureMtrx_H[m][n]*predict_P[n][j];
                }
            }
            update_P[i][j]+=sum;
        }
    }

    memcpy(lst_x,update_x,2*sizeof(double));
    memcpy(lst_P,*update_P,4*sizeof(double));
    memcpy(pre_x,predict_x,2*sizeof(double));
}
