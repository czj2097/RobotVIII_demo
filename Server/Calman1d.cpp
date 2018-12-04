#include "Calman1d.h"

void Calman1d::init(double mtrx_F, double noise_Q, double mtrx_H, double noise_R, double input_u, double input_z)
{
    predictMtrx_F=mtrx_F;
    predictNoise_Q=noise_Q;
    measureMtrx_H=mtrx_H;
    measureNoise_R=noise_R;
    predictInput_u=input_u;
    measureInput_z=input_z;
}

void Calman1d::init(double mtrx_F, double noise_Q, double mtrx_H, double noise_R)
{
    predictMtrx_F=mtrx_F;
    predictNoise_Q=noise_Q;
    measureMtrx_H=mtrx_H;
    measureNoise_R=noise_R;
}

void Calman1d::SetCalmanParam(double input_u, double input_z)
{
    predictInput_u=input_u;
    measureInput_z=input_z;
}

void Calman1d::doFilter(double &lst_x, double &pre_x, double &lst_P)
{
    predict_x=predictInput_u+predictMtrx_F*lst_x;
    predict_P=predictNoise_Q+predictMtrx_F*lst_P*predictMtrx_F;

    double gain_tmp {0};
    gain_tmp=measureNoise_R+measureMtrx_H*predict_P*measureMtrx_H;
    calmanGain_K=predict_P*measureMtrx_H/gain_tmp;

    double z_tmp {0};
    z_tmp=measureInput_z-measureMtrx_H*predict_x;
    update_x=predict_x+calmanGain_K*z_tmp;
    update_P=predict_P-calmanGain_K*measureMtrx_H*predict_P;

    lst_x=update_x;
    lst_P=update_P;
    pre_x=predict_x;
}
