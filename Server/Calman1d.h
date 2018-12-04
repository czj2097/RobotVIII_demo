#ifndef CALMAN_1D_H
#define CALMAN_1D_H

class Calman1d
{
public:
    void init(double mtrx_F, double noise_Q, double mtrx_H, double noise_R, double input_u, double input_z);
    void init(double mtrx_F, double noise_Q, double mtrx_H, double noise_R);
    void SetCalmanParam(double input_u, double input_z);
    void doFilter(double &lst_x, double &pre_x, double &lst_P);

private:
    double predictMtrx_F;
    double predictInput_u;
    double predictNoise_Q;

    double predict_x;
    double predict_P;

    double measureMtrx_H;
    double measureInput_z;
    double measureNoise_R;

    double calmanGain_K;
    double update_x;
    double update_P;
};

#endif
