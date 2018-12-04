#ifndef CALMAN_2D_H
#define CALMAN_2D_H

class Calman2d
{
public:
    void init(double *mtrx_F, double *noise_Q, double *mtrx_H, double *noise_R, double *input_u, double *input_z);
    void init(double *mtrx_F, double *noise_Q, double *mtrx_H, double *noise_R);
    void SetCalmanParam(double *input_u, double *input_z);
    void doFilter(double *lst_x, double *pre_x, double *lst_P);

private:
    double predictMtrx_F[2][2];
    double predictInput_u[2];
    double predictNoise_Q[2][2];

    double predict_x[2];
    double predict_P[2][2];

    double measureMtrx_H[2][2];
    double measureInput_z[2];
    double measureNoise_R[2][2];

    double calmanGain_K[2][2];
    double update_x[2];
    double update_P[2][2];
};

#endif
