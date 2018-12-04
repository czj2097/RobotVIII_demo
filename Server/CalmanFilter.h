#ifndef CALMAN_FILTER_H
#define CALMAN_FILTER_H

#include <cstddef>
#include <Eigen/LU>

template <size_t WIDTH>
class CalmanFilter
{
public:
    void SetCalmanParam(double *mtrx_F, double *input_u, double *noise_Q, double *mtrx_H, double *input_z, double *noise_R);
    void doFilter();

private:
    double predictMtrx_F[WIDTH][WIDTH];
    double predictInput_u[WIDTH];
    double predictNoise_Q[WIDTH][WIDTH];

    double predict_x[WIDTH];
    double predict_P[WIDTH][WIDTH];

    double measureMtrx_H[WIDTH][WIDTH];
    double measureInput_z[WIDTH];
    double measureNoise_R[WIDTH][WIDTH];

    double calmanGain_K[WIDTH][WIDTH];
    double calmanGain_tmp[WIDTH][WIDTH];
    double update_x[WIDTH];
    double update_P[WIDTH][WIDTH];

    double * inverse(const Eigen::MatrixBase<double>& mtrx);
};

template <size_t WIDTH>
double * CalmanFilter<WIDTH>::inverse(const Eigen::MatrixBase<double>& mtrx)
{
    return mtrx.inverse();
}

template <size_t WIDTH>
void CalmanFilter<WIDTH>::SetCalmanParam(double *mtrx_F, double *input_u, double *noise_Q, double *mtrx_H, double *input_z, double *noise_R)
{
    memcpy(*predictMtrx_F,mtrx_F,WIDTH*WIDTH*sizeof(double));
    memcpy(predictInput_u,input_u,WIDTH*sizeof(double));
    memcpy(*predictNoise_Q,noise_Q,WIDTH*WIDTH*sizeof(double));
    memcpy(*measureMtrx_H,mtrx_H,WIDTH*WIDTH*sizeof(double));
    memcpy(measureInput_z,input_z,WIDTH*sizeof(double));
    memcpy(*measureNoise_R,noise_R,WIDTH*WIDTH*sizeof(double));
}

template <size_t WIDTH>
void CalmanFilter<WIDTH>::doFilter(double *lst_x, double *lst_P)
{
    for(int i=0;i<WIDTH;i++)
    {
        predict_x[i]=predictInput_u[i];
        for(int j=0;j<WIDTH;j++)
        {
            predict_x[i]+=predictMtrx_F[i][j]*lst_x[j];

            predict_P[i][j]=predictNoise_Q[i][j];
            double sum {0};
            for(int m=0;m<WIDTH;m++)
            {
                for(int n=0;n<WIDTH;n++)
                {
                    sum+=predictMtrx_F[i][m]*lst_P[m][n]*predictMtrx_F[j][n];
                }
            }
            predict_P[i][j]+=sum;
        }
    }

    for(int i=0;i<WIDTH;i++)
    {
        for(int j=0;j<WIDTH;j++)
        {
            calmanGain_tmp[i][j]=measureNoise_R[i][j];
            double sum {0};
            for(int m=0;m<WIDTH;m++)
            {
                for(int n=0;n<WIDTH;n++)
                {
                    sum+=measureMtrx_H[i][m]*predict_P[m][n]*measureMtrx_H[j][n];
                }
            }
            calmanGain_tmp[i][j]+=sum;
        }
    }

}

#endif
