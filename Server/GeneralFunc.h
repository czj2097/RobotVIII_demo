#ifndef GENERAL_FUNC_H
#define GENERAL_FUNC_H

namespace GeneralFunc
{
    void GetPmOfU(double angle, double *axis, double *pm);
    void Get5thPolynomial(double startP, double startV, double startA, double endP, double endV, double endA, double endT, double *c);
    void inv3(double * mtrx,double * inv_mtrx);
    double norm(double *vec_in);
    void forceInit(int count, const double* forceRaw_in, double* forceInF_out);
    void FitCycle2D(double *pnts, int pntsNum, double *cycle);
}

namespace Controller
{
    struct lstPIDparam
    {
        double lstInt;
        double lstErr;
    };
    double doPID(lstPIDparam &param, double err, double kp, double ki, double kd, double delta_t);

    struct lstLagParam
    {
        double lstFstInt;
        double lstSndInt;
    };
    double SndOrderLag(lstLagParam &param,double startP, double input, double wn, double damping, double delta_t);
}

#endif
