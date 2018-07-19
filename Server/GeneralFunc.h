#ifndef GENERAL_FUNC_H
#define GENERAL_FUNC_H

namespace GeneralFunc
{
    void GetPmOfU(double angle, double *axis, double *pm);
    void Get5thPolynomial(double startP, double startV, double startA, double endP, double endV, double endA, double endT, double *c);
    void inv3(double * mtrx,double * inv_mtrx);
    double norm(double *vec_in);

}

namespace Controller
{
    double ApplyPID(double err, double kp, double ki, double kd, double delta_t);
    double SndOrderLag(double startP, double input, double wn, double damping, double delta_t);
}

#endif
