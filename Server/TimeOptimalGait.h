#ifndef TIMEOPTIMALGAIT_H
#define TIMEOPTIMALGAIT_H

class TimeOptimalGait
{
public:
    TimeOptimalGait();
    ~TimeOptimalGait();

    void SetInitPos(double *pEB, double *pEE);
    void SetBodyTraj(double *traj);
    void SetSwingLegTraj(int legID, double *traj);

private:
    const int totalCount {1800};
    const double vLmt {0.9};
    const double aLmt {3.2};

    double initPeb[6];
    double initPee[18];
    double bodyTraj[totalCount+1];
    double swingLegTraj[totalCount+1][6];


};

#endif
