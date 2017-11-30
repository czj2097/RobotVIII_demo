#ifndef REAL_TIME_OPTIMAL_H
#define REAL_TIME_OPTIMAL_H

#include <aris.h>
#include <Robot_Type_I.h>
#include <sys/time.h>
#include "Move_Gait.h"

class RealTimeOptimal
{
public:
    RealTimeOptimal();
    ~RealTimeOptimal();
    void screwInterpolationTraj();

private:
    const double vLmt {1.0};
    const double aLmt {3.2};
    const double stepH {0.05};
    const double stepD {0.5};
    const double initPeb[6] {0};
    const double initVeb[6] {0};
    //const double initPee[18]{ -0.3, -0.85, -0.65,
    //                   -0.45, -0.85, 0,
    //                    -0.3, -0.85, 0.65,
    //                     0.3, -0.85, -0.65,
    //                    0.45, -0.85, 0,
    //                     0.3, -0.85, 0.65 };

    const double initPee[18] { -0.30, -0.58, -0.52,
                         -0.60, -0.58,  0,
                         -0.30, -0.58,  0.52,
                          0.30, -0.58, -0.52,
                          0.60, -0.58,  0,
                          0.30, -0.58,  0.52 };
    const double initVee[18] {0};
    double vIn[3000][18];
    double pIn[3000][18];
    void GetSpline(double startP, double endP, double startV, double endV, double *a, double *t);
    void GetTraj(int screwID, int startCount, int totalCount, double startP, double endP, double startV, double endV, double *a, double *t);
    void ScalingTraj(double startP, double endP, double startV, double endV, double *a, double *t);

    void GetBezier(int screwID, double startP, double endP, double startV, double endV);
};

struct JointSpaceWalkParam final:public aris::server::GaitParamBase
{
};

class JointSpaceWalk
{
public:
    JointSpaceWalk();
    ~JointSpaceWalk();
    static void parseJointSpaceFastWalk(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg);
    static int jointSpaceFastWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

private:
    static void swingLegTg(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in, int legID);
    static void stanceLegTg(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in, int legID);

    static double bodyAcc;
    static double bodyDec;
    static int totalCount;
    static double height;
    static double beta;
    static NormalGait::WalkState walkState;

    static double beginPee[18];
    static double beginVel;
    static double endPee[18];
    static double endVel;
    static double distance;

    static NormalGait::GaitPhase gaitPhase[6];//swing true, stance false
    static bool constFlag;
};


#endif
