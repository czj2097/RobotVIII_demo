#ifndef REAL_TIME_OPTIMAL_H
#define REAL_TIME_OPTIMAL_H

#include <aris.h>
#include <Robot_Type_I.h>
#include <sys/time.h>
#include "Move_Gait.h"

void screwInterpolationTraj();
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
