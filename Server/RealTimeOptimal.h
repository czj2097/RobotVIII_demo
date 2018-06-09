#ifndef REAL_TIME_OPTIMAL_H
#define REAL_TIME_OPTIMAL_H

#include <aris.h>
#include <Robot_Type_I.h>
#include <sys/time.h>
#include "Move_Gait.h"

enum VelType
{
    AccDec,
    DecAcc,
    AccConstDec,
    DecConstAcc,
};
struct P2PMotionParam
{
    VelType velType;
    int inopInterNum;//inoperative interval
    double minTime;
    double inopInter[2];
    double trajAcc[3];//trapezoidal trajectory
    double trajTime[3];
    double trajPos[3];
    double maxVel;
};
enum PolylineType //polyline type
{
    I,
    II,
    III,
};
struct PolylineParam
{
    PolylineType lineType;
    double pnts[1000];
    int pntsNum;
    double startV;
    double endV;
    double nxtPntTime;
    int fstBrkPntNum;
    int lstBrkPntNum;
    double fstBrkPntVel;
    double lstBrkPntVel;
};



class RTOptimal
{
public:
    RTOptimal();
    ~RTOptimal();
    void screwInterpolationTraj();

private:
    const double vLmt {1.0};
    const double aLmt {3.2};
    const double initPeb[6] {0};
    const double initVeb[6] {0};
    //const double initPee[18]{ -0.3, -0.85, -0.65,
    //                   -0.45, -0.85, 0,
    //                    -0.3, -0.85, 0.65,
    //                     0.3, -0.85, -0.65,
    //                    0.45, -0.85, 0,
    //                     0.3, -0.85, 0.65 };//Robot_VIII

    const double initPee[18] { -0.30, -0.58, -0.52,
                         -0.60, -0.58,  0,
                         -0.30, -0.58,  0.52,
                          0.30, -0.58, -0.52,
                          0.60, -0.58,  0,
                          0.30, -0.58,  0.52 };//Robot EDU2
    const double initVee[18] {0};
    double vIn[3000][18];
    double pIn[3000][18];
    void GetTraj(int screwID, int startCount, int totalCount, double startP, double endP, double startV, double endV, double *a, double *t);
    void ScalingTraj(double startP, double endP, double startV, double endV, double *a, double *t);

    Robots::RobotTypeI rbt;
    PolylineParam lineParam[18];
    P2PMotionParam curP2PParam[18];

    double curPnt[3];
    double curSeg[3];
    double curPntVel[3];
    double curSegVel[3];
    double curMinTime;

    void GetParamInFromParamEE(double *pnts, int pntsNum, double *startV, double *endV, int legID);
    void JudgeLineType(PolylineParam &lnParam);
    void GetNxtPntTime(PolylineParam &lnParam, P2PMotionParam &p2pParam);
    void GetNxtPntTimeTypeI(PolylineParam &lnParam, P2PMotionParam &p2pParam);
    void GetNxtPntTimeTypeII(PolylineParam &lnParam, P2PMotionParam &p2pParam);
    void GetNxtPntTimeTypeIII(PolylineParam &lnParam, P2PMotionParam &p2pParam);

    void GetTwoTimeThreeDof();
    void GetCurSegOptParam(P2PMotionParam *p2pParam);

    void GetP2PMotionParam(double startP, double endP, double startV, double endV, double midV, P2PMotionParam &p2pParam);
    void GetOptimalP2PMotionAcc(double startP, double endP, double startV, double endV, P2PMotionParam &p2pParam);
    void GetOptimalP2PMotionJerk(double startP, double endP, double startV, double endV);

    void GetTrajOneLeg(double *pnts, int pntsNum, double *startV, double *endV, int legID);
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
