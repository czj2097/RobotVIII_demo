#ifndef REAL_TIME_OPTIMAL_H
#define REAL_TIME_OPTIMAL_H

#include <aris.h>
#include <Robot_Type_I.h>
#include <sys/time.h>
#include "Move_Gait.h"

struct P2PMotionParam
{
    VelType velType;
    int inopInterNum;//inoperative interval
    double minT;
    double inopInter[2];
    double trajAcc[4];//trapezoidal trajectory
    double trajTim[4];
    double trajPos[4];
};
struct PolylineParam
{
    int fstBrkPntNum;
    int lstBrkPntNum;
    double fstBrkPntVel;
    double lstBrkPntVel;

    double pntsI[2];
    double pntsII[3];
    double pntsIII[4];

    double nxtPntTime;
    double remainTime;
};
enum VelType
{
    Acc,
    Dec,
    AccDec,
    DecAcc,
    AccPosDec,
    DecNegAcc,
};
enum LineType //polyline type
{
    I,
    II,
    III,
    IV,
};

class RealTimeOptimal
{
public:
    RealTimeOptimal();
    ~RealTimeOptimal();
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


    LineType lineType;
    P2PMotionParam fstP2PParam;
    P2PMotionParam lstP2PParam;
    //P2PMotionParam curP2PParam;
    PolylineParam lineParam;

    double curPos;
    double nxtPos;

    int segStartNum;
    int segEndNum;
    double segStartVel;
    double segEndVel;
    bool isCurSegFnshed;

    double curPnt[3];
    double curSeg[3];
    double curPntVel[3];
    double curSegVel[3];
    double curMinTime;
    P2PMotionParam curP2PParam[3];

    void Init();
    void JudgeLineType(double *pnts, int pntsNum);
    void GetInitSegParam(double *pnts, int pntsNum, double startV, double endV);
    void GetTwoTimeThreeDof();
    void GetCurSegOptParam(P2PMotionParam *p2pParam);

    void GetParamOfTypeI(double *pnts, int pntsNum, double startV, double endV);
    void GetCurSegParam(double segStartP, double segEndP, double segStartV, double segEndV, P2PMotionParam &p2pParam);
    void GetParamOfTypeIII(double *pnts, int pntsNum, double startV, double endV);

    void InitP2PParam(P2PMotionParam &p2pParam);
    void GetP2PMotionParam(double startP, double endP, double startV, double endV, double midV, P2PMotionParam &p2pParam);
    void GetOptimalP2PMotionAcc(double startP, double endP, double startV, double endV, P2PMotionParam &p2pParam);
    void GetOptimalP2PMotionJerk(double startP, double endP, double startV, double endV);

    void tg(double *pnts, int pntsNum, double startV, double endV, int count);
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
