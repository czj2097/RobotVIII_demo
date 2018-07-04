#ifndef REAL_TIME_OPTIMAL_H
#define REAL_TIME_OPTIMAL_H

#include <aris.h>
#include <Robot_Type_I.h>
#include <sys/time.h>

void dlmwrite(const char *filename, const double *mtx, const int m, const int n);

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
    double startPos;
    double endPos;
    double startVel;
    double endVel;
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
    int fstBrkPntNum;
    int lstBrkPntNum;
    double fstBrkPntVel;
    double lstBrkPntVel;

    int passNxtPntNum;
    double passNxtPntTime[4];//most 3 times, use the fourth ele to deal with error impossible to happen
    double totalTime;
    double nxtPntTime;
};
struct NextPointParam
{
    double min_t1;
    double min_t2;
    double min_vm;
    double max_t1;
    double max_t2;
    double max_vm;
    double rch_t1;
    double rch_t2;
    double rch_vm;
    double reachableVel[2];
    double optVel;

    double record_t1;
    double record_t2;
};



class RTOptimal
{
public:
    RTOptimal();
    ~RTOptimal();
    void GetParamInFromParamEE(double *pnts, int pntsNum, double *startV, double *endV, int legID);
    bool GetTrajOneLeg(double *pnts, int pntsNum, double *startV, double *endV, int legID, int count, double *pIn, double *pEE);


    double nxtPntTime[3000][9];
    void RecordNxtPntTime(NextPointParam &nxtParam, int count, int screwID);

private:
    const double stepD {0.5};
    const double stepH {0.05};
    const double vLmt {1.0};
    const double aLmt {3.2};
    const double initPeb[6] {0};
    const double initVeb[6] {0};
    const double initPee[18]{ -0.3, -0.85, -0.65,
                       -0.45, -0.85, 0,
                        -0.3, -0.85, 0.65,
                         0.3, -0.85, -0.65,
                        0.45, -0.85, 0,
                         0.3, -0.85, 0.65 };//Robot_VIII

//    const double initPee[18] { -0.30, -0.58, -0.52,
//                         -0.60, -0.58,  0,
//                         -0.30, -0.58,  0.52,
//                          0.30, -0.58, -0.52,
//                          0.60, -0.58,  0,
//                          0.30, -0.58,  0.52 };//Robot EDU2
    const double initVee[18] {0};

    Robots::RobotTypeI rbt;
    PolylineParam lineParam[18];
    P2PMotionParam curP2PParam[18];
    NextPointParam nxtPntParam[18];

    double nxtPntMinTime;
    int actScrewID;
    bool isCurPntPassed;
    double lstPntTime;

    //void GetParamInFromParamEE(double *pnts, int pntsNum, double *startV, double *endV, int legID);
    void JudgeLineType(PolylineParam &lnParam);

    void GetTwoTimeTypeI(PolylineParam &lnParam, P2PMotionParam &p2pParam);
    void GetTwoTimeTypeII(PolylineParam &lnParam, P2PMotionParam &p2pParam);
    void GetTwoTimeTypeIII(PolylineParam &lnParam, P2PMotionParam &p2pParam);
    double GetBrkPntSndPassVel(double startP, double endP, double startV, double endV, double pnt);
    void GetThreeNxtPntTime(PolylineParam &lnParam, P2PMotionParam &p2pParam);
    double LocatePntInLmtInterval(P2PMotionParam &p2pParam, double pnt, double *limitedTimeInterval);
    void DecideNxtPntTime(PolylineParam *lnParam);

    void GetNxtPntOptVel(PolylineParam &lnParam, NextPointParam &nxtParam);
    bool GetNxtPntReachableVel(PolylineParam &lnParam, NextPointParam &nxtParam, int screwID);
    double GetNxtPin(PolylineParam &lnParam, NextPointParam &nxtParam, int count);

    void GetP2PMotionParam(P2PMotionParam &p2pParam, double midV);
    void GetOptimalP2PMotionAcc(double startP, double endP, double startV, double endV, P2PMotionParam &p2pParam);
    void GetOptimalP2PMotionJerk(double startP, double endP, double startV, double endV);

    //void GetTrajOneLeg(double *pnts, int pntsNum, double *startV, double *endV, int legID, double *pIn);
};

#endif
