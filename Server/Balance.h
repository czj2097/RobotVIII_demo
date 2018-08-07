#ifndef BALANCE_H
#define BALANCE_H

#include <aris.h>
#include <aris_control_pipe.h>
#include <Robot_Type_I.h>
#include "GeneralFunc.h"
#include "NormalGait.h"
#include "EmergencyStop.h"

using namespace aris::control;

void GetTargetEulFromAcc(double *planAccInG, double *reqAccInG, double *targetEul);
double GetAngleFromAcc(double acc);

struct BalanceParam final :public aris::server::GaitParamBase
{
    int preCount;
    int totalCount;
    int n;
    double d;
};

enum BalanceWalkState
{
    Init,
    Balance,
    ForwardAcc,
    TurnAcc,
    Const,
    Stop,
};
enum GaitPhase
{
    Swing,
    Stance,
    Touch,
};

struct BallBalanceParam
{
    double fceInB[6];
    double fceInB_filtered[6];
    double bodyPos[3];
    double ballPos[2];
    double tarAcc[2];
    double tarEul[2];
    double actEul[3];
    double imuEul[3];
};

class BallBalance
{
public:
    static void parseBallBalance(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
    static int ballBalance(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

    static void bodyPosTg(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);
    static void bodyEulTg(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);
    static void recordData();

private:
    static double realPeb213[6];
    static BalanceWalkState workState;
    static int countIter;

    static BallBalanceParam bbParam;
    static Pipe<BallBalanceParam> bbPipe;
    static std::thread bbThread;
};

struct BalanceWalkParam:public BallBalanceParam
{
    double legForce[6];
    double footPos[18];
};
class BalanceWalk
{
public:
    static void parseBalanceWalk(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
    static int balanceWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);
    static void recordData();

private:
    static double forwardAcc;
    static double turnAcc;
    static int totalCount;
    static int totalCount_tmp;
    static double height;
    static double height_tmp;
    static double alpha;
    static double alpha_tmp;

    static BalanceWalkState walkState;
    static bool constFlag;
    static double beginXVel;
    static double endXVel;
    static double beginZVel;
    static double endZVel;
    static double beginOmega;
    static double endOmega;

    static double beginPeb[6];
    static double pEB[6];
    static GaitPhase gaitPhase[6];
    static double swingPee[18];
    static double swingBeginPee[18];
    static double swingEndPee[18];
    static double stancePee[18];
    static double stanceBeginPee[18];
    static double stanceEndPee[18];
    static double followPee[18];
    static double followRecordPee[18];
    static double followBeginPee[18];
    static double followBeginVee[18];
    static int followBeginCount[6];
    static bool followFlag[6];
    static bool filterFlag[6];
    static int filterCount[6];

    static double initPee[18];
    static double avgRealH;
    static double planH;

    static double bodyEul213[3];
    static double sumEul[3];
    static double avgEul[3];
    static double targetEul[3];
    static double inputEul[3];
    static double inputEul_tmp[3];

    static double forceRaw[42];
    static double forceInF[42];
    static double forceSum[42];
    static double forceAvg[42];

    static void forceInit(int count, int legID);
    static void swingLegTg(const aris::dynamic::PlanParamBase &param_in, int legID);
    static void stanceLegTg(const aris::dynamic::PlanParamBase &param_in, int legID);
    static void followLegTg(const aris::dynamic::PlanParamBase &param_in, int legID);
    static void bodyEulTg(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

    static BalanceWalkParam bwParam;
    static Pipe<BalanceWalkParam> bwPipe;
    static std::thread bwThread;
};

#endif
