#ifndef BALANCE_H
#define BALANCE_H

#include <aris.h>
#include <Robot_Type_I.h>
#include "GeneralFunc.h"
#include "Move_Gait.h"
#include "EmergencyStop.h"

void GetTargetEulFromAcc(double *planAccInG, double *reqAccInG, double *targetEul);
void GetReqAccInBFromFce(double *fce, double *reqAccInB);

struct BalanceParam final :public aris::server::GaitParamBase
{
    int preCount;
    int totalCount;
    int n;
    double d;
};
void parseBalance(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
int balance(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);


enum BalanceState
{
    Init,
    Balance,
    Stop,
};
struct BallBalanceParam
{
    double fceInB[6];
    double fceInB_filtered[6];
    double bodyPos[3];
    double ballPos[2];
    double tarAcc[2];
    double tarEul[2];
    double actEul[2];
    double imuEul[2];
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
    static BalanceState workState;
    static int countIter;
    static double GetAngleFromAcc(double acc);

    static BallBalanceParam bbParam;
    static Pipe<BallBalanceParam> bbPipe;
    static std::thread bbThread;
};

#endif
