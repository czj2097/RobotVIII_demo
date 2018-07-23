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
struct BallBalanceParam final:public aris::server::GaitParamBase{};

class BallBalance
{
public:
    static void parseBallBalance(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
    static int ballBalance(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

private:
    static BalanceState workState;
    static int countIter;
    static double GetAngleFromAcc(double acc);
};

#endif
