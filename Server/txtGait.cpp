#include "txtGait.h"

double TxtGait::Peb[4000][6];
double TxtGait::Pee[4000][18];
double TxtGait::Pin[4000][18];

void TxtGait::parseTxtGait(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
{
    TxtGaitParam param;
    param.totalCount=4000;

    aris::dynamic::dlmread("./Peb.txt",*Peb);
    aris::dynamic::dlmread("./Pee.txt",*Pee);
    aris::dynamic::dlmread("./Pin.txt",*Pin);

    msg.copyStruct(param);

    std::cout<<"finished parse"<<std::endl;
}
int TxtGait::txtGait(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
{
    auto &robot = static_cast<Robots::RobotBase &>(model);
    auto &param = static_cast<const TxtGaitParam &>(param_in);

    robot.SetPeb(*Peb+6*param.count);
    robot.SetPee(*Pee+18*param.count);
    robot.SetPin(*Pin+18*param.count);

    return param.count-param.totalCount+1;
}
