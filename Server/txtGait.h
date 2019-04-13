#ifndef TXT_GAIT_H
#define TXT_GAIT_H

#include <aris.h>
#include <Robot_Gait.h>
#include <Robot_Type_I.h>

struct TxtGaitParam final:public aris::server::GaitParamBase
{
    int totalCount;
};

class TxtGait
{
public:
    static void parseTxtGait(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
    static int txtGait(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

private:
    static double Peb[4000][6];
    static double Pee[4000][18];
    static double Pin[4000][18];
};


#endif
