#ifndef NORMAL_GAIT_H
#define NORMAL_GAIT_H

#include <aris.h>
#include <Robot_Gait.h>
#include <Robot_Type_I.h>

namespace NormalGait
{
    struct MoveRotateParam final :public aris::server::GaitParamBase
    {
        double targetBodyPE213[6]{0};
        std::int32_t totalCount;
    };
    void parseMoveWithRotate(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
    int moveWithRotate(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

    struct AdjustRcParam final :public aris::server::GaitParamBase
    {
        double distance {0};
        bool isAll {false};
        int legID {0};
        bool actLegs[6];
        std::int32_t totalCount {1500};
    };
    void parseAdjustRc(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
    int adjustRc(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

    struct P2PWalkParam final :public aris::server::GaitParamBase
    {
        double targetPointInB[3] {0};
        std::int32_t totalCount {1000};
        double height {0.05};
    };
    void parseP2PWalk(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
    int p2pWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

    struct CircleWalkParam final :public aris::server::GaitParamBase
    {
        double startAngle {0};
        double radius {1};
        std::int32_t totalCount {1000};
        std::int32_t n {1};
        double h {0.05};
        double beta {0};
        std::int32_t direction {1};//-1 clockwise, 1 anticlockwise
    };
    void parseCircleWalk(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
    int circleWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

    struct LegWorkParam final :public aris::server::GaitParamBase
    {
        double distance {0};
        std::int32_t totalCount {1000};
        double h {0.05};
        double beta {0};
    };
    void parseLegWork(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
    int legWork(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);
}

#endif
