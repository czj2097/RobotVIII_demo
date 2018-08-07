#ifndef OPENDOOR_H
#define OPENDOOR_H

#include <thread>
#include <functional>
#include <cstdint>
#include <map>

#include <aris.h>
#include <aris_control_pipe.h>
#include <Robot_Gait.h>
#include <Robot_Type_I.h>
#include <sys/time.h>

#include "LowPassFilter.h"
#include "GeneralFunc.h"
#include "NormalGait.h"
#include "EmergencyStop.h"

using namespace aris::control;

struct OpenDoorParam
{
    double PebInG[6];
    double PeeInG[18];
    double PebInMak[6];
    double PeeInMak[18];
    double forceInB[6];
    double forceInB_filtered[6];
    std::int32_t count {0};
    std::int32_t countIter {0};
};

//extern Pipe<OpenDoorParam> openDoorPipe;
enum OpenDoorState
{
    Init,
    WaitCmd,
    DoorLocate,
    AdjustD2H,
    HandleLocate,
    TurnHandle,
    PushHandle,
    PullHandle,
    AdjustPush2W,
    AdjustPull2W,
    WalkThrough,
    PullWork,
    Stop,
};

enum DoorLocateState
{
    Point1,
    Point2,
    Point3,
    PrepareD2H,
};

enum HandleLocateState
{
    Forward,
    Backward,
    Leftward,
    RightWard,
    Follow,
    Jump,
    NonJump,
    Downward,
};

enum PullWorkState
{
    CircleWalk,
    LegWork,
    Rotate,
    PullWalk,
};

enum WalkThroughState
{
    AlignWalk,
    AdjustWidth,
    PushWalk,
};

class OpenDoor
{
public:
    static void parseOpenDoor(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
    static int openDoor(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);
    static void recordData();

private:
    static double Fbody[6];
    static double ForceRange[2];
    static double bodyLastVel[6];//derivative of 213 euler angle, not spatial velocity
    static double planeVerticalInB[3];
    static double ODBeginPeb[6];
    static double ODBeginPmb[4][4];
    static double adjustBeginPeb[6];
    static double jumpBeginPeb[6];
    static double nowPeb[6];
    static double nowPmb[4][4];
    static double nowPee[18];
    static double now2startPebDistInG[6];
    static double now2startPeeDistInG[3];

    static int moveDir[6];
    static bool isPull;
    static bool isLeft;
    static bool isJump;

    static OpenDoorState ODState;
    static DoorLocateState DLState;
    static HandleLocateState HLState;
    static PullWorkState PWState;
    static WalkThroughState WTState;

    static OpenDoorParam ODP;
    static Robots::WalkParam walkParam;
    static NormalGait::AdjustRcParam arcParam;
    static NormalGait::CircleWalkParam circleWalkParam;
    static LowpassFilter<6> filter;
    static Pipe<OpenDoorParam> odPipe;
    static std::thread odThread;

    static void initialize(Robots::RobotBase &robot);
    static int locateDoor(Robots::RobotBase &robot, int count);
    static void GetAdjustD2HParam(Robots::RobotBase &robot);
    static void GenerateBodyMotionByFce(Robots::RobotBase &robot);
    static int adjustRobotD2H(Robots::RobotBase &robot, int count);
    static int locateHandle(Robots::RobotBase &robot, int count);
    static int turnHandle(Robots::RobotBase &robot, int count);
    static int pushHandle(Robots::RobotBase &robot, int count);
    static int pullHandle(Robots::RobotBase &robot, int count);
    static void GetAdjustPush2WParam(Robots::RobotBase &robot);
    static void GetAdjustPull2WParam(Robots::RobotBase &robot);
    static int adjustRobotPush2W(Robots::RobotBase &robot, int count);
    static int adjustRobotPull2W(Robots::RobotBase &robot, int count);
    static void GetPushWalkParam();
    static void GetPullWalkParam();
    static int pushWalkThrough(Robots::RobotBase &robot, int count);
    static int pullDoorWork(Robots::RobotBase &robot, int count);
};

#endif
