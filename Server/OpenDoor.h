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

#include "Move_Gait.h"

void parseContinueMoveBegin(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
void parseContinueMoveJudge(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
int continueMove(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

void recordOpenDoorData();

enum MoveState
{
    None,
    PointLocate1,
    PointLocate2,
    LocateAjust,
    Forward,
    Backward,
    Rightward,
    Leftward,
    Follow,
    Jump,
    NonJump,
    Downward,
    Pullhandle,
    Pushhandle,
    PrePull,
    PrePush,
    Pull,
    Push,
};

enum PushState
{
    push2Start,
    alignWalk,
    pushWalk,
};

enum PullState
{
    pull2Start,
    circleWalk,
    legWork,
    rotate,
    pullWalk,
};

struct ContinueMoveParam final :public aris::server::GaitParamBase
{
    std::int32_t move_direction;
};

struct ForceTaskParamBase
{
    double bodyPE_last[6];//213 euler angle
    double bodyVel_last[6];//derivative of 213 euler angle, not spatial velocity
    double pEE_last[18];

    double forceInB[6];
    double forceInB_filtered[6];
};

struct OpenDoorParam :public ForceTaskParamBase
{
    double beginBodyPE[6];
    double beginBodyPm[4][4];
    MoveState moveState;
    PushState pushState;
    PullState pullState;
    int ret {0};
    std::int32_t count {0};
    std::int32_t countIter {0};
    Robots::WalkParam walkParam;
    NormalGait::CircleWalkParam circleWalkParam;
    LowpassFilter<6> filter;

    const double toolInR[3] {0,0.08,-0.385};
    double toolInG[3] {0};
    double toolVel[3] {0};

    //MoveState: PointLocation


    //Door Location
    double planeYPR[3] {0};
    double handleLocation[3] {0};

    bool isLastFollow;

    //startPE
    double beginPE[6] {0};


    //now2Start used twice
    double startPE[6] {0};
    const int now2StartCount {1000};

    //MoveState: Follow
    double startPeeInB[18] {0};
    double endPeeInB[18] {0};

    //MoveState: Downward
    bool downwardFlag;
    int downwardCount {0};

    //PushState
    double xNowInG[3] {0};
    double yNowInG[3] {0};
    double now2startDistance[3] {0};
    double now2startDistanceInB[6] {0};
    double now2startDistanceInG[6] {0};
    double now2startPeeDistInG[3] {0};

    //PullState
    const int legWorkCount {1000};

    //pause
    MoveState moveState_last;
    int pauseCount{0};
    bool isPause;
};

extern Pipe<OpenDoorParam> openDoorPipe;
static std::thread openDoorThread;

enum OpenDoorState
{
    Init,
    WaitCmd,
    DoorLocate,
    AdjustD2H,
    HandleLocate,
    TurnHandle,
    PushDoor,
    PullDoor,
    AdjustP2W,
    WalkThrough,
    Quit,
};

enum DoorLocateState
{
    Point1,
    Point2,
    Point3,
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

class OpenDoor
{
public:
    OpenDoor();
    ~OpenDoor();
    static void parseOpenDoor(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
    static int openDoor(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

private:
    static double Fbody[6];
    static double ForceRange[2];
    static double planeVerticalInB[3];
    static double ODBeginPeb[6];
    static double ODBeginPmb[4][4];
    static double adjustBeginPeb[6];
    static double forwardBeginPeb[6];
    static double jumpBeginPeb[6];
    static double nowPeb[6];
    static double nowPmb[4][4];
    static double nowPee[18];

    static bool isContinue;
    static int moveDir[6];
    static bool isPull;
    static bool isLeft;
    static bool isConfirm;
    static bool isJump;
    static OpenDoorState ODState;
    static DoorLocateState DLState;
    static HandleLocateState HLState;
    static OpenDoorParam ODP;

    static void initialize(Robots::RobotBase &robot);
    static int locateDoor(Robots::RobotBase &robot, int count);
    static void GetAdjustParam(Robots::RobotBase &robot);
    static void GenerateBodyMotionByFce(Robots::RobotBase &robot);
    static int adjustRobotD2H(Robots::RobotBase &robot, int count);
    static int locateHandle(Robots::RobotBase &robot, int count);
    static int turnHandle(Robots::RobotBase &robot, int count);


};

#endif
