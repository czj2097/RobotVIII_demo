#ifndef MOVE_GAIT_H
#define MOVE_GAIT_H

#include <iostream>
#include <memory>
#include <cstring>
#include <iomanip>
#include <bitset>
#include <map>
#include <string>
#include <stdlib.h>

#include <Aris_Pipe.h>
#include <Aris_Core.h>
#include <Aris_Message.h>
#include <Aris_DynKer.h>
#include <Aris_Motion.h>
#include <Robot_Server.h>
#include <Robot_Gait.h>
#include <Robot_Base.h>

using namespace Aris::Core;
using namespace std;
using namespace Aris::Control;


//static const double PI = 3.141592653589793;

struct MOVES_PARAM :public Robots::GAIT_PARAM_BASE
{
    double targetPee[18]{0};
    double targetBodyPE[6]{0};
    std::int32_t periodCount;
    int comID; //移动的部件（component）序号
    bool isAbsolute{false}; //用于判断移动命令是绝对坐标还是相对坐标
};

struct SWING_PARAM :public Robots::GAIT_PARAM_BASE
{
    double centreP[3]{0};
    double swingRad;
    std::int32_t periodCount;
};

struct MOVEWITHROTATE_PARAM :public Robots::GAIT_PARAM_BASE
{
	double targetBodyPE213[6]{0};
	std::int32_t totalCount;
};

struct CONTINUEMOVE_PARAM :public Robots::GAIT_PARAM_BASE
{
    std::int8_t move_direction[6];
};

//CWF
enum WALK_DIRECTION
{
    STOP,
    FORWARD,
    BACKWARD,
    RIGHTWARD,
    LEFTWARD
};

enum MoveState
{
	None=0,
	PointLocate1=1,
	PointLocate2=2,
	LocateAjust=3,
	Forward=4,
	Backward=5,
	Rightward=6,
	Leftward=7,
	Follow=8,
	Downward=9,
	Upward=10,
	Pullhandle=11,
	Pushhandle=12,
    PrePush=13,
    Push=14,

};

enum PushState
{
    now2Start=0,
    rightWalk=1,
    forwardWalk=2,
};

struct CM_LAST_PARAM
{
	MoveState moveState;
    MoveState moveState_last;
    PushState pushState;

    std::int32_t count;
    std::int32_t countIter{0};
    int pauseCount{0};
	double forceSum[3];
	double forceAvg[3]{0,0,0};
	double force[3];
	double forceYZ;

	const double posLimit[6]{0.3,0.25,0.5,0.524,0.349,0.349};

	double bodyPE_last[6];
	double bodyVel_last[6];

    bool downwardFlag;
    bool pauseFlag;
    int forwardWalkFlag{0};

    double startPE[6];
    double startPm[4][4];
    double handlePE[6];
    double nowPE[6];
    double nowPm[4][4];
    double nowPee[18];
    double realPE[6];
    double now2startDistance[3];
    double now2startDistanceModified[6]{0,0,0,0,0,0};
    double handle2startDistance[3];
    double xNowInG[3];

    const int now2StartCount{2000};
    int downwardCount;
    int ret{0};

    Robots::WALK_PARAM walkParam;

    double pointLocation1[6];
    double pointLocation2[6];
    double pointLocation3[6];
    double planeParam[3];
    double planeVertical[3];
    double planeVerticalInB[3];
    double planeYPR[3]{0,0,0};
    //double forwardDistance;

    double horizontal[3];
    double horizontalInB[3];

    double startPeeInB[18];//for case Follow
    double endPeeInB[18];
    double realPeeInB[18];
    const int followCount{2000};
};

extern PIPE<CM_LAST_PARAM> gaitDataPipe;
static std::thread gaitDataThread;

/*parse function*/
Aris::Core::MSG parseMove2(const std::string &cmd, const map<std::string, std::string> &params);
Aris::Core::MSG parseSwing(const std::string &cmd, const map<std::string, std::string> &params);
Aris::Core::MSG parseMoveWithRotate(const std::string &cmd, const map<std::string, std::string> &params);
Aris::Core::MSG parseContinueMoveBegin(const std::string &cmd, const map<std::string, std::string> &params);
Aris::Core::MSG parseContinueMoveJudge(const std::string &cmd, const map<std::string, std::string> &params);
Aris::Core::MSG parseOpenDoorBegin(const std::string &cmd, const map<std::string, std::string> &params);
Aris::Core::MSG parseOpenDoorJudge(const std::string &cmd, const map<std::string, std::string> &params);
Aris::Core::MSG parseCWF(const std::string &cmd, const std::map<std::string, std::string> &params);
Aris::Core::MSG parseCWFStop(const std::string &cmd, const std::map<std::string, std::string> &params);

/*operation function*/
int move2(Robots::ROBOT_BASE * pRobot, const Robots::GAIT_PARAM_BASE * pParam);
int swing(Robots::ROBOT_BASE * pRobot, const Robots::GAIT_PARAM_BASE * pParam);
int moveWithRotate(Robots::ROBOT_BASE * pRobot, const Robots::GAIT_PARAM_BASE * pParam);
//int continueMoveBegin(Robots::ROBOT_BASE * pRobot, const Robots::GAIT_PARAM_BASE * pParam);
int continueMove(Robots::ROBOT_BASE * pRobot, const Robots::GAIT_PARAM_BASE * pParam);
int openDoor(Robots::ROBOT_BASE * pRobot, const Robots::GAIT_PARAM_BASE * pParam);
int continuousWalkWithForce(Robots::ROBOT_BASE * pRobot, const Robots::GAIT_PARAM_BASE * pParam);
WALK_DIRECTION forceJudge(const double *force, const double *threshold);

#endif // MOVE_GAIT_H
