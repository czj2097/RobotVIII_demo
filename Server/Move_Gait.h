#ifndef MOVE_GAIT_H
#define MOVE_GAIT_H

#include <thread>
#include <functional>
#include <cstdint>
#include <map>

#include <aris.h>
#include <aris_control_pipe.h>
#include <Robot_Gait.h>
#include <Robot_Type_I.h>
#include <sys/time.h>

using namespace aris::control;

//static const double PI = 3.141592653589793;

namespace ForceTask
{
	void parseContinueMoveBegin(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg);
	void parseContinueMoveJudge(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg);
	void parseOpenDoorBegin(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg);
	void parseOpenDoorJudge(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg);
	void forceInit(int count, const double* forceRaw_in, const double* forcePm, double* forceInB_out);
	int continueMove(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);
	int openDoor(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);
	int forceWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

	void StartRecordData();
	void inv3(double * matrix,double * invmatrix);
	void crossMultiply(double * vector_in1, double *vector_in2, double * vector_out);
	void s_admdm(int m, int n, int k, double alpha, double *matrix_in1, double *matrix_in2, double *matrix_out);
	double dotMultiply(double *vector_in1, double *vector_in2);
	double norm(double * vector_in);

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
		Downward,
		Upward,
		Pullhandle,
		Pushhandle,
		PrePush,
		Push,
	};

	enum PushState
	{
		now2Start,
		leftWalk,
		forwardWalk,
	};

	struct maxVelParam
	{
		double bodyVel_last_spatial[6]{0,0,0,0,0,0};
		bool legState{false};
	};
	void getMaxPin(double* maxPin, aris::dynamic::Model &model, ForceTask::maxVelParam &param_in);

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
	};

	struct OpenDoorParam :public ForceTaskParamBase
	{
		MoveState moveState;
		PushState pushState;
		int ret{0};
		std::int32_t count;
		std::int32_t countIter{0};
		Robots::WalkParam walkParam;

		const double toolInR[3]{0,0.08,-0.385};
		double toolInG[3];
		double toolVel[3]{0,0,0};

		//MoveState: PointLocation
		double pointLocation1[6];
		double pointLocation2[6];
		double pointLocation3[6];
		double location[3][3];

		//Door Location
		double planeYPR[3]{0,0,0};
		double handleLocation[3]{0,0,0};
		//startPE
		double beginPE[6];
		double vector0[3];
		double vector1[3];
		double vector2[3];

		//now2Start used twice
		double nowPE[6]; //used in Follow again
		double nowPee[18];
		double startPE[6];
		const int now2StartCount{2000};

		//MoveState: Follow
		double startPeeInB[18];
		double endPeeInB[18];
		const int followCount{2000};

		//MoveState: Downward
		bool downwardFlag;
		int downwardCount;

		//PushState
		double handlePE[6];
		double nowPm[4][4];
		double xNowInG[3];
		double yNowInG[3];
		double now2startDistance[3];
		double now2startDistanceInB[6]{0,0,0,0,0,0};
		double now2startDistanceInG[6]{0,0,0,0,0,0};
		double handle2startDistance[3];

		//pause
		MoveState moveState_last;
		int pauseCount{0};
		bool pauseFlag;
	};
};

extern Pipe<ForceTask::OpenDoorParam> openDoorPipe;
static std::thread openDoorThread;


struct MoveRotateParam final :public aris::server::GaitParamBase
{
	double targetBodyPE213[6]{0};
	std::int32_t totalCount;
};
struct FastWalkByScrewParam final:public aris::server::GaitParamBase
{
	double swingPee[1000][18];
	double bodyForwardAcc{0.5};
	double bodyForwardDec{0.5};
	double bodyForwardVel{0.5};
	std::int32_t n{2};
};

void parseMoveWithRotate(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg);
int moveWithRotate(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);
void parseFastWalkByScrew(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg);
int fastWalkByScrew(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);
void fastTgByScrew();
void fastTg();
void ellipseTrajAnalyse();
void screwInterpolationTraj();
void fastTgByPYAnalyse();
void wkByPYAnalyse();
void maxCal(aris::dynamic::Model &model, double alpha, int legID, double *maxVel, double *maxAcc);
void maxVel();

#endif // MOVE_GAIT_H
