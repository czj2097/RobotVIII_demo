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
#include "rtdk.h"
#include "GeneralFunc.h"

using namespace aris::control;

//static const double PI = 3.141592653589793;
namespace NormalGait
{
	enum WalkState
	{
		Init,
        ForwardAcc,
        TurnAcc,
        Const,
		Stop,
	};
	enum GaitPhase
	{
		Swing,
		Stance,
		Follow,
	};

	struct MoveRotateParam final :public aris::server::GaitParamBase
	{
		double targetBodyPE213[6]{0};
		std::int32_t totalCount;
	};
	void parseMoveWithRotate(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg);
	int moveWithRotate(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

	void StartRecordData();
}

namespace ForceTask
{
	void parseContinueMoveBegin(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg);
	void parseContinueMoveJudge(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg);
    int continueMove(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

    void parseOpenDoorBegin(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg);
	void parseOpenDoorJudge(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg);
    int openDoor(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

    void forceInit(int count, const double* forceRaw_in, double* forceInF_out);
	int forceForward(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

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
        adjustLeg,
		forwardWalk,
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


	struct ForceWalkParam final:public aris::server::GaitParamBase
	{
	};
	class ForceWalk
	{
	public:
		ForceWalk();
		~ForceWalk();
		static void parseForceWalk(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg);
		static int forceWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

	private:
        static double forwardAcc;
        static double turnAcc;
		static int totalCount;
        static int totalCount_tmp;
		static double height;
        static double height_tmp;
        static double alpha;
        static double alpha_tmp;

        static NormalGait::WalkState walkState;
        static bool constFlag;
        static double beginXVel;
        static double endXVel;
        static double beginZVel;
        static double endZVel;
        static double beginOmega;
        static double endOmega;

        static double beginPeb[6];
        static double pEB[6];
        static NormalGait::GaitPhase gaitPhase[6];
        static double swingPee[18];
        static double swingBeginPee[18];
        static double swingEndPee[18];
        static double stancePee[18];
        static double stanceBeginPee[18];
        static double stanceEndPee[18];
        static double followBeginPee[18]; 
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

        static double forceRaw[36];
        static double forceInF[36];
        static double forceSum[36];
        static double forceAvg[36];

        static void forceInit(int count, int legID);
        static void swingLegTg(const aris::dynamic::PlanParamBase &param_in, int legID);
        static void stanceLegTg(const aris::dynamic::PlanParamBase &param_in, int legID);
    };
}

namespace FastWalk
{

    void TestGetdJacOverPee();
    void TimeOptimalGait3();//3 swing legs calculated together
    void TimeOptimalGait1by1();//3 swing legs calculated one by one

    //calculate max vel & cal at any position in workspace, vel finished, cal unfinished
    void maxCal(aris::dynamic::Model &model, double *direction, int legID, double vLmt, double *maxVel);
    void maxVel();
	struct maxVelParam
	{
		double bodyVel_last_spatial[6]{0,0,0,0,0,0};
		bool legState{false};
	};
    void getMaxPin(double* maxPin, aris::dynamic::Model &model, maxVelParam &param_in);//unfinished, useless


    void FastWalkPYAnalyse();
    void WalkPYAnalyse();
	struct FastWalkByPYParam final:public aris::server::GaitParamBase
	{
		std::int32_t n{2};
		int totalCount{900};
	};
	class FastWalkPY
	{
	public:
		FastWalkPY();
		~FastWalkPY();
		static void parseFastWalkByPY(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg);
		static int fastWalkByPY(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

	private:
		static double pIn_acc[900][18];
		static double pIn_const[1800][18];
		static double pIn_dec[900][18];
    };

    struct outputParam
    {
        double outputPin[18];
        double outputPee[18];
    };

}

extern Pipe<ForceTask::OpenDoorParam> openDoorPipe;
static std::thread openDoorThread;
extern Pipe<FastWalk::outputParam> fastWalkPipe;
static std::thread fastWalkThread;

#endif // MOVE_GAIT_H
