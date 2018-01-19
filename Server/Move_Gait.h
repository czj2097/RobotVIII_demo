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

#include "LowPassFilter.h"

using namespace aris::control;

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
    void parseMoveWithRotate(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
	int moveWithRotate(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

	void StartRecordData();
	void inv3(double * matrix,double * invmatrix);
	double norm(double * vector_in);

    void testWorkSpace();

    struct AdjustRcParam final :public aris::server::GaitParamBase
    {
        double distance {0};
        bool isAll {false};
        int legID {0};
        std::int32_t totalCount {1500};
    };
    void parseAdjustRc(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
    int adjustRc(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

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

namespace ForceTask
{
    void forceInit(int count, const double* forceRaw_in, double* forceInF_out);
	int forceForward(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

	struct ForceWalkParam final:public aris::server::GaitParamBase
	{
	};
	class ForceWalk
	{
	public:
		ForceWalk();
		~ForceWalk();
        static void parseForceWalk(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
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
	void maxCal(aris::dynamic::Model &model, double alpha, int legID, double *maxVel, double *maxAcc);
    void maxVel();
	struct maxVelParam
	{
		double bodyVel_last_spatial[6]{0,0,0,0,0,0};
		bool legState{false};
	};
    void getMaxPin(double* maxPin, aris::dynamic::Model &model, maxVelParam &param_in);//unfinished, useless


    void FastWalkPYAnalyse();
    void WalkPYAnalyse();
    struct OfflineGaitParam final:public aris::server::GaitParamBase
	{
        std::int32_t n {2};
        int totalCountPY {900};
        int totalCountCZJ {3258};
	};
    class OfflineGait
	{
	public:
        OfflineGait();
        ~OfflineGait();
        static void parseFastWalkByPY(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
		static int fastWalkByPY(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

        static void parseFastWalkByCZJ(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
        static int fastWalkByCZJ(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

	private:
		static double pIn_acc[900][18];
		static double pIn_const[1800][18];
		static double pIn_dec[900][18];

        static double pIn_entire[3258][18];
	};


    void screwInterpolationTraj();
	struct JointSpaceWalkParam final:public aris::server::GaitParamBase
	{
	};
	struct outputParam
	{
		double outputPin[18];
		double outputPee[18];
	};
	class JointSpaceWalk
	{
		public:
			JointSpaceWalk();
			~JointSpaceWalk();
            static void parseJointSpaceFastWalk(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
			static int jointSpaceFastWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

		private:
			static void swingLegTg(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in, int legID);
			static void stanceLegTg(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in, int legID);

			static double bodyAcc;
			static double bodyDec;
			static int totalCount;
			static double height;
			static double beta;
			static NormalGait::WalkState walkState;

			static double beginPee[18];
			static double beginVel;
			static double endPee[18];
			static double endVel;
			static double distance;

			static NormalGait::GaitPhase gaitPhase[6];//swing true, stance false
			static bool constFlag;
	};
}

extern Pipe<FastWalk::outputParam> fastWalkPipe;
static std::thread fastWalkThread;

#endif // MOVE_GAIT_H
