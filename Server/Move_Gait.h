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
#include "NormalGait.h"

using namespace aris::control;

//static const double PI = 3.141592653589793;

namespace ForceTask
{
    struct ContinueMoveParam final :public aris::server::GaitParamBase
    {
        std::int32_t move_direction;
    };
	void parseContinueMoveBegin(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
	void parseContinueMoveJudge(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
    int continueMove(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

    void forceInit(int count, const double* forceRaw_in, double* forceInF_out);
	int forceForward(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);

    struct ForceTaskParamBase
    {
        double bodyPE_last[6];//213 euler angle
        double bodyVel_last[6];//derivative of 213 euler angle, not spatial velocity
        double pEE_last[18];

        double forceInB[6];
    };
}

namespace FastWalk
{
    void TestGetdJacOverPee();

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
		static void parseFastWalkByPY(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
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

extern Pipe<FastWalk::outputParam> fastWalkPipe;
static std::thread fastWalkThread;

#endif // MOVE_GAIT_H
