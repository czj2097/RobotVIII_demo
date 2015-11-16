#include <Platform.h>

#include <iostream>
#include <cstring>
#include <iomanip> 
#include <bitset>
#include <cstring>
#include <map>
#include <string>

using namespace std;


#ifdef PLATFORM_IS_WINDOWS
#define rt_printf printf
#endif



#include <Aris_Core.h>
#include <Aris_Message.h>
#include <Aris_IMU.h>
#include <Robot_Server.h>
#include <Robot_Gait.h>
#include <Robot_Type_I.h>

using namespace Aris::Core;
#include "Move_Gait.h"


int main()
{
	auto rs = Robots::ROBOT_SERVER::GetInstance();
	rs->CreateRobot<Robots::ROBOT_TYPE_I>();
 	rs->LoadXml("/home/hex/Desktop/RobotVIII_demo/resource/RobotVIII_exhibition.xml");


	rs->AddGait("wk", Robots::walk, Robots::parseWalk);
	rs->AddGait("ad", Robots::adjust, Robots::parseAdjust);
	rs->AddGait("fw", Robots::fastWalk, Robots::parseFastWalk);
	rs->AddGait("ro", Robots::resetOrigin, Robots::parseResetOrigin);

    rs->AddGait("move2",move2,parseMove2);
    rs->AddGait("sw",swing,parseSwing);

    rs->AddGait("mr",moveWithRotate,parseMoveWithRotate);
    //rs->AddGait("cmb",continueMoveBegin,parseContinueMoveBegin);
    rs->AddGait("cmj",continueMoveJudge,parseContinueMoveJudge);

	rs->Start();
	std::cout<<"started"<<std::endl;

	
	
	Aris::Core::RunMsgLoop();

	return 0;
}
