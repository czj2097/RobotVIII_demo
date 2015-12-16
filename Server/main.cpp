#include <Platform.h>

#include <iostream>
#include <cstring>
#include <iomanip> 
#include <bitset>
#include <cstring>
#include <map>
#include <string>
#include <thread>

using namespace std;


#ifdef PLATFORM_IS_WINDOWS
#define rt_printf printf
#endif



#include <Aris_Core.h>
#include <Aris_Pipe.h>
#include <Aris_Message.h>
#include <Aris_IMU.h>
#include <Robot_Server.h>
#include <Robot_Gait.h>
#include <Robot_Type_I.h>
#include <memory>

using namespace Aris::Core;
using namespace Aris::Control;
#include "Move_Gait.h"

void startRecordGaitData();

int main()
{
    startRecordGaitData();

	auto rs = Robots::ROBOT_SERVER::GetInstance();
	rs->CreateRobot<Robots::ROBOT_TYPE_I>();
    rs->LoadXml("/home/hex/Desktop/mygit/RobotVIII_demo/resource/Robot_III.xml");

	rs->AddGait("wk", Robots::walk, Robots::parseWalk);
	rs->AddGait("ad", Robots::adjust, Robots::parseAdjust);
	//rs->AddGait("fw", Robots::fastWalk, Robots::parseFastWalk);
	rs->AddGait("ro", Robots::resetOrigin, Robots::parseResetOrigin);

    //rs->AddGait("move2",move2,parseMove2);
    //rs->AddGait("sw",swing,parseSwing);

    //rs->AddGait("mr",moveWithRotate,parseMoveWithRotate);
    rs->AddGait("cmb",continueMove,parseContinueMoveBegin);
    rs->AddGait("cmj",continueMove,parseContinueMoveJudge);
    //rs->AddGait("odb",openDoor,parseOpenDoorBegin);
    //rs->AddGait("odj",openDoor,parseOpenDoorJudge);
    rs->AddGait("cwf",continuousWalkWithForce,parseCWF);
    rs->AddGait("cwfs",continuousWalkWithForce,parseCWFStop);

	rs->Start();
	
	std::cout<<"started"<<std::endl;
	
	Aris::Core::RunMsgLoop();

	std::cout << "Program exited!" << std::endl;
	return 0;
}

void startRecordGaitData()
{
	gaitDataThread = std::thread([&]()
	{
		struct CM_LAST_PARAM param;
		static std::fstream fileGait;
		std::string name = Aris::Core::logFileName();
		name.replace(name.rfind("log.txt"), std::strlen("gait.txt"), "gait.txt");
		fileGait.open(name.c_str(), std::ios::out | std::ios::trunc);

		long long count = -1;
		while (1)
		{
			gaitDataPipe.RecvInNRT(param);

			//fileGait << ++count << " ";
			fileGait << param.count << "  ";
			fileGait << param.countIter << "  ";
			fileGait << param.moveState << "  ";
			fileGait << param.planeYPR[0] << "  ";
			fileGait << param.walkParam.beta << "  ";


			for (int i = 0; i < 6; i++)
			{
				fileGait << param.bodyPE_last[i] << "  ";
			}
			for (int i = 0; i < 6; i++)
			{
				fileGait << param.walkParam.beginBodyPE[i] << "  ";
			}
			for (int i = 0; i < 3; i++)
			{
				fileGait << param.force[i] << "  ";
			}
			fileGait << std::endl;
		}

		fileGait.close();
	});

}

