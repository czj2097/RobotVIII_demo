#include <Platform.h>

#include <iostream>
#include <cstring>
#include <iomanip> 
#include <bitset>
#include <map>
#include <string>
#include <memory>
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


using namespace Aris::Core;
using namespace Aris::Control;

#include "Move_Gait.h"

int main()
{
	ForceTask::StartRecordData();

	auto rs = Robots::ROBOT_SERVER::GetInstance();
	rs->CreateRobot<Robots::ROBOT_TYPE_I>();
    rs->LoadXml("/home/hex/Desktop/mygit/RobotVIII_demo/resource/Robot_VIII.xml");

	rs->AddGait("wk", Robots::walk, Robots::parseWalk);
	rs->AddGait("ad", Robots::adjust, Robots::parseAdjust);
	rs->AddGait("ro", Robots::resetOrigin, Robots::parseResetOrigin);
    rs->AddGait("move2",move2,parseMove2);
    rs->AddGait("cmb",ForceTask::continueMove,ForceTask::parseContinueMoveBegin);
    rs->AddGait("cmj",ForceTask::continueMove,ForceTask::parseContinueMoveJudge);
    rs->AddGait("cwf",continuousWalkWithForce,parseCWF);
    rs->AddGait("cwfs",continuousWalkWithForce,parseCWFStop);
    rs->AddGait("fw", Robots::fastWalk, Robots::parseFastWalk);
    //rs->AddGait("sw",swing,parseSwing);
    //rs->AddGait("mr",moveWithRotate,parseMoveWithRotate);
    rs->AddGait("odb",ForceTask::openDoor,ForceTask::parseOpenDoorBegin);
    rs->AddGait("odj",ForceTask::openDoor,ForceTask::parseOpenDoorJudge);

	rs->Start();
	
	std::cout<<"started"<<std::endl;
	
	Aris::Core::RunMsgLoop();

	std::cout << "Program exited!" << std::endl;
	return 0;
}

<<<<<<< HEAD
void startRecordGaitData()
{
	openDoorThread = std::thread([&]()
	{
		struct CM_LAST_PARAM param;
		static std::fstream fileGait;
		std::string name = Aris::Core::logFileName();
		name.replace(name.rfind("log.txt"), std::strlen("openDoor.txt"), "openDoor.txt");
		fileGait.open(name.c_str(), std::ios::out | std::ios::trunc);

		long long count = -1;
		while (1)
		{
			openDoorPipe.RecvInNRT(param);

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

=======
>>>>>>> dev
