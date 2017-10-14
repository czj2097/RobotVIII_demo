#include <iostream>
#include <cstring>
#include <iomanip> 
#include <bitset>
#include <cstring>
#include <map>
#include <string>

using namespace std;

#include <aris.h>
#include <Robot_Type_I.h>

#ifdef WIN32
#define rt_printf printf
#endif
#ifdef UNIX
#include "rtdk.h"
#include "unistd.h"
#endif

#include "Move_Gait.h"
#include "TimeOptimalGait.h"


int main(int argc, char *argv[])
{
    /*
    //test IMU//
    double peIMU2Body[6] {0};//{0,0,0,0,-PI/2,PI};
    double pmIMU2Body[4][4] {0};
    double pmBody2IMU[4][4] {0};
    pmIMU2Body[0][0] = -1;
    pmIMU2Body[1][2] = -1;
    pmIMU2Body[2][1] = -1;
    pmIMU2Body[3][3] = 1;
    aris::dynamic::s_pm2pe(*pmIMU2Body,peIMU2Body,"313");
    aris::dynamic::s_inv_pm(*pmIMU2Body,*pmBody2IMU);
    printf("peIMU2Body:%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",peIMU2Body[0],peIMU2Body[1],peIMU2Body[2],peIMU2Body[3],peIMU2Body[4],peIMU2Body[5]);

    double peIMU2IMUGround[6]{0,0,0,PI,0,PI};
    double pmIMU2IMUGround[4][4];
    aris::dynamic::s_pe2pm(peIMU2IMUGround,*pmIMU2IMUGround,"321");

    double peIMUGround2BodyGround[6] {0};//{0,0,0,0,-PI/2,PI};
    double pmIMUGround2BodyGround[4][4] {0};
    pmIMUGround2BodyGround[0][0] = 1;
    pmIMUGround2BodyGround[1][2] = 1;
    pmIMUGround2BodyGround[2][1] = -1;
    pmIMUGround2BodyGround[3][3] = 1;
    aris::dynamic::s_pm2pe(*pmIMUGround2BodyGround,peIMUGround2BodyGround);
    printf("peIMUGround2BodyGround:%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",peIMUGround2BodyGround[0],peIMUGround2BodyGround[1],peIMUGround2BodyGround[2],peIMUGround2BodyGround[3],peIMUGround2BodyGround[4],peIMUGround2BodyGround[5]);

    double pm[4][4] {0};
    aris::dynamic::s_pm_dot_pm(*pmIMUGround2BodyGround,*pmIMU2IMUGround,*pmBody2IMU,*pm);
    printf("pmBody2BodyGround:\n");
    for (int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            printf("%.4f,",pm[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    printf("pmImuGround2BodyGround:\n");
    for (int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            printf("%.4f,",pmIMUGround2BodyGround[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    printf("pmImu2ImuGround:\n");
    for (int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            printf("%.4f,",pmIMU2IMUGround[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    printf("pmBodu2Imu:\n");
    for (int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            printf("%.4f,",pmBody2IMU[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    */

    //FastWalk::TimeOptimalGait1by1();
//    FastWalk::WalkPYAnalyse();
    TimeOptimalGait optimalgait;
    optimalgait.GetOptimalDsByMajorIteration();
    optimalgait.GetOptimalGait2s();
    optimalgait.GetOptimalGait2t();
    optimalgait.OutputData();
    optimalgait.GetNormalGait();

    //FastWalk::JointSpaceWalk jointspacewalker;
    //FastWalk::FastWalkPY pyfastwalker;
    //ForceTask::ForceWalk forcewalker;
    //NormalGait::StartRecordData();

	std::string xml_address;

	if (argc <= 1)
	{
		std::cout << "you did not type in robot name, in this case ROBOT-VIII will start" << std::endl;
		xml_address = "/usr/Robots/resource/Robot_Type_I/Robot_VIII/Robot_VIII.xml";
	}
	else if (std::string(argv[1]) == "III")
	{
		xml_address = "/usr/Robots/resource/Robot_Type_I/Robot_III/Robot_III.xml";
	}
	else if (std::string(argv[1]) == "VIII")
	{
		xml_address = "/home/hex/Desktop/mygit/RobotVIII_demo/resource/Robot_VIII.xml";
	}
	else
	{
		throw std::runtime_error("invalid robot name, please type in III or VIII");
	}

	auto &rs = aris::server::ControlServer::instance();

	rs.createModel<Robots::RobotTypeI>();
	rs.loadXml(xml_address.c_str());
	rs.addCmd("en", Robots::basicParse, nullptr);
	rs.addCmd("ds", Robots::basicParse, nullptr);
	rs.addCmd("hm", Robots::basicParse, nullptr);
	rs.addCmd("rc", Robots::recoverParse, Robots::recoverGait);
	rs.addCmd("wk", Robots::walkParse, Robots::walkGait);
	rs.addCmd("ro", Robots::resetOriginParse, Robots::resetOriginGait);

	rs.addCmd("mwr",NormalGait::parseMoveWithRotate,NormalGait::moveWithRotate);

    rs.addCmd("cmb",ForceTask::parseContinueMoveBegin,ForceTask::continueMove);
    rs.addCmd("cmj",ForceTask::parseContinueMoveJudge,ForceTask::continueMove);
    rs.addCmd("odb",ForceTask::parseOpenDoorBegin,ForceTask::openDoor);
    rs.addCmd("odj",ForceTask::parseOpenDoorJudge,ForceTask::openDoor);
    rs.addCmd("ffd",Robots::walkParse, ForceTask::forceForward);

    //rs.addCmd("jfw",jointspacewalker.parseJointSpaceFastWalk,jointspacewalker.jointSpaceFastWalk);
    //rs.addCmd("fsw",pyfastwalker.parseFastWalkByPY,pyfastwalker.fastWalkByPY);
    //rs.addCmd("fcw",forcewalker.parseForceWalk,forcewalker.forceWalk);

	rs.open();

	rs.setOnExit([&]()
	{
		aris::core::XmlDocument xml_doc;
		xml_doc.LoadFile(xml_address.c_str());
		auto model_xml_ele = xml_doc.RootElement()->FirstChildElement("Model");
		if (!model_xml_ele)throw std::runtime_error("can't find Model element in xml file");
		rs.model().saveXml(*model_xml_ele);

		aris::core::stopMsgLoop();
	});
	aris::core::runMsgLoop();


	return 0;
}
