#include <iostream>
#include <cstring>
#include <iomanip> 
#include <bitset>
#include <cstring>
#include <map>
#include <string>

using namespace std;

#include <aris.h>
#include <Robot_Type_III.h>
#include <Basic_Gait.h>

#ifdef WIN32
#define rt_printf printf
#endif
#ifdef UNIX
#include "rtdk.h"
#include "unistd.h"
#endif

#include "Move_Gait.h"


int main(int argc, char *argv[])
{
    /*
    //test IMU//
    double peIMU2Body[6];//{0,0,0,0,-PI/2,PI};
    double pmIMU2Body[4][4];
    double pmBody2IMU[4][4];
    std::fill_n(&pmIMU2Body[0][0], 16, 0);
    pmIMU2Body[0][1] = -1;
    pmIMU2Body[1][2] = -1;
    pmIMU2Body[2][0] = 1;
    pmIMU2Body[3][3] = 1;
    aris::dynamic::s_pm2pe(*pmIMU2Body,peIMU2Body,"313");
    aris::dynamic::s_inv_pm(*pmIMU2Body,*pmBody2IMU);
    printf("%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",peIMU2Body[0],peIMU2Body[1],peIMU2Body[2],peIMU2Body[3],peIMU2Body[4],peIMU2Body[5]);

    double peIMU2IMUGround[6]{0,0,0,PI,0,PI};
    double pmIMU2IMUGround[4][4];
    aris::dynamic::s_pe2pm(peIMU2IMUGround,*pmIMU2IMUGround,"321");

    double peIMUGround2BodyGround[6];//{0,0,0,0,-PI/2,PI};
    double pmIMUGround2BodyGround[4][4];
    std::fill_n(&pmIMUGround2BodyGround[0][0], 16, 0);
    pmIMUGround2BodyGround[0][1] = -1;
    pmIMUGround2BodyGround[1][2] = 1;
    pmIMUGround2BodyGround[2][0] = -1;
    pmIMUGround2BodyGround[3][3] = 1;
    aris::dynamic::s_pm2pe(*pmIMUGround2BodyGround,peIMUGround2BodyGround);
    printf("%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",peIMUGround2BodyGround[0],peIMUGround2BodyGround[1],peIMUGround2BodyGround[2],peIMUGround2BodyGround[3],peIMUGround2BodyGround[4],peIMUGround2BodyGround[5]);

    double pm[4][4];
    aris::dynamic::s_pm_dot_pm(*pmIMUGround2BodyGround,*pmIMU2IMUGround,*pmBody2IMU,*pm);
    for (int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            printf("%.4f,",pm[i][j]);
        }
        printf("\n");
    }

    for (int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            printf("%.4f,",pmIMUGround2BodyGround[i][j]);
        }
        printf("\n");
    }

    for (int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            printf("%.4f,",pmIMU2IMUGround[i][j]);
        }
        printf("\n");
    }

    for (int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            printf("%.4f,",pmBody2IMU[i][j]);
        }
        printf("\n");
    }*/

    //FastWalk::wkByPYAnalyse();

    FastWalk::JointSpaceWalk jointspacewalker;
	FastWalk::FastWalkPY pyfastwalker;
    ForceTask::ForceWalk forcewalker;
	NormalGait::StartRecordData();
	std::string xml_address;

    if (argc <= 1)
    {
        std::cout << "you did not type in robot name, in this case ROBOT-IX will start" << std::endl;
        xml_address = "/usr/Robots/resource/Robot_Type_III/Robot_IX/Robot_IX.xml";
    }
    else if (std::string(argv[1]) == "IX")
    {
        xml_address = "/home/hex/Desktop/mygit/RobotVIII_demo/resource/Robot_IX.xml";
    }
    else if (std::string(argv[1]) == "VIII")
    {
        xml_address = "/home/hex/Desktop/mygit/RobotVIII_demo/resource/Robot_VIII.xml";
    }
    else if (std::string(argv[1]) == "III")
    {
        xml_address = "/home/hex/Desktop/mygit/RobotVIII_demo/resource/Robot_III.xml";
    }
	else
	{
        throw std::runtime_error("invalid robot name, please type in III, VIII or IX");
	}

	auto &rs = aris::server::ControlServer::instance();

    rs.createModel<Robots::RobotTypeIII>();
	rs.loadXml(xml_address.c_str());
	rs.addCmd("en", Robots::basicParse, nullptr);
	rs.addCmd("ds", Robots::basicParse, nullptr);
	rs.addCmd("hm", Robots::basicParse, nullptr);
	rs.addCmd("hmsw", Robots::basicParse, nullptr);
	rs.addCmd("rc", Robots::recoverParse, Robots::recoverGait);
	rs.addCmd("wk", Robots::walkParse, Robots::walkGait);
	rs.addCmd("ro", Robots::resetOriginParse, Robots::resetOriginGait);
    rs.addCmd("ec", Robots::Gait::extendChainParse, Robots::Gait::extendChainGait);
    //waist cmd
    rs.addCmd("rcw", Robots::Gait::recoverWaistParse, Robots::Gait::recoverWaistGait);
    rs.addCmd("aw", Robots::Gait::adjustWaistParse, Robots::Gait::adjustWaistGait);


	rs.addCmd("mwr",NormalGait::parseMoveWithRotate,NormalGait::moveWithRotate);
    //rs.addCmd("swk",NormalGait::parseSpecialWalk,NormalGait::specialWalk);

    rs.addCmd("cmb",ForceTask::parseContinueMoveBegin,ForceTask::continueMove);
    rs.addCmd("cmj",ForceTask::parseContinueMoveJudge,ForceTask::continueMove);
    rs.addCmd("odb",ForceTask::parseOpenDoorBegin,ForceTask::openDoor);
    rs.addCmd("odj",ForceTask::parseOpenDoorJudge,ForceTask::openDoor);
    rs.addCmd("ffd",Robots::walkParse, ForceTask::forceForward);

    rs.addCmd("jfw",jointspacewalker.parseJointSpaceFastWalk,jointspacewalker.jointSpaceFastWalk);
    rs.addCmd("fsw",pyfastwalker.parseFastWalkByPY,pyfastwalker.fastWalkByPY);
    rs.addCmd("fcw",forcewalker.parseForceWalk,forcewalker.forceWalk);

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
