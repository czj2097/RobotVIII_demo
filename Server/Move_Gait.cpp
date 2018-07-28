#include "Move_Gait.h"

#ifdef UNIX
#include "rtdk.h"
#endif
#ifdef WIN32
#define rt_printf printf
#endif

#include <cstring>
#include <cmath>
#include <algorithm>
#include <memory>

Pipe<FastWalk::outputParam> fastWalkPipe(true);

namespace NormalGait
{
    void parseMoveWithRotate(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
	{
		MoveRotateParam param;

		for(auto &i:params)
		{
			if(i.first=="u")
			{
                param.targetBodyPE213[0]=std::stod(i.second);
			}
			else if(i.first=="v")
			{
                param.targetBodyPE213[1]=std::stod(i.second);
			}
			else if(i.first=="w")
			{
                param.targetBodyPE213[2]=std::stod(i.second);
			}
			else if(i.first=="yaw")
			{
                param.targetBodyPE213[3]=std::stod(i.second)*PI/180;
			}
			else if(i.first=="pitch")
			{
                param.targetBodyPE213[4]=std::stod(i.second)*PI/180;
			}
			else if(i.first=="roll")
			{
                param.targetBodyPE213[5]=std::stod(i.second)*PI/180;
			}
			else if(i.first=="totalCount")
			{
                param.totalCount=std::stoi(i.second);
			}
			else
			{
				std::cout<<"parse failed"<<std::endl;
			}
		}

		msg.copyStruct(param);

		std::cout<<"finished parse"<<std::endl;
	}

	int moveWithRotate(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
    {
		auto &robot = static_cast<Robots::RobotBase &>(model);
		auto &param = static_cast<const MoveRotateParam &>(param_in);

        /*****test and calibrate imu*****/
        double peImuGrnd2BodyGrnd[6] {0,0,0,PI/6,-PI/2,0};
        double pmImuGrnd2BodyGrnd[16];
        double pmBody2Imu[16];
        double peImu2ImuGrnd[6] {0};
        peImu2ImuGrnd[3]=param.imu_data->at(0).euler[2];
        peImu2ImuGrnd[4]=param.imu_data->at(0).euler[1];
        peImu2ImuGrnd[5]=param.imu_data->at(0).euler[0];
        double pmImu2ImuGrnd[16];
        double bodyPe[6];//213
        double bodyPm[16];
        aris::dynamic::s_pe2pm(peImuGrnd2BodyGrnd,pmImuGrnd2BodyGrnd,"213");
        aris::dynamic::s_pe2pm(peImu2ImuGrnd,pmImu2ImuGrnd,"321");
        aris::dynamic::s_inv_pm(pmImuGrnd2BodyGrnd,pmBody2Imu);
        aris::dynamic::s_pm_dot_pm(pmImuGrnd2BodyGrnd,pmImu2ImuGrnd,pmBody2Imu,bodyPm);
        aris::dynamic::s_pm2pe(bodyPm,bodyPe,"213");
        printf("bodyPe:%f,%f,%f\n",bodyPe[3],bodyPe[4],bodyPe[5]);
        /*************************************/

		static double beginBodyPE213[6];
		static double pEE[18];
		if(param.count==0)
		{
			robot.GetPeb(beginBodyPE213,"213");
			robot.GetPee(pEE);
		}

		double realBodyPE213[6];
		for(int i=0;i<6;i++)
		{
			double s = -(param.targetBodyPE213[i] / 2)*cos(PI * (param.count + 1) / param.totalCount ) + param.targetBodyPE213[i] / 2;
			realBodyPE213[i]=beginBodyPE213[i]+s; //target of current ms
		}

		robot.SetPeb(realBodyPE213,"213");
		robot.SetPee(pEE);

		return param.totalCount - param.count - 1;
	}

    void parseAdjustRc(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
    {
        AdjustRcParam param;

        for(auto &i:params)
        {
            if(i.first=="distance")
            {
                param.distance=std::stod(i.second);
            }
            else if(i.first=="all")
            {
                param.isAll=true;
            }
            else if(i.first=="leg")
            {
                param.legID=std::stoi(i.second);
            }
            else if(i.first=="totalCount")
            {
                param.totalCount=std::stoi(i.second);
            }
            else
            {
                std::cout<<"parse failed"<<std::endl;
            }
        }

        msg.copyStruct(param);

        std::cout<<"finished parse"<<std::endl;
    }

    int adjustRc(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
    {
        auto &robot = static_cast<Robots::RobotBase &>(model);
        auto &param = static_cast<const AdjustRcParam &>(param_in);

        static double StartPeb[6];
        static double StartPee[18];
        double currentPee[18];
        if(param.count==0)
        {
            robot.GetPeb(StartPeb);
            robot.GetPee(StartPee);
        }

        if(param.isAll==true)
        {
            for(int i=0;i<3;i++)//0 2 4
            {
                if(param.count<param.totalCount)
                {
                    currentPee[6*i]=StartPee[6*i]+param.distance/2*sin(PI/6+2*i*PI/3)*(1-cos(param.count*PI/param.totalCount));
                    currentPee[6*i+1]=StartPee[6*i+1]+0.025*(1-cos(param.count*2*PI/param.totalCount));
                    currentPee[6*i+2]=StartPee[6*i+2]+param.distance/2*cos(PI/6+2*i*PI/3)*(1-cos(param.count*PI/param.totalCount));

                    currentPee[6*i+3]=StartPee[6*i+3];
                    currentPee[6*i+4]=StartPee[6*i+4];
                    currentPee[6*i+5]=StartPee[6*i+5];
                }
                else
                {
                    currentPee[6*i]=StartPee[6*i]+param.distance*sin(PI/6+2*i*PI/3);
                    currentPee[6*i+1]=StartPee[6*i+1];
                    currentPee[6*i+2]=StartPee[6*i+2]+param.distance*cos(PI/6+2*i*PI/3);

                    currentPee[6*i+3]=StartPee[6*i+3]+param.distance/2*sin(PI/2-2*i*PI/3)*(1-cos((param.count-param.totalCount)*PI/param.totalCount));
                    currentPee[6*i+4]=StartPee[6*i+4]+0.025*(1-cos((param.count-param.totalCount)*2*PI/param.totalCount));
                    currentPee[6*i+5]=StartPee[6*i+5]+param.distance/2*cos(PI/2-2*i*PI/3)*(1-cos((param.count-param.totalCount)*PI/param.totalCount));
                }
            }
        }
        else
        {
            double k[6] {PI/6,PI/2,5*PI/6,-PI/6,3*PI/2,-5*PI/6};
            for(int i=0;i<6;i++)
            {
                if(i==param.legID)
                {
                    currentPee[3*i]=StartPee[3*i]+param.distance/2*sin(k[i])*(1-cos(param.count*PI/param.totalCount));
                    currentPee[3*i+1]=StartPee[3*i+1]+0.025*(1-cos(param.count*2*PI/param.totalCount));
                    currentPee[3*i+2]=StartPee[3*i+2]+param.distance/2*cos(k[i])*(1-cos(param.count*PI/param.totalCount));
                }
                else
                {
                    currentPee[3*i]=StartPee[3*i];
                    currentPee[3*i+1]=StartPee[3*i+1];
                    currentPee[3*i+2]=StartPee[3*i+2];
                }
            }
        }

        robot.SetPeb(StartPeb);
        robot.SetPee(currentPee);
        if(param.isAll==true)
        {
            return 2*param.totalCount - param.count - 1;
        }
        else
        {
            return param.totalCount - param.count - 1;
        }
    }

    void parseP2PWalk(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
    {
        P2PWalkParam param;

        for(auto &i:params)
        {
            if(i.first=="x")
            {
                param.targetPointInB[0]=std::stod(i.second);
            }
            else if(i.first=="y")
            {
                param.targetPointInB[1]=std::stod(i.second);
            }
            else if(i.first=="z")
            {
                param.targetPointInB[2]=std::stod(i.second);
            }
            else if(i.first=="totalCount")
            {
                param.totalCount=std::stoi(i.second);
            }
            else if(i.first=="height")
            {
                param.height=std::stod(i.second);
            }
            else
            {
                std::cout<<"parse failed"<<std::endl;
            }
        }

        msg.copyStruct(param);

        std::cout<<"finished parse"<<std::endl;
    }

    int p2pWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
    {
        auto &robot = static_cast<Robots::RobotBase &>(model);
        auto &param = static_cast<const P2PWalkParam &>(param_in);

        static double beginPeb[6] {0};
        if(param.count==0)
        {
            robot.GetPeb(beginPeb);
        }
        double Peb[6] {0};
        double dist = sqrt(param.targetPointInB[0]*param.targetPointInB[0]+param.targetPointInB[2]*param.targetPointInB[2]);

        Robots::WalkParam walkParam;
        walkParam.alpha=atan2(param.targetPointInB[0],param.targetPointInB[2])+PI;
        walkParam.beta=0;
        walkParam.h=param.height;
        if(dist<=0.2)
        {
            walkParam.n=1;
        }
        else
        {
            walkParam.n=(int)((dist-0.2)/0.4)+2;
        }
        walkParam.d=dist/(walkParam.n-0.5);
        walkParam.count=param.count;
        walkParam.totalCount=param.totalCount;

        Robots::walkGait(robot,walkParam);
        robot.GetPeb(Peb);
        Peb[1]=beginPeb[1]+param.targetPointInB[1]*(-cos(param.count/(2*walkParam.n*param.totalCount)*PI)+1)/2;

        robot.SetPeb(Peb);

        return 2 * walkParam.n * param.totalCount - param.count - 1;
    }

	void StartRecordData()
	{
        fastWalkThread = std::thread([&]()
		{
            struct FastWalk::outputParam param;
			static std::fstream fileGait;
			std::string name = aris::core::logFileName();
            name.replace(name.rfind("log.txt"), std::strlen("fastWalk.txt"), "fastWalk.txt");
			fileGait.open(name.c_str(), std::ios::out | std::ios::trunc);

			long long count = -1;
			while (1)
			{
                fastWalkPipe.recvInNrt(param);

				//fileGait << ++count << " ";
				for (int i = 0; i < 18; i++)
				{
                    fileGait << param.outputPee[i] << "  ";
				}
				for (int i = 0; i < 18; i++)
				{
                    fileGait << param.outputPin[i] << "  ";
				}
				fileGait << std::endl;
			}

			fileGait.close();
		});
	}

    void testWorkSpace()
    {
        Robots::RobotTypeI rbt;
        rbt.loadXml("../../resource/Robot_VIII.xml");

        double Peb[6] {0};
        double initPee[18] {-0.6*sin(PI/6), -0.58, -0.6*cos(PI/6),
                        -0.6,           -0.58,  0,
                        -0.6*sin(PI/6), -0.58,  0.6*cos(PI/6),
                         0.6*sin(PI/6), -0.58, -0.6*cos(PI/6),
                         0.6,           -0.58,  0,
                         0.6*sin(PI/6), -0.58,  0.6*cos(PI/6)};

        double Pee[18] {-0.6*sin(PI/6), -0.58, -0.6*cos(PI/6),
                        -0.6,           -0.58,  0,
                        -0.6*sin(PI/6), -0.58,  0.6*cos(PI/6),
                         -0.1,          -0.58, -0.5,
                         0.6,           -0.58,  0,
                         0.6*sin(PI/6), -0.58,  0.6*cos(PI/6)};
        double Pin[3];

        rbt.SetPeb(Peb);
        rbt.SetPee(Pee);

        rbt.pLegs[3]->GetPin(Pin);

        printf("Pin:%.4f, %.4f, %.4f\n", Pin[0],Pin[1],Pin[2]);
    }

    void printEverything(Robots::RobotBase &robot)
    {
        double printPmb[16];
        robot.GetPmb(printPmb);
        rt_printf("Pmb[16]={%.4f,%.4f,%.4f,%.4f,\n%.4f,%.4f,%.4f,%.4f\n%.4f,%.4f,%.4f,%.4f\n%.4f,%.4f,%.4f,%.4f}\n",
                  printPmb[0],printPmb[1],printPmb[2],printPmb[3],
                printPmb[4],printPmb[5],printPmb[6],printPmb[7],
                printPmb[8],printPmb[9],printPmb[10],printPmb[11],
                  printPmb[12],printPmb[13],printPmb[14],printPmb[15]);
        double printPeb[6];
        double printPee[18];
        robot.GetPeb(printPeb);
        robot.GetPee(printPee);
        rt_printf("Peb[6]={%.4f,%.4f,%.4f,%.4f,%.4f,%.4f}\nPee[18]={%.4f,%.4f,%.4f,\n%.4f,%.4f,%.4f,\n%.4f,%.4f,%.4f,\n%.4f,%.4f,%.4f,\n%.4f,%.4f,%.4f,\n%.4f,%.4f,%.4f}\n",
                  printPeb[0],printPeb[1],printPeb[2],printPeb[3],printPeb[4],printPeb[5],
                printPee[0],printPee[1],printPee[2],printPee[3],printPee[4],printPee[5],
                printPee[6],printPee[7],printPee[8],printPee[9],printPee[10],printPee[11],
                printPee[12],printPee[13],printPee[14],printPee[15],printPee[16],printPee[17]);
        double printPin[18];
        robot.GetPin(printPin);
        rt_printf("Pin[18]={%.4f,%.4f,%.4f,\n%.4f,%.4f,%.4f,\n%.4f,%.4f,%.4f,\n%.4f,%.4f,%.4f,\n%.4f,%.4f,%.4f,\n%.4f,%.4f,%.4f}\n",
                  printPin[0],printPin[1],printPin[2],printPin[3],printPin[4],printPin[5],
                  printPin[6],printPin[7],printPin[8],printPin[9],printPin[10],printPin[11],
                  printPin[12],printPin[13],printPin[14],printPin[15],printPin[16],printPin[17]);
    }

    void parseCircleWalk(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
    {
        CircleWalkParam param;

        for(auto &i:params)
        {
            if(i.first=="startangle")
            {
                param.startAngle=std::stod(i.second);
            }
            else if(i.first=="radius")
            {
                param.radius=std::stod(i.second);
            }
            else if(i.first=="totalCount")
            {
                param.totalCount=std::stoi(i.second);
            }
            else if(i.first=="n")
            {
                param.n=std::stoi(i.second);
            }
            else if(i.first=="h")
            {
                param.h=std::stod(i.second);
            }
            else if(i.first=="beta")
            {
                param.beta=std::stod(i.second);
            }
            else if(i.first=="direction")
            {
                param.direction=std::stoi(i.second);
            }
            else
            {
                std::cout<<"parse failed"<<std::endl;
            }
        }

        msg.copyStruct(param);

        std::cout<<"finished parse"<<std::endl;
    }

    int circleWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
    {
        auto &robot = static_cast<Robots::RobotBase &>(model);
        auto &param = static_cast<const CircleWalkParam &>(param_in);

        static aris::dynamic::FloatMarker beginMak{ robot.ground() };
        static double beginPee[18];

        if (param.count%param.totalCount == 0)
        {
            beginMak.setPrtPm(*robot.body().pm());
            beginMak.update();
            robot.GetPee(beginPee, beginMak);
        }

        double Peb[6] {0};
        double Pee[18] {0};
        memcpy(Pee,beginPee,18*sizeof(double));

        int period_count = param.count%param.totalCount;
        int half_period = param.count/param.totalCount;
        int period = param.count/(2*param.totalCount);
        double bodyStartAngle=param.startAngle+(0.5*half_period-0.25)*param.beta*param.direction;
        double legStartAngle0=param.startAngle+(period-0.5)*param.beta*param.direction;
        double legStartAngle1=param.startAngle+period*param.beta*param.direction;

        if(param.count<param.totalCount)//acc
        {
            Peb[0]=param.radius*(sin(param.startAngle)-sin(param.startAngle+0.25*param.direction*param.beta*param.count*param.count/param.totalCount/param.totalCount));
            Peb[2]=param.radius*(cos(param.startAngle)-cos(param.startAngle+0.25*param.direction*param.beta*param.count*param.count/param.totalCount/param.totalCount));

            for(int i=0;i<3;i++)
            {
                if(param.direction==-1)
                {
                    Pee[6*i]=beginPee[6*i]+param.radius*(sin(param.startAngle)-sin(param.startAngle+0.5*param.direction*param.beta))/2*(1-cos(param.count*PI/param.totalCount));
                    Pee[6*i+1]=beginPee[6*i+1]+0.05*(1-cos(param.count*2*PI/param.totalCount));
                    Pee[6*i+2]=beginPee[6*i+2]+param.radius*(cos(param.startAngle)-cos(param.startAngle+0.5*param.direction*param.beta))/2*(1-cos(param.count*PI/param.totalCount));
                }
                else
                {
                    Pee[6*i+3]=beginPee[6*i+3]+param.radius*(sin(param.startAngle)-sin(param.startAngle+0.5*param.direction*param.beta))/2*(1-cos(param.count*PI/param.totalCount));
                    Pee[6*i+4]=beginPee[6*i+4]+0.05*(1-cos(param.count*2*PI/param.totalCount));
                    Pee[6*i+5]=beginPee[6*i+5]+param.radius*(cos(param.startAngle)-cos(param.startAngle+0.5*param.direction*param.beta))/2*(1-cos(param.count*PI/param.totalCount));
                }
            }
        }
        else if(param.count<(2*param.n-1)*param.totalCount)//const
        {
            Peb[0]=param.radius*(sin(bodyStartAngle)-sin(bodyStartAngle+0.5*param.direction*param.beta*period_count/param.totalCount));
            Peb[2]=param.radius*(cos(bodyStartAngle)-cos(bodyStartAngle+0.5*param.direction*param.beta*period_count/param.totalCount));

            for(int i=0;i<3;i++)
            {
                if(half_period%2==1)
                {
                    if(param.direction==-1)
                    {
                        Pee[6*i+3]=beginPee[6*i+3]+param.radius*(sin(legStartAngle1)-sin(legStartAngle1+param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                        Pee[6*i+4]=beginPee[6*i+4]+0.05*(1-cos(period_count*2*PI/param.totalCount));
                        Pee[6*i+5]=beginPee[6*i+5]+param.radius*(cos(legStartAngle1)-cos(legStartAngle1+param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                    }
                    else
                    {
                        Pee[6*i]=beginPee[6*i]+param.radius*(sin(legStartAngle1)-sin(legStartAngle1+param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                        Pee[6*i+1]=beginPee[6*i+1]+0.05*(1-cos(period_count*2*PI/param.totalCount));
                        Pee[6*i+2]=beginPee[6*i+2]+param.radius*(cos(legStartAngle1)-cos(legStartAngle1+param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                    }
                }
                else
                {
                    if(param.direction==-1)
                    {
                        Pee[6*i]=beginPee[6*i]+param.radius*(sin(legStartAngle0)-sin(legStartAngle0+param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                        Pee[6*i+1]=beginPee[6*i+1]+0.05*(1-cos(period_count*2*PI/param.totalCount));
                        Pee[6*i+2]=beginPee[6*i+2]+param.radius*(cos(legStartAngle0)-cos(legStartAngle0+param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                    }
                    else
                    {
                        Pee[6*i+3]=beginPee[6*i+3]+param.radius*(sin(legStartAngle0)-sin(legStartAngle0+param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                        Pee[6*i+4]=beginPee[6*i+4]+0.05*(1-cos(period_count*2*PI/param.totalCount));
                        Pee[6*i+5]=beginPee[6*i+5]+param.radius*(cos(legStartAngle0)-cos(legStartAngle0+param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                    }
                }
            }
        }
        else//dec
        {
            Peb[0]=param.radius*(sin(bodyStartAngle)-sin(bodyStartAngle+0.5*param.direction*param.beta*period_count/param.totalCount-0.25*param.direction*param.beta*param.count*param.count/param.totalCount/param.totalCount));
            Peb[2]=param.radius*(cos(bodyStartAngle)-cos(bodyStartAngle+0.5*param.direction*param.beta*period_count/param.totalCount-0.25*param.direction*param.beta*param.count*param.count/param.totalCount/param.totalCount));

            for(int i=0;i<3;i++)
            {
                if(param.direction==-1)
                {
                    Pee[6*i+3]=beginPee[6*i+3]+param.radius*(sin(legStartAngle1)-sin(legStartAngle1+0.5*param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                    Pee[6*i+4]=beginPee[6*i+4]+0.05*(1-cos(period_count*2*PI/param.totalCount));
                    Pee[6*i+5]=beginPee[6*i+5]+param.radius*(cos(legStartAngle1)-cos(legStartAngle1+0.5*param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                }
                else
                {
                    Pee[6*i]=beginPee[6*i]+param.radius*(sin(legStartAngle1)-sin(legStartAngle1+0.5*param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                    Pee[6*i+1]=beginPee[6*i+1]+0.05*(1-cos(period_count*2*PI/param.totalCount));
                    Pee[6*i+2]=beginPee[6*i+2]+param.radius*(cos(legStartAngle1)-cos(legStartAngle1+0.5*param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                }
            }
        }
        return param.count-2*param.n*param.totalCount-1;
    }
}

namespace FastWalk
{
    void TestGetdJacOverPee()//work well in LegBase and Body Frame
	{
		timeval tpstart,tpend;
		float tused;
		gettimeofday(&tpstart,NULL);


		Robots::RobotTypeI rbt;
        rbt.loadXml("../../resource/Robot_VIII.xml");

        double f_s[18]{0};
        double df_s[18]{0};
        double ddf_s[18]{0};
        double s_t[2001]{0};
        double ds_t[2001]{0};
        double dds_t[2001]{0};

        double Jvi[3][3]{0};
        double dJvi_t[3][3]{0};
        double dJvi_t_in[3][3]{0};
        double outputJvi[2001][9]{0};
        double outputdJvi_t[2001][9]{0};
        double outputdJvi_t_in[2001][9]{0};

        double dJvi_x[3][3]{0};
        double dJvi_y[3][3]{0};
        double dJvi_z[3][3]{0};
		double stepH=0.05;
		double stepD=0.15;
        double initPeb[6]{0};
        double initVeb[6]{0};
        double initAeb[6]{0};
		double initPee[18]{ -0.3, -0.85, -0.65,
						   -0.45, -0.85, 0,
							-0.3, -0.85, 0.65,
							 0.3, -0.85, -0.65,
							0.45, -0.85, 0,
							 0.3, -0.85, 0.65 };
        double vEE[18]{0};
        double aEE[18]{0};
        double vEE_B[3]{0};

        for (int i=0;i<2001;i++)//i is the time(ms)
		{
			s_t[i]=PI*sin(PI*i/4000)-PI;//[-PI,0]rad
			ds_t[i]=PI*PI/4*cos(PI*i/4000);//rad/s
			dds_t[i]=-PI*PI*PI/16*sin(PI*i/4000);//rad/s^2

			for(int j=0;j<6;j++)
			{
				f_s[3*j]=stepD/2*cos(s_t[i])+stepD/2+initPee[3*j];
				f_s[3*j+1]=-stepH*sin(s_t[i])+initPee[3*j+1];
				f_s[3*j+2]=0+initPee[3*j+2];

				df_s[3*j]=-stepD/2*sin(s_t[i]);
				df_s[3*j+1]=-stepH*cos(s_t[i]);
				df_s[3*j+2]=0;

				vEE[3*j]=df_s[3*j]*ds_t[i];
				vEE[3*j+1]=df_s[3*j+1]*ds_t[i];
				vEE[3*j+2]=df_s[3*j+2]*ds_t[i];

				ddf_s[3*j]=-stepD/2*cos(s_t[i]);
				ddf_s[3*j+1]=stepH*sin(s_t[i]);
				ddf_s[3*j+2]=0;

                aEE[3*j]=ddf_s[3*j]*ds_t[i]*ds_t[i]+df_s[3*j]*dds_t[i];
                aEE[3*j+1]=ddf_s[3*j+1]*ds_t[i]*ds_t[i]+df_s[3*j+1]*dds_t[i];
                aEE[3*j+2]=ddf_s[3*j+2]*ds_t[i]*ds_t[i]+df_s[3*j+2]*dds_t[i];
			}

			rbt.SetPeb(initPeb);
			rbt.SetPee(f_s);

            //to GetDifJvi & GetVee, vel must be set, acc is not essential
            rbt.SetVb(initVeb);
            rbt.SetVee(vEE);
            //rbt.SetAb(initAeb);
            //rbt.SetAee(aEE);

			//for(int j=0;j<6;j++)
			//{
                //only related to position, not related to vel
                rbt.pLegs[0]->GetJvi(*Jvi,rbt.body());
                rbt.pLegs[0]->GetdJacOverPee(*dJvi_x,*dJvi_y,*dJvi_z,"B");

                rbt.pLegs[0]->GetDifJvi(*dJvi_t,rbt.body());
                rbt.pLegs[0]->GetVee(vEE_B,rbt.body());

				std::fill_n(*dJvi_t_in,9,0);

				aris::dynamic::s_daxpy(9,vEE_B[0],*dJvi_x,1,*dJvi_t_in,1);//for t
				aris::dynamic::s_daxpy(9,vEE_B[1],*dJvi_y,1,*dJvi_t_in,1);//for t
				aris::dynamic::s_daxpy(9,vEE_B[2],*dJvi_z,1,*dJvi_t_in,1);//for t

				//if(j==0)//to check is dJvi right by diff Jvi in matlab
				//{
					memcpy(*outputJvi+9*i,*Jvi,9*sizeof(double));
                    memcpy(*outputdJvi_t+9*i,*dJvi_t,9*sizeof(double));
                    memcpy(*outputdJvi_t_in+9*i,*dJvi_t_in,9*sizeof(double));
				//}
			//}
		}

        aris::dynamic::dlmwrite("./Jvi.txt",*outputJvi,2001,9);
        aris::dynamic::dlmwrite("./dJvi_t.txt",*outputdJvi_t,2001,9);
        aris::dynamic::dlmwrite("./dJvi_t_in.txt",*outputdJvi_t_in,2001,9);



		gettimeofday(&tpend,NULL);
		tused=tpend.tv_sec-tpstart.tv_sec+(double)(tpend.tv_usec-tpstart.tv_usec)/1000000;
		printf("UsedTime:%f\n",tused);
	}

    void TimeOptimalGait3()
    {
        timeval tpstart,tpend;
        float tused;
        gettimeofday(&tpstart,NULL);



        Robots::RobotTypeI rbt;
        rbt.loadXml("../../resource/Robot_VIII.xml");

        double initPeb[6] {0};
        double initVeb[6] {0};
        double initAeb[6] {0};
        double initPee[18] { -0.3, -0.85, -0.65,
                            -0.45, -0.85, 0,
                             -0.3, -0.85, 0.65,
                              0.3, -0.85, -0.65,
                             0.45, -0.85, 0,
                              0.3, -0.85, 0.65 };
        double pEB[6] {0};
        double pEE[18] {0};
        double stepH {0.05};
        double stepD {1.1};

        double s {0};
        double b_s {0};
        double db_s {0};
        double ddb_s {0};
        double f_s[3] {0};
        double df_s[3] {0};
        double ddf_s[3] {0};
        double f_s_B[3] {0};
        double df_s_B[3] {0};
        double ddf_s_B[3] {0};

        double Jvi[9] {0};
        double dJvi_x[9] {0};
        double dJvi_y[9] {0};
        double dJvi_z[9] {0};
        double dJvi[9] {0};

        double param_dsds[18] {0};
        double param_dsds1[18] {0};
        double param_dsds2[18] {0};
        double param_dds[18] {0};
        double abs_param_dds[9] {0};
        double param_fsB[18] {0};
        double param_Lmt[18] {0};//param_dsds divided by param_dds
        double param_ConstL[18] {0};
        double param_ConstH[18] {0};

        double output_dsds[1800][18] {0};
        double output_dsds1[1800][18] {0};
        double output_dsds2[1800][18] {0};
        double output_dds[1800][18] {0};
        double output_fsB[1800][18] {0};
        double output_Lmt[1800][18] {0};
        double output_ConstL[1800][18] {0};
        double output_ConstH[1800][18] {0};

        double ds_max_aLmt[1800][4] {0};
        double ds_max_vLmt[1800][4] {0};
        double ds_max[1800][4] {0};
        double output_ValueL[1800][4] {0};
        double output_ValueH[1800][4] {0};

        double vLmt {0.9};
        double aLmt {3.2};

        for (int i=0;i<1800;i++)
        {
            /**********Trajectory design & generation**********/
            s=0.1*i * PI/180;//degree to rad

            f_s[0]=0;
            f_s[1]=stepH*sin(PI/2*(1-cos(s)));
            f_s[2]=stepD/2*cos(PI/2*(1-cos(s)));
            b_s=stepD/4-stepD/2*(s/PI);

            df_s[0]=0;
            df_s[1]=stepH*cos(PI/2*(1-cos(s)))*PI/2*sin(s);
            df_s[2]=-stepD/2*sin(PI/2*(1-cos(s)))*PI/2*sin(s);
            db_s=-stepD/2/PI;

            ddf_s[0]=0;
            ddf_s[1]=-stepH*sin(PI/2*(1-cos(s)))*PI/2*sin(s)*PI/2*sin(s)+stepH*cos(PI/2*(1-cos(s)))*PI/2*cos(s);
            ddf_s[2]=-stepD/2*cos(PI/2*(1-cos(s)))*PI/2*sin(s)*PI/2*sin(s)-stepD/2*sin(PI/2*(1-cos(s)))*PI/2*cos(s);
            ddb_s=0;

            memcpy(df_s_B,df_s,3*sizeof(double));
            df_s_B[2]=df_s[2]-db_s;
            memcpy(ddf_s_B,ddf_s,3*sizeof(double));

            pEB[2]=initPeb[2]+b_s;
            for (int j=0;j<3;j++)
            {
                //swing leg
                pEE[6*j]=initPee[6*j]+f_s[0];
                pEE[6*j+1]=initPee[6*j+1]+f_s[1];
                pEE[6*j+2]=initPee[6*j+2]+f_s[2];

                //stance leg
                pEE[6*j+3]=initPee[6*j+3];
                pEE[6*j+4]=initPee[6*j+4];
                pEE[6*j+5]=initPee[6*j+5];
            }

            rbt.SetPeb(pEB);
            rbt.SetPee(pEE);

            /**********param of all legs calculation**********/
            for (int j=0;j<3;j++)
            {
                //swing leg
                rbt.pLegs[2*j]->GetJvi(Jvi,rbt.body());
                rbt.pLegs[2*j]->GetdJacOverPee(dJvi_x,dJvi_y,dJvi_z,"B");
                rbt.pLegs[2*j]->GetPee(f_s_B,rbt.body());
                memcpy(param_fsB+6*j,f_s_B,3*sizeof(double));

                std::fill_n(dJvi,9,0);
                aris::dynamic::s_daxpy(9,df_s[0],dJvi_x,1,dJvi,1);//for s
                aris::dynamic::s_daxpy(9,df_s[1],dJvi_y,1,dJvi,1);
                aris::dynamic::s_daxpy(9,df_s[2],dJvi_z,1,dJvi,1);

                std::fill_n(param_dsds1+6*j,3,0);
                aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,ddf_s_B,1,1,param_dsds1+6*j,1);
                std::fill_n(param_dsds2+6*j,3,0);
                aris::dynamic::s_dgemm(3,1,3,1,dJvi,3,df_s_B,1,1,param_dsds2+6*j,1);
                std::fill_n(param_dds+6*j,3,0);
                aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,df_s_B,1,1,param_dds+6*j,1);
                std::fill_n(param_dsds+6*j,3,0);
                for (int k=0;k<3;k++)
                {
                    param_dsds[6*j+k]=param_dsds1[6*j+k]+param_dsds2[6*j+k];
                    abs_param_dds[3*j+k]=fabs(param_dds[6*j+k]);

                    if(param_dds[6*j+k]>0)
                    {
                        param_Lmt[6*j+k]=-param_dsds[6*j+k]/param_dds[6*j+k];
                        param_ConstL[6*j+k]=-aLmt/param_dds[6*j+k];
                        param_ConstH[6*j+k]=aLmt/param_dds[6*j+k];
                    }
                    else if(param_dds[6*j+k]<0)
                    {
                        param_Lmt[6*j+k]=-param_dsds[6*j+k]/param_dds[6*j+k];
                        param_ConstL[6*j+k]=aLmt/param_dds[6*j+k];
                        param_ConstH[6*j+k]=-aLmt/param_dds[6*j+k];
                    }
                    else
                    {
                        printf("WARNING!!! param_dds equals zero!!! Swing : i=%d \n",i);
                    }
                }


                //stance leg, only vel in z direction
                rbt.pLegs[2*j+1]->GetJvi(Jvi,rbt.body());
                rbt.pLegs[2*j+1]->GetdJacOverPee(dJvi_x,dJvi_y,dJvi_z,"B");
                rbt.pLegs[2*j+1]->GetPee(f_s_B,rbt.body());
                memcpy(param_fsB+6*j+3,f_s_B,3*sizeof(double));

                std::fill_n(dJvi,9,0);
                aris::dynamic::s_daxpy(9,0    ,dJvi_x,1,dJvi,1);//for s
                aris::dynamic::s_daxpy(9,0    ,dJvi_y,1,dJvi,1);
                aris::dynamic::s_daxpy(9,-db_s,dJvi_z,1,dJvi,1);

                std::fill_n(param_dsds1+6*j+3,3,0);
                aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,ddf_s_B,1,1,param_dsds1+6*j+3,1);
                std::fill_n(param_dsds2+6*j+3,3,0);
                aris::dynamic::s_dgemm(3,1,3,1,dJvi,3,df_s_B,1,1,param_dsds2+6*j+3,1);
                std::fill_n(param_dds+6*j+3,3,0);
                aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,df_s_B,1,1,param_dds+6*j+3,1);
                std::fill_n(param_dsds+6*j+3,3,0);
                for (int k=0;k<3;k++)
                {
                    param_dsds[6*j+3+k]=param_dsds1[6*j+3+k]+param_dsds2[6*j+3+k];

                    if(param_dds[6*j+3+k]>0)
                    {
                        param_Lmt[6*j+3+k]=-param_dsds[6*j+3+k]/param_dds[6*j+3+k];
                        param_ConstL[6*j+3+k]=-aLmt/param_dds[6*j+3+k];
                        param_ConstH[6*j+3+k]=aLmt/param_dds[6*j+3+k];
                    }
                    else if(param_dds[6*j+3+k]<0)
                    {
                        param_Lmt[6*j+3+k]=-param_dsds[6*j+3+k]/param_dds[6*j+3+k];
                        param_ConstL[6*j+3+k]=aLmt/param_dds[6*j+3+k];
                        param_ConstH[6*j+3+k]=-aLmt/param_dds[6*j+3+k];
                    }
                    else
                    {
                        printf("WARNING!!! param_dds equals zero!!! Stance : i=%d \n",i);
                    }
                }
            }

            /**********maxds of 3 swing legs calculation**********/
            //3 swing leg, together
            int kw {0};
            bool stopFlag=false;
            while (stopFlag==false && kw<5000)
            {
                double ds=0.001*kw;
                double value_FuncL[9]{0};
                double value_FuncH[9]{0};
                double max_ValueL {0};
                double min_ValueH {0};
                for (int j=0;j<3;j++)
                {
                    for (int k=0;k<3;k++)
                    {
                        value_FuncL[3*j+k]=param_Lmt[6*j+k]*ds*ds+param_ConstL[6*j+k];
                        value_FuncH[3*j+k]=param_Lmt[6*j+k]*ds*ds+param_ConstH[6*j+k];
                    }
                }
                max_ValueL=*std::max_element(value_FuncL,value_FuncL+9);
                min_ValueH=*std::min_element(value_FuncH,value_FuncH+9);

                kw++;

                if(min_ValueH<max_ValueL)
                {
                    stopFlag=true;
                    if(i%18==0)
                    {
                        //printf("minValueH=%.4f,maxValueL=%.4f,kw=%d\n",min_ValueH,max_ValueL,kw);
                    }
                    if(kw==5000 || kw==1)
                    {
                        printf("WARNING!!! Error with itration count!!!");
                    }
                }
                else
                {
                    output_ValueL[i][0]=max_ValueL;
                    output_ValueH[i][0]=min_ValueH;
                    ds_max_aLmt[i][0]=ds;
                }
            }

            ds_max_vLmt[i][0]=vLmt/(*std::max_element(abs_param_dds,abs_param_dds+9));
            ds_max[i][0]=std::min(ds_max_aLmt[i][0],ds_max_vLmt[i][0]);

            //3 swing leg, one by one
            for (int j=0;j<3;j++)
            {
                kw=0;
                stopFlag=false;
                while (stopFlag==false && kw<5000)
                {
                    double ds=0.002*kw;
                    double value_FuncL[9]{0};
                    double value_FuncH[9]{0};
                    double max_ValueL {0};
                    double min_ValueH {0};
                    for (int k=0;k<3;k++)
                    {
                        value_FuncL[3*j+k]=param_Lmt[3*2*j+k]*ds*ds+param_ConstL[3*2*j+k];
                        value_FuncH[3*j+k]=param_Lmt[3*2*j+k]*ds*ds+param_ConstH[3*2*j+k];
                    }
                    max_ValueL=*std::max_element(value_FuncL+3*j,value_FuncL+3*j+3);
                    min_ValueH=*std::min_element(value_FuncH+3*j,value_FuncH+3*j+3);

                    kw++;

                    if(min_ValueH<max_ValueL)
                    {
                        stopFlag=true;
                        if(kw==5000 || kw==1)
                        {
                            printf("WARNING!!! Error with itration count!!! Leg:%d\n",2*j);
                        }
                        if(i%18==0)
                        {
                            //printf("minValueH=%.4f,maxValueL=%.4f,kw=%d\n",min_ValueH,max_ValueL,kw);
                        }
                    }
                    else
                    {
                        output_ValueL[i][j+1]=max_ValueL;
                        output_ValueH[i][j+1]=min_ValueH;
                        ds_max_aLmt[i][j+1]=ds;
                    }
                }

                ds_max_vLmt[i][j+1]=vLmt/(*std::max_element(abs_param_dds+3*j,abs_param_dds+3*j+3));
                ds_max[i][j+1]=std::min(ds_max_aLmt[i][j+1],ds_max_vLmt[i][j+1]);
            }

            memcpy(*output_dsds1+18*i,param_dsds1,18*sizeof(double));
            memcpy(*output_dsds2+18*i,param_dsds2,18*sizeof(double));
            memcpy(*output_dsds+18*i,param_dsds,18*sizeof(double));
            memcpy(*output_dds+18*i,param_dds,18*sizeof(double));
            memcpy(*output_fsB+18*i,param_fsB,18*sizeof(double));
            memcpy(*output_Lmt+18*i,param_Lmt,18*sizeof(double));
            memcpy(*output_ConstL+18*i,param_ConstL,18*sizeof(double));
            memcpy(*output_ConstH+18*i,param_ConstH,18*sizeof(double));
        }

        /**********Iteration to calculate ds**********/
        double ds_forward[1800] {ds_max_aLmt[0][0]};
        double ds_backward[1800] {0};
        ds_backward[1799]=ds_max_aLmt[1799][0];
        double dds_forward[1800] {0};
        double dds_backward[1800] {0};
        double delta_s {PI/1800};
        int ki_back {1799};
        int stop_back {0};
        int ki_for {0};
        bool stop_Iter {false};
        bool switch_Flag {true};//true acc, false dec
        int dec_start {0};
        int dec_end {0};
        double real_ds[1800] {0};
        double real_dds[1800] {0};
        double real_ddsMax[1800] {0};
        double real_ddsMin[1800] {0};

        double min_dist[1800];
        std::fill_n(min_dist,1800,1);

        //backward
        while (stop_Iter==false && ki_back>=0)
        {
            double dec[9] {0};
            for (int j=0;j<3;j++)
            {
                for (int k=0;k<3;k++)
                {
                    dec[3*j+k]=output_Lmt[ki_back][3*2*j+k]*ds_backward[ki_back]*ds_backward[ki_back]+output_ConstL[ki_back][3*2*j+k];
                }
            }
            dds_backward[ki_back]=*std::max_element(dec,dec+9);
            ds_backward[ki_back-1]=ds_backward[ki_back]-dds_backward[ki_back]*delta_s/ds_backward[ki_back];

            if (ds_backward[ki_back-1]>ds_max_aLmt[ki_back-1][0])
            {
                stop_Iter=true;
                stop_back=ki_back;
                printf("Backward Iteration ends at k=%d, ds_backward:%.4f\n",ki_back,ds_backward[ki_back]);
            }
            else
            {
                ki_back--;
            }
        }

        //forward
        unsigned int ki {0};
        stop_Iter=false;
        while (stop_Iter==false)
        {
            if (switch_Flag==true)
            {
                double acc[9] {0};
                for (int j=0;j<3;j++)
                {
                    for (int k=0;k<3;k++)
                    {
                        acc[3*j+k]=output_Lmt[ki_for][3*2*j+k]*ds_forward[ki_for]*ds_forward[ki_for]+output_ConstH[ki_for][3*2*j+k];
                    }
                }
                dds_forward[ki_for]=*std::min_element(acc,acc+9);
                ds_forward[ki_for+1]=ds_forward[ki_for]+dds_forward[ki_for]*delta_s/ds_forward[ki_for];

                if (ds_forward[ki_for+1]>ds_max_aLmt[ki_for+1][0])
                {
                    switch_Flag=false;
                    dec_start=ki_for;
                    printf("acc touching at k=%d, ds_forward=%.4f; ",ki_for,ds_forward[ki_for]);
                }
                else
                {
                    ki_for++;
                }
            }
            else
            {
                double dec[9] {0};
                for (int j=0;j<3;j++)
                {
                    for (int k=0;k<3;k++)
                    {
                        dec[3*j+k]=output_Lmt[ki_for][3*2*j+k]*ds_forward[ki_for]*ds_forward[ki_for]+output_ConstL[ki_for][3*2*j+k];
                    }
                }
                dds_forward[ki_for]=*std::max_element(dec,dec+9);
                ds_forward[ki_for+1]=ds_forward[ki_for]+dds_forward[ki_for]*delta_s/ds_forward[ki_for];

                if (ds_forward[ki_for+1]>ds_max_aLmt[ki_for+1][0])
                {
                    //printf("dec trying\n");
                    if(dec_start>0)
                    {
                        dec_start--;
                        ki_for=dec_start;
                    }
                    else
                    {
                        printf("dec_start=%d\n",dec_start);
                        ki_for=dec_start;
                        ds_forward[0]-=ds_max_aLmt[0][0]/1000;
                    }
                }
                else
                {
                    //printf("dec ending,ki_for=%d, ds_forward=%.4f\n",ki_for,ds_forward[ki_for+1]);
                    if (ds_forward[ki_for+1]<1)//min_maxds
                    {
                        switch_Flag=true;
                        for(int k=dec_start;k<(ki_for+2);k++)
                        {
                            min_dist[k]=ds_max_aLmt[k][0]-ds_forward[k];
                        }
                        dec_end=std::min_element(min_dist+dec_start+1,min_dist+ki_for+2)-min_dist;
                        //dec_start must be ignored, if dec_start is the min_dist, the calculation will cycle between dec_start & dec_start+1
                        ki_for=dec_end-1;
                        printf("dec finished, start at k=%d, end at k=%d, ds_forward=%.4f\n",dec_start,dec_end,ds_forward[ki_for+1]);
                    }
                    ki_for++;
                }
            }

            if(ki_for==1799)
            {
                stop_Iter=true;
                memcpy(real_ds,ds_forward,(ki_for+1)*sizeof(double));
                memcpy(real_dds,dds_forward,(ki_for+1)*sizeof(double));
                printf("forward reach the end, and never encounter with the backward, ki=%u < %u\n",ki,0xFFFFFFFF);
            }
            if(ki_for>=stop_back && ds_forward[ki_for]>=ds_backward[ki_for])
            {
                stop_Iter=true;
                memcpy(real_ds,ds_forward,ki_for*sizeof(double));
                memcpy(real_ds+ki_for,ds_backward+ki_for,(1800-ki_for)*sizeof(double));
                memcpy(real_dds,dds_forward,ki_for*sizeof(double));
                memcpy(real_dds+ki_for,dds_backward+ki_for,(1800-ki_for)*sizeof(double));
                printf("forward & backward encounters at k=%d\n",ki_for);
            }
            ki++;
            if(ki==0xFFFFFFFF)
            {
                stop_Iter=true;
                printf("WARNING!!! Iteration takes too long, force stop!!! ki=%d\n",ki);
            }

        }

        for (int i=0;i<1800;i++)
        {
            double acc[9] {0};
            double dec[9] {0};
            for (int j=0;j<3;j++)
            {
                for (int k=0;k<3;k++)
                {
                    acc[3*j+k]=output_Lmt[i][3*2*j+k]*real_ds[i]*real_ds[i]+output_ConstH[i][3*2*j+k];
                    dec[3*j+k]=output_Lmt[i][3*2*j+k]*real_ds[i]*real_ds[i]+output_ConstL[i][3*2*j+k];
                }
            }
            real_ddsMax[i]=*std::min_element(acc,acc+9);
            real_ddsMin[i]=*std::max_element(dec,dec+9);
        }

        /**********test forward curve, output some more info**********/
        //dec from k=215,k=218,k=230
        double output_ds_k[1800] {0};
        double dds_forward_test[1800] {0};
        double ds_forward_test[1800] {0};

        int ki_test {1203};
        ds_forward_test[ki_test]=real_ds[ki_test];
        output_ds_k[ki_test]=ds_forward_test[ki_test];
        double acc_test[9] {0};
        for (int j=0;j<3;j++)
        {
            for (int k=0;k<3;k++)
            {
                acc_test[3*j+k]=output_Lmt[ki_test][3*2*j+k]*ds_forward_test[ki_test]*ds_forward_test[ki_test]+output_ConstH[ki_test][3*2*j+k];
            }
        }
        dds_forward_test[ki_test]=*std::min_element(acc_test,acc_test+9);
        ds_forward_test[ki_test+1]=ds_forward_test[ki_test]+dds_forward_test[ki_test]*delta_s/ds_forward_test[ki_test];
        ki_test++;
        output_ds_k[ki_test]=ds_forward_test[ki_test];
        while(output_ds_k[ki_test]<ds_max[ki_test][0])
        {
            double dec_test[9] {0};
            for (int j=0;j<3;j++)
            {
                for (int k=0;k<3;k++)
                {
                    dec_test[3*j+k]=output_Lmt[ki_test][3*2*j+k]*ds_forward_test[ki_test]*ds_forward_test[ki_test]+output_ConstL[ki_test][3*2*j+k];
                }
            }
            dds_forward_test[ki_test]=*std::max_element(dec_test,dec_test+9);
            ds_forward_test[ki_test+1]=ds_forward_test[ki_test]+dds_forward_test[ki_test]*delta_s/ds_forward_test[ki_test];

            output_ds_k[ki_test+1]=ds_forward_test[ki_test+1];
            ki_test++;
        }

        aris::dynamic::dlmwrite("./test_ds_k.txt",output_ds_k,1800,1);


        /**********calculate final trajectory**********/
        //for s
        double vEE[18] {0};
        double aEE[18] {0};
        double output_Pee[1800][9] {0};
        double output_Pin[1800][9] {0};
        double output_Vin[1800][9] {0};
        double output_Ain[1800][9] {0};
        for (int i=0;i<1800;i++)
        {
            s=0.1*i * PI/180;//degree to rad

            f_s[0]=0;
            f_s[1]=stepH*sin(PI/2*(1-cos(s)));
            f_s[2]=stepD/2*cos(PI/2*(1-cos(s)));
            b_s=stepD/4-stepD/2*(s/PI);

            df_s[0]=0;
            df_s[1]=stepH*cos(PI/2*(1-cos(s)))*PI/2*sin(s);
            df_s[2]=-stepD/2*sin(PI/2*(1-cos(s)))*PI/2*sin(s);
            db_s=-stepD/2/PI;

            ddf_s[0]=0;
            ddf_s[1]=-stepH*sin(PI/2*(1-cos(s)))*PI/2*sin(s)*PI/2*sin(s)+stepH*cos(PI/2*(1-cos(s)))*PI/2*cos(s);
            ddf_s[2]=-stepD/2*cos(PI/2*(1-cos(s)))*PI/2*sin(s)*PI/2*sin(s)-stepD/2*sin(PI/2*(1-cos(s)))*PI/2*cos(s);
            ddb_s=0;

            memcpy(df_s_B,df_s,3*sizeof(double));
            df_s_B[2]=df_s[2]-db_s;
            memcpy(ddf_s_B,ddf_s,3*sizeof(double));

            rbt.SetPeb(initPeb);
            rbt.SetVb(initVeb);
            rbt.SetAb(initAeb);
            for (int j=0;j<3;j++)
            {
                //swing leg
                pEE[6*j]=initPee[6*j]+f_s[0];
                pEE[6*j+1]=initPee[6*j+1]+f_s[1];
                pEE[6*j+2]=initPee[6*j+2]+f_s[2]-b_s;

                vEE[6*j]=df_s[0]*real_ds[i];
                vEE[6*j+1]=df_s[1]*real_ds[i];
                vEE[6*j+2]=(df_s[2]-db_s)*real_ds[i];

                aEE[6*j]=ddf_s[0]*real_ds[i]*real_ds[i]+df_s[0]*real_dds[i];
                aEE[6*j+1]=ddf_s[1]*real_ds[i]*real_ds[i]+df_s[1]*real_dds[i];
                aEE[6*j+2]=(ddf_s[2]-ddb_s)*real_ds[i]*real_ds[i]+(df_s[2]-db_s)*real_dds[i];

                rbt.pLegs[2*j]->SetPee(pEE+6*j,rbt.body());
                rbt.pLegs[2*j]->SetVee(vEE+6*j,rbt.body());
                rbt.pLegs[2*j]->SetAee(aEE+6*j,rbt.body());

                rbt.pLegs[2*j]->GetPee(*output_Pee+9*i+3*j,rbt.body());
                rbt.pLegs[2*j]->GetPin(*output_Pin+9*i+3*j);
                rbt.pLegs[2*j]->GetVin(*output_Vin+9*i+3*j);
                rbt.pLegs[2*j]->GetAin(*output_Ain+9*i+3*j);
            }
        }

        //fot t
        double totalTime {0};
        int totalCount {0};
        double v0 {0};
        double vm {0};
        double vt {stepD/2/PI*real_ds[0]};
        double stance_begin_s;

        for (int i=0;i<1800;i++)
        {
            totalTime+=delta_s/real_ds[i];
        }
        totalCount=(int)(totalTime*1000)+1;
        printf("totalTime is %.4f, totalCount is %d\n",totalTime,totalCount);

        double * real_s=new double [totalCount];
        double * real_Pee=new double [9*2*totalCount];
        double * real_Pin=new double [9*2*totalCount];
        real_s[0]=0;
        for (int i=1;i<totalCount;i++)
        {
            double ds=0.5*(real_ds[(int)(real_s[i-1]/(PI/1800))]+real_ds[(int)(real_s[i-1]/(PI/1800))+1]);
            real_s[i]=real_s[i-1]+ds*0.001;
            if (i==totalCount-1)
            {
                double dds=0.5*(real_dds[(int)(real_s[i-1]/(PI/1800))]+real_dds[(int)(real_s[i-1]/(PI/1800))+1]);
                v0=stepD/2/PI*(ds+dds*0.001);
                stance_begin_s=real_s[totalCount-1]+ds*0.001;
            }
        }
        vm=stepD/(0.001*totalCount)-(v0+vt)/2;

        for (int i=0;i<2*totalCount;i++)
        {
            //swing phase
            if(i<totalCount)
            {
                f_s[0]=0;
                f_s[1]=stepH*sin(PI/2*(1-cos(real_s[i])));
                f_s[2]=stepD/2*cos(PI/2*(1-cos(real_s[i])))-(stepD/4-stepD/2*(real_s[i]/PI));
            }
            //stance phase
            else
            {
                f_s[0]=0;
                f_s[1]=0;
                if((i-totalCount)<(double)totalCount/2)
                {
                    f_s[2]=stepD/2*cos(PI/2*(1-cos(stance_begin_s)))-(stepD/4-stepD/2*(stance_begin_s/PI))
                            +v0*(0.001*(i-totalCount))+0.5*(vm-v0)/(0.001*totalCount/2)*0.001*(i-totalCount)*0.001*(i-totalCount);
                }
                else
                {
                    f_s[2]=stepD/2*cos(PI/2*(1-cos(stance_begin_s)))-(stepD/4-stepD/2*(stance_begin_s/PI))
                            +v0*(0.001*totalCount/2)+0.5*(vm-v0)/(0.001*totalCount/2)*0.001*totalCount/2*0.001*totalCount/2
                            +vm*(0.001*(i-1.5*totalCount))+0.5*(vt-vm)/(0.001*totalCount/2)*0.001*(i-1.5*totalCount)*0.001*(i-1.5*totalCount);
                }
            }


            rbt.SetPeb(initPeb);

            for (int j=0;j<3;j++)
            {
                //swing leg
                real_Pee[9*i+3*j]=initPee[6*j]+f_s[0];
                real_Pee[9*i+3*j+1]=initPee[6*j+1]+f_s[1];
                real_Pee[9*i+3*j+2]=initPee[6*j+2]+f_s[2];

                rbt.pLegs[2*j]->SetPee(real_Pee+9*i+3*j,rbt.pLegs[2*j]->body());
                rbt.pLegs[2*j]->GetPin(real_Pin+9*i+3*j);
            }
        }

        aris::dynamic::dlmwrite("./real_s.txt",real_s,totalCount,1);
        aris::dynamic::dlmwrite("./real_Pee.txt",real_Pee,2*totalCount,9);
        aris::dynamic::dlmwrite("./real_Pin.txt",real_Pin,2*totalCount,9);
        delete [] real_s;
        delete [] real_Pee;
        delete [] real_Pin;

        aris::dynamic::dlmwrite("./param_dsds1.txt",*output_dsds1,1800,18);
        aris::dynamic::dlmwrite("./param_dsds2.txt",*output_dsds2,1800,18);
        aris::dynamic::dlmwrite("./param_dsds.txt",*output_dsds,1800,18);
        aris::dynamic::dlmwrite("./param_dds.txt",*output_dds,1800,18);
        aris::dynamic::dlmwrite("./param_fsB.txt",*output_fsB,1800,18);
        aris::dynamic::dlmwrite("./param_Lmt.txt",*output_Lmt,1800,18);
        aris::dynamic::dlmwrite("./param_ConstL.txt",*output_ConstL,1800,18);
        aris::dynamic::dlmwrite("./param_ConstH.txt",*output_ConstH,1800,18);
        aris::dynamic::dlmwrite("./max_ValueL.txt",*output_ValueL,1800,4);
        aris::dynamic::dlmwrite("./min_ValueH.txt",*output_ValueH,1800,4);
        aris::dynamic::dlmwrite("./ds_max_aLmt.txt",*ds_max_aLmt,1800,4);
        aris::dynamic::dlmwrite("./ds_max_vLmt.txt",*ds_max_vLmt,1800,4);
        aris::dynamic::dlmwrite("./ds_max.txt",*ds_max,1800,4);

        aris::dynamic::dlmwrite("./dds_max.txt",real_ddsMax,1800,1);
        aris::dynamic::dlmwrite("./dds_min.txt",real_ddsMin,1800,1);
        aris::dynamic::dlmwrite("./ds_forward.txt",ds_forward,1800,1);
        aris::dynamic::dlmwrite("./ds_backward.txt",ds_backward,1800,1);
        aris::dynamic::dlmwrite("./dds_forward.txt",dds_forward,1800,1);
        aris::dynamic::dlmwrite("./dds_backward.txt",dds_backward,1800,1);

        aris::dynamic::dlmwrite("./Pee.txt",*output_Pee,1800,9);
        aris::dynamic::dlmwrite("./Pin.txt",*output_Pin,1800,9);
        aris::dynamic::dlmwrite("./Vin.txt",*output_Vin,1800,9);
        aris::dynamic::dlmwrite("./Ain.txt",*output_Ain,1800,9);


        gettimeofday(&tpend,NULL);
        tused=tpend.tv_sec-tpstart.tv_sec+(double)(tpend.tv_usec-tpstart.tv_usec)/1000000;
        printf("UsedTime:%f\n",tused);
    }

    void TimeOptimalGait1by1()
    {
        timeval tpstart,tpend;
        float tused;
        gettimeofday(&tpstart,NULL);


        Robots::RobotTypeI rbt;
        rbt.loadXml("../../resource/Robot_VIII.xml");

        double stepH {0.1};
        double stepD {0.8};
        double vLmt {0.9};
        double aLmt {3.2};
        const int sTotalCount {1801};
        const double delta_s {PI/(sTotalCount-1)};

        double initPeb[6] {0};
        double initPee[18] {  -0.3, -0.85, -0.65,
                             -0.45, -0.85,     0,
                              -0.3, -0.85,  0.65,
                               0.3, -0.85, -0.65,
                              0.45, -0.85,     0,
                               0.3, -0.85,  0.65 };//024 swing, 135 stance
        double pEB[6] {0};
        double pEE[18] {0};

        double s_t {0};
        double b_st[sTotalCount] {0};
        double db_st[sTotalCount] {0};
        double ddb_st[sTotalCount] {0};
        double pb_sw[sTotalCount] {0};
        double vb_sw[sTotalCount] {0};
        double ab_sw[sTotalCount] {0};

        double s_w {0};
        double f_sw[3] {0};
        double df_sw[3] {0};
        double ddf_sw[3] {0};

        double Jvi[9] {0};
        double dJvi_x[9] {0};
        double dJvi_y[9] {0};
        double dJvi_z[9] {0};
        double dJvi[9] {0};//used in stanceLeg
        double dJvi_dot_f[9] {0};//used in swingLeg
        double dJvi_dot_vb[9] {0};//used in swingLeg

        double param_dsds[18] {0};
        double param_dsds1[18] {0};
        double param_dsds2[18] {0};
        double param_dds[18] {0};
        double abs_param_dds[9] {0};//for vLmt of stanceLeg
        double param_ds[18] {0};//for aLmt of swingLeg
        double param_ds1[18] {0};
        double param_ds2[18] {0};
        double param_const[18] {0};//for aLmt of swingLeg
        double param_const1[18] {0};
        double param_const2[18] {0};

        double param_a2[18] {0};//coefficient of ds*ds, param_dsds divided by param_dds
        double param_a1[18] {0};//coefficient of ds, param_ds divided by param_dds
        double param_a0L[18] {0};
        double param_a0H[18] {0};

        double output_PeeB[sTotalCount][18] {0};
        double output_dsds[sTotalCount][18] {0};
        double output_dsds1[sTotalCount][18] {0};
        double output_dsds2[sTotalCount][18] {0};
        double output_ds[sTotalCount][18] {0};
        double output_ds1[sTotalCount][18] {0};
        double output_ds2[sTotalCount][18] {0};
        double output_const[sTotalCount][18] {0};
        double output_const1[sTotalCount][18] {0};
        double output_const2[sTotalCount][18] {0};
        double output_dds[sTotalCount][18] {0};
        double output_a2[sTotalCount][18] {0};
        double output_a1[sTotalCount][18] {0};
        double output_a0L[sTotalCount][18] {0};
        double output_a0H[sTotalCount][18] {0};

        double ds_upBound_aLmt[sTotalCount][4] {0};
        double ds_lowBound_aLmt[sTotalCount][4] {0};
        double ds_upBound_vLmt[sTotalCount][4] {0};
        double ds_upBound[sTotalCount][4] {0};
        double ds_lowBound[sTotalCount][4] {0};
        double dds_lowBound[sTotalCount][4] {0};
        double dds_upBound[sTotalCount][4] {0};

        double slopedsBound[sTotalCount][4] {0};
        double slopeDelta[sTotalCount][4] {0};
        double paramdds0Point[sTotalCount][18] {0};
        double tangentPoint[sTotalCount][4] {0};
        double switchPoint[sTotalCount][4] {0};
        int paramdds0Count[18] {0};
        int tangentCount[4] {0};
        int switchCount[4] {0};

        bool stopFlag {false};
        bool accFlag {true};//true acc, false dec
        int dec_start {0};
        int dec_end {0};
        double min_dist[sTotalCount];
        unsigned int cycleCount {0};//used to limit the cycleCount of forward integration
        int stop_back {0};
        int ki_back {0};
        int ki_for {0};

        double ds_forward[sTotalCount][4] {0};
        double ds_backward[sTotalCount][4] {0};
        double dds_forward[sTotalCount][4] {0};
        double dds_backward[sTotalCount][4] {0};
        double timeArray[sTotalCount][4] {0};
        double timeArray_tmp[sTotalCount][4] {0};
        double totalTime[4] {0};
        double maxTime {0};
        int maxTotalCount {0};
        int maxTotalCount_last {1};
        int maxTimeID {0};

        double real_ds[sTotalCount][4] {0};
        double real_dds[sTotalCount][4] {0};
        double real_ddsMax[sTotalCount][4] {0};
        double real_ddsMin[sTotalCount][4] {0};

        /*******************************************************************************/
        /********** StanceLeg : generate the traj & calculate the bound of ds **********/
        for (int i=0;i<sTotalCount;i++)
        {
            s_t=i*PI/(sTotalCount-1);//degree to rad
            b_st[i]=stepD/4-stepD/2*(s_t/PI);
            db_st[i]=-stepD/2/PI;
            ddb_st[i]=0;

            pEB[2]=initPeb[2]+b_st[i];
            rbt.SetPeb(pEB);
            for(int j=0;j<3;j++)//pEE of stance leg
            {
                pEE[6*j+3]=initPee[6*j+3];
                pEE[6*j+4]=initPee[6*j+4];
                pEE[6*j+5]=initPee[6*j+5];
                rbt.pLegs[2*j+1]->SetPee(pEE+6*j+3);
            }

            //calculate param of stance leg, only vel in z direction
            for(int j=0;j<3;j++)
            {
                rbt.pLegs[2*j+1]->GetJvi(Jvi,rbt.body());
                rbt.pLegs[2*j+1]->GetdJacOverPee(dJvi_x,dJvi_y,dJvi_z,"B");
                rbt.pLegs[2*j+1]->GetPee(*output_PeeB+18*i+6*j+3,rbt.body());

                double db_st_tmp[3] {0};
                db_st_tmp[2]=-db_st[i];
                double ddb_st_tmp[3] {0};
                ddb_st_tmp[2]=-ddb_st[i];

                std::fill_n(dJvi,9,0);
                aris::dynamic::s_daxpy(9,db_st_tmp[0],dJvi_x,1,dJvi,1);//for s_t
                aris::dynamic::s_daxpy(9,db_st_tmp[1],dJvi_y,1,dJvi,1);
                aris::dynamic::s_daxpy(9,db_st_tmp[2],dJvi_z,1,dJvi,1);

                std::fill_n(param_dds+6*j+3,3,0);
                aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,db_st_tmp,1,1,param_dds+6*j+3,1);

                std::fill_n(param_dsds1+6*j+3,3,0);
                aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,ddb_st_tmp,1,1,param_dsds1+6*j+3,1);
                std::fill_n(param_dsds2+6*j+3,3,0);
                aris::dynamic::s_dgemm(3,1,3,1,dJvi,3,db_st_tmp,1,1,param_dsds2+6*j+3,1);
                std::fill_n(param_dsds+6*j+3,3,0);
                for (int k=0;k<3;k++)
                {
                    param_dsds[6*j+3+k]=param_dsds1[6*j+3+k]+param_dsds2[6*j+3+k];
                    abs_param_dds[3*j+k]=fabs(param_dds[6*j+3+k]);

                    if(param_dds[6*j+3+k]>0)
                    {
                        param_a2[6*j+3+k]=-param_dsds[6*j+3+k]/param_dds[6*j+3+k];
                        param_a0L[6*j+3+k]=(-param_const[6*j+k+3]-aLmt)/param_dds[6*j+3+k];
                        param_a0H[6*j+3+k]=(-param_const[6*j+k+3]+aLmt)/param_dds[6*j+3+k];
                    }
                    else if(param_dds[6*j+3+k]<0)
                    {
                        param_a2[6*j+3+k]=-param_dsds[6*j+3+k]/param_dds[6*j+3+k];
                        param_a0L[6*j+3+k]=(-param_const[6*j+k+3]+aLmt)/param_dds[6*j+3+k];
                        param_a0H[6*j+3+k]=(-param_const[6*j+k+3]-aLmt)/param_dds[6*j+3+k];
                    }
                    else
                    {
                        printf("WARNING!!! param_dds equals zero!!! StanceLeg : i=%d \n",i);
                    }
                }
            }

            //calculate bound of ds of 3 stance legs
            int k_st {0};
            bool dsBoundFlag_st {false};
            const int kstCount {15000};
            while (dsBoundFlag_st==false)
            {
                double ds=0.001*k_st;
                double dec[9]{0};
                double acc[9]{0};
                double max_dec {0};
                double min_acc {0};
                for (int j=0;j<3;j++)
                {
                    for (int k=0;k<3;k++)
                    {
                        dec[3*j+k]=param_a2[6*j+3+k]*ds*ds+param_a0L[6*j+3+k];
                        acc[3*j+k]=param_a2[6*j+3+k]*ds*ds+param_a0H[6*j+3+k];
                    }
                }

                max_dec=*std::max_element(dec,dec+9);
                min_acc=*std::min_element(acc,acc+9);

                k_st++;
                if(k_st==kstCount)
                {
                    dsBoundFlag_st=true;
                    printf("WARNING!!! kstCount=%d is too small!!!\n",kstCount);
                }
                else
                {
                    if(min_acc<max_dec)
                    {
                        for(int k_st2=0;k_st2<10000;k_st2++)
                        {
                            ds=ds_upBound_aLmt[i][0]+0.001*0.0001*k_st2;
                            for (int j=0;j<3;j++)
                            {
                                for (int k=0;k<3;k++)
                                {
                                    dec[3*j+k]=param_a2[6*j+3+k]*ds*ds+param_a0L[6*j+3+k];
                                    acc[3*j+k]=param_a2[6*j+3+k]*ds*ds+param_a0H[6*j+3+k];
                                }
                            }
                            max_dec=*std::max_element(dec,dec+9);
                            min_acc=*std::min_element(acc,acc+9);

                            if(min_acc<max_dec)
                            {
                                //printf("stance ds_upBound_aLmt=%.6f\t",ds_upBound_aLmt[i][0]);
                                dsBoundFlag_st=true;
                                break;
                            }
                            else
                            {
                                dds_lowBound[i][0]=max_dec;
                                dds_upBound[i][0]=min_acc;
                                ds_upBound_aLmt[i][0]=ds;
                            }
                        }
                    }
                    else
                    {
                        ds_upBound_aLmt[i][0]=ds;
                    }
                }
            }

            ds_upBound_vLmt[i][0]=vLmt/(*std::max_element(abs_param_dds,abs_param_dds+9));
            ds_upBound[i][0]=std::min(ds_upBound_aLmt[i][0],ds_upBound_vLmt[i][0]);
            ds_lowBound[i][0]=ds_lowBound_aLmt[i][0];

            for (int j=0;j<3;j++)
            {
                memcpy(*output_dsds1+18*i+6*j+3, param_dsds1+6*j+3, 3*sizeof(double));
                memcpy(*output_dsds2+18*i+6*j+3, param_dsds2+6*j+3, 3*sizeof(double));
                memcpy( *output_dsds+18*i+6*j+3,  param_dsds+6*j+3, 3*sizeof(double));
                memcpy(  *output_dds+18*i+6*j+3,   param_dds+6*j+3, 3*sizeof(double));
                memcpy(   *output_a2+18*i+6*j+3,    param_a2+6*j+3, 3*sizeof(double));
                memcpy(  *output_a0L+18*i+6*j+3,   param_a0L+6*j+3, 3*sizeof(double));
                memcpy(  *output_a0H+18*i+6*j+3,   param_a0H+6*j+3, 3*sizeof(double));
            }
        }

        /********** StanceLeg : find switch point - param_dds=0 point & tangent point **********/
        //initialize
        for(int i=0;i<sTotalCount;i++)
        {
            tangentPoint[i][0]=-1;
            switchPoint[i][0]=-1;
            for(int j=0;j<3;j++)
            {
                for(int k=0;k<3;k++)
                {
                    paramdds0Point[i][6*j+k+3]=-1;
                }
            }
        }

        //calculate the slope of ds_upBound
        slopedsBound[0][0]=(ds_upBound[1][0]-ds_upBound[0][0])/delta_s;
        for(int i=1;i<sTotalCount-1;i++)
        {
            slopedsBound[i][0]=(ds_upBound[i+1][0]-ds_upBound[i-1][0])/2/delta_s;
        }
        slopedsBound[sTotalCount-1][0]=(ds_upBound[sTotalCount-1][0]-ds_upBound[sTotalCount-2][0])/delta_s;

        for(int i=0;i<sTotalCount;i++)
        {
            slopeDelta[i][0]=slopedsBound[i][0]-(dds_upBound[i][0]+dds_lowBound[i][0])/2;
        }

        for(int i=0;i<sTotalCount-1;i++)
        {
            for(int j=0;j<3;j++)
            {
                for(int k=0;k<3;k++)
                {
                    if(output_dds[i+1][6*j+k+3]*output_dds[i][6*j+k+3]<0 || output_dds[i][6*j+k+3]==0)
                    {
                        paramdds0Point[paramdds0Count[6*j+k+3]][6*j+k+3]=i+fabs(output_dds[i][6*j+k+3])/(fabs(output_dds[i][6*j+k+3])+fabs(output_dds[i+1][6*j+k+3]));
                        paramdds0Count[6*j+k+3]++;
                    }
                }
            }

            if(slopeDelta[i+1][0]*slopeDelta[i][0]<0 || slopeDelta[i][0]==0)
            {
                tangentPoint[tangentCount[0]][0]=i+fabs(slopeDelta[i][0])/(fabs(slopeDelta[i][0])+fabs(slopeDelta[i+1][0]));
                tangentCount[0]++;
            }
        }

        //merge tangentPoint & paramdds0Point into switchPoint
        for(int i=0;i<tangentCount[0];i++)
        {
            switchPoint[i][0]=tangentPoint[i][0];
        }
        switchCount[0]=tangentCount[0];
        for(int j=0;j<3;j++)
        {
            for(int k=0;k<3;k++)
            {
                for(int i=0;i<paramdds0Count[6*j+k+3];i++)
                {
                    switchPoint[i+switchCount[0]][0]=paramdds0Point[i][6*j+k+3];
                }
                switchCount[0]+=paramdds0Count[6*j+k+3];
            }
        }

        //filtering the same point & sorting by the value
        for(int i=0;i<switchCount[0];i++)
        {
            for(int j=i+1;j<switchCount[0];j++)
            {
                if(switchPoint[j][0]<switchPoint[i][0])
                {
                    auto tmp=switchPoint[i][0];
                    switchPoint[i][0]=switchPoint[j][0];
                    switchPoint[j][0]=tmp;
                }
                else if(switchPoint[j][0]==switchPoint[i][0])
                {
                    switchPoint[j][0]=switchPoint[switchCount[0]-1][0];
                    switchPoint[switchCount[0]-1][0]=-1;
                    j--;
                    switchCount[0]--;
                }
            }
        }

        printf("\nStanceLeg Switch Point:");
        for(int i=0;i<switchCount[0]+1;i++)
        {
            printf("%.4f,",switchPoint[i][0]);
        }
        printf("\n");


        /********** StanceLeg : numerical integration to calculate ds **********/
        /*
        for(int i=0;i<sTotalCount;i++)
        {
            real_ds[i][0]=ds_upBound[i][0];
            real_dds[i][0]=dds_upBound[i][0];
        }
        for(int m=0;m<switchCount[0];m++)
        {
            int k_st {(int)switchPoint[m][0]};
            stopFlag=false;
            ds_backward[k_st][0]=ds_upBound[k_st][0];
            while(stopFlag==false)
            {
                double dec[9] {0};
                for (int j=0;j<3;j++)
                {
                    for (int k=0;k<3;k++)
                    {
                        dec[3*j+k]=output_a2[k_st][6*j+3+k]*ds_backward[k_st][0]*ds_backward[k_st][0]+output_a0L[k_st][6*j+3+k];
                    }
                }
                dds_backward[k_st][0]=*std::max_element(dec,dec+9);
                ds_backward[k_st-1][0]=sqrt(ds_backward[k_st][0]*ds_backward[k_st][0]-2*dds_backward[k_st][0]*delta_s);

                if(ds_backward[k_st-1][0]>ds_upBound[k_st-1][0])
                {
                    stopFlag=true;
                    printf("StanceLeg backward touching upBound at %d, quit switchPoint %.4f\n",k_st-1,switchPoint[m][0]);
                }
                else if(k_st==1)
                {
                    for (int j=0;j<3;j++)
                    {
                        for (int k=0;k<3;k++)
                        {
                            dec[3*j+k]=output_a2[k_st-1][6*j+3+k]*ds_backward[k_st-1][0]*ds_backward[k_st-1][0]+output_a0L[k_st-1][6*j+3+k];
                        }
                    }
                    dds_backward[k_st-1][0]=*std::max_element(dec,dec+9);
                    for(int i=k_st-1;i<switchPoint[m][0]+1;i++)
                    {
                        real_ds[i][0]=ds_backward[i][0];
                        real_dds[i][0]=dds_backward[i][0];
                    }
                    stopFlag=true;
                    printf("StanceLeg backward touching 0, from switchPoint %.4f\n",switchPoint[m][0]);
                }
                else if(ds_backward[k_st-1][0]>=real_ds[k_st-1][0])
                {
                    for (int j=0;j<3;j++)
                    {
                        for (int k=0;k<3;k++)
                        {
                            dec[3*j+k]=output_a2[k_st-1][6*j+3+k]*ds_backward[k_st-1][0]*ds_backward[k_st-1][0]+output_a0L[k_st-1][6*j+3+k];
                        }
                    }
                    dds_backward[k_st-1][0]=*std::max_element(dec,dec+9);
                    for(int i=k_st-1;i<switchPoint[m][0]+1;i++)
                    {
                        real_ds[i][0]=ds_backward[i][0];
                        real_dds[i][0]=dds_backward[i][0];
                    }
                    stopFlag=true;
                    printf("StanceLeg backward touching last curve at %d, from switchPoint %.4f\n",k_st-1,switchPoint[m][0]);
                }
                else
                {
                    k_st--;
                }
            }
            if(ds_backward[k_st-1][0]>ds_upBound[k_st-1][0])
            {
                continue;
            }

            bool isEqual0Point {false};
            for(int j=0;j<3;j++)
            {
                for(int k=0;k<3;k++)
                {
                    for(int i=0;i<paramdds0Count[6*j+k+3];i++)
                    {
                        if(output_dds[(int)switchPoint[i][0]][6*j+k+3]==0 || slopeDelta[(int)switchPoint[i][0]][0]==0)
                        {
                            isEqual0Point=true;
                        }
                    }
                }
            }
            if(isEqual0Point==false)
            {
                k_st=switchPoint[m][0]+1;
            }
            else
            {
                k_st=switchPoint[m][0];
            }
            stopFlag=false;
            ds_forward[k_st][0]=ds_upBound[k_st][0];
            while(stopFlag==false)
            {
                double acc[9] {0};
                for (int j=0;j<3;j++)
                {
                    for (int k=0;k<3;k++)
                    {
                        acc[3*j+k]=output_a2[k_st][6*j+3+k]*ds_forward[k_st][0]*ds_forward[k_st][0]+output_a0H[k_st][6*j+3+k];
                    }
                }
                dds_forward[k_st][0]=*std::min_element(acc,acc+9);
                ds_forward[k_st+1][0]=sqrt(ds_forward[k_st][0]*ds_forward[k_st][0]+2*dds_forward[k_st][0]*delta_s);

                if(ds_forward[k_st+1][0]>ds_upBound[k_st+1][0] || k_st==sTotalCount-2)
                {
                    for(int i=switchPoint[m][0]+1;i<k_st+1;i++)
                    {
                        real_ds[i][0]=ds_forward[i][0];
                        real_dds[i][0]=dds_forward[i][0];
                    }
                    stopFlag=true;
                    printf("StanceLeg forward touching upBound at %d, from switchPoint %.4f\n",k_st,switchPoint[m][0]);
                }
                else
                {
                    k_st++;
                }
            }
        }*/
        //backward integration
        stopFlag=false;
        ki_back=sTotalCount-1;
        ds_backward[ki_back][0]=ds_upBound[ki_back][0];
        while (stopFlag==false && ki_back>=0)
        {
            double dec[9] {0};
            for (int j=0;j<3;j++)
            {
                for (int k=0;k<3;k++)
                {
                    dec[3*j+k]=output_a2[ki_back][6*j+3+k]*ds_backward[ki_back][0]*ds_backward[ki_back][0]
                              +output_a0L[ki_back][6*j+3+k];
                }
            }
            dds_backward[ki_back][0]=*std::max_element(dec,dec+9);
            ds_backward[ki_back-1][0]=sqrt(ds_backward[ki_back][0]*ds_backward[ki_back][0]-2*dds_backward[ki_back][0]*delta_s);

            if (ds_backward[ki_back-1][0]>ds_upBound[ki_back-1][0])
            {
                stopFlag=true;
                stop_back=ki_back;
                printf("StanceLeg Backward Integration ends at k=%d, ds_backward:%.4f; ds_backward[ki_back-1][0]=%.4f > ds_upBound[ki_back-1][0]=%.4f\n",
                       ki_back,ds_backward[ki_back][0],ds_backward[ki_back-1][0],ds_upBound[ki_back-1][0]);
            }
            else
            {
                ki_back--;
            }
        }

        //forward integration
        stopFlag=false;
        accFlag=true;
        ki_for=0;
        cycleCount=0;
        ds_forward[ki_for][0]=ds_upBound[ki_for][0];
        std::fill_n(min_dist,sTotalCount,1);
        while (stopFlag==false)
        {
            if (accFlag==true)
            {
                double acc[9] {0};
                for (int j=0;j<3;j++)
                {
                    for (int k=0;k<3;k++)
                    {
                        acc[3*j+k]=output_a2[ki_for][6*j+3+k]*ds_forward[ki_for][0]*ds_forward[ki_for][0]
                                  +output_a0H[ki_for][6*j+3+k];
                    }
                }
                dds_forward[ki_for][0]=*std::min_element(acc,acc+9);
                ds_forward[ki_for+1][0]=sqrt(ds_forward[ki_for][0]*ds_forward[ki_for][0]+2*dds_forward[ki_for][0]*delta_s);
                //printf("acc,sqrt:%.4f\n",ds_forward[ki_for][0]*ds_forward[ki_for][0]+2*dds_forward[ki_for][0]*delta_s);

                if (ds_forward[ki_for+1][0]>ds_upBound[ki_for+1][0])
                {
                    accFlag=false;
                    dec_start=ki_for;
                    printf("StanceLeg acc reach bound at k=%d, ds_forward=%.4f; ",ki_for,ds_forward[ki_for][0]);
                }
                else
                {
                    ki_for++;
                }
            }
            else
            {
                double dec[9] {0};
                for (int j=0;j<3;j++)
                {
                    for (int k=0;k<3;k++)
                    {
                        dec[3*j+k]=output_a2[ki_for][6*j+3+k]*ds_forward[ki_for][0]*ds_forward[ki_for][0]
                                  +output_a0L[ki_for][6*j+3+k];
                    }
                }
                dds_forward[ki_for][0]=*std::max_element(dec,dec+9);
                ds_forward[ki_for+1][0]=sqrt(ds_forward[ki_for][0]*ds_forward[ki_for][0]+2*dds_forward[ki_for][0]*delta_s);

                if (ds_forward[ki_for+1][0]>ds_upBound[ki_for+1][0])
                {
                    //printf("dec trying\n");
                    if(dec_start>0)
                    {
                        dec_start--;
                        ki_for=dec_start;
                    }
                    else
                    {
                        printf("dec_start=%d\n",dec_start);
                        ki_for=dec_start;
                        ds_forward[0][0]-=ds_upBound[0][0]/1000;
                    }
                }
                else
                {
                    if ((ds_forward[ki_for][0]*ds_forward[ki_for][0]+2*dds_forward[ki_for][0]*delta_s)<=(ds_lowBound[ki_for+1][0]*ds_lowBound[ki_for+1][0]) || ki_for==(sTotalCount-2))
                    {

                        accFlag=true;
                        for(int k=dec_start;k<(ki_for+2);k++)
                        {
                            min_dist[k]=ds_upBound[k][0]-ds_forward[k][0];
                        }
                        dec_end=std::min_element(min_dist+dec_start+1,min_dist+ki_for+2)-min_dist;
                        //dec_start must be ignored, if dec_start is the min_dist, the calculation will cycle between dec_start & dec_start+1
                        ki_for=dec_end-1;
                        printf("dec finished, start at k=%d, end at k=%d, ds_forward=%.4f\n",dec_start,dec_end,ds_forward[ki_for+1][0]);
                    }
                    ki_for++;
                }
            }

            //loop end condition
            if(ki_for==sTotalCount-1)
            {
                stopFlag=true;
                for(int i=0;i<sTotalCount;i++)
                {
                    real_ds[i][0]=ds_forward[i][0];
                    real_dds[i][0]=dds_forward[i][0];
                    //printf("i=%d, ds_forward:%.4f, dds_forward:%.4f\n",i,ds_forward[i][0],dds_forward[i][0]);
                }
                printf("StanceLeg forward reach the end, and never encounter with the backward, cycleCount=%u < %u\n",cycleCount,0x0FFFFFFF);
            }
            if(ki_for>=stop_back && ds_forward[ki_for][0]>=ds_backward[ki_for][0])
            {
                stopFlag=true;
                for(int i=0;i<ki_for;i++)
                {
                    real_ds[i][0]=ds_forward[i][0];
                    real_dds[i][0]=dds_forward[i][0];
                }
                for(int i=ki_for;i<sTotalCount;i++)
                {
                    real_ds[i][0]=ds_backward[i][0];
                    real_dds[i][0]=dds_backward[i][0];
                }
                printf("StanceLeg forward & backward encounters at k=%d\n",ki_for);
            }
            cycleCount++;
            if(cycleCount==0x0FFFFFFF)
            {
                stopFlag=true;
                printf("WARNING!!! StanceLeg integration takes too long, force stop!!! ki=%d\n",cycleCount);
            }
        }


        //pb_sw, vb_sw, ab_sw initialized here, and need to be updated during iteration
        totalTime[0]=0;
        for (int i=0;i<sTotalCount;i++)
        {
            if(i!=0)
            {
                totalTime[0]+=2*delta_s/(real_ds[i-1][0]+real_ds[i][0]);
                timeArray[i][0]=totalTime[0];
            }
            pb_sw[i]=b_st[i];
            vb_sw[i]=db_st[i]*real_ds[i][0];
            ab_sw[i]=ddb_st[i]*real_ds[i][0]*real_ds[i][0]+db_st[i]*real_dds[i][0];
        }


        /******************************* Iteration ************************************/
        int iterCount {0};
        while(maxTotalCount!=maxTotalCount_last)
        {
            iterCount++;
            maxTotalCount_last=maxTotalCount;
            /********** SwingLeg : generate the traj & calculate the bound of ds **********/
            for (int i=0;i<sTotalCount;i++)
            {
                s_w=i*PI/(sTotalCount-1);
                f_sw[0]=0;
                f_sw[1]=stepH*sin(s_w);
                f_sw[2]=stepD/2*cos(s_w);
                df_sw[0]=0;
                df_sw[1]=stepH*cos(s_w);
                df_sw[2]=-stepD/2*sin(s_w);
                ddf_sw[0]=0;
                ddf_sw[1]=-stepH*sin(s_w);
                ddf_sw[2]=-stepD/2*cos(s_w);

                pEB[2]=initPeb[2]+pb_sw[i];
                rbt.SetPeb(pEB);
                for (int j=0;j<3;j++)//pEE of swing leg
                {
                    pEE[6*j]=initPee[6*j]+f_sw[0];
                    pEE[6*j+1]=initPee[6*j+1]+f_sw[1];
                    pEE[6*j+2]=initPee[6*j+2]+f_sw[2];
                    rbt.pLegs[2*j]->SetPee(pEE+6*j);
                }

                //calculate param of swing leg
                for (int j=0;j<3;j++)
                {
                    rbt.pLegs[2*j]->GetJvi(Jvi,rbt.body());
                    rbt.pLegs[2*j]->GetdJacOverPee(dJvi_x,dJvi_y,dJvi_z,"B");
                    rbt.pLegs[2*j]->GetPee(*output_PeeB+18*i+6*j,rbt.body());

                    double vb_sw_tmp[3] {0};
                    vb_sw_tmp[2]=vb_sw[i];
                    double ab_sw_tmp[3] {0};
                    ab_sw_tmp[2]=ab_sw[i];

                    std::fill_n(dJvi_dot_f,9,0);
                    aris::dynamic::s_daxpy(9,df_sw[0],dJvi_x,1,dJvi_dot_f,1);
                    aris::dynamic::s_daxpy(9,df_sw[1],dJvi_y,1,dJvi_dot_f,1);
                    aris::dynamic::s_daxpy(9,df_sw[2],dJvi_z,1,dJvi_dot_f,1);

                    std::fill_n(dJvi_dot_vb,9,0);
                    aris::dynamic::s_daxpy(9,vb_sw_tmp[0],dJvi_x,1,dJvi_dot_vb,1);
                    aris::dynamic::s_daxpy(9,vb_sw_tmp[1],dJvi_y,1,dJvi_dot_vb,1);
                    aris::dynamic::s_daxpy(9,vb_sw_tmp[2],dJvi_z,1,dJvi_dot_vb,1);

                    std::fill_n(param_dds+6*j,3,0);
                    aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,df_sw,1,1,param_dds+6*j,1);

                    std::fill_n(param_dsds1+6*j,3,0);
                    aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,ddf_sw,1,1,param_dsds1+6*j,1);
                    std::fill_n(param_dsds2+6*j,3,0);
                    aris::dynamic::s_dgemm(3,1,3,1,dJvi_dot_f,3,df_sw,1,1,param_dsds2+6*j,1);

                    std::fill_n(param_ds1+6*j,3,0);
                    aris::dynamic::s_dgemm(3,1,3,1,dJvi_dot_vb,3,df_sw,1,1,param_ds1+6*j,1);
                    std::fill_n(param_ds2+6*j,3,0);
                    aris::dynamic::s_dgemm(3,1,3,1,dJvi_dot_f,3,vb_sw_tmp,1,1,param_ds2+6*j,1);

                    std::fill_n(param_const1+6*j,3,0);
                    aris::dynamic::s_dgemm(3,1,3,1,dJvi_dot_vb,3,vb_sw_tmp,1,1,param_const1+6*j,1);
                    std::fill_n(param_const2+6*j,3,0);
                    aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,ab_sw_tmp,1,1,param_const2+6*j,1);

                    std::fill_n(param_dsds+6*j,3,0);
                    std::fill_n(param_ds+6*j,3,0);
                    for (int k=0;k<3;k++)
                    {
                        param_dsds[6*j+k]=param_dsds1[6*j+k]+param_dsds2[6*j+k];
                        param_ds[6*j+k]=-param_ds1[6*j+k]-param_ds2[6*j+k];
                        param_const[6*j+k]=param_const1[6*j+k]-param_const2[6*j+k];

                        if(param_dds[6*j+k]>0)
                        {
                            param_a2[6*j+k]=-param_dsds[6*j+k]/param_dds[6*j+k];
                            param_a1[6*j+k]=-param_ds[6*j+k]/param_dds[6*j+k];
                            param_a0L[6*j+k]=(-aLmt-param_const[6*j+k])/param_dds[6*j+k];
                            param_a0H[6*j+k]=(aLmt-param_const[6*j+k])/param_dds[6*j+k];
                        }
                        else if(param_dds[6*j+k]<0)
                        {
                            param_a2[6*j+k]=-param_dsds[6*j+k]/param_dds[6*j+k];
                            param_a1[6*j+k]=-param_ds[6*j+k]/param_dds[6*j+k];
                            param_a0L[6*j+k]=(aLmt-param_const[6*j+k])/param_dds[6*j+k];
                            param_a0H[6*j+k]=(-aLmt-param_const[6*j+k])/param_dds[6*j+k];
                        }
                        else
                        {
                            printf("WARNING!!! param_dds equals zero!!! SwingLeg : i=%d \n",i);
                        }
                    }
                }

                /********** calculate ds bound of swingLegs 1 by 1 **********/
                for (int j=0;j<3;j++)
                {
                    bool ds_lowBoundFlag_sw {false};
                    bool ds_upBoundFlag_sw {false};
                    int k_sw {0};
                    const int kswCount {30000};
                    while (ds_upBoundFlag_sw==false)
                    {
                        double ds=0.001*k_sw;
                        double dec[3] {0};
                        double acc[3] {0};
                        double max_dec {0};
                        double min_acc {0};
                        for (int k=0;k<3;k++)
                        {
                            dec[k]=param_a2[6*j+k]*ds*ds+param_a1[6*j+k]*ds+param_a0L[6*j+k];
                            acc[k]=param_a2[6*j+k]*ds*ds+param_a1[6*j+k]*ds+param_a0H[6*j+k];
                        }
                        max_dec=*std::max_element(dec,dec+3);
                        min_acc=*std::min_element(acc,acc+3);

                        k_sw++;
                        if(k_sw==kswCount)
                        {
                            ds_upBoundFlag_sw=true;
                            printf("WARNING!!! kswCount=%d is too small!!! Leg:%d, i=%d\n",kswCount,2*j,i);
                        }
                        else
                        {
                            if(ds_lowBoundFlag_sw==false && ds_upBoundFlag_sw==false && min_acc>max_dec)
                            {
                                ds_lowBound_aLmt[i][j+1]=ds;

                                if(k_sw==1)
                                {
                                    //printf("swing ds_lowBound_aLmt=%.6f\n",ds_lowBound_aLmt[i][j+1]);
                                    ds_lowBoundFlag_sw=true;
                                }
                                else
                                {
                                    for(int k_sw2=0;k_sw2<10000;k_sw2++)
                                    {
                                        ds=ds_lowBound_aLmt[i][j+1]-0.001*0.0001*k_sw2;
                                        for (int k=0;k<3;k++)
                                        {
                                            dec[k]=param_a2[6*j+k]*ds*ds+param_a1[6*j+k]*ds+param_a0L[6*j+k];
                                            acc[k]=param_a2[6*j+k]*ds*ds+param_a1[6*j+k]*ds+param_a0H[6*j+k];
                                        }
                                        max_dec=*std::max_element(dec,dec+3);
                                        min_acc=*std::min_element(acc,acc+3);

                                        if(min_acc<max_dec)
                                        {
                                            //printf("swing ds_lowBound_aLmt=%.6f\n",ds_lowBound_aLmt[i][j+1]);
                                            ds_lowBoundFlag_sw=true;
                                            break;
                                        }
                                        else
                                        {
                                            ds_lowBound_aLmt[i][j+1]=ds;
                                        }
                                    }
                                }
                            }
                            else if(ds_lowBoundFlag_sw==true && ds_upBoundFlag_sw==false && min_acc>=max_dec)
                            {
                                dds_lowBound[i][j+1]=max_dec;
                                dds_upBound[i][j+1]=min_acc;
                                ds_upBound_aLmt[i][j+1]=ds;
                            }
                            else if(ds_lowBoundFlag_sw==true && ds_upBoundFlag_sw==false && min_acc<max_dec)
                            {
                                for(int k_sw2=0;k_sw2<10000;k_sw2++)
                                {
                                    ds=ds_upBound_aLmt[i][j+1]+0.001*0.0001*k_sw2;
                                    for (int k=0;k<3;k++)
                                    {
                                        dec[k]=param_a2[6*j+k]*ds*ds+param_a1[6*j+k]*ds+param_a0L[6*j+k];
                                        acc[k]=param_a2[6*j+k]*ds*ds+param_a1[6*j+k]*ds+param_a0H[6*j+k];
                                    }
                                    max_dec=*std::max_element(dec,dec+3);
                                    min_acc=*std::min_element(acc,acc+3);

                                    if(min_acc<max_dec)
                                    {
                                        //printf("swing ds_upBound_aLmt=%.6f,i=%d,j=%d\n",ds_upBound_aLmt[i][j+1],i,j);
                                        ds_upBoundFlag_sw=true;
                                        break;
                                    }
                                    else
                                    {
                                        dds_lowBound[i][j+1]=max_dec;
                                        dds_upBound[i][j+1]=min_acc;
                                        ds_upBound_aLmt[i][j+1]=ds;
                                    }
                                }
                            }
                        }
                    }

                    double vLmt_value[3] {0};
                    double Jvi_dot_vb[3] {0};
                    double vb_sw_tmp[3] {0};
                    vb_sw_tmp[2]=vb_sw[i];
                    aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,vb_sw_tmp,1,1,Jvi_dot_vb,1);
                    for (int k=0;k<3;k++)
                    {
                        if(param_dds[6*j+k]>0)
                        {
                            vLmt_value[k]=(Jvi_dot_vb[k]+vLmt)/param_dds[6*j+k];
                        }
                        else if(param_dds[6*j+k]<0)
                        {
                            vLmt_value[k]=(Jvi_dot_vb[k]-vLmt)/param_dds[6*j+k];
                        }
                        else
                        {
                            printf("WARNING!!! param_dds equals zero!!! SwingLeg : i=%d \n",i);
                        }

                    }
                    ds_upBound_vLmt[i][j+1]=*std::max_element(vLmt_value,vLmt_value+3);
                    ds_upBound[i][j+1]=std::min(ds_upBound_aLmt[i][j+1],ds_upBound_vLmt[i][j+1]);
                    ds_lowBound[i][j+1]=ds_lowBound_aLmt[i][j+1];

                    memcpy(*output_dsds1+18*i+6*j, param_dsds1+6*j, 3*sizeof(double));
                    memcpy(*output_dsds2+18*i+6*j, param_dsds2+6*j, 3*sizeof(double));
                    memcpy(*output_dsds+18*i+6*j,  param_dsds+6*j,  3*sizeof(double));
                    memcpy(*output_ds1+18*i+6*j, param_ds1+6*j, 3*sizeof(double));
                    memcpy(*output_ds2+18*i+6*j, param_ds2+6*j, 3*sizeof(double));
                    memcpy(*output_ds+18*i+6*j,  param_ds+6*j,  3*sizeof(double));
                    memcpy(*output_const1+18*i+6*j, param_const1+6*j, 3*sizeof(double));
                    memcpy(*output_const2+18*i+6*j, param_const2+6*j, 3*sizeof(double));
                    memcpy(*output_const+18*i+6*j,  param_const+6*j,  3*sizeof(double));
                    memcpy(*output_dds+18*i+6*j, param_dds+6*j, 3*sizeof(double));
                    memcpy(*output_a2+18*i+6*j,  param_a2+6*j,  3*sizeof(double));
                    memcpy(*output_a1+18*i+6*j,  param_a1+6*j,  3*sizeof(double));
                    memcpy(*output_a0L+18*i+6*j, param_a0L+6*j, 3*sizeof(double));
                    memcpy(*output_a0H+18*i+6*j, param_a0H+6*j, 3*sizeof(double));
                }
            }

            /********** SwingLeg : find switch point - param_dds=0 point & tangent point **********/
            for(int j=0;j<3;j++)
            {
                //initialize
                for(int k=0;k<3;k++)
                {
                    paramdds0Count[6*j+k]=0;
                }
                tangentCount[j+1]=0;
                switchCount[j+1]=0;
                for(int i=0;i<sTotalCount;i++)
                {
                    tangentPoint[i][j+1]=-1;
                    switchPoint[i][j+1]=-1;
                    for(int k=0;k<3;k++)
                    {
                        paramdds0Point[i][6*j+k]=-1;
                    }
                }

                slopedsBound[0][j+1]=(ds_upBound[1][j+1]-ds_upBound[0][j+1])/delta_s;
                for(int i=1;i<sTotalCount-1;i++)
                {
                    slopedsBound[i][j+1]=(ds_upBound[i+1][j+1]-ds_upBound[i-1][j+1])/2/delta_s;
                }
                slopedsBound[sTotalCount-1][j+1]=(ds_upBound[sTotalCount-1][j+1]-ds_upBound[sTotalCount-2][j+1])/delta_s;

                for(int i=0;i<sTotalCount;i++)
                {
                    slopeDelta[i][j+1]=slopedsBound[i][j+1]-(dds_upBound[i][j+1]+dds_lowBound[i][j+1])/2;
                }

                for(int i=0;i<sTotalCount-1;i++)
                {
                    for(int k=0;k<3;k++)
                    {
                        if(output_dds[i+1][6*j+k]*output_dds[i][6*j+k]<0 || output_dds[i][6*j+k]==0)
                        {
                            paramdds0Point[paramdds0Count[6*j+k]][6*j+k]=i+fabs(output_dds[i][6*j+k])/(fabs(output_dds[i][6*j+k])+fabs(output_dds[i+1][6*j+k]));;
                            paramdds0Count[6*j+k]++;
                        }
                    }

                    if(slopeDelta[i+1][j+1]*slopeDelta[i][j+1]<0 || slopeDelta[i][j+1]==0)
                    {
                        tangentPoint[tangentCount[j+1]][j+1]=i+fabs(slopeDelta[i][j+1])/(fabs(slopeDelta[i][j+1])+fabs(slopeDelta[i+1][j+1]));
                        tangentCount[j+1]++;
                    }
                }

                //merge tangentPoint & paramdds0Point into switchPoint
                for(int i=0;i<tangentCount[j+1];i++)
                {
                    switchPoint[i][j+1]=tangentPoint[i][j+1];
                }
                switchCount[j+1]=tangentCount[j+1];
                for(int k=0;k<3;k++)
                {
                    for(int i=0;i<paramdds0Count[6*j+k];i++)
                    {
                        switchPoint[i+switchCount[j+1]][j+1]=paramdds0Point[i][6*j+k];
                    }
                    switchCount[j+1]+=paramdds0Count[6*j+k];
                }

                //filtering the same point & sorting by the value
                for(int i=0;i<switchCount[j+1];i++)
                {
                    for(int k=i+1;k<switchCount[j+1];k++)
                    {
                        if(switchPoint[k][j+1]<switchPoint[i][j+1])
                        {
                            auto tmp=switchPoint[i][j+1];
                            switchPoint[i][j+1]=switchPoint[k][j+1];
                            switchPoint[k][j+1]=tmp;
                        }
                        else if(switchPoint[k][j+1]==switchPoint[i][j+1])
                        {
                            switchPoint[k][j+1]=switchPoint[switchCount[j+1]-1][j+1];
                            switchPoint[switchCount[j+1]-1][j+1]=-1;
                            k--;
                            switchCount[j+1]--;
                        }
                    }
                }

                printf("\nSwingLeg Switch Point:");
                for(int i=0;i<switchCount[j+1]+1;i++)
                {
                    printf("%.4f,",switchPoint[i][j+1]);
                }
                printf("\n");
            }

            /********** SwingLeg : numerical integration to calculate ds**********/

            for(int j=0;j<3;j++)
            {/*
                for(int i=0;i<sTotalCount;i++)
                {
                    real_ds[i][j+1]=ds_upBound[i][j+1];
                    real_dds[i][j+1]=dds_upBound[i][j+1];
                }
                for(int m=0;m<switchCount[j+1];m++)
                {
                    int k_sw {(int)switchPoint[m][j+1]};
                    stopFlag=false;
                    ds_backward[k_sw][j+1]=ds_upBound[k_sw][j+1];
                    while(stopFlag==false)
                    {
                        double dec[3] {0};
                        for (int k=0;k<3;k++)
                        {
                            dec[k]=output_a2[k_sw][6*j+k]*ds_backward[k_sw][j+1]*ds_backward[k_sw][j+1]
                                  +output_a1[k_sw][6*j+k]*ds_backward[k_sw][j+1]
                                  +output_a0L[k_sw][6*j+k];
                        }
                        dds_backward[k_sw][j+1]=*std::max_element(dec,dec+3);
                        ds_backward[k_sw-1][j+1]=sqrt(ds_backward[k_sw][j+1]*ds_backward[k_sw][j+1]-2*dds_backward[k_sw][j+1]*delta_s);

                        if(ds_backward[k_sw-1][j+1]>ds_upBound[k_sw-1][j+1])
                        {
                            stopFlag=true;
                            printf("SwingLeg backward touching upBound at %d, quit switchPoint %.4f\n",k_sw-1,switchPoint[m][j+1]);
                        }
                        else if(k_sw==1)
                        {
                            for (int k=0;k<3;k++)
                            {
                                dec[k]=output_a2[k_sw][6*j+k]*ds_backward[k_sw][j+1]*ds_backward[k_sw][j+1]
                                      +output_a1[k_sw][6*j+k]*ds_backward[k_sw][j+1]
                                      +output_a0L[k_sw][6*j+k];
                            }
                            dds_backward[k_sw-1][j+1]=*std::max_element(dec,dec+3);
                            for(int i=k_sw-1;i<switchPoint[m][j+1]+1;i++)
                            {
                                real_ds[i][j+1]=ds_backward[i][j+1];
                                real_dds[i][j+1]=dds_backward[i][j+1];
                            }
                            stopFlag=true;
                            printf("SwingLeg backward touching 0, from switchPoint %.4f\n",switchPoint[m][j+1]);
                        }
                        else if(ds_backward[k_sw-1][j+1]>=real_ds[k_sw-1][j+1])
                        {
                            for (int k=0;k<3;k++)
                            {
                                dec[k]=output_a2[k_sw-1][6*j+k]*ds_backward[k_sw-1][j+1]*ds_backward[k_sw-1][j+1]
                                      +output_a1[k_sw-1][6*j+k]*ds_backward[k_sw-1][j+1]
                                      +output_a0L[k_sw-1][6*j+k];
                            }
                            dds_backward[k_sw-1][j+1]=*std::max_element(dec,dec+3);
                            for(int i=k_sw-1;i<switchPoint[m][j+1]+1;i++)
                            {
                                real_ds[i][j+1]=ds_backward[i][j+1];
                                real_dds[i][j+1]=dds_backward[i][j+1];
                            }
                            stopFlag=true;
                            printf("SwingLeg backward touching last curve at %d, from switchPoint %.4f\n",k_sw-1,switchPoint[m][j+1]);
                        }
                        else
                        {
                            k_sw--;
                        }
                    }
                    if(ds_backward[k_sw-1][j+1]>ds_upBound[k_sw-1][j+1])
                    {
                        continue;
                    }

                    bool isEqual0Point {false};
                    for(int k=0;k<3;k++)
                    {
                        for(int i=0;i<paramdds0Count[6*j+k];i++)
                        {
                            if(output_dds[(int)switchPoint[i][j+1]][6*j+k]==0 || slopeDelta[(int)switchPoint[i][j+1]][j+1]==0)
                            {
                                isEqual0Point=true;
                            }
                        }
                    }
                    if(isEqual0Point==false)
                    {
                        k_sw=switchPoint[m][j+1]+1;
                    }
                    else
                    {
                        k_sw=switchPoint[m][j+1];
                    }
                    stopFlag=false;
                    ds_forward[k_sw][j+1]=ds_upBound[k_sw][j+1];
                    while(stopFlag==false)
                    {
                        double acc[3] {0};
                        for (int k=0;k<3;k++)
                        {
                            acc[k]=output_a2[k_sw][6*j+k]*ds_forward[k_sw][j+1]*ds_forward[k_sw][j+1]
                                  +output_a1[k_sw][6*j+k]*ds_forward[k_sw][j+1]
                                  +output_a0H[k_sw][6*j+k];
                        }
                        dds_forward[k_sw][j+1]=*std::min_element(acc,acc+3);
                        ds_forward[k_sw+1][j+1]=sqrt(ds_forward[k_sw][j+1]*ds_forward[k_sw][j+1]+2*dds_forward[k_sw][j+1]*delta_s);

                        if(ds_forward[k_sw+1][j+1]>ds_upBound[k_sw+1][j+1] || k_sw==sTotalCount-2)
                        {
                            for(int i=switchPoint[m][j+1]+1;i<k_sw+1;i++)
                            {
                                real_ds[i][j+1]=ds_forward[i][j+1];
                                real_dds[i][j+1]=dds_forward[i][j+1];
                            }
                            stopFlag=true;
                            printf("SwingLeg forward touching upBound at %d, from switchPoint %.4f\n",k_sw,switchPoint[m][j+1]);
                        }
                        else
                        {
                            k_sw++;
                        }
                    }
                }*/
                //backward integration
                stopFlag=false;
                ki_back=sTotalCount-1;
                ds_backward[ki_back][j+1]=0;
                while (stopFlag==false && ki_back>=0)
                {
                    double dec[3] {0};
                    for (int k=0;k<3;k++)
                    {
                        dec[k]=output_a2[ki_back][6*j+k]*ds_backward[ki_back][j+1]*ds_backward[ki_back][j+1]
                              +output_a1[ki_back][6*j+k]*ds_backward[ki_back][j+1]
                              +output_a0L[ki_back][6*j+k];
                    }
                    dds_backward[ki_back][j+1]=*std::max_element(dec,dec+3);
                    ds_backward[ki_back-1][j+1]=sqrt(ds_backward[ki_back][j+1]*ds_backward[ki_back][j+1]-2*dds_backward[ki_back][j+1]*delta_s);

                    if (ds_backward[ki_back-1][j+1]>ds_upBound[ki_back-1][j+1])
                    {
                        stopFlag=true;
                        stop_back=ki_back;
                        //printf("SwingLeg Backward Iteration ends at k=%d, ds_backward:%.4f\n",ki_back,ds_backward[ki_back][j+1]);
                        //printf("ds_backward[ki_back-1]:%.4f; ds_bound[ki_back-1]:%.4f\n",ds_backward[ki_back-1][j+1],ds_upBound[ki_back-1][j+1]);
                    }
                    else
                    {
                        ki_back--;
                    }
                }

                //forward integration
                stopFlag=false;
                accFlag=true;
                ki_for=0;
                cycleCount=0;
                ds_forward[ki_for][j+1]=0;
                std::fill_n(min_dist,sTotalCount,1);
                while (stopFlag==false)
                {
                    if (accFlag==true)
                    {
                        double acc[3] {0};
                        for (int k=0;k<3;k++)
                        {
                            acc[k]=output_a2[ki_for][6*j+k]*ds_forward[ki_for][j+1]*ds_forward[ki_for][j+1]
                                  +output_a1[ki_for][6*j+k]*ds_forward[ki_for][j+1]
                                  +output_a0H[ki_for][6*j+k];
                        }
                        dds_forward[ki_for][j+1]=*std::min_element(acc,acc+3);
                        ds_forward[ki_for+1][j+1]=sqrt(ds_forward[ki_for][j+1]*ds_forward[ki_for][j+1]+2*dds_forward[ki_for][j+1]*delta_s);

                        if (ds_forward[ki_for+1][j+1]>ds_upBound[ki_for+1][j+1])
                        {
                            accFlag=false;
                            dec_start=ki_for;
                            //printf("SwingLeg acc touching at k=%d, ds_forward=%.4f; ",ki_for,ds_forward[ki_for][j+1]);
                        }
                        else
                        {
                            ki_for++;
                        }
                    }
                    else
                    {
                        double dec[3] {0};
                        for (int k=0;k<3;k++)
                        {
                            dec[k]=output_a2[ki_for][6*j+k]*ds_forward[ki_for][j+1]*ds_forward[ki_for][j+1]
                                  +output_a1[ki_for][6*j+k]*ds_forward[ki_for][j+1]
                                  +output_a0L[ki_for][6*j+k];
                        }
                        dds_forward[ki_for][j+1]=*std::max_element(dec,dec+3);
                        ds_forward[ki_for+1][j+1]=sqrt(ds_forward[ki_for][j+1]*ds_forward[ki_for][j+1]+2*dds_forward[ki_for][j+1]*delta_s);

                        if (ds_forward[ki_for+1][j+1]>ds_upBound[ki_for+1][j+1])
                        {
                            //printf("dec trying\n");
                            if(dec_start>0)
                            {
                                dec_start--;
                                ki_for=dec_start;
                            }
                            else
                            {
                                //printf("dec_start=%d, ki_for=%d, ds_forward=%.4f, ds_bound=%.4f\n",dec_start,ki_for+1,ds_forward[ki_for+1][j+1],ds_upBound[ki_for+1][j+1]);
                                ki_for=dec_start;
                                ds_forward[0][j+1]-=ds_upBound[0][j+1]/1000;
                            }
                        }
                        else
                        {
                            //printf("dec ending,ki_for=%d, ds_forward=%.4f\n",ki_for,ds_forward[ki_for+1][j+1]);
                            if ((ds_forward[ki_for][j+1]*ds_forward[ki_for][j+1]+2*dds_forward[ki_for][j+1]*delta_s)<=(ds_lowBound[ki_for+1][j+1]*ds_lowBound[ki_for+1][j+1]) || ki_for==(sTotalCount-2))
                            {
                                accFlag=true;
                                for(int k=dec_start;k<(ki_for+2);k++)
                                {
                                    min_dist[k]=ds_upBound[k][j+1]-ds_forward[k][j+1];
                                }
                                dec_end=std::min_element(min_dist+dec_start+1,min_dist+ki_for+2)-min_dist;
                                //dec_start must be ignored, if dec_start is the min_dist, the calculation will cycle between dec_start & dec_start+1
                                ki_for=dec_end-1;
                                //printf("dec finished, start at k=%d, end at k=%d, ds_forward=%.4f\n",dec_start,dec_end,ds_forward[ki_for+1][j+1]);
                            }
                            ki_for++;
                        }
                    }

                    if(ki_for==sTotalCount-1)
                    {
                        stopFlag=true;
                        for(int i=0;i<sTotalCount;i++)
                        {
                            real_ds[i][j+1]=ds_forward[i][j+1];
                            real_dds[i][j+1]=dds_forward[i][j+1];
                        }
                        //printf("SwingLeg forward reach the end, and never encounter with the backward, ki=%u < %u\n",cycleCount,0x0FFFFFFF);
                    }
                    if(ki_for>=stop_back && ds_forward[ki_for][j+1]>=ds_backward[ki_for][j+1])
                    {
                        stopFlag=true;
                        for(int i=0;i<ki_for;i++)
                        {
                            real_ds[i][j+1]=ds_forward[i][j+1];
                            real_dds[i][j+1]=dds_forward[i][j+1];
                        }
                        for(int i=ki_for;i<sTotalCount;i++)
                        {
                            real_ds[i][j+1]=ds_backward[i][j+1];
                            real_dds[i][j+1]=dds_backward[i][j+1];
                        }
                        //printf("SwingLeg forward & backward encounters at k=%d\n",ki_for);
                    }
                    cycleCount++;
                    if(cycleCount==0x0FFFFFFF)
                    {
                        stopFlag=true;
                        printf("WARNING!!! SwingLeg integration takes too long, force stop!!! ki=%d\n",cycleCount);
                    }
                }

                totalTime[j+1]=0;
                for (int i=1;i<sTotalCount;i++)
                {
                    totalTime[j+1]+=2*delta_s/(real_ds[i-1][j+1]+real_ds[i][j+1]);
                    timeArray[i][j+1]=totalTime[j+1];
                }
            }

            maxTime=*std::max_element(totalTime,totalTime+4);
            maxTotalCount=(int)(maxTime*1000)+1;
            maxTime=0.001*maxTotalCount;
            maxTimeID=std::max_element(totalTime,totalTime+4)-totalTime;
            printf("totalTime: %.4f, %.4f, %.4f, %.4f; maxTime:%.4f\n",totalTime[0],totalTime[1],totalTime[2],totalTime[3],maxTime);

            //update pb_sw, vb_sw, ab_sw here
            //double timeArray_tmp[sTotalCount][4] {0};
            for(int i=0;i<sTotalCount;i++)
            {
                memcpy(*timeArray_tmp+4*i,*timeArray+4*i,4*sizeof(double));
                for(int j=0;j<4;j++)
                {
                    timeArray_tmp[i][j]*=maxTime/totalTime[j];
                }
            }
            int j_start {0};
            double pb_sw_tmp[sTotalCount] {0};
            double vb_sw_tmp[sTotalCount] {0};
            double ab_sw_tmp[sTotalCount] {0};
            for(int i=0;i<(sTotalCount-1);i++)
            {
                //if(i%100==99)
                    //printf("j_start=%d, ",j_start);
                for(int j=j_start;j<(sTotalCount-1);j++)
                {
                    if(timeArray_tmp[i][maxTimeID]>=timeArray_tmp[j][0] && timeArray_tmp[i][maxTimeID]<timeArray_tmp[j+1][0])
                    {
                        j_start=j;
                        pb_sw_tmp[i]=(pb_sw[j+1]-pb_sw[j])/(timeArray_tmp[j+1][0]-timeArray_tmp[j][0])*(timeArray_tmp[i][maxTimeID]-timeArray_tmp[j][0])+pb_sw[j];
                        vb_sw_tmp[i]=(vb_sw[j+1]-vb_sw[j])/(timeArray_tmp[j+1][0]-timeArray_tmp[j][0])*(timeArray_tmp[i][maxTimeID]-timeArray_tmp[j][0])+vb_sw[j];
                        vb_sw_tmp[i]*=totalTime[0]/maxTime;
                        ab_sw_tmp[i]=(ab_sw[j+1]-ab_sw[j])/(timeArray_tmp[j+1][0]-timeArray_tmp[j][0])*(timeArray_tmp[i][maxTimeID]-timeArray_tmp[j][0])+ab_sw[j];
                        ab_sw_tmp[i]*=totalTime[0]/maxTime*totalTime[0]/maxTime;
                        break;
                    }
                }
                //if(i%100==99)
                    //printf("j_end=%d\n",j_start);
            }
            pb_sw_tmp[sTotalCount-1]=pb_sw[sTotalCount-1];
            vb_sw_tmp[sTotalCount-1]=vb_sw[sTotalCount-1]*totalTime[0]/maxTime;
            ab_sw_tmp[sTotalCount-1]=ab_sw[sTotalCount-1]*totalTime[0]/maxTime*totalTime[0]/maxTime;

            memcpy(pb_sw,pb_sw_tmp,sTotalCount*sizeof(double));
            memcpy(vb_sw,vb_sw_tmp,sTotalCount*sizeof(double));
            memcpy(ab_sw,ab_sw_tmp,sTotalCount*sizeof(double));
        }

        printf("Iteration finished, iterCount=%d, maxTime=%.4f, timeArray:%.4f,%.4f\n\n",iterCount,maxTime,timeArray[0][0],timeArray[1][0]);

        aris::dynamic::dlmwrite("./param_PeeB.txt",*output_PeeB,sTotalCount,18);
        aris::dynamic::dlmwrite("./param_dsds1.txt",*output_dsds1,sTotalCount,18);
        aris::dynamic::dlmwrite("./param_dsds2.txt",*output_dsds2,sTotalCount,18);
        aris::dynamic::dlmwrite("./param_dsds.txt", *output_dsds, sTotalCount,18);
        aris::dynamic::dlmwrite("./param_ds1.txt",*output_ds1,sTotalCount,18);
        aris::dynamic::dlmwrite("./param_ds2.txt",*output_ds2,sTotalCount,18);
        aris::dynamic::dlmwrite("./param_ds.txt", *output_ds, sTotalCount,18);
        aris::dynamic::dlmwrite("./param_const1.txt",*output_const1,sTotalCount,18);
        aris::dynamic::dlmwrite("./param_const2.txt",*output_const2,sTotalCount,18);
        aris::dynamic::dlmwrite("./param_const.txt", *output_const, sTotalCount,18);
        aris::dynamic::dlmwrite("./param_dds.txt",*output_dds,sTotalCount,18);

        aris::dynamic::dlmwrite("./param_a2.txt",*output_a2,sTotalCount,18);
        aris::dynamic::dlmwrite("./param_a1.txt",*output_a1,sTotalCount,18);
        aris::dynamic::dlmwrite("./param_a0L.txt",*output_a0L,sTotalCount,18);
        aris::dynamic::dlmwrite("./param_a0H.txt",*output_a0H,sTotalCount,18);

        aris::dynamic::dlmwrite("./max_ValueL.txt",*dds_lowBound,sTotalCount,4);
        aris::dynamic::dlmwrite("./min_ValueH.txt",*dds_upBound,sTotalCount,4);
        aris::dynamic::dlmwrite("./ds_upBound_aLmt.txt",*ds_upBound_aLmt,sTotalCount,4);
        aris::dynamic::dlmwrite("./ds_lowBound_aLmt.txt",*ds_lowBound_aLmt,sTotalCount,4);
        aris::dynamic::dlmwrite("./ds_bound_vLmt.txt",*ds_upBound_vLmt,sTotalCount,4);
        aris::dynamic::dlmwrite("./ds_lowBound.txt",*ds_lowBound,sTotalCount,4);
        aris::dynamic::dlmwrite("./ds_upBound.txt",*ds_upBound,sTotalCount,4);
        aris::dynamic::dlmwrite("./dds_max.txt",*real_ddsMax,sTotalCount,4);
        aris::dynamic::dlmwrite("./dds_min.txt",*real_ddsMin,sTotalCount,4);
        aris::dynamic::dlmwrite("./ds_forward.txt",*ds_forward,sTotalCount,4);
        aris::dynamic::dlmwrite("./ds_backward.txt",*ds_backward,sTotalCount,4);
        aris::dynamic::dlmwrite("./dds_forward.txt",*dds_forward,sTotalCount,4);
        aris::dynamic::dlmwrite("./dds_backward.txt",*dds_backward,sTotalCount,4);

        aris::dynamic::dlmwrite("./real_ds.txt",*real_ds,sTotalCount,4);
        aris::dynamic::dlmwrite("./real_dds.txt",*real_dds,sTotalCount,4);

        aris::dynamic::dlmwrite("./timeArray.txt",*timeArray_tmp,sTotalCount,4);
        aris::dynamic::dlmwrite("./pb_sw.txt",pb_sw,sTotalCount,1);
        aris::dynamic::dlmwrite("./vb_sw.txt",vb_sw,sTotalCount,1);
        aris::dynamic::dlmwrite("./ab_sw.txt",ab_sw,sTotalCount,1);

        //aris::dynamic::dlmwrite("./Pee.txt",*output_Pee,sTotalCount,9);
        //aris::dynamic::dlmwrite("./Pin.txt",*output_Pin,sTotalCount,9);
        //aris::dynamic::dlmwrite("./Vin.txt",*output_Vin,sTotalCount,9);
        //aris::dynamic::dlmwrite("./Ain.txt",*output_Ain,sTotalCount,9);


        gettimeofday(&tpend,NULL);
        tused=tpend.tv_sec-tpstart.tv_sec+(double)(tpend.tv_usec-tpstart.tv_usec)/1000000;
        printf("UsedTime:%f\n",tused);
    }


    double OfflineGait::pIn_acc[900][18];
    double OfflineGait::pIn_const[1800][18];
    double OfflineGait::pIn_dec[900][18];
    double OfflineGait::pIn_entire[3258][18];
    OfflineGait::OfflineGait()
	{
	}
    OfflineGait::~OfflineGait()
	{
	}
    void OfflineGait::parseFastWalkByPY(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
	{
        OfflineGaitParam param;

		aris::dynamic::dlmread("/home/hex/Desktop/mygit/Robots/src/Robot_Type_I/resource/Robot_VIII/pIn_acc.txt",*pIn_acc);
		aris::dynamic::dlmread("/home/hex/Desktop/mygit/Robots/src/Robot_Type_I/resource/Robot_VIII/pIn_const.txt",*pIn_const);
		aris::dynamic::dlmread("/home/hex/Desktop/mygit/Robots/src/Robot_Type_I/resource/Robot_VIII/pIn_dec.txt",*pIn_dec);
		for(auto &i:params)
		{
			if(i.first=="n")
			{
                param.n=std::stoi(i.second);
			}
			else
			{
				std::cout<<"parse failed"<<std::endl;
			}
		}
		msg.copyStruct(param);
		std::cout<<"finished parse"<<std::endl;
	}
    int OfflineGait::fastWalkByPY(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
	{
		auto &robot = static_cast<Robots::RobotBase &>(model);
        auto &param = static_cast<const OfflineGaitParam &>(param_in);

        if (param.count<param.totalCountPY)
		{
			robot.SetPin(*pIn_acc+18*param.count);
		}
        else if(param.count<param.totalCountPY*(2*param.n-1))
		{
            robot.SetPin(*pIn_const+18*((param.count-param.totalCountPY)%(2*param.totalCountPY)));
		}
		else
		{
            robot.SetPin(*pIn_dec+18*(param.count-param.totalCountPY*(2*param.n-1)));
		}
        return param.count-param.totalCountPY*2*param.n+1;
	}
    void OfflineGait::parseFastWalkByCZJ(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
    {
        OfflineGaitParam param;
        aris::dynamic::dlmread("./entirePin.txt",*pIn_entire);
        msg.copyStruct(param);
    }
    int OfflineGait::fastWalkByCZJ(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
    {
        auto &robot = static_cast<Robots::RobotBase &>(model);
        auto &param = static_cast<const OfflineGaitParam &>(param_in);

        robot.SetPin(*pIn_entire+18*param.count);

        return param.count-param.totalCountCZJ+1;
    }

    void FastWalkPYAnalyse()
	{
		Robots::RobotTypeI rbt;
        rbt.loadXml("../../resource/Robot_VIII.xml");

		double pIn_acc[900][19];
		double pIn_const[1800][19];
		double pIn_dec[900][19];
		double pEB[6]{0,0,0,0,0,0};

		double pEE_acc[900][18];
		double pEE_const[1800][18];
		double pEE_dec[900][18];

		aris::dynamic::dlmread("/home/hex/Desktop/mygit/Robots/src/Robot_Type_I/resource/Robot_VIII/pIn_acc.txt",*pIn_acc);
		aris::dynamic::dlmread("/home/hex/Desktop/mygit/Robots/src/Robot_Type_I/resource/Robot_VIII/pIn_const.txt",*pIn_const);
		aris::dynamic::dlmread("/home/hex/Desktop/mygit/Robots/src/Robot_Type_I/resource/Robot_VIII/pIn_dec.txt",*pIn_dec);

		for (int i=0;i<900;i++)
		{
			rbt.SetPeb(pEB);
			rbt.SetPin(*pIn_acc+18*i);
			rbt.GetPee(*pEE_acc+18*i);
		}
		for (int i=0;i<1800;i++)
		{
			rbt.SetPeb(pEB);
			rbt.SetPin(*pIn_const+18*i);
			rbt.GetPee(*pEE_const+18*i);
		}
		for (int i=0;i<900;i++)
		{
			rbt.SetPeb(pEB);
			rbt.SetPin(*pIn_dec+18*i);
			rbt.GetPee(*pEE_dec+18*i);
		}

		aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/Robots/src/Robot_Type_I/resource/Robot_VIII/pEE_acc.txt",*pEE_acc,900,18);
		aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/Robots/src/Robot_Type_I/resource/Robot_VIII/pEE_const.txt",*pEE_const,1800,18);
		aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/Robots/src/Robot_Type_I/resource/Robot_VIII/pEE_dec.txt",*pEE_dec,900,18);
	}

    void WalkPYAnalyse()//ellipse traj
	{
		Robots::RobotTypeI rbt;
        rbt.loadXml("../../resource/Robot_VIII.xml");

        printf("wk analyse begin\n");
		Robots::WalkParam param;
        param.alpha=0;
        param.beta=0;
        param.d=0.8;
        param.h=0.1;
        param.totalCount=1335;//Ain reach the limit when t=1335, d=0.8, h=0.1
        param.n=3;

        double pEE[2*param.n*param.totalCount][18];
        double pIn[2*param.n*param.totalCount][18];
        double pEE_B[2*param.n*param.totalCount][18];
        double pEB[2*param.n*param.totalCount][6];

        for(param.count=0;param.count<2*param.n*param.totalCount;param.count++)
		{
			Robots::walkGait(rbt,param);
            rbt.GetPeb(*pEB+6*param.count);
            rbt.GetPee(*pEE+18*param.count);
            rbt.GetPee(*pEE_B+18*param.count,rbt.body());
			rbt.GetPin(*pIn+18*param.count);
		}

        aris::dynamic::dlmwrite("./wk_pEB.txt",*pEB,2*param.n*param.totalCount,6);
        aris::dynamic::dlmwrite("./wk_pEE.txt",*pEE,2*param.n*param.totalCount,18);
        aris::dynamic::dlmwrite("./wk_pEE_B.txt",*pEE_B,2*param.n*param.totalCount,18);
        aris::dynamic::dlmwrite("./wk_pIn.txt",*pIn,2*param.n*param.totalCount,18);
	}


	/*find gait with maxVel in joint space by iteration*/
	void screwInterpolationTraj()
	{
		Robots::RobotTypeI rbt;
        rbt.loadXml("../../resource/Robot_VIII.xml");

		timeval tpstart,tpend;
		float tused;
		gettimeofday(&tpstart,NULL);

		double totalTime{1.5};

		double ratio{0.7};//control point, from middle to edge [0,1]

		double stepH=0.04;
		double stepD=1.1;
		double initPeb[6]{0,0,0,0,0,0};
		double initVeb[6]{0,0,0,0,0,0};
		double initPee[18]{ -0.3, -0.85, -0.65,
						   -0.45, -0.85, 0,
							-0.3, -0.85, 0.65,
							 0.3, -0.85, -0.65,
							0.45, -0.85, 0,
							 0.3, -0.85, 0.65 };
		double initVee[18]{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		double pEE[18];
		double pEE_B[18];

		double aEE[18];
		double vLmt{0.9};
		double aLmt{3.2*2};

		int keyPointNum{3};
		double keyPee_B[keyPointNum][18];
		double keyPin[keyPointNum][18];
		double keyVin[keyPointNum][18];

		rbt.SetPeb(initPeb);
		rbt.SetVb(initVeb);

		for (int i=0;i<6;i++)
		{
			keyPee_B[0][3*i]=initPee[3*i];
			keyPee_B[0][3*i+1]=initPee[3*i+1];
			keyPee_B[0][3*i+2]=initPee[3*i+2]+stepD/4;

			keyPee_B[1][3*i]=initPee[3*i];
			keyPee_B[1][3*i+1]=initPee[3*i+1]+stepH;
			if(i==0 || i==3)
			{
				keyPee_B[1][3*i+2]=initPee[3*i+2]+stepD/4*ratio;
			}
			else if (i==2 || i==5)
			{
				keyPee_B[1][3*i+2]=initPee[3*i+2]-stepD/4*ratio;
			}
			else
			{
				keyPee_B[1][3*i+2]=initPee[3*i+2];
			}

			keyPee_B[2][3*i]=initPee[3*i];
			keyPee_B[2][3*i+1]=initPee[3*i+1];
			keyPee_B[2][3*i+2]=initPee[3*i+2]-stepD/4;
		}

		double e{1};
		int k{0};
		double totalTmax{0},totalT[18];
		double t1[18],t2[18],t3[18],t4[18],t5[18],t6[18];
		double s1[18],s2[18],s3[18],s4[18],s5[18],s6[18];

		while(e>0.001)
		{
			for(int i=0;i<6;i++)
			{
				initVee[3*i+2]=stepD/2/totalTime;
			}

			for (int i=0;i<keyPointNum;i++)
			{
				rbt.SetPee(*keyPee_B+18*i);
				rbt.SetVee(initVee);

				rbt.GetPin(*keyPin+18*i);
				rbt.GetVin(*keyVin+18*i);
			}

			totalTmax=0;

			for (int i=0;i<18;i++)
			{
				s1[i]=(vLmt*vLmt-keyVin[0][i]*keyVin[0][i])/2/aLmt;
				s3[i]=vLmt*vLmt/2/aLmt;
				s2[i]=keyPin[0]-keyPin[1]-s1[i]-s3[i];
				if(s2[i]>=0)
				{
					//printf("s2_const exist, the s2 is %.4f\n",s2[i]);
					t1[i]=(vLmt+keyVin[0][i])/aLmt;//dec
					t2[i]=s2[i]/vLmt;//const
					t3[i]=vLmt/aLmt;//acc
				}
				else
				{
					//printf("s2_const dont exist,the max vel is %f\n",sqrt(aLmt*(keyPin[0][i]-keyPin[1][i])+0.5*keyVin[0][i]*keyVin[0][i]));
					t1[i]=(sqrt(aLmt*(keyPin[0][i]-keyPin[1][i])+0.5*keyVin[0][i]*keyVin[0][i])+keyVin[0][i])/aLmt;//dec
					t2[i]=0;
					t3[i]=sqrt(aLmt*(keyPin[0][i]-keyPin[1][i])+0.5*keyVin[0][i]*keyVin[0][i])/aLmt;//acc
				}

				s4[i]=vLmt*vLmt/2/aLmt;
				s6[i]=(vLmt*vLmt-keyVin[2][i]*keyVin[2][i])/2/aLmt;
				s5[i]=keyPin[2][i]-keyPin[1][i]-s4[i]-s6[i];

				if(s5[i]>=0)
				{
					//printf("s5_const exist, the s5 is %.4f\n",s5[i]);
					t4[i]=vLmt/aLmt;
					t5[i]=s5[i]/vLmt;
					t6[i]=(vLmt-keyVin[2][i])/aLmt;
				}
				else
				{
					//printf("s5_const dont exist,the max vel is %f\n",sqrt(aLmt*(keyPin[2][i]-keyPin[1][i])+0.5*keyVin[2][i]*keyVin[2][i]));
					t4[i]=sqrt(aLmt*(keyPin[2][i]-keyPin[1][i])+0.5*keyVin[2][i]*keyVin[2][i])/aLmt;
					t5[i]=0;
					t6[i]=(sqrt(aLmt*(keyPin[2][i]-keyPin[1][i])+0.5*keyVin[2][i]*keyVin[2][i])-keyVin[2][i])/aLmt;
				}

				totalT[i]=t1[i]+t2[i]+t3[i]+t4[i]+t5[i]+t6[i];
				if(totalT[i]>totalTmax)
				{
					totalTmax=totalT[i];
					k=i;
				}
			}

			e=fabs(totalTmax-totalTime);

			totalTime=totalTmax;

			printf("%.4f,screwID:%d\n",totalTmax,k);
			printf("\n");
		}

		gettimeofday(&tpend,NULL);
		tused=tpend.tv_sec-tpstart.tv_sec+(double)(tpend.tv_usec-tpstart.tv_usec)/1000000;
		printf("UsedTime:%f\n",tused);

		totalTime=1;
		int totalCount=round(totalTime*1000);
		double vIn[3000][18];
		double pIn[3000][18];
		double vIn_adjust[totalCount][18];
		double pIn_adjust[totalCount][18];
		double pEE_adjust[totalCount][18];

		for (int i=0;i<18;i++)
		{
			s1[i]=keyVin[0][i]*t1[i]-0.5*aLmt*t1[i]*t1[i];//shift in phase 1 ( vector with direction)
			s3[i]=-0.5*aLmt*t3[i]*t3[i];
			s2[i]=keyPin[1][i]-keyPin[0][i]-s1[i]-s3[i];
			s4[i]=0.5*aLmt*t4[i]*t4[i];
			s6[i]=keyVin[2][i]*t6[i]+0.5*aLmt*t6[i]*t6[i];
			s5[i]=keyPin[2][i]-keyPin[1][i]-s4[i]-s6[i];
		}

		for (int i=0;i<18;i++)
		{
			//printf("t:%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",t1[i],t2[i],t3[i],t4[i],t5[i],t6[i],totalT[i]);
			//printf("s:%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",s1[i],s2[i],s3[i],s4[i],s5[i],s6[i]);
		}

		for (int i=0;i<3000;i++)
		{
			for (int j=0;j<18;j++)
			{
				if(((double)i/1000)<t1[j])
				{
					vIn[i][j]=keyVin[0][j]-aLmt*i/1000;
					pIn[i][j]=keyPin[0][j]+keyVin[0][j]*i/1000-0.5*aLmt*i/1000*i/1000;
				}
				else if(((double)i/1000)<(t1[j]+t2[j]))
				{
					vIn[i][j]=-vLmt;
					pIn[i][j]=keyPin[0][j]+s1[j]-vLmt*((double)i/1000-t1[j]);
				}
				else if(((double)i/1000)<(t1[j]+t2[j]+t3[j]+t4[j]))
				{
					vIn[i][j]=keyVin[0][j]-aLmt*t1[j]+aLmt*((double)i/1000-t1[j]-t2[j]);
					pIn[i][j]=keyPin[0][j]+s1[j]+s2[j]+(keyVin[0][j]-aLmt*t1[j])*((double)i/1000-t1[j]-t2[j])+0.5*aLmt*((double)i/1000-t1[j]-t2[j])*((double)i/1000-t1[j]-t2[j]);
				}
				else if(((double)i/1000)<(t1[j]+t2[j]+t3[j]+t4[j]+t5[j]))
				{
					vIn[i][j]=vLmt;
					pIn[i][j]=keyPin[0][j]+s1[j]+s2[j]+s3[j]+s4[j]+vLmt*((double)i/1000-(t1[j]+t2[j]+t3[j]+t4[j]));
				}
				else if(((double)i/1000)<(t1[j]+t2[j]+t3[j]+t4[j]+t5[j]+t6[j]))
				{
					vIn[i][j]=keyVin[0][j]-aLmt*t1[j]+aLmt*(t3[j]+t4[j])-aLmt*((double)i/1000-t1[j]-t3[j]-t4[j]);
					pIn[i][j]=keyPin[0][j]+s1[j]+s2[j]+s3[j]+s4[j]+s5[j]+(keyVin[0][j]-aLmt*t1[j]+aLmt*(t3[j]+t4[j]))*((double)i/1000-t1[j]-t2[j]-t3[j]-t4[j]-t5[j])-0.5*aLmt*((double)i/1000-t1[j]-t2[j]-t3[j]-t4[j]-t5[j])*((double)i/1000-t1[j]-t2[j]-t3[j]-t4[j]-t5[j]);
				}
				else
				{
					vIn[i][j]=keyVin[2][j];
					pIn[i][j]=keyPin[2][j];
				}
			}
		}

		for (int i=0;i<totalCount;i++)
		{
			for (int j=0;j<18;j++)
			{
				if(((double)i/1000)<(t1[j]*totalTime/totalT[j]))
				{
					pIn_adjust[i][j]=keyPin[0][j]+keyVin[0][j]*((double)i/1000*totalT[j]/totalTime)-0.5*aLmt*((double)i/1000*totalT[j]/totalTime)*((double)i/1000*totalT[j]/totalTime);
				}
				else if(((double)i/1000)<((t1[j]+t2[j])*totalTime/totalT[j]))
				{
					pIn_adjust[i][j]=keyPin[0][j]+s1[j]-vLmt*((double)i/1000*totalT[j]/totalTime-t1[j]);
				}
				else if(((double)i/1000)<((t1[j]+t2[j]+t3[j]+t4[j])*totalTime/totalT[j]))
				{
					pIn_adjust[i][j]=keyPin[0][j]+s1[j]+s2[j]+(keyVin[0][j]-aLmt*t1[j])*((double)i/1000*totalT[j]/totalTime-t1[j]-t2[j])+0.5*aLmt*((double)i/1000*totalT[j]/totalTime-t1[j]-t2[j])*((double)i/1000*totalT[j]/totalTime-t1[j]-t2[j]);
				}
				else if(((double)i/1000)<((t1[j]+t2[j]+t3[j]+t4[j]+t5[j])*totalTime/totalT[j]))
				{
					pIn_adjust[i][j]=keyPin[0][j]+s1[j]+s2[j]+s3[j]+s4[j]+vLmt*((double)i/1000*totalT[j]/totalTime-t1[j]-t2[j]-t3[j]-t4[j]);
				}
				else if(((double)i/1000)<((t1[j]+t2[j]+t3[j]+t4[j]+t5[j]+t6[j])*totalTime/totalT[j]))
				{
					pIn_adjust[i][j]=keyPin[0][j]+s1[j]+s2[j]+s3[j]+s4[j]+s5[j]+(keyVin[0][j]-aLmt*t1[j]+aLmt*(t3[j]+t4[j]))*((double)i/1000*totalT[j]/totalTime-t1[j]-t2[j]-t3[j]-t4[j]-t5[j])-0.5*aLmt*((double)i/1000*totalT[j]/totalTime-t1[j]-t2[j]-t3[j]-t4[j]-t5[j])*((double)i/1000*totalT[j]/totalTime-t1[j]-t2[j]-t3[j]-t4[j]-t5[j]);
				}
				else
				{
					pIn_adjust[i][j]=keyPin[2][j];
				}
			}

			rbt.SetPin(*pIn_adjust+18*i);
			rbt.GetPee(*pEE_adjust+18*i,rbt.body());
		}

        aris::dynamic::dlmwrite("../../Server/vIn.txt",*vIn,3000,18);
        aris::dynamic::dlmwrite("../../Server/pIn.txt",*pIn,3000,18);
        aris::dynamic::dlmwrite("../../Server/pIn_adjust.txt",*pIn_adjust,totalCount,18);
        aris::dynamic::dlmwrite("../../Server/pEE_adjust.txt",*pEE_adjust,totalCount,18);

	}


	NormalGait::WalkState JointSpaceWalk::walkState;
	double JointSpaceWalk::bodyAcc;
	double JointSpaceWalk::bodyDec;
	int JointSpaceWalk::totalCount;
	double JointSpaceWalk::height;
	double JointSpaceWalk::beta;
	double JointSpaceWalk::beginPee[18];
	double JointSpaceWalk::beginVel;
	double JointSpaceWalk::endPee[18];
	double JointSpaceWalk::endVel;
	double JointSpaceWalk::distance;
	NormalGait::GaitPhase JointSpaceWalk::gaitPhase[6];//swing true, stance false
	bool JointSpaceWalk::constFlag;

	JointSpaceWalk::JointSpaceWalk()
	{
	}
	JointSpaceWalk::~JointSpaceWalk()
	{
	}

    void JointSpaceWalk::parseJointSpaceFastWalk(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
	{
		JointSpaceWalkParam param;

		for(auto &i:params)
		{
			if(i.first=="init")
			{
				walkState=NormalGait::WalkState::Init;
				msg.copyStruct(param);
			}
			else if(i.first=="acc")
			{
                bodyAcc=std::stod(i.second);
                walkState=NormalGait::WalkState::ForwardAcc;
			}
			else if(i.first=="stop")
			{
				walkState=NormalGait::WalkState::Stop;
			}
			else if(i.first=="totalCount")
			{
                totalCount=std::stoi(i.second);
			}
			else if(i.first=="height")
			{
                height=std::stod(i.second);
			}
			else if(i.first=="beta")
			{
                beta=std::stod(i.second);
			}
			else
			{
				std::cout<<"parse failed"<<std::endl;
			}
		}

		std::cout<<"finished parse"<<std::endl;
	}

	void JointSpaceWalk::swingLegTg(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in, int legID)
	{
		auto &robot = static_cast<Robots::RobotBase &>(model);
		auto &param = static_cast<const JointSpaceWalkParam &>(param_in);

		double ratio{0.7};//control point, from middle to edge [0,1]
		double stepH=height;
		double stepD=distance;
		double vLmt{0.9};
		double aLmt{3.2};

		int keyPointNum{3};
		double keyPee_B[keyPointNum][3];
		double keyVee_B[keyPointNum][3];
		std::fill_n(*keyVee_B,keyPointNum*3,0);
		double keyPin[keyPointNum][3];
		double keyVin[keyPointNum][3];

		double halfGaitTime=(double)(param.count%totalCount)/1000;
		double totalTime=(double)totalCount/1000;
		double totalT[3];
		double t1[3],t2[3],t3[3],t4[3],t5[3],t6[3];
		double s1[3],s2[3],s3[3],s4[3],s5[3],s6[3];
		double pIn_adjust[3];

		//keyPee_B
		keyPee_B[0][0]=beginPee[3*legID];
		keyPee_B[0][1]=beginPee[3*legID+1];
		keyPee_B[0][2]=beginPee[3*legID+2];

		keyPee_B[1][0]=beginPee[3*legID];
		keyPee_B[1][1]=beginPee[3*legID+1]+stepH;
		if(legID==0 || legID==3)
		{
			keyPee_B[1][2]=beginPee[3*legID+2]-stepD/2*(1-ratio);
		}
		else if (legID==2 || legID==5)
		{
			keyPee_B[1][2]=beginPee[3*legID+2]-stepD/2*(1+ratio);
		}
		else
		{
			keyPee_B[1][2]=beginPee[3*legID+2]-stepD/2;
		}

		keyPee_B[2][0]=beginPee[3*legID];
		keyPee_B[2][1]=beginPee[3*legID+1];
		keyPee_B[2][2]=beginPee[3*legID+2]-stepD;

		//keyVee_B
		keyVee_B[0][2]=beginVel;
		keyVee_B[2][2]=endVel;

		//keyPin & keyVin
		robot.pLegs[legID]->SetPee(*keyPee_B);
		robot.pLegs[legID]->SetVee(*keyVee_B);
		robot.pLegs[legID]->GetPin(*keyPin);
		robot.pLegs[legID]->GetVin(*keyVin);

		robot.pLegs[legID]->SetPee(*keyPee_B+3);
		robot.pLegs[legID]->GetPin(*keyPin+3);

		robot.pLegs[legID]->SetPee(*keyPee_B+6);
		robot.pLegs[legID]->SetVee(*keyVee_B+6);
		robot.pLegs[legID]->GetPin(*keyPin+6);
		robot.pLegs[legID]->GetVin(*keyVin+6);

		for (int i=0;i<3;i++)
		{
			//cal t1[i]~t6[i]
			s1[i]=(vLmt*vLmt-keyVin[0][i]*keyVin[0][i])/2/aLmt;
			s3[i]=vLmt*vLmt/2/aLmt;
			s2[i]=keyPin[0]-keyPin[1]-s1[i]-s3[i];
			if(s2[i]>=0)
			{
				//printf("s2_const exist, the s2 is %.4f\n",s2[i]);
				t1[i]=(vLmt+keyVin[0][i])/aLmt;//dec
				t2[i]=s2[i]/vLmt;//const
				t3[i]=vLmt/aLmt;//acc
			}
			else
			{
				//printf("s2_const dont exist,the max vel is %f\n",sqrt(aLmt*(keyPin[0][i]-keyPin[1][i])+0.5*keyVin[0][i]*keyVin[0][i]));
				t1[i]=(sqrt(aLmt*(keyPin[0][i]-keyPin[1][i])+0.5*keyVin[0][i]*keyVin[0][i])+keyVin[0][i])/aLmt;//dec
				t2[i]=0;
				t3[i]=sqrt(aLmt*(keyPin[0][i]-keyPin[1][i])+0.5*keyVin[0][i]*keyVin[0][i])/aLmt;//acc
			}

			s4[i]=vLmt*vLmt/2/aLmt;
			s6[i]=(vLmt*vLmt-keyVin[2][i]*keyVin[2][i])/2/aLmt;
			s5[i]=keyPin[2][i]-keyPin[1][i]-s4[i]-s6[i];

			if(s5[i]>=0)
			{
				//printf("s5_const exist, the s5 is %.4f\n",s5[i]);
				t4[i]=vLmt/aLmt;
				t5[i]=s5[i]/vLmt;
				t6[i]=(vLmt-keyVin[2][i])/aLmt;
			}
			else
			{
				//printf("s5_const dont exist,the max vel is %f\n",sqrt(aLmt*(keyPin[2][i]-keyPin[1][i])+0.5*keyVin[2][i]*keyVin[2][i]));
				t4[i]=sqrt(aLmt*(keyPin[2][i]-keyPin[1][i])+0.5*keyVin[2][i]*keyVin[2][i])/aLmt;
				t5[i]=0;
				t6[i]=(sqrt(aLmt*(keyPin[2][i]-keyPin[1][i])+0.5*keyVin[2][i]*keyVin[2][i])-keyVin[2][i])/aLmt;
			}

			totalT[i]=t1[i]+t2[i]+t3[i]+t4[i]+t5[i]+t6[i];

			//cal s1[i]~s6[i]
			s1[i]=keyVin[0][i]*t1[i]-0.5*aLmt*t1[i]*t1[i];//shift in phase 1 ( vector with direction)
			s3[i]=-0.5*aLmt*t3[i]*t3[i];
			s2[i]=keyPin[1][i]-keyPin[0][i]-s1[i]-s3[i];
			s4[i]=0.5*aLmt*t4[i]*t4[i];
			s6[i]=keyVin[2][i]*t6[i]+0.5*aLmt*t6[i]*t6[i];
			s5[i]=keyPin[2][i]-keyPin[1][i]-s4[i]-s6[i];

			//scaling to adjust to the totalCount
			if(halfGaitTime<(t1[i]*totalTime/totalT[i]))
			{
				pIn_adjust[i]=keyPin[0][i]+keyVin[0][i]*(halfGaitTime*totalT[i]/totalTime)-0.5*aLmt*(halfGaitTime*totalT[i]/totalTime)*(halfGaitTime*totalT[i]/totalTime);
			}
			else if(halfGaitTime<((t1[i]+t2[i])*totalTime/totalT[i]))
			{
				pIn_adjust[i]=keyPin[0][i]+s1[i]-vLmt*(halfGaitTime*totalT[i]/totalTime-t1[i]);
			}
			else if(halfGaitTime<((t1[i]+t2[i]+t3[i]+t4[i])*totalTime/totalT[i]))
			{
				pIn_adjust[i]=keyPin[0][i]+s1[i]+s2[i]+(keyVin[0][i]-aLmt*t1[i])*(halfGaitTime*totalT[i]/totalTime-t1[i]-t2[i])+0.5*aLmt*(halfGaitTime*totalT[i]/totalTime-t1[i]-t2[i])*(halfGaitTime*totalT[i]/totalTime-t1[i]-t2[i]);
			}
			else if(halfGaitTime<((t1[i]+t2[i]+t3[i]+t4[i]+t5[i])*totalTime/totalT[i]))
			{
				pIn_adjust[i]=keyPin[0][i]+s1[i]+s2[i]+s3[i]+s4[i]+vLmt*(halfGaitTime*totalT[i]/totalTime-t1[i]-t2[i]-t3[i]-t4[i]);
			}
			else if(halfGaitTime<((t1[i]+t2[i]+t3[i]+t4[i]+t5[i]+t6[i])*totalTime/totalT[i]))
			{
				pIn_adjust[i]=keyPin[0][i]+s1[i]+s2[i]+s3[i]+s4[i]+s5[i]+(keyVin[0][i]-aLmt*t1[i]+aLmt*(t3[i]+t4[i]))*(halfGaitTime*totalT[i]/totalTime-t1[i]-t2[i]-t3[i]-t4[i]-t5[i])-0.5*aLmt*(halfGaitTime*totalT[i]/totalTime-t1[i]-t2[i]-t3[i]-t4[i]-t5[i])*(halfGaitTime*totalT[i]/totalTime-t1[i]-t2[i]-t3[i]-t4[i]-t5[i]);
			}
			else
			{
				pIn_adjust[i]=keyPin[2][i];
			}
		}

		robot.pLegs[legID]->SetPin(pIn_adjust);
	}

	void JointSpaceWalk::stanceLegTg(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in, int legID)
	{
		auto &robot = static_cast<Robots::RobotBase &>(model);
		auto &param = static_cast<const JointSpaceWalkParam &>(param_in);

		double pEE_B[3];
		double halfGaitTime=(double)(param.count%totalCount)/1000;

		pEE_B[0]=beginPee[3*legID];
		pEE_B[1]=beginPee[3*legID+1];
		pEE_B[2]=beginPee[3*legID+2]+(beginVel*halfGaitTime+0.5*(endVel-beginVel)/totalCount*1000*halfGaitTime*halfGaitTime);

		robot.pLegs[legID]->SetPee(pEE_B,robot.body());
	}

	int JointSpaceWalk::jointSpaceFastWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
	{
		auto &robot = static_cast<Robots::RobotBase &>(model);
		auto &param = static_cast<const JointSpaceWalkParam &>(param_in);

		outputParam OPP;

		if(param.count%(2*totalCount)<totalCount)
		{
			for (int i=0;i<3;i++)
			{
				gaitPhase[2*i]=NormalGait::GaitPhase::Swing;//swing
				gaitPhase[2*i+1]=NormalGait::GaitPhase::Stance;//stance
			}
		}
		else
		{
			for (int i=0;i<3;i++)
			{
				gaitPhase[2*i]=NormalGait::GaitPhase::Stance;//stance
				gaitPhase[2*i+1]=NormalGait::GaitPhase::Swing;//swing
			}
		}

		if(param.count%totalCount==(totalCount-1))
		{
            if(walkState==NormalGait::WalkState::ForwardAcc && constFlag==true)
			{
				walkState=NormalGait::WalkState::Const;
				constFlag=false;
				rt_printf("Acc finished, the coming Const Vel is: %.4f\n",endVel);
			}
		}
		if(param.count%totalCount==0)
		{
			beginVel=endVel;
			memcpy(beginPee,endPee,sizeof(endPee));

			switch(walkState)
			{
			case NormalGait::WalkState::Init:
				robot.GetPee(beginPee,robot.body());
				robot.GetPee(endPee,robot.body());
				beginVel=0;
				endVel=0;
				constFlag=false;
				break;

            case NormalGait::WalkState::ForwardAcc:
				endVel=beginVel+bodyAcc*totalCount/1000;
				constFlag=true;
				break;

			case NormalGait::WalkState::Const:
				endVel=beginVel;
				break;
			}
			distance=(beginVel+endVel)/2*totalCount/1000;
			for (int i=0;i<6;i++)
			{
				endPee[3*i]=beginPee[3*i];
				endPee[3*i+1]=beginPee[3*i+1];
                if(gaitPhase[i]==NormalGait::GaitPhase::Swing)
				{
					endPee[3*i+2]=beginPee[3*i+2]-distance;
				}
				else
				{
					endPee[3*i+2]=beginPee[3*i+2]+distance;
				}
			}
		}

		double initPeb[6]{0,0,0,0,0,0};
		double initVeb[6]{0,0,0,0,0,0};
		robot.SetPeb(initPeb);
		robot.SetVb(initVeb);

		for (int i=0;i<6;i++)
		{
			if(gaitPhase[i]==NormalGait::GaitPhase::Swing)
			{
				swingLegTg(robot,param,i);
			}
			else if(gaitPhase[i]==NormalGait::GaitPhase::Stance)
			{
				stanceLegTg(robot,param,i);
			}
		}

		robot.GetPin(OPP.outputPin);
		robot.GetPee(OPP.outputPee,robot.body());
		fastWalkPipe.sendToNrt(OPP);

		if(walkState==NormalGait::WalkState::Stop && param.count%totalCount==(totalCount-1))
		{
			return 0;
		}
		else
		{
			return 1;
		}
	}


	/*analyse the maxVel of Pee in workspace*/
	void maxCal(aris::dynamic::Model &model, double alpha, int legID, double *maxVel, double *maxAcc)
	{
		auto &robot = static_cast<Robots::RobotBase &>(model);

		double vLmt[2]{-0.9,0.9};
		double aLmt[2]{-3.2,3.2};

		double direction[3]{0,cos(alpha-PI/2),sin(alpha-PI/2)};
		double Kv{0};
		double Ka{0};

		double Jvi[3][3];
		double dJvi[3][3];
		double tem[3];

		//maxVel
		robot.pLegs[legID]->GetJvi(*Jvi,robot.body());
		aris::dynamic::s_dgemm(3,1,3,1,*Jvi,3,direction,1,0,tem,1);
		for (int i=0;i<3;i++)
		{
			if(i==0)
			{
				Kv=vLmt[1]/fabs(tem[i]);
			}
			else
			{
				if(Kv>=vLmt[1]/fabs(tem[i]))
				{
					Kv=vLmt[1]/fabs(tem[i]);
				}
			}
		}
		for (int i=0;i<3;i++)
		{
			maxVel[i]=Kv*direction[i];
		}

		//printf("alpha:%f ; tem:%f,%f,%f ; Kv:%f\n",alpha,tem[0],tem[1],tem[2],Kv);

		//maxAcc
		robot.pLegs[legID]->GetDifJvi(*dJvi,robot.body());
	}

	void maxVel()
	{
		timeval tpstart,tpend;
		float tused;
		gettimeofday(&tpstart,NULL);

		Robots::RobotTypeI rbt;
        rbt.loadXml("../../resource/Robot_VIII.xml");

		double initPeb[6]{0,0,0,0,0,0};
		double initPee[18]{ -0.3, -0.85, -0.65,
						   -0.45, -0.85, 0,
							-0.3, -0.85, 0.65,
							 0.3, -0.85, -0.65,
							0.45, -0.85, 0,
							 0.3, -0.85, 0.65 };
		double pEE[18];
		int stepHNum=11;
		int stepDNum=21;
		int angleNum=36;
		double alpha;
		double maxVel[stepDNum*stepHNum*angleNum][18];
		double maxAcc[stepDNum*stepHNum*angleNum][18];

		for (int i=0;i<stepDNum;i++)
		{
			for (int j=0;j<stepHNum;j++)
			{
				for (int k=0;k<6;k++)
				{
					pEE[3*k]=initPee[3*k];
					pEE[3*k+1]=initPee[3*k+1]-0.1+0.4*j/(stepHNum-1);
					pEE[3*k+2]=initPee[3*k+2]-0.4+0.8*i/(stepDNum-1);
				}
				rbt.SetPeb(initPeb);
				rbt.SetPee(pEE);

				for(int a=0;a<angleNum;a++)
				{
					alpha=(double)a/angleNum*PI;
					for (int k=0;k<6;k++)
					{
						maxCal(rbt,alpha,k,*maxVel+((i*stepHNum+j)*angleNum+a)*18+3*k,*maxAcc+((i*stepHNum+j)*angleNum+a)*18+3*k);
					}
				}
			}
		}

        aris::dynamic::dlmwrite("../../Server/maxVel.txt",*maxVel,stepDNum*stepHNum*angleNum,18);

		gettimeofday(&tpend,NULL);
		tused=tpend.tv_sec-tpstart.tv_sec+(double)(tpend.tv_usec-tpstart.tv_usec)/1000000;
		printf("UsedTime:%f\n",tused);
	}

	//when this function called, make sure that the robot model is at the state set in last ms.
	void getMaxPin(double* maxPin,aris::dynamic::Model &model,maxVelParam &param_in)
	{
		auto &robot = static_cast<Robots::RobotBase &>(model);
		//param not casted? problems may appear

		double Jvi[9][6];
		double difJvi[9][6];
		double bodyVel_spatial[6];
		double bodyAcc_spatial[6];
		double aIn[18]{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		double vIn[18]{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		//if(param_in.legState==false)

		robot.GetJvi(*Jvi,"000111000111000111");
		robot.GetDifJvi(*difJvi,"000111000111000111");
		robot.GetVb(bodyVel_spatial);//cannot get
		robot.GetAb(bodyAcc_spatial);//cannot get
	}
}


namespace ForceTask
{
    void forceInit(int count, const double* forceRaw_in, double* forceInF_out)
	{
		static double forceSum[6];
		static double forceAvg[6]{0,0,0,0,0,0};
        //double forceInF[6]{0,0,0,0,0,0};
		if(count==0)
		{
            std::fill_n(forceSum,6,0);
		}
        if(count<200)
        {
            for(int i=0;i<6;i++)
            {
                forceSum[i]+=forceRaw_in[i];
                forceAvg[i]=forceSum[i]/(count+1);
            }
        }

		for(int i=0;i<6;i++)
		{
            forceInF_out[i]=forceRaw_in[i]-forceAvg[i];
		}
        //aris::dynamic::s_f2f(forcePm,forceInF,forceInB_out);
	}

	int forceForward(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
	{
		auto &robot = static_cast<Robots::RobotBase &>(model);
		auto &param = static_cast<const Robots::WalkParam &>(param_in);

        double forceInF[6]{0,0,0,0,0,0};
        double forceInB[6]{0,0,0,0,0,0};
		double maxForce{200};
		static bool actFlag=false;
		static int countIter;
		static Robots::WalkParam realParam=param;
		static double pIn1[18];
		static double pIn2[18];
		double pIn[18];

        forceInit(param.count, param.force_data->at(0).fce, forceInF);
        aris::dynamic::s_f2f(*robot.forceSensorMak().prtPm(),forceInF,forceInB);

		if((fabs(forceInB[0])>maxForce || fabs(forceInB[1])>maxForce || fabs(forceInB[2])>maxForce) && actFlag==false)
		{
			rt_printf("force exceed limit at count %d, start force adjust\n",param.count);
			actFlag=true;
			realParam.count=param.count;

			countIter=0;
			Robots::walkGait(robot,param);
			robot.GetPin(pIn2);
			param.count--;
			Robots::walkGait(robot,param);
			robot.GetPin(pIn1);
		}

		if(actFlag==false)
		{
			return Robots::walkGait(robot,param);
		}
		else
		{
			if(countIter<=20)
			{
				for(int i=0;i<18;i++)
				{
					pIn[i]=pIn1[i]+(pIn2[i]-pIn1[i])*countIter*(1-(double)countIter/20);
				}
				countIter++;
				robot.SetPin(pIn);
				return 1;
			}
			else
			{
				realParam.count--;
				Robots::walkGait(robot,realParam);

				if(realParam.count%100==0)
				{
					rt_printf("%d,%d,%f,%f,%f\n",param.count,realParam.count,forceInB[0],forceInB[1],forceInB[2]);
				}
				if(realParam.count%param.totalCount==0)
				{
					actFlag=false;
				}
				return realParam.count%param.totalCount;
			}
		}
	}

//    std::atomic_bool isForce;
//    std::atomic_bool isContinue;

//    void parseContinueMoveBegin(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
//    {
//        ContinueMoveParam param;

//        for(auto &i:params)
//        {
//            if(i.first=="u")
//            {
//                moveDir[0]=std::stoi(i.second);
//            }
//            else if(i.first=="v")
//            {
//                moveDir[1]=std::stoi(i.second);
//            }
//            else if(i.first=="w")
//            {
//                moveDir[2]=std::stoi(i.second);
//            }
//            else if(i.first=="yaw")
//            {
//                moveDir[3]=std::stoi(i.second);
//            }
//            else if(i.first=="pitch")
//            {
//                moveDir[4]=std::stoi(i.second);
//            }
//            else if(i.first=="roll")
//            {
//                moveDir[5]=std::stoi(i.second);
//            }
//            else
//            {
//                std::cout<<"parse failed"<<std::endl;
//            }
//        }

//        isContinue=true;
//        isForce=true;

//        msg.copyStruct(param);

//        std::cout<<"finished parse"<<std::endl;
//    }

//    void parseContinueMoveJudge(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
//    {
//        for(auto &i:params)
//        {
//            if(i.first=="isStop")
//            {
//                if(i.second=="0")
//                    isContinue=true;
//                else
//                    isContinue=false;
//            }
//            else if(i.first=="isForce")
//            {
//                if(i.second=="1")
//                    isForce=true;
//                else
//                    isForce=false;
//            }
//            else if(i.first=="u")
//            {
//                moveDir[0]=std::stoi(i.second);
//            }
//            else if(i.first=="v")
//            {
//                moveDir[1]=std::stoi(i.second);
//            }
//            else if(i.first=="w")
//            {
//                moveDir[2]=std::stoi(i.second);
//            }
//            else if(i.first=="yaw")
//            {
//                moveDir[3]=std::stoi(i.second);
//            }
//            else if(i.first=="pitch")
//            {
//                moveDir[4]=std::stoi(i.second);
//            }
//            else if(i.first=="roll")
//            {
//                moveDir[5]=std::stoi(i.second);
//            }
//            else
//            {
//                std::cout<<"parse failed"<<std::endl;
//            }
//        }

//        std::cout<<"finished parse"<<std::endl;
//    }

//    /*****C must be adjusted when used on different robots*****/
//    int continueMove(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
//    {
//        auto &robot = static_cast<Robots::RobotBase &>(model);
//        auto &param = static_cast<const ContinueMoveParam &>(param_in);

//        double bodyVel[6];
//        double bodyAcc[6];
//        double bodyPm[4][4];
//        double deltaPE[6];
//        double deltaPm[4][4];
//        double realPE[6];
//        double realPm[4][4];
//        double nowPee[18];

//        double Fbody[6]{0,0,0,0,0,0};
//        double C[6]{30,30,30,30,30,30};
//        double M[6]{1,1,1,1,1,1};
//        double deltaT{0.001};
//        double forceRange[6]{30,30,30,20,20,20};

//        static ForceTaskParamBase FTP;

//        if(isContinue==true)
//        {
//            //rt_printf("gait continuing\n");
//            if(isForce==false)
//            {
//                for (int i=0;i<6;i++)
//                {
//                    Fbody[i]=moveDir[i];
//                }
//            }
//            else
//            {
//                //initialize
//                if (param.count==0)
//                {
//                    for(int i=0;i<6;i++)
//                    {
//                        FTP.bodyVel_last[i]=0;
//                    }
//                }
//                double forceInF[6];
//                ForceTask::forceInit(param.count,param.force_data->at(0).fce,forceInF);
//                aris::dynamic::s_f2f(*robot.forceSensorMak().prtPm(),forceInF,FTP.forceInB);

//                //Find the max force direction & Set the direction as the move dircetion
//                int num1;
//                int num2;
//                double fmax{0};
//                double mmax{0};
//                for (int i=0;i<3;i++)
//                {
//                    if (fabs(FTP.forceInB[i])>fabs(fmax))
//                    {
//                        fmax=FTP.forceInB[i];
//                        num1=i;
//                    }
//                    if (fabs(FTP.forceInB[i+3])>fabs(mmax))
//                    {
//                        mmax=FTP.forceInB[i+3];
//                        num2=i+3;
//                    }
//                }

//                //Judge whether to rotate in priority. If not, then judge whether to translate
//                if(FTP.forceInB[num2]>forceRange[num2])
//                {
//                    Fbody[num2]=1;
//                }
//                else if(FTP.forceInB[num2]<-forceRange[num2])
//                {
//                    Fbody[num2]=-1;
//                }
//                else
//                {
//                    if(FTP.forceInB[num1]>forceRange[num1])
//                    {
//                        Fbody[num1]=1;
//                    }
//                    else if(FTP.forceInB[num1]<-forceRange[num1])
//                    {
//                        Fbody[num1]=-1;
//                    }
//                }
//            }

//            for (int i=0;i<6;i++)
//            {
//                bodyAcc[i]=(Fbody[i]-C[i]*FTP.bodyVel_last[i])/M[i];
//                bodyVel[i]=FTP.bodyVel_last[i]+bodyAcc[i]*deltaT;
//                deltaPE[i]=bodyVel[i]*deltaT;
//            }

//            robot.GetPmb(*bodyPm);
//            robot.GetPee(nowPee);
//            aris::dynamic::s_pe2pm(deltaPE,*deltaPm,"213");
//            aris::dynamic::s_pm_dot_pm(*bodyPm,*deltaPm,*realPm);
//            aris::dynamic::s_pm2pe(*realPm,realPE,"313");

//            robot.SetPeb(realPE);
//            robot.SetPee(nowPee);

//            if(param.count%1000==0)
//            {
//                rt_printf("count:%d\n",param.count);
//                //rt_printf("bodyPm:%f,%f,%f,%f\n%f,%f,%f,%f\n%f,%f,%f,%f\n%f,%f,%f,%f\n",bodyPm[0][0],bodyPm[0][1],bodyPm[0][2],bodyPm[0][3],bodyPm[1][0],bodyPm[1][1],bodyPm[1][2],bodyPm[1][3],bodyPm[2][0],bodyPm[2][1],bodyPm[2][2],bodyPm[2][3],bodyPm[3][0],bodyPm[3][1],bodyPm[3][2],bodyPm[3][3]);
//                rt_printf("RawForce:%f,%f,%f\n",param.force_data->at(0).Fx,param.force_data->at(0).Fy,param.force_data->at(0).Fz);
//                rt_printf("Force:%f,%f,%f,%f,%f,%f\n",FTP.forceInB[0],FTP.forceInB[1],FTP.forceInB[2],FTP.forceInB[3],FTP.forceInB[4],FTP.forceInB[5]);
//                rt_printf("realPE:%f,%f,%f,%f,%f,%f\n\n",realPE[0],realPE[1],realPE[2],realPE[3],realPE[4],realPE[5]);
//            }

//            memcpy(FTP.bodyVel_last,bodyVel,sizeof(double)*6);
//            return 1;
//        }

//        else
//        {
//            //rt_printf("gait stopping\n");
//            for (int i=0;i<6;i++)
//            {
//                bodyAcc[i]=(Fbody[i]-C[i]*FTP.bodyVel_last[i])/M[i];
//                bodyVel[i]=FTP.bodyVel_last[i]+bodyAcc[i]*deltaT;
//                deltaPE[i]=bodyVel[i]*deltaT;
//            }

//            robot.GetPmb(*bodyPm);
//            robot.GetPee(nowPee);
//            aris::dynamic::s_pe2pm(deltaPE,*deltaPm,"213");
//            aris::dynamic::s_pm_dot_pm(*bodyPm,*deltaPm,*realPm);
//            aris::dynamic::s_pm2pe(*realPm,realPE,"313");

//            robot.SetPeb(realPE);
//            robot.SetPee(nowPee);

//            memcpy(FTP.bodyVel_last,bodyVel,sizeof(double)*6);

//            if ( fabs(bodyVel[0])<1e-10 && fabs(bodyVel[1])<1e-10 && fabs(bodyVel[2])<1e-10 && fabs(bodyVel[3])<1e-10 && fabs(bodyVel[4])<1e-10 && fabs(bodyVel[5])<1e-10)
//                return 0;
//            else
//                return 1;
//        }
//    }

    double ForceWalk::forwardAcc;
    double ForceWalk::turnAcc;
    int ForceWalk::totalCount;
    int ForceWalk::totalCount_tmp;
    double ForceWalk::height;
    double ForceWalk::height_tmp;
    double ForceWalk::alpha;
    double ForceWalk::alpha_tmp;

    NormalGait::WalkState ForceWalk::walkState;
    bool ForceWalk::constFlag;
    double ForceWalk::beginXVel;
    double ForceWalk::endXVel;
    double ForceWalk::beginZVel;
    double ForceWalk::endZVel;
    double ForceWalk::beginOmega;
    double ForceWalk::endOmega;

    double ForceWalk::beginPeb[6];
    double ForceWalk::pEB[6];
    NormalGait::GaitPhase ForceWalk::gaitPhase[6];
    double ForceWalk::swingPee[18];
    double ForceWalk::swingBeginPee[18];
    double ForceWalk::swingEndPee[18];
    double ForceWalk::stancePee[18];
    double ForceWalk::stanceBeginPee[18];
    double ForceWalk::stanceEndPee[18];
    double ForceWalk::followBeginPee[18];
    bool ForceWalk::followFlag[6];
    bool ForceWalk::filterFlag[6];
    int ForceWalk::filterCount[6];

    double ForceWalk::initPee[18];
    double ForceWalk::avgRealH;
    double ForceWalk::planH;

    double ForceWalk::bodyEul213[3];
    double ForceWalk::sumEul[3];
    double ForceWalk::avgEul[3];
    double ForceWalk::targetEul[3];
    double ForceWalk::inputEul[3];
    double ForceWalk::inputEul_tmp[3];

    double ForceWalk::forceRaw[36];
    double ForceWalk::forceInF[36];
    double ForceWalk::forceSum[36];
    double ForceWalk::forceAvg[36];

	ForceWalk::ForceWalk()
	{
	}
	ForceWalk::~ForceWalk()
	{
	}

    void ForceWalk::parseForceWalk(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
	{
		ForceWalkParam param;

		for(auto &i:params)
		{
			if(i.first=="init")
			{
				walkState=NormalGait::WalkState::Init;
				msg.copyStruct(param);
			}
            else if(i.first=="forwardAcc")
			{
                forwardAcc=std::stod(i.second);
                walkState=NormalGait::WalkState::ForwardAcc;
			}
            else if(i.first=="betaAcc")
            {
                turnAcc=std::stod(i.second);
                walkState=NormalGait::WalkState::TurnAcc;
            }
			else if(i.first=="stop")
			{
				walkState=NormalGait::WalkState::Stop;
			}
			else if(i.first=="totalCount")
			{
                totalCount=std::stoi(i.second);
			}
			else if(i.first=="height")
			{
                height_tmp=std::stod(i.second);
			}
            else if(i.first=="roll")
			{
                inputEul_tmp[2]=std::stod(i.second);
			}
            else if(i.first=="pitch")
            {
                inputEul_tmp[1]=std::stod(i.second);
            }
            else if(i.first=="alpha")
            {
                alpha_tmp=std::stod(i.second);
            }
			else
			{
				std::cout<<"parse failed"<<std::endl;
			}
		}

		std::cout<<"finished parse"<<std::endl;
	}

    void ForceWalk::forceInit(int count,int legID)
    {
        if(count==0)
        {
            std::fill(forceSum+6*legID,forceSum+6*legID+5,0);
        }
        if(count<100)
        {
            for(int i=6*legID;i<6*legID+6;i++)
            {
                forceSum[i]+=forceRaw[i];
                forceAvg[i]=forceSum[i]/(count+1);
            }
        }
        for(int i=6*legID;i<6*legID+6;i++)
        {
            forceInF[i]=forceRaw[i]-forceAvg[i];
        }
    }

    void ForceWalk::swingLegTg(const aris::dynamic::PlanParamBase &param_in, int legID)
    {
        auto &param = static_cast<const ForceWalkParam &>(param_in);

        int period_count = param.count%totalCount;
        if(period_count==0)
        {
            double bodyXDist=(beginXVel+endXVel)*totalCount/2*0.001;
            double bodyZDist=(beginZVel+endZVel)*totalCount/2*0.001;
            double bodyAngle=(beginOmega+endOmega)*totalCount/2*0.001;
            double peOmega[6]{0,0,0,PI/2,endOmega*totalCount/2*0.001+bodyAngle,-PI/2};
            double pmOmega[4][4];
            aris::dynamic::s_pe2pm(peOmega,*pmOmega);
            aris::dynamic::s_pm_dot_pnt(*pmOmega,initPee+3*legID,swingEndPee+3*legID);//rotate
            swingEndPee[3*legID]=swingEndPee[3*legID]-(endXVel*totalCount/2*0.001+bodyXDist);//translate
            swingEndPee[3*legID+2]=swingEndPee[3*legID+2]-(endZVel*totalCount/2*0.001+bodyZDist);//translate
        }

        if(period_count>=(totalCount/2-100))
        {
            memcpy(forceRaw+6*legID,param.force_data->at(legID).fce,6*sizeof(double));
            forceInit(period_count-totalCount/2+100,legID);
        }
        //if((legID==1 || legID==4) && period_count>=(totalCount/2-100) && period_count<(totalCount/2+100))
        //rt_printf("count:%d,leg:%d,forceInF:%.4f\n",period_count,legID,forceInF[6*legID+2]);

        const double delayTouch=asin(0/height);
        double s=-(PI/2+delayTouch/2)*cos(PI*(period_count+1)/totalCount)+PI/2+delayTouch/2;//0-PI
        swingPee[3*legID]=(swingBeginPee[3*legID]+swingEndPee[3*legID])/2-(swingEndPee[3*legID]-swingBeginPee[3*legID])/2*cos(s);
        swingPee[3*legID+2]=(swingBeginPee[3*legID+2]+swingEndPee[3*legID+2])/2-(swingEndPee[3*legID+2]-swingBeginPee[3*legID+2])/2*cos(s);
        if(s<PI/2)
        {
            swingPee[3*legID+1]=swingBeginPee[3*legID+1]+(height+swingEndPee[3*legID+1]-swingBeginPee[3*legID+1])*sin(s);
        }
        else
        {
            swingPee[3*legID+1]=swingEndPee[3*legID+1]+height*sin(s);
        }

        if(period_count==totalCount-1)
        {
            //rt_printf("leg:%d, swingBegin:%.4f,%.4f,%.4f\n",legID,swingBeginPee[3*legID],swingBeginPee[3*legID+1],swingBeginPee[3*legID+2]);
            //rt_printf("leg %d, swingPee:  %.4f,%.4f,%.4f\n",legID,swingPee[3*legID],swingPee[3*legID+1],swingPee[3*legID+2]);
        }
    }

    void ForceWalk::stanceLegTg(const aris::dynamic::PlanParamBase &param_in, int legID)
    {
        auto &param = static_cast<const ForceWalkParam &>(param_in);

        int period_count = param.count%totalCount;
        double pe[6]{0,0,0,0,0,0};
        double pm[4][4];
        double stanceEul[3];
        double stancePee_tmp[3];
        memcpy(stancePee,stanceBeginPee,sizeof(stanceBeginPee));

        if(period_count==0)
        {
            std::fill_n(sumEul,3,0);
        }
        else if(period_count>=(totalCount/4-50) && period_count<totalCount/4)
        {
            if(legID==0 || legID==1)
            {
                //param.imu_data->toEulBody2Ground(bodyEul213,PI,"213");
                for (int i=0;i<3;i++)
                {
                    if(bodyEul213[i]>PI)
                    {
                        bodyEul213[i]-=2*PI;
                    }
                    sumEul[i]+=bodyEul213[i];
                    if(period_count==(totalCount/4-1))
                    {
                        avgEul[i]=sumEul[i]/50;
                        //inputEul[i]=inputEul_tmp[i];
                        if(param.count/totalCount==0)
                        {
                            targetEul[i]=avgEul[i];
                        }
                    }
                }
            }
        }    
        else if(period_count>=totalCount/4 && period_count<3*totalCount/4)
        {
            stancePee[3*legID+1]=stanceBeginPee[3*legID+1]+(planH-avgRealH)/2*(1-cos(PI*(period_count-totalCount/4)/(totalCount/2)));

            for(int i=0;i<3;i++)
            {
                stanceEul[i]=(avgEul[i]-targetEul[i]-inputEul[i])/2*(1-cos(PI*(period_count-totalCount/4)/(totalCount/2)));
            }
            memcpy(pe+3,stanceEul,sizeof(stanceEul));
            aris::dynamic::s_pe2pm(pe,*pm,"213");
            memcpy(stancePee_tmp,stancePee+3*legID,sizeof(stancePee_tmp));
            aris::dynamic::s_pm_dot_pnt(*pm,stancePee_tmp,stancePee+3*legID);
        }
        else if(period_count>=3*totalCount/4)
        {
            stancePee[3*legID+1]=stanceBeginPee[3*legID+1]+(planH-avgRealH);

            for(int i=0;i<3;i++)
            {
                stanceEul[i]=avgEul[i]-targetEul[i]-inputEul[i];
            }
            memcpy(pe+3,stanceEul,sizeof(stanceEul));
            aris::dynamic::s_pe2pm(pe,*pm,"213");
            memcpy(stancePee_tmp,stancePee+3*legID,sizeof(stancePee_tmp));
            aris::dynamic::s_pm_dot_pnt(*pm,stancePee_tmp,stancePee+3*legID);
        }
        if(period_count==totalCount-1)
        {
            //rt_printf("leg:%d, stancePee:%.4f,%.4f,%.4f\n",legID,stancePee[3*legID],stancePee[3*legID+1],stancePee[3*legID+2]);
        }
    }

	int ForceWalk::forceWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
	{
		auto &robot = static_cast<Robots::RobotBase &>(model);
		auto &param = static_cast<const ForceWalkParam &>(param_in);

        const double frcRange[6]{-100,-100,-100,-100,-100,-100};

        static aris::dynamic::FloatMarker beginMak{robot.ground()};
        if (param.count == 0)
        {
            beginMak.setPrtPm(*robot.body().pm());
            beginMak.update();
            robot.GetPee(initPee,beginMak);
            robot.GetPee(followBeginPee,beginMak);
            robot.GetPee(swingEndPee,beginMak);
            std::fill_n(followFlag,6,false);
            std::fill_n(filterFlag,6,false);
            std::fill_n(filterCount,6,0);
            inputEul_tmp[0]=0;//set yaw const 0
            //totalCount=totalCount_tmp;//param.count%totalCount illegal at count 0

            planH=initPee[1]-0.005;
        }

        if(param.count%totalCount==0)
        {
            beginMak.setPrtPm(*robot.body().pm());
            beginMak.update();
            robot.GetPeb(beginPeb,beginMak,"213");

            //for test
            double pEBinG[6];
            double pEEinG[18];
            robot.GetPeb(pEBinG,"213");
            robot.GetPee(pEEinG);
            rt_printf("count:%d,pEBinG:%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",param.count,pEBinG[0],pEBinG[1],pEBinG[2],pEBinG[3],pEBinG[4],pEBinG[5]);
            rt_printf("count:%d,pEEinG:%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",param.count,pEEinG[0],pEEinG[1],pEEinG[2],pEEinG[3],pEEinG[4],pEEinG[5]);

            //totalCount=totalCount_tmp;
            height=height_tmp;
            alpha=alpha_tmp;
            memcpy(inputEul,inputEul_tmp,sizeof(inputEul_tmp));

            beginXVel=endXVel;
            beginZVel=endZVel;
            beginOmega=endOmega;
            switch(walkState)
            {
            case NormalGait::WalkState::Init:
                beginXVel=0;
                endXVel=0;
                beginZVel=0;
                endZVel=0;
                beginOmega=0;
                endOmega=0;
                constFlag=false;
                break;

            case NormalGait::WalkState::ForwardAcc:
                endXVel=beginXVel+sin(alpha)*forwardAcc*totalCount/1000;
                endZVel=beginZVel+cos(alpha)*forwardAcc*totalCount/1000;
                endOmega=beginOmega;
                constFlag=true;
                break;

            case NormalGait::WalkState::TurnAcc:
                endXVel=beginXVel;
                endZVel=beginZVel;
                endOmega=beginOmega+turnAcc*totalCount/1000;
                constFlag=true;
                break;

            case NormalGait::WalkState::Const:
                endXVel=beginXVel;
                endZVel=beginZVel;
                endOmega=beginOmega;
                break;

            default:
                rt_printf("Error: unknown walkstate!\n");
                break;
            }
        }
        if(param.count%totalCount==(totalCount-1))
        {
            if(walkState==NormalGait::WalkState::ForwardAcc && constFlag==true)
            {
                walkState=NormalGait::WalkState::Const;
                constFlag=false;
                rt_printf("ForwardAcc finished, the coming Const Vel is:(-x)%.4f,(-z)%.4f,Omega is:%.4f\n",endXVel,endZVel,endOmega);
            }
            if(walkState==NormalGait::WalkState::TurnAcc && constFlag==true)
            {
                walkState=NormalGait::WalkState::Const;
                constFlag=false;
                rt_printf("TurnAcc finished, the coming Const Vel is:(-x)%.4f,(-z)%.4f,Omega is:%.4f\n",endXVel,endZVel,endOmega);
            }
        }

        if(param.count%(2*totalCount)==0)
		{
            double min{0};
			for (int i=0;i<3;i++)
			{
                gaitPhase[2*i]=NormalGait::GaitPhase::Swing;
                gaitPhase[2*i+1]=NormalGait::GaitPhase::Stance;
                robot.pLegs[2*i]->GetPee(swingBeginPee+6*i,beginMak);
                robot.pLegs[2*i+1]->GetPee(stanceBeginPee+6*i+3,beginMak);
                if(stanceBeginPee[6*i+4]<min)
                {
                    min=stanceBeginPee[6*i+4];
                }
			}
            avgRealH=min;
            //avgRealH=(stanceBeginPee[3*1+1]+stanceBeginPee[3*3+1]+stanceBeginPee[3*5+1])/3;
		}
        else if(param.count%(2*totalCount)==totalCount)
		{
            double min{0};
			for (int i=0;i<3;i++)
			{
                gaitPhase[2*i]=NormalGait::GaitPhase::Stance;
                gaitPhase[2*i+1]=NormalGait::GaitPhase::Swing;
                robot.pLegs[2*i]->GetPee(stanceBeginPee+6*i,beginMak);
                robot.pLegs[2*i+1]->GetPee(swingBeginPee+6*i+3,beginMak);
                if(stanceBeginPee[6*i+1]<min)
                {
                    min=stanceBeginPee[6*i+1];
                }
			}
            avgRealH=min;
            //avgRealH=(stanceBeginPee[3*0+1]+stanceBeginPee[3*2+1]+stanceBeginPee[3*4+1])/3;
		}

        for(int i=0;i<6;i++)
        {
            if(gaitPhase[i]==NormalGait::GaitPhase::Swing  && param.count%totalCount>(3*totalCount/4) && followFlag[i]==false)
            {
                //detect 5 points to confirm touching ground
                if(forceInF[6*i+2]<frcRange[i] && filterFlag[i]==false)//param.force_data->at(leg2frc[i]).Fz<frcRange[i]
                {
                    filterFlag[i]=true;
                    filterCount[i]=param.count;
                    rt_printf("leg %d detects force:%.4f, going into Follow in 5 ms after count %d\n",i,forceInF[6*i+2],filterCount[i]);
                }
                if(forceInF[6*i+2]>frcRange[i] && filterFlag[i]==true &&
                   (param.count==filterCount[i]+1 || param.count==filterCount[i]+2 || param.count==filterCount[i]+3 || param.count==filterCount[i]+4))
                {
                    filterFlag[i]=false;
                    rt_printf("leg %d gets fake touching signal at count %d\n",i,filterCount[i]);
                    filterCount[i]=0;
                }

                if(filterCount[i]!=0 && param.count>(filterCount[i]+4))
                {
                    rt_printf("leg %d detects force:%.4f, transfer into Follow at count %d\n",i,param.count);
                    gaitPhase[i]=NormalGait::GaitPhase::Follow;
                    followFlag[i]=true;
                    robot.pLegs[i]->GetPee(followBeginPee+3*i,beginMak);

                    filterFlag[i]=false;
                    filterCount[i]=0;
                }
            }
            if(gaitPhase[i]==NormalGait::GaitPhase::Stance && param.count%totalCount==0)
            {
                if(followFlag[i]==false)
                {
                    robot.pLegs[i]->GetPee(followBeginPee+3*i,beginMak);
                }
                else
                {
                    followFlag[i]=false;
                }
            }
        }
        memcpy(pEB,beginPeb,sizeof(beginPeb));
        pEB[0]=beginPeb[0]-(beginXVel*(param.count%totalCount+1)*0.001
               +0.5*(endXVel-beginXVel)/totalCount*(param.count%totalCount+1)*(param.count%totalCount+1)*0.001);
        pEB[2]=beginPeb[2]-(beginZVel*(param.count%totalCount+1)*0.001
               +0.5*(endZVel-beginZVel)/totalCount*(param.count%totalCount+1)*(param.count%totalCount+1)*0.001);
        pEB[3]=beginPeb[3]+(beginOmega*(param.count%totalCount+1)*0.001
               +0.5*(endOmega-beginOmega)/totalCount*(param.count%totalCount+1)*(param.count%totalCount+1)*0.001);;
        robot.SetPeb(pEB,beginMak,"213");
		for (int i=0;i<6;i++)
		{
			if(gaitPhase[i]==NormalGait::GaitPhase::Swing)
			{
                swingLegTg(param,i);
                robot.pLegs[i]->SetPee(swingPee+3*i,beginMak);
			}
            else if(gaitPhase[i]==NormalGait::GaitPhase::Follow)
			{
                robot.pLegs[i]->SetPee(followBeginPee+3*i,beginMak);
			}
            else if(gaitPhase[i]==NormalGait::GaitPhase::Stance)
            {
                stanceLegTg(param,i);
                robot.pLegs[i]->SetPee(stancePee+3*i,beginMak);
            }
            else
            {
                rt_printf("Error: unknown gaitphase!\n");
            }
		}

		if(walkState==NormalGait::WalkState::Stop && param.count%totalCount==(totalCount-1))
        {
			return 0;
		}
		else
		{
			return 1;
		}
    }
}
