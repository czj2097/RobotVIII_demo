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

namespace FastWalk
{
    void TestGetdJacOverPee()//work well in LegBase and Body Frame
	{
		timeval tpstart,tpend;
		float tused;
		gettimeofday(&tpstart,NULL);


		Robots::RobotTypeI rbt;
		rbt.loadXml("/home/hex/Desktop/mygit/RobotVIII_demo/resource/Robot_VIII.xml");

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

        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Jvi.txt",*outputJvi,2001,9);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dJvi_t.txt",*outputdJvi_t,2001,9);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dJvi_t_in.txt",*outputdJvi_t_in,2001,9);



		gettimeofday(&tpend,NULL);
		tused=tpend.tv_sec-tpstart.tv_sec+(double)(tpend.tv_usec-tpstart.tv_usec)/1000000;
		printf("UsedTime:%f\n",tused);
	}

	double FastWalkPY::pIn_acc[900][18];
	double FastWalkPY::pIn_const[1800][18];
	double FastWalkPY::pIn_dec[900][18];
	FastWalkPY::FastWalkPY()
	{
	}
	FastWalkPY::~FastWalkPY()
	{
	}
	void FastWalkPY::parseFastWalkByPY(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
	{
		FastWalkByPYParam param;

		aris::dynamic::dlmread("/home/hex/Desktop/mygit/Robots/src/Robot_Type_I/resource/Robot_VIII/pIn_acc.txt",*pIn_acc);
		aris::dynamic::dlmread("/home/hex/Desktop/mygit/Robots/src/Robot_Type_I/resource/Robot_VIII/pIn_const.txt",*pIn_const);
		aris::dynamic::dlmread("/home/hex/Desktop/mygit/Robots/src/Robot_Type_I/resource/Robot_VIII/pIn_dec.txt",*pIn_dec);
		for(auto &i:params)
		{
			if(i.first=="n")
			{
				param.n=stoi(i.second);
			}
			else
			{
				std::cout<<"parse failed"<<std::endl;
			}
		}
		msg.copyStruct(param);
		std::cout<<"finished parse"<<std::endl;
	}
	int FastWalkPY::fastWalkByPY(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
	{
		auto &robot = static_cast<Robots::RobotBase &>(model);
		auto &param = static_cast<const FastWalkByPYParam &>(param_in);

		if (param.count<param.totalCount)
		{
			robot.SetPin(*pIn_acc+18*param.count);
		}
		else if(param.count<param.totalCount*(2*param.n-1))
		{
			robot.SetPin(*pIn_const+18*((param.count-param.totalCount)%(2*param.totalCount)));
		}
		else
		{
			robot.SetPin(*pIn_dec+18*(param.count-param.totalCount*(2*param.n-1)));
		}
		return param.count-param.totalCount*2*param.n+1;
	}

    void FastWalkPYAnalyse()
	{
		Robots::RobotTypeI rbt;
		rbt.loadXml("/home/hex/Desktop/mygit/RobotVIII_demo/resource/Robot_VIII.xml");

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
		rbt.loadXml("/home/hex/Desktop/mygit/RobotVIII_demo/resource/Robot_VIII.xml");

        printf("wk analyse begin\n");
		Robots::WalkParam param;
        param.alpha=0;
        param.beta=0;
        param.d=0.8;
        param.h=0.1;
        param.totalCount=1335;//Ain reach the limit when t=1335, d=0.8, h=0.1 for RobotVIII
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

        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/wk_pEB.txt",*pEB,2*param.n*param.totalCount,6);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/wk_pEE.txt",*pEE,2*param.n*param.totalCount,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/wk_pEE_B.txt",*pEE_B,2*param.n*param.totalCount,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/wk_pIn.txt",*pIn,2*param.n*param.totalCount,18);
	}

	/*analyse the maxVel of Pee in workspace*/
    void maxCal(aris::dynamic::Model &model, double * direction, int legID, double vLmt, double *maxVel)
	{
		auto &robot = static_cast<Robots::RobotBase &>(model);
		double Kv{0};
		double Jvi[3][3];
		double tem[3];

		//maxVel
		robot.pLegs[legID]->GetJvi(*Jvi,robot.body());
		aris::dynamic::s_dgemm(3,1,3,1,*Jvi,3,direction,1,0,tem,1);
		for (int i=0;i<3;i++)
		{
			if(i==0)
			{
                Kv=vLmt/fabs(tem[i]);
			}
			else
			{
                if(Kv>=vLmt/fabs(tem[i]))
				{
                    Kv=vLmt/fabs(tem[i]);
				}
			}
		}
		for (int i=0;i<3;i++)
		{
			maxVel[i]=Kv*direction[i];
		}

        printf("maxVel:%f,%f,%f\n",maxVel[0],maxVel[1],maxVel[2]);

//		//maxAcc
//        double aLmt[2]{-3.2,3.2};
//        double Ka{0};
//        double dJvi[3][3];
//        robot.pLegs[legID]->SetVee(*maxVel,robot.body());
//		robot.pLegs[legID]->GetDifJvi(*dJvi,robot.body());
	}

	void maxVel()
	{
		timeval tpstart,tpend;
		float tused;
		gettimeofday(&tpstart,NULL);

		Robots::RobotTypeI rbt;
		rbt.loadXml("/home/hex/Desktop/mygit/RobotVIII_demo/resource/Robot_VIII.xml");

		double initPeb[6]{0,0,0,0,0,0};
		double initPee[18]{ -0.3, -0.85, -0.65,
						   -0.45, -0.85, 0,
							-0.3, -0.85, 0.65,
							 0.3, -0.85, -0.65,
							0.45, -0.85, 0,
							 0.3, -0.85, 0.65 };
		double pEE[18];
        double vLmt {0.9};
		int stepHNum=11;
		int stepDNum=21;
		int angleNum=36;
		double maxVel[stepDNum*stepHNum*angleNum][18];
        //double maxAcc[stepDNum*stepHNum*angleNum][18];

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
                    double alpha=(double)a/angleNum*PI;
                    double direction[3] {0,cos(alpha-PI/2),sin(alpha-PI/2)};
					for (int k=0;k<6;k++)
					{
                        maxCal(rbt,direction,k,vLmt,*maxVel+((i*stepHNum+j)*angleNum+a)*18+3*k);
					}
				}
			}
		}

		aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/Server/maxVel.txt",*maxVel,stepDNum*stepHNum*angleNum,18);

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
		if(count<100)
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

	std::atomic_bool isForce;
	std::atomic_bool isContinue;
	std::atomic_int moveDir[6];
	std::atomic_bool isPull;
    std::atomic_bool isLeft;
	std::atomic_bool isConfirm;

	void parseContinueMoveBegin(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
	{
		ContinueMoveParam param;

		for(auto &i:params)
		{
			if(i.first=="u")
			{
				moveDir[0]=stoi(i.second);
			}
			else if(i.first=="v")
			{
				moveDir[1]=stoi(i.second);
			}
			else if(i.first=="w")
			{
				moveDir[2]=stoi(i.second);
			}
			else if(i.first=="yaw")
			{
				moveDir[3]=stoi(i.second);
			}
			else if(i.first=="pitch")
			{
				moveDir[4]=stoi(i.second);
			}
			else if(i.first=="roll")
			{
				moveDir[5]=stoi(i.second);
			}
			else
			{
				std::cout<<"parse failed"<<std::endl;
			}
		}

		isContinue=true;
		isForce=true;

		msg.copyStruct(param);

		std::cout<<"finished parse"<<std::endl;
	}

	void parseContinueMoveJudge(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
	{
		for(auto &i:params)
		{
            if(i.first=="isStop")
			{
				if(i.second=="0")
					isContinue=true;
				else
					isContinue=false;
			}
			else if(i.first=="isForce")
			{
				if(i.second=="1")
					isForce=true;
				else
					isForce=false;
			}
			else if(i.first=="u")
			{
				moveDir[0]=stoi(i.second);
			}
			else if(i.first=="v")
			{
				moveDir[1]=stoi(i.second);
			}
			else if(i.first=="w")
			{
				moveDir[2]=stoi(i.second);
			}
			else if(i.first=="yaw")
			{
				moveDir[3]=stoi(i.second);
			}
			else if(i.first=="pitch")
			{
				moveDir[4]=stoi(i.second);
			}
			else if(i.first=="roll")
			{
				moveDir[5]=stoi(i.second);
			}
			else
			{
				std::cout<<"parse failed"<<std::endl;
			}
		}

		std::cout<<"finished parse"<<std::endl;
	}

	/*****C must be adjusted when used on different robots*****/
	int continueMove(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
	{
		auto &robot = static_cast<Robots::RobotBase &>(model);
		auto &param = static_cast<const ContinueMoveParam &>(param_in);

		double bodyVel[6];
		double bodyAcc[6];
		double bodyPm[4][4];
		double deltaPE[6];
		double deltaPm[4][4];
		double realPE[6];
		double realPm[4][4];
		double nowPee[18];

		double Fbody[6]{0,0,0,0,0,0};
		double C[6]{30,30,30,30,30,30};
		double M[6]{1,1,1,1,1,1};
		double deltaT{0.001};
		double forceRange[6]{30,30,30,20,20,20};

		static ForceTaskParamBase FTP;

		if(isContinue==true)
		{
			//rt_printf("gait continuing\n");
			if(isForce==false)
			{
				for (int i=0;i<6;i++)
				{
					Fbody[i]=moveDir[i];
				}
			}
			else
			{
				//initialize
				if (param.count==0)
				{
					for(int i=0;i<6;i++)
					{
						FTP.bodyVel_last[i]=0;
					}
				}
                double forceInF[6];
                forceInit(param.count,param.force_data->at(0).fce,forceInF);
                aris::dynamic::s_f2f(*robot.forceSensorMak().prtPm(),forceInF,FTP.forceInB);

				//Find the max force direction & Set the direction as the move dircetion
				int num1;
				int num2;
				double fmax{0};
				double mmax{0};
				for (int i=0;i<3;i++)
				{
					if (fabs(FTP.forceInB[i])>fabs(fmax))
					{
						fmax=FTP.forceInB[i];
						num1=i;
					}
					if (fabs(FTP.forceInB[i+3])>fabs(mmax))
					{
						mmax=FTP.forceInB[i+3];
						num2=i+3;
					}
				}

				if(FTP.forceInB[num2]>forceRange[num2])
				{
					Fbody[num2]=1;
				}
				else if(FTP.forceInB[num2]<-forceRange[num2])
				{
					Fbody[num2]=-1;
				}
				else
				{
					if(FTP.forceInB[num1]>forceRange[num1])
					{
						Fbody[num1]=1;
					}
					else if(FTP.forceInB[num1]<-forceRange[num1])
					{
						Fbody[num1]=-1;
					}
				}
			}

			for (int i=0;i<6;i++)
			{
				bodyAcc[i]=(Fbody[i]-C[i]*FTP.bodyVel_last[i])/M[i];
				bodyVel[i]=FTP.bodyVel_last[i]+bodyAcc[i]*deltaT;
				deltaPE[i]=bodyVel[i]*deltaT;
			}

			robot.GetPmb(*bodyPm);
			robot.GetPee(nowPee);
			double pBody[6];
			aris::dynamic::s_pe2pm(deltaPE,*deltaPm,"213");
			aris::dynamic::s_pm_dot_pm(*bodyPm,*deltaPm,*realPm);
			aris::dynamic::s_pm2pe(*realPm,realPE,"313");


			robot.SetPeb(realPE);
			robot.SetPee(nowPee);

			if(param.count%1000==0)
			{
				rt_printf("count:%d\n",param.count);
				//rt_printf("bodyPm:%f,%f,%f,%f\n%f,%f,%f,%f\n%f,%f,%f,%f\n%f,%f,%f,%f\n",bodyPm[0][0],bodyPm[0][1],bodyPm[0][2],bodyPm[0][3],bodyPm[1][0],bodyPm[1][1],bodyPm[1][2],bodyPm[1][3],bodyPm[2][0],bodyPm[2][1],bodyPm[2][2],bodyPm[2][3],bodyPm[3][0],bodyPm[3][1],bodyPm[3][2],bodyPm[3][3]);
				rt_printf("RawForce:%f,%f,%f\n",param.force_data->at(0).Fx,param.force_data->at(0).Fy,param.force_data->at(0).Fz);
				rt_printf("Force:%f,%f,%f,%f,%f,%f\n",FTP.forceInB[0],FTP.forceInB[1],FTP.forceInB[2],FTP.forceInB[3],FTP.forceInB[4],FTP.forceInB[5]);
				rt_printf("realPE:%f,%f,%f,%f,%f,%f\n\n",realPE[0],realPE[1],realPE[2],realPE[3],realPE[4],realPE[5]);
			}

			memcpy(FTP.bodyVel_last,bodyVel,sizeof(double)*6);
			return 1;
		}

		else
		{
			//rt_printf("gait stopping\n");
			for (int i=0;i<6;i++)
			{
				bodyAcc[i]=(Fbody[i]-C[i]*FTP.bodyVel_last[i])/M[i];
				bodyVel[i]=FTP.bodyVel_last[i]+bodyAcc[i]*deltaT;
				deltaPE[i]=bodyVel[i]*deltaT;
			}

			robot.GetPmb(*bodyPm);
			robot.GetPee(nowPee);
			aris::dynamic::s_pe2pm(deltaPE,*deltaPm,"213");
			aris::dynamic::s_pm_dot_pm(*bodyPm,*deltaPm,*realPm);
			aris::dynamic::s_pm2pe(*realPm,realPE,"313");

			robot.SetPeb(realPE);
			robot.SetPee(nowPee);

			memcpy(FTP.bodyVel_last,bodyVel,sizeof(double)*6);

			if ( fabs(bodyVel[0])<1e-10 && fabs(bodyVel[1])<1e-10 && fabs(bodyVel[2])<1e-10 && fabs(bodyVel[3])<1e-10 && fabs(bodyVel[4])<1e-10 && fabs(bodyVel[5])<1e-10)
				return 0;
			else
				return 1;
		}
	}
}
