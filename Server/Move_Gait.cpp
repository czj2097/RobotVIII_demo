#include "Move_Gait.h"
#include "rtdk.h"

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

Pipe<ForceTask::OpenDoorParam> openDoorPipe(true);

void parseMoveWithRotate(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg)
{
	MoveRotateParam param;

	for(auto &i:params)
    {
        if(i.first=="u")
        {
            param.targetBodyPE213[0]=stod(i.second);
        }
        else if(i.first=="v")
		{
			param.targetBodyPE213[1]=stod(i.second);
		}
        else if(i.first=="w")
		{
			param.targetBodyPE213[2]=stod(i.second);
		}
        else if(i.first=="yaw")
        {
            param.targetBodyPE213[3]=stod(i.second)*PI/180;
        }
        else if(i.first=="pitch")
		{
			param.targetBodyPE213[4]=stod(i.second)*PI/180;
		}
        else if(i.first=="roll")
		{
			param.targetBodyPE213[5]=stod(i.second)*PI/180;
		}
        else if(i.first=="totalCount")
        {
        	param.totalCount=stoi(i.second);
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

	double pBody[6];

	robot.SetPeb(realBodyPE213,"213");
	robot.SetPee(pEE);

	return param.totalCount - param.count - 1;
}

namespace FastWalk
{

	void fastTg()
	{
		timeval tpstart,tpend;
		float tused;
		gettimeofday(&tpstart,NULL);


		Robots::RobotTypeI rbt;
		rbt.loadXml("/home/hex/Desktop/mygit/RobotVIII_demo/resource/Robot_VIII.xml");

		double f_s[18];
		double df_s[18];
		double ddf_s[18];
		double s_t[2001];
		double ds_t[2001];
		double dds_t[2001];

		double Jvi[3][3];
		double dJvi_t[3][3]{0,0,0,0,0,0,0,0,0};
		double dJvi_t_tem[3][3]{0,0,0,0,0,0,0,0,0};
		double dJvi_t_in[3][3]{0,0,0,0,0,0,0,0,0};
		double outputJvi[2001][9];
		double outputdJvi_t[2001][9];
		double outputdJvi_t_in[2001][9];

		double dJvi_x[3][3];
		double dJvi_y[3][3];
		double dJvi_z[3][3];
		double stepH=0.05;
		double stepD=0.15;
		double initPeb[6]{0,0,0,0,0,0};
		double initVeb[6]{0,0,0,0,0,0};
		double initAeb[6]{0,0,0,0,0,0};
		double initPee[18]{ -0.3, -0.85, -0.65,
						   -0.45, -0.85, 0,
							-0.3, -0.85, 0.65,
							 0.3, -0.85, -0.65,
							0.45, -0.85, 0,
							 0.3, -0.85, 0.65 };
		double vEE[18];
		double aEE[18];
		double vEE_B[3];
		double vLmt[2]{-0.9,0.9};
		double aLmt[2]{-3.2,3.2};


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

				aEE[3*j]=ddf_s[3*j]*ds_t[i]+df_s[3*j]*dds_t[i];
				aEE[3*j+1]=ddf_s[3*j+1]*ds_t[i]+df_s[3*j+1]*dds_t[i];
				aEE[3*j+2]=ddf_s[3*j+2]*ds_t[i]+df_s[3*j+2]*dds_t[i];
			}

			rbt.SetPeb(initPeb);
			rbt.SetPee(f_s);
			rbt.SetVb(initVeb);
			rbt.SetVee(vEE);
			rbt.SetAb(initAeb);
			rbt.SetAee(aEE);

			//for(int j=0;j<6;j++)
			//{
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

		aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/Server/Jvi.txt",*outputJvi,2001,9);
		aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/Server/dJvi_t.txt",*outputdJvi_t,2001,9);
		aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/Server/dJvi_t_in.txt",*outputdJvi_t_in,2001,9);



		gettimeofday(&tpend,NULL);
		tused=tpend.tv_sec-tpstart.tv_sec+(double)(tpend.tv_usec-tpstart.tv_usec)/1000000;
		printf("UsedTime:%f\n",tused);
	}

	void fastTgByPYAnalyse()
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

	void wkByPYAnalyse()
	{
		Robots::RobotTypeI rbt;
		rbt.loadXml("/home/hex/Desktop/mygit/RobotVIII_demo/resource/Robot_VIII.xml");

		Robots::WalkParam param;
		param.alpha=0;
		param.beta=0;
		param.d=1.0;
		param.h=0.1;
		param.totalCount=914;
		param.n=2;

		double pEE[4*param.totalCount][18];
		double pIn[4*param.totalCount][18];

		for(param.count=0;param.count<4*param.totalCount;param.count++)
		{
			Robots::walkGait(rbt,param);
			rbt.GetPee(*pEE+18*param.count,rbt.body());
			rbt.GetPin(*pIn+18*param.count);
		}

		aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/Server/pEE_WK.txt",*pEE,4*param.totalCount,18);
		aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/Server/pIn_WK.txt",*pIn,4*param.totalCount,18);
	}


	/*find gait with maxVel in joint space by iteration*/
	void screwInterpolationTraj()
	{
		Robots::RobotTypeI rbt;
		rbt.loadXml("/home/hex/Desktop/mygit/RobotVIII_demo/resource/Robot_VIII.xml");

		timeval tpstart,tpend;
		float tused;
		gettimeofday(&tpstart,NULL);

		double totalTime{3};

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

		double pEB[6];
		double pEE[18];
		double pEE_B[18];

		double aEE[18];
		double vLmt{0.9};
		double aLmt{3.2};

		int keyPointNum{3};
		double keyPee_B[keyPointNum][18];
		double keyPin[keyPointNum][18];
		double keyVin[keyPointNum][18];

		rbt.SetPeb(pEB);
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

		aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/Server/vIn.txt",*vIn,3000,18);
		aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/Server/pIn.txt",*pIn,3000,18);
		aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/Server/pIn_adjust.txt",*pIn_adjust,totalCount,18);
		aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/Server/pEE_adjust.txt",*pEE_adjust,totalCount,18);

	}


	JointSpaceWalk::JointSpaceWalk()
	{
	}
	JointSpaceWalk::~JointSpaceWalk()
	{
	}

	WalkState JointSpaceWalk::walkState=WalkState::None;

	void JointSpaceWalk::parseJointSpaceFastWalk(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg)
	{
		JointSpaceWalkParam param;

		for(auto &i:params)
		{
			if(i.first=="init")
			{
				walkState=WalkState::Init;
				msg.copyStruct(param);
				break;
			}
			else if(i.first=="acc")
			{
				bodyAcc=stod(i.second);
				walkState=WalkState::Acc;
				break;
			}
			else if(i.first=="dec")
			{
				bodyDec=stod(i.second);
				walkState=WalkState::Dec;
				break;
			}
			else if(i.first=="stop")
			{
				walkState=WalkState::Stop;
				break;
			}
			else
			{
				std::cout<<"parse failed"<<std::endl;
			}
		}

		std::cout<<"finished parse"<<std::endl;
	}

	int JointSpaceWalk::jointSpaceFastWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
	{
		auto &robot = static_cast<Robots::RobotBase &>(model);
		auto &param = static_cast<const JointSpaceWalkParam &>(param_in);

		int totalCount=1000;
		if(param.count%(2*totalCount)<totalCount)
		{
			for (int i=0;i<3;i++)
			{
				gaitPhase[2*i]=true;//swing
				gaitPhase[2*i+1]=false;//stance
			}
		}
		else
		{
			for (int i=0;i<3;i++)
			{
				gaitPhase[2*i]=false;//stance
				gaitPhase[2*i+1]=true;//swing
			}
		}

		switch(walkState)
		{
		case WalkState::None:
			return 1;
			break;

		case WalkState::Init:
			robot.GetPeb(beginPeb);
			robot.GetPee(beginPee);
			std::fill_n(beginVb,6,0);
			std::fill_n(beginVee,18,0);

			return 1;
			break;

		case WalkState::Acc:

			break;

		case WalkState::Const:

			break;

		case WalkState::Dec:

			break;

		case WalkState::Stop:

			break;
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
		rbt.loadXml("/home/hex/Desktop/mygit/RobotVIII_demo/resource/Robot_VIII.xml");

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

		aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/Server/maxVel.txt",*maxVel,stepDNum*stepHNum*angleNum,18);

		gettimeofday(&tpend,NULL);
		tused=tpend.tv_sec-tpstart.tv_sec+(double)(tpend.tv_usec-tpstart.tv_usec)/1000000;
		printf("UsedTime:%f\n",tused);
	}

	/*realtime gait of any given stepD by scaling the Pee, seems not to be reasonalbe*/
	void parseFastWalkByPeeScaling(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg)
	{
		FastWalkByScrewParam param;

		aris::dynamic::dlmread("/home/hex/Desktop/mygit/RobotVIII_demo/Server/pEE_adjust.txt",*param.swingPee);
		for(auto &i:params)
		{
			if(i.first=="bodyForwardAcc")
			{
				param.bodyForwardAcc=stod(i.second);
			}
			else if(i.first=="bodyForwardDec")
			{
				param.bodyForwardDec=stod(i.second);
			}
			else if(i.first=="bodyForwardVel")
			{
				param.bodyForwardVel=stod(i.second);
			}
			else if(i.first=="n")
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

	int fastWalkByPeeScaling(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
	{
		auto &robot = static_cast<Robots::RobotBase &>(model);
		auto &param = static_cast<const FastWalkByScrewParam &>(param_in);

		double maxVel=0.5;//m/s
		double bodyVel;
		int totalCount=1000;
		int accCount=round(param.bodyForwardVel/param.bodyForwardAcc*1000);
		int decCount=round(param.bodyForwardVel/param.bodyForwardDec*1000);
		double initPeb[6]{0,0,0,0,0,0};
		double initPee[18]{ -0.3, -0.85, -0.65,
						   -0.45, -0.85, 0,
							-0.3, -0.85, 0.65,
							 0.3, -0.85, -0.65,
							0.45, -0.85, 0,
							 0.3, -0.85, 0.65 };

		if(param.count<accCount)
		{
			bodyVel=param.bodyForwardAcc*param.count/1000;
		}
		else if(param.count<param.n*2*totalCount-decCount)
		{
			bodyVel=param.bodyForwardVel;
		}
		else
		{
			bodyVel=param.bodyForwardVel-param.bodyForwardDec*(param.count-param.n*2*totalCount+decCount)/1000;
		}

		double ratio=bodyVel/maxVel;
		double pEB[6]{0,0,0,0,0,0};
		double pEE[18];
		double pEB_last[6];

		robot.GetPeb(pEB_last);
		memcpy(pEB,pEB_last,sizeof(pEB_last));
		pEB[2]=pEB_last[2]-bodyVel*param.count/1000;
		if(param.count%(2*totalCount)<totalCount)
		{
			for(int i=0;i<3;i++)
			{
				//swing leg in B
				pEE[6*i]=param.swingPee[param.count%totalCount][6*i];
				pEE[6*i+1]=param.swingPee[param.count%totalCount][6*i+1];
				pEE[6*i+2]=(param.swingPee[param.count%totalCount][6*i+2]-initPee[6*i+2])*ratio+initPee[6*i+2];
				//stance leg in B
				pEE[6*i+3]=initPee[6*i+3];
				pEE[6*i+4]=initPee[6*i+4];
				pEE[6*i+5]=-bodyVel*(param.count%totalCount)/totalCount/2*ratio+initPee[6*i+5];
			}
		}
		else
		{
			for(int i=0;i<3;i++)
			{
				//swing leg in B
				pEE[6*i+3]=param.swingPee[param.count%totalCount-totalCount][6*i+3];
				pEE[6*i+4]=param.swingPee[param.count%totalCount-totalCount][6*i+4];
				pEE[6*i+5]=(param.swingPee[param.count%totalCount-totalCount][6*i+5]-initPee[6*i+5])*ratio+initPee[6*i+5];
				//stance leg in B
				pEE[6*i]=initPee[6*i];
				pEE[6*i+1]=initPee[6*i+1];
				pEE[6*i+2]=-bodyVel*(param.count%totalCount-totalCount)/totalCount/2*ratio+initPee[6*i+2];
			}
		}

		robot.SetPeb(pEB);
		robot.SetPee(pEE,robot.body());

		return param.n*2*totalCount - param.count - 1;
	}

	/*analyse the Pee after scaling*/
	void fastTgByPeeScaling()
	{
		Robots::RobotTypeI robot;
		robot.loadXml("/home/hex/Desktop/mygit/RobotVIII_demo/resource/Robot_VIII.xml");
		FastWalkByScrewParam param;
		aris::dynamic::dlmread("/home/hex/Desktop/mygit/RobotVIII_demo/Server/pEE_adjust.txt",*param.swingPee);

		double maxVel=0.5;//m/s
		double bodyVel;
		int totalCount=1000;
		int accCount=round(param.bodyForwardVel/param.bodyForwardAcc*1000);
		int decCount=round(param.bodyForwardVel/param.bodyForwardDec*1000);
		double initPeb[6]{0,0,0,0,0,0};
		double initPee[18]{ -0.3, -0.85, -0.65,
						   -0.45, -0.85, 0,
							-0.3, -0.85, 0.65,
							 0.3, -0.85, -0.65,
							0.45, -0.85, 0,
							 0.3, -0.85, 0.65 };
		double pEB[6]{0,0,0,0,0,0};
		double ratio;
		double pEE[18];
		double pIn[18];
		double pEB_last[6];
		double pEE_last[6];
		double pEE_output[2*param.n*totalCount][18];
		double pIn_output[2*param.n*totalCount][18];

		for(int count=0;count<totalCount*2*param.n;count++)
		{
			if(count<accCount)
			{
				bodyVel=param.bodyForwardAcc*count/1000;
			}
			else if(count<param.n*2*totalCount-decCount)
			{
				bodyVel=param.bodyForwardVel;
			}
			else
			{
				bodyVel=param.bodyForwardVel-param.bodyForwardDec*(count-param.n*2*totalCount+decCount)/1000;
			}

			ratio=bodyVel/maxVel;
			robot.GetPeb(pEB_last);

			if(count%100==0)
			{
				printf("Ratio:%.4f,Peb:%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",ratio,pEB[0],pEB[1],pEB[2],pEB[3],pEB[4],pEB[5]);
			}
			memcpy(pEB,pEB_last,sizeof(pEB_last));
			pEB[2]=pEB_last[2]-bodyVel/1000;
			if(count%(2*totalCount)<totalCount)
			{
				for(int i=0;i<3;i++)
				{
					//swing leg in B
					pEE[6*i]=param.swingPee[count%totalCount][6*i];
					pEE[6*i+1]=param.swingPee[count%totalCount][6*i+1];
					pEE[6*i+2]=(param.swingPee[count%totalCount][6*i+2]-initPee[6*i+2])*ratio+initPee[6*i+2];
					//stance leg in B
					pEE[6*i+3]=initPee[6*i+3];
					pEE[6*i+4]=initPee[6*i+4];
					pEE[6*i+5]=maxVel*(count%totalCount-totalCount/2)/1000*ratio+initPee[6*i+5];
				}
			}
			else
			{
				for(int i=0;i<3;i++)
				{
					//swing leg in B
					pEE[6*i+3]=param.swingPee[count%totalCount][6*i+3];
					pEE[6*i+4]=param.swingPee[count%totalCount][6*i+4];
					pEE[6*i+5]=(param.swingPee[count%totalCount][6*i+5]-initPee[6*i+5])*ratio+initPee[6*i+5];
					//stance leg in B
					pEE[6*i]=initPee[6*i];
					pEE[6*i+1]=initPee[6*i+1];
					pEE[6*i+2]=maxVel*(count%totalCount-totalCount/2)/1000*ratio+initPee[6*i+2];
				}
			}

			robot.SetPeb(pEB);
			robot.SetPee(pEE,robot.body());
			robot.GetPin(pIn);
			memcpy(*pEE_output+18*count,pEE,sizeof(pEE));
			memcpy(*pIn_output+18*count,pIn,sizeof(pIn));
		}
		aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/Server/pEE_output.txt",*pEE_output,param.n*2*totalCount,18);
		aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/Server/pIn_output.txt",*pIn_output,param.n*2*totalCount,18);
	}

}



namespace ForceTask
{
	void forceInit(int count, const double* forceRaw_in, const double* forcePm, double* forceInB_out)
	{
		static double forceSum[6];
		static double forceAvg[6]{0,0,0,0,0,0};
		double forceInF[6]{0,0,0,0,0,0};
		if(count==0)
		{
			for(int i=0;i<6;i++)
			{
				forceSum[i]=0;
			}
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
			forceInF[i]=forceRaw_in[i]-forceAvg[i];
		}
		aris::dynamic::s_f2f(forcePm,forceInF,forceInB_out);
	}

	int forceWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
	{
		auto &robot = static_cast<Robots::RobotBase &>(model);
		auto &param = static_cast<const Robots::WalkParam &>(param_in);

		double forceInB[6]{0,0,0,0,0,0};
		double maxForce{200};
		static bool actFlag=false;
		static int countIter;
		static Robots::WalkParam realParam=param;
		static double pIn1[18];
		static double pIn2[18];
		double pIn[18];

		forceInit(param.count, param.force_data->at(0).fce, *robot.forceSensorMak().prtPm(), forceInB);

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

	//when this function called, make sure that the robot model is at the state set in last ms.
	void getMaxPin(double* maxPin,aris::dynamic::Model &model,ForceTask::maxVelParam &param_in)
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

	std::atomic_bool isForce;
	std::atomic_bool isContinue;
	std::atomic_int moveDir[6];
	std::atomic_bool isPull;
	std::atomic_bool isConfirm;

	void parseContinueMoveBegin(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg)
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

	void parseContinueMoveJudge(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg)
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
				forceInit(param.count,param.force_data->at(0).fce,*robot.forceSensorMak().prtPm(),FTP.forceInB);

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

	void parseOpenDoorBegin(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg)
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
		isPull=true;
		isContinue=true;
		isConfirm=false;

		msg.copyStruct(param);

		std::cout<<"finished parse"<<std::endl;
	}

	void parseOpenDoorJudge(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg)
	{
		for(auto &i:params)
		{
			if(i.first=="isPull")
			{
				if(i.second=="1")
					isPull=true;
				else
					isPull=false;
			}
			else if(i.first=="isStop")
			{
				if(i.second=="0")
					isContinue=true;
				else
					isContinue=false;
			}
			else if(i.first=="isConfirm")
			{
				if(i.second=="1")
					isConfirm=true;
				else
					isConfirm=false;
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

	void inv3(double * matrix,double * invmatrix)
	{
		double a1=matrix[0];
		double b1=matrix[1];
		double c1=matrix[2];
		double a2=matrix[3];
		double b2=matrix[4];
		double c2=matrix[5];
		double a3=matrix[6];
		double b3=matrix[7];
		double c3=matrix[8];

		double value_matrix=a1*(b2*c3-c2*b3)-a2*(b1*c3-c1*b3)+a3*(b1*c2-c1*b2);

		invmatrix[0]=(b2*c3-c2*b3)/value_matrix;
		invmatrix[1]=(c1*b3-b1*c3)/value_matrix;
		invmatrix[2]=(b1*c2-c1*b2)/value_matrix;
		invmatrix[3]=(c2*a3-a2*c3)/value_matrix;
		invmatrix[4]=(a1*c3-c1*a3)/value_matrix;
		invmatrix[5]=(a2*c1-a1*c2)/value_matrix;
		invmatrix[6]=(a2*b3-a3*b2)/value_matrix;
		invmatrix[7]=(b1*a3-a1*b3)/value_matrix;
		invmatrix[8]=(a1*b2-b1*a2)/value_matrix;
	}

	void crossMultiply(double * vector_in1, double *vector_in2, double * vector_out)
	{
		vector_out[0]=vector_in1[1]*vector_in2[2]-vector_in1[2]*vector_in2[1];
		vector_out[1]=vector_in1[2]*vector_in2[0]-vector_in1[0]*vector_in2[2];
		vector_out[2]=vector_in1[0]*vector_in2[1]-vector_in1[1]*vector_in2[0];
	}

	double dotMultiply(double *vector_in1, double *vector_in2)
	{
		double sum{0};
		for(int i=0;i<3;i++)
		{
			sum+=vector_in1[i]*vector_in2[i];
		}
		return sum;
	}

	double norm(double * vector_in)
	{
		return	sqrt(vector_in[0]*vector_in[0]+vector_in[1]*vector_in[1]+vector_in[2]*vector_in[2]);
	}


	//*****only for Robot VIII*****
	int openDoor(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
	{
		auto &robot = static_cast<Robots::RobotBase &>(model);
		auto &param = static_cast<const ContinueMoveParam &>(param_in);

		//for data record
		double pEE[18];
		double bodyVel123[6]{0,0,0,0,0,0};

		//MoveState: PointLocation
		static double beginBodyPE[6];
		static double beginBodyPm[4][4];
		double invBeginBodyPm[4][4];
		//double location[3][3];
		double invLocation[3][3];
		double planeConst[3]{1,1,1};
		double planeVertical[3]{0,0,0};
		double planeVerticalInB[3]{0,0,0};

		double touchPE[6];

		//PushState
		double followPeeInB[18];
		double now2StartPE[6];
		double xBodyInB[3]{1,0,0};
		double yBodyInB[3]{0,1,0};
		double pushBodyPE313[6];//for pause
		double pushPee[18];//for pause
		double d0=0.43;//distance from the handle to the middle of the door
		double h0=0.05;//height from the start position to the walk through position

		//Force Control
		double Fbody[6]{0,0,0,0,0,0};
		double C[6]{50,50,50,50,50,50};
		double M[6]{1,1,1,1,1,1};
		double deltaT{0.001};
		double ForceRange[2]{10,90};
		double forceRatio{1};//1 on RobotIII, 1000 on RobotVIII & single motor

		//Position Generetion From Force
		double bodyVel[6]{0,0,0,0,0,0};
		double bodyAcc[6];
		double bodyPE[6];//delta PE every ms
		double bodyPm[4][4];//bodyPm to Ground every ms
		double deltaPE[6];
		double deltaPm[4][4];
		double realPm[4][4];
		double realPE[6];

		static OpenDoorParam ODP;
		ODP.count=param.count;

		if(param.count==0)
		{
			robot.GetPeb(beginBodyPE);
			robot.GetPmb(*beginBodyPm);
		}

		if(isContinue==true)
		{
			if (param.count==0)
			{
				//initialize
				for(int i=0;i<6;i++)
				{
					ODP.bodyVel_last[i]=0;
				}
				aris::dynamic::s_pe2pe("313",beginBodyPE,"213",ODP.bodyPE_last);
				ODP.pauseFlag=false;
				ODP.moveState_last=MoveState::None;

				//start param of now2start in LocateAjust
				robot.GetPeb(ODP.startPE);
			}
			forceInit(param.count,param.force_data->at(0).fce,*robot.forceSensorMak().prtPm(),ODP.forceInB);

			switch(ODP.moveState)
			{
			case MoveState::None:
				for (int i=0;i<6;i++)
				{
					Fbody[i]=moveDir[i];
				}

				if (fabs(ODP.forceInB[2])>ForceRange[0] && param.count>100 && ODP.pauseFlag==false)
				{
					ODP.countIter=param.count;
					robot.GetPeb(ODP.pointLocation1);
					memcpy(*ODP.location,ODP.pointLocation1,sizeof(double)*3);
					ODP.moveState=MoveState::PointLocate1;
				}

				break;

			case MoveState::PointLocate1:
				if (param.count-ODP.countIter<5000)
				{
					Fbody[2]=1;
					if(isPull==true)
						Fbody[0]=1;
					else
						Fbody[0]=-1;
				}
				else
				{
					Fbody[2]=-1;
				}

				if (fabs(ODP.forceInB[2])>ForceRange[0] && (param.count-ODP.countIter)>5000)
				{
					ODP.countIter=param.count;
					robot.GetPeb(ODP.pointLocation2);
					memcpy(*ODP.location+3,ODP.pointLocation2,sizeof(double)*3);
					ODP.moveState=MoveState::PointLocate2;
				}

				break;

			case MoveState::PointLocate2:
				if (param.count-ODP.countIter<5000)
				{
					Fbody[2]=1;
					Fbody[1]=1;
				}
				else
				{
					Fbody[2]=-1;
				}

				if (fabs(ODP.forceInB[2])>ForceRange[0] && (param.count-ODP.countIter)>5000)
				{
					ODP.countIter=param.count+1;
					robot.GetPeb(ODP.pointLocation3);
					memcpy(*ODP.location+6,ODP.pointLocation3,sizeof(double)*3);
					ODP.moveState=MoveState::LocateAjust;

					//calculate the plane of the door. ax+by+cz=1,(a,b,c) is the vertical vector of the plane
					inv3(*ODP.location,*invLocation);
					aris::dynamic::s_dgemm(3,1,3,1,*invLocation,3,planeConst,1,1,planeVertical,1);
					aris::dynamic::s_inv_pm(*beginBodyPm,*invBeginBodyPm);//have not rotate, beginBodyPm is right here
					aris::dynamic::s_pm_dot_v3(*invBeginBodyPm,planeVertical,planeVerticalInB);
					ODP.planeYPR[0]=atan(planeVerticalInB[0]/planeVerticalInB[2]);
					ODP.planeYPR[1]=-asin(planeVerticalInB[1]/norm(planeVerticalInB));
					ODP.planeYPR[2]=0;

					//Set now param of now2start in LocateAjust
					robot.GetPeb(ODP.nowPE);
					robot.GetPee(ODP.nowPee);

					//Set rotate param
					if(fabs(ODP.planeYPR[0])>PI/9)
					{
						ODP.walkParam.n=2;
						ODP.walkParam.beta=ODP.planeYPR[0]*0.88/3*2;
					}
					else
					{
						ODP.walkParam.n=1;
						ODP.walkParam.beta=ODP.planeYPR[0]*0.88*2;
					}
					ODP.walkParam.d=0;
					ODP.walkParam.alpha=0;
					ODP.walkParam.totalCount=2000;

					rt_printf("yaw:%f,pitch:%f,roll:%f\n",ODP.planeYPR[0],ODP.planeYPR[1],ODP.planeYPR[2]);
				}

				break;

			case MoveState::LocateAjust:
				//1.now2start
				if(param.count-ODP.countIter<ODP.now2StartCount)
				{
					for (int i=0;i<6;i++)
					{
						now2StartPE[i]=ODP.nowPE[i]+(ODP.startPE[i]-ODP.nowPE[i])/2*(1-cos((param.count-ODP.countIter)*PI/ODP.now2StartCount));
					}

					robot.SetPeb(now2StartPE);
					robot.SetPee(ODP.nowPee);
				}
				//2.rotate
				else
				{
					ODP.walkParam.count=param.count-ODP.countIter-ODP.now2StartCount;
					ODP.ret=Robots::walkGait(robot,ODP.walkParam);

					if(ODP.ret==0)
					{
						ODP.countIter=param.count;
						ODP.moveState=MoveState::Forward;
					}
				}

				robot.GetPmb(*bodyPm);
				robot.GetPeb(bodyPE,"213");
				for(int i=0;i<6;i++)
				{
					bodyVel[i]=bodyPE[i]-ODP.bodyPE_last[i];
				}

				break;

			case MoveState::Forward:
				Fbody[2]=-1;

				if(param.count==ODP.countIter+1)
				{
					robot.GetPeb(ODP.startPE);//Used in PrePush for now2Start
				}

				if(fabs(ODP.forceInB[2])>ForceRange[0])
				{
					ODP.countIter=param.count;
					ODP.moveState=MoveState::Backward;
					robot.GetPeb(touchPE);//For calculation of DoorLocation

					for (int i=0;i<3;i++)
					{
						ODP.vector0[i]=touchPE[i]-ODP.startPE[i];
					}
				}

				break;

			case MoveState::Backward:
				Fbody[2]=1;

				if (param.count-ODP.countIter>750)
				{
					ODP.countIter=param.count;
					robot.GetPee(ODP.endPeeInB,robot.body());
					robot.GetPeb(ODP.beginPE);//For calculation of DoorLocation
					if(isPull==false)
						ODP.moveState=MoveState::Rightward;
					else
						ODP.moveState=MoveState::Leftward;
				}

				break;

			case MoveState::Leftward:

				Fbody[0]=-1;

				if(isPull==true)
				{
					if (fabs(ODP.forceInB[0])>ForceRange[0])
					{
						ODP.countIter=param.count;
						robot.GetPeb(ODP.handlePE);//for PushState: leftWalk
						ODP.moveState=MoveState::Rightward;

						for (int i=0;i<3;i++)
						{
							ODP.vector1[i]=ODP.handlePE[i]-ODP.beginPE[i];
						}
					}
					else if (fabs(ODP.forceInB[0])<=ForceRange[0] && param.count-ODP.countIter>7500)
					{
						ODP.countIter=param.count+1;
						robot.GetPee(ODP.startPeeInB,robot.body());
						robot.GetPeb(ODP.nowPE);
						ODP.moveState=MoveState::Follow;
					}
				}
				else//Push
				{
					if (param.count-ODP.countIter>3000)
					{
						ODP.countIter=param.count;
						ODP.moveState=MoveState::Downward;
						ODP.downwardCount=(int)(PI/2*2500);
						ODP.downwardFlag = false;
						robot.GetPeb(ODP.beginPE);//For calculation of DoorLocation
					}
				}

				break;

			case MoveState::Rightward:

				Fbody[0]=1;

				if(isPull==false)
				{
					if (fabs(ODP.forceInB[0])>ForceRange[0])
					{
						ODP.countIter=param.count;
						robot.GetPeb(ODP.handlePE);
						ODP.moveState=MoveState::Leftward;
					}
					else if (fabs(ODP.forceInB[0])<ForceRange[0] && param.count-ODP.countIter>7500)
					{
						ODP.countIter=param.count+1;
						robot.GetPee(ODP.startPeeInB,robot.body());
						robot.GetPeb(ODP.nowPE);
						ODP.moveState=MoveState::Follow;
					}
				}
				else//Pull
				{
					if (param.count-ODP.countIter>3000)
					{
						ODP.countIter=param.count;
						ODP.moveState=MoveState::Downward;
						ODP.downwardCount=(int)(PI/2*2500);
						ODP.downwardFlag = false;
					}
				}

				break;

			case MoveState::Follow:
				if(param.count-ODP.countIter<ODP.followCount)
				{
					for(int i=0;i<3;i++)
					{
						//leg 0,2,4
						followPeeInB[6*i]=ODP.startPeeInB[6*i]+(ODP.endPeeInB[6*i]-ODP.startPeeInB[6*i])/2*(1-cos((param.count-ODP.countIter)*PI/ODP.followCount));
						followPeeInB[6*i+1]=ODP.startPeeInB[6*i+1]+0.05*(1-cos((param.count-ODP.countIter)*2*PI/ODP.followCount));
						followPeeInB[6*i+2]=ODP.startPeeInB[6*i+2];
						//leg 1,3,5
						followPeeInB[6*i+3]=ODP.startPeeInB[6*i+3];
						followPeeInB[6*i+4]=ODP.startPeeInB[6*i+4];
						followPeeInB[6*i+5]=ODP.startPeeInB[6*i+5];
					}
				}
				else if(param.count-ODP.countIter<2*ODP.followCount)
				{
					for(int i=0;i<3;i++)
					{
						//leg 1,3,5
						followPeeInB[6*i+3]=ODP.startPeeInB[6*i+3]+(ODP.endPeeInB[6*i+3]-ODP.startPeeInB[6*i+3])/2*(1-cos((param.count-ODP.countIter-ODP.followCount)*PI/ODP.followCount));
						followPeeInB[6*i+4]=ODP.startPeeInB[6*i+4]+0.05*(1-cos((param.count-ODP.countIter-ODP.followCount)*2*PI/ODP.followCount));
						followPeeInB[6*i+5]=ODP.startPeeInB[6*i+5];
						//leg 0,2,4
						followPeeInB[6*i]=ODP.endPeeInB[6*i];
						followPeeInB[6*i+1]=ODP.endPeeInB[6*i+1];
						followPeeInB[6*i+2]=ODP.endPeeInB[6*i+2];
					}
				}
				robot.SetPeb(ODP.nowPE);
				robot.SetPee(followPeeInB,robot.body());

				if(param.count-ODP.countIter+1==2*ODP.followCount)
				{
					ODP.countIter=param.count;
					ODP.moveState=MoveState::Forward;
				}

				robot.GetPmb(*bodyPm);
				robot.GetPeb(bodyPE,"213");
				for(int i=0;i<6;i++)
				{
					bodyVel[i]=bodyPE[i]-ODP.bodyPE_last[i];
				}

				break;

			case MoveState::Downward:
				Fbody[1]=-1;

				if (ODP.downwardFlag)
				{
					if(isPull==false)
					{
						Fbody[1]=-cos((param.count-ODP.countIter)*PI/ODP.downwardCount/2);
						Fbody[0]=sin((param.count-ODP.countIter)*PI/ODP.downwardCount/2);
					}
					else
					{
						Fbody[1]=-cos((param.count-ODP.countIter)*PI/ODP.downwardCount/2);
						Fbody[0]=-sin((param.count-ODP.countIter)*PI/ODP.downwardCount/2);
					}
				}

				if(ODP.downwardFlag==false && fabs(ODP.forceInB[1])>20)
				{
					ODP.downwardFlag = true;
					ODP.countIter=param.count;
					robot.GetPeb(touchPE);//For calculation of DoorLocation

					for (int i=0;i<3;i++)
					{
						ODP.vector2[i]=touchPE[i]-ODP.beginPE[i];
					}

					ODP.handleLocation[2]=-norm(ODP.vector0);
					ODP.handleLocation[0]=norm(ODP.vector1);
					ODP.handleLocation[1]=-norm(ODP.vector2);

					rt_printf("handleLocation:%f,%f,%f\n",ODP.handleLocation[0],ODP.handleLocation[1],ODP.handleLocation[2]);
				}

				if (fabs(ODP.forceInB[0])>ForceRange[1] && param.count-ODP.countIter>ODP.downwardCount/2)
				{
					ODP.countIter=param.count;
					if(isPull==true)
						ODP.moveState=MoveState::Pullhandle;
					else
						ODP.moveState=MoveState::Pushhandle;
				}

				break;

			case MoveState::Pullhandle:

				Fbody[2]=1;

				if (param.count-ODP.countIter>2500)
					ODP.moveState=MoveState::Upward;

				break;

			case MoveState::Upward:

				Fbody[1]=1;
				if(isPull==true)
				{
					if (param.count-ODP.countIter>5000)
					{
						isContinue=false;//Stop here when Pull
					}
				}

				break;

			case MoveState::Pushhandle:

				Fbody[2]=-1;

				if (param.count-ODP.countIter>2500)
				{
					ODP.moveState=MoveState::PrePush;
				}
				break;

			case MoveState::PrePush:

				if ( fabs(bodyVel[0])<1e-10 && fabs(bodyVel[1])<1e-10 && fabs(bodyVel[2])<1e-10 && fabs(bodyVel[3])<1e-10 && fabs(bodyVel[4])<1e-10 && fabs(bodyVel[5])<1e-10)
				{
					ODP.countIter=param.count+1;
					ODP.moveState=MoveState::Push;
					ODP.pushState=PushState::now2Start;

					//now2start param
					robot.GetPeb(ODP.nowPE);
					robot.GetPmb(*ODP.nowPm);
					robot.GetPee(ODP.nowPee);
					for(int i=0;i<3;i++)
					{
						ODP.now2startDistance[i]=ODP.startPE[i]-ODP.nowPE[i]; //startPE recorded at the last forward
						ODP.handle2startDistance[i]=ODP.startPE[i]-ODP.handlePE[i];
					}
					aris::dynamic::s_pm_dot_v3(*ODP.nowPm,xBodyInB,ODP.xNowInG);
					aris::dynamic::s_pm_dot_v3(*ODP.nowPm,yBodyInB,ODP.yNowInG);

					ODP.now2startDistanceInB[0]=dotMultiply(ODP.now2startDistance,ODP.xNowInG);
					ODP.now2startDistanceInB[1]=dotMultiply(ODP.now2startDistance,ODP.yNowInG)+h0;
					ODP.now2startDistanceInB[2]=0;

					aris::dynamic::s_pm_dot_v3(*ODP.nowPm,ODP.now2startDistanceInB,ODP.now2startDistanceInG);
				}

				break;

			case MoveState::Push:
				//Three Step using position control
				switch(ODP.pushState)
				{
				case PushState::now2Start:
					//1.Move back to startPE;
					for (int i=0;i<6;i++)
					{
						now2StartPE[i]=ODP.nowPE[i]+(ODP.now2startDistanceInG[i])/2*(1-cos((param.count-ODP.countIter)*PI/ODP.now2StartCount));
					}

					robot.SetPeb(now2StartPE);
					robot.SetPee(ODP.nowPee);

					if(param.count-ODP.countIter+1==ODP.now2StartCount)
					{
						if(isConfirm==true)
						{
							ODP.pushState=PushState::leftWalk;
							ODP.countIter=param.count+1;

							ODP.walkParam.n=2;
							ODP.walkParam.alpha=PI/2;
							ODP.walkParam.beta=0;
							ODP.walkParam.totalCount=2000;
							ODP.walkParam.d=(d0-fabs(dotMultiply(ODP.handle2startDistance,ODP.xNowInG)))/3*2;
						}
						else//pause, tested useless
						{
							robot.SetPeb(ODP.startPE);
							robot.SetPee(ODP.nowPee);
						}
					}

					break;

				case PushState::leftWalk:
					//2.Move rightward to align with the door;
					ODP.walkParam.count=param.count-ODP.countIter;

					ODP.ret=Robots::walkGait(robot,ODP.walkParam);

					if(ODP.ret==0)
					{
						if(isConfirm==true)
						{
							ODP.pushState=PushState::forwardWalk;
							ODP.countIter=param.count+1;

							ODP.walkParam.totalCount=2000;
							ODP.walkParam.n=4;
							ODP.walkParam.alpha=0;
							ODP.walkParam.beta=0;
							ODP.walkParam.d=0.5;
						}
						else//for pause, teseted useless
						{
							robot.GetPee(pushPee);
							robot.GetPeb(pushBodyPE313);
							robot.SetPeb(pushBodyPE313);
							robot.SetPee(pushPee);
						}
					}

					break;

				case PushState::forwardWalk:
					//3.Move through the door
					ODP.walkParam.count=param.count-ODP.countIter;

					ODP.ret=Robots::walkGait(robot,ODP.walkParam);

					if(ODP.ret==0)
					{
						if(isConfirm==true)
						{
							return 0;
						}
						else//for pause, teseted useless
						{
							robot.GetPee(pushPee);
							robot.GetPeb(pushBodyPE313);
							robot.SetPeb(pushBodyPE313);
							robot.SetPee(pushPee);
						}
					}

					break;

				default:
					break;
				}

				robot.GetPmb(*bodyPm);
				robot.GetPeb(bodyPE,"213");
				for(int i=0;i<6;i++)
				{
					bodyVel[i]=bodyPE[i]-ODP.bodyPE_last[i];
				}

				break;

			default:
				break;
			}

			if(ODP.moveState!=MoveState::Push && ODP.moveState!=MoveState::LocateAjust && ODP.moveState!=MoveState::Follow)
			{
			   if (isConfirm==false && ODP.pauseFlag==false)
				{
					ODP.moveState_last=ODP.moveState;
					ODP.moveState=MoveState::None;
					ODP.pauseFlag=true;
				}
				else if(isConfirm==false && ODP.pauseFlag==true)
				{
					ODP.moveState=MoveState::None;
					ODP.pauseCount++;
					//Do not change F to zeros here, so that we can move the robot manually when paused
				}
				else if(isConfirm==true && ODP.pauseFlag==true)
				{
					ODP.moveState=ODP.moveState_last;
					ODP.pauseFlag=false;
					ODP.pauseCount++;
				}
				else//isConfirm=true && pauseFlag=false
				{
					ODP.countIter+=ODP.pauseCount;
					ODP.pauseCount=0;
				}


				for (int i=0;i<6;i++)
				{
					bodyAcc[i]=(Fbody[i]-C[i]*ODP.bodyVel_last[i])/M[i];
					bodyVel[i]=ODP.bodyVel_last[i]+bodyAcc[i]*deltaT;
					deltaPE[i]=bodyVel[i]*deltaT;
				}

				robot.GetPmb(*bodyPm);
				robot.GetPeb(bodyPE);
				double pBody[6];
				aris::dynamic::s_pe2pm(deltaPE,*deltaPm,"213");
				aris::dynamic::s_pm_dot_pm(*bodyPm,*deltaPm,*realPm);
				aris::dynamic::s_pm2pe(*realPm,realPE,"313");

				robot.GetPee(pEE);
				robot.SetPeb(realPE);
				robot.SetPee(pEE);
			}
			if (param.count%100==0)
			{
				rt_printf("moveState:%d,forceRaw:%f,%f,%f,force:%f,%f,%f\n",ODP.moveState,param.force_data->at(0).Fx,param.force_data->at(0).Fy,param.force_data->at(0).Fz,ODP.forceInB[0],ODP.forceInB[1],ODP.forceInB[2]);
			}

			//Aris::Dynamic::s_pm_dot_pnt(*bodyPm,ODP.toolInR,ODP.toolInG);
			bodyVel123[0]=bodyVel[0];
			bodyVel123[1]=bodyVel[1];
			bodyVel123[2]=bodyVel[2];
			bodyVel123[3]=bodyVel[4];
			bodyVel123[4]=bodyVel[3];
			bodyVel123[5]=bodyVel[5];
			aris::dynamic::s_vp2vp(*bodyPm,bodyVel123,ODP.toolInR,nullptr,ODP.toolInG,ODP.toolVel);

			robot.GetPee(pEE);
			memcpy(ODP.pEE_last,pEE,sizeof(double)*18);
			memcpy(ODP.bodyPE_last,bodyPE,sizeof(double)*6);
			memcpy(ODP.bodyVel_last,bodyVel,sizeof(double)*6);

			openDoorPipe.sendToNrt(ODP);

			return 1;
		}

		else
		{
			//rt_printf("gait stopping\n");
			for (int i=0;i<6;i++)
			{
				bodyAcc[i]=(Fbody[i]-C[i]*ODP.bodyVel_last[i])/M[i];
				bodyVel[i]=ODP.bodyVel_last[i]+bodyAcc[i]*deltaT;
				deltaPE[i]=bodyVel[i]*deltaT;
			}

			robot.GetPmb(*bodyPm);
			robot.GetPeb(bodyPE);
			double pBody[6];
			aris::dynamic::s_pe2pm(deltaPE,*deltaPm,"213");
			aris::dynamic::s_pm_dot_pm(*bodyPm,*deltaPm,*realPm);
			aris::dynamic::s_pm2pe(*realPm,realPE,"313");
			double nowPee[18];
			robot.GetPee(nowPee);
			robot.SetPeb(realPE);
			robot.SetPee(nowPee);

			memcpy(ODP.pEE_last,nowPee,sizeof(double)*18);
			memcpy(ODP.bodyPE_last,bodyPE,sizeof(double)*6);
			memcpy(ODP.bodyVel_last,bodyVel,sizeof(double)*6);

			if ( fabs(bodyVel[0])<1e-10 && fabs(bodyVel[1])<1e-10 && fabs(bodyVel[2])<1e-10 && fabs(bodyVel[3])<1e-10 && fabs(bodyVel[4])<1e-10 && fabs(bodyVel[5])<1e-10)
				return 0;
			else
				return 1;
		}
	}

	void StartRecordData()
	{
		openDoorThread = std::thread([&]()
		{
			struct OpenDoorParam param;
			static std::fstream fileGait;
			std::string name = aris::core::logFileName();
			name.replace(name.rfind("log.txt"), std::strlen("openDoor.txt"), "openDoor.txt");
			fileGait.open(name.c_str(), std::ios::out | std::ios::trunc);

			long long count = -1;
			while (1)
			{
				openDoorPipe.recvInNrt(param);

				//fileGait << ++count << " ";
				fileGait << param.walkParam.count << "  ";
				/*
				fileGait << param.moveState << "  ";

				for (int i = 0; i < 18; i++)
				{
					fileGait << param.pEE_last[i] << "  ";
				}
				for (int i = 0; i < 3; i++)
				{
					fileGait << param.toolInG[i] << "  ";
				}
				for (int i = 0; i < 3; i++)
				{
					fileGait << param.toolVel[i] << "  ";
				}
				for (int i = 0; i < 3; i++)
				{
					fileGait << param.forceInB[i] << "  ";
				}*/
				for (int i=0;i<18;i++)
				{
					fileGait << param.walkParam.motion_raw_data->at(i).feedback_cur << "  ";
				}
				fileGait << std::endl;
			}

			fileGait.close();
		});
	}
}
