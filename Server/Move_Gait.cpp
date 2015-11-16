#include "Move_Gait.h"
#include "rtdk.h"

Aris::Core::MSG parseMove2(const std::string &cmd, const map<std::string, std::string> &params)
{
    MOVES_PARAM  param;

    double targetPos[3]; //移动目标位置

    for(auto &i:params)
    {
        if(i.first=="component")
        {
            if(i.second=="lf")
            {
                param.comID=0;
            }
            else if(i.second=="lm")
            {
                param.comID=1;
            }
            else if(i.second=="lr")
            {
                param.comID=2;
            }
            else if(i.second=="rf")
            {
                param.comID=3;
            }
            else if(i.second=="rm")
            {
                param.comID=4;
            }
            else if(i.second=="rr")
            {
                param.comID=5;
            }
            else if(i.second=="bd")
            {
                param.comID=6;
            }
            else
            {
                std::cout<<"parse failed"<<std::endl;
                return MSG{};
            }
        }
        //绝对坐标移动
        else if(i.first=="x")
        {
            targetPos[0]=std::stod(i.second);
            param.isAbsolute=true;
        }
        else if(i.first=="y")
        {
            targetPos[1]=stod(i.second);
            param.isAbsolute=true;
        }
        else if(i.first=="z")
        {
            targetPos[2]=stod(i.second);
            param.isAbsolute=true;
        }
        //相对坐标移动
        else if(i.first=="u")
        {
            targetPos[0]=stod(i.second);
        }
        else if(i.first=="v")
        {
            targetPos[1]=stod(i.second);
        }
        else if(i.first=="w")
        {
            targetPos[2]=stod(i.second);
        }
        else
        {
            std::cout<<"parse failed"<<std::endl;
            return MSG{};
        }
    }

    if(param.comID==6)
    {
        std::copy_n(targetPos, 3, param.targetBodyPE);
    }
    else
    {
        std::copy_n(targetPos, 3, &param.targetPee[3*param.comID]);
    }

    param.periodCount=3000;

    Aris::Core::MSG msg;
    msg.CopyStruct(param);

    std::cout<<"finished parse"<<std::endl;

    return msg;
}

int move2(Robots::ROBOT_BASE * pRobot, const Robots::GAIT_PARAM_BASE * pParam)
{
    const MOVES_PARAM *pMP = static_cast<const MOVES_PARAM *>(pParam);

    double realTargetPee[18];
    double realTargetPbody[6];
    std::copy_n(pMP->beginPee, 18, realTargetPee);
    std::copy_n(pMP->beginBodyPE, 6, realTargetPbody);

    //绝对坐标
    if (pMP->isAbsolute)
    {
        if(pMP->comID==6)
        {
            std::copy_n(pMP->targetBodyPE, 3, realTargetPbody);
        }
        else
        {
            std::copy_n(&(pMP->targetPee[pMP->comID*3]), 3, &realTargetPee[pMP->comID*3]);
        }
    }
    //相对坐标
    else
    {
        if(pMP->comID==6)
        {
            for(int i=0;i<3;i++)
            {
                realTargetPbody[i]+=pMP->targetBodyPE[i];
            }
        }
        else
        {
            for(int i=0;i<18;i++)
            {
                realTargetPee[i]+=pMP->targetPee[i];
            }
        }
    }

    double s = -(PI / 2)*cos(PI * (pMP->count  + 1) / pMP->periodCount ) + PI / 2;

    /*插值当前的末端和身体位置*/
    double pEE[18], pBody[6];

    for (int i = 0; i < 18; ++i)
    {
        pEE[i] = pMP->beginPee[i] * (cos(s) + 1) / 2 + realTargetPee[i] * (1 - cos(s)) / 2;
    }
    for (int i = 0; i < 6; ++i)
    {
        pBody[i] = pMP->beginBodyPE[i] * (cos(s) + 1) / 2 + realTargetPbody[i] * (1 - cos(s)) / 2;
    }

    pRobot->SetPee(pEE, pBody);

    /*test*/
//    if(pMP->count%500==0)
//    {
//        rt_printf("rr: %f %f %f\n"
//                        , pEE[15], pEE[16], pEE[17]);
//    }

    /*返回剩余的count数*/

    return pMP->periodCount - pMP->count - 1;

}

Aris::Core::MSG parseSwing(const std::string &cmd, const map<std::string, std::string> &params)
{
    SWING_PARAM  param;

    for(auto &i:params)
    {
        if(i.first=="y")
        {
            param.centreP[1]=stod(i.second);
        }
        else if(i.first=="z")
        {
            param.centreP[2]=stod(i.second);
        }
        else if(i.first=="deg")
        {
            param.swingRad=stod(i.second)/180*PI;//计算身体摆动的弧度
        }
        else
        {
            std::cout<<"parse failed"<<std::endl;
            return MSG{};
        }
    }

    param.periodCount=3000;

    Aris::Core::MSG msg;
    msg.CopyStruct(param);

    std::cout<<"finished parse"<<std::endl;

    return msg;
}

int swing(Robots::ROBOT_BASE * pRobot, const Robots::GAIT_PARAM_BASE * pParam)
{
    const SWING_PARAM *pSWP = static_cast<const SWING_PARAM *>(pParam);

    //计算圆弧半径和起始俯仰角
    double radius;
    double beginRad;
    double beginBodyPE123[6];
    Aris::DynKer::s_pe2pe("313", pSWP->beginBodyPE, "123", beginBodyPE123); //将起始身体位姿转换到123欧拉角
    radius = sqrt(pow((pSWP->centreP[1] - pSWP->beginBodyPE[1]),2) + pow((pSWP->centreP[2] - pSWP->beginBodyPE[2]),2));
    beginRad = atan2((pSWP->beginBodyPE[1] - pSWP->centreP[1]) , -(pSWP->beginBodyPE[2] - pSWP->centreP[2]));

    double s = -(pSWP->swingRad / 2)*cos(PI * (pSWP->count  + 1) / pSWP->periodCount ) + pSWP->swingRad / 2;

    /*插值当前的身体位置*/
    double pBody[6];
    double pBody123[6];
    double currentRad;
    currentRad = beginRad + s;
    std::copy_n(beginBodyPE123, 6, pBody123);

    pBody123[1] = pSWP->centreP[1] + radius*sin(currentRad);
    pBody123[2] = pSWP->centreP[2] - radius*cos(currentRad);
    pBody123[3] = beginBodyPE123[3] + s;//当前俯仰角
    Aris::DynKer::s_pe2pe("123", pBody123, "313", pBody);

    pRobot->SetPee(pSWP->beginPee, pBody);

    /*test*/
    if(pSWP->count%500==0)
    {
        rt_printf("beginRad: %f\n"
                  , beginRad);
        rt_printf("pBody: %f %f %f\n"
                        , pBody[3], pBody[4], pBody[5]);
    }

    /*返回剩余的count数*/

    return pSWP->periodCount - pSWP->count - 1;

}

Aris::Core::MSG parseMoveWithRotate(const std::string &cmd, const map<std::string, std::string> &params)
{
	MOVEWITHROTATE_PARAM param;

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
            return MSG{};
        }
    }

    Aris::Core::MSG msg;
    msg.CopyStruct(param);

    std::cout<<"finished parse"<<std::endl;

    return msg;
}

int moveWithRotate(Robots::ROBOT_BASE * pRobot, const Robots::GAIT_PARAM_BASE * pParam)
{
	const MOVEWITHROTATE_PARAM *pMRP = static_cast<const MOVEWITHROTATE_PARAM *>(pParam);

	double beginBodyPE213[6];
	Aris::DynKer::s_pe2pe("313",pMRP->beginBodyPE,"213",beginBodyPE213);

	double realBodyPE213[6];
	for(int i=0;i<6;i++)
	{
		double s = -(pMRP->targetBodyPE213[i] / 2)*cos(PI * (pMRP->count  + 1) / pMRP->totalCount ) + pMRP->targetBodyPE213[i] / 2;
		realBodyPE213[i]=beginBodyPE213[i]+s; //target of current ms
	}

	double pBody[6];
	Aris::DynKer::s_pe2pe("213",realBodyPE213,"313",pBody);

	pRobot->SetPee(pMRP->beginPee,pBody);

	return pMRP->totalCount - pMRP->count - 1;
}

std::atomic_bool isContinue;
std::atomic_bool moveDir[12];
/*Aris::Core::MSG parseContinueMoveBegin(const std::string &cmd, const map<std::string, std::string> &params)
{
	CONTINUEMOVE_PARAM param;

	for(auto &i:params)
	{
		if(i.first=="u")
		{
			param.moveDir[0]=stod(i.second);
		}
		else if(i.first=="v")
		{
			param.moveDir[1]=stod(i.second);
		}
		else if(i.first=="w")
		{
			param.moveDir[2]=stod(i.second);
		}
		else if(i.first=="yaw")
		{
			param.moveDir[3]=stod(i.second);
		}
		else if(i.first=="pitch")
		{
			param.moveDir[4]=stod(i.second);
		}
		else if(i.first=="roll")
		{
			param.moveDir[5]=stod(i.second);
		}
		else
		{
			std::cout<<"parse failed"<<std::endl;
			return MSG{};
		}
	}

	Aris::Core::MSG msg;
	msg.CopyStruct(param);

	std::cout<<"finished parse"<<std::endl;

	return msg;
}

int continueMoveBegin(Robots::ROBOT_BASE * pRobot, const Robots::GAIT_PARAM_BASE * pParam)
{
	const CONTINUEMOVE_PARAM * pCMP = static_cast<const CONTINUEMOVE_PARAM *>(pParam);

	static CM_LAST_PARAM pCMLP;

	double bodyPE[6];
	double bodyVel[6];
	double bodyAcc[6];

	double F[6];
	double C{2.5};
	double M[6]{0.2,0.2,0.2,0.2,0.2,0.2};
	double deltaT{0.001};

	memcpy(F,pCMP->moveDir,sizeof(double)*6);
	//rt_printf("%f,%f,%f,%f,%f,%f\n",F[0],F[1],F[2],F[3],F[4],F[5]);
	//rt_printf("%f,%f,%f,%f,%f,%f\n",C[0],C[1],C[2],C[3],C[4],C[5]);
	//rt_printf("%f,%f,%f,%f,%f,%f\n",M[0],M[1],M[2],M[3],M[4],M[5]);
	for (int i=0;i<6;i++)
	{
		bodyAcc[i]=(F[i]-C*pCMLP.bodyVel_last[i])/M[i];
		bodyVel[i]=pCMLP.bodyVel_last[i]+bodyAcc[i]*deltaT;
		bodyPE[i]=pCMLP.bodyPE_last[i]+bodyVel[i]*deltaT;
	}

	//rt_printf("%d,%f,%f,%f,%f,%f,%f\n",pCMP->count,bodyPE[0],bodyPE[1],bodyPE[2],bodyPE[3],bodyPE[4],bodyPE[5]);

	double pBody[6];
	Aris::DynKer::s_pe2pe("213",bodyPE,"313",pBody);

	pRobot->SetPee(pParam->beginPee,pBody);

	memcpy(pCMLP.bodyPE_last,bodyPE,sizeof(double)*6);
	memcpy(pCMLP.bodyVel_last,bodyVel,sizeof(double)*6);

}*/

Aris::Core::MSG parseContinueMoveJudge(const std::string &cmd, const map<std::string, std::string> &params)
{
    //   cmd -u -v -w -r -p -y -iscontinue -isstart
	CONTINUEMOVE_PARAM param;

	for(auto &i:params)
	{
		if(i.first=="isContinue")
		{
			if(i.second=="1")
				isContinue=true;
            else
                isContinue=false;
		}
    	/*else if(i.first=="isStart")
		{
			if(i.second=="1")
				param.isStart=true;
            else
                param.isStart=false;
        }*/
		else if(i.first=="u")
		{
            param.move_direction[0]=stoi(i.second);

			if(i.second=="1")
			{
				moveDir[0]=true;
				moveDir[1]=false;
			}
			else if(i.second=="-1")
			{
				moveDir[0]=false;
				moveDir[1]=true;
			}
			else
			{
				moveDir[0]=false;
				moveDir[1]=false;
			}
		}
		else if(i.first=="v")
		{
            param.move_direction[1]=stoi(i.second);

			if(i.second=="1")
			{
				moveDir[2]=true;
				moveDir[3]=false;
			}
			else if(i.second=="-1")
			{
				moveDir[2]=false;
				moveDir[3]=true;
			}
			else
			{
				moveDir[2]=false;
				moveDir[3]=false;
			}
		}
		else if(i.first=="w")
		{
            param.move_direction[2]=stoi(i.second);

			if(i.second=="1")
			{
				moveDir[4]=true;
				moveDir[5]=false;
			}
			else if(i.second=="-1")
			{
				moveDir[4]=false;
				moveDir[5]=true;
			}
			else
			{
				moveDir[4]=false;
				moveDir[5]=false;
			}
		}
		else if(i.first=="yaw")
		{
            param.move_direction[3]=stoi(i.second);

			if(i.second=="1")
			{
				moveDir[6]=true;
				moveDir[7]=false;
			}
			else if(i.second=="-1")
			{
				moveDir[6]=false;
				moveDir[7]=true;
			}
			else
			{
				moveDir[6]=false;
				moveDir[7]=false;
			}
		}
		else if(i.first=="pitch")
		{
            param.move_direction[4]=stoi(i.second);

			if(i.second=="1")
			{
				moveDir[8]=true;
				moveDir[9]=false;
			}
			else if(i.second=="-1")
			{
				moveDir[8]=false;
				moveDir[9]=true;
			}
			else
			{
				moveDir[8]=false;
				moveDir[9]=false;
			}
		}
		else if(i.first=="roll")
		{
            param.move_direction[5]=stoi(i.second);

			if(i.second=="1")
			{
				moveDir[10]=true;
				moveDir[11]=false;
			}
			else if(i.second=="-1")
			{
				moveDir[10]=false;
				moveDir[11]=true;
			}
			else
			{
				moveDir[10]=false;
				moveDir[11]=false;
			}
		}
		else
		{
			std::cout<<"parse failed"<<std::endl;
			return MSG{};
		}
	}

	Aris::Core::MSG msg;
	msg.CopyStruct(param);

	std::cout<<"finished parse"<<std::endl;

	return msg;
}

int continueMoveJudge(Robots::ROBOT_BASE * pRobot, const Robots::GAIT_PARAM_BASE * pParam)
{
    double bodyVel[6];
    double bodyAcc[6];
    double bodyPE[6];

    double beginBodyPE213[6];
    double realBodyPE[6];

    double zeros[6];

    double F[6];
    double C[6]{10,10,10,10,10,10};
    double M[6]{0.2,0.2,0.2,0.2,0.2,0.2};
    double deltaT{0.001};

	const CONTINUEMOVE_PARAM * pCMP = static_cast<const CONTINUEMOVE_PARAM *>(pParam);

	static CM_LAST_PARAM pCMLP;

    Aris::DynKer::s_pe2pe("313",pCMP->beginBodyPE,"213",beginBodyPE213);

    if(isContinue==true)
	{
		rt_printf("gait continuing\n");
		for (int i=0;i<6;i++)
		{
			if (moveDir[2*i+0]==true && moveDir[2*i+1]==false)
				F[i]=1;
			else if (moveDir[2*i+0]==false && moveDir[2*i+1]==true)
				F[i]=-1;
			else if (moveDir[2*i+0]==false && moveDir[2*i+1]==false)
				F[i]=0;

			bodyAcc[i]=(F[i]-C[i]*pCMLP.bodyVel_last[i])/M[i];
			bodyVel[i]=pCMLP.bodyVel_last[i]+bodyAcc[i]*deltaT;
			bodyPE[i]=pCMLP.bodyPE_last[i]+bodyVel[i]*deltaT;

			realBodyPE[i]=beginBodyPE213[i]+bodyPE[i];
		}

		double pBody[6];
		Aris::DynKer::s_pe2pe("213",realBodyPE,"313",pBody);

		pRobot->SetPee(pParam->beginPee,pBody);

		memcpy(pCMLP.bodyPE_last,bodyPE,sizeof(double)*6);
		memcpy(pCMLP.bodyVel_last,bodyVel,sizeof(double)*6);
		rt_printf("%d,%f,%f,%f,%f,%f,%f\n",pCMP->count,beginBodyPE213[0],beginBodyPE213[1],beginBodyPE213[2],beginBodyPE213[3],beginBodyPE213[4],beginBodyPE213[5]);
		rt_printf("%d,%f,%f,%f,%f,%f,%f\n",pCMP->count,bodyPE[0],bodyPE[1],bodyPE[2],bodyPE[3],bodyPE[4],bodyPE[5]);
	}
    else if(isContinue==false)
	{
		rt_printf("gait stopping\n");

		memcpy(F,zeros,sizeof(double)*6);
		for (int i=0;i<6;i++)
		{
			bodyAcc[i]=(F[i]-C[i]*pCMLP.bodyVel_last[i])/M[i];
			bodyVel[i]=pCMLP.bodyVel_last[i]+bodyAcc[i]*deltaT;
			bodyPE[i]=pCMLP.bodyPE_last[i]+bodyVel[i]*deltaT;

			realBodyPE[i]=beginBodyPE213[i]+bodyPE[i];
		}

		double pBody[6];
		Aris::DynKer::s_pe2pe("213",realBodyPE,"313",pBody);
		pRobot->SetPee(pParam->beginPee,pBody);

		memcpy(pCMLP.bodyPE_last,bodyPE,sizeof(double)*6);
		memcpy(pCMLP.bodyVel_last,bodyVel,sizeof(double)*6);

		//rt_printf("%d,%f,%f,%f,%f,%f,%f\n",pCMP->count,beginBodyPE213[0],beginBodyPE213[1],beginBodyPE213[2],beginBodyPE213[3],beginBodyPE213[4],beginBodyPE213[5]);
		//rt_printf("%d,%f,%f,%f,%f,%f,%f\n",pCMP->count,bodyPE[0],bodyPE[1],bodyPE[2],bodyPE[3],bodyPE[4],bodyPE[5]);

		if ( fabs(bodyVel[0])<1e-10 && fabs(bodyVel[1])<1e-10 && fabs(bodyVel[2])<1e-10 && fabs(bodyVel[3])<1e-10 && fabs(bodyVel[4])<1e-10 && fabs(bodyVel[5])<1e-10)
			return 0;
	}
}
