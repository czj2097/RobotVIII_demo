#include "Move_Gait.h"
#include "rtdk.h"

PIPE<MOVES_PARAM> move2Pipe(11,true);
PIPE<ForceTask::OPENDOOR_PARAM> openDoorPipe(16,true);

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

    move2Pipe.SendToNRT(*pMP);

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

std::atomic_bool isStoppingCWF;

Aris::Core::MSG parseCWF(const std::string &cmd, const std::map<std::string, std::string> &params)
{
    Robots::WALK_PARAM  param;

    for (auto &i : params)
    {
        if (i.first == "totalCount")
        {
            param.totalCount = std::stoi(i.second);
        }
        else if (i.first == "walkDirection")
        {
            param.walkDirection = std::stoi(i.second);
        }
        else if (i.first == "upDirection")
        {
            param.upDirection = std::stoi(i.second);
        }
        else if (i.first == "distance")
        {
            param.d = std::stod(i.second);
        }
        else if (i.first == "height")
        {
            param.h = std::stod(i.second);
        }
    }

    isStoppingCWF=false;

    Aris::Core::MSG msg;

    msg.CopyStruct(param);

    return msg;
}

Aris::Core::MSG parseCWFStop(const std::string &cmd, const std::map<std::string, std::string> &params)
{
    isStoppingCWF = true;

    Aris::Core::MSG msg;

    return msg;
}

int continuousWalkWithForce(Robots::ROBOT_BASE * pRobot, const Robots::GAIT_PARAM_BASE * pParam)
{
    const Robots::WALK_PARAM *pCWFP = static_cast<const Robots::WALK_PARAM *>(pParam);

    static bool isWalking=false;
    static int walkBeginCount{0};
    //static double walkBeginPee[18]{0};
    //static double walkBeginBodyPE[6]{0};
    static double forceOffsetSum[3]{0};

    double forceOffsetAvg[3]{0};
    double realForceData[3]{0};
    const double forceThreshold[3]{40,40,40};//力传感器的触发阈值,单位N或Nm
    const double forceRatio{1};

    //力传感器手动清零
    if (pCWFP->count<100)
    {
        if(pCWFP->count==0)
        {
            for(int i=0;i<3;i++)
            {
                forceOffsetSum[i]=0;
            }
        }
        forceOffsetSum[0]+=pCWFP->pForceData->at(0).Fx;
        forceOffsetSum[1]+=pCWFP->pForceData->at(0).Fy;
        forceOffsetSum[2]+=pCWFP->pForceData->at(0).Mz;
    }
    else
    {
        for(int i=0;i<3;i++)
        {
            forceOffsetAvg[i]=forceOffsetSum[i]/100;
        }
        if(pCWFP->count==100)
        {
            rt_printf("forceOffsetAvg: %f %f %f\n",forceOffsetAvg[0],forceOffsetAvg[1],forceOffsetAvg[2]);
        }
        realForceData[0]=(pCWFP->pForceData->at(0).Fx-forceOffsetAvg[0])/forceRatio;
        realForceData[1]=(pCWFP->pForceData->at(0).Fy-forceOffsetAvg[1])/forceRatio;
        realForceData[2]=(pCWFP->pForceData->at(0).Mz-forceOffsetAvg[2])/forceRatio;

        static Robots::WALK_PARAM realParam = *pCWFP;

        if(!isWalking)
        {
            WALK_DIRECTION walkDir=forceJudge(realForceData, forceThreshold);
            if(walkDir!=STOP)
            {
                switch (walkDir)
                {
                case FORWARD:
                    realParam.d=pCWFP->d;
                    realParam.alpha=0;
                    realParam.beta=0;
                    rt_printf("Walking Forward\n");
                    break;
                case BACKWARD:
                    realParam.d=-1*pCWFP->d;
                    realParam.alpha=0;
                    realParam.beta=0;
                    rt_printf("Walking Backward\n");
                    break;
                case LEFTWARD:
                    realParam.d=pCWFP->d/2;
                    realParam.alpha=PI/2;
                    realParam.beta=0;
                    rt_printf("Walking Leftward\n");
                    break;
                case RIGHTWARD:
                    realParam.d=pCWFP->d/2;
                    realParam.alpha=-PI/2;
                    realParam.beta=0;
                    rt_printf("Walking Rightward\n");
                    break;
                case TURNLEFT:
                    realParam.d=0;
                    realParam.alpha=0;
                    realParam.beta=PI/12;
                    rt_printf("Turning Left\n");
                    break;
                case TURNRIGHT:
                    realParam.d=0;
                    realParam.alpha=0;
                    realParam.beta=-PI/12;
                    rt_printf("Turning Right\n");
                    break;
                case FAST_TURNLEFT:
                    realParam.d=0;
                    realParam.alpha=0;
                    realParam.beta=PI/6;
                    rt_printf("Fast Turning Left\n");
                    break;
                case FAST_TURNRIGHT:
                    realParam.d=0;
                    realParam.alpha=0;
                    realParam.beta=-PI/6;
                    rt_printf("Fast Turning Right\n");
                    break;
                default:
                    break;
                }
                isWalking=true;
                walkBeginCount=pCWFP->count;
                pRobot->GetPee(realParam.beginPee);
                pRobot->GetBodyPe(realParam.beginBodyPE);

                rt_printf("realForceData: %f %f %f\n",realForceData[0],realForceData[1],realForceData[2]);
                rt_printf("beginBodyPE: %f %f %f\n",realParam.beginBodyPE[0],realParam.beginBodyPE[1],realParam.beginBodyPE[2]);
            }
        }
        else
        {
            realParam.count=pCWFP->count-walkBeginCount;
            int ret=Robots::walk(pRobot, &realParam);
            if(ret==0)
            {
                rt_printf("Finish One Walking Step\n");
                isWalking=false;
            }
        }
    }

    if(isStoppingCWF && (!isWalking))
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

WALK_DIRECTION forceJudge(const double *force, const double *threshold)
{
    WALK_DIRECTION walkDir{STOP};
    if(std::fabs(force[2]) > threshold[2])
    {
        if(force[2] < -2*threshold[2])
            walkDir=FAST_TURNRIGHT;
        else if(force[2] < -threshold[2])
            walkDir=TURNRIGHT;
        else if(force[2] > 2*threshold[2])
            walkDir=FAST_TURNLEFT;
        else if(force[2] > threshold[2])
            walkDir=TURNLEFT;
    }
    else if(std::fabs(std::fabs(force[0]) - threshold[0]) > std::fabs(std::fabs(force[1]) - threshold[1]))
    {
        if(force[0] < -threshold[0])
            walkDir=FORWARD;
        else if(force[0] > threshold[0])
            walkDir=BACKWARD;
    }
    else
    {
        if(force[1] < -threshold[1])
            walkDir=LEFTWARD;
        else if(force[1] > threshold[1])
            walkDir=RIGHTWARD;
    }
    return walkDir;
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



std::atomic_bool isForce;
std::atomic_bool isContinue;
std::atomic_int moveDir[6];
std::atomic_bool isPull;
std::atomic_bool isConfirm;

Aris::Core::MSG ForceTask::parseContinueMoveBegin(const std::string &cmd, const map<std::string, std::string> &params)
{
	CONTINUEMOVE_PARAM param;

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
			return MSG{};
		}
	}

	isContinue=true;
	isForce=false;

	Aris::Core::MSG msg;
	msg.CopyStruct(param);

	std::cout<<"finished parse"<<std::endl;

	return msg;
}

Aris::Core::MSG ForceTask::parseContinueMoveJudge(const std::string &cmd, const map<std::string, std::string> &params)
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
			return MSG{};
		}
	}

	Aris::Core::MSG msg;
	//msg.CopyStruct(param);

	std::cout<<"finished parse"<<std::endl;

	return msg;
}

/*****C & forceRatio must be adjusted when used on different robots*****/
int ForceTask::continueMove(Robots::ROBOT_BASE * pRobot, const Robots::GAIT_PARAM_BASE * pParam)
{
    double bodyVel[6];
    double bodyAcc[6];
    double bodyPm[4][4];
    double deltaPE[6];
    double deltaPm[4][4];
    double realPE[6];
    double realPm[4][4];

    double Fbody[6]{0,0,0,0,0,0};
    double C[6]{30,30,30,30,30,30};
    double M[6]{1,1,1,1,1,1};
    double deltaT{0.001};
    double forceRange[6]{30,30,30,20,20,20};
    double forceRatio{1};//1 on RobotIII, 1000 on RobotVIII & single motor

	const CONTINUEMOVE_PARAM * pCMP = static_cast<const CONTINUEMOVE_PARAM *>(pParam);

	static ForceTask::CM_RECORD_PARAM CMRP;

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
    		if(pCMP->count<100)
    		{
    			//initialize
                if (pCMP->count==0)
                {
                	for(int i=0;i<6;i++)
                	{
                		CMRP.forceSum[i]=0;
                		CMRP.bodyVel_last[i]=0;
                	}
                }

                CMRP.forceSum[0]+=pCMP->pForceData->at(0).Fx;
                CMRP.forceSum[1]+=pCMP->pForceData->at(0).Fy;
                CMRP.forceSum[2]+=pCMP->pForceData->at(0).Fz;
                CMRP.forceSum[3]+=pCMP->pForceData->at(0).Mx;
                CMRP.forceSum[4]+=pCMP->pForceData->at(0).My;
                CMRP.forceSum[5]+=pCMP->pForceData->at(0).Mz;
        	}
        	else if(pCMP->count==100)
        	{
        		for(int i=0;i<6;i++)
        		{
        			CMRP.forceAvg[i]=CMRP.forceSum[i]/100;
        		}
        	}
        	else
        	{
        		CMRP.force[0]=(pParam->pForceData->at(0).Fx-CMRP.forceAvg[0])/forceRatio;
        		CMRP.force[1]=(pParam->pForceData->at(0).Fy-CMRP.forceAvg[1])/forceRatio;
        		CMRP.force[2]=(pParam->pForceData->at(0).Fz-CMRP.forceAvg[2])/forceRatio;
        		CMRP.force[3]=(pParam->pForceData->at(0).Mx-CMRP.forceAvg[3])/forceRatio;
        		CMRP.force[4]=(pParam->pForceData->at(0).My-CMRP.forceAvg[4])/forceRatio;
        		CMRP.force[5]=(pParam->pForceData->at(0).Mz-CMRP.forceAvg[5])/forceRatio;

				int force2robot[6]={2,0,1,5,4,3};
				int num1;
				int num2;
				double fmax{0};
				double mmax{0};
				for (int i=0;i<3;i++)
				{
					if (fabs(CMRP.force[i])>fabs(fmax))
					{
						fmax=CMRP.force[i];
						num1=i;
					}
					if (fabs(CMRP.force[i+3])>fabs(mmax))
					{
						mmax=CMRP.force[i+3];
						num2=i+3;
					}
				}

				if(CMRP.force[num2]>forceRange[num2])
				{
					Fbody[force2robot[num2]]=1;
				}
				else if(CMRP.force[num2]<-forceRange[num2])
				{
					Fbody[force2robot[num2]]=-1;
				}
				else
				{
					if(CMRP.force[num1]>forceRange[num1])
					{
						Fbody[force2robot[num1]]=1;
					}
					else if(CMRP.force[num1]<-forceRange[num1])
					{
						Fbody[force2robot[num1]]=-1;
					}
				}
        	}
    	}

		for (int i=0;i<6;i++)
		{
			bodyAcc[i]=(Fbody[i]-C[i]*CMRP.bodyVel_last[i])/M[i];
			bodyVel[i]=CMRP.bodyVel_last[i]+bodyAcc[i]*deltaT;
			deltaPE[i]=bodyVel[i]*deltaT;
		}

		pRobot->GetBodyPm(*bodyPm);
		double pBody[6];
		Aris::DynKer::s_pe2pm(deltaPE,*deltaPm,"213");
		Aris::DynKer::s_pm_dot_pm(*bodyPm,*deltaPm,*realPm);
		Aris::DynKer::s_pm2pe(*realPm,realPE,"313");
		double nowPee[18];
		pRobot->GetPee(nowPee);
		pRobot->SetPee(nowPee,realPE);

		if(pCMP->count%100==0)
		{
			rt_printf("bodyPm:%f,%f,%f,%f\n%f,%f,%f,%f\n%f,%f,%f,%f\n%f,%f,%f,%f\n",bodyPm[0][0],bodyPm[0][1],bodyPm[0][2],bodyPm[0][3],bodyPm[1][0],bodyPm[1][1],bodyPm[1][2],bodyPm[1][3],bodyPm[2][0],bodyPm[2][1],bodyPm[2][2],bodyPm[2][3],bodyPm[3][0],bodyPm[3][1],bodyPm[3][2],bodyPm[3][3]);
			rt_printf("Fbody:%f,%f,%f,%f,%f,%f\n",Fbody[0],Fbody[1],Fbody[2],Fbody[3],Fbody[4],Fbody[5]);
			rt_printf("realPE:%f,%f,%f,%f,%f,%f\n\n",realPE[0],realPE[1],realPE[2],realPE[3],realPE[4],realPE[5]);
		}

		memcpy(CMRP.bodyVel_last,bodyVel,sizeof(double)*6);
		return 1;
	}

    else
	{
		//rt_printf("gait stopping\n");
		for (int i=0;i<6;i++)
		{
			bodyAcc[i]=(Fbody[i]-C[i]*CMRP.bodyVel_last[i])/M[i];
			bodyVel[i]=CMRP.bodyVel_last[i]+bodyAcc[i]*deltaT;
			deltaPE[i]=bodyVel[i]*deltaT;
		}

		pRobot->GetBodyPm(*bodyPm);
		double pBody[6];
		Aris::DynKer::s_pe2pm(deltaPE,*deltaPm,"213");
		Aris::DynKer::s_pm_dot_pm(*bodyPm,*deltaPm,*realPm);
		Aris::DynKer::s_pm2pe(*realPm,realPE,"313");
		double nowPee[18];
		pRobot->GetPee(nowPee);
		pRobot->SetPee(nowPee,realPE);

		memcpy(CMRP.bodyVel_last,bodyVel,sizeof(double)*6);

		if ( fabs(bodyVel[0])<1e-10 && fabs(bodyVel[1])<1e-10 && fabs(bodyVel[2])<1e-10 && fabs(bodyVel[3])<1e-10 && fabs(bodyVel[4])<1e-10 && fabs(bodyVel[5])<1e-10)
			return 0;
		else
			return 1;
	}
}

Aris::Core::MSG ForceTask::parseOpenDoorBegin(const std::string &cmd, const map<std::string, std::string> &params)
{
	CONTINUEMOVE_PARAM param;

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
			return MSG{};
		}
	}
    isPull=true;
	isContinue=true;
	isConfirm=false;

	Aris::Core::MSG msg;
	msg.CopyStruct(param);

	std::cout<<"finished parse"<<std::endl;

	return msg;
}

Aris::Core::MSG ForceTask::parseOpenDoorJudge(const std::string &cmd, const map<std::string, std::string> &params)
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
			return MSG{};
		}
	}

	Aris::Core::MSG msg;
	//msg.CopyStruct(param);

	std::cout<<"finished parse"<<std::endl;

	return msg;
}

void ForceTask::inv3(double * matrix,double * invmatrix)
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

void ForceTask::crossMultiply(double * vector_in1, double *vector_in2, double * vector_out)
{
	vector_out[0]=vector_in1[1]*vector_in2[2]-vector_in1[2]*vector_in2[1];
	vector_out[1]=vector_in1[2]*vector_in2[0]-vector_in1[0]*vector_in2[2];
	vector_out[2]=vector_in1[0]*vector_in2[1]-vector_in1[1]*vector_in2[0];
}

double ForceTask::dotMultiply(double *vector_in1, double *vector_in2)
{
	double sum{0};
	for(int i=0;i<3;i++)
	{
		sum+=vector_in1[i]*vector_in2[i];
	}
	return sum;
}

double ForceTask::norm(double * vector_in)
{
	return	sqrt(vector_in[0]*vector_in[0]+vector_in[1]*vector_in[1]+vector_in[2]*vector_in[2]);
}



//*****only for Robot VIII*****
int ForceTask::openDoor(Robots::ROBOT_BASE * pRobot, const Robots::GAIT_PARAM_BASE * pParam)
{
    //MoveState: PointLocation
    double beginBodyPm[4][4];
    double invBeginBodyPm[4][4];
    double location[3][3];
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
    double d0=0.42;//distance from the handle to the middle of the door
    double h0=0.05;//height from the start position to the walk through position

    //Force Control
    double Fbody[6]{0,0,0,0,0,0};
	double C[6]{50,50,50,50,50,50};
	double M[6]{1,1,1,1,1,1};
	double deltaT{0.001};
    double ForceRange[2]{10,90};
    double forceRatio{1000};//1 on RobotIII, 1000 on RobotVIII & single motor

    //Position Generetion From Force
    double bodyVel[6];
    double bodyAcc[6];
    double bodyPE[6];//delta PE every ms
    double bodyPm[4][4];//bodyPm to Ground every ms
    double deltaPE[6];
    double deltaPm[4][4];
    double realPm[4][4];
    double realPE[6];

	const CONTINUEMOVE_PARAM * pCMP = static_cast<const CONTINUEMOVE_PARAM *>(pParam);

	static ForceTask::OPENDOOR_PARAM ODP;
	ODP.count=pCMP->count;

    Aris::DynKer::s_pe2pm(pCMP->beginBodyPE,*beginBodyPm,"313");

    if(isContinue==true)
	{
    	if (pCMP->count<100)
    	{
            if (pCMP->count==0)
            {
            	//initialize
            	for(int i=0;i<6;i++)
            	{
            		ODP.forceSum[i]=0;
            		ODP.bodyVel_last[i]=0;
            	}
            	Aris::DynKer::s_pe2pe("313",pCMP->beginBodyPE,"213",ODP.bodyPE_last);
            	ODP.pauseFlag=false;
            	ODP.moveState_last=MoveState::None;

            	//start param of now2start in LocateAjust
    			pRobot->GetBodyPe(ODP.startPE);
            }

            ODP.moveState=MoveState::None;

            ODP.forceSum[0]+=pCMP->pForceData->at(0).Fx;
            ODP.forceSum[1]+=pCMP->pForceData->at(0).Fy;
            ODP.forceSum[2]+=pCMP->pForceData->at(0).Fz;
            ODP.forceSum[3]+=pCMP->pForceData->at(0).Mx;
            ODP.forceSum[4]+=pCMP->pForceData->at(0).My;
            ODP.forceSum[5]+=pCMP->pForceData->at(0).Mz;

    	}
    	else if(pCMP->count==100)
    	{
    		for(int i=0;i<6;i++)
    		{
    			ODP.forceAvg[i]=ODP.forceSum[i]/100;
    		}
    	}
    	else
    	{
    		ODP.force[0]=(pParam->pForceData->at(0).Fx-ODP.forceAvg[0])/forceRatio;
    		ODP.force[1]=(pParam->pForceData->at(0).Fy-ODP.forceAvg[1])/forceRatio;
    		ODP.force[2]=(pParam->pForceData->at(0).Fz-ODP.forceAvg[2])/forceRatio;
    		ODP.force[3]=(pParam->pForceData->at(0).Mx-ODP.forceAvg[3])/forceRatio;
    		ODP.force[4]=(pParam->pForceData->at(0).My-ODP.forceAvg[4])/forceRatio;
    		ODP.force[5]=(pParam->pForceData->at(0).Mz-ODP.forceAvg[5])/forceRatio;
    	}

    	switch(ODP.moveState)
    	{
    	case MoveState::None:
            for (int i=0;i<6;i++)
            {
            	Fbody[i]=moveDir[i];
            }

            if (fabs(ODP.force[0])>ForceRange[0] && pCMP->count>200 && ODP.pauseFlag==false)
            {
                ODP.countIter=pCMP->count;
                pRobot->GetBodyPe(ODP.pointLocation1);
                memcpy(*location,ODP.pointLocation1,sizeof(double)*3);
                ODP.moveState=MoveState::PointLocate1;
            }

    		break;

    	case MoveState::PointLocate1:
    		if (pCMP->count-ODP.countIter<5000)
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

    		if (fabs(ODP.force[0])>ForceRange[0] && (pCMP->count-ODP.countIter)>5000)
			{
				ODP.countIter=pCMP->count;
				pRobot->GetBodyPe(ODP.pointLocation2);
				memcpy(*location+3,ODP.pointLocation2,sizeof(double)*3);
				ODP.moveState=MoveState::PointLocate2;
			}

    		break;

    	case MoveState::PointLocate2:
    		if (pCMP->count-ODP.countIter<5000)
			{
				Fbody[2]=1;
				Fbody[1]=1;
			}
			else
			{
				Fbody[2]=-1;
			}

			if (fabs(ODP.force[0])>ForceRange[0] && (pCMP->count-ODP.countIter)>5000)
			{
				ODP.countIter=pCMP->count+1;
				pRobot->GetBodyPe(ODP.pointLocation3);
				memcpy(*location+6,ODP.pointLocation3,sizeof(double)*3);
				ODP.moveState=MoveState::LocateAjust;

				//calculate the plane of the door. ax+by+cz=1,(a,b,c) is the vertical vector of the plane
				ForceTask::inv3(*location,*invLocation);
				Aris::DynKer::s_dgemm(3,1,3,1,*invLocation,3,planeConst,1,1,planeVertical,1);
				Aris::DynKer::s_inv_pm(*beginBodyPm,*invBeginBodyPm);//have not rotate, beginBodyPm is right here
				Aris::DynKer::s_pm_dot_v3(*invBeginBodyPm,planeVertical,planeVerticalInB);
				ODP.planeYPR[0]=atan(planeVerticalInB[0]/planeVerticalInB[2]);
				ODP.planeYPR[1]=-asin(planeVerticalInB[1]/ForceTask::norm(planeVerticalInB));
				ODP.planeYPR[2]=0;

				//Set now param of now2start in LocateAjust
				pRobot->GetBodyPe(ODP.nowPE);
				pRobot->GetPee(ODP.nowPee);

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
    		if(pCMP->count-ODP.countIter<ODP.now2StartCount)
    		{
				for (int i=0;i<6;i++)
				{
					now2StartPE[i]=ODP.nowPE[i]+(ODP.startPE[i]-ODP.nowPE[i])/2*(1-cos((pCMP->count-ODP.countIter)*PI/ODP.now2StartCount));
				}

				pRobot->SetPee(ODP.nowPee,now2StartPE);

				if(pCMP->count-ODP.countIter+1==ODP.now2StartCount)
				{
					pRobot->GetPee(ODP.walkParam.beginPee);
					pRobot->GetBodyPe(ODP.walkParam.beginBodyPE);
				}
    		}
    		//2.rotate
    		else
    		{
				ODP.walkParam.count=pCMP->count-ODP.countIter-ODP.now2StartCount;
				ODP.ret=Robots::walk(pRobot,& ODP.walkParam);

				if(ODP.ret==0)
				{
					ODP.countIter=pCMP->count;
					ODP.moveState=MoveState::Forward;
				}
    		}

    		pRobot->GetBodyPm(*bodyPm);
    		pRobot->GetBodyPe(bodyPE,"213");
			for(int i=0;i<6;i++)
			{
				bodyVel[i]=bodyPE[i]-ODP.bodyPE_last[i];
			}

    		break;

    	case MoveState::Forward:
    		Fbody[2]=-1;

    		if(pCMP->count==ODP.countIter+1)
    		{
    			pRobot->GetBodyPe(ODP.startPE);//Used in PrePush for now2Start
    		}

    		if(fabs(ODP.force[0])>ForceRange[0])
    		{
    			ODP.countIter=pCMP->count;
    			ODP.moveState=MoveState::Backward;
    			pRobot->GetBodyPe(touchPE);//For calculation of DoorLocation

    			for (int i=0;i<3;i++)
    			{
    				ODP.vector0[i]=touchPE[i]-ODP.startPE[i];
    			}
    		}

    		break;

    	case MoveState::Backward:
    		Fbody[2]=1;

    		if (pCMP->count-ODP.countIter>750)
    		{
    			ODP.countIter=pCMP->count;
    			pRobot->GetPee(ODP.endPeeInB,"B");
    			pRobot->GetBodyPe(ODP.beginPE);//For calculation of DoorLocation
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
				if (fabs(ODP.force[1])>ForceRange[0])
				{
					ODP.countIter=pCMP->count;
					pRobot->GetBodyPe(ODP.handlePE);//for PushState: leftWalk
					ODP.moveState=MoveState::Rightward;

					for (int i=0;i<3;i++)
					{
						ODP.vector1[i]=ODP.handlePE[i]-ODP.beginPE[i];
					}
				}
				else if (fabs(ODP.force[1])<=ForceRange[0] && pCMP->count-ODP.countIter>7500)
				{
					ODP.countIter=pCMP->count+1;
					pRobot->GetPee(ODP.startPeeInB,"B");
					pRobot->GetBodyPe(ODP.nowPE);
					ODP.moveState=MoveState::Follow;
				}
			}
			else//Push
			{
				if (pCMP->count-ODP.countIter>3000)
				{
					ODP.countIter=pCMP->count;
					ODP.moveState=MoveState::Downward;
					ODP.downwardCount=(int)(PI/2*2500);
				    ODP.downwardFlag = false;
				    pRobot->GetBodyPe(ODP.beginPE);//For calculation of DoorLocation
				}
			}

    		break;

    	case MoveState::Rightward:

    		Fbody[0]=1;

    		if(isPull==false)
    		{
    			if (fabs(ODP.force[1])>ForceRange[0])
    			{
    				ODP.countIter=pCMP->count;
    				pRobot->GetBodyPe(ODP.handlePE);
    				ODP.moveState=MoveState::Leftward;
    			}
    			else if (fabs(ODP.force[1])<ForceRange[0] && pCMP->count-ODP.countIter>7500)
				{
					ODP.countIter=pCMP->count+1;
					pRobot->GetPee(ODP.startPeeInB,"B");
					pRobot->GetBodyPe(ODP.nowPE);
					ODP.moveState=MoveState::Follow;
				}
    		}
    		else//Pull
    		{
    			if (pCMP->count-ODP.countIter>3000)
    			{
    				ODP.countIter=pCMP->count;
    				ODP.moveState=MoveState::Downward;
    				ODP.downwardCount=(int)(PI/2*2500);
    				ODP.downwardFlag = false;
    			}
    		}

    		break;

    	case MoveState::Follow:
    		if(pCMP->count-ODP.countIter<ODP.followCount)
			{
    			for(int i=0;i<3;i++)
    			{
    				//leg 0,2,4
    				followPeeInB[6*i]=ODP.startPeeInB[6*i]+(ODP.endPeeInB[6*i]-ODP.startPeeInB[6*i])/2*(1-cos((pCMP->count-ODP.countIter)*PI/ODP.followCount));
    				followPeeInB[6*i+1]=ODP.startPeeInB[6*i+1]+0.05*(1-cos((pCMP->count-ODP.countIter)*2*PI/ODP.followCount));
					followPeeInB[6*i+2]=ODP.startPeeInB[6*i+2];
					//leg 1,3,5
					followPeeInB[6*i+3]=ODP.startPeeInB[6*i+3];
					followPeeInB[6*i+4]=ODP.startPeeInB[6*i+4];
					followPeeInB[6*i+5]=ODP.startPeeInB[6*i+5];
    			}
			}
    		else if(pCMP->count-ODP.countIter<2*ODP.followCount)
    		{
    			for(int i=0;i<3;i++)
				{
					//leg 1,3,5
					followPeeInB[6*i+3]=ODP.startPeeInB[6*i+3]+(ODP.endPeeInB[6*i+3]-ODP.startPeeInB[6*i+3])/2*(1-cos((pCMP->count-ODP.countIter-ODP.followCount)*PI/ODP.followCount));
					followPeeInB[6*i+4]=ODP.startPeeInB[6*i+4]+0.05*(1-cos((pCMP->count-ODP.countIter-ODP.followCount)*2*PI/ODP.followCount));
					followPeeInB[6*i+5]=ODP.startPeeInB[6*i+5];
					//leg 0,2,4
					followPeeInB[6*i]=ODP.endPeeInB[6*i];
					followPeeInB[6*i+1]=ODP.endPeeInB[6*i+1];
					followPeeInB[6*i+2]=ODP.endPeeInB[6*i+2];
				}
    		}
    		pRobot->SetPee(followPeeInB,ODP.nowPE,"B");

    		if(pCMP->count-ODP.countIter+1==2*ODP.followCount)
    		{
    			ODP.countIter=pCMP->count;
    			ODP.moveState=MoveState::Forward;
    		}

    		pRobot->GetBodyPm(*bodyPm);
    		pRobot->GetBodyPe(bodyPE,"213");
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
    				Fbody[1]=-cos((pCMP->count-ODP.countIter)*PI/ODP.downwardCount/2);
    				Fbody[0]=sin((pCMP->count-ODP.countIter)*PI/ODP.downwardCount/2);
    			}
    			else
    			{
    				Fbody[1]=-cos((pCMP->count-ODP.countIter)*PI/ODP.downwardCount/2);
    				Fbody[0]=-sin((pCMP->count-ODP.countIter)*PI/ODP.downwardCount/2);
    			}
    			Fbody[2]=(fabs(ODP.force[0])-50)/50;
            }

            if(ODP.downwardFlag==false && fabs(ODP.force[2])>ForceRange[0])
            {
            	ODP.downwardFlag = true;
            	ODP.countIter=pCMP->count;
            	pRobot->GetBodyPe(touchPE);//For calculation of DoorLocation

            	for (int i=0;i<3;i++)
				{
					ODP.vector2[i]=touchPE[i]-ODP.beginPE[i];
				}
            	pRobot->GetBodyPm(* ODP.nowPm);

            	ODP.handleLocation[2]=-ForceTask::norm(ODP.vector0);
            	ODP.handleLocation[0]=ForceTask::norm(ODP.vector1);
            	ODP.handleLocation[1]=-ForceTask::norm(ODP.vector2);

            	rt_printf("handleLocation:%f,%f,%f\n",ODP.handleLocation[0],ODP.handleLocation[1],ODP.handleLocation[2]);
            }

    		if (fabs(ODP.force[1])>ForceRange[1] && pCMP->count-ODP.countIter>ODP.downwardCount/2)
    		{
    			ODP.countIter=pCMP->count;
    			if(isPull==true)
    				ODP.moveState=MoveState::Pullhandle;
    			else
    				ODP.moveState=MoveState::Pushhandle;
    		}

    		break;

    	case MoveState::Pullhandle:

    		Fbody[2]=1;

    		if (pCMP->count-ODP.countIter>2500)
    			ODP.moveState=MoveState::Upward;

    		break;

    	case MoveState::Upward:

			Fbody[1]=1;
			if(isPull==true)
			{
				if (pCMP->count-ODP.countIter>5000)
				{
					isContinue=false;//Stop here when Pull
				}
			}

			break;

    	case MoveState::Pushhandle:

    		Fbody[2]=-1;

    		if (pCMP->count-ODP.countIter>2500)
            {
                ODP.moveState=MoveState::PrePush;
            }
    		break;

        case MoveState::PrePush:

            if ( fabs(bodyVel[0])<1e-10 && fabs(bodyVel[1])<1e-10 && fabs(bodyVel[2])<1e-10 && fabs(bodyVel[3])<1e-10 && fabs(bodyVel[4])<1e-10 && fabs(bodyVel[5])<1e-10)
            {
                ODP.countIter=pCMP->count+1;
                ODP.moveState=MoveState::Push;
                ODP.pushState=PushState::now2Start;

                //now2start param
                pRobot->GetBodyPe(ODP.nowPE);
				pRobot->GetPee(ODP.nowPee);
				for(int i=0;i<3;i++)
				{
					ODP.now2startDistance[i]=ODP.startPE[i]-ODP.nowPE[i]; //startPE recorded at the last forward
					ODP.handle2startDistance[i]=ODP.startPE[i]-ODP.handlePE[i];
				}
				Aris::DynKer::s_pm_dot_v3(*ODP.nowPm,xBodyInB,ODP.xNowInG);
				Aris::DynKer::s_pm_dot_v3(*ODP.nowPm,yBodyInB,ODP.yNowInG);

				ODP.now2startDistanceModified[0]=ForceTask::dotMultiply(ODP.now2startDistance,ODP.xNowInG);
				ODP.now2startDistanceModified[1]=ForceTask::dotMultiply(ODP.now2startDistance,ODP.yNowInG)+h0;
				ODP.now2startDistanceModified[2]=0;

				Aris::DynKer::s_inv_pm_dot_v3(*ODP.nowPm,ODP.now2startDistanceModified,ODP.now2startDistanceReal);
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
                    now2StartPE[i]=ODP.nowPE[i]+(ODP.now2startDistanceReal[i])/2*(1-cos((pCMP->count-ODP.countIter)*PI/ODP.now2StartCount));
                }

                pRobot->SetPee(ODP.nowPee,now2StartPE);

                if(pCMP->count-ODP.countIter+1==ODP.now2StartCount)
                {
                    if(isConfirm==true)
                    {
                        ODP.pushState=PushState::leftWalk;
                        ODP.countIter=pCMP->count+1;

                        ODP.walkParam.n=2;
                        ODP.walkParam.alpha=PI/2;
                        ODP.walkParam.beta=0;
                        ODP.walkParam.totalCount=2000;
                        ODP.walkParam.d=(d0-fabs(ForceTask::dotMultiply(ODP.handle2startDistance,ODP.xNowInG)))/3*2;
                        pRobot->GetPee(ODP.walkParam.beginPee);
                        pRobot->GetBodyPe(ODP.walkParam.beginBodyPE);
                    }
                    else//pause
                    {
                        pRobot->SetPee(ODP.nowPee,ODP.startPE);
                    }
                }

                break;

            case PushState::leftWalk:
                //2.Move rightward to align with the door;
                ODP.walkParam.count=pCMP->count-ODP.countIter;

                ODP.ret=Robots::walk(pRobot,& ODP.walkParam);

                if(ODP.ret==0)
                {
                    if(isConfirm==true)
                    {
                        ODP.pushState=PushState::forwardWalk;
                        ODP.countIter=pCMP->count+1;

                        ODP.walkParam.totalCount=2000;
                        ODP.walkParam.n=4;
                        ODP.walkParam.alpha=0;
                        ODP.walkParam.beta=0;
                        ODP.walkParam.d=0.5;
                        pRobot->GetPee(ODP.walkParam.beginPee);
                        pRobot->GetBodyPe(ODP.walkParam.beginBodyPE);
                    }
                    else
                    {
                        pRobot->GetPee(pushPee);
                        pRobot->GetBodyPe(pushBodyPE313);
                        pRobot->SetPee(pushPee,pushBodyPE313);
                    }
                }

                break;

            case PushState::forwardWalk:
                //3.Move through the door
                ODP.walkParam.count=pCMP->count-ODP.countIter;

                ODP.ret=Robots::walk(pRobot,& ODP.walkParam);

                if(ODP.ret==0)
                {
                    if(isConfirm==true)
                    {
						return 0;
                    }
                    else
                    {
                        pRobot->GetPee(pushPee);
                        pRobot->GetBodyPe(pushBodyPE313);
                        pRobot->SetPee(pushPee,pushBodyPE313);
                    }
                }

                break;

            default:
                break;
            }

            pRobot->GetBodyPm(*bodyPm);
            pRobot->GetBodyPe(bodyPE,"213");
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

            pRobot->GetBodyPm(*bodyPm);
            pRobot->GetBodyPe(bodyPE);
			double pBody[6];
			Aris::DynKer::s_pe2pm(deltaPE,*deltaPm,"213");
			Aris::DynKer::s_pm_dot_pm(*bodyPm,*deltaPm,*realPm);
			Aris::DynKer::s_pm2pe(*realPm,realPE,"313");
			double nowPee[18];
			pRobot->GetPee(nowPee);
			pRobot->SetPee(nowPee,realPE);
        }

        rt_printf("moveState:%d,force:%f,%f,%f\n",ODP.moveState,ODP.force[0],ODP.force[1],ODP.force[2]);

        Aris::DynKer::s_pm_dot_pnt(*bodyPm,ODP.toolInR,ODP.toolInG);

        memcpy(ODP.bodyPE_last,bodyPE,sizeof(double)*6);
        memcpy(ODP.bodyVel_last,bodyVel,sizeof(double)*6);

		openDoorPipe.SendToNRT(ODP);

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

		pRobot->GetBodyPm(*bodyPm);
		pRobot->GetBodyPe(bodyPE);
		double pBody[6];
		Aris::DynKer::s_pe2pm(deltaPE,*deltaPm,"213");
		Aris::DynKer::s_pm_dot_pm(*bodyPm,*deltaPm,*realPm);
		Aris::DynKer::s_pm2pe(*realPm,realPE,"313");
		double nowPee[18];
		pRobot->GetPee(nowPee);
		pRobot->SetPee(nowPee,realPE);

		memcpy(ODP.bodyPE_last,bodyPE,sizeof(double)*6);
		memcpy(ODP.bodyVel_last,bodyVel,sizeof(double)*6);

		if ( fabs(bodyVel[0])<1e-10 && fabs(bodyVel[1])<1e-10 && fabs(bodyVel[2])<1e-10 && fabs(bodyVel[3])<1e-10 && fabs(bodyVel[4])<1e-10 && fabs(bodyVel[5])<1e-10)
			return 0;
		else
			return 1;
	}
}

void ForceTask::StartRecordData()
{
	openDoorThread = std::thread([&]()
	{
		struct OPENDOOR_PARAM param;
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
			fileGait << param.moveState << "  ";

			for (int i = 0; i < 6; i++)
			{
				fileGait << param.bodyPE_last[i] << "  ";
			}
			for (int i = 0; i < 3; i++)
			{
				fileGait << param.toolInG[i] << "  ";
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
