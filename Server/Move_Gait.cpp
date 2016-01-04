#include "Move_Gait.h"
#include "rtdk.h"

PIPE<MOVES_PARAM> move2Pipe(11,true);
PIPE<CM_LAST_PARAM> openDoorPipe(16,true);

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
std::atomic_bool moveDir[12];
std::atomic_bool isPull;
std::atomic_bool isConfirm;

Aris::Core::MSG parseContinueMoveBegin(const std::string &cmd, const map<std::string, std::string> &params)
{
	CONTINUEMOVE_PARAM param;

	for(auto &i:params)
	{
		if(i.first=="u")
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

	isContinue=true;
	isForce=true;

	Aris::Core::MSG msg;
	msg.CopyStruct(param);

	std::cout<<"finished parse"<<std::endl;

	return msg;
}

Aris::Core::MSG parseContinueMoveJudge(const std::string &cmd, const map<std::string, std::string> &params)
{
    //   cmd -u -v -w -r -p -y -iscontinue -ispull -isconfirm
	//CONTINUEMOVE_PARAM param;

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
            //param.move_direction[0]=stoi(i.second);

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
            //param.move_direction[1]=stoi(i.second);

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
            //param.move_direction[2]=stoi(i.second);

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
            //param.move_direction[3]=stoi(i.second);

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
            //param.move_direction[4]=stoi(i.second);

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
            //param.move_direction[5]=stoi(i.second);

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
	//msg.CopyStruct(param);

	std::cout<<"finished parse"<<std::endl;

	return msg;
}

int continueMove(Robots::ROBOT_BASE * pRobot, const Robots::GAIT_PARAM_BASE * pParam)
{
    double bodyVel[6];
    double bodyAcc[6];
    double bodyPE[6];

    double realBodyPE[6];
    double bodyPm[4][4];
    double deltaPm[4][4];
    double realPm[4][4];

    double zeros[6];

    double Fbody[6]{0,0,0,0,0,0};
    double C[6]{30,30,30,30,30,30};
    double M[6]{1,1,1,1,1,1};
    double deltaT{0.001};
    double forceRange[6]{30,30,30,20,20,20};
    double forceRatio{1};

	const CONTINUEMOVE_PARAM * pCMP = static_cast<const CONTINUEMOVE_PARAM *>(pParam);

	static CM_LAST_PARAM pCMLP;

    if(pCMP->count % 200 ==0 ){
        //rt_printf("Force x=%.3lf y=%.3lf z=%.3lf\n", pCMP->pForceData->at(0).Fx,pCMP->pForceData->at(0).Fy,pCMP->pForceData->at(0).Fz);
        //rt_printf("FoAct x=%.3lf y=%.3lf z=%.3lf\n", pCMLP.force[0],pCMLP.force[1], pCMLP.force[2]);
    }

    if(isContinue==true)
	{
    	//rt_printf("gait continuing\n");
    	if(isForce==false)
    	{
			for (int i=0;i<6;i++)
			{
				if (moveDir[2*i+0]==true && moveDir[2*i+1]==false)
					Fbody[i]=1;
				else if (moveDir[2*i+0]==false && moveDir[2*i+1]==true)
					Fbody[i]=-1;
				else if (moveDir[2*i+0]==false && moveDir[2*i+1]==false)
					Fbody[i]=0;
				else
					rt_printf("parse move diretion wrong!!!\n");
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
                		pCMLP.forceSum[i]=0;
                		pCMLP.bodyVel_last[i]=0;
                	}
                	Aris::DynKer::s_pe2pe("313",pCMP->beginBodyPE,"213",pCMLP.bodyPE_last);
                }

        		pCMLP.forceSum[0]+=pCMP->pForceData->at(0).Fx;
        		pCMLP.forceSum[1]+=pCMP->pForceData->at(0).Fy;
        		pCMLP.forceSum[2]+=pCMP->pForceData->at(0).Fz;
        		pCMLP.forceSum[3]+=pCMP->pForceData->at(0).Mx;
				pCMLP.forceSum[4]+=pCMP->pForceData->at(0).My;
				pCMLP.forceSum[5]+=pCMP->pForceData->at(0).Mz;
        	}
        	else if(pCMP->count==100)
        	{
        		for(int i=0;i<6;i++)
        		{
        			pCMLP.forceAvg[i]=pCMLP.forceSum[i]/100;
        		}
        	}
        	else
        	{
    			pCMLP.force[0]=(pParam->pForceData->at(0).Fx-pCMLP.forceAvg[0])/forceRatio;
    			pCMLP.force[1]=(pParam->pForceData->at(0).Fy-pCMLP.forceAvg[1])/forceRatio;
    			pCMLP.force[2]=(pParam->pForceData->at(0).Fz-pCMLP.forceAvg[2])/forceRatio;
    			pCMLP.force[3]=(pParam->pForceData->at(0).Mx-pCMLP.forceAvg[3])/forceRatio;
				pCMLP.force[4]=(pParam->pForceData->at(0).My-pCMLP.forceAvg[4])/forceRatio;
				pCMLP.force[5]=(pParam->pForceData->at(0).Mz-pCMLP.forceAvg[5])/forceRatio;

				int force2robot[6]={2,0,1,5,4,3};
				int num1;
				int num2;
				double fmax{0};
				double mmax{0};
				for (int i=0;i<3;i++)
				{
					if (fabs(pCMLP.force[i])>fabs(fmax))
					{
						fmax=pCMLP.force[i];
						num1=i;
					}
					if (fabs(pCMLP.force[i+3])>fabs(mmax))
					{
						mmax=pCMLP.force[i+3];
						num2=i+3;
					}
				}
				//rt_printf("The max force dirction is %d\n",num);

				if(pCMLP.force[num2]>forceRange[num2])
				{
					Fbody[force2robot[num2]]=1;
				}
				else if(pCMLP.force[num2]<-forceRange[num2])
				{
					Fbody[force2robot[num2]]=-1;
				}
				else
				{
					if(pCMLP.force[num1]>forceRange[num1])
					{
						Fbody[force2robot[num1]]=1;
					}
					else if(pCMLP.force[num1]<-forceRange[num1])
					{
						Fbody[force2robot[num1]]=-1;
					}
				}
        	}
    	}

		for (int i=0;i<6;i++)
		{
			// Set Position Limit
			/*
			if (fabs(pCMLP.bodyPE_last[i])>pCMLP.posLimit[i])
			{
				Fground[0]=0;
				rt_printf("position of direction %d is out of limit\n",i);
			}*/

			bodyAcc[i]=(Fbody[i]-C[i]*pCMLP.bodyVel_last[i])/M[i];
			bodyVel[i]=pCMLP.bodyVel_last[i]+bodyAcc[i]*deltaT;
			bodyPE[i]=bodyVel[i]*deltaT;
		}

		pRobot->GetBodyPm(*bodyPm);
		double pBody[6];
		Aris::DynKer::s_pe2pm(bodyPE,*deltaPm,"213");
		Aris::DynKer::s_pm_dot_pm(*bodyPm,*deltaPm,*realPm);
		Aris::DynKer::s_pm2pe(*realPm,realBodyPE,"313");
		double nowPee[18];
		pRobot->GetPee(nowPee);
		pRobot->SetPee(nowPee,realBodyPE);

		if(pCMP->count%100==0)
		{
			rt_printf("bodyPm:%f,%f,%f,%f\n%f,%f,%f,%f\n%f,%f,%f,%f\n%f,%f,%f,%f\n",bodyPm[0][0],bodyPm[0][1],bodyPm[0][2],bodyPm[0][3],bodyPm[1][0],bodyPm[1][1],bodyPm[1][2],bodyPm[1][3],bodyPm[2][0],bodyPm[2][1],bodyPm[2][2],bodyPm[2][3],bodyPm[3][0],bodyPm[3][1],bodyPm[3][2],bodyPm[3][3]);
			rt_printf("Fbody:%f,%f,%f,%f,%f,%f\n",Fbody[0],Fbody[1],Fbody[2],Fbody[3],Fbody[4],Fbody[5]);
			rt_printf("realBodyPE:%f,%f,%f,%f,%f,%f\n\n",realBodyPE[0],realBodyPE[1],realBodyPE[2],realBodyPE[3],realBodyPE[4],realBodyPE[5]);
		}

		memcpy(pCMLP.bodyPE_last,bodyPE,sizeof(double)*6);
		memcpy(pCMLP.bodyVel_last,bodyVel,sizeof(double)*6);
		return 1;
	}

    else
	{
		//rt_printf("gait stopping\n");
		for (int i=0;i<6;i++)
		{
			bodyAcc[i]=(Fbody[i]-C[i]*pCMLP.bodyVel_last[i])/M[i];
			bodyVel[i]=pCMLP.bodyVel_last[i]+bodyAcc[i]*deltaT;
			bodyPE[i]=bodyVel[i]*deltaT;
		}

		double pBody[6];
		Aris::DynKer::s_pe2pm(bodyPE,*deltaPm,"213");
		Aris::DynKer::s_pm_dot_pm(*bodyPm,*deltaPm,*realPm);
		Aris::DynKer::s_pm2pe(*realPm,realBodyPE,"313");
		double nowPee[18];
		pRobot->GetPee(nowPee);
		pRobot->SetPee(nowPee,realBodyPE);

		memcpy(pCMLP.bodyPE_last,bodyPE,sizeof(double)*6);
		memcpy(pCMLP.bodyVel_last,bodyVel,sizeof(double)*6);

		if ( fabs(bodyVel[0])<1e-10 && fabs(bodyVel[1])<1e-10 && fabs(bodyVel[2])<1e-10 && fabs(bodyVel[3])<1e-10 && fabs(bodyVel[4])<1e-10 && fabs(bodyVel[5])<1e-10)
			return 0;
		else
			return 1;
	}
}

Aris::Core::MSG parseOpenDoorBegin(const std::string &cmd, const map<std::string, std::string> &params)
{
	CONTINUEMOVE_PARAM param;

	for(auto &i:params)
	{

        if(i.first=="u")
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
    isPull=true;
	isContinue=true;
	isConfirm=false;

	Aris::Core::MSG msg;
	msg.CopyStruct(param);

	std::cout<<"finished parse"<<std::endl;

	return msg;
}

Aris::Core::MSG parseOpenDoorJudge(const std::string &cmd, const map<std::string, std::string> &params)
{
    //   cmd -u -v -w -r -p -y -iscontinue -isstart
	//CONTINUEMOVE_PARAM param;

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
            //param.move_direction[0]=stoi(i.second);

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
            //param.move_direction[1]=stoi(i.second);

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
            //param.move_direction[2]=stoi(i.second);

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
            //param.move_direction[3]=stoi(i.second);

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
            //param.move_direction[4]=stoi(i.second);

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
            //param.move_direction[5]=stoi(i.second);

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
	//msg.CopyStruct(param);

	std::cout<<"finished parse"<<std::endl;

	return msg;
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

void crossMultiply(double * matrix_in1, double *matrix_in2, double * matrix_out)
{
	matrix_out[0]=matrix_in1[1]*matrix_in2[2]-matrix_in1[2]*matrix_in2[1];
	matrix_out[1]=matrix_in1[2]*matrix_in2[0]-matrix_in1[0]*matrix_in2[2];
	matrix_out[2]=matrix_in1[0]*matrix_in2[1]-matrix_in1[1]*matrix_in2[0];
}

double dotMultiply(double *matrix_in1, double *matrix_in2)
{
	double sum{0};

	for(int i=0;i<3;i++)
	{
		sum+=matrix_in1[i]*matrix_in2[i];
	}

	return sum;
}

double norm(double * matrix_in)
{
	return	sqrt(matrix_in[0]*matrix_in[0]+matrix_in[1]*matrix_in[1]+matrix_in[2]*matrix_in[2]);
}

int openDoor(Robots::ROBOT_BASE * pRobot, const Robots::GAIT_PARAM_BASE * pParam)
{

    double bodyVel[6];
    double bodyAcc[6];
    double bodyPE[6];
    double bodyPm[4][4];//bodyPm to Ground every ms

    double beginBodyPm[4][4];//in force control to calculate realBodyPE
    double invBeginBodyPm[4][4];
    double realBodyPm[4][4];

    double location[3][3];
    double invLocation[3][3];
    double planeConst[3]{1,1,1};

    double xBodyInB[3]{1,0,0};
    double yBodyInB[3]{0,1,0};
    double zBodyInB[3]{0,0,1};

    double pushBodyPE313[6];//in position control during push
    double pushPee[18];

    double zeros[6]{0,0,0,0,0,0};

    double Fbody[6]{0,0,0,0,0,0};
    double Fground[6]{0,0,0,0,0,0};
    double C[6]{50,50,50,50,50,50};
    double M[6]{1,1,1,1,1,1};
    double deltaT{0.001};

    double ForceRange[2]{10,90};

    double d0=0.42;//distance from the handle to the middle of the door

	const CONTINUEMOVE_PARAM * pCMP = static_cast<const CONTINUEMOVE_PARAM *>(pParam);

	static CM_LAST_PARAM pCMLP;
	pCMLP.count=pCMP->count;

    Aris::DynKer::s_pe2pm(pCMP->beginBodyPE,*beginBodyPm,"313");
   // rt_printf("beginBodyPe %f %f %f\n",beginBodyPE213[0],beginBodyPE213[1],beginBodyPE213[2]);

    if(isContinue==true)
	{
    	if (pCMP->count<100)
    	{
            if (pCMP->count==0)
            {
            	//initialize
            	for(int i=0;i<3;i++)
            	{
            		pCMLP.forceSum[i]=0;
            		pCMLP.bodyVel_last[i]=0;
            		pCMLP.bodyVel_last[i+3]=0;
            	}
            	Aris::DynKer::s_pe2pe("313",pCMP->beginBodyPE,"213",pCMLP.bodyPE_last);
            	pCMLP.pauseFlag=false;
            	pCMLP.moveState_last=MoveState::None;

            	//start param of now2start in LocateAjust
    			pRobot->GetBodyPe(pCMLP.startPE);
            }

    		pCMLP.moveState=MoveState::None;

    		pCMLP.forceSum[0]+=pCMP->pForceData->at(0).Fx;
    		pCMLP.forceSum[1]+=pCMP->pForceData->at(0).Fy;
    		pCMLP.forceSum[2]+=pCMP->pForceData->at(0).Fz;

    	}
    	else if(pCMP->count==100)
    	{
    		for(int i=0;i<3;i++)
    		{
    			pCMLP.forceAvg[i]=pCMLP.forceSum[i]/100/1000;
    		}
    	}
    	else
    	{
			pCMLP.force[0]=pParam->pForceData->at(0).Fx/1000-pCMLP.forceAvg[0];
			pCMLP.force[1]=pParam->pForceData->at(0).Fy/1000-pCMLP.forceAvg[1];
			pCMLP.force[2]=pParam->pForceData->at(0).Fz/1000-pCMLP.forceAvg[2];
    	}

    	switch(pCMLP.moveState)
    	{
    	case MoveState::None:
            for (int i=0;i<6;i++)
            {
                if (moveDir[2*i+0]==true && moveDir[2*i+1]==false)
                    Fbody[i]=1;
                else if (moveDir[2*i+0]==false && moveDir[2*i+1]==true)
                    Fbody[i]=-1;
                else if (moveDir[2*i+0]==false && moveDir[2*i+1]==false)
                    Fbody[i]=0;
                else
                    rt_printf("parse move diretion wrong!!!\n");
            }

            if (fabs(pCMLP.force[0])>ForceRange[0] && pCMP->count>200 && pCMLP.pauseFlag==false)
            {
                pCMLP.countIter=pCMP->count;
                pRobot->GetBodyPe(pCMLP.pointLocation1);
                memcpy(*location,pCMLP.pointLocation1,sizeof(double)*3);
                pCMLP.moveState=MoveState::PointLocate1;
            }

    		break;

    	case MoveState::PointLocate1:
    		if (pCMP->count-pCMLP.countIter<5000)
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

    		if (fabs(pCMLP.force[0])>ForceRange[0] && (pCMP->count-pCMLP.countIter)>5000)
			{
				pCMLP.countIter=pCMP->count;
				pRobot->GetBodyPe(pCMLP.pointLocation2);
				memcpy(*location+3,pCMLP.pointLocation2,sizeof(double)*3);
				pCMLP.moveState=MoveState::PointLocate2;
			}

    		break;

    	case MoveState::PointLocate2:
    		if (pCMP->count-pCMLP.countIter<5000)
			{
				Fbody[2]=1;
				Fbody[1]=1;
			}
			else
			{
				Fbody[2]=-1;
			}

			if (fabs(pCMLP.force[0])>ForceRange[0] && (pCMP->count-pCMLP.countIter)>5000)
			{
				pCMLP.countIter=pCMP->count+1;
				pRobot->GetBodyPe(pCMLP.pointLocation3);
				memcpy(*location+6,pCMLP.pointLocation3,sizeof(double)*3);
				pCMLP.moveState=MoveState::LocateAjust;

				//calculate the plane of the door. ax+by+cz=1,(a,b,c) is the vertical vector of the plane
				inv3(*location,*invLocation);
				Aris::DynKer::s_dgemm(3,1,3,1,*invLocation,3,planeConst,1,1,pCMLP.planeVertical,1);
				Aris::DynKer::s_inv_pm(*beginBodyPm,*invBeginBodyPm);//have not rotate, beginBodyPm is right here
				Aris::DynKer::s_pm_dot_v3(*invBeginBodyPm,pCMLP.planeVertical,pCMLP.planeVerticalInB);
				if(pCMLP.planeVerticalInB[2]<0)
				{
					for (int i=0;i<3;i++)
					{
						pCMLP.planeVertical[i]=-pCMLP.planeVertical[i];
						pCMLP.planeVerticalInB[i]=-pCMLP.planeVerticalInB[i];
					}
				}
				pCMLP.planeVerticalInB[1]=0;
				pCMLP.planeYPR[0]=acos(dotMultiply(zBodyInB,pCMLP.planeVerticalInB)/norm(pCMLP.planeVerticalInB));
				if(pCMLP.planeVerticalInB[0]<0)
				{
					pCMLP.planeYPR[0]=-pCMLP.planeYPR[0];
				}

				//only use location1 & location2 to calculate horizontal, thus get yaw, stored in YPR[1]&YPR[2]
				for (int i=0;i<3;i++)
				{
					pCMLP.horizontal[i]=pCMLP.pointLocation2[i]-pCMLP.pointLocation1[i];
				}
				Aris::DynKer::s_pm_dot_v3(*invBeginBodyPm,pCMLP.horizontal,pCMLP.horizontalInB);
				pCMLP.planeYPR[1]=-atan(pCMLP.horizontalInB[2]/pCMLP.horizontalInB[0]);
				pCMLP.planeYPR[2]=acos(dotMultiply(xBodyInB,pCMLP.horizontalInB)/norm(pCMLP.horizontalInB));
				if(pCMLP.planeYPR[2]>PI/2)
				{
					pCMLP.planeYPR[2]=-(PI-pCMLP.planeYPR[2]);
				}

				//Set now param of now2start in LocateAjust
				pRobot->GetBodyPe(pCMLP.nowPE);
				pRobot->GetPee(pCMLP.nowPee);

				//Set rotate param
				if(fabs(pCMLP.planeYPR[0])>PI/9)
				{
					pCMLP.walkParam.n=2;
					pCMLP.walkParam.beta=pCMLP.planeYPR[0]*0.88/3*2;
				}
				else
				{
					pCMLP.walkParam.n=1;
					pCMLP.walkParam.beta=pCMLP.planeYPR[0]*0.88*2;
				}
				pCMLP.walkParam.d=0;
				pCMLP.walkParam.alpha=0;
				pCMLP.walkParam.totalCount=2000;

				rt_printf("location1:%f,%f,%f\n",pCMLP.pointLocation1[0],pCMLP.pointLocation1[1],pCMLP.pointLocation1[2]);
				rt_printf("location2:%f,%f,%f\n",pCMLP.pointLocation2[0],pCMLP.pointLocation2[1],pCMLP.pointLocation2[2]);
				rt_printf("location3:%f,%f,%f\n",pCMLP.pointLocation3[0],pCMLP.pointLocation3[1],pCMLP.pointLocation3[2]);
				rt_printf("Vertical:%f,%f,%f\n",pCMLP.planeVertical[0],pCMLP.planeVertical[1],pCMLP.planeVertical[2]);
				rt_printf("VerticalInB:%f,%f,%f\n",pCMLP.planeVerticalInB[0],pCMLP.planeVerticalInB[1],pCMLP.planeVerticalInB[2]);
				rt_printf("horizontalInB:%f,%f,%f\n",pCMLP.horizontalInB[0],pCMLP.horizontalInB[1],pCMLP.horizontalInB[2]);
				rt_printf("yaw:%f,%f,%f\n",pCMLP.planeYPR[0],pCMLP.planeYPR[1],pCMLP.planeYPR[2]);
			}

    		break;

    	case MoveState::LocateAjust:
    		//1.now2start
    		if(pCMP->count-pCMLP.countIter<pCMLP.now2StartCount)
    		{
				for (int i=0;i<6;i++)
				{
					pCMLP.realPE[i]=pCMLP.nowPE[i]+(pCMLP.startPE[i]-pCMLP.nowPE[i])/2*(1-cos((pCMP->count-pCMLP.countIter)*PI/pCMLP.now2StartCount));
				}

				pRobot->SetPee(pCMLP.nowPee,pCMLP.realPE);

				if(pCMP->count-pCMLP.countIter+1==pCMLP.now2StartCount)
				{
					pRobot->GetPee(pCMLP.walkParam.beginPee);
					pRobot->GetBodyPe(pCMLP.walkParam.beginBodyPE);
				}
    		}
    		//2.rotate
    		else
    		{
				pCMLP.walkParam.count=pCMP->count-pCMLP.countIter-pCMLP.now2StartCount;
				pCMLP.ret=Robots::walk(pRobot,& pCMLP.walkParam);

				if(pCMLP.ret==0)
				{
					pCMLP.countIter=pCMP->count;
					pCMLP.moveState=MoveState::Forward;
					//pCMLP.forwardDistance=(0.9/cos(pCMLP.planeYPR[0])-0.9)/0.02*1000;

				}
    		}

    		pRobot->GetBodyPe(bodyPE,"213");
			for(int i=0;i<6;i++)
			{
				bodyVel[i]=bodyPE[i]-pCMLP.bodyPE_last[i];
			}

    		break;

    	case MoveState::Forward:
    		Fbody[2]=-1;
    		/*
		if(pCMP->count-pCMLP.countIter>(int)pCMLP.forwardDistance && fabs(pCMLP.force[0])<ForceRange[0])
    		{
    			pCMLP.countIter=pCMP->count;
    			pRobot->GetPee(pCMLP.endPeeInB,"B");
    			if(isPull==true)
					pCMLP.moveState=MoveState::Rightward;
				else
					pCMLP.moveState=MoveState::Leftward;
    		}*/
    		if(pCMP->count==pCMLP.countIter+1)
    		{
    			pRobot->GetBodyPe(pCMLP.startPE);
    		}

    		if(fabs(pCMLP.force[0])>ForceRange[0])
    		{
    			pCMLP.countIter=pCMP->count;
    			pCMLP.moveState=MoveState::Backward;
    		}

    		break;

    	case MoveState::Backward:
    		Fbody[2]=1;

    		if (pCMP->count-pCMLP.countIter>750)
    		{
    			pCMLP.countIter=pCMP->count;
    			pRobot->GetPee(pCMLP.endPeeInB,"B");
    			if(isPull==true)
    				pCMLP.moveState=MoveState::Rightward;
    			else
    				pCMLP.moveState=MoveState::Leftward;
    		}

    		break;

    	case MoveState::Leftward:

    		Fbody[0]=-1;

    		if(isPull==false)
			{
				if (fabs(pCMLP.force[1])>ForceRange[0])
				{
					pCMLP.countIter=pCMP->count;
					pRobot->GetBodyPe(pCMLP.handlePE);
					pCMLP.moveState=MoveState::Rightward;
				}
				else if (fabs(pCMLP.force[1])<=ForceRange[0] && pCMP->count-pCMLP.countIter>7500)
				{
					pCMLP.countIter=pCMP->count+1;
					pRobot->GetPee(pCMLP.startPeeInB,"B");
					pRobot->GetBodyPe(pCMLP.nowPE);
					pCMLP.moveState=MoveState::Follow;
				}
			}
			else
			{
				if (pCMP->count-pCMLP.countIter>3000)
				{
					pCMLP.countIter=pCMP->count;
					pCMLP.moveState=MoveState::Downward;
					pCMLP.downwardCount=(int)1250*PI;
				    pCMLP.downwardFlag = false;
				}
			}

    		break;

    	case MoveState::Rightward:

    		Fbody[0]=1;

    		if(isPull==true)
    		{
    			if (fabs(pCMLP.force[1])>ForceRange[0])
    			{
    				pCMLP.countIter=pCMP->count;
    				pRobot->GetBodyPe(pCMLP.handlePE);
    				pCMLP.moveState=MoveState::Leftward;
    			}
    			else if (fabs(pCMLP.force[1])<ForceRange[0] && pCMP->count-pCMLP.countIter>7500)
				{
					pCMLP.countIter=pCMP->count+1;
					pRobot->GetPee(pCMLP.startPeeInB,"B");
					pRobot->GetBodyPe(pCMLP.nowPE);
					pCMLP.moveState=MoveState::Follow;
				}
    		}
    		else
    		{
    			if (pCMP->count-pCMLP.countIter>3000)
    			{
    				pCMLP.countIter=pCMP->count;
    				pCMLP.moveState=MoveState::Downward;
    				pCMLP.downwardCount=(int)(PI*1250);
    				pCMLP.downwardFlag = false;
    			}
    		}

    		break;

    	case MoveState::Follow:
    		if(pCMP->count-pCMLP.countIter<pCMLP.followCount)
			{
    			for(int i=0;i<3;i++)
    			{
    				//leg 0,2,4
    				pCMLP.realPeeInB[6*i]=pCMLP.startPeeInB[6*i]+(pCMLP.endPeeInB[6*i]-pCMLP.startPeeInB[6*i])/2*(1-cos((pCMP->count-pCMLP.countIter)*PI/pCMLP.followCount));
					pCMLP.realPeeInB[6*i+1]=pCMLP.startPeeInB[6*i+1]+0.05*(1-cos((pCMP->count-pCMLP.countIter)*2*PI/pCMLP.followCount));
					pCMLP.realPeeInB[6*i+2]=pCMLP.startPeeInB[6*i+2];
					//leg 1,3,5
					pCMLP.realPeeInB[6*i+3]=pCMLP.startPeeInB[6*i+3];
					pCMLP.realPeeInB[6*i+4]=pCMLP.startPeeInB[6*i+4];
					pCMLP.realPeeInB[6*i+5]=pCMLP.startPeeInB[6*i+5];
    			}
			}
    		else if(pCMP->count-pCMLP.countIter<2*pCMLP.followCount)
    		{
    			for(int i=0;i<3;i++)
				{
					//leg 1,3,5
					pCMLP.realPeeInB[6*i+3]=pCMLP.startPeeInB[6*i+3]+(pCMLP.endPeeInB[6*i+3]-pCMLP.startPeeInB[6*i+3])/2*(1-cos((pCMP->count-pCMLP.countIter-pCMLP.followCount)*PI/pCMLP.followCount));
					pCMLP.realPeeInB[6*i+4]=pCMLP.startPeeInB[6*i+4]+0.05*(1-cos((pCMP->count-pCMLP.countIter-pCMLP.followCount)*2*PI/pCMLP.followCount));
					pCMLP.realPeeInB[6*i+5]=pCMLP.startPeeInB[6*i+5];
					//leg 0,2,4
					pCMLP.realPeeInB[6*i]=pCMLP.endPeeInB[6*i];
					pCMLP.realPeeInB[6*i+1]=pCMLP.endPeeInB[6*i+1];
					pCMLP.realPeeInB[6*i+2]=pCMLP.endPeeInB[6*i+2];
				}
    		}
    		pRobot->SetPee(pCMLP.realPeeInB,pCMLP.nowPE,"B");

    		if(pCMP->count-pCMLP.countIter+1==2*pCMLP.followCount)
    		{
    			pCMLP.countIter=pCMP->count;
    			pCMLP.moveState=MoveState::Forward;
    			/*
    			if(isPull==true)
					pCMLP.moveState=MoveState::Rightward;
				else
					pCMLP.moveState=MoveState::Leftward;
					*/
    		}

    		pRobot->GetBodyPe(bodyPE,"213");
			for(int i=0;i<6;i++)
			{
				bodyVel[i]=bodyPE[i]-pCMLP.bodyPE_last[i];
			}

    		break;

    	case MoveState::Downward:

    		Fbody[1]=-1;
    		if (pCMLP.downwardFlag)
            {
    			if(isPull==true)
    			{
    				Fbody[1]=-cos((pCMP->count-pCMLP.countIter)*PI/pCMLP.downwardCount/2);
    				Fbody[0]=sin((pCMP->count-pCMLP.countIter)*PI/pCMLP.downwardCount/2);
    			}
    			else
    			{
    				Fbody[1]=-cos((pCMP->count-pCMLP.countIter)*PI/pCMLP.downwardCount/2);
    				Fbody[0]=-sin((pCMP->count-pCMLP.countIter)*PI/pCMLP.downwardCount/2);
    			}
            }

            if(pCMLP.downwardFlag==false && fabs(pCMP->pForceData->at(0).Fz/1000-pCMLP.forceAvg[2])>ForceRange[0])
            {
            	pCMLP.downwardFlag = true;
            	pCMLP.countIter=pCMP->count;
            }

            //pCMLP.forceYZ=sqrt((pCMLP.force[1])*(pCMLP.force[1])+(pCMLP.force[2])*(pCMLP.force[2]));

    		if (fabs(pCMLP.force[1])>ForceRange[1] && pCMP->count-pCMLP.countIter>pCMLP.downwardCount/2)
    		{
    			pCMLP.countIter=pCMP->count;
    			if(isPull==true)
    				pCMLP.moveState=MoveState::Pullhandle;
    			else
    				pCMLP.moveState=MoveState::Pushhandle;
    		}

    		break;

    	case MoveState::Pullhandle:

    		Fbody[2]=1;

    		if (pCMP->count-pCMLP.countIter>2500)
    			pCMLP.moveState=MoveState::Upward;

    		break;

    	case MoveState::Pushhandle:

    		Fbody[2]=-1;

    		if (pCMP->count-pCMLP.countIter>2500)
            {
                pCMLP.moveState=MoveState::Upward;
            }
    		break;

    	case MoveState::Upward:

			Fbody[1]=1;
			if(isPull==true)
			{
				if (pCMP->count-pCMLP.countIter>5000)
				{
					isContinue=false;
				}
			}
			else
			{
				if (pCMP->count-pCMLP.countIter>10500)
				{
					pCMLP.moveState=MoveState::PrePush;
				}
			}
			break;

        case MoveState::PrePush:

            if ( fabs(bodyVel[0])<1e-10 && fabs(bodyVel[1])<1e-10 && fabs(bodyVel[2])<1e-10 && fabs(bodyVel[3])<1e-10 && fabs(bodyVel[4])<1e-10 && fabs(bodyVel[5])<1e-10)
            {
                pCMLP.countIter=pCMP->count+1;
                pCMLP.moveState=MoveState::Push;
                pCMLP.pushState=PushState::now2Start;

                //now2start param
                pRobot->GetBodyPe(pCMLP.nowPE);
				pRobot->GetBodyPm(* pCMLP.nowPm);
				pRobot->GetPee(pCMLP.nowPee);
				for(int i=0;i<3;i++)
				{
					pCMLP.now2startDistance[i]=pCMLP.startPE[i]-pCMLP.nowPE[i];
					pCMLP.handle2startDistance[i]=pCMLP.startPE[i]-pCMLP.handlePE[i];
				}
				Aris::DynKer::s_pm_dot_v3(*pCMLP.nowPm,xBodyInB,pCMLP.xNowInG);
				for(int i=0;i<3;i++)
				{
					pCMLP.now2startDistanceModified[i]=dotMultiply(pCMLP.now2startDistance,pCMLP.xNowInG)*pCMLP.xNowInG[i];
					pCMLP.now2startDistanceModified[i+3]=pCMLP.startPE[i+3]-pCMLP.nowPE[i+3];
				}
            }

            break;

    	case MoveState::Push:
            //Three Step using position control
            switch(pCMLP.pushState)
            {
            case PushState::now2Start:
                //1.Move back to startPE;
                for (int i=0;i<6;i++)
                {
                    pCMLP.realPE[i]=pCMLP.nowPE[i]+(pCMLP.now2startDistanceModified[i])/2*(1-cos((pCMP->count-pCMLP.countIter)*PI/pCMLP.now2StartCount));
                }

                pRobot->SetPee(pCMLP.nowPee,pCMLP.realPE);

                if(pCMP->count-pCMLP.countIter+1==pCMLP.now2StartCount)
                {
                    if(isContinue==true)
                    {
                        pCMLP.pushState=PushState::rightWalk;
                        pCMLP.countIter=pCMP->count+1;

                        pCMLP.walkParam.n=2;
                        pCMLP.walkParam.alpha=-PI/2;
                        pCMLP.walkParam.beta=0;
                        pCMLP.walkParam.totalCount=2000;
                        pCMLP.walkParam.d=(d0-fabs(dotMultiply(pCMLP.handle2startDistance,pCMLP.xNowInG)))/3*2;
                        pRobot->GetPee(pCMLP.walkParam.beginPee);
                        pRobot->GetBodyPe(pCMLP.walkParam.beginBodyPE);
                    }
                    else
                    {
                        pRobot->SetPee(pCMLP.nowPee,pCMLP.startPE);
                    }
                }

                break;

            case PushState::rightWalk:
                //2.Move rightward to align with the door;
                pCMLP.walkParam.count=pCMP->count-pCMLP.countIter;

                pCMLP.ret=Robots::walk(pRobot,& pCMLP.walkParam);

                if(pCMLP.ret==0)
                {
                    if(isContinue==true)
                    {
                        pCMLP.pushState=PushState::forwardWalk;
                        pCMLP.countIter=pCMP->count+1;
                        pCMLP.forwardWalkFlag=0;

                        pCMLP.walkParam.totalCount=2000;
                        pCMLP.walkParam.n=4;
                        pCMLP.walkParam.alpha=0;
                        pCMLP.walkParam.beta=0;
                        pCMLP.walkParam.d=0.5;
                        pRobot->GetPee(pCMLP.walkParam.beginPee);
                        pRobot->GetBodyPe(pCMLP.walkParam.beginBodyPE);
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
                pCMLP.walkParam.count=pCMP->count-pCMLP.countIter;

                pCMLP.ret=Robots::walk(pRobot,& pCMLP.walkParam);

                if(pCMLP.ret==0)
                {
                    if(isContinue==true)
                    {
                     //   if(pCMLP.forwardWalkFlag<5)
                       // {
                         //   pCMLP.forwardWalkFlag++;
                           // pCMLP.countIter=pCMP->count+1;
                           // pRobot->GetPee(pCMLP.walkParam.beginPee);
                           // pRobot->GetBodyPe(pCMLP.walkParam.beginBodyPE);
                       // }
                       // else
                       // {
                            return 0;
                       // }
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

            pRobot->GetBodyPe(bodyPE,"213");
            for(int i=0;i<6;i++)
            {
                bodyVel[i]=bodyPE[i]-pCMLP.bodyPE_last[i];
            }

    		break;

    	default:
    		break;
    	}

        if(pCMLP.moveState!=MoveState::Push && pCMLP.moveState!=MoveState::LocateAjust && pCMLP.moveState!=MoveState::Follow)
        {
            if (isConfirm==false && pCMLP.pauseFlag==false)
            {
                pCMLP.moveState_last=pCMLP.moveState;
                pCMLP.moveState=MoveState::None;
                pCMLP.pauseFlag=true;
            }
            else if(isConfirm==false && pCMLP.pauseFlag==true)
            {
                pCMLP.moveState=MoveState::None;
                pCMLP.pauseCount++;
                //Do not change F to zeros here, so that we can move the robot manually when paused
            }
            else if(isConfirm==true && pCMLP.pauseFlag==true)
            {
                pCMLP.moveState=pCMLP.moveState_last;
                pCMLP.pauseFlag=false;
                pCMLP.pauseCount++;
            }
            else//isConfirm=true && pauseFlag=false
            {
                pCMLP.countIter+=pCMLP.pauseCount;
                pCMLP.pauseCount=0;
            }
/*
            for (int i=0;i<3;i++)
            {
                if (fabs(pCMLP.bodyPE_last[i])>pCMLP.posLimit[i])
                {
                    Fbody[0]=0;
                    if(pCMP->count%300==0)
                    	rt_printf("position of direction %d is out of limit\n",i);
                }
            }
*/
            pRobot->GetBodyPm(*bodyPm);
            Aris::DynKer::s_pm_dot_v3(*bodyPm,Fbody,Fground);

            for (int i=0;i<6;i++)
            {
                bodyAcc[i]=(Fground[i]-C[i]*pCMLP.bodyVel_last[i])/M[i];
                bodyVel[i]=pCMLP.bodyVel_last[i]+bodyAcc[i]*deltaT;
                bodyPE[i]=pCMLP.bodyPE_last[i]+bodyVel[i]*deltaT;
            }

            double pBody[6];
            Aris::DynKer::s_pe2pe("213",bodyPE,"313",pBody);
            double nowPee[18];
            pRobot->GetPee(nowPee);
            pRobot->SetPee(nowPee,pBody);
        }

        memcpy(pCMLP.bodyPE_last,bodyPE,sizeof(double)*6);
        memcpy(pCMLP.bodyVel_last,bodyVel,sizeof(double)*6);

        if(pCMP->count%1000==0)
        {
            rt_printf("count:%d,countIter:%d,state:%d,Fz:%f,Fy:%f,Fx:%f,Fyz:%f\n",pCMP->count,pCMLP.countIter,pCMLP.moveState,fabs(pCMP->pForceData->at(0).Fz/1000-pCMLP.forceAvg[2]),fabs(pCMP->pForceData->at(0).Fy/1000-pCMLP.forceAvg[1]),fabs(pCMP->pForceData->at(0).Fx/1000-pCMLP.forceAvg[0]),fabs(pCMLP.forceYZ));
        }

		openDoorPipe.SendToNRT(pCMLP);

		return 1;
	}

    else
	{
		//rt_printf("gait stopping\n");

		memcpy(Fground,zeros,sizeof(double)*6);
		for (int i=0;i<6;i++)
		{
			bodyAcc[i]=(Fground[i]-C[i]*pCMLP.bodyVel_last[i])/M[i];
			bodyVel[i]=pCMLP.bodyVel_last[i]+bodyAcc[i]*deltaT;
			bodyPE[i]=pCMLP.bodyPE_last[i]+bodyVel[i]*deltaT;
		}

		double nowPee[18];
		pRobot->GetPee(nowPee);
		pRobot->SetPee(nowPee,bodyPE,"G","213");

		memcpy(pCMLP.bodyPE_last,bodyPE,sizeof(double)*6);
		memcpy(pCMLP.bodyVel_last,bodyVel,sizeof(double)*6);

		//rt_printf("%d,%f,%f,%f,%f,%f,%f\n",pCMP->count,beginBodyPE213[0],beginBodyPE213[1],beginBodyPE213[2],beginBodyPE213[3],beginBodyPE213[4],beginBodyPE213[5]);
		//rt_printf("%d,%f,%f,%f,%f,%f,%f\n",pCMP->count,bodyPE[0],bodyPE[1],bodyPE[2],bodyPE[3],bodyPE[4],bodyPE[5]);

		if ( fabs(bodyVel[0])<1e-10 && fabs(bodyVel[1])<1e-10 && fabs(bodyVel[2])<1e-10 && fabs(bodyVel[3])<1e-10 && fabs(bodyVel[4])<1e-10 && fabs(bodyVel[5])<1e-10)
			return 0;
		else
			return 1;
	}
}
