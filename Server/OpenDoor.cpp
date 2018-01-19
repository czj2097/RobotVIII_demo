#include "OpenDoor.h"

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

Pipe<OpenDoorParam> openDoorPipe(true);
std::atomic_bool isForce;
std::atomic_bool isContinue;
std::atomic_int moveDir[6];
std::atomic_bool isPull;
std::atomic_bool isLeft;
std::atomic_bool isConfirm;
std::atomic_bool isJump;

void parseContinueMoveBegin(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
{
    ContinueMoveParam param;

    for(auto &i:params)
    {
        if(i.first=="u")
        {
            moveDir[0]=std::stoi(i.second);
        }
        else if(i.first=="v")
        {
            moveDir[1]=std::stoi(i.second);
        }
        else if(i.first=="w")
        {
            moveDir[2]=std::stoi(i.second);
        }
        else if(i.first=="yaw")
        {
            moveDir[3]=std::stoi(i.second);
        }
        else if(i.first=="pitch")
        {
            moveDir[4]=std::stoi(i.second);
        }
        else if(i.first=="roll")
        {
            moveDir[5]=std::stoi(i.second);
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
            moveDir[0]=std::stoi(i.second);
        }
        else if(i.first=="v")
        {
            moveDir[1]=std::stoi(i.second);
        }
        else if(i.first=="w")
        {
            moveDir[2]=std::stoi(i.second);
        }
        else if(i.first=="yaw")
        {
            moveDir[3]=std::stoi(i.second);
        }
        else if(i.first=="pitch")
        {
            moveDir[4]=std::stoi(i.second);
        }
        else if(i.first=="roll")
        {
            moveDir[5]=std::stoi(i.second);
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
            ForceTask::forceInit(param.count,param.force_data->at(0).fce,forceInF);
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

            //Judge whether to rotate in priority. If not, then judge whether to translate
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



void OpenDoor::parseOpenDoor(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
{
    aris::server::GaitParamBase param;

    for(auto &i:params)
    {
        if(i.first=="init")
        {
            ODState=OpenDoorState::Init;
            msg.copyStruct(param);
        }
        else if(i.first=="doorLocate")
        {
            ODState=OpenDoorState::DoorLocate;
            DLState=DoorLocateState::Point1;
        }
        else if(i.first=="handleLocate")
        {
            ODState=OpenDoorState::HandleLocate;
        }
        else if(i.first=="turnHandle")
        {
            ODState=OpenDoorState::TurnHandle;
        }
        else if(i.first=="pushDoor")
        {
            ODState=OpenDoorState::PushDoor;
        }
        else if(i.first=="quit")
        {
            ODState=OpenDoorState::Quit;
        }
        else if(i.first=="isForward")
        {
            if(i.second=="1")
                isPull=false;
            else
                isPull=true;
        }
        else if(i.first=="isLeft")
        {
            if(i.second=="1")
                isLeft=true;
            else
                isLeft=false;
        }
        else if(i.first=="isJump")
        {
            if(i.second=="1")
                isJump=true;
            else
                isJump=false;
        }
        else if(i.first=="u")
        {
            moveDir[0]=std::stoi(i.second);
        }
        else if(i.first=="v")
        {
            moveDir[1]=std::stoi(i.second);
        }
        else if(i.first=="w")
        {
            moveDir[2]=std::stoi(i.second);
        }
        else if(i.first=="yaw")
        {
            moveDir[3]=std::stoi(i.second);
        }
        else if(i.first=="pitch")
        {
            moveDir[4]=std::stoi(i.second);
        }
        else if(i.first=="roll")
        {
            moveDir[5]=std::stoi(i.second);
        }
        else
        {
            std::cout<<"parse failed"<<std::endl;
        }
    }
    isContinue=true;
    isConfirm=false;

    std::cout<<"finished parse"<<std::endl;
}

void parseOpenDoorJudge(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
{
    for(auto &i:params)
    {
        if(i.first=="isQuit")
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
            moveDir[0]=std::stoi(i.second);
        }
        else if(i.first=="v")
        {
            moveDir[1]=std::stoi(i.second);
        }
        else if(i.first=="w")
        {
            moveDir[2]=std::stoi(i.second);
        }
        else if(i.first=="yaw")
        {
            moveDir[3]=std::stoi(i.second);
        }
        else if(i.first=="pitch")
        {
            moveDir[4]=std::stoi(i.second);
        }
        else if(i.first=="roll")
        {
            moveDir[5]=std::stoi(i.second);
        }
        else
        {
            std::cout<<"parse failed"<<std::endl;
        }
    }

    std::cout<<"finished parse"<<std::endl;
}

//*****only for Robot VIII*****
int OpenDoor::openDoor(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
{
    auto &robot = static_cast<Robots::RobotBase &>(model);
    auto &param = static_cast<const aris::server::GaitParamBase &>(param_in);

    //for data record
    double pEE[18];
    double bodyVel123[6]{0,0,0,0,0,0};

    //PushState & PullState
    double currentPe[6] {0};
    double currentPee[18] {0};
    double xBodyInB[3]{1,0,0};
    double yBodyInB[3]{0,1,0};
    double pushBodyPE313[6];//for pause
    double pushPee[18];//for pause
    double handle2MiddleDist=0.55;//distance from the handle to the middle of the door

    double handleLen=0.1;//for pull

    //Position Generetion From Force
    double bodyVel[6]{0,0,0,0,0,0};
    double bodyAcc[6];
    double bodyPE[6];//delta PE every ms
    double bodyPm[4][4];//bodyPm to Ground every ms
    double deltaPE[6];
    double deltaPm[4][4];
    double realPm[4][4];
    double realPE[6];

    ODP.count=param.count;
    if(param.count>200)
    {
        double forceInF[6];
        ForceTask::forceInit(param.count,param.ruicong_data->at(0).force[0].fce,forceInF);
        aris::dynamic::s_f2f(*robot.forceSensorMak().prtPm(),forceInF,ODP.forceInB);
        ODP.filter.DoFilter(ODP.forceInB,ODP.forceInB_filtered);
    }

    int ret;
    switch(ODState)
    {
    case OpenDoorState::Init:
        initialize(robot);
        ODState=OpenDoorState::WaitCmd;
        break;
    case OpenDoorState::WaitCmd:
        //Do nothing
        break;
    case OpenDoorState::DoorLocate:
        ret=locateDoor(robot,param.count);
        if(ret==0)
        {
            GetAdjustParam(robot);
            ODState=OpenDoorState::AdjustD2H;
        }
        break;
    case OpenDoorState::AdjustD2H:
        ret=adjustRobotD2H(robot,param.count);
        if(ret==0)
        {
            ODState=OpenDoorState::HandleLocate;
            HLState=HandleLocateState::Forward;
        }
        break;
    case OpenDoorState::HandleLocate:
        ret=locateHandle(robot,param.count);
        if(ret==0)
        {
            ODState=OpenDoorState::TurnHandle;
        }
        break;
    case OpenDoorState::TurnHandle:
        ret=turnHandle(robot,param.count);
        if(ret==0)
        {
            if(isPull==true)
                ODState=OpenDoorState::PullDoor;
            else
                ODState=OpenDoorState::PushDoor;
        }
        break;
    case OpenDoorState::PushDoor:
        break;
    case OpenDoorState::PullDoor:
        break;
    case OpenDoorState::AdjustP2W:
        break;
    case OpenDoorState::WalkThrough:
        break;
    case OpenDoorState::Quit:
        break;
    default:
        break;
    }

    if(isContinue==true)
    {
        switch(ODP.moveState)
        {
        case MoveState::None:
            for (int i=0;i<6;i++)
            {
                Fbody[i]=moveDir[i];
            }

            if (fabs(ODP.forceInB_filtered[2])>ForceRange[0] && param.count>200 && ODP.isPause==false)
            {
                ODP.countIter=param.count+1;
                robot.GetPeb(ODP.pointLocation1);
                memcpy(*ODP.location,ODP.pointLocation1,sizeof(double)*3);
                ODP.moveState=MoveState::PointLocate1;
            }

            break;

        case MoveState::PointLocate1:
            if (param.count-ODP.countIter<2500)
            {
                Fbody[2]=1;
                if(isLeft==true)
                    Fbody[0]=1;
                else
                    Fbody[0]=-1;
            }
            else
            {
                Fbody[2]=-1;
            }

            if (fabs(ODP.forceInB_filtered[2])>ForceRange[0] && (param.count-ODP.countIter)>=2500)
            {
                ODP.countIter=param.count+1;
                robot.GetPeb(ODP.pointLocation2);
                memcpy(*ODP.location+3,ODP.pointLocation2,sizeof(double)*3);
                ODP.moveState=MoveState::PointLocate2;
            }

            break;

        case MoveState::PointLocate2:
            if (param.count-ODP.countIter<1500)
            {
                Fbody[2]=1;
                Fbody[1]=1;
            }
            else
            {
                Fbody[2]=-1;
            }

            if (fabs(ODP.forceInB_filtered[2])>ForceRange[0] && (param.count-ODP.countIter)>=1500)
            {
                ODP.countIter=param.count+1;
                robot.GetPeb(ODP.pointLocation3);
                memcpy(*ODP.location+6,ODP.pointLocation3,sizeof(double)*3);
                ODP.moveState=MoveState::LocateAjust;

                //calculate the plane of the door. ax+by+cz=1,(a,b,c) is the vertical vector of the plane
                NormalGait::inv3(*ODP.location,*invLocation);
                aris::dynamic::s_dgemm(3,1,3,1,*invLocation,3,planeConst,1,1,planeVertical,1);
                aris::dynamic::s_inv_pm(*beginBodyPm,*invBeginBodyPm);//have not rotate, beginBodyPm is right here
                aris::dynamic::s_pm_dot_v3(*invBeginBodyPm,planeVertical,planeVerticalInB);

                ODP.planeYPR[0]=atan(planeVerticalInB[0]/planeVerticalInB[2]);
                ODP.planeYPR[1]=-asin(planeVerticalInB[1]/NormalGait::norm(planeVerticalInB));
                ODP.planeYPR[2]=0;

                //Set now param of now2start in LocateAjust
                robot.GetPeb(ODP.nowPE);
                robot.GetPee(ODP.nowPee);

                //Set rotate param
                if(fabs(ODP.planeYPR[0])>=PI/4)
                {
                    rt_printf("WARNING!!! Too large angle between the door and the robot!");
                    return 1;
                }
                else if(fabs(ODP.planeYPR[0])>PI/9)
                {
                    ODP.walkParam.n=2;
                    ODP.walkParam.beta=ODP.planeYPR[0]*rotateRatio/3*2;
                }
                else
                {
                    ODP.walkParam.n=1;
                    ODP.walkParam.beta=ODP.planeYPR[0]*rotateRatio*2;
                }
                ODP.walkParam.d=0;
                ODP.walkParam.alpha=0;
                ODP.walkParam.totalCount=1000;

                rt_printf("yaw:%f,pitch:%f,roll:%f\n",ODP.planeYPR[0],ODP.planeYPR[1],ODP.planeYPR[2]);
            }

            break;

        case MoveState::LocateAjust:
            //now2start & rotate
            ODP.walkParam.count=param.count-ODP.countIter;
            ODP.ret=Robots::walkGait(robot,ODP.walkParam);
            robot.GetPeb(currentPe);
            robot.GetPee(currentPee);

            for (int i=0;i<3;i++)
            {
                currentPe[i]=ODP.nowPE[i]+(ODP.startPE[i]-ODP.nowPE[i])/2*(1-cos((param.count-ODP.countIter)*PI/(2*ODP.walkParam.n*ODP.walkParam.totalCount)));
            }
            robot.SetPeb(currentPe);
            robot.SetPee(currentPee);

            if(ODP.ret==0)
            {
                ODP.countIter=param.count+1;
                ODP.moveState=MoveState::Forward;
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
            if(param.count==ODP.countIter)
            {
                robot.GetPeb(ODP.beginPE);//For calculation of DoorLocation
            }

            if(fabs(ODP.forceInB_filtered[2])>ForceRange[0])
            {
                ODP.countIter=param.count+1;
                ODP.moveState=MoveState::Backward;

                robot.GetPeb(ODP.nowPE);//For calculation of DoorLocation
                for (int i=0;i<3;i++)
                {
                    ODP.vector0[i]=ODP.nowPE[i]-ODP.beginPE[i];//calculation of DoorLocation
                }
            }

            break;

        case MoveState::Backward:
            Fbody[2]=1;

            if (param.count-ODP.countIter>=1000)
            {
                ODP.countIter=param.count+1;
                robot.GetPee(ODP.endPeeInB,robot.body());
                if(isLeft==false)
                    ODP.moveState=MoveState::Rightward;
                else
                    ODP.moveState=MoveState::Leftward;
            }

            break;

        case MoveState::Leftward:

            Fbody[0]=-1;

            if(param.count==ODP.countIter)
            {
                robot.GetPeb(ODP.beginPE);//For calculation of DoorLocation
            }

            if (fabs(ODP.forceInB_filtered[0])>ForceRange[0] || param.count-ODP.countIter>3000)
            {
                ODP.countIter=param.count+1;
                robot.GetPee(ODP.startPeeInB,robot.body());
                robot.GetPeb(ODP.nowPE);
                ODP.moveState=MoveState::Follow;
                if(fabs(ODP.forceInB_filtered[0])>ForceRange[0])
                {
                    ODP.isLastFollow=true;

                    for (int i=0;i<3;i++)
                    {
                        ODP.vector1[i]=ODP.nowPE[i]-ODP.beginPE[i];//calculation of DoorLocation
                    }
                }
            }

            break;

        case MoveState::Rightward:

            Fbody[0]=1;

            if(param.count==ODP.countIter)
            {
                robot.GetPeb(ODP.beginPE);//For calculation of DoorLocation
            }

            if (fabs(ODP.forceInB_filtered[0])>ForceRange[0] || param.count-ODP.countIter>3000)
            {
                ODP.countIter=param.count+1;
                robot.GetPee(ODP.startPeeInB,robot.body());
                robot.GetPeb(ODP.nowPE);
                ODP.moveState=MoveState::Follow;
                if(fabs(ODP.forceInB_filtered[0])>ForceRange[0])
                {
                    ODP.isLastFollow=true;

                    for (int i=0;i<3;i++)
                    {
                        ODP.vector1[i]=ODP.nowPE[i]-ODP.beginPE[i];//calculation of DoorLocation
                    }
                }
            }

            break;

        case MoveState::Follow:
            if(param.count-ODP.countIter<ODP.followCount)
            {
                for(int i=0;i<3;i++)
                {
                    //leg 0,2,4
                    currentPee[6*i]=ODP.startPeeInB[6*i]+(ODP.endPeeInB[6*i]-ODP.startPeeInB[6*i])/2*(1-cos((param.count-ODP.countIter)*PI/ODP.followCount));
                    currentPee[6*i+1]=ODP.startPeeInB[6*i+1]+0.025*(1-cos((param.count-ODP.countIter)*2*PI/ODP.followCount));
                    currentPee[6*i+2]=ODP.startPeeInB[6*i+2];
                    //leg 1,3,5
                    currentPee[6*i+3]=ODP.startPeeInB[6*i+3];
                    currentPee[6*i+4]=ODP.startPeeInB[6*i+4];
                    currentPee[6*i+5]=ODP.startPeeInB[6*i+5];
                }
            }
            else if(param.count-ODP.countIter<2*ODP.followCount)
            {
                for(int i=0;i<3;i++)
                {
                    //leg 1,3,5
                    currentPee[6*i+3]=ODP.startPeeInB[6*i+3]+(ODP.endPeeInB[6*i+3]-ODP.startPeeInB[6*i+3])/2*(1-cos((param.count-ODP.countIter-ODP.followCount)*PI/ODP.followCount));
                    currentPee[6*i+4]=ODP.startPeeInB[6*i+4]+0.025*(1-cos((param.count-ODP.countIter-ODP.followCount)*2*PI/ODP.followCount));
                    currentPee[6*i+5]=ODP.startPeeInB[6*i+5];
                    //leg 0,2,4
                    currentPee[6*i]=ODP.endPeeInB[6*i];
                    currentPee[6*i+1]=ODP.endPeeInB[6*i+1];
                    currentPee[6*i+2]=ODP.endPeeInB[6*i+2];
                }
            }
            robot.SetPeb(ODP.nowPE);
            robot.SetPee(currentPee,robot.body());

            if(param.count-ODP.countIter+1==2*ODP.followCount)
            {
                ODP.countIter=param.count+1;
                if(ODP.isLastFollow==false)
                {
                    ODP.moveState=MoveState::Forward;
                }
                else
                {
                    robot.GetPeb(ODP.jumpStartPe);
                    robot.GetPmb(*ODP.nowPm);
                    robot.GetPee(ODP.nowPee);
                    if(isJump==false)
                    {
                        ODP.moveState=MoveState::NonJump;
                        robot.GetPeb(ODP.startPE);//For now2start in push & pull
                    }
                    else
                    {
                        ODP.moveState=MoveState::Jump;
                    }
                }
            }

            robot.GetPmb(*bodyPm);
            robot.GetPeb(bodyPE,"213");
            for(int i=0;i<6;i++)
            {
                bodyVel[i]=bodyPE[i]-ODP.bodyPE_last[i];
            }

            break;

        case MoveState::NonJump:
            nonJumpDistInB[0]=(isLeft==true ? nonJumpDist : -nonJumpDist);
            aris::dynamic::s_pm_dot_v3(*ODP.nowPm,nonJumpDistInB,nonJumpDistInG);

            memcpy(currentPe,ODP.jumpStartPe,6*sizeof(double));
            currentPe[0]=ODP.jumpStartPe[0]+nonJumpDistInG[0]/2*(1-cos((param.count-ODP.countIter)*PI/ODP.jumpCount));
            currentPe[2]=ODP.jumpStartPe[2]+nonJumpDistInG[2]/2*(1-cos((param.count-ODP.countIter)*PI/ODP.jumpCount));

            robot.SetPeb(currentPe);
            robot.SetPee(ODP.nowPee);

            if (param.count-ODP.countIter==ODP.jumpCount-1)
            {
                ODP.moveState=MoveState::Downward;
                ODP.downwardCount=(int)(PI/2*2500);
                ODP.downwardFlag = false;
            }

            robot.GetPmb(*bodyPm);
            robot.GetPeb(bodyPE,"213");
            for(int i=0;i<6;i++)
            {
                bodyVel[i]=bodyPE[i]-ODP.bodyPE_last[i];
            }
            break;

        case MoveState::Jump:
            jumpDistInB[0]=(isLeft==true ? -jumpDist : jumpDist);
            aris::dynamic::s_pm_dot_v3(*ODP.nowPm,jumpDistInB,jumpDistInG);

            if(param.count-ODP.countIter<ODP.jumpCount)
            {
                memcpy(currentPe,ODP.jumpStartPe,6*sizeof(double));
                memcpy(currentPee,ODP.nowPee,18*sizeof(double));
                currentPe[1]=ODP.jumpStartPe[1]+jumpH/2*(1-cos((param.count-ODP.countIter)*PI/ODP.jumpCount));
                for(int i=0;i<3;i++)
                {
                    currentPee[6*i]=ODP.nowPee[6*i]+jumpDistInG[0]/2*(1-cos((param.count-ODP.countIter)*PI/ODP.jumpCount));
                    currentPee[6*i+1]=ODP.nowPee[6*i+1]+0.025*(1-cos((param.count-ODP.countIter)*2*PI/ODP.followCount));
                    currentPee[6*i+2]=ODP.nowPee[6*i+2]+jumpDistInG[2]/2*(1-cos((param.count-ODP.countIter)*PI/ODP.jumpCount));
                }
                if(param.count-ODP.countIter==ODP.jumpCount-1)
                {
                    ODP.jumpStartPe[1]=currentPe[1];
                }
            }
            else if(param.count-ODP.countIter<2*ODP.jumpCount)
            {
                memcpy(currentPe,ODP.jumpStartPe,6*sizeof(double));
                currentPe[0]=ODP.jumpStartPe[0]+jumpDistInG[0]/2*(1-cos((param.count-ODP.jumpCount-ODP.countIter)*PI/ODP.jumpCount));
                currentPe[2]=ODP.jumpStartPe[2]+jumpDistInG[2]/2*(1-cos((param.count-ODP.jumpCount-ODP.countIter)*PI/ODP.jumpCount));
                for(int i=0;i<3;i++)
                {
                    currentPee[6*i]=ODP.nowPee[6*i]+jumpDistInG[0];
                    currentPee[6*i+1]=ODP.nowPee[6*i+1];
                    currentPee[6*i+2]=ODP.nowPee[6*i+2]+jumpDistInG[2];

                    currentPee[6*i+3]=ODP.nowPee[6*i+3]+jumpDistInG[0]/2*(1-cos((param.count-ODP.jumpCount-ODP.countIter)*PI/ODP.jumpCount));
                    currentPee[6*i+4]=ODP.nowPee[6*i+4]+0.025*(1-cos((param.count-ODP.jumpCount-ODP.countIter)*2*PI/ODP.jumpCount));
                    currentPee[6*i+5]=ODP.nowPee[6*i+5]+jumpDistInG[2]/2*(1-cos((param.count-ODP.jumpCount-ODP.countIter)*PI/ODP.jumpCount));
                }
                if(param.count-ODP.countIter==2*ODP.jumpCount-1)
                {
                    ODP.moveState=MoveState::Downward;
                    ODP.downwardCount=(int)(PI/2*2500);
                    ODP.downwardFlag = false;
                    robot.GetPeb(ODP.startPE);//For now2start in push & pull
                }
            }
            robot.SetPeb(currentPe);
            robot.SetPee(currentPee);

            robot.GetPmb(*bodyPm);
            robot.GetPeb(bodyPE,"213");
            for(int i=0;i<6;i++)
            {
                bodyVel[i]=bodyPE[i]-ODP.bodyPE_last[i];
            }

            break;

        case MoveState::Downward:
            Fbody[1]=-1;

            if(param.count==ODP.countIter)
            {
                robot.GetPeb(ODP.beginPE);//For calculation of DoorLocation
            }

            if (ODP.downwardFlag)
            {
                if(isLeft==false)
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

            if(ODP.downwardFlag==false && fabs(ODP.forceInB_filtered[1])>ForceRange[0])
            {
                ODP.downwardFlag = true;
                ODP.countIter=param.count+1;

                robot.GetPeb(ODP.nowPE);//For calculation of DoorLocation
                for (int i=0;i<3;i++)
                {
                    ODP.vector2[i]=ODP.nowPE[i]-ODP.beginPE[i];
                }
                ODP.handleLocation[2]=-NormalGait::norm(ODP.vector0);
                ODP.handleLocation[0]=NormalGait::norm(ODP.vector1);
                ODP.handleLocation[1]=-NormalGait::norm(ODP.vector2);
                rt_printf("handleLocation:%f,%f,%f\n",ODP.handleLocation[0],ODP.handleLocation[1],ODP.handleLocation[2]);
            }

            if (fabs(ODP.forceInB_filtered[0])>ForceRange[1] || fabs(ODP.forceInB_filtered[1])>ForceRange[1])
            {
                ODP.countIter=param.count+1;
                if(isPull==true)
                    ODP.moveState=MoveState::Pullhandle;
                else
                    ODP.moveState=MoveState::Pushhandle;
            }

            break;

        case MoveState::Pullhandle:

            Fbody[2]=1;

            if (param.count-ODP.countIter>2000)
                ODP.moveState=MoveState::PrePull;

            break;

        case MoveState::Pushhandle:

            Fbody[2]=-1;

            if (param.count-ODP.countIter>2000)
            {
                ODP.moveState=MoveState::PrePush;
            }
            break;

        case MoveState::PrePull:
            if ( fabs(bodyVel[0])<1e-10 && fabs(bodyVel[1])<1e-10 && fabs(bodyVel[2])<1e-10
              && fabs(bodyVel[3])<1e-10 && fabs(bodyVel[4])<1e-10 && fabs(bodyVel[5])<1e-10)
            {
                ODP.countIter=param.count+1;
                ODP.moveState=MoveState::Pull;
                ODP.pullState=PullState::pull2Start;

                //now2start param
                robot.GetPeb(ODP.nowPE);
                robot.GetPmb(*ODP.nowPm);
                robot.GetPee(ODP.nowPee);
                for(int i=0;i<3;i++)
                {
                    ODP.now2startDistance[i]=ODP.startPE[i]-ODP.nowPE[i]; //startPE recorded at the last forward
                }

                aris::dynamic::s_pm_dot_v3(*ODP.nowPm,xBodyInB,ODP.xNowInG);
                aris::dynamic::s_pm_dot_v3(*ODP.nowPm,yBodyInB,ODP.yNowInG);

                if((isJump==false && isLeft==false) || (isJump==true && isLeft==true))
                {
                    ODP.now2startDistanceInB[0]=aris::dynamic::s_vn_dot_vn(3,ODP.now2startDistance,ODP.xNowInG)-handleLen/2;
                    for(int i=0;i<3;i++)
                    {
                        ODP.now2startPeeDistInG[i]=-handleLen/2*ODP.xNowInG[i];
                    }
                }
                else
                {
                    ODP.now2startDistanceInB[0]=aris::dynamic::s_vn_dot_vn(3,ODP.now2startDistance,ODP.xNowInG)+handleLen/2;
                    for(int i=0;i<3;i++)
                    {
                        ODP.now2startPeeDistInG[i]=handleLen/2*ODP.xNowInG[i];
                    }
                }
                ODP.now2startDistanceInB[1]=aris::dynamic::s_vn_dot_vn(3,ODP.now2startDistance,ODP.yNowInG);
                ODP.now2startDistanceInB[2]=0;

                aris::dynamic::s_pm_dot_v3(*ODP.nowPm,ODP.now2startDistanceInB,ODP.now2startDistanceInG);
            }

            break;

        case MoveState::Pull:
            if (param.count%100==0)
            {
                rt_printf("pushState:%d\n",ODP.pushState);
            }
            switch(ODP.pullState)
            {
            case PullState::pull2Start:
                for (int i=0;i<6;i++)
                {
                    currentPe[i]=ODP.nowPE[i]+(ODP.now2startDistanceInG[i])/2*(1-cos((param.count-ODP.countIter)*PI/(2*ODP.now2StartCount)));
                }
                if(param.count-ODP.countIter<ODP.now2StartCount)
                {
                    for(int i=0;i<3;i++)
                    {
                        //leg 0,2,4
                        currentPee[6*i]=ODP.nowPee[6*i]+ODP.now2startPeeDistInG[0]/2*(1-cos((param.count-ODP.countIter)*PI/ODP.now2StartCount));
                        currentPee[6*i+1]=ODP.nowPee[6*i+1]+0.025*(1-cos((param.count-ODP.countIter)*2*PI/ODP.now2StartCount));
                        currentPee[6*i+2]=ODP.nowPee[6*i+2]+ODP.now2startPeeDistInG[2]/2*(1-cos((param.count-ODP.countIter)*PI/ODP.now2StartCount));
                        //leg 1,3,5
                        currentPee[6*i+3]=ODP.nowPee[6*i+3];
                        currentPee[6*i+4]=ODP.nowPee[6*i+4];
                        currentPee[6*i+5]=ODP.nowPee[6*i+5];
                    }
                }
                else if(param.count-ODP.countIter<2*ODP.now2StartCount)
                {
                    for(int i=0;i<3;i++)
                    {
                        //leg 1,3,5
                        currentPee[6*i+3]=ODP.nowPee[6*i+3]+ODP.now2startPeeDistInG[0]/2*(1-cos((param.count-ODP.countIter-ODP.now2StartCount)*PI/ODP.now2StartCount));
                        currentPee[6*i+4]=ODP.nowPee[6*i+4]+0.025*(1-cos((param.count-ODP.countIter-ODP.now2StartCount)*2*PI/ODP.now2StartCount));
                        currentPee[6*i+5]=ODP.nowPee[6*i+5]+ODP.now2startPeeDistInG[2]/2*(1-cos((param.count-ODP.countIter-ODP.now2StartCount)*PI/ODP.now2StartCount));
                        //leg 0,2,4
                        currentPee[6*i]=ODP.nowPee[6*i]+ODP.now2startPeeDistInG[0];
                        currentPee[6*i+1]=ODP.nowPee[6*i+1];
                        currentPee[6*i+2]=ODP.nowPee[6*i+2]+ODP.now2startPeeDistInG[2];
                    }
                }

                robot.SetPeb(currentPe);
                robot.SetPee(currentPee);

                if(param.count-ODP.countIter==2*ODP.now2StartCount-1)
                {
                    if(isConfirm==true)
                    {
                        ODP.pullState=PullState::circleWalk;
                        ODP.countIter=param.count+1;

                        ODP.circleWalkParam.radius=1.1;
                        ODP.circleWalkParam.n=4;
                        ODP.circleWalkParam.beta=0.1;
                        ODP.circleWalkParam.h=0.05;
                        ODP.circleWalkParam.totalCount=1000;
                        ODP.circleWalkParam.direction=(isLeft==false ? -1 : 1);
                        ODP.circleWalkParam.startAngle=(isLeft==false ? -PI/2 : PI/2);
                    }
                }
                break;

            case PullState::circleWalk:
                ODP.circleWalkParam.count=param.count-ODP.countIter;
                ODP.ret=NormalGait::circleWalk(robot,ODP.circleWalkParam);

                if(ODP.ret==0)
                {
                    if(isConfirm==true)
                    {
                        ODP.pullState=PullState::legWork;
                        ODP.countIter=param.count+1;
                        robot.GetPeb(ODP.nowPE);
                        robot.GetPee(ODP.nowPee,robot.body());
                    }
                }
                break;

            case PullState::legWork:
                memcpy(currentPe,ODP.nowPE,6*sizeof(double));
                memcpy(currentPee,ODP.nowPee,18*sizeof(double));
                if(param.count-ODP.countIter<ODP.legWorkCount)
                {
                    if(isLeft==false)//leg 2
                    {
                        currentPee[9]=ODP.nowPee[9];
                        currentPee[10]=ODP.nowPee[10]+0.025*(1-cos((param.count-ODP.countIter)*2*PI/ODP.legWorkCount));
                        currentPee[11]=ODP.nowPee[11]-0.1/2*(1-cos((param.count-ODP.countIter)*PI/ODP.legWorkCount));
                    }
                    else//leg 0
                    {
                        currentPee[0]=ODP.nowPee[0];
                        currentPee[1]=ODP.nowPee[1]+0.025*(1-cos((param.count-ODP.countIter)*2*PI/ODP.legWorkCount));
                        currentPee[2]=ODP.nowPee[2]-0.1/2*(1-cos((param.count-ODP.countIter)*PI/ODP.legWorkCount));
                    }
                    if(param.count-ODP.countIter==ODP.legWorkCount-1)
                    {
                        memcpy(ODP.nowPee,currentPee,18*sizeof(double));
                    }
                }
                else if(param.count-ODP.countIter<2*ODP.legWorkCount)
                {
                    currentPe[1]=ODP.nowPE[1]+jumpH/2*(1-cos((param.count-ODP.countIter-ODP.legWorkCount)*PI/ODP.legWorkCount));
                    if(param.count-ODP.countIter==2*ODP.legWorkCount-1)
                    {
                        memcpy(ODP.nowPE,currentPe,18*sizeof(double));
                        aris::dynamic::s_pe2pm(ODP.nowPE,*ODP.nowPm,"313");
                    }
                }
                else
                {
                    double backDistInB[3] {0,0,0.1};
                    double backDistInG[3] {0};
                    aris::dynamic::s_pm_dot_v3(*ODP.nowPm,backDistInB,backDistInG);
                    currentPe[0]=ODP.nowPE[0]+backDistInG[0]/2*(1-cos((param.count-ODP.countIter-2*ODP.legWorkCount)*PI/ODP.legWorkCount));
                    currentPe[2]=ODP.nowPE[2]+backDistInG[2]/2*(1-cos((param.count-ODP.countIter-2*ODP.legWorkCount)*PI/ODP.legWorkCount));
                    if(param.count-ODP.countIter==3*ODP.legWorkCount-1)
                    {
                        ODP.pullState=PullState::rotate;
                        ODP.countIter=param.count+1;

                        ODP.circleWalkParam.radius=(isLeft==false ? sqrt(currentPee[9]*currentPee[9]+currentPee[11]*currentPee[11]) : sqrt(currentPee[0]*currentPee[0]+currentPee[2]*currentPee[2]));
                        ODP.circleWalkParam.n=6;
                        ODP.circleWalkParam.beta=PI/2/(ODP.circleWalkParam.n-0.5);
                        ODP.circleWalkParam.h=0.05;
                        ODP.circleWalkParam.totalCount=1000;
                        ODP.circleWalkParam.startAngle=0;
                        ODP.circleWalkParam.direction=(isLeft==false ? 1 : -1);
                    }
                }

                robot.SetPeb(currentPe);
                robot.SetPee(currentPee,robot.body());

                break;

            case PullState::rotate:
                ODP.circleWalkParam.count=param.count-ODP.countIter;
                ODP.ret=NormalGait::circleWalk(robot,ODP.circleWalkParam);

                if(ODP.ret==0)
                {
                    if(isConfirm==true)
                    {
                        ODP.pullState=PullState::pullWalk;
                        ODP.countIter=param.count+1;
                        robot.GetPeb(ODP.nowPE);
                        robot.GetPee(ODP.nowPee,robot.body());
                    }
                }
                break;

            case PullState::pullWalk:

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

        case MoveState::PrePush:

            if ( fabs(bodyVel[0])<1e-10 && fabs(bodyVel[1])<1e-10 && fabs(bodyVel[2])<1e-10
              && fabs(bodyVel[3])<1e-10 && fabs(bodyVel[4])<1e-10 && fabs(bodyVel[5])<1e-10)
            {
                ODP.countIter=param.count+1;
                ODP.moveState=MoveState::Push;
                ODP.pushState=PushState::push2Start;

                //now2start param
                robot.GetPeb(ODP.nowPE);
                robot.GetPmb(*ODP.nowPm);
                robot.GetPee(ODP.nowPee);
                for(int i=0;i<3;i++)
                {
                    ODP.now2startDistance[i]=ODP.startPE[i]-ODP.nowPE[i]; //startPE recorded at the last forward
                }
                aris::dynamic::s_pm_dot_v3(*ODP.nowPm,xBodyInB,ODP.xNowInG);
                aris::dynamic::s_pm_dot_v3(*ODP.nowPm,yBodyInB,ODP.yNowInG);

                ODP.now2startDistanceInB[0]=aris::dynamic::s_vn_dot_vn(3,ODP.now2startDistance,ODP.xNowInG);
                if(isJump==false)
                {
                    ODP.now2startDistanceInB[1]=aris::dynamic::s_vn_dot_vn(3,ODP.now2startDistance,ODP.yNowInG);
                }
                else
                {
                    ODP.now2startDistanceInB[1]=aris::dynamic::s_vn_dot_vn(3,ODP.now2startDistance,ODP.yNowInG);
                }
                ODP.now2startDistanceInB[2]=0;

                aris::dynamic::s_pm_dot_v3(*ODP.nowPm,ODP.now2startDistanceInB,ODP.now2startDistanceInG);
            }

            break;

        case MoveState::Push:
            //Three Step using position control
            if (param.count%100==0)
            {
                rt_printf("pushState:%d\n",ODP.pushState);
            }
            switch(ODP.pushState)
            {
            case PushState::push2Start:
                //1.Move back to startPE;
                for (int i=0;i<6;i++)
                {
                    currentPe[i]=ODP.nowPE[i]+(ODP.now2startDistanceInG[i])/2*(1-cos((param.count-ODP.countIter)*PI/ODP.now2StartCount));
                }

                robot.SetPeb(currentPe);
                robot.SetPee(ODP.nowPee);

                if(param.count-ODP.countIter+1==ODP.now2StartCount)
                {
                    if(isConfirm==true)
                    {
                        ODP.pushState=PushState::alignWalk;
                        ODP.countIter=param.count+1;

                        ODP.walkParam.n=2;
                        ODP.walkParam.alpha=(isLeft==true ? -PI/2 : PI/2);
                        ODP.walkParam.beta=0;
                        ODP.walkParam.totalCount=1000;
                        ODP.walkParam.d=(isJump==false ? (handle2MiddleDist-0.03)/3*2 : (handle2MiddleDist-0.1)/3*2);
                    }
                    else//pause, tested useless
                    {
                        robot.SetPeb(ODP.startPE);
                        robot.SetPee(ODP.nowPee);
                    }
                }

                break;

            case PushState::alignWalk:
                //2.Move rightward to align with the door;
                ODP.walkParam.count=param.count-ODP.countIter;

                ODP.ret=Robots::walkGait(robot,ODP.walkParam);

                if(ODP.ret==0)
                {
                    if(isConfirm==true)
                    {
                        ODP.pushState=PushState::pushWalk;
                        ODP.countIter=param.count+1;

                        ODP.walkParam.totalCount=1500;
                        ODP.walkParam.n=4;
                        ODP.walkParam.alpha=0;
                        ODP.walkParam.beta=0;
                        ODP.walkParam.d=0.4;
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

            case PushState::pushWalk:
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

        if(ODP.moveState!=MoveState::Push && ODP.moveState!=MoveState::LocateAjust && ODP.moveState!=MoveState::Follow
                && ODP.moveState!=MoveState::Jump && ODP.moveState!=MoveState::NonJump)
        {
           if (isConfirm==false && ODP.isPause==false)
            {
                ODP.moveState_last=ODP.moveState;
                ODP.moveState=MoveState::None;
                ODP.isPause=true;
            }
            else if(isConfirm==false && ODP.isPause==true)
            {
                ODP.moveState=MoveState::None;
                ODP.pauseCount++;
                //Do not change F to zeros here, so that we can move the robot manually when paused
            }
            else if(isConfirm==true && ODP.isPause==true)
            {
                ODP.moveState=ODP.moveState_last;
                ODP.isPause=false;
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
            aris::dynamic::s_pe2pm(deltaPE,*deltaPm,"213");
            aris::dynamic::s_pm_dot_pm(*bodyPm,*deltaPm,*realPm);
            aris::dynamic::s_pm2pe(*realPm,realPE,"313");

            robot.GetPee(pEE);
            robot.SetPeb(realPE);
            robot.SetPee(pEE);
        }
        if (param.count%100==0)
        {
            rt_printf("moveState:%d,forceRaw:%.2f,%.2f,%.2f,forceInB:%.2f,%.2f,%.2f,forceFiltered:%.2f,%.2f,%.2f\n",ODP.moveState,param.ruicong_data->at(0).force[0].Fx,param.ruicong_data->at(0).force[0].Fy,param.ruicong_data->at(0).force[0].Fz,ODP.forceInB[0],ODP.forceInB[1],ODP.forceInB[2],ODP.forceInB_filtered[0],ODP.forceInB_filtered[1],ODP.forceInB_filtered[2]);
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
        memcpy(ODP.bodyPE_last,bodyPE,sizeof(double)*6);//used to record data
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
        aris::dynamic::s_pe2pm(deltaPE,*deltaPm,"213");
        aris::dynamic::s_pm_dot_pm(*bodyPm,*deltaPm,*realPm);
        aris::dynamic::s_pm2pe(*realPm,realPE,"313");
        double nowPee[18];
        robot.GetPee(nowPee);
        robot.SetPeb(realPE);
        robot.SetPee(nowPee);

        memcpy(ODP.pEE_last,nowPee,sizeof(double)*18);
        memcpy(ODP.bodyPE_last,bodyPE,sizeof(double)*6);//used to record data
        memcpy(ODP.bodyVel_last,bodyVel,sizeof(double)*6);

        if ( fabs(bodyVel[0])<1e-10 && fabs(bodyVel[1])<1e-10 && fabs(bodyVel[2])<1e-10
          && fabs(bodyVel[3])<1e-10 && fabs(bodyVel[4])<1e-10 && fabs(bodyVel[5])<1e-10)
            return 0;
        else
            return 1;
    }
}

void OpenDoor::initialize(Robots::RobotBase &robot)
{
    std::fill_n(ODP.bodyVel_last,6,0);
    ODP.isPause=false;
    ODP.moveState_last=MoveState::None;
    ODP.moveState=MoveState::None;
    ODP.isLastFollow=false;
    robot.GetPeb(ODBeginPeb);
    robot.GetPmb(*ODBeginPmb);
    aris::dynamic::s_pe2pe("313",ODBeginPeb,"213",ODP.bodyPE_last);

    robot.GetPeb(adjustBeginPeb);
    ODP.filter.Initialize();
    ODP.filter.SetCutFrequency(0.03,1000);
}

int OpenDoor::locateDoor(Robots::RobotBase &robot, int count)
{
    double point1[6] {0};
    double point2[6] {0};
    double point3[6] {0};
    static double Location[3][3] {0};
    double invLocation[3][3] {0};
    double planeConst[3]{1,1,1};
    double planeVertical[3] {0};

    std::fill_n(Fbody,6,0);
    switch(DLState)
    {
    case DoorLocateState::Point1:
        Fbody[2]=-1;
        GenerateBodyMotionByFce(robot);
        if (fabs(ODP.forceInB_filtered[2])>ForceRange[0])
        {
            ODP.countIter=count+1;
            robot.GetPeb(point1);
            memcpy(*Location,point1,sizeof(double)*3);
            DLState=DoorLocateState::Point2;
        }
        break;
    case DoorLocateState::Point2:
        if (count-ODP.countIter<2500)
        {
            Fbody[2]=1;
            if(isLeft==true)
                Fbody[0]=1;
            else
                Fbody[0]=-1;
        }
        else
        {
            Fbody[2]=-1;
        }
        GenerateBodyMotionByFce(robot);
        if (fabs(ODP.forceInB_filtered[2])>ForceRange[0] && (count-ODP.countIter)>=2500)
        {
            ODP.countIter=count+1;
            robot.GetPeb(point2);
            memcpy(*Location+3,point2,sizeof(double)*3);
            DLState=DoorLocateState::Point3;
        }
        break;
    case DoorLocateState::Point3:
        if (count-ODP.countIter<1500)
        {
            Fbody[2]=1;
            Fbody[1]=1;
        }
        else
        {
            Fbody[2]=-1;
        }
        GenerateBodyMotionByFce(robot);
        if (fabs(ODP.forceInB_filtered[2])>ForceRange[0] && (count-ODP.countIter)>=1500)
        {
            ODP.countIter=count+1;
            robot.GetPeb(point3);
            memcpy(*Location+6,point3,sizeof(double)*3);

            //calculate the plane of the door. ax+by+cz=1,(a,b,c) is the vertical vector of the plane
            NormalGait::inv3(*Location,*invLocation);
            aris::dynamic::s_dgemm(3,1,3,1,*invLocation,3,planeConst,1,1,planeVertical,1);
            aris::dynamic::s_inv_pm_dot_v3(*ODBeginPmb,planeVertical,planeVerticalInB);

            return 0;
        }
        break;
    default:
        break;
    }

    return 1;
}

void OpenDoor::GetAdjustParam(Robots::RobotBase &robot)
{
    double rotateRatio {1};//0.88

    ODP.planeYPR[0]=atan(planeVerticalInB[0]/planeVerticalInB[2]);
    ODP.planeYPR[1]=-asin(planeVerticalInB[1]/NormalGait::norm(planeVerticalInB));
    ODP.planeYPR[2]=0;

    //Set now param of now2start in ajustRobot
    robot.GetPeb(nowPeb);
    robot.GetPee(nowPee);

    //Set rotate param
    if(fabs(ODP.planeYPR[0])>=PI/4)
    {
        rt_printf("WARNING!!! Too large angle between the door and the robot!");
    }
    else if(fabs(ODP.planeYPR[0])>PI/9)
    {
        ODP.walkParam.n=2;
        ODP.walkParam.beta=ODP.planeYPR[0]*rotateRatio/3*2;
    }
    else
    {
        ODP.walkParam.n=1;
        ODP.walkParam.beta=ODP.planeYPR[0]*rotateRatio*2;
    }
    ODP.walkParam.d=0;
    ODP.walkParam.alpha=0;
    ODP.walkParam.totalCount=1000;

    rt_printf("yaw:%f,pitch:%f,roll:%f\n",ODP.planeYPR[0],ODP.planeYPR[1],ODP.planeYPR[2]);
}

void OpenDoor::GenerateBodyMotionByFce(Robots::RobotBase &robot)
{
    double C[6] {50,50,50,50,50,50};
    double M[6] {1,1,1,1,1,1};
    double deltaT {0.001};
    double bodyPm[4][4] {0};
    double deltaPm[4][4] {0};
    double deltaPE[6] {0};
    double realPm[4][4] {0};

    for (int i=0;i<6;i++)
    {
        bodyAcc[i]=(Fbody[i]-C[i]*ODP.bodyVel_last[i])/M[i];
        bodyVel[i]=ODP.bodyVel_last[i]+bodyAcc[i]*deltaT;
        deltaPE[i]=bodyVel[i]*deltaT;
    }

    robot.GetPmb(*bodyPm);
    aris::dynamic::s_pe2pm(deltaPE,*deltaPm,"213");
    aris::dynamic::s_pm_dot_pm(*bodyPm,*deltaPm,*realPm);
    robot.SetPmb(*realPm);
}

int OpenDoor::adjustRobotD2H(Robots::RobotBase &robot, int count)
{
    double currentPe[6] {0};

    ODP.walkParam.count=count-ODP.countIter;
    int ret=Robots::walkGait(robot,ODP.walkParam);
    robot.GetPeb(currentPe);

    for (int i=0;i<3;i++)
    {
        currentPe[i]=nowPeb[i]+(adjustBeginPeb[i]-nowPeb[i])/2*(1-cos((param.count-ODP.countIter)*PI/(2*ODP.walkParam.n*ODP.walkParam.totalCount)));
    }
    robot.SetPeb(currentPe);

    if(ret==0)
    {
        ODP.countIter=count+1;
        std::fill_n(ODP.bodyVel_last,6,0);
    }

    return ret;
}

int OpenDoor::locateHandle(Robots::RobotBase &robot, int count)
{
    static double beginPeb[6] {0};
    static double handlePebRefer[6] {0};
    static double vector0[3] {0};
    static double vector1[3] {0};
    static double vector2[3] {0};
    double currentPeb[6] {0};
    double currentPee[18] {0};
    std::fill_n(Fbody,6,0);
    switch(HLState)
    {
    case HandleLocateState::Forward:
        if(count==ODP.countIter)
        {
            robot.GetPeb(beginPeb);//To calculate vector0
            robot.GetPeb(handlePosRefer);
        }
        Fbody[2]=-1;
        GenerateBodyMotionByFce(robot);

        if(fabs(ODP.forceInB_filtered[2])>ForceRange[0])
        {
            ODP.countIter=count+1;
            HLState=HandleLocateState::Backward;
            robot.GetPeb(nowPeb);
            for (int i=0;i<3;i++)
            {
                vector0[i]=nowPeb[i]-beginPeb[i];
            }
        }
        break;
    case HandleLocateState::Backward:
        Fbody[2]=1;
        GenerateBodyMotionByFce(robot);

        if (param.count-ODP.countIter>=1000)
        {
            ODP.countIter=count+1;
            robot.GetPee(ODP.endPeeInB,robot.body());
            if(isLeft==false)
                HLState=HandleLocateState::Rightward;
            else
                HLState=HandleLocateState::Leftward;
        }
        break;
    case HandleLocateState::Leftward:
        if(count==ODP.countIter)
        {
            robot.GetPeb(beginPeb);//To calculate vector1
        }
        Fbody[0]=-1;
        GenerateBodyMotionByFce(robot);

        if (fabs(ODP.forceInB_filtered[0])>ForceRange[0] || count-ODP.countIter>3000)
        {
            ODP.countIter=count+1;
            robot.GetPee(ODP.startPeeInB,robot.body());
            robot.GetPeb(nowPeb);
            HLState=HandleLocateState::Follow;
            if(fabs(ODP.forceInB_filtered[0])>ForceRange[0])
            {
                ODP.isLastFollow=true;
                for (int i=0;i<3;i++)
                {
                    vector1[i]=nowPeb[i]-beginPeb[i];
                }
            }
        }
        break;
    case HandleLocateState::RightWard:
        if(count==ODP.countIter)
        {
            robot.GetPeb(beginPeb);//To calculate vector1
        }
        Fbody[0]=1;
        GenerateBodyMotionByFce(robot);

        if (fabs(ODP.forceInB_filtered[0])>ForceRange[0] || count-ODP.countIter>3000)
        {
            ODP.countIter=count+1;
            robot.GetPee(ODP.startPeeInB,robot.body());
            robot.GetPeb(nowPeb);
            HLState=HandleLocateState::Follow;
            if(fabs(ODP.forceInB_filtered[0])>ForceRange[0])
            {
                ODP.isLastFollow=true;
                for (int i=0;i<3;i++)
                {
                    vector1[i]=nowPeb[i]-beginPeb[i];
                }
            }
        }
        break;
    case HandleLocateState::Follow:
        const int followCount {1000};
        if(count-ODP.countIter<followCount)
        {
            for(int i=0;i<3;i++)
            {
                //leg 0,2,4
                currentPee[6*i]=ODP.startPeeInB[6*i]+(ODP.endPeeInB[6*i]-ODP.startPeeInB[6*i])/2*(1-cos((count-ODP.countIter)*PI/followCount));
                currentPee[6*i+1]=ODP.startPeeInB[6*i+1]+0.025*(1-cos((count-ODP.countIter)*2*PI/followCount));
                currentPee[6*i+2]=ODP.startPeeInB[6*i+2];
                //leg 1,3,5
                currentPee[6*i+3]=ODP.startPeeInB[6*i+3];
                currentPee[6*i+4]=ODP.startPeeInB[6*i+4];
                currentPee[6*i+5]=ODP.startPeeInB[6*i+5];
            }
        }
        else if(count-ODP.countIter<2*followCount)
        {
            for(int i=0;i<3;i++)
            {
                //leg 1,3,5
                currentPee[6*i+3]=ODP.startPeeInB[6*i+3]+(ODP.endPeeInB[6*i+3]-ODP.startPeeInB[6*i+3])/2*(1-cos((count-ODP.countIter-followCount)*PI/followCount));
                currentPee[6*i+4]=ODP.startPeeInB[6*i+4]+0.025*(1-cos((count-ODP.countIter-followCount)*2*PI/followCount));
                currentPee[6*i+5]=ODP.startPeeInB[6*i+5];
                //leg 0,2,4
                currentPee[6*i]=ODP.endPeeInB[6*i];
                currentPee[6*i+1]=ODP.endPeeInB[6*i+1];
                currentPee[6*i+2]=ODP.endPeeInB[6*i+2];
            }
        }
        robot.SetPeb(nowPeb);
        robot.SetPee(currentPee,robot.body());

        if(count-ODP.countIter+1==2*followCount)
        {
            ODP.countIter=count+1;
            if(ODP.isLastFollow==false)
            {
                HLState=HandleLocateState::Forward;
            }
            else
            {
                robot.GetPeb(jumpBeginPeb);
                robot.GetPmb(*nowPmb);
                robot.GetPee(nowPee);
                if(isJump==false)
                {
                    HLState=HandleLocateState::NonJump;
                    robot.GetPeb(adjustBeginPeb);//For adjustRobotT2P
                }
                else
                {
                    HLState=HandleLocateState::Jump;
                }
            }
        }
        break;
    case HandleLocateState::Jump:
        double jumpH {0.04};
        double jumpDist {0.05};
        const int jumpCount {1000};
        double jumpDistInB[3] {0};
        double jumpDistInG[3] {0};

        jumpDistInB[0]=(isLeft==true ? -jumpDist : jumpDist);
        aris::dynamic::s_pm_dot_v3(*nowPmb,jumpDistInB,jumpDistInG);

        if(count-ODP.countIter<jumpCount)
        {
            memcpy(currentPeb,jumpBeginPeb,6*sizeof(double));
            memcpy(currentPee,nowPee,18*sizeof(double));
            currentPeb[1]=jumpBeginPeb[1]+jumpH/2*(1-cos((count-ODP.countIter)*PI/jumpCount));
            for(int i=0;i<3;i++)
            {
                currentPee[6*i]=nowPee[6*i]+jumpDistInG[0]/2*(1-cos((count-ODP.countIter)*PI/jumpCount));
                currentPee[6*i+1]=nowPee[6*i+1]+0.025*(1-cos((count-ODP.countIter)*2*PI/ODP.followCount));
                currentPee[6*i+2]=nowPee[6*i+2]+jumpDistInG[2]/2*(1-cos((count-ODP.countIter)*PI/jumpCount));
            }
            if(count-ODP.countIter==jumpCount-1)
            {
                jumpBeginPeb[1]=currentPeb[1];
            }
        }
        else if(count-ODP.countIter<2*jumpCount)
        {
            memcpy(currentPeb,jumpBeginPeb,6*sizeof(double));
            currentPeb[0]=jumpBeginPeb[0]+jumpDistInG[0]/2*(1-cos((count-jumpCount-ODP.countIter)*PI/jumpCount));
            currentPeb[2]=jumpBeginPeb[2]+jumpDistInG[2]/2*(1-cos((count-jumpCount-ODP.countIter)*PI/jumpCount));
            for(int i=0;i<3;i++)
            {
                currentPee[6*i]=nowPee[6*i]+jumpDistInG[0];
                currentPee[6*i+1]=nowPee[6*i+1];
                currentPee[6*i+2]=nowPee[6*i+2]+jumpDistInG[2];

                currentPee[6*i+3]=nowPee[6*i+3]+jumpDistInG[0]/2*(1-cos((count-jumpCount-ODP.countIter)*PI/jumpCount));
                currentPee[6*i+4]=nowPee[6*i+4]+0.025*(1-cos((count-jumpCount-ODP.countIter)*2*PI/jumpCount));
                currentPee[6*i+5]=nowPee[6*i+5]+jumpDistInG[2]/2*(1-cos((count-jumpCount-ODP.countIter)*PI/jumpCount));
            }
            if(count-ODP.countIter==2*jumpCount-1)
            {
                HLState=HandleLocateState::Downward;
                robot.GetPeb(adjustBeginPeb);//For adjustRobotT2P
            }
        }
        robot.SetPeb(currentPeb);
        robot.SetPee(currentPee);
        break;
    case HandleLocateState::NonJump:
        double nonJumpDist=0.05;
        const int nonJumpCount {1000};
        double nonJumpDistInB[3] {0};
        double nonJumpDistInG[3] {0};
        nonJumpDistInB[0]=(isLeft==true ? nonJumpDist : -nonJumpDist);
        aris::dynamic::s_pm_dot_v3(*nowPmb,nonJumpDistInB,nonJumpDistInG);

        memcpy(currentPeb,jumpBeginPeb,6*sizeof(double));
        currentPeb[0]=jumpBeginPeb[0]+nonJumpDistInG[0]/2*(1-cos((count-ODP.countIter)*PI/nonJumpCount));
        currentPeb[2]=jumpBeginPeb[2]+nonJumpDistInG[2]/2*(1-cos((count-ODP.countIter)*PI/nonJumpCount));

        robot.SetPeb(currentPeb);
        robot.SetPee(nowPee);

        if (count-ODP.countIter==nonJumpCount-1)
        {
            HLState=HandleLocateState::Downward;
        }
        break;
    case HandleLocateState::Downward:
        if(count==ODP.countIter)
        {
            robot.GetPeb(beginPeb);//To calculate vector2
        }
        Fbody[1]=-1;
        GenerateBodyMotionByFce(robot);

        if(fabs(ODP.forceInB_filtered[1])>ForceRange[0])
        {
            ODP.countIter=count+1;
            robot.GetPeb(nowPeb);
            for (int i=0;i<3;i++)
            {
                vector2[i]=nowPeb[i]-beginPeb[i];
            }

            double handleLocation[3];
            handleLocation[2]=-NormalGait::norm(vector0);
            handleLocation[0]=NormalGait::norm(vector1);
            handleLocation[1]=-NormalGait::norm(vector2);
            rt_printf("HandleLocate Finished! (%.4f,%.4f,%.4f,0,0,0) to PebRefer (%.4f,%.4f,%.4f,%.4f,%.4f,%.4f).\n",
                      handleLocation[0],handleLocation[1],handleLocation[2],handlePebRefer[0],handlePebRefer[1],
                      handlePebRefer[2],handlePebRefer[3],handlePebRefer[4],handlePebRefer[5]);

            return 0;
        }
        break;
    default:
        break;
    }

    return 1;
}

int OpenDoor::turnHandle(Robots::RobotBase &robot, int count)
{
    const int downwardCount=(int)(PI/2*2500);
    if(isLeft==false)
    {
        Fbody[1]=-cos((count-ODP.countIter)*PI/downwardCount/2);
        Fbody[0]=sin((count-ODP.countIter)*PI/downwardCount/2);
    }
    else
    {
        Fbody[1]=-cos((count-ODP.countIter)*PI/downwardCount/2);
        Fbody[0]=-sin((count-ODP.countIter)*PI/downwardCount/2);
    }
    GenerateBodyMotionByFce(robot);

    if(fabs(ODP.forceInB_filtered[0])>ForceRange[1] || fabs(ODP.forceInB_filtered[1])>ForceRange[1])
    {
        ODP.countIter=count+1;

        return 0;
    }

    return 1;
}

void recordOpenDoorData()
{
    openDoorThread = std::thread([&]()
    {
        struct OpenDoorParam param;
        static std::fstream fileGait;
        std::string name = aris::core::logFileName();
        name.replace(name.rfind("log.txt"), std::strlen("openDoorData.txt"), "openDoorData.txt");
        fileGait.open(name.c_str(), std::ios::out | std::ios::trunc);

        long long count = -1;
        while (1)
        {
            openDoorPipe.recvInNrt(param);

            //fileGait << ++count << " ";
            fileGait << param.count << "  ";
            fileGait << param.moveState << "  ";
            fileGait << param.pushState << "  ";
            for (int i = 0; i < 18; i++)
            {
                fileGait << param.pEE_last[i] << "  ";
            }
            for (int i = 0; i < 18; i++)
            {
                fileGait << param.bodyPE_last[i] << "  ";
            }
            for (int i = 0; i < 18; i++)
            {
                fileGait << param.bodyVel_last[i] << "  ";
            }
            for (int i=0;i<6;i++)
            {
                fileGait << param.forceInB[i] << "  ";
            }
            fileGait << std::endl;
        }

        fileGait.close();
    });
}
