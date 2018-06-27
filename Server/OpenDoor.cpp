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

//Pipe<OpenDoorParam> openDoorPipe(true);
double OpenDoor::Fbody[6];
double OpenDoor::ForceRange[2] {1,6};
double OpenDoor::bodyLastVel[6];//derivative of 213 euler angle, not spatial velocity
double OpenDoor::planeVerticalInB[3];
double OpenDoor::ODBeginPeb[6];
double OpenDoor::ODBeginPmb[4][4];
double OpenDoor::adjustBeginPeb[6];
double OpenDoor::jumpBeginPeb[6];
double OpenDoor::nowPeb[6];
double OpenDoor::nowPmb[4][4];
double OpenDoor::nowPee[18];
double OpenDoor::now2startPebDistInG[6];
double OpenDoor::now2startPeeDistInG[3];

int OpenDoor::moveDir[6];
bool OpenDoor::isPull;
bool OpenDoor::isLeft;
bool OpenDoor::isJump;

OpenDoorState OpenDoor::ODState;
DoorLocateState OpenDoor::DLState;
HandleLocateState OpenDoor::HLState;
PullWorkState OpenDoor::PWState;
WalkThroughState OpenDoor::WTState;

OpenDoorParam OpenDoor::ODP;
Robots::WalkParam OpenDoor::walkParam;
NormalGait::CircleWalkParam OpenDoor::circleWalkParam;
LowpassFilter<6> OpenDoor::filter;
Pipe<OpenDoorParam> OpenDoor::openDoorPipe;
std::thread OpenDoor::openDoorThread;

OpenDoor::OpenDoor(){}
OpenDoor::~OpenDoor(){}

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

    std::cout<<"finished parse"<<std::endl;
}

//*****only for Robot VIII*****
int OpenDoor::openDoor(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
{
    auto &robot = static_cast<Robots::RobotBase &>(model);
    auto &param = static_cast<const aris::server::GaitParamBase &>(param_in);

    double forceInF[6];
    ForceTask::forceInit(param.count,param.ruicong_data->at(0).force[0].fce,forceInF);
    aris::dynamic::s_f2f(*robot.forceSensorMak().prtPm(),forceInF,ODP.forceInB);
    filter.DoFilter(ODP.forceInB,ODP.forceInB_filtered);

    EmergencyStop stoper;
    stoper.doRecord(robot);

    ODP.count=param.count;
    std::fill_n(Fbody,6,0);
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
            GetAdjustD2HParam(robot);
            ODState=OpenDoorState::AdjustD2H;
            rt_printf("ODState changed to AdjustD2H at count %d\n",param.count);
        }
        break;
    case OpenDoorState::AdjustD2H:
        ret=adjustRobotD2H(robot,param.count);
        if(ret==0)
        {
            ODState=OpenDoorState::HandleLocate;
            HLState=HandleLocateState::Forward;
            rt_printf("ODState changed to HandleLocate at count %d\n",param.count);
        }
        break;
    case OpenDoorState::HandleLocate:
        ret=locateHandle(robot,param.count);
        if(ret==0)
        {
            ODState=OpenDoorState::TurnHandle;
            rt_printf("ODState changed to TurnHandle at count %d\n",param.count);
        }
        break;
    case OpenDoorState::TurnHandle:
        ret=turnHandle(robot,param.count);
        if(ret==0)
        {
            if(isPull==true)
            {
                ODState=OpenDoorState::PullHandle;
                rt_printf("ODState changed to PullHandle at count %d\n",param.count);
            }
            else
            {
                ODState=OpenDoorState::PushHandle;
                rt_printf("ODState changed to PushHandle at count %d\n",param.count);
            }
        }
        break;
    case OpenDoorState::PushHandle:
        ret=pushHandle(robot,param.count);
        if(ret==0)
        {
            ODState=OpenDoorState::AdjustPush2W;
            GetAdjustPush2WParam(robot);
            rt_printf("ODState changed to AdjustPush2W at count %d\n",param.count);
        }
        break;
    case OpenDoorState::PullHandle:
        ret=pullHandle(robot,param.count);
        if(ret==0)
        {
            ODState=OpenDoorState::AdjustPull2W;
            GetAdjustPull2WParam(robot);
            rt_printf("ODState changed to AdjustPull2W at count %d\n",param.count);
        }
        break;
    case OpenDoorState::AdjustPush2W:
        ret=adjustRobotPush2W(robot,param.count);
        if(ret==0)
        {
            ODState=OpenDoorState::WalkThrough;
            WTState=WalkThroughState::AlignWalk;
            GetPushWalkParam();
            rt_printf("ODState changed to WalkThrough at count %d\n",param.count);
        }
        break;
    case OpenDoorState::AdjustPull2W:
        ret=adjustRobotPull2W(robot,param.count);
        if(ret==0)
        {
            return 0;
            ODState=OpenDoorState::PullWork;
            PWState=PullWorkState::CircleWalk;
            GetPullWalkParam();
            rt_printf("ODState changed to PullWork at count %d\n",param.count);
        }
        break;
    case OpenDoorState::WalkThrough:
        ret=pushWalkThrough(robot,param.count);
        if(ret==0)
        {
            rt_printf("Push Door Finished at count %d\n", param.count);
            return 0;
        }
        break;
    case OpenDoorState::PullWork:
        ret=pullDoorWork(robot,param.count);
        if(ret==0)\
        {
            rt_printf("Pull Door Finished at count %d\n", param.count);
            return 0;
        }
        break;
    case OpenDoorState::Quit:
        ret=stoper.stop(robot,param.count,3.2);
        return ret;
        break;
    default:
        break;
    }

    if (param.count%100==0)
    {
        rt_printf("ODState:%d,DLState:%d,HLState:%d,WTState:%d\nforceRaw:%.2f,%.2f,%.2f,forceInB:%.2f,%.2f,%.2f,forceFiltered:%.2f,%.2f,%.2f\n",
                  ODState,DLState,HLState,WTState,param.ruicong_data->at(0).force[0].Fx,param.ruicong_data->at(0).force[0].Fy,param.ruicong_data->at(0).force[0].Fz,
                ODP.forceInB[0],ODP.forceInB[1],ODP.forceInB[2],ODP.forceInB_filtered[0],ODP.forceInB_filtered[1],ODP.forceInB_filtered[2]);
    }

    openDoorPipe.sendToNrt(ODP);

    return 1;
}

void OpenDoor::initialize(Robots::RobotBase &robot)
{
    std::fill_n(bodyLastVel,6,0);
    robot.GetPmb(*ODBeginPmb);
    robot.GetPeb(adjustBeginPeb);
    filter.Initialize();
    filter.SetCutFrequency(0.03,1000);
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
            rt_printf("DLState changed to Point2 at count %d\n",count);
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
            rt_printf("DLState changed to Point3 at count %d\n",count);
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

            rt_printf("DoorLocate Finished. Going to adjustRobotD2H\n");

            return 0;
        }
        break;
    default:
        break;
    }

    return 1;
}

void OpenDoor::GenerateBodyMotionByFce(Robots::RobotBase &robot)
{
    double C[6] {50,50,50,50,50,50};
    double M[6] {1,1,1,1,1,1};
    double deltaT {0.001};
    double bodyPm[4][4] {0};
    double deltaPm[4][4] {0};
    double deltaPe[6] {0};
    double nextPm[4][4] {0};
    double Pee[18] {0};

    double bodyAcc[6] {0};
    double bodyVel[6] {0};

    for (int i=0;i<6;i++)
    {
        bodyAcc[i]=(Fbody[i]-C[i]*bodyLastVel[i])/M[i];
        bodyVel[i]=bodyLastVel[i]+bodyAcc[i]*deltaT;
        deltaPe[i]=bodyVel[i]*deltaT;
    }

    robot.GetPmb(*bodyPm);
    robot.GetPee(Pee);
    aris::dynamic::s_pe2pm(deltaPe,*deltaPm,"213");
    aris::dynamic::s_pm_dot_pm(*bodyPm,*deltaPm,*nextPm);
    robot.SetPmb(*nextPm);
    robot.SetPee(Pee);

    memcpy(bodyLastVel,bodyVel,6*sizeof(double));
}

void OpenDoor::GetAdjustD2HParam(Robots::RobotBase &robot)
{
    double rotateRatio {1};//0.88
    double planeYPR[3] {0};

    planeYPR[0]=atan(planeVerticalInB[0]/planeVerticalInB[2]);
    planeYPR[1]=-asin(planeVerticalInB[1]/NormalGait::norm(planeVerticalInB));
    planeYPR[2]=0;

    //Set now param of now2start in ajustRobot
    robot.GetPeb(nowPeb);
    robot.GetPee(nowPee);

    //Set rotate param
    if(fabs(planeYPR[0])>=PI/4)
    {
        rt_printf("WARNING!!! Too large angle between the door and the robot!");
    }
    else if(fabs(planeYPR[0])>PI/9)
    {
        walkParam.n=2;
        walkParam.beta=planeYPR[0]*rotateRatio/3*2;
    }
    else
    {
        walkParam.n=1;
        walkParam.beta=planeYPR[0]*rotateRatio*2;
    }
    walkParam.d=0;
    walkParam.alpha=0;
    walkParam.totalCount=1000;

    rt_printf("yaw:%f,pitch:%f,roll:%f\n",planeYPR[0],planeYPR[1],planeYPR[2]);
}

int OpenDoor::adjustRobotD2H(Robots::RobotBase &robot, int count)
{
    double currentPeb[6] {0};
    double currentPee[18] {0};

    walkParam.count=count-ODP.countIter;
    int ret=Robots::walkGait(robot,walkParam);
    robot.GetPeb(currentPeb);
    robot.GetPee(currentPee);

    for (int i=0;i<3;i++)
    {
        currentPeb[i]=nowPeb[i]+(adjustBeginPeb[i]-nowPeb[i])/2*(1-cos((count-ODP.countIter)*PI/(2*walkParam.n*walkParam.totalCount)));
    }
    robot.SetPeb(currentPeb);
    robot.SetPee(currentPee);

    if(ret==0)
    {
        ODP.countIter=count+1;
        std::fill_n(bodyLastVel,6,0);
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
    static double startPeeInB[18] {0};
    static double endPeeInB[18] {0};
    static bool isLastFollow {false};
    const int followCount {1000};
    double currentPeb[6] {0};
    double currentPee[18] {0};
    double jumpH {0.04};
    double jumpDist {0.05};
    const int jumpCount {1000};
    double jumpDistInB[3] {0};
    double jumpDistInG[3] {0};
    double nonJumpDist=0.05;
    const int nonJumpCount {1000};
    double nonJumpDistInB[3] {0};
    double nonJumpDistInG[3] {0};

    switch(HLState)
    {
    case HandleLocateState::Forward:
        if(count==ODP.countIter)
        {
            robot.GetPeb(beginPeb);//To calculate vector0
            robot.GetPeb(handlePebRefer);
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

        if (count-ODP.countIter>=1000)
        {
            ODP.countIter=count+1;
            robot.GetPee(endPeeInB,robot.body());
            if(isLeft==false)
                HLState=HandleLocateState::RightWard;
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
            robot.GetPee(startPeeInB,robot.body());
            robot.GetPeb(nowPeb);
            HLState=HandleLocateState::Follow;
            if(fabs(ODP.forceInB_filtered[0])>ForceRange[0])
            {
                isLastFollow=true;
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
            robot.GetPee(startPeeInB,robot.body());
            robot.GetPeb(nowPeb);
            HLState=HandleLocateState::Follow;
            if(fabs(ODP.forceInB_filtered[0])>ForceRange[0])
            {
                isLastFollow=true;
                for (int i=0;i<3;i++)
                {
                    vector1[i]=nowPeb[i]-beginPeb[i];
                }
            }
        }
        break;
    case HandleLocateState::Follow:
        if(count-ODP.countIter<followCount)
        {
            for(int i=0;i<3;i++)
            {
                //leg 0,2,4
                currentPee[6*i]=startPeeInB[6*i]+(endPeeInB[6*i]-startPeeInB[6*i])/2*(1-cos((count-ODP.countIter)*PI/followCount));
                currentPee[6*i+1]=startPeeInB[6*i+1]+0.025*(1-cos((count-ODP.countIter)*2*PI/followCount));
                currentPee[6*i+2]=startPeeInB[6*i+2];
                //leg 1,3,5
                currentPee[6*i+3]=startPeeInB[6*i+3];
                currentPee[6*i+4]=startPeeInB[6*i+4];
                currentPee[6*i+5]=startPeeInB[6*i+5];
            }
        }
        else if(count-ODP.countIter<2*followCount)
        {
            for(int i=0;i<3;i++)
            {
                //leg 1,3,5
                currentPee[6*i+3]=startPeeInB[6*i+3]+(endPeeInB[6*i+3]-startPeeInB[6*i+3])/2*(1-cos((count-ODP.countIter-followCount)*PI/followCount));
                currentPee[6*i+4]=startPeeInB[6*i+4]+0.025*(1-cos((count-ODP.countIter-followCount)*2*PI/followCount));
                currentPee[6*i+5]=startPeeInB[6*i+5];
                //leg 0,2,4
                currentPee[6*i]=endPeeInB[6*i];
                currentPee[6*i+1]=endPeeInB[6*i+1];
                currentPee[6*i+2]=endPeeInB[6*i+2];
            }
        }
        robot.SetPeb(nowPeb);
        robot.SetPee(currentPee,robot.body());

        if(count-ODP.countIter+1==2*followCount)
        {
            ODP.countIter=count+1;
            if(isLastFollow==false)
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
                    robot.GetPeb(adjustBeginPeb);//For adjustRobotP2W
                }
                else
                {
                    HLState=HandleLocateState::Jump;
                }
            }
        }
        break;
    case HandleLocateState::Jump:
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
                currentPee[6*i+1]=nowPee[6*i+1]+0.025*(1-cos((count-ODP.countIter)*2*PI/followCount));
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
                std::fill_n(bodyLastVel,6,0);
                robot.GetPeb(adjustBeginPeb);//For adjustRobotP2W
            }
        }
        robot.SetPeb(currentPeb);
        robot.SetPee(currentPee);
        break;
    case HandleLocateState::NonJump:
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
            std::fill_n(bodyLastVel,6,0);
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
            isLastFollow=false;
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

int OpenDoor::pushHandle(Robots::RobotBase &robot, int count)
{
    if(count-ODP.countIter<=2500)
    {
        Fbody[2]=-1;
    }
    GenerateBodyMotionByFce(robot);
    if (count-ODP.countIter>2500 && fabs(bodyLastVel[2])<1e-10)
    {
        ODP.countIter=count+1;
        return 0;
    }
    return 1;
}

int OpenDoor::pullHandle(Robots::RobotBase &robot, int count)
{
    if(count-ODP.countIter<=2000)
    {
        Fbody[2]=1;
    }
    GenerateBodyMotionByFce(robot);
    if(count-ODP.countIter>2000 && fabs(bodyLastVel[2])<1e-10)
    {
        ODP.countIter=count+1;
        return 0;
    }
    return 1;
}

void OpenDoor::GetAdjustPush2WParam(Robots::RobotBase &robot)
{
    const double xBodyInB[3] {1,0,0};
    const double yBodyInB[3] {0,1,0};
    double now2startDistance[3] {0};
    double now2startPebDistInB[6] {0};
    double xNowInG[3] {0};
    double yNowInG[3] {0};

    robot.GetPeb(nowPeb);
    robot.GetPmb(*nowPmb);
    robot.GetPee(nowPee);
    for(int i=0;i<3;i++)
    {
        now2startDistance[i]=adjustBeginPeb[i]-nowPeb[i];
    }
    aris::dynamic::s_pm_dot_v3(*nowPmb,xBodyInB,xNowInG);
    aris::dynamic::s_pm_dot_v3(*nowPmb,yBodyInB,yNowInG);

    now2startPebDistInB[0]=aris::dynamic::s_vn_dot_vn(3,now2startDistance,xNowInG);
    now2startPebDistInB[1]=aris::dynamic::s_vn_dot_vn(3,now2startDistance,yNowInG);
    now2startPebDistInB[2]=0;
    now2startPebDistInB[3]=nowPeb[3];
    now2startPebDistInB[4]=nowPeb[4];
    now2startPebDistInB[5]=nowPeb[5];

    aris::dynamic::s_pm_dot_v3(*nowPmb,now2startPebDistInB,now2startPebDistInG);
}

void OpenDoor::GetAdjustPull2WParam(Robots::RobotBase &robot)
{
    const double xBodyInB[3] {1,0,0};
    const double yBodyInB[3] {0,1,0};
    const double handleLen {0.1};
    double now2startDistance[3] {0};
    double now2startPebDistInB[6] {0};
    double xNowInG[3] {0};
    double yNowInG[3] {0};

    robot.GetPeb(nowPeb);
    robot.GetPmb(*nowPmb);
    robot.GetPee(nowPee);
    for(int i=0;i<3;i++)
    {
        now2startDistance[i]=adjustBeginPeb[i]-nowPeb[i];
    }

    aris::dynamic::s_pm_dot_v3(*nowPmb,xBodyInB,xNowInG);
    aris::dynamic::s_pm_dot_v3(*nowPmb,yBodyInB,yNowInG);

    if((isJump==false && isLeft==false) || (isJump==true && isLeft==true))
    {
        now2startPebDistInB[0]=aris::dynamic::s_vn_dot_vn(3,now2startDistance,xNowInG)-handleLen/2;
        for(int i=0;i<3;i++)
        {
            now2startPeeDistInG[i]=-handleLen/2*xNowInG[i];
        }
    }
    else
    {
        now2startPebDistInB[0]=aris::dynamic::s_vn_dot_vn(3,now2startDistance,xNowInG)+handleLen/2;
        for(int i=0;i<3;i++)
        {
            now2startPeeDistInG[i]=handleLen/2*xNowInG[i];
        }
    }
    now2startPebDistInB[1]=aris::dynamic::s_vn_dot_vn(3,now2startDistance,yNowInG);
    now2startPebDistInB[2]=0;

    aris::dynamic::s_pm_dot_v3(*nowPmb,now2startPebDistInB,now2startPebDistInG);
}

int OpenDoor::adjustRobotPush2W(Robots::RobotBase &robot, int count)
{
    double currentPe[6] {0};
    const int now2startCount {1000};

    for (int i=0;i<6;i++)
    {
        currentPe[i]=nowPeb[i]+(now2startPebDistInG[i])/2*(1-cos((count-ODP.countIter)*PI/now2startCount));
    }

    robot.SetPeb(currentPe);
    robot.SetPee(nowPee);

    if(count-ODP.countIter+1==now2startCount)
    {
        ODP.countIter=count+1;
        return 0;
    }
    return 1;
}

int OpenDoor::adjustRobotPull2W(Robots::RobotBase &robot, int count)
{
    double currentPeb[6] {0};
    double currentPee[18] {0};
    const int now2startCount {1000};
    if(count-ODP.countIter<now2startCount)
    {
        for (int i=0;i<6;i++)
        {
            currentPeb[i]=nowPeb[i]+(now2startPebDistInG[i])/2*(1-cos((count-ODP.countIter)*PI/now2startCount));
        }
        memcpy(currentPee,nowPee,18*sizeof(double));
    }
    else if(count-ODP.countIter<2*now2startCount)
    {
        for (int i=0;i<6;i++)
        {
            currentPeb[i]=nowPeb[i]+now2startPebDistInG[i];
        }
        for(int i=0;i<3;i++)
        {
            //leg 0,2,4
            currentPee[6*i]=nowPee[6*i]+now2startPeeDistInG[0]/2*(1-cos((count-ODP.countIter-now2startCount)*PI/now2startCount));
            currentPee[6*i+1]=nowPee[6*i+1]+0.025*(1-cos((count-ODP.countIter-now2startCount)*2*PI/now2startCount));
            currentPee[6*i+2]=nowPee[6*i+2]+now2startPeeDistInG[2]/2*(1-cos((count-ODP.countIter-now2startCount)*PI/now2startCount));
            //leg 1,3,5
            currentPee[6*i+3]=nowPee[6*i+3];
            currentPee[6*i+4]=nowPee[6*i+4];
            currentPee[6*i+5]=nowPee[6*i+5];
        }
    }
    else if(count-ODP.countIter<3*now2startCount)
    {
        for (int i=0;i<6;i++)
        {
            currentPeb[i]=nowPeb[i]+now2startPebDistInG[i];
        }
        for(int i=0;i<3;i++)
        {
            //leg 1,3,5
            currentPee[6*i+3]=nowPee[6*i+3]+now2startPeeDistInG[0]/2*(1-cos((count-ODP.countIter-2*now2startCount)*PI/now2startCount));
            currentPee[6*i+4]=nowPee[6*i+4]+0.025*(1-cos((count-ODP.countIter-2*now2startCount)*2*PI/now2startCount));
            currentPee[6*i+5]=nowPee[6*i+5]+now2startPeeDistInG[2]/2*(1-cos((count-ODP.countIter-2*now2startCount)*PI/now2startCount));
            //leg 0,2,4
            currentPee[6*i]=nowPee[6*i]+now2startPeeDistInG[0];
            currentPee[6*i+1]=nowPee[6*i+1];
            currentPee[6*i+2]=nowPee[6*i+2]+now2startPeeDistInG[2];
        }
    }

    robot.SetPeb(currentPeb);
    robot.SetPee(currentPee);

    if(count-ODP.countIter==3*now2startCount-1)
    {
        ODP.countIter=count+1;
        return 0;
    }
    return 1;
}

void OpenDoor::GetPushWalkParam()
{
    const double handle2MiddleDist {0.50};
    walkParam.n=2;
    walkParam.alpha=(isLeft==true ? -PI/2 : PI/2);
    walkParam.beta=0;
    walkParam.totalCount=1000;
    walkParam.d=(isJump==false ? (handle2MiddleDist-0.03)/3*2 : (handle2MiddleDist-0.1)/3*2);
}

void OpenDoor::GetPullWalkParam()
{
    circleWalkParam.radius=1.1;
    circleWalkParam.n=4;
    circleWalkParam.beta=0.1;
    circleWalkParam.h=0.05;
    circleWalkParam.totalCount=1000;
    circleWalkParam.direction=(isLeft==false ? -1 : 1);
    circleWalkParam.startAngle=(isLeft==false ? -PI/2 : PI/2);
}

int OpenDoor::pushWalkThrough(Robots::RobotBase &robot, int count)
{
    int ret {0};
    switch(WTState)
    {
    case WalkThroughState::AlignWalk:
        walkParam.count=count-ODP.countIter;
        ret=Robots::walkGait(robot,walkParam);

        if(ret==0)
        {
            WTState=WalkThroughState::PushWalk;
            ODP.countIter=count+1;

            walkParam.totalCount=1500;
            walkParam.n=6;
            walkParam.alpha=0;
            walkParam.beta=0;
            walkParam.d=0.2;
        }
        break;

    case WalkThroughState::PushWalk:
        walkParam.count=count-ODP.countIter;
        ret=Robots::walkGait(robot,walkParam);

        if(ret==0)
        {
            return 0;
        }
        break;

    default:
        break;
    }
    return 1;
}

int OpenDoor::pullDoorWork(Robots::RobotBase &robot, int count)
{
    double currentPeb[6] {0};
    double currentPee[18] {0};
    double jumpH {0.04};
    const int legWorkCount {1000};
    int ret {0};

    switch(PWState)
    {
    case PullWorkState::CircleWalk:
        circleWalkParam.count=count-ODP.countIter;
        ret=NormalGait::circleWalk(robot,circleWalkParam);

        if(ret==0)
        {
            PWState=PullWorkState::LegWork;
            ODP.countIter=count+1;
            robot.GetPeb(nowPeb);
            robot.GetPee(nowPee,robot.body());
        }
        break;

    case PullWorkState::LegWork:
        memcpy(currentPeb,nowPeb,6*sizeof(double));
        memcpy(currentPee,nowPee,18*sizeof(double));
        if(count-ODP.countIter<legWorkCount)
        {
            if(isLeft==false)//leg 2
            {
                currentPee[9]=nowPee[9];
                currentPee[10]=nowPee[10]+0.025*(1-cos((count-ODP.countIter)*2*PI/legWorkCount));
                currentPee[11]=nowPee[11]-0.1/2*(1-cos((count-ODP.countIter)*PI/legWorkCount));
            }
            else//leg 0
            {
                currentPee[0]=nowPee[0];
                currentPee[1]=nowPee[1]+0.025*(1-cos((count-ODP.countIter)*2*PI/legWorkCount));
                currentPee[2]=nowPee[2]-0.1/2*(1-cos((count-ODP.countIter)*PI/legWorkCount));
            }
            if(count-ODP.countIter==legWorkCount-1)
            {
                memcpy(nowPee,currentPee,18*sizeof(double));
            }
        }
        else if(count-ODP.countIter<2*legWorkCount)
        {
            currentPeb[1]=nowPeb[1]+jumpH/2*(1-cos((count-ODP.countIter-legWorkCount)*PI/legWorkCount));
            if(count-ODP.countIter==2*legWorkCount-1)
            {
                memcpy(nowPeb,currentPeb,6*sizeof(double));
                aris::dynamic::s_pe2pm(nowPeb,*nowPmb,"313");
            }
        }
        else
        {
            double backDistInB[3] {0,0,0.1};
            double backDistInG[3] {0};
            aris::dynamic::s_pm_dot_v3(*nowPmb,backDistInB,backDistInG);
            currentPeb[0]=nowPeb[0]+backDistInG[0]/2*(1-cos((count-ODP.countIter-2*legWorkCount)*PI/legWorkCount));
            currentPeb[2]=nowPeb[2]+backDistInG[2]/2*(1-cos((count-ODP.countIter-2*legWorkCount)*PI/legWorkCount));
            if(count-ODP.countIter==3*legWorkCount-1)
            {
                PWState=PullWorkState::Rotate;
                ODP.countIter=count+1;

                circleWalkParam.radius=(isLeft==false ? sqrt(currentPee[9]*currentPee[9]+currentPee[11]*currentPee[11]) : sqrt(currentPee[0]*currentPee[0]+currentPee[2]*currentPee[2]));
                circleWalkParam.n=6;
                circleWalkParam.beta=PI/2/(circleWalkParam.n-0.5);
                circleWalkParam.h=0.05;
                circleWalkParam.totalCount=1000;
                circleWalkParam.startAngle=0;
                circleWalkParam.direction=(isLeft==false ? 1 : -1);
            }
        }

        robot.SetPeb(currentPeb);
        robot.SetPee(currentPee,robot.body());

        break;

    case PullWorkState::Rotate:
        circleWalkParam.count=count-ODP.countIter;
        ret=NormalGait::circleWalk(robot,circleWalkParam);

        if(ret==0)
        {
            PWState=PullWorkState::PullWalk;
            ODP.countIter=count+1;
            robot.GetPeb(nowPeb);
            robot.GetPee(nowPee,robot.body());
        }
        break;

    case PullWorkState::PullWalk:

        break;

    default:
        break;
    }

    return 1;
}

void OpenDoor::recordOpenDoorData()
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
            for (int i=0;i<6;i++)
            {
                fileGait << param.forceInB[i] << "  ";
            }
            fileGait << std::endl;
        }

        fileGait.close();
    });
}
