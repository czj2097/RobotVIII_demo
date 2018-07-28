#include "Balance.h"
#include "rtdk.h"
#include "LowPassFilter.h"

void GetReqAccInBFromFce(double *fce, double *reqAccInB)
{
    static double lstDistX {0};
    static double lstDistZ {0};

    double ballDistX=fce[5]/fce[1];
    double ballDistZ=-fce[3]/fce[1];

    double ballVx=ballDistX-lstDistX;
    double ballVz=ballDistZ-lstDistZ;

    if(ballDistX==0)
    {
        reqAccInB[0]=0;
    }
    else
    {
        reqAccInB[0]=-ballVx*ballVx/(2*ballDistX);
    }
    reqAccInB[1]=0;
    if(ballDistZ==0)
    {
        reqAccInB[2]=0;
    }
    else
    {
        reqAccInB[2]=-ballVz*ballVz/(2*ballDistZ);
    }

    lstDistX=ballDistX;
    lstDistZ=ballDistZ;
}

void GetTargetEulFromAcc(double *planAccInG, double *reqAccInG, double *targetEul)
{
    const double g=9.80665;
    double FnInG[3];//Fn/m, direction is the same, m can be not taken into account
    double ballG[3] {0,-g,0};
    for(int i=0;i<3;i++)
    {
        FnInG[i]=reqAccInG[i]+planAccInG[i]-ballG[i];
    }
    targetEul[0]=0;
    targetEul[1]=atan2(FnInG[2],FnInG[1]);
    targetEul[2]=-asin(FnInG[0]/GeneralFunc::norm(FnInG));

//    double tiltAngle=atan((g-acc[1])/sqrt(acc[0]*acc[0]+acc[2]*acc[2]));
//    double tiltAxis[3];

//    double yAxis[3] {0,1,0};
//    double acc_tmp[3] {0};
//    acc_tmp[0]=acc[0];
//    acc_tmp[2]=acc[2];
//    aris::dynamic::s_cro3(yAxis,acc_tmp,tiltAxis);

//    double Pmb[4][4];
//    GeneralFunc::GetPmOfU(tiltAngle,tiltAxis,*Pmb);
//    aris::dynamic::s_pm2pe(*Pmb,eul213,"213");
}

void parseBalance(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
{
    BalanceParam param;

    for(auto &i:params)
    {
        if(i.first=="distance")
        {
            param.d=stod(i.second);
        }
        else if(i.first=="n")
        {
            param.n=stoi(i.second);
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
    param.preCount=120;

    msg.copyStruct(param);

    std::cout<<"finished parse"<<std::endl;
}

int balance(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
{
    auto &robot = static_cast<Robots::RobotBase &>(model);
    auto &param = static_cast<const BalanceParam &>(param_in);

    static double beginPeb213[6];
    static double beginPee[18];

    if(param.count==0 || (param.count-param.preCount)%param.totalCount==0)
    {
        robot.GetPeb(beginPeb213,"213");
        robot.GetPee(beginPee);
    }

    //fce init
    double fceInF[6];
    double fceInB[6];
    double fceInB_filtered[6];
    static double fceMtx[10][6] {0};
    ForceTask::forceInit(param.count,param.ruicong_data->at(0).force[0].fce,fceInF);
    aris::dynamic::s_f2f(*robot.forceSensorMak().prtPm(),fceInF,fceInB);

    for(int i=0;i<9;i++)
    {
        memcpy(*fceMtx+6*i,*fceMtx+6*i+6,6*sizeof(double));
    }
    memcpy(*fceMtx+6*9,fceInB,6*sizeof(double));
    double fceSum[6] {0};
    for(int i=0;i<6;i++)
    {
        for(int j=0;j<10;j++)
        {
            fceSum[i]+=fceMtx[j][i];
        }
        fceInB_filtered[i]=fceSum[i]/10;
    }

    if(param.count<param.preCount)
    {
        robot.SetPeb(beginPeb213,"213");
        robot.SetPee(beginPee);
        return 1;
    }

    double realPeb213[6];
    memcpy(realPeb213,beginPeb213,6*sizeof(double));
    double t=(double)((param.count-param.preCount)%param.totalCount)/param.totalCount;
    double c[6];
    if((param.count-param.preCount)/param.totalCount==0)
    {
        GeneralFunc::Get5thPolynomial(beginPeb213[2],0,0,beginPeb213[2]-param.d,0,0,0.001*param.totalCount,c);
    }
    else if((param.count-param.preCount)/param.totalCount==2*param.n-1)
    {
        GeneralFunc::Get5thPolynomial(beginPeb213[2],0,0,beginPeb213[2]+param.d,0,0,0.001*param.totalCount,c);
    }
    else if(((param.count-param.preCount)/param.totalCount)%2==0)
    {
        GeneralFunc::Get5thPolynomial(beginPeb213[2],0,0,beginPeb213[2]-2*param.d,0,0,0.001*param.totalCount,c);
    }
    else
    {
        GeneralFunc::Get5thPolynomial(beginPeb213[2],0,0,beginPeb213[2]+2*param.d,0,0,0.001*param.totalCount,c);
    }
    realPeb213[2]=c[0]*(1-t)*(1-t)*(1-t)*(1-t)*(1-t)
               +5*c[1]*(1-t)*(1-t)*(1-t)*(1-t)*t
              +10*c[2]*(1-t)*(1-t)*(1-t)*t*t
              +10*c[3]*(1-t)*(1-t)*t*t*t
               +5*c[4]*(1-t)*t*t*t*t
                 +c[5]*t*t*t*t*t;

    double fce[6] {0,1,0,0,0,0};
    double reqAccInB[3] {0};
    double reqAccInG[3] {0};
    double lstPmb[4][4];
    robot.GetPmb(*lstPmb);
    GetReqAccInBFromFce(fce,reqAccInB);
    aris::dynamic::s_pm_dot_v3(*lstPmb,reqAccInB,reqAccInG);
    //rt_printf("reqAccInB:%f,%f,%f\n",reqAccInB[0],reqAccInB[1],reqAccInB[2]);
    //rt_printf("reqAccInG:%f,%f,%f\n",reqAccInG[0],reqAccInG[1],reqAccInG[2]);

    double planAccInG[3] {0};
    planAccInG[2]=(20*c[0]-40*c[1]+20*c[2])*(1-t)*(1-t)*(1-t)
                 +(60*c[1]-120*c[2]+60*c[3])*(1-t)*(1-t)*t
                 +(60*c[2]-120*c[3]+60*c[4])*(1-t)*t*t
                 +(20*c[3]-40*c[4]+20*c[5])*t*t*t;
    rt_printf("planAcc:%f,%f,%f\n",planAccInG[0],planAccInG[1],planAccInG[2]);
    GetTargetEulFromAcc(planAccInG,reqAccInG,realPeb213+3);

    rt_printf("peb:%f,%f,%f,%f,%f,%f\n",realPeb213[0],realPeb213[1],realPeb213[2],realPeb213[3],realPeb213[4],realPeb213[5]);
    //rt_printf("pee:%f,%f,%f,%f,%f,%f\n",beginPee[0],beginPee[1],beginPee[2],beginPee[3],beginPee[4],beginPee[5]);
    robot.SetPeb(realPeb213,"213");
    robot.SetPee(beginPee);

    return 2*param.n*param.totalCount+param.preCount-param.count-1;
}

/*
<bl default="bl_param">
    <bl_param type="group">
        <distance abbreviation="d" type="double" default="0.1"/>
        <totalCount abbreviation="t" type="int" default="1000"/>
        <n abbreviation="n" type="int" default="2"/>
    </bl_param>
</bl>

rs.addCmd("bl",parseBalance,balance);
 */


BalanceState BallBalance::workState;
int BallBalance::countIter;
Pipe<BallBalanceParam> BallBalance::bbPipe(true);
std::thread BallBalance::bbThread;


double BallBalance::GetAngleFromAcc(double acc)
{
    const double g=9.80665;
    const double angleLmt=10;//degree
    double accLmt=g*sin(angleLmt/180*PI);

    if(acc>accLmt)
    {
        acc=accLmt;
    }
    else if(acc<-accLmt)
    {
        acc=-accLmt;
    }

    return asin(acc/g);
}

void BallBalance::parseBallBalance(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
{
    aris::server::GaitParamBase param;

    for(auto &i:params)
    {
        if(i.first=="init")
        {
            workState=BalanceState::Init;
            msg.copyStruct(param);
        }
        else if(i.first=="balance")
        {
            workState=BalanceState::Balance;
        }
        else if(i.first=="stop")
        {
            workState=BalanceState::Stop;
        }
        else
        {
            std::cout<<"parse failed"<<std::endl;
        }
    }

    std::cout<<"finished parse"<<std::endl;
}

int BallBalance::ballBalance(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
{
    auto &robot = static_cast<Robots::RobotBase &>(model);
    auto &param = static_cast<const aris::server::GaitParamBase &>(param_in);

    BallBalanceParam bbParam;

    EmergencyStop stoper;
    stoper.doRecord(robot);

    /*****IMU*****/
    double peImuGrnd2BodyGrnd[6] {0,0,0,PI/6,-PI/2,0};//213
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

    const double m {400};
    double mInG[3] {0,-m,0};
    double mInB[3];
    aris::dynamic::s_inv_pm_dot_v3(bodyPm,mInG,mInB);

    /*****FSR*****/
    double fceInF[6];
    double fceInB[6];
    ForceTask::forceInit(param.count,param.ruicong_data->at(0).force[6].fce,fceInF);
    aris::dynamic::s_f2f(*robot.forceSensorMak().prtPm(),fceInF,fceInB);
    //rt_printf("fceInF:%f,%f,%f,fceInB:%f,%f,%f\n",fceInF[0],fceInF[1],fceInF[2],fceInB[0],fceInB[1],fceInB[2]);

    double fceInB_filtered[6];
    static LowpassFilter<6> filter;
    if(param.count==0)
    {
        filter.Initialize();
        filter.SetCutFrequency(0.010,1000);
    }
    filter.DoFilter(fceInB,fceInB_filtered);

//    double ballPosZ=-fceInB_filtered[3]/fceInB_filtered[1];
//    double ballPosX=fceInB_filtered[5]/fceInB_filtered[1];

    double ballPosZ=-fceInB_filtered[3]/mInB[1];
    double ballPosX=fceInB_filtered[5]/mInB[1];

    static double beginPeb213[6];
    static double beginPee[18];
    static Controller::lstPIDparam zPIDparam;
    static Controller::lstPIDparam xPIDparam;
    static Controller::lstLagParam zLagParam;
    static Controller::lstLagParam xLagParam;

    double tarAccZ {0};
    double tarAccX {0};
    double lagStartEul[3] {0};
    double tarEul[3] {0};//213 eul
    double actEul[3] {0};//213 eul
    double fx { 0.8 };
    double fz { 0.8 };

    switch(workState)
    {
    case BalanceState::Init:
        robot.GetPeb(beginPeb213,"213");
        robot.GetPee(beginPee);

        zPIDparam.lstErr=ballPosZ;
        zPIDparam.lstInt=0;
        xPIDparam.lstErr=ballPosX;
        xPIDparam.lstInt=0;

        zLagParam.lstFstInt=0;
        zLagParam.lstSndInt=0;
        xLagParam.lstFstInt=0;
        xLagParam.lstSndInt=0;

        countIter=param.count+1;
	
        memcpy(bbParam.fceInB,fceInB,6*sizeof(double));
        memcpy(bbParam.fceInB_filtered,fceInB_filtered,6*sizeof(double));
        memcpy(bbParam.bodyPos,bodyPe,3*sizeof(double));
        bbParam.ballPos[0]=ballPosZ;
        bbParam.ballPos[1]=ballPosX;
        bbParam.tarAcc[0]=tarAccZ;
        bbParam.tarAcc[1]=tarAccX;
        bbParam.tarEul[0]=tarEul[1];
        bbParam.tarEul[1]=tarEul[2];
        bbParam.actEul[0]=tarEul[1];
        bbParam.actEul[1]=tarEul[2];
        bbParam.imuEul[0]=bodyPe[4];
        bbParam.imuEul[1]=bodyPe[5];
        bbPipe.sendToNrt(bbParam);

        return 1;
        break;

    case BalanceState::Balance:
        if((param.count-countIter)%1000==0)
        {
            memcpy(lagStartEul,bodyPe+3,3*sizeof(double));
        }

        tarAccZ=Controller::doPID(zPIDparam,ballPosZ,3,0,3,0.001);
        tarAccX=Controller::doPID(xPIDparam,ballPosX,3,0,3,0.001);

        tarEul[1]=-GetAngleFromAcc(tarAccZ);
        tarEul[2]=GetAngleFromAcc(tarAccX);

        actEul[1]=Controller::SndOrderLag(zLagParam,0,tarEul[1]-lagStartEul[1],2*PI*fz,1,0.001);
        actEul[2]=Controller::SndOrderLag(xLagParam,0,tarEul[2]-lagStartEul[2],2*PI*fx,1,0.001);

        rt_printf("pos:%f, acc:%f, tarEul:%f, actEul:%f\n",ballPosZ,tarAccZ,tarEul[1],actEul[1]);

        beginPeb213[3+1]=actEul[1];
        beginPeb213[3+2]=actEul[2];

        robot.SetPeb(beginPeb213,"213");
        robot.SetPee(beginPee);

        memcpy(bbParam.fceInB,fceInB,6*sizeof(double));
        memcpy(bbParam.fceInB_filtered,fceInB_filtered,6*sizeof(double));
        memcpy(bbParam.bodyPos,bodyPe,3*sizeof(double));
        bbParam.ballPos[0]=ballPosZ;
        bbParam.ballPos[1]=ballPosX;
        bbParam.tarAcc[0]=tarAccZ;
        bbParam.tarAcc[1]=tarAccX;
        bbParam.tarEul[0]=tarEul[1];
        bbParam.tarEul[1]=tarEul[2];
        bbParam.actEul[0]=tarEul[1];
        bbParam.actEul[1]=tarEul[2];
        bbParam.imuEul[0]=bodyPe[4];
        bbParam.imuEul[1]=bodyPe[5];
        bbPipe.sendToNrt(bbParam);

        return 1;
        break;

    case BalanceState::Stop:
        return stoper.stop(robot,param.count,3.2);
        break;

    default:
        break;
    }

    return 0;
}

void BallBalance::recordData()
{
    bbThread = std::thread([&]()
    {
        struct BallBalanceParam param;
        static std::fstream fileGait;
        std::string name = aris::core::logFileName();
        name.replace(name.rfind("log.txt"), std::strlen("ballBalanceData.txt"), "ballBalanceData.txt");
        fileGait.open(name.c_str(), std::ios::out | std::ios::trunc);

        long long count = -1;
        while (1)
        {
            bbPipe.recvInNrt(param);

            fileGait << ++count << " ";
            for (int i=0;i<6;i++)
            {
                fileGait << param.fceInB[i] << "  ";
            }
            for (int i=0;i<6;i++)
            {
                fileGait << param.fceInB_filtered[i] << "  ";
            }
            for (int i=0;i<3;i++)
            {
                fileGait << param.bodyPos[i] << "  ";
            }
            for(int i=0;i<2;i++)
            {
                fileGait << param.ballPos[i] << "  ";
                fileGait << param.tarAcc[i] << "  ";
                fileGait << param.tarEul[i] << "  ";
                fileGait << param.actEul[i] << "  ";
                fileGait << param.imuEul[i] << "  ";
            }
            fileGait << std::endl;
        }

        fileGait.close();
    });
}

/*
<bb default="bb_param">
    <bb_param type="unique">
        <init abbreviation="i"/>
        <balance abbreviation="b"/>
        <stop abbreviation="s"/>
    </bb_param>
</bb>

rs.addCmd("bb",BallBalance::parseBallBalance,BallBalance::ballBalance);
BallBalance::recordData();
 */
