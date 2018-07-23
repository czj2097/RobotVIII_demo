#include "Balance.h"
#include "rtdk.h"

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
    BallBalanceParam param;

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
    auto &param = static_cast<const BallBalanceParam &>(param_in);

    EmergencyStop stoper;
    stoper.doRecord(robot);

    const double m {0.4};
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

    double ballPosZ=-fceInB_filtered[3]/fceInB_filtered[1];
    double ballPosX=fceInB_filtered[5]/fceInB_filtered[1];

    static double beginPeb213[6];
    static double beginPee[18];
    static Controller::lstPIDparam zPIDparam;
    static Controller::lstPIDparam xPIDparam;
    static Controller::lstLagParam zLagParam;
    static Controller::lstLagParam xLagParam;

    double tarAcc[3] {0};
    double tarEul[3] {0};//213 eul
    double actEul[3] {0};//213 eul
    double f {2};

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
        zLagParam.lstSndInt=0;

        return 1;
        break;

    case BalanceState::Balance:
        tarAcc[2]=Controller::doPID(zPIDparam,ballPosZ,3,0.5,3,0.001);
        tarAcc[0]=Controller::doPID(xPIDparam,ballPosX,3,0.5,3,0.001);

        tarEul[1]=-GetAngleFromAcc(tarAcc[2]);
        tarEul[2]=GetAngleFromAcc(tarAcc[0]);

        actEul[1]=Controller::SndOrderLag(zLagParam,0,tarEul[1],2*PI*f,1,0.001);
        actEul[2]=Controller::SndOrderLag(xLagParam,0,tarEul[2],2*PI*f,1,0.001);

        beginPeb213[3+1]=actEul[1];
        beginPeb213[3+2]=actEul[2];

        robot.SetPeb(beginPeb213,"213");
        robot.SetPee(beginPee);
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

/*
<bb default="bb_param">
    <bb_param type="unique">
        <init abbreviation="i"/>
        <balance abbreviation="b"/>
        <stop abbreviation="s"/>
    </bb_param>
</bb>

rs.addCmd("bb",BallBalance::parseBallBalance,BallBalance::ballBalance);
 */
