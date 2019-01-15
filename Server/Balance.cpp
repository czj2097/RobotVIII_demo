#include "Balance.h"
#include "rtdk.h"
#include "LowPassFilter.h"
#include "LeadCompensator.h"
#include "Calman2d.h"
#include "Calman1d.h"

void GetTarEulFromAcc(double *planAccInG, double *reqAccInG, double *targetEul)
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

void GetTarEulFromAcc(double *reqAccInB, double *tarEul)
{
    const double g=9.80665;
    double FnInG[3];//Fn/m, direction is the same, m can be not taken into account
    double ballG[3] {0,-g,0};
    for(int i=0;i<3;i++)
    {
        FnInG[i]=reqAccInB[i]-ballG[i];
    }
    tarEul[0]=0;
    tarEul[1]=atan2(FnInG[2],FnInG[1]);
    tarEul[2]=-asin(FnInG[0]/GeneralFunc::norm(FnInG));
    for(int i=0;i<3;i++)
    {
        if(tarEul[i]>(10*PI/180))
        {
            tarEul[i]=10*PI/180;
        }
        if(tarEul[i]<(-10*PI/180))
        {
            tarEul[i]=-10*PI/180;
        }
    }
//    rt_printf("%f,%f,%f\n",tarEul[0],tarEul[1],tarEul[2]);
}

void GetAccFromActEul(double *actAccInB, double *actEul)
{
    const double g=9.80665;
    double ballG[3] {0,-g,0};
    double actPm[16];
    aris::dynamic::s_pe2pm(actEul,actPm,"213");
    aris::dynamic::s_inv_pm_dot_v3(actPm,ballG,actAccInB);
}

void GetAlphaFromEul(double *eul, double *w, double *alpha)
{
    static double lst_u[3] {0};
    static double lst_w[3] {0};
    double cur_u[3] {eul[2],eul[1],eul[0]};//312
    double cur_w[3] {0};
    double d_u[3] {0};

    double mtrx[9] {0, cos(cur_u[0]),sin(cur_u[0])*cos(cur_u[1]),
                    0,-sin(cur_u[0]),cos(cur_u[0])*cos(cur_u[1]),
                    1, 0            ,-sin(cur_u[1]) };
    for(int i=0;i<3;i++)//312
    {
        d_u[i]=(cur_u[i]-lst_u[i])/0.001;
    }
    aris::dynamic::s_dgemm(3,1,3,1,mtrx,3,d_u,1,0,cur_w,1);

    for(int i=0;i<3;i++)//123
    {
        w[i]=cur_w[i];
        alpha[i]=(cur_w[i]-lst_w[i])/0.001;
    }

    memcpy(lst_u,cur_u,3*sizeof(double));
    memcpy(lst_w,cur_w,3*sizeof(double));
}

double BallBalance::realPeb213[6];
BalanceWalkState BallBalance::workState;
int BallBalance::countIter;
BallBalanceParam BallBalance::bbParam;
Pipe<BallBalanceParam> BallBalance::bbPipe(true);
std::thread BallBalance::bbThread;


void BallBalance::bodyPosTg(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
{
    auto &robot = static_cast<Robots::RobotBase &>(model);
    auto &param = static_cast<const BalanceParam &>(param_in);
	
	double beginPeb213[6];
	if((param.count-countIter)%param.totalCount==0)
    {
        robot.GetPeb(beginPeb213,"213");
    }

    double c[6];
    double t=(double)((param.count-countIter)%param.totalCount)/param.totalCount;
    if(param.count-countIter<2*param.n*param.totalCount)
    {
        if((param.count-countIter)/param.totalCount==0)
        {
            GeneralFunc::Get5thPolynomial(beginPeb213[2],0,0,beginPeb213[2]-param.d,0,0,0.001*param.totalCount,c);
        }
        else if((param.count-countIter)/param.totalCount==2*param.n-1)
        {
            GeneralFunc::Get5thPolynomial(beginPeb213[2],0,0,beginPeb213[2]+param.d,0,0,0.001*param.totalCount,c);
        }
        else if(((param.count-countIter)/param.totalCount)%2==0)
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
    }
}

void BallBalance::bodyEulTg(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
{
    auto &robot = static_cast<Robots::RobotBase &>(model);
    auto &param = static_cast<const BalanceParam &>(param_in);

    static Controller::lstPIDparam zPIDparam;
    static Controller::lstPIDparam xPIDparam;
    static Controller::lstLagParam zLagParam;
    static Controller::lstLagParam xLagParam;

    double lagStartEul[2] {0};
    double fx { 0.8 };
    double fz { 0.8 };

    if(param.count==countIter)
    {
        zPIDparam.lstErr=bbParam.ballPos[0];//pos on z-axis
        zPIDparam.lstInt=0;
        xPIDparam.lstErr=bbParam.ballPos[1];//pos on x-axis
        xPIDparam.lstInt=0;

        zLagParam.lstFstInt=0;
        zLagParam.lstSndInt=0;
        xLagParam.lstFstInt=0;
        xLagParam.lstSndInt=0;

        std::fill_n(bbParam.ballPos_filtered,2,0);
    }
//    if((param.count-countIter)%1000==0)
    {
        lagStartEul[0]=bbParam.actEul[1];
        lagStartEul[1]=bbParam.actEul[2];
    }

    static double lst_calmanPos[2] {0};
    bbParam.ballVel[0]=(bbParam.ballPos_filtered[0]-lst_calmanPos[0])/0.001;
    bbParam.ballVel[1]=(bbParam.ballPos_filtered[1]-lst_calmanPos[1])/0.001;
    memcpy(lst_calmanPos,bbParam.ballPos_filtered,2*sizeof(double));
    doCalman2d();
    bbParam.tarAcc[0]=1.5*bbParam.ballPos_filtered[0]+1.5*bbParam.ballVel_filtered[0];
    bbParam.tarAcc[1]=1.5*bbParam.ballPos_filtered[1]+1.5*bbParam.ballVel_filtered[1];

//    bbParam.tarAcc[0]=Controller::doPID(zPIDparam,bbParam.ballPos[0],1.5,0,1.5,0.001);//acc on z-axis
//    bbParam.tarAcc[1]=Controller::doPID(xPIDparam,bbParam.ballPos[1],1.5,0,1.5,0.001);//acc on x-axis

//    double ratio[2];//change from 1/3 to 1 as ballPos change from 0 to 0.2m;
//    for(int i=0;i<2;i++)
//    {
//        if(fabs(bbParam.ballPos[i])>0.01)
//        {
//            ratio[i]=1;
//        }
//        else
//        {
//            ratio[i]=1.0;
//        }
//        bbParam.tarAcc[i]*=ratio[i];
//    }

//    bbParam.tarEul[0]=-GetAngleFromAcc(bbParam.tarAcc[0]);//angle around x-axis
//    bbParam.tarEul[1]=GetAngleFromAcc(bbParam.tarAcc[1]);//angle around z-axis
    double tarAcc_tmp[3] {bbParam.tarAcc[1],0,bbParam.tarAcc[0]};
    double tarEul_tmp[3];
    GetTarEulFromAcc(tarAcc_tmp,tarEul_tmp);
    bbParam.tarEul[0]=-tarEul_tmp[1]-0.013;
    bbParam.tarEul[1]=-tarEul_tmp[2]+0.006;

    bbParam.actEul[1]=Controller::SndOrderLag(zLagParam,0,bbParam.tarEul[0]-lagStartEul[0],2*PI*fz,1,0.001);
    bbParam.actEul[2]=Controller::SndOrderLag(xLagParam,0,bbParam.tarEul[1]-lagStartEul[1],2*PI*fx,1,0.001);

    realPeb213[3+1]=bbParam.actEul[1];
    realPeb213[3+2]=bbParam.actEul[2];
}

void BallBalance::doCalman1d()
{
    double acc[3];
    GetAccFromActEul(acc,bbParam.imuEul);

    const double delta_t {0.001};
    static double lst_x1 {0};
    static double lst_P1 {0};
    static double lst_x2 {0};
    static double lst_P2 {0};
    static double pre_x1 {0};
    static double pre_x2 {0};

    double mtrx_F {1};
    double input_u1 {delta_t*acc[0]};//vel on x
    double input_u2 {delta_t*acc[2]};//vel on z
    double noise_Q {9e-4*9e-4};
    double mtrx_H {1};
    double noise_R {0.4*0.4};
    double input_z1 {bbParam.ballVel[1]};//vel on x
    double input_z2 {bbParam.ballVel[0]};//vel on z

    Calman1d filter1;
    filter1.init(mtrx_F,noise_Q,mtrx_H,noise_R,input_u1,input_z1);
    filter1.doFilter(lst_x1,pre_x1,lst_P1);

    Calman1d filter2;
    filter2.init(mtrx_F,noise_Q,mtrx_H,noise_R,input_u2,input_z2);
    filter2.doFilter(lst_x2,pre_x2,lst_P2);

    bbParam.predictVel[1]=pre_x1;
    bbParam.predictVel[0]=pre_x2;
    bbParam.ballVel_filtered[1]=lst_x1;//pos on x
    bbParam.ballVel_filtered[0]=lst_x2;//pos on z
}

void BallBalance::doCalman2d()
{
    double acc[3];
    GetAccFromActEul(acc,bbParam.imuEul);

    const double delta_t {0.001};
    static double lst_x1[2] {0};
    static double lst_P1[2][2] {0};
    static double lst_x2[2] {0};
    static double lst_P2[2][2] {0};
    static double pre_x1[2] {0};
    static double pre_x2[2] {0};

    double mtrx_F[2][2] {1, delta_t, 0, 1};
    double input_u1[2] {delta_t*delta_t/2*acc[0], delta_t*acc[0]};//pos on x
    double input_u2[2] {delta_t*delta_t/2*acc[2], delta_t*acc[2]};//pos on z
    double noise_Q[2][2] {5e-4*5e-4, 0, 0, 9e-3*9e-3};
    double mtrx_H[2][2] {1, 0, 0, 1};
    double noise_R[2][2] {0.08*0.08, 0, 0, 0.4*0.4};
    double input_z1[2] {bbParam.ballPos[1], bbParam.ballVel[1]};
    double input_z2[2] {bbParam.ballPos[0], bbParam.ballVel[0]};

    Calman2d filter1;
    filter1.init(*mtrx_F,*noise_Q,*mtrx_H,*noise_R,input_u1,input_z1);
    filter1.doFilter(lst_x1,pre_x1,*lst_P1);

    Calman2d filter2;
    filter2.init(*mtrx_F,*noise_Q,*mtrx_H,*noise_R,input_u2,input_z2);
    filter2.doFilter(lst_x2,pre_x2,*lst_P2);


    bbParam.ballPos_filtered[1]=lst_x1[0];//pos on x
    bbParam.ballVel_filtered[1]=lst_x1[1];
    bbParam.ballPos_filtered[0]=lst_x2[0];//pos on z
    bbParam.ballVel_filtered[0]=lst_x2[1];

    bbParam.predictPos[1]=pre_x1[0];
    bbParam.predictVel[1]=pre_x1[1];
    bbParam.predictPos[0]=pre_x2[0];
    bbParam.predictVel[0]=pre_x2[1];
}

void BallBalance::parseBallBalance(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
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
        else if(i.first=="init")
        {
            workState=BalanceWalkState::Init;
        }
        else if(i.first=="balance")
        {
            workState=BalanceWalkState::Balance;
        }
        else if(i.first=="stop")
        {
            workState=BalanceWalkState::Stop;
        }
        else
        {
            std::cout<<"parse failed"<<std::endl;
        }
    }
    if(workState==BalanceWalkState::Init)
    {
	printf("totalCount=%d,n=%d,d=%f\n",param.totalCount,param.n,param.d);
        msg.copyStruct(param);
    }
    std::cout<<"finished parse"<<std::endl;
}

int BallBalance::ballBalance(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
{
    auto &robot = static_cast<Robots::RobotBase &>(model);
    auto &param = static_cast<const BalanceParam &>(param_in);

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
    memcpy(bbParam.imuEul,bodyPe+3,3*sizeof(double));
    for(int i=0;i<3;i++)
    {
        if(bbParam.imuEul[i]>PI)
        {
            bbParam.imuEul[i]-=2*PI;
        }
    }

    const double m {400};
    double mInG[3] {0,-m,0};
    double mInB[3];
    aris::dynamic::s_inv_pm_dot_v3(bodyPm,mInG,mInB);

    /*****FSR*****/
    double fceInF[6];
    GeneralFunc::forceInit(param.count,param.ruicong_data->at(0).force[6].fce,fceInF);
    aris::dynamic::s_f2f(*robot.forceSensorMak().prtPm(),fceInF,bbParam.fceInB);
    //rt_printf("fceInF:%f,%f,%f,fceInB:%f,%f,%f\n",fceInF[0],fceInF[1],fceInF[2],bbParam.fceInB[0],bbParam.fceInB[1],bbParam.fceInB[2]);

    //inertial force of the table
//    double J=0.045238;
//    GetAlphaFromEul(bbParam.actEul,bbParam.w,bbParam.alpha);
//    memcpy(bbParam.fceInB_inertia,bbParam.fceInB,6*sizeof(double));
//    bbParam.fceInB_inertia[3]-=J*bbParam.alpha[0];
//    bbParam.fceInB_inertia[5]-=J*bbParam.alpha[2];

    static LowpassFilter<6> filter;
    if(param.count==0)
    {
        filter.Initialize();
        filter.SetCutFrequency(0.020,1000);
    }
    filter.DoFilter(bbParam.fceInB,bbParam.fceInB_filtered);

    static LeadCompensator<6> compensator;
    if(param.count==0)
    {
        compensator.Initialize();
        compensator.SetLeadParam(1.6,1.6,0.3,1000);
    }
    compensator.Compensate(bbParam.fceInB_filtered,bbParam.fceInB_lead);

//    bbParam.ballPos[0]=-bbParam.fceInB_filtered[3]/bbParam.fceInB_filtered[1];
//    bbParam.ballPos[1]=bbParam.fceInB_filtered[5]/bbParam.fceInB_filtered[1];

    bbParam.ballPos[0]=-bbParam.fceInB_lead[3]/mInB[1];
    bbParam.ballPos[1]=bbParam.fceInB_lead[5]/mInB[1];

    static double beginPeb213[6];
    static double beginPee[18];

    switch(workState)
    {
    case BalanceWalkState::Init:
        robot.GetPeb(beginPeb213,"213");
        robot.GetPee(beginPee);

        countIter=param.count+1;
	
        memcpy(bbParam.bodyPos,beginPeb213,3*sizeof(double));
        std::fill_n(bbParam.tarAcc,2,0);
        std::fill_n(bbParam.tarEul,2,0);
        std::fill_n(bbParam.actEul,3,0);
        bbPipe.sendToNrt(bbParam);

        return 1;
        break;

    case BalanceWalkState::Balance:
        if((param.count-countIter)%param.totalCount==0)
        {
            robot.GetPeb(beginPeb213,"213");
            robot.GetPee(beginPee);
        }
        memcpy(realPeb213,beginPeb213,6*sizeof(double));

        /*****pos*****/
        //bodyPosTg(robot,param);

        /*****eul*****/
        bodyEulTg(robot,param);

        robot.SetPeb(realPeb213,"213");
        robot.SetPee(beginPee);

        memcpy(bbParam.bodyPos,realPeb213,3*sizeof(double));
        bbPipe.sendToNrt(bbParam);

        return 1;
        break;

    case BalanceWalkState::Stop:
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
//            for (int i=0;i<6;i++)
//            {
//                fileGait << param.fceInB_inertia[i] << "  ";
//            }
            for (int i=0;i<6;i++)
            {
                fileGait << param.fceInB_filtered[i] << "  ";
            }
            for (int i=0;i<6;i++)
            {
                fileGait << param.fceInB_lead[i] << "  ";
            }
            for (int i=0;i<3;i++)
            {
                fileGait << param.bodyPos[i] << "  ";
            }
            for(int i=0;i<2;i++)
            {
                fileGait << param.ballPos[i] << "  ";
                fileGait << param.predictPos[i] << "  ";
                fileGait << param.ballPos_filtered[i] << "  ";
                fileGait << param.ballVel[i] << "  ";
                fileGait << param.predictVel[i] << "  ";
                fileGait << param.ballVel_filtered[i] << "  ";
            }
            for(int i=0;i<2;i++)
            {
                fileGait << param.tarAcc[i] << "  ";
                fileGait << param.tarEul[i] << "  ";
                fileGait << param.actEul[i+1] << "  ";
                fileGait << param.imuEul[i+1] << "  ";
            }
            fileGait << std::endl;
        }

        fileGait.close();
    });
}

/*
{ 0,0,0,PI,PI/2,PI*3/2 }

<bb default="bb_param">
        <bb_param type="group">
            <state_param type="unique"/>
                <init abbreviation="i"/>
                <balance abbreviation="b"/>
                <stop abbreviation="s"/>
            </state_param>
            <distance abbreviation="d" type="double" default="0.1"/>
            <totalCount abbreviation="t" type="int" default="1000"/>
            <n abbreviation="n" type="int" default="2"/>
    </bb_param>
</bb>

rs.addCmd("bb",BallBalance::parseBallBalance,BallBalance::ballBalance);
BallBalance::recordData();
 */


double BalanceWalk::forwardAcc;
double BalanceWalk::turnAcc;
int BalanceWalk::totalCount;
int BalanceWalk::totalCount_tmp;
double BalanceWalk::height;
double BalanceWalk::height_tmp;
double BalanceWalk::alpha;
double BalanceWalk::alpha_tmp;

BalanceWalkState BalanceWalk::walkState;
bool BalanceWalk::constFlag;
double BalanceWalk::beginXVel;
double BalanceWalk::endXVel;
double BalanceWalk::endXVelinB;
double BalanceWalk::beginZVel;
double BalanceWalk::endZVel;
double BalanceWalk::endZVelinB;
double BalanceWalk::beginOmega;
double BalanceWalk::endOmega;
double BalanceWalk::endGama;

double BalanceWalk::beginPeb[6];
double BalanceWalk::pEB[6];
GaitPhase BalanceWalk::gaitPhase[6];
double BalanceWalk::swingPee[18];
double BalanceWalk::swingBeginPee[18];
double BalanceWalk::swingEndPee[18];
double BalanceWalk::stancePee[18];
double BalanceWalk::stanceBeginPee[18];
double BalanceWalk::stanceEndPee[18];
double BalanceWalk::followPee[18];
double BalanceWalk::followRecordPee[18];
double BalanceWalk::followBeginPee[18];
double BalanceWalk::followBeginVee[18];
int BalanceWalk::followBeginCount[6];
bool BalanceWalk::isTrans2Follow[6];
bool BalanceWalk::isPreTrans2Follow[6];
bool BalanceWalk::isFollowEnd[6];
bool BalanceWalk::isSndHalf[6];
int BalanceWalk::preTrans2FollowCount[6];

double BalanceWalk::initPee[18];
double BalanceWalk::avgRealH;
double BalanceWalk::planH;

double BalanceWalk::bodyEul213[3];
double BalanceWalk::sumEul[3];
double BalanceWalk::avgEul[3];
double BalanceWalk::targetEul[3];
double BalanceWalk::inputEul[3];
double BalanceWalk::inputEul_tmp[3];

double BalanceWalk::forceRaw[42];
double BalanceWalk::forceInF[42];
double BalanceWalk::forceSum[42];
double BalanceWalk::forceAvg[42];

double BalanceWalk::imuRaw[3];
double BalanceWalk::imuInB[3];
double BalanceWalk::imuSum[3];
double BalanceWalk::imuAvg[3];

BalanceWalkParam BalanceWalk::bwParam;
Pipe<BalanceWalkParam> BalanceWalk::bwPipe(true);
std::thread BalanceWalk::bwThread;


void BalanceWalk::parseBalanceWalk(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
{
    aris::server::GaitParamBase param;

    for(auto &i:params)
    {
        if(i.first=="init")
        {
            walkState=BalanceWalkState::Init;
            msg.copyStruct(param);
        }
        else if(i.first=="forwardAcc")
        {
            forwardAcc=stod(i.second);
            walkState=BalanceWalkState::ForwardAcc;
        }
        else if(i.first=="betaAcc")
        {
            turnAcc=stod(i.second);
            walkState=BalanceWalkState::TurnAcc;
        }
        else if(i.first=="stop")
        {
            walkState=BalanceWalkState::Stop;
        }
        else if(i.first=="totalCount")
        {
            totalCount=stoi(i.second);
        }
        else if(i.first=="height")
        {
            height_tmp=stod(i.second);
        }
        else if(i.first=="roll")
        {
            inputEul_tmp[2]=stod(i.second);
        }
        else if(i.first=="pitch")
        {
            inputEul_tmp[1]=stod(i.second);
        }
        else if(i.first=="alpha")
        {
            alpha_tmp=stod(i.second);
        }
        else
        {
            std::cout<<"parse failed"<<std::endl;
        }
    }

    std::cout<<"finished parse"<<std::endl;
}

void BalanceWalk::forceInit(int count,int legID)
{
    if(count==0)
    {
        std::fill(forceSum+6*legID,forceSum+6*legID+6,0);
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

void BalanceWalk::imuInit(int count)
{
    double peImuGrnd2BodyGrnd[6] {0,0,0,PI/6,-PI/2,0};//213
    double pmImuGrnd2BodyGrnd[16];
    double pmBody2Imu[16];
    double peImu2ImuGrnd[6] {0};
    peImu2ImuGrnd[3]=0;
    peImu2ImuGrnd[4]=imuRaw[1];
    peImu2ImuGrnd[5]=imuRaw[0];
    double pmImu2ImuGrnd[16];
    double bodyPe[6];//213
    double bodyPm[16];
    aris::dynamic::s_pe2pm(peImuGrnd2BodyGrnd,pmImuGrnd2BodyGrnd,"213");
    aris::dynamic::s_pe2pm(peImu2ImuGrnd,pmImu2ImuGrnd,"321");
    aris::dynamic::s_inv_pm(pmImuGrnd2BodyGrnd,pmBody2Imu);
    aris::dynamic::s_pm_dot_pm(pmImuGrnd2BodyGrnd,pmImu2ImuGrnd,pmBody2Imu,bodyPm);
    aris::dynamic::s_pm2pe(bodyPm,bodyPe,"213");
    for(int i=0;i<3;i++)
    {
        if(bodyPe[i+3]>PI)
        {
            bodyPe[i+3]-=2*PI;
        }
    }
 
    if(count==0)
    {
        std::fill(imuSum,imuSum+3,0);
    }
    if(count<500)
    {
        for(int i=0;i<3;i++)
        {
            imuSum[i]+=bodyPe[i+3];
            imuAvg[i]=imuSum[i]/(count+1);
        }
    }
    for(int i=0;i<3;i++)
    {
        imuInB[i]=bodyPe[i+3]-imuAvg[i];
        while(imuInB[i]>PI || imuInB[i]<-PI)
        {
            rt_printf("The situation is thought to be not happen!\n");
            if(imuInB[i]>PI)
            {
                imuInB[i]-=2*PI;
            }
            else
            {
                imuInB[i]+=2*PI;
            }
        }
    }
    //imuInB[1]+=0.005;
    imuInB[2]-=0.010;


    for(int i=0;i<2;i++)
    {
        if(imuInB[i+1]>10*PI/180)
        {
            rt_printf("ERROR! imuInB[%d]=%f too large!\n",i+1,imuInB[i+1]);
            imuInB[i+1]=10*PI/180;
        }
        else if(imuInB[i+1]<-10*PI/180)
        {
            rt_printf("ERROR! imuInB[%d]=%f too small!\n",i+1,imuInB[i+1]);
            imuInB[i+1]=-10*PI/180;
        }
    }


}

void BalanceWalk::GetAngleFromAcc(double *acc, double *angle)
{
    const double g=9.80665;
    const double angleLmt=5;//degree
    double accLmt=g*sin(angleLmt/180*PI);

    for(int i=0;i<2;i++)
    {
        if(acc[i]>accLmt)
        {
            acc[i]=accLmt;
        }
        else if(acc[i]<-accLmt)
        {
            acc[i]=-accLmt;
        }
    }

    angle[0]=-(7/5*acc[0])/g;//about z
    angle[1]=(7/5*acc[1])/g;//about x
}

void BalanceWalk::bodyEulTg(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
{
    auto &robot = static_cast<Robots::RobotBase &>(model);
    auto &param = static_cast<const aris::server::GaitParamBase &>(param_in);

    /*****IMU*****/
//    static LowpassFilter<3> imuFilter;
//    if(param.count==0)
//    {
//        imuFilter.Initialize();
//        imuFilter.SetCutFrequency(0.010,1000);
//    }
//    imuFilter.DoFilter(param.imu_data->at(0).euler,imuRaw);
    memcpy(imuRaw,param.imu_data->at(0).euler,3*sizeof(double));
    BalanceWalk::imuInit(param.count);
    memcpy(bwParam.imuRaw,imuRaw,3*sizeof(double));
    memcpy(bwParam.imuEul,imuInB,3*sizeof(double));

    double bodyPm[16];
    double bodyPe[6] {0};
    memcpy(bodyPe,imuInB,3*sizeof(double));
    aris::dynamic::s_pe2pm(bodyPe,bodyPm,"213");
    const double m {400};
    double mInG[3] {0,-m,0};
    double mInB[3];
    aris::dynamic::s_inv_pm_dot_v3(bodyPm,mInG,mInB);

    /*****FSR*****/
    memcpy(forceRaw+6*6,param.ruicong_data->at(0).force[6].fce,6*sizeof(double));
    BalanceWalk::forceInit(param.count,6);
    aris::dynamic::s_f2f(*robot.forceSensorMak().prtPm(),forceInF+6*6,bwParam.fceInB);
    //rt_printf("fceInF:%f,%f,%f,fceInB:%f,%f,%f\n",forceInF[0],forceInF[1],forceInF[2],bwParam.fceInB[0],bwParam.fceInB[1],bwParam.fceInB[2]);

    static LowpassFilter<6> filter;
    if(param.count==0)
    {
        filter.Initialize();
        filter.SetCutFrequency(0.020,1000);
    }
    filter.DoFilter(bwParam.fceInB,bwParam.fceInB_filtered);

    static LeadCompensator<6> compensator;
    if(param.count==0)
    {
        compensator.Initialize();
        compensator.SetLeadParam(1.6,1.6,0.4,1000);
    }
    compensator.Compensate(bwParam.fceInB_filtered,bwParam.fceInB_lead);

    bwParam.ballPos[0]=-bwParam.fceInB_lead[3]/mInB[1];
    bwParam.ballPos[1]=bwParam.fceInB_lead[5]/mInB[1];

    if(param.count<totalCount)
    {
        std::fill_n(bwParam.tarAcc,2,0);
        std::fill_n(bwParam.tarEul,2,0);
        std::fill_n(bwParam.actEul,3,0);
    }
    else
    {
        /*****eul*****/
        static Controller::lstPIDparam zPIDparam;
        static Controller::lstPIDparam xPIDparam;
        static Controller::lstLagParam zLagParam;
        static Controller::lstLagParam xLagParam;
        double lagStartEul[2] {0};
        double fx { 0.8 };
        double fz { 0.8 };

        if(param.count==totalCount)
        {
            zPIDparam.lstErr=bwParam.ballPos[0];
            zPIDparam.lstInt=0;
            xPIDparam.lstErr=bwParam.ballPos[1];
            xPIDparam.lstInt=0;

            zLagParam.lstFstInt=0;
            zLagParam.lstSndInt=0;
            xLagParam.lstFstInt=0;
            xLagParam.lstSndInt=0;
        }
        //if(param.count%totalCount==0)//must do this
        //{
            lagStartEul[0]=bwParam.actEul[1];
            lagStartEul[1]=bwParam.actEul[2];
        //}

        static double lst_calmanPos[2] {0};
        bwParam.ballVel[0]=(bwParam.ballPos_filtered[0]-lst_calmanPos[0])/0.001;
        bwParam.ballVel[1]=(bwParam.ballPos_filtered[1]-lst_calmanPos[1])/0.001;
        memcpy(lst_calmanPos,bwParam.ballPos_filtered,2*sizeof(double));
        doCalman();
        bwParam.tarAcc[0]=1.5*bwParam.ballPos_filtered[0]+1.5*bwParam.ballVel_filtered[0];//x
        bwParam.tarAcc[1]=1.5*bwParam.ballPos_filtered[1]+1.5*bwParam.ballVel_filtered[1];//z

//        bwParam.tarAcc[0]=Controller::doPID(zPIDparam,bwParam.ballPos[0],3,0,3,0.001);
//        bwParam.tarAcc[1]=Controller::doPID(xPIDparam,bwParam.ballPos[1],3,0,3,0.001);

//        double ratio[2];//change from 1/3 to 1 as ballPos change from 0 to 0.2m;
//        for(int i=0;i<2;i++)
//        {
//            if(fabs(bwParam.ballPos[i])>0.01)
//	    {
//                ratio[i]=1;
//            }
//            else
//            {
//                ratio[i]=0.3;
//            }
//            bwParam.tarAcc[i]*=ratio[i];
//        }

        GetAngleFromAcc(bwParam.tarAcc,bwParam.tarEul);//angle around x-axis
        bwParam.tarEul[0]-=bwParam.imuEul[1];
        bwParam.tarEul[1]-=bwParam.imuEul[2];
 
//        double tarAcc_tmp[3] {bwParam.tarAcc[1],0,bwParam.tarAcc[0]};
//        double tarEul_tmp[3];
//        GetTarEulFromAcc(tarAcc_tmp,tarEul_tmp);
//        bwParam.tarEul[0]=-tarEul_tmp[1]-(bwParam.imuEul[1]-bwParam.actEul[1]);
//        bwParam.tarEul[1]=-tarEul_tmp[2]-(bwParam.imuEul[2]-bwParam.actEul[2]);
        for(int i=0;i<2;i++)
        {
            if(bwParam.tarEul[i]>10*PI/180)
            {
                rt_printf("ERROR! tarEul[%d]=%f too large!\n",i,bwParam.tarEul[i]);
                bwParam.tarEul[i]=10*PI/180;
            }
            else if(bwParam.tarEul[i]<-10*PI/180)
            {
                rt_printf("ERROR! tarEul[%d]=%f too small!\n",i,bwParam.tarEul[i]);
                bwParam.tarEul[i]=-10*PI/180;

            }
        }

        bwParam.input[3]=1.5*bwParam.tarEul[0];
        bwParam.input[5]=1.5*bwParam.tarEul[1];

        bwParam.actEul[1]=Controller::SndOrderLag(zLagParam,0,bwParam.tarEul[0]+lagStartEul[0],2*PI*fz,0.707,0.001);
        bwParam.actEul[2]=Controller::SndOrderLag(xLagParam,0,bwParam.tarEul[1]+lagStartEul[1],2*PI*fx,0.707,0.001);

        pEB[4]=bwParam.actEul[1];
        pEB[5]=bwParam.actEul[2];
    }
}

void BalanceWalk::doCalman()
{
    double acc[3];
    GetAccFromActEul(acc,bwParam.actEul);

    const double delta_t {0.001};
    static double lst_x1[2] {0};
    static double lst_P1[2][2] {0};
    static double lst_x2[2] {0};
    static double lst_P2[2][2] {0};
    static double pre_x1[2] {0};
    static double pre_x2[2] {0};

    double mtrx_F[2][2] {1, delta_t, 0, 1};
    double input_u1[2] {delta_t*delta_t/2*acc[0], delta_t*acc[0]};//pos on x
    double input_u2[2] {delta_t*delta_t/2*acc[2], delta_t*acc[2]};//pos on z
    double noise_Q[2][2] {5e-2*5e-2, 0, 0, 9e-2*9e-2};
    double mtrx_H[2][2] {1, 0, 0, 1};
    double noise_R[2][2] {0.20*0.20, 0, 0, 0.8*0.8};
    double input_z1[2] {bwParam.ballPos[1], bwParam.ballVel[1]};
    double input_z2[2] {bwParam.ballPos[0], bwParam.ballVel[0]};

    Calman2d filter1;
    filter1.init(*mtrx_F,*noise_Q,*mtrx_H,*noise_R,input_u1,input_z1);
    filter1.doFilter(lst_x1,pre_x1,*lst_P1);

    Calman2d filter2;
    filter2.init(*mtrx_F,*noise_Q,*mtrx_H,*noise_R,input_u2,input_z2);
    filter2.doFilter(lst_x2,pre_x2,*lst_P2);

    bwParam.ballPos_filtered[1]=lst_x1[0];//pos on x
    bwParam.ballVel_filtered[1]=lst_x1[1];
    bwParam.ballPos_filtered[0]=lst_x2[0];//pos on z
    bwParam.ballVel_filtered[0]=lst_x2[1];

    bwParam.predictPos[1]=pre_x1[0];
    bwParam.predictVel[1]=pre_x1[1];
    bwParam.predictPos[0]=pre_x2[0];
    bwParam.predictVel[0]=pre_x2[1];
}

void BalanceWalk::swingLegTg(const aris::dynamic::PlanParamBase &param_in, int legID)
{
    auto &param = static_cast<const aris::server::GaitParamBase &>(param_in);

    int period_count = param.count%totalCount;
    if(isSndHalf[legID]==false)
    {
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

        //initialize the force on each leg in every cycle
        if(period_count>=(totalCount/2-100))
        {
            forceInit(period_count-totalCount/2+100,legID);
        }
        else
        {
            forceInit(1000,legID);
        }
    }
    else
    {
        forceInit(1000,legID);
    }

    if(isSndHalf[legID]==false && period_count<(int)(0.9186*totalCount)+1)
    {
        double s;
        if(period_count<(int)(5.0/6*totalCount)+1)
        {
            s=PI/2-PI/2*cos(PI*period_count/totalCount);//0-5/6*PI
        }
        else
        {
            s=(2+sqrt(3))*PI/4+PI*PI/4/totalCount*(period_count-5.0/6*totalCount);//5/6*PI-PI
        }
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
    }
    else
    {
        if(period_count<(int)(0.2226*totalCount)+1)
        {
            swingPee[3*legID+1]=swingEndPee[3*legID+1]-PI*PI*height/4*(period_count+totalCount-0.9186*totalCount)/totalCount;
        }
        else if(period_count<(int)(0.4252*totalCount)+1)
        {
            swingPee[3*legID+1]=swingEndPee[3*legID+1]-0.75*height-PI*PI*height/4*(period_count+totalCount-1.2226*totalCount)/totalCount
                    +PI*PI*PI*PI*height/16*(period_count+totalCount-1.2226*totalCount)/totalCount*(period_count+totalCount-1.2226*totalCount)/totalCount;
        }
        else if(period_count>=(int)(0.9186*totalCount)+1)
        {
            swingPee[3*legID+1]=swingEndPee[3*legID+1]-PI*PI*height/4*(period_count-0.9186*totalCount)/totalCount;
        }
        else
        {
            rt_printf("Error! Impossible period_count appears in swing phase!");
        }
    }
}

void BalanceWalk::followLegTg(const aris::dynamic::PlanParamBase &param_in, int legID)
{
    auto &param = static_cast<const aris::server::GaitParamBase &>(param_in);
    forceInit(1000,legID);
    memcpy(followPee+3*legID,followBeginPee+3*legID,3*sizeof(double));

    double accLmt {2.0};
    int decCount[3];
    int maxDecCount {0};
    double acc[3];
    for(int i=0;i<3;i++)
    {
        decCount[i]=(int)fabs(followBeginVee[3*legID+i]/accLmt*1000)+1;
        if(decCount[i]>maxDecCount)
        {
            maxDecCount=decCount[i];
        }
    }
    //if(param.count==followBeginCount[legID])
    //{
        //rt_printf("leg %d going to dec within count %d\n",legID,maxDecCount);
        //rt_printf("followRecordPee:%f, %f, %f\n",followRecordPee[3*legID],followRecordPee[3*legID+1],followRecordPee[3*legID+2]);
        //rt_printf("followBeginPee:%f, %f, %f\n",followBeginPee[3*legID],followBeginPee[3*legID+1],followBeginPee[3*legID+2]);
        //rt_printf("followBeginVee:%f, %f, %f\n",followBeginVee[3*legID],followBeginVee[3*legID+1],followBeginVee[3*legID+2]);
    //}

    for(int i=0;i<3;i++)
    {
        acc[i]=-followBeginVee[3*legID+i]/(maxDecCount*0.001);
        followPee[3*legID+i]=followBeginPee[3*legID+i]+followBeginVee[3*legID+i]*(param.count-followBeginCount[legID]+1)*0.001
                +0.5*acc[i]*(param.count-followBeginCount[legID]+1)*0.001*(param.count-followBeginCount[legID]+1)*0.001;
    }

    if(param.count-followBeginCount[legID]==maxDecCount-1)
    {
        isFollowEnd[legID]=true;
    }
}
/*
void BalanceWalk::stanceLegTg(const aris::dynamic::PlanParamBase &param_in, int legID)
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
}*/

void BalanceWalk::stanceLegTg(const aris::dynamic::PlanParamBase &param_in, int legID)
{
    auto &param = static_cast<const aris::server::GaitParamBase &>(param_in);
    forceInit(1000,legID);

    int period_count = param.count%totalCount;
    memcpy(stancePee+3*legID,stanceBeginPee+3*legID,3*sizeof(double));
    if(isSndHalf[legID]==true && period_count>=totalCount/2)
    {
        stancePee[3*legID+1]=stanceBeginPee[3*legID+1]+(planH-avgRealH)/2*(1-cos(2*PI*(period_count-totalCount/2)/totalCount));
    }

    //if(legID==3)
    {
        //rt_printf("leg:%d, stancePee:%.4f,%.4f,%.4f\n",legID,stancePee[3*legID],stancePee[3*legID+1],stancePee[3*legID+2]);
    //rt_printf("planH=%f,avgH=%f\n",planH,avgRealH);
    }
}

int BalanceWalk::balanceWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
{
    auto &robot = static_cast<Robots::RobotBase &>(model);
    auto &param = static_cast<const aris::server::GaitParamBase &>(param_in);

    EmergencyStop stoper;
    stoper.doRecord(robot);
    
    const double frcRange[6]{-10,-10,-10,-10,-10,-10};

    static aris::dynamic::FloatMarker beginMakLG0{robot.ground()};
    static aris::dynamic::FloatMarker beginMakLG1{robot.ground()};
    static aris::dynamic::FloatMarker beginMak{robot.ground()};
    if (param.count == 0)
    {
        beginMakLG0.setPrtPm(*robot.body().pm());
        beginMakLG0.update();
        beginMakLG1.setPrtPm(*robot.body().pm());
        beginMakLG1.update();
        beginMak.setPrtPm(*robot.body().pm());
        beginMak.update();
        robot.GetPee(initPee,beginMak);
        robot.GetPee(stanceBeginPee,beginMak);
        //robot.GetPee(swingEndPee,beginMak);
        std::fill_n(gaitPhase,6,GaitPhase::Stance);
        std::fill_n(isSndHalf,6,false);
        std::fill_n(isFollowEnd,6,false);
        std::fill_n(isTrans2Follow,6,false);
        std::fill_n(isPreTrans2Follow,6,false);
        std::fill_n(preTrans2FollowCount,6,0);
        std::fill_n(bwParam.input,6,0);
        inputEul_tmp[0]=0;//set yaw const 0
        //totalCount=totalCount_tmp;//param.count%totalCount illegal at count 0

        //planH=initPee[1]-0.005;
        planH=initPee[1]-0.000;
    }

    if(param.count%totalCount==0)
    {
        double bodyPe_tmp[6];
        double bodyPm_tmp[4][4];
        robot.GetPeb(bodyPe_tmp,"213");
        bodyPe_tmp[4]=0;
        bodyPe_tmp[5]=0;
        aris::dynamic::s_pe2pm(bodyPe_tmp,*bodyPm_tmp,"213");
        beginMak.setPrtPm(*bodyPm_tmp);
        beginMak.update();
        robot.GetPeb(beginPeb,beginMak,"213");

        //totalCount=totalCount_tmp;
        height=height_tmp;
        alpha=alpha_tmp;
        memcpy(inputEul,inputEul_tmp,sizeof(inputEul_tmp));

        beginXVel=endXVelinB;
        beginZVel=endZVelinB;
        beginOmega=endOmega;
        switch(walkState)
        {
        case BalanceWalkState::Init:
            beginXVel=0;
            endXVel=0;
            endXVelinB=0;
            beginZVel=0;
            endZVel=0;
            endZVelinB=0;
            beginOmega=0;
            endOmega=0;
            constFlag=false;
            break;

        case BalanceWalkState::ForwardAcc:
            endOmega=beginOmega;
            endGama=(endOmega+beginOmega)/2*totalCount*0.001;
            endXVelinB=beginXVel+sin(alpha)*forwardAcc*totalCount*0.001;
            endZVelinB=beginZVel+cos(alpha)*forwardAcc*totalCount*0.001;
            endXVel=cos(endGama)*endXVelinB+sin(endGama)*endZVelinB;
            endZVel=-sin(endGama)*endXVelinB+cos(endGama)*endZVelinB;
            constFlag=true;

            bwParam.input[0]=(endXVel-beginXVel)/(0.001*totalCount);
            bwParam.input[2]=(endZVel-beginXVel)/(0.001*totalCount);
            bwParam.input[4]=(endOmega-beginOmega)/(0.001*totalCount);
            break;

        case BalanceWalkState::TurnAcc:
            endOmega=beginOmega+turnAcc*totalCount/1000;
            endGama=(endOmega+beginOmega)/2*totalCount*0.001;
            endXVelinB=beginXVel;
            endZVelinB=beginZVel;
            endXVel=cos(endGama)*endXVelinB+sin(endGama)*endZVelinB;
            endZVel=-sin(endGama)*endXVelinB+cos(endGama)*endZVelinB;
            constFlag=true;

            bwParam.input[0]=(endXVel-beginXVel)/(0.001*totalCount);
            bwParam.input[2]=(endZVel-beginXVel)/(0.001*totalCount);
            bwParam.input[4]=(endOmega-beginOmega)/(0.001*totalCount);
            break;

        case BalanceWalkState::Const:
            endOmega=beginOmega;
            endGama=(endOmega+beginOmega)/2*totalCount*0.001;
            endXVelinB=beginXVel;
            endZVelinB=beginZVel;
            endXVel=cos(endGama)*endXVelinB+sin(endGama)*endZVelinB;
            endZVel=-sin(endGama)*endXVelinB+cos(endGama)*endZVelinB;

            bwParam.input[0]=(endXVel-beginXVel)/(0.001*totalCount);
            bwParam.input[2]=(endZVel-beginXVel)/(0.001*totalCount);
            bwParam.input[4]=(endOmega-beginOmega)/(0.001*totalCount);
            break;

        default:
            rt_printf("Error: unknown walkstate!\n");
            break;
        }
    }
    if(param.count%totalCount==(totalCount-1))
    {
        if(walkState==BalanceWalkState::ForwardAcc && constFlag==true)
        {
            walkState=BalanceWalkState::Const;
            constFlag=false;
            rt_printf("ForwardAcc finished, the coming Const Vel is:(-x)%.4f,(-z)%.4f,Omega is:%.4f\n",endXVel,endZVel,endOmega);
        }
        if(walkState==BalanceWalkState::TurnAcc && constFlag==true)
        {
            walkState=BalanceWalkState::Const;
            constFlag=false;
            rt_printf("TurnAcc finished, the coming Const Vel is:(-x)%.4f,(-z)%.4f,Omega is:%.4f\n",endXVel,endZVel,endOmega);
        }
    }

    //initialize gaitPhase to Swing at the beginning of every cycle
    if(param.count%(2*totalCount)==0)
    {
        double bodyPe_tmp[6];
        double bodyPm_tmp[4][4];
        robot.GetPeb(bodyPe_tmp,"213");
        bodyPe_tmp[4]-=bwParam.imuEul[1];
        bodyPe_tmp[5]-=bwParam.imuEul[2];
        aris::dynamic::s_pe2pm(bodyPe_tmp,*bodyPm_tmp,"213");
        beginMakLG0.setPrtPm(*bodyPm_tmp);
        beginMakLG0.update();
        for (int i=0;i<3;i++)
        {
            gaitPhase[2*i]=GaitPhase::Swing;
            robot.pLegs[2*i]->GetPee(swingBeginPee+6*i,beginMakLG0);
            isSndHalf[2*i]=false;
            isSndHalf[2*i+1]=true;
        }
    }
    else if(param.count%(2*totalCount)==totalCount)
    {
        double bodyPe_tmp[6];
        double bodyPm_tmp[4][4];
        robot.GetPeb(bodyPe_tmp,"213");
        bodyPe_tmp[4]-=bwParam.imuEul[1];
        bodyPe_tmp[5]-=bwParam.imuEul[2];
        aris::dynamic::s_pe2pm(bodyPe_tmp,*bodyPm_tmp,"213");
        beginMakLG1.setPrtPm(*bodyPm_tmp);
        beginMakLG1.update();
        for (int i=0;i<3;i++)
        {
            gaitPhase[2*i+1]=GaitPhase::Swing;
            robot.pLegs[2*i+1]->GetPee(swingBeginPee+6*i+3,beginMakLG1);
            isSndHalf[2*i]=true;
            isSndHalf[2*i+1]=false;
        }
    }

    if(param.count%totalCount==totalCount/2)
    {
        double min{0};
        double sum{0};
        int num{0};
        for(int i=0;i<6;i++)
        {
            if(gaitPhase[i]==GaitPhase::Stance)
            {
                if(stanceBeginPee[3*i+1]<min)
                {
                    min=stanceBeginPee[3*i+1];
                }
                sum+=stanceBeginPee[3*i+1];
                num++;
            }
        }
        rt_printf("num=%d\n",num);
        //avgRealH=min;
        avgRealH=sum/3;
    }

    memcpy(pEB,beginPeb,sizeof(beginPeb));
    pEB[0]=beginPeb[0]-(beginXVel*(param.count%totalCount+1)*0.001
           +0.5*(endXVel-beginXVel)/totalCount*(param.count%totalCount+1)*(param.count%totalCount+1)*0.001);
    pEB[1]=0;
    pEB[2]=beginPeb[2]-(beginZVel*(param.count%totalCount+1)*0.001
           +0.5*(endZVel-beginZVel)/totalCount*(param.count%totalCount+1)*(param.count%totalCount+1)*0.001);
    pEB[3]=bwParam.actEul[0]=beginPeb[3]+(beginOmega*(param.count%totalCount+1)*0.001
           +0.5*(endOmega-beginOmega)/totalCount*(param.count%totalCount+1)*(param.count%totalCount+1)*0.001);
    bodyEulTg(robot,param);
    for(int i=0;i<3;i++)
    {
        if(pEB[i+3]>PI)
        {
            pEB[i+3]-=2*PI;
        }
    }
    robot.SetPeb(pEB,beginMak,"213");
    memcpy(bwParam.bodyPos,pEB,3*sizeof(double));

    const double leg2fce[6] {0,1,2,5,4,3};
    for(int i=0;i<6;i++)
    {
        memcpy(forceRaw+6*i,param.ruicong_data->at(0).force[leg2fce[i]].fce,6*sizeof(double));
        switch(gaitPhase[i])
        {
        case GaitPhase::Swing:
            swingLegTg(param,i);
            if(i%2==0)
            {
                robot.pLegs[i]->SetPee(swingPee+3*i,beginMakLG0);
            }
            else
            {
                robot.pLegs[i]->SetPee(swingPee+3*i,beginMakLG1);
            }

            if(isTrans2Follow[i]==false && ((isSndHalf[i]==false && param.count%totalCount>totalCount/2) || isSndHalf[i]==true))
            {
                //detect 5 points to confirm touching ground
                int pntsNum {5};
                if(forceInF[6*i+2]<frcRange[i] && isPreTrans2Follow[i]==false)//param.force_data->at(leg2frc[i]).Fz<frcRange[i]
                {
                    isPreTrans2Follow[i]=true;
                    preTrans2FollowCount[i]=param.count;
                    //rt_printf("leg %d detects force:%.4f, going into Follow in 5 ms after count %d\n",i,forceInF[6*i+2],preTrans2FollowCount[i]);
                }
                if(forceInF[6*i+2]>frcRange[i] && isPreTrans2Follow[i]==true && param.count<preTrans2FollowCount[i]+pntsNum)
                {
                    isPreTrans2Follow[i]=false;
                    //rt_printf("leg %d gets fake touching signal at count %d\n",i,preTrans2FollowCount[i]);
                    preTrans2FollowCount[i]=0;
                }

                //touching ground confirmed, transfer into follow
                if(isPreTrans2Follow[i]==true && param.count==preTrans2FollowCount[i]+pntsNum-1)
                {
                    if(i%2==0)
                    {
                        robot.pLegs[i]->GetPee(followRecordPee+3*i,beginMakLG0);
                    }
                    else
                    {
                        robot.pLegs[i]->GetPee(followRecordPee+3*i,beginMakLG1);
                    }
                }
                if(isPreTrans2Follow[i]==true && param.count>(preTrans2FollowCount[i]+pntsNum-1))
                {
                    //rt_printf("leg %d detects force:%.4f, transfer into Follow at count %d\n",i,forceInF[6*i+2],param.count);
                    isTrans2Follow[i]=true;
                    gaitPhase[i]=GaitPhase::Follow;
                    if(i%2==0)
                    {
                        robot.pLegs[i]->GetPee(followBeginPee+3*i,beginMakLG0);
                    }
                    else
                    {
                        robot.pLegs[i]->GetPee(followBeginPee+3*i,beginMakLG1);
                    }
                    for(int j=0;j<3;j++)
                    {
                        followBeginVee[3*i+j]=(followBeginPee[3*i+j]-followRecordPee[3*i+j])*1000;
                    }
                    followBeginCount[i]=param.count+1;
                }
            }

            if(isSndHalf[i]==true && param.count%totalCount==(int)(0.4252*totalCount) && isTrans2Follow[i]==false)
            {
                gaitPhase[i]=GaitPhase::Stance;
                if(i%2==0)
                {
                    robot.pLegs[i]->GetPee(stanceBeginPee+3*i,beginMakLG0);
                }
                else
                {
                    robot.pLegs[i]->GetPee(stanceBeginPee+3*i,beginMakLG1);
                }
            }

            break;

        case GaitPhase::Follow:
            followLegTg(param,i);
            if(i%2==0)
            {
                robot.pLegs[i]->SetPee(followPee+3*i,beginMakLG0);
            }
            else
            {
                robot.pLegs[i]->SetPee(followPee+3*i,beginMakLG1);
            }

            if(isFollowEnd[i]==true)
            {
                gaitPhase[i]=GaitPhase::Stance;
                if(i%2==0)
                {
                    robot.pLegs[i]->GetPee(stanceBeginPee+3*i,beginMakLG0);
                }
                else
                {
                    robot.pLegs[i]->GetPee(stanceBeginPee+3*i,beginMakLG1);
                }
            }

            break;

        case GaitPhase::Stance:
            stanceLegTg(param,i);
            if(i%2==0)
            {
                robot.pLegs[i]->SetPee(stancePee+3*i,beginMakLG0);
            }
            else
            {
                robot.pLegs[i]->SetPee(stancePee+3*i,beginMakLG1);
            }

            if(param.count%totalCount==totalCount-1 && isSndHalf[i]==true)
            {
                isPreTrans2Follow[i]=false;
                isTrans2Follow[i]=false;
                isFollowEnd[i]=false;
            }

            break;

        default:
            rt_printf("Error: unknown gaitphase!\n");
            break;
        }
    }

    for(int i=0;i<6;i++)
    {
        bwParam.phase[i]=gaitPhase[i];
        bwParam.legForce[i]=forceInF[6*i+2];
    }
    robot.GetPeb(bwParam.PebInG,robot.ground(),"213");
    for(int i=0;i<3;i++)
    {
        if(bwParam.PebInG[i+3]>PI)
        {
            bwParam.PebInG[i+3]-=2*PI;
        }
    }
    robot.GetPee(bwParam.footPos,beginMak);
    robot.GetPin(bwParam.Pin);
    bwPipe.sendToNrt(bwParam);

    if(walkState==BalanceWalkState::Stop)
    {
        int ret=stoper.stop(robot,param.count,3.2);
        return ret;
    }
    else
    {
        return 1;
    }
}

void BalanceWalk::recordData()
{
    bwThread = std::thread([&]()
    {
        struct BalanceWalkParam param;
        static std::fstream fileGait;
        std::string name = aris::core::logFileName();
        name.replace(name.rfind("log.txt"), std::strlen("BalanceWalkData.txt"), "BalanceWalkData.txt");
        fileGait.open(name.c_str(), std::ios::out | std::ios::trunc);

        long long count = -1;
        while (1)
        {
            bwPipe.recvInNrt(param);

            fileGait << ++count << " ";
            for (int i=0;i<6;i++)
            {
                fileGait << param.legForce[i] << "  ";
            }
            for (int i=0;i<6;i++)
            {
                fileGait << param.fceInB[i] << "  ";
            }
            for (int i=0;i<6;i++)
            {
                fileGait << param.fceInB_filtered[i] << "  ";
            }
            for (int i=0;i<6;i++)
            {
                fileGait << param.fceInB_lead[i] << "  ";
            }
            for (int i=0;i<6;i++)
            {
                fileGait << param.phase[i] << "  ";
            }
            for (int i=0;i<18;i++)
            {
                fileGait << param.footPos[i] << "  ";
            }
            for (int i=0;i<18;i++)
            {
                fileGait << param.Pin[i] << "  ";
            }
            for (int i=0;i<6;i++)
            {
                fileGait << param.PebInG[i] << "  ";
            }
            for (int i=0;i<6;i++)
            {
                fileGait << param.input[i] << "  ";
            }
            for(int i=0;i<2;i++)
            {
                fileGait << param.ballPos[i] << "  ";
                fileGait << param.ballPos_filtered[i] << "  ";
                fileGait << param.ballVel[i] << "  ";
                fileGait << param.ballVel_filtered[i] << "  ";
            }
            for(int i=0;i<2;i++)
            {
                fileGait << param.tarAcc[i] << "  ";
            }
            for(int i=0;i<2;i++)
            {
                fileGait << param.tarEul[i] << "  ";
            }
            for(int i=0;i<3;i++)
            {
                fileGait << param.actEul[i] << "  ";
            }
            for(int i=0;i<3;i++)
            {
                fileGait << param.imuRaw[i] << "  ";
            }
            for(int i=0;i<3;i++)
            {
                fileGait << param.imuEul[i] << "  ";
            }
            fileGait << std::endl;
        }

        fileGait.close();
    });
}

/*
{ 0,0,0,PI,PI/2,PI*3/2 }

<bw>
    <bw_param type="group">
        <state_param type="unique">
            <init abbreviation="i"/>
            <forwardAcc abbreviation="f" type="double" default="0"/>
            <betaAcc abbreviation="b" type="double" default="0"/>
            <stop abbreviation="s"/>
        </state_param>
        <totalCount abbreviation="t" type="int" default="1000"/>
        <height abbreviation="h" type="double" default="0.05"/>
        <roll abbreviation="r" type="double" default="0"/>
        <pitch abbreviation="p" type="double" default="0"/>
        <alpha abbreviation="a" type="double" default="0"/>
    </bw_param>
</bw>

rs.addCmd("bw",BalanceWalk::parseBalanceWalk,BalanceWalk::balanceWalk);
BalanceWalk::recordData();
*/
