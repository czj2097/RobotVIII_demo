#include "RealTimeOptimal.h"

RealTimeOptimal::RealTimeOptimal(){}
RealTimeOptimal::~RealTimeOptimal(){}

/*find gait with maxVel in joint space by iteration*/
void RealTimeOptimal::GetOptimalPTPMotion(double startP, double endP, double startV, double endV, double *a, double *t)
{
    int k1 = endP > startP ? 1 : -1;
    double v0 {0};
    double vm {0};
    double vt {0};
    if(startV*k1<0)
    {
        v0=-startV;
        a[0]=k1*aLmt;
        t[0]=-2*startV/(k1*aLmt);
    }
    else
    {
        v0=startV;
        a[0]=0;
        t[0]=0;
    }
    if(endV*k1<0)
    {
        vt=-endV;
        a[4]=-k1*aLmt;
        t[4]=2*endV/(-k1*aLmt);
    }
    else
    {
        vt=endV;
        a[4]=0;
        t[4]=0;
    }

    printf("startP=%.4f, endP=%.4f, startV=%.4f, endV=%.4f\n",startP,endP,startV,endV);

    double s1,s2,s3;
    s1=(vLmt*vLmt-v0*v0)/(2*k1*aLmt);
    s3=(vLmt*vLmt-vt*vt)/(2*k1*aLmt);
    s2=endP-startP-s1-s3;

    if(s2*k1>0)
    {
        printf("s2_const exist, the s2 is %.4f\n",s2);
        a[1]=k1*aLmt;
        t[1]=(k1*vLmt-v0)/(k1*aLmt);//dec
        a[2]=0;
        t[2]=s2/(k1*vLmt);//const
        a[3]=-k1*aLmt;
        t[3]=(k1*vLmt-vt)/(k1*aLmt);//acc
    }
    else
    {
        if(k1*aLmt*(endP-startP)+0.5*v0*v0+0.5*vt*vt<0)
        {
            printf("WARNING!!! sqrt() applied to a negative data.\n");
        }
        double judgeP=startP+fabs(vt*vt-v0*v0)/(2*k1*aLmt);
        vm=k1*sqrt(k1*aLmt*(endP-startP)+0.5*v0*v0+0.5*vt*vt);
        if((judgeP-endP)*(endP-startP)>0)
        {
            printf("acc in current dir is not right, acc in the opposite dir\n");
            vm=k1*sqrt(-k1*aLmt*(endP-startP)+0.5*v0*v0+0.5*vt*vt);
            printf("s2_const dont exist,the max vel is %.4f\n",vm);
            if((vm-v0)*k1>0 || (vm-vt)*k1>0)
            {
                vm=-k1*sqrt(-k1*aLmt*(endP-startP)+0.5*v0*v0+0.5*vt*vt);
                printf("vm not right, modified to %.4f\n",vm);
            }
            a[1]=-k1*aLmt;
            t[1]=(vm-v0)/(-k1*aLmt);
            a[2]=0;
            t[2]=0;
            a[3]=k1*aLmt;
            t[3]=(vm-vt)/(-k1*aLmt);//acc
        }
        else
        {
            vm=k1*sqrt(k1*aLmt*(endP-startP)+0.5*v0*v0+0.5*vt*vt);
            printf("s2_const dont exist,the max vel is %f\n",vm);
            a[1]=k1*aLmt;
            t[1]=(vm-v0)/(k1*aLmt);
            a[2]=0;
            t[2]=0;
            a[3]=-k1*aLmt;
            t[3]=(vm-vt)/(k1*aLmt);//acc
        }
    }

    if(t[0]<-1e-6 || t[1]<-1e-6 || t[2]<-1e-6 || t[3]<-1e-6 || t[4]<-1e-6)
    {
        printf("\nWARNING!!! t<0 : v0=%.4f, vm=%.4f, vt=%.4f\n\n",v0,vm,vt);
    }

//    printf("s:%.4f,%.4f,%.4f\n",s1,s2,s3);
}

void RealTimeOptimal::GetTraj(int screwID, int startCount, int totalCount, double startP, double endP, double startV, double endV, double *a, double *t)
{
    double v[5] {startV};
    double s[5] {startP};
    for(int i=1;i<5;i++)
    {
        v[i]=v[i-1]+a[i-1]*t[i-1];
        s[i]=s[i-1]+v[i-1]*t[i-1]+0.5*a[i-1]*t[i-1]*t[i-1];
    }
    printf("t:%.4f,%.4f,%.4f,%.4f,%.4f\n",t[0],t[1],t[2],t[3],t[4]);
    printf("a:%.4f,%.4f,%.4f,%.4f,%.4f\n",a[0],a[1],a[2],a[3],a[4]);
    printf("v:%.4f,%.4f,%.4f,%.4f,%.4f\n",v[0],v[1],v[2],v[3],v[4]);
    printf("s:%.4f,%.4f,%.4f,%.4f,%.4f\n\n",s[0],s[1],s[2],s[3],s[4]);

    for (int i=0;i<totalCount;i++)
    {
        if(i*0.001<t[0])
        {
            vIn[i+startCount][screwID]=v[0]+a[0]*i*0.001;
            pIn[i+startCount][screwID]=s[0]+v[0]*i*0.001+0.5*a[0]*i*0.001*i*0.001;
        }
        else if(i*0.001<(t[0]+t[1]))
        {
            vIn[i+startCount][screwID]=v[1]+a[1]*(i*0.001-t[0]);
            pIn[i+startCount][screwID]=s[1]+v[1]*(i*0.001-t[0])+0.5*a[1]*(i*0.001-t[0])*(i*0.001-t[0]);
        }
        else if(i*0.001<(t[0]+t[1]+t[2]))
        {
            vIn[i+startCount][screwID]=v[2]+a[2]*(i*0.001-t[0]-t[1]);
            pIn[i+startCount][screwID]=s[2]+v[2]*(i*0.001-t[0]-t[1])+0.5*a[2]*(i*0.001-t[0]-t[1])*(i*0.001-t[0]-t[1]);
        }
        else if(i*0.001<(t[0]+t[1]+t[2]+t[3]))
        {
            vIn[i+startCount][screwID]=v[3]+a[3]*(i*0.001-t[0]-t[1]-t[2]);
            pIn[i+startCount][screwID]=s[3]+v[3]*(i*0.001-t[0]-t[1]-t[2])+0.5*a[3]*(i*0.001-t[0]-t[1]-t[2])*(i*0.001-t[0]-t[1]-t[2]);
        }
        else if(i*0.001<(t[0]+t[1]+t[2]+t[3]+t[4]))
        {
            vIn[i+startCount][screwID]=v[4]+a[4]*(i*0.001-t[0]-t[1]-t[2]-t[3]);
            pIn[i+startCount][screwID]=s[4]+v[4]*(i*0.001-t[0]-t[1]-t[2]-t[3])+0.5*a[4]*(i*0.001-t[0]-t[1]-t[2]-t[3])*(i*0.001-t[0]-t[1]-t[2]-t[3]);
        }
        else if(i*0.001<(t[0]+t[1]+t[2]+t[3]+t[4]+t[5]))
        {
            vIn[i+startCount][screwID]=v[5]+a[5]*(i*0.001-t[0]-t[1]-t[2]-t[3]-t[4]);
            pIn[i+startCount][screwID]=s[5]+v[5]*(i*0.001-t[0]-t[1]-t[2]-t[3]-t[4])+0.5*a[5]*(i*0.001-t[0]-t[1]-t[2]-t[3]-t[4])*(i*0.001-t[0]-t[1]-t[2]-t[3]-t[4]);
        }
//        else
//        {
//            vIn[i+startCount][screwID]=endV;
//            pIn[i+startCount][screwID]=endP;
//        }
    }
}

void RealTimeOptimal::GetBezier(int screwID, double startP, double endP, double startV, double endV)
{
    double c0,c1,c2,c3,c4;

}

void RealTimeOptimal::screwInterpolationTraj()
{
    Robots::RobotTypeI rbt;
    rbt.loadXml("/home/hex/Desktop/mygit/RobotVIII_demo/resource/RobotEDU2.xml");

    timeval tpstart,tpend;
    float tused;
    gettimeofday(&tpstart,NULL);

    double totalTime{1.5};
//    double pEE[18];
//    double pEE_B[18];
//    double aEE[18];
    double ratio {0};//control point, from middle to edge [0,1]

//    int keyPointNum{3};
    double keyPee_B[3][18];
    double keyPin[3][18];
    double keyVin[3][18];
    double dir[3] {1,-1,1};

    rbt.SetPeb(initPeb);
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
    double totalTmax {0};
    double totalT[18] {0};
    double totalT1[18] {0};
    double totalT2[18] {0};
    int totalCount1[18] {0};
    int totalCount2[18] {0};
    double t1[18][5] {0};
    double t2[18][5] {0};
    double a1[18][5] {0};
    double a2[18][5] {0};

    while(e>=0.0001)
    {
        for (int i=0;i<3;i++)
        {
            double maxVel[18];
            double Vbody[3] {0,0,stepD/totalTime/2};
            rbt.SetPee(*keyPee_B+18*i,rbt.body());
            for(int j=0;j<6;j++)
            {
                double direction[3] {0,0,dir[i]};
                FastWalk::maxCal(rbt,direction,j,vLmt,maxVel+3*j);
            }
            rbt.SetVee(maxVel,rbt.body());
            if(i!=1)
            {
                for(int j=0;j<6;j++)
                {
                    rbt.pLegs[j]->SetVee(Vbody,rbt.body());
                }
            }
            rbt.GetPin(*keyPin+18*i);
            rbt.GetVin(*keyVin+18*i);

            printf("keyPin: %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f",
                   keyPin[i][0],keyPin[i][1],keyPin[i][2],keyPin[i][3],keyPin[i][4],keyPin[i][5],keyPin[i][6],keyPin[i][7],keyPin[i][8],
                   keyPin[i][9],keyPin[i][10],keyPin[i][11],keyPin[i][12],keyPin[i][13],keyPin[i][14],keyPin[i][15],keyPin[i][16],keyPin[i][17]);
        }

        totalTmax=0;
        std::fill_n(*pIn,3000*18,0);
        std::fill_n(*vIn,3000*18,0);

        for (int i=0;i<18;i++)
        {
            GetOptimalPTPMotion(keyPin[0][i],keyPin[1][i],keyVin[0][i],keyVin[1][i],*a1+5*i,*t1+5*i);
            GetOptimalPTPMotion(keyPin[1][i],keyPin[2][i],keyVin[1][i],keyVin[2][i],*a2+5*i,*t2+5*i);

            std::fill_n(totalT1,18,0);
            std::fill_n(totalT2,18,0);
            for(int j=0;j<5;j++)
            {
                totalT1[i]+=t1[i][j];
                totalT2[i]+=t2[i][j];
            }
            totalCount1[i]=(int)(totalT1[i]*1000)+1;
            totalCount2[i]=(int)(totalT2[i]*1000)+1;

            totalT[i]=totalT1[i]+totalT2[i];

            if(totalT[i]>totalTmax)
            {
                totalTmax=totalT[i];
                k=i;
            }

            GetTraj(i,0,totalCount1[i],keyPin[0][i],keyPin[1][i],keyVin[0][i],keyVin[1][i],*a1+5*i,*t1+5*i);
            GetTraj(i,totalCount1[i],totalCount2[i],keyPin[1][i],keyPin[2][i],keyVin[1][i],keyVin[2][i],*a2+5*i,*t2+5*i);
        }

        e=fabs(totalTmax-totalTime);
        totalTime=totalTmax;

        printf("%.4f,screwID:%d\n",totalTmax,k);
        printf("\n");
    }

    gettimeofday(&tpend,NULL);
    tused=tpend.tv_sec-tpstart.tv_sec+(double)(tpend.tv_usec-tpstart.tv_usec)/1000000;
    printf("UsedTime:%f\n",tused);

//    totalTime=1;
//    int totalCount=round(totalTime*1000);
//    double vIn[3000][18];
//    double pIn[3000][18];
//    double vIn_adjust[totalCount][18];
//    double pIn_adjust[totalCount][18];
//    double pEE_adjust[totalCount][18];

//    for (int i=0;i<18;i++)
//    {
//        s1[i]=keyVin[0][i]*t1[i]-0.5*aLmt*t1[i]*t1[i];//shift in phase 1 (vector with direction)
//        s3[i]=-0.5*aLmt*t3[i]*t3[i];
//        s2[i]=keyPin[1][i]-keyPin[0][i]-s1[i]-s3[i];
//        s4[i]=0.5*aLmt*t4[i]*t4[i];
//        s6[i]=keyVin[2][i]*t6[i]+0.5*aLmt*t6[i]*t6[i];
//        s5[i]=keyPin[2][i]-keyPin[1][i]-s4[i]-s6[i];
//    }

//    for (int i=0;i<18;i++)
//    {
//        printf("t:%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",t1[i],t2[i],t3[i],t4[i],t5[i],t6[i],totalT[i]);
//        printf("s:%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",s1[i],s2[i],s3[i],s4[i],s5[i],s6[i]);
//    }

//    for (int i=0;i<3000;i++)
//    {
//        for (int j=0;j<18;j++)
//        {
//            if(((double)i/1000)<t1[0][j])
//            {
//                vIn[i][j]=keyVin[0][j]+a1[0][j]*i/1000;
//                pIn[i][j]=keyPin[0][j]+keyVin[0][j]*i/1000+0.5*a1[0][j]*i/1000*i/1000;
//            }
//            else if(((double)i/1000)<(t1[0][j]+t1[1][j]))
//            {
//                vIn[i][j]=-vLmt;
//                pIn[i][j]=keyPin[0][j]+s1[j]-vLmt*((double)i/1000-t1[j]);
//            }
//            else if(((double)i/1000)<(t1[j]+t2[j]+t3[j]+t4[j]))
//            {
//                vIn[i][j]=keyVin[0][j]-aLmt*t1[j]+aLmt*((double)i/1000-t1[j]-t2[j]);
//                pIn[i][j]=keyPin[0][j]+s1[j]+s2[j]+(keyVin[0][j]-aLmt*t1[j])*((double)i/1000-t1[j]-t2[j])+0.5*aLmt*((double)i/1000-t1[j]-t2[j])*((double)i/1000-t1[j]-t2[j]);
//            }
//            else if(((double)i/1000)<(t1[j]+t2[j]+t3[j]+t4[j]+t5[j]))
//            {
//                vIn[i][j]=vLmt;
//                pIn[i][j]=keyPin[0][j]+s1[j]+s2[j]+s3[j]+s4[j]+vLmt*((double)i/1000-(t1[j]+t2[j]+t3[j]+t4[j]));
//            }
//            else if(((double)i/1000)<(t1[j]+t2[j]+t3[j]+t4[j]+t5[j]+t6[j]))
//            {
//                vIn[i][j]=keyVin[0][j]-aLmt*t1[j]+aLmt*(t3[j]+t4[j])-aLmt*((double)i/1000-t1[j]-t3[j]-t4[j]);
//                pIn[i][j]=keyPin[0][j]+s1[j]+s2[j]+s3[j]+s4[j]+s5[j]+(keyVin[0][j]-aLmt*t1[j]+aLmt*(t3[j]+t4[j]))*((double)i/1000-t1[j]-t2[j]-t3[j]-t4[j]-t5[j])-0.5*aLmt*((double)i/1000-t1[j]-t2[j]-t3[j]-t4[j]-t5[j])*((double)i/1000-t1[j]-t2[j]-t3[j]-t4[j]-t5[j]);
//            }
//            else
//            {
//                vIn[i][j]=keyVin[2][j];
//                pIn[i][j]=keyPin[2][j];
//            }
//        }
//    }

//    for (int i=0;i<totalCount;i++)
//    {
//        for (int j=0;j<18;j++)
//        {
//            if(((double)i/1000)<(t1[j]*totalTime/totalT[j]))
//            {
//                pIn_adjust[i][j]=keyPin[0][j]+keyVin[0][j]*((double)i/1000*totalT[j]/totalTime)-0.5*aLmt*((double)i/1000*totalT[j]/totalTime)*((double)i/1000*totalT[j]/totalTime);
//            }
//            else if(((double)i/1000)<((t1[j]+t2[j])*totalTime/totalT[j]))
//            {
//                pIn_adjust[i][j]=keyPin[0][j]+s1[j]-vLmt*((double)i/1000*totalT[j]/totalTime-t1[j]);
//            }
//            else if(((double)i/1000)<((t1[j]+t2[j]+t3[j]+t4[j])*totalTime/totalT[j]))
//            {
//                pIn_adjust[i][j]=keyPin[0][j]+s1[j]+s2[j]+(keyVin[0][j]-aLmt*t1[j])*((double)i/1000*totalT[j]/totalTime-t1[j]-t2[j])+0.5*aLmt*((double)i/1000*totalT[j]/totalTime-t1[j]-t2[j])*((double)i/1000*totalT[j]/totalTime-t1[j]-t2[j]);
//            }
//            else if(((double)i/1000)<((t1[j]+t2[j]+t3[j]+t4[j]+t5[j])*totalTime/totalT[j]))
//            {
//                pIn_adjust[i][j]=keyPin[0][j]+s1[j]+s2[j]+s3[j]+s4[j]+vLmt*((double)i/1000*totalT[j]/totalTime-t1[j]-t2[j]-t3[j]-t4[j]);
//            }
//            else if(((double)i/1000)<((t1[j]+t2[j]+t3[j]+t4[j]+t5[j]+t6[j])*totalTime/totalT[j]))
//            {
//                pIn_adjust[i][j]=keyPin[0][j]+s1[j]+s2[j]+s3[j]+s4[j]+s5[j]+(keyVin[0][j]-aLmt*t1[j]+aLmt*(t3[j]+t4[j]))*((double)i/1000*totalT[j]/totalTime-t1[j]-t2[j]-t3[j]-t4[j]-t5[j])-0.5*aLmt*((double)i/1000*totalT[j]/totalTime-t1[j]-t2[j]-t3[j]-t4[j]-t5[j])*((double)i/1000*totalT[j]/totalTime-t1[j]-t2[j]-t3[j]-t4[j]-t5[j]);
//            }
//            else
//            {
//                pIn_adjust[i][j]=keyPin[2][j];
//            }
//        }

//        rbt.SetPin(*pIn_adjust+18*i);
//        rbt.GetPee(*pEE_adjust+18*i,rbt.body());
//    }

    aris::dynamic::dlmwrite("./vIn.txt",*vIn,3000,18);
    aris::dynamic::dlmwrite("./pIn.txt",*pIn,3000,18);
//    aris::dynamic::dlmwrite("./pIn_adjust.txt",*pIn_adjust,totalCount,18);
//    aris::dynamic::dlmwrite("./pEE_adjust.txt",*pEE_adjust,totalCount,18);

}

NormalGait::WalkState JointSpaceWalk::walkState;
double JointSpaceWalk::bodyAcc;
double JointSpaceWalk::bodyDec;
int JointSpaceWalk::totalCount;
double JointSpaceWalk::height;
double JointSpaceWalk::beta;
double JointSpaceWalk::beginPee[18];
double JointSpaceWalk::beginVel;
double JointSpaceWalk::endPee[18];
double JointSpaceWalk::endVel;
double JointSpaceWalk::distance;
NormalGait::GaitPhase JointSpaceWalk::gaitPhase[6];//swing true, stance false
bool JointSpaceWalk::constFlag;

JointSpaceWalk::JointSpaceWalk()
{
}
JointSpaceWalk::~JointSpaceWalk()
{
}

void JointSpaceWalk::parseJointSpaceFastWalk(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg)
{
    JointSpaceWalkParam param;

    for(auto &i:params)
    {
        if(i.first=="init")
        {
            walkState=NormalGait::WalkState::Init;
            msg.copyStruct(param);
        }
        else if(i.first=="acc")
        {
            bodyAcc=stod(i.second);
            walkState=NormalGait::WalkState::ForwardAcc;
        }
        else if(i.first=="stop")
        {
            walkState=NormalGait::WalkState::Stop;
        }
        else if(i.first=="totalCount")
        {
            totalCount=stoi(i.second);
        }
        else if(i.first=="height")
        {
            height=stod(i.second);
        }
        else if(i.first=="beta")
        {
            beta=stod(i.second);
        }
        else
        {
            std::cout<<"parse failed"<<std::endl;
        }
    }

    std::cout<<"finished parse"<<std::endl;
}

void JointSpaceWalk::swingLegTg(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in, int legID)
{
    auto &robot = static_cast<Robots::RobotBase &>(model);
    auto &param = static_cast<const JointSpaceWalkParam &>(param_in);

    double ratio{0.7};//control point, from middle to edge [0,1]
    double stepH=height;
    double stepD=distance;
    double vLmt{0.9};
    double aLmt{3.2};

    int keyPointNum{3};
    double keyPee_B[keyPointNum][3];
    double keyVee_B[keyPointNum][3];
    std::fill_n(*keyVee_B,keyPointNum*3,0);
    double keyPin[keyPointNum][3];
    double keyVin[keyPointNum][3];

    double halfGaitTime=(double)(param.count%totalCount)/1000;
    double totalTime=(double)totalCount/1000;
    double totalT[3];
    double t1[3],t2[3],t3[3],t4[3],t5[3],t6[3];
    double s1[3],s2[3],s3[3],s4[3],s5[3],s6[3];
    double pIn_adjust[3];

    //keyPee_B
    keyPee_B[0][0]=beginPee[3*legID];
    keyPee_B[0][1]=beginPee[3*legID+1];
    keyPee_B[0][2]=beginPee[3*legID+2];

    keyPee_B[1][0]=beginPee[3*legID];
    keyPee_B[1][1]=beginPee[3*legID+1]+stepH;
    if(legID==0 || legID==3)
    {
        keyPee_B[1][2]=beginPee[3*legID+2]-stepD/2*(1-ratio);
    }
    else if (legID==2 || legID==5)
    {
        keyPee_B[1][2]=beginPee[3*legID+2]-stepD/2*(1+ratio);
    }
    else
    {
        keyPee_B[1][2]=beginPee[3*legID+2]-stepD/2;
    }

    keyPee_B[2][0]=beginPee[3*legID];
    keyPee_B[2][1]=beginPee[3*legID+1];
    keyPee_B[2][2]=beginPee[3*legID+2]-stepD;

    //keyVee_B
    keyVee_B[0][2]=beginVel;
    keyVee_B[2][2]=endVel;

    //keyPin & keyVin
    robot.pLegs[legID]->SetPee(*keyPee_B);
    robot.pLegs[legID]->SetVee(*keyVee_B);
    robot.pLegs[legID]->GetPin(*keyPin);
    robot.pLegs[legID]->GetVin(*keyVin);

    robot.pLegs[legID]->SetPee(*keyPee_B+3);
    robot.pLegs[legID]->GetPin(*keyPin+3);

    robot.pLegs[legID]->SetPee(*keyPee_B+6);
    robot.pLegs[legID]->SetVee(*keyVee_B+6);
    robot.pLegs[legID]->GetPin(*keyPin+6);
    robot.pLegs[legID]->GetVin(*keyVin+6);

    for (int i=0;i<3;i++)
    {
        //cal t1[i]~t6[i]
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

        //cal s1[i]~s6[i]
        s1[i]=keyVin[0][i]*t1[i]-0.5*aLmt*t1[i]*t1[i];//shift in phase 1 ( vector with direction)
        s3[i]=-0.5*aLmt*t3[i]*t3[i];
        s2[i]=keyPin[1][i]-keyPin[0][i]-s1[i]-s3[i];
        s4[i]=0.5*aLmt*t4[i]*t4[i];
        s6[i]=keyVin[2][i]*t6[i]+0.5*aLmt*t6[i]*t6[i];
        s5[i]=keyPin[2][i]-keyPin[1][i]-s4[i]-s6[i];

        //scaling to adjust to the totalCount
        if(halfGaitTime<(t1[i]*totalTime/totalT[i]))
        {
            pIn_adjust[i]=keyPin[0][i]+keyVin[0][i]*(halfGaitTime*totalT[i]/totalTime)-0.5*aLmt*(halfGaitTime*totalT[i]/totalTime)*(halfGaitTime*totalT[i]/totalTime);
        }
        else if(halfGaitTime<((t1[i]+t2[i])*totalTime/totalT[i]))
        {
            pIn_adjust[i]=keyPin[0][i]+s1[i]-vLmt*(halfGaitTime*totalT[i]/totalTime-t1[i]);
        }
        else if(halfGaitTime<((t1[i]+t2[i]+t3[i]+t4[i])*totalTime/totalT[i]))
        {
            pIn_adjust[i]=keyPin[0][i]+s1[i]+s2[i]+(keyVin[0][i]-aLmt*t1[i])*(halfGaitTime*totalT[i]/totalTime-t1[i]-t2[i])+0.5*aLmt*(halfGaitTime*totalT[i]/totalTime-t1[i]-t2[i])*(halfGaitTime*totalT[i]/totalTime-t1[i]-t2[i]);
        }
        else if(halfGaitTime<((t1[i]+t2[i]+t3[i]+t4[i]+t5[i])*totalTime/totalT[i]))
        {
            pIn_adjust[i]=keyPin[0][i]+s1[i]+s2[i]+s3[i]+s4[i]+vLmt*(halfGaitTime*totalT[i]/totalTime-t1[i]-t2[i]-t3[i]-t4[i]);
        }
        else if(halfGaitTime<((t1[i]+t2[i]+t3[i]+t4[i]+t5[i]+t6[i])*totalTime/totalT[i]))
        {
            pIn_adjust[i]=keyPin[0][i]+s1[i]+s2[i]+s3[i]+s4[i]+s5[i]+(keyVin[0][i]-aLmt*t1[i]+aLmt*(t3[i]+t4[i]))*(halfGaitTime*totalT[i]/totalTime-t1[i]-t2[i]-t3[i]-t4[i]-t5[i])-0.5*aLmt*(halfGaitTime*totalT[i]/totalTime-t1[i]-t2[i]-t3[i]-t4[i]-t5[i])*(halfGaitTime*totalT[i]/totalTime-t1[i]-t2[i]-t3[i]-t4[i]-t5[i]);
        }
        else
        {
            pIn_adjust[i]=keyPin[2][i];
        }
    }

    robot.pLegs[legID]->SetPin(pIn_adjust);
}

void JointSpaceWalk::stanceLegTg(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in, int legID)
{
    auto &robot = static_cast<Robots::RobotBase &>(model);
    auto &param = static_cast<const JointSpaceWalkParam &>(param_in);

    double pEE_B[3];
    double halfGaitTime=(double)(param.count%totalCount)/1000;

    pEE_B[0]=beginPee[3*legID];
    pEE_B[1]=beginPee[3*legID+1];
    pEE_B[2]=beginPee[3*legID+2]+(beginVel*halfGaitTime+0.5*(endVel-beginVel)/totalCount*1000*halfGaitTime*halfGaitTime);

    robot.pLegs[legID]->SetPee(pEE_B,robot.body());
}

int JointSpaceWalk::jointSpaceFastWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
{
    auto &robot = static_cast<Robots::RobotBase &>(model);
    auto &param = static_cast<const JointSpaceWalkParam &>(param_in);

    FastWalk::outputParam OPP;

    if(param.count%(2*totalCount)<totalCount)
    {
        for (int i=0;i<3;i++)
        {
            gaitPhase[2*i]=NormalGait::GaitPhase::Swing;//swing
            gaitPhase[2*i+1]=NormalGait::GaitPhase::Stance;//stance
        }
    }
    else
    {
        for (int i=0;i<3;i++)
        {
            gaitPhase[2*i]=NormalGait::GaitPhase::Stance;//stance
            gaitPhase[2*i+1]=NormalGait::GaitPhase::Swing;//swing
        }
    }

    if(param.count%totalCount==(totalCount-1))
    {
        if(walkState==NormalGait::WalkState::ForwardAcc && constFlag==true)
        {
            walkState=NormalGait::WalkState::Const;
            constFlag=false;
            rt_printf("Acc finished, the coming Const Vel is: %.4f\n",endVel);
        }
    }
    if(param.count%totalCount==0)
    {
        beginVel=endVel;
        memcpy(beginPee,endPee,sizeof(endPee));

        switch(walkState)
        {
        case NormalGait::WalkState::Init:
            robot.GetPee(beginPee,robot.body());
            robot.GetPee(endPee,robot.body());
            beginVel=0;
            endVel=0;
            constFlag=false;
            break;

        case NormalGait::WalkState::ForwardAcc:
            endVel=beginVel+bodyAcc*totalCount/1000;
            constFlag=true;
            break;

        case NormalGait::WalkState::Const:
            endVel=beginVel;
            break;
        }
        distance=(beginVel+endVel)/2*totalCount/1000;
        for (int i=0;i<6;i++)
        {
            endPee[3*i]=beginPee[3*i];
            endPee[3*i+1]=beginPee[3*i+1];
            if(gaitPhase[i]==NormalGait::GaitPhase::Swing)
            {
                endPee[3*i+2]=beginPee[3*i+2]-distance;
            }
            else
            {
                endPee[3*i+2]=beginPee[3*i+2]+distance;
            }
        }
    }

    double initPeb[6]{0,0,0,0,0,0};
    double initVeb[6]{0,0,0,0,0,0};
    robot.SetPeb(initPeb);
    robot.SetVb(initVeb);

    for (int i=0;i<6;i++)
    {
        if(gaitPhase[i]==NormalGait::GaitPhase::Swing)
        {
            swingLegTg(robot,param,i);
        }
        else if(gaitPhase[i]==NormalGait::GaitPhase::Stance)
        {
            stanceLegTg(robot,param,i);
        }
    }

    /*
    //test IMU//
    double eul[3];
    param.imu_data->toEulBody2Ground(eul,"213");
    rt_printf("eul:%.4f,%.4f,%.4f\n",eul[0],eul[1],eul[2]);*/

    robot.GetPin(OPP.outputPin);
    robot.GetPee(OPP.outputPee,robot.body());
    fastWalkPipe.sendToNrt(OPP);

    if(walkState==NormalGait::WalkState::Stop && param.count%totalCount==(totalCount-1))
    {
        return 0;
    }
    else
    {
        return 1;
    }
}
