#include "RealTimeOptimal.h"

/*find gait with maxVel in joint space by iteration*/
void screwInterpolationTraj()
{
    Robots::RobotTypeI rbt;
    rbt.loadXml("/home/hex/Desktop/mygit/RobotVIII_demo/resource/RobotEDU2.xml");

    timeval tpstart,tpend;
    float tused;
    gettimeofday(&tpstart,NULL);

    double totalTime{1.5};

    double ratio{0};//control point, from middle to edge [0,1]

    double stepH=0.05;
    double stepD=0.5;
    double initPeb[6]{0,0,0,0,0,0};
    double initVeb[6]{0,0,0,0,0,0};
//    double initPee[18]{ -0.3, -0.85, -0.65,
//                       -0.45, -0.85, 0,
//                        -0.3, -0.85, 0.65,
//                         0.3, -0.85, -0.65,
//                        0.45, -0.85, 0,
//                         0.3, -0.85, 0.65 };

    double initPee[18] { -0.30, -0.58, -0.52,
                         -0.60, -0.58,  0,
                         -0.30, -0.58,  0.52,
                          0.30, -0.58, -0.52,
                          0.60, -0.58,  0,
                          0.30, -0.58,  0.52 };

    double initVee[18]{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    double pEE[18];
    double pEE_B[18];

    double aEE[18];
    double vLmt{1.0};
    double aLmt{3.2};

    int keyPointNum{3};
    double keyPee_B[keyPointNum][18];
    double keyPin[keyPointNum][18];
    double keyVin[keyPointNum][18];

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
    double totalTmax{0},totalT[18];
    double t1[18],t2[18],t3[18],t4[18],t5[18],t6[18];
    double s1[18],s2[18],s3[18],s4[18],s5[18],s6[18];

    while(e>=0.0001)
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
        s1[i]=keyVin[0][i]*t1[i]-0.5*aLmt*t1[i]*t1[i];//shift in phase 1 (vector with direction)
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

    aris::dynamic::dlmwrite("./vIn.txt",*vIn,3000,18);
    aris::dynamic::dlmwrite("./pIn.txt",*pIn,3000,18);
    aris::dynamic::dlmwrite("./pIn_adjust.txt",*pIn_adjust,totalCount,18);
    aris::dynamic::dlmwrite("./pEE_adjust.txt",*pEE_adjust,totalCount,18);

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
