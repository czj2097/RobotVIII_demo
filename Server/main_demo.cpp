#include <aris.h>
#include "Robot_Type_I.h"
#include "plan.h"
#include "NonRTOptimalGCSRotate.h"
#include "NonRTOptimalGCS360.h"

void BoPdataProcessing();
void twoLevelPlanner();
void proveModel();
void analyseRobot360Ability();
void analyseRobotRotAbility();

int main()
{
    //proveModel();
    //twoLevelPlanner();
    analyseRobot360Ability();


    printf("Finished!\n");
	char aaa;
	std::cin >> aaa;
	return 0;
}

void analyseRobot360Ability()
{
    double v_mtx[361][12];
    double t_mtx[361][12];
    double dinB=0.2;
    double dinG;
    double h=0.05;
    double a=0;
    double b=0;
    double dutyCycle=0.50001;
    double aLmt=6.4;
    double vLmt=1.5;
    double initPee[18]
    { -0.3,   -0.58,  -0.52,
      -0.6,   -0.58,   0,
      -0.3,   -0.58,   0.52,
       0.3,   -0.58,  -0.52,
       0.6,   -0.58,   0,
       0.3,   -0.58,   0.52 };
    double out_TipPos[3000][18];
    double bodyConstVel;
    double out_opt_period;
    double out_nml_period;
    TimeOptimal::NonRTOptimalGCS360 planner;

    for(int j=0;j<6;j++)
    {
        for(int i=0;i<361;i++)
        {
            dutyCycle=0.5+0.05*j;
            dinG=dinB/dutyCycle;
            a=2*PI*i/360;

            planner.GetTimeOptimalGait(dinG,h,a,b,dutyCycle,aLmt,vLmt,initPee,*out_TipPos,bodyConstVel,out_opt_period);
            planner.GetNormalGait(out_nml_period);
            t_mtx[i][2*j+0]=out_opt_period;
            t_mtx[i][2*j+1]=out_nml_period;
            v_mtx[i][2*j+0]=dinG/out_opt_period;
            v_mtx[i][2*j+1]=dinG/out_nml_period;

            printf("\n\n\nOptimal gait of i=%d Finished!\n\n\n",i);
        }
    }
    aris::dynamic::dlmwrite("./twotime.txt",*t_mtx,361,12);
    aris::dynamic::dlmwrite("./twovel.txt",*v_mtx,361,12);


}

void analyseRobotRotAbility()
{
    double v_mtx[51][2];
    double t_mtx[51][2];
    double d=0;
    double h=0.05;
    double a=0;
    double binB=0.2;
    double binG;
    double dutyCycle=0.5;
    double aLmt=6.4;
    double vLmt=1.5;
    double initPee[18]
    { -0.3,   -0.58,  -0.52,
      -0.6,   -0.58,   0,
      -0.3,   -0.58,   0.52,
       0.3,   -0.58,  -0.52,
       0.6,   -0.58,   0,
       0.3,   -0.58,   0.52 };
    double out_TipPos[3000][18];
    double bodyConstVel;
    double out_opt_period;
    double out_nml_period;
    TimeOptimal::NonRTOptimalGCSRotate planner;
    for(int j=0;j<51;j++)
    {
        dutyCycle=0.5+0.005*j;
        binG=binB/dutyCycle;
        planner.GetTimeOptimalGait(d,h,a,binG,dutyCycle,aLmt,vLmt,initPee,*out_TipPos,bodyConstVel,out_opt_period);
        planner.GetNormalGait(out_nml_period);
        t_mtx[j][0]=out_opt_period;
        t_mtx[j][1]=out_nml_period;
        v_mtx[j][0]=binG/out_opt_period;
        v_mtx[j][1]=binG/out_nml_period;

        printf("\n\n\nOptimal gait of j=%d Finished!\n\n\n",j);
    }
    aris::dynamic::dlmwrite("./rottime.txt",*t_mtx,51,2);
    aris::dynamic::dlmwrite("./rotvel.txt",*v_mtx,51,2);
}

void printPm(const double *pm)
{
    printf("%f,%f,%f,%f\n%f,%f,%f,%f\n%f,%f,%f,%f\n%f,%f,%f,%f\n",
           pm[0],pm[1],pm[2],pm[3],
           pm[4],pm[5],pm[6],pm[7],
           pm[8],pm[9],pm[10],pm[11],
           pm[12],pm[13],pm[14],pm[15]);
}

void proveModel()
{
    Robots::RobotTypeI robot;
    robot.loadXml("/home/hex/Desktop/mygit/RobotVIII_demo/resource/RobotEDU2_re.xml");

    double initPeb[6] {0};
    double initPee[18]
    { -0.3,   -0.58,  -0.52,
      -0.6,   -0.58,   0,
      -0.3,   -0.58,   0.52,
       0.3,   -0.58,  -0.52,
       0.6,   -0.58,   0,
       0.3,   -0.58,   0.52 };
    double initPin[18];
    robot.SetPeb(initPeb,"313");
    robot.SetPee(initPee);
    robot.GetPin(initPin);

    double L2 {0.05};
    double L1 {0.05};
    int totalCount {1000};

    double output[1001][24];
    double Pee[3];
    double Vee[3];
    double Aee[3];
    double VeeL[3];
    double AeeL[3];
    double inv_prtPm[16];
    double Pin[3];
    double VinL[3];
    double AinL[3];
    double Vin[3];
    double Ain[3];
    double dJi_x_L[9];
    double dJi_y_L[9];
    double dJi_z_L[9];
    double dJi_L[9];
    double Ji_L[9];
    double dJi_x_B[9];
    double dJi_y_B[9];
    double dJi_z_B[9];
    double dJi_B[9];
    double Ji_B[9];
    memcpy(Pee,initPee+3,3*sizeof(double));


    printPm(*robot.pLegs[1]->base().prtPm());
    aris::dynamic::s_inv_pm(*robot.pLegs[1]->base().prtPm(),inv_prtPm);

    printf("Peb:%f,%f,%f,%f,%f,%f\n",initPeb[0],initPeb[1],initPeb[2],initPeb[3],initPeb[4],initPeb[5]);
    printf("Pee:%f,%f,%f\n",Pee[0],Pee[1],Pee[2]);

    for(int i=0;i<1001;i++)
    {
        Pin[0]=initPin[3]-0.5*L1*(1-cos(PI*i/totalCount));
        Pin[1]=initPin[4]-0.5*L2*(1-cos(2*PI*i/totalCount));
        Pin[2]=initPin[5]+0.5*L1*(1-cos(PI*i/totalCount));

        Pee[0]=initPee[3]-0.5*L1*(1-cos(PI*i/totalCount));
        Pee[1]=initPee[4]+0.5*L2*(1-cos(2*PI*i/totalCount));
        Pee[2]=initPee[5]+0.5*L1*(1-cos(PI*i/totalCount));

        Vee[0]=-0.5*L1*sin(PI*i/totalCount)*PI/totalCount*1000;
        Vee[1]=0.5*L2*sin(2*PI*i/totalCount)*2*PI/totalCount*1000;
        Vee[2]=0.5*L1*sin(PI*i/totalCount)*PI/totalCount*1000;

        Aee[0]=-0.5*L1*cos(PI*i/totalCount)*PI/totalCount*1000*PI/totalCount*1000;
        Aee[1]=0.5*L2*cos(2*PI*i/totalCount)*2*PI/totalCount*1000*2*PI/totalCount*1000;
        Aee[2]=0.5*L1*cos(PI*i/totalCount)*PI/totalCount*1000*PI/totalCount*1000;

        //prove kinematics
        robot.pLegs[1]->SetPee(Pee);
        robot.pLegs[1]->GetPee(*output+24*i);
        robot.pLegs[1]->GetPin(*output+24*i+3);

        robot.pLegs[1]->SetPin(Pin);
        robot.pLegs[1]->GetPee(*output+24*i+6);
        robot.pLegs[1]->GetPin(*output+24*i+9);

        //prove influence coefficient in "L"
        aris::dynamic::s_dgemm(3, 1, 3, 1, inv_prtPm, 4, Vee, 1, 0, VeeL, 1);
        aris::dynamic::s_dgemm(3, 1, 3, 1, inv_prtPm, 4, Aee, 1, 0, AeeL, 1);

        robot.pLegs[1]->SetPee(Pee);
        robot.pLegs[1]->GetJvi(Ji_L,robot.pLegs[1]->base());
        robot.pLegs[1]->GetdJacOverPee(dJi_x_L,dJi_y_L,dJi_z_L,"L");
        std::fill_n(dJi_L,9,0);
        aris::dynamic::s_daxpy(9, VeeL[0], dJi_x_L, 1, dJi_L, 1);
        aris::dynamic::s_daxpy(9, VeeL[1], dJi_y_L, 1, dJi_L, 1);
        aris::dynamic::s_daxpy(9, VeeL[2], dJi_z_L, 1, dJi_L, 1);

        double tmp1[3];
        double tmp2[3];
        aris::dynamic::s_dgemm(3, 1, 3, 1, dJi_L, 3, VeeL, 1, 0, tmp1, 1);
        aris::dynamic::s_dgemm(3, 1, 3, 1, Ji_L, 3, AeeL, 1, 0, tmp2, 1);
        for(int j=0;j<3;j++)
        {
            AinL[j]=tmp1[j]+tmp2[j];
        }

        aris::dynamic::s_dgemm(3, 1, 3, 1, Ji_L, 3, VeeL, 1, 0, VinL, 1);
        memcpy(*output+24*i+12,VinL,3*sizeof(double));
        memcpy(*output+24*i+15,AinL,3*sizeof(double));

        //prove influence coefficient in "B"
        robot.pLegs[1]->SetPee(Pee);
        robot.pLegs[1]->GetJvi(Ji_B,robot.body());
        robot.pLegs[1]->GetdJacOverPee(dJi_x_B,dJi_y_B,dJi_z_B,"B");
        std::fill_n(dJi_B,9,0);
        aris::dynamic::s_daxpy(9, Vee[0], dJi_x_B, 1, dJi_B, 1);
        aris::dynamic::s_daxpy(9, Vee[1], dJi_y_B, 1, dJi_B, 1);
        aris::dynamic::s_daxpy(9, Vee[2], dJi_z_B, 1, dJi_B, 1);

        double tmp3[3];
        double tmp4[3];
        aris::dynamic::s_dgemm(3, 1, 3, 1, dJi_B, 3, Vee, 1, 0, tmp3, 1);
        aris::dynamic::s_dgemm(3, 1, 3, 1, Ji_B, 3, Aee, 1, 0, tmp4, 1);
        for(int j=0;j<3;j++)
        {
            Ain[j]=tmp3[j]+tmp4[j];
        }

        aris::dynamic::s_dgemm(3, 1, 3, 1, Ji_B, 3, Vee, 1, 0, Vin, 1);
        memcpy(*output+24*i+18,Vin,3*sizeof(double));
        memcpy(*output+24*i+21,Ain,3*sizeof(double));
    }

    aris::dynamic::dlmwrite("./proveModel.txt",*output,1001,24);
}

void twoLevelPlanner()
{
    Robots::RobotTypeI robot;
    robot.loadXml("/home/hex/Desktop/mygit/RobotVIII_demo/resource/RobotEDU2_re.xml");
    double initPee[18]
    { -0.3,   -0.58,  -0.52,
      -0.6,   -0.58,   0,
      -0.3,   -0.58,   0.52,
       0.3,   -0.58,  -0.52,
       0.6,   -0.58,   0,
       0.3,   -0.58,   0.52 };
    double initPeb[6] {0};
    double L {0.3};
    double H {0.05};
    int totalCount {1000};
    double s1;
    double s2;

    double beginPeb[6];
    double beginPee1[18];
    double beginPee2[18];
    double output[4001][78];
    double Peb[6];
    double Pee1[18];
    double Pee2[18];
    memcpy(Peb,initPeb,6*sizeof(double));
    memcpy(Pee1,initPee,18*sizeof(double));
    memcpy(Pee2,initPee,18*sizeof(double));

    for(int i=0;i<4001;i++)
    {
        s1=PI/2*(1-cos(PI*(i%totalCount)/totalCount));
        if((i%totalCount)<totalCount/2)
        {
            s2=2*PI/totalCount/totalCount*(i%totalCount)*(i%totalCount);
        }
        else
        {
            s2=-2*PI/totalCount/totalCount*(i%totalCount)*(i%totalCount)+4*PI/totalCount*(i%totalCount)-PI;
        }

        if(i%1000==0)
        {
            memcpy(beginPeb,Peb,6*sizeof(double));
            memcpy(beginPee1,Pee1,18*sizeof(double));
            memcpy(beginPee2,Pee2,18*sizeof(double));
        }
        if(i<1000)
        {
            memcpy(Peb,beginPeb,6*sizeof(double));
            Peb[2]=beginPeb[2]-L/4*(i%totalCount)/totalCount*(i%totalCount)/totalCount;

            for(int j=0;j<6;j++)
            {
                if(j%2==0)
                {
                    Pee1[3*j]  =beginPee1[3*j];
                    Pee1[3*j+1]=beginPee1[3*j+1]+H*sin(s1);
                    Pee1[3*j+2]=beginPee1[3*j+2]-L/4+L/4*cos(s1);

                    Pee2[3*j]  =beginPee2[3*j];
                    Pee2[3*j+1]=beginPee2[3*j+1]+H*sin(s2);
                    Pee2[3*j+2]=beginPee2[3*j+2]-L/4+L/4*cos(s2);
                }
                else
                {
                    Pee1[3*j]  =beginPee1[3*j];
                    Pee1[3*j+1]=beginPee1[3*j+1];
                    Pee1[3*j+2]=beginPee1[3*j+2];

                    Pee2[3*j]  =beginPee2[3*j];
                    Pee2[3*j+1]=beginPee2[3*j+1];
                    Pee2[3*j+2]=beginPee2[3*j+2];
                }
            }
        }
        else if(i>=3000)
        {
            memcpy(Peb,beginPeb,6*sizeof(double));
            Peb[2]=beginPeb[2]-L/2*(i%totalCount)/totalCount+L/4*(i%totalCount)/totalCount*(i%totalCount)/totalCount;

            for(int j=0;j<6;j++)
            {
                if(j%2==1)
                {
                    Pee1[3*j]  =beginPee1[3*j];
                    Pee1[3*j+1]=beginPee1[3*j+1]+H*sin(s1);
                    Pee1[3*j+2]=beginPee1[3*j+2]-L/4+L/4*cos(s1);

                    Pee2[3*j]  =beginPee2[3*j];
                    Pee2[3*j+1]=beginPee2[3*j+1]+H*sin(s2);
                    Pee2[3*j+2]=beginPee2[3*j+2]-L/4+L/4*cos(s2);
                }
                else
                {
                    Pee1[3*j]  =beginPee1[3*j];
                    Pee1[3*j+1]=beginPee1[3*j+1];
                    Pee1[3*j+2]=beginPee1[3*j+2];

                    Pee2[3*j]  =beginPee2[3*j];
                    Pee2[3*j+1]=beginPee2[3*j+1];
                    Pee2[3*j+2]=beginPee2[3*j+2];
                }
            }
        }
        else if(i%(2*totalCount)<1000)
        {
            memcpy(Peb,beginPeb,6*sizeof(double));
            Peb[2]=beginPeb[2]-L/2*(i%totalCount)/totalCount;

            for(int j=0;j<6;j++)
            {
                if(j%2==0)
                {
                    Pee1[3*j]  =beginPee1[3*j];
                    Pee1[3*j+1]=beginPee1[3*j+1]+H*sin(s1);
                    Pee1[3*j+2]=beginPee1[3*j+2]-L/2+L/2*cos(s1);

                    Pee2[3*j]  =beginPee2[3*j];
                    Pee2[3*j+1]=beginPee2[3*j+1]+H*sin(s2);
                    Pee2[3*j+2]=beginPee2[3*j+2]-L/2+L/2*cos(s2);
                }
                else
                {
                    Pee1[3*j]  =beginPee1[3*j];
                    Pee1[3*j+1]=beginPee1[3*j+1];
                    Pee1[3*j+2]=beginPee1[3*j+2];

                    Pee2[3*j]  =beginPee2[3*j];
                    Pee2[3*j+1]=beginPee2[3*j+1];
                    Pee2[3*j+2]=beginPee2[3*j+2];
                }
            }
        }
        else
        {
            memcpy(Peb,beginPeb,6*sizeof(double));
            Peb[2]=beginPeb[2]-L/2*(i%totalCount)/totalCount;

            for(int j=0;j<6;j++)
            {
                if(j%2==1)
                {
                    Pee1[3*j]  =beginPee1[3*j];
                    Pee1[3*j+1]=beginPee1[3*j+1]+H*sin(s1);
                    Pee1[3*j+2]=beginPee1[3*j+2]-L/2+L/2*cos(s1);

                    Pee2[3*j]  =beginPee2[3*j];
                    Pee2[3*j+1]=beginPee2[3*j+1]+H*sin(s2);
                    Pee2[3*j+2]=beginPee2[3*j+2]-L/2+L/2*cos(s2);
                }
                else
                {
                    Pee1[3*j]  =beginPee1[3*j];
                    Pee1[3*j+1]=beginPee1[3*j+1];
                    Pee1[3*j+2]=beginPee1[3*j+2];

                    Pee2[3*j]  =beginPee2[3*j];
                    Pee2[3*j+1]=beginPee2[3*j+1];
                    Pee2[3*j+2]=beginPee2[3*j+2];
                }
            }
        }

        *(*output+78*i+18*4)=s1;
        *(*output+78*i+18*4+1)=s2;
        *(*output+78*i+18*4+2)=Peb[2];
        robot.SetPeb(Peb,"313");
        robot.SetPee(Pee1);
        robot.GetPee(*output+78*i);
        robot.GetPin(*output+78*i+18);

        robot.SetPee(Pee2);
        robot.GetPee(*output+78*i+18*2);
        robot.GetPin(*output+78*i+18*3);
    }
    aris::dynamic::dlmwrite("./twoleveltripod.txt",*output,4001,78);
}

void BoPdataProcessing()
{
    Robots::RobotTypeI robot;
    robot.loadXml("/home/hex/Desktop/mygit/RobotVIII_demo/resource/RobotEDU2_re.xml");

    double input[16739][26];
    double output[16739][24];
    aris::dynamic::dlmread("./8.txt",*input);
//    for(int i=0;i<16739*2;i++)
//    {
//        for(int j=0;j<6;j++)
//        {
//            printf("%f\t",input[i][j+3]);
//        }
//        printf("\n");
//    }
    double delta[2];
    double bodyPe[6];
    double beginPe[6];
    double beginPm[16];
    for(int i=0;i<16739;i++)
    {
        static aris::dynamic::FloatMarker beginMak{robot.ground()};

        for(int j=0;j<6;j++)
        {
            bodyPe[j]=input[i][j];
            beginPe[j]=input[i][j];
        }

        if(i%1000==0)
        {
            beginPe[4]=0;
            beginPe[5]=0;
            aris::dynamic::s_pe2pm(beginPe,beginPm,"213");
            beginMak.setPrtPm(beginPm);
            beginMak.update();

            delta[0]=input[i][4]-input[i][6];
            delta[1]=input[i][5]-input[i][7];
        }
        bodyPe[4]-=delta[0];
        bodyPe[5]-=delta[1];

        robot.SetPeb(bodyPe,robot.ground(),"213");
        robot.SetPin(*input+26*i+8);
        if(i<100)
        {
        for(int j=0;j<18;j++)
        {
            printf("%.4f, ",*(*input+26*i+8+j));
        }
        printf("\n");
        }


        robot.GetPeb(*output+24*i,robot.ground(),"213");
        robot.GetPee(*output+24*i+6,beginMak);
    }
    aris::dynamic::dlmwrite("./8_out.txt",*output,16739,24);
}
