#include "TimeOptimalGait.h"
#include <Robot_Type_I.h>

TimeOptimalGait::TimeOptimalGait(){}
TimeOptimalGait::~TimeOptimalGait(){}

void TimeOptimalGait::SetInitPos(double *pEB, double *pEE)
{
    memcpy(initPeb,pEB,6*sizeof(double));
    memcpy(initPee,pEE,18*sizeof(double));
}

void TimeOptimalGait::SetBodyTraj(double *traj)
{
    memcpy(bodyTraj,traj,(totalCount+1)*sizeof(double));
}

void TimeOptimalGait::SetSwingLegTraj(int legID, double *traj)
{
    for(int i=0;i<totalCount+1;i++)
    {
        swingLegTraj[i][legID]=traj[i];
    }
}



void TimeOptimalGait1by1()
{
    timeval tpstart,tpend;
    float tused;
    gettimeofday(&tpstart,NULL);


    Robots::RobotTypeI rbt;
    rbt.loadXml("/home/hex/Desktop/mygit/RobotVIII_demo/resource/Robot_VIII.xml");

    double stepH {0.1};
    double stepD {0.8};
    double vLmt {0.9};
    double aLmt {3.2};
    const int sTotalCount {1801};
    const double delta_s {PI/(sTotalCount-1)};

    double initPeb[6] {0};
    double initPee[18] {  -0.3, -0.85, -0.65,
                         -0.45, -0.85,     0,
                          -0.3, -0.85,  0.65,
                           0.3, -0.85, -0.65,
                          0.45, -0.85,     0,
                           0.3, -0.85,  0.65 };//024 swing, 135 stance
    double pEB[6] {0};
    double pEE[18] {0};

    double s_t {0};
    double b_st[sTotalCount] {0};
    double db_st[sTotalCount] {0};
    double ddb_st[sTotalCount] {0};
    double pb_sw[sTotalCount] {0};
    double vb_sw[sTotalCount] {0};
    double ab_sw[sTotalCount] {0};

    double s_w {0};
    double f_sw[3] {0};
    double df_sw[3] {0};
    double ddf_sw[3] {0};

    double Jvi[9] {0};
    double dJvi_x[9] {0};
    double dJvi_y[9] {0};
    double dJvi_z[9] {0};
    double dJvi[9] {0};//used in stanceLeg
    double dJvi_dot_f[9] {0};//used in swingLeg
    double dJvi_dot_vb[9] {0};//used in swingLeg

    double param_dsds[18] {0};
    double param_dsds1[18] {0};
    double param_dsds2[18] {0};
    double param_dds[18] {0};
    double abs_param_dds[9] {0};//for vLmt of stanceLeg
    double param_ds[18] {0};//for aLmt of swingLeg
    double param_ds1[18] {0};
    double param_ds2[18] {0};
    double param_const[18] {0};//for aLmt of swingLeg
    double param_const1[18] {0};
    double param_const2[18] {0};

    double param_a2[18] {0};//coefficient of ds*ds, param_dsds divided by param_dds
    double param_a1[18] {0};//coefficient of ds, param_ds divided by param_dds
    double param_a0L[18] {0};
    double param_a0H[18] {0};

    double output_PeeB[sTotalCount][18] {0};
    double output_dsds[sTotalCount][18] {0};
    double output_dsds1[sTotalCount][18] {0};
    double output_dsds2[sTotalCount][18] {0};
    double output_ds[sTotalCount][18] {0};
    double output_ds1[sTotalCount][18] {0};
    double output_ds2[sTotalCount][18] {0};
    double output_const[sTotalCount][18] {0};
    double output_const1[sTotalCount][18] {0};
    double output_const2[sTotalCount][18] {0};
    double output_dds[sTotalCount][18] {0};
    double output_a2[sTotalCount][18] {0};
    double output_a1[sTotalCount][18] {0};
    double output_a0L[sTotalCount][18] {0};
    double output_a0H[sTotalCount][18] {0};

    double ds_upBound_aLmt[sTotalCount][4] {0};
    double ds_lowBound_aLmt[sTotalCount][4] {0};
    double ds_upBound_vLmt[sTotalCount][4] {0};
    double ds_upBound[sTotalCount][4] {0};
    double ds_lowBound[sTotalCount][4] {0};
    double dds_lowBound[sTotalCount][4] {0};
    double dds_upBound[sTotalCount][4] {0};

    double slopedsBound[sTotalCount][4] {0};
    double slopeDelta[sTotalCount][4] {0};
    int paramdds0Point[sTotalCount][18] {0};
    int tangentPoint[sTotalCount][4] {0};
    int switchPoint[sTotalCount][4] {0};
    int paramdds0Count[18] {0};
    int tangentCount[4] {0};
    int switchCount[4] {0};

    bool stopFlag {false};
    bool accFlag {true};//true acc, false dec
    int dec_start {0};
    int dec_end {0};
    double min_dist[sTotalCount];
    unsigned int cycleCount {0};//used to limit the cycleCount of forward integration
    int stop_back {0};
    int ki_back {0};
    int ki_for {0};

    double ds_forward[sTotalCount][4] {0};
    double ds_backward[sTotalCount][4] {0};
    double dds_forward[sTotalCount][4] {0};
    double dds_backward[sTotalCount][4] {0};
    double timeArray[sTotalCount][4] {0};
    double timeArray_tmp[sTotalCount][4] {0};
    double totalTime[4] {0};
    double maxTime {0};
    int maxTotalCount {0};
    int maxTotalCount_last {1};
    int maxTimeID {0};

    double real_ds[sTotalCount][4] {0};
    double real_dds[sTotalCount][4] {0};
    double real_ddsMax[sTotalCount][4] {0};
    double real_ddsMin[sTotalCount][4] {0};

    /*******************************************************************************/
    /********** StanceLeg : generate the traj & calculate the bound of ds **********/
    for (int i=0;i<sTotalCount;i++)
    {
        s_t=i*PI/(sTotalCount-1);//degree to rad
        b_st[i]=stepD/4-stepD/2*(s_t/PI);
        db_st[i]=-stepD/2/PI;
        ddb_st[i]=0;

        pEB[2]=initPeb[2]+b_st[i];
        rbt.SetPeb(pEB);
        for(int j=0;j<3;j++)//pEE of stance leg
        {
            pEE[6*j+3]=initPee[6*j+3];
            pEE[6*j+4]=initPee[6*j+4];
            pEE[6*j+5]=initPee[6*j+5];
            rbt.pLegs[2*j+1]->SetPee(pEE+6*j+3);
        }

        //calculate param of stance leg, only vel in z direction
        for(int j=0;j<3;j++)
        {
            rbt.pLegs[2*j+1]->GetJvi(Jvi,rbt.body());
            rbt.pLegs[2*j+1]->GetdJacOverPee(dJvi_x,dJvi_y,dJvi_z,"B");
            rbt.pLegs[2*j+1]->GetPee(*output_PeeB+18*i+6*j+3,rbt.body());

            double db_st_tmp[3] {0};
            db_st_tmp[2]=-db_st[i];
            double ddb_st_tmp[3] {0};
            ddb_st_tmp[2]=-ddb_st[i];

            std::fill_n(dJvi,9,0);
            aris::dynamic::s_daxpy(9,db_st_tmp[0],dJvi_x,1,dJvi,1);//for s_t
            aris::dynamic::s_daxpy(9,db_st_tmp[1],dJvi_y,1,dJvi,1);
            aris::dynamic::s_daxpy(9,db_st_tmp[2],dJvi_z,1,dJvi,1);

            std::fill_n(param_dds+6*j+3,3,0);
            aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,db_st_tmp,1,1,param_dds+6*j+3,1);

            std::fill_n(param_dsds1+6*j+3,3,0);
            aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,ddb_st_tmp,1,1,param_dsds1+6*j+3,1);
            std::fill_n(param_dsds2+6*j+3,3,0);
            aris::dynamic::s_dgemm(3,1,3,1,dJvi,3,db_st_tmp,1,1,param_dsds2+6*j+3,1);
            std::fill_n(param_dsds+6*j+3,3,0);
            for (int k=0;k<3;k++)
            {
                param_dsds[6*j+3+k]=param_dsds1[6*j+3+k]+param_dsds2[6*j+3+k];
                abs_param_dds[3*j+k]=fabs(param_dds[6*j+3+k]);

                if(param_dds[6*j+3+k]>0)
                {
                    param_a2[6*j+3+k]=-param_dsds[6*j+3+k]/param_dds[6*j+3+k];
                    param_a0L[6*j+3+k]=(-param_const[6*j+k+3]-aLmt)/param_dds[6*j+3+k];
                    param_a0H[6*j+3+k]=(-param_const[6*j+k+3]+aLmt)/param_dds[6*j+3+k];
                }
                else if(param_dds[6*j+3+k]<0)
                {
                    param_a2[6*j+3+k]=-param_dsds[6*j+3+k]/param_dds[6*j+3+k];
                    param_a0L[6*j+3+k]=(-param_const[6*j+k+3]+aLmt)/param_dds[6*j+3+k];
                    param_a0H[6*j+3+k]=(-param_const[6*j+k+3]-aLmt)/param_dds[6*j+3+k];
                }
                else
                {
                    printf("WARNING!!! param_dds equals zero!!! StanceLeg : i=%d \n",i);
                }
            }
        }

        //calculate bound of ds of 3 stance legs
        int k_st {0};
        bool dsBoundFlag_st {false};
        const int kstCount {15000};
        while (dsBoundFlag_st==false)
        {
            double ds=0.001*k_st;
            double dec[9]{0};
            double acc[9]{0};
            double max_dec {0};
            double min_acc {0};
            for (int j=0;j<3;j++)
            {
                for (int k=0;k<3;k++)
                {
                    if(param_dds[6*j+3+k]>0)
                    {
                        dec[3*j+k]=(-param_dsds[6*j+k+3]*ds*ds-param_const[6*j+k+3]-aLmt)/param_dds[6*j+3+k];
                        acc[3*j+k]=(-param_dsds[6*j+k+3]*ds*ds-param_const[6*j+k+3]+aLmt)/param_dds[6*j+3+k];
                    }
                    else if(param_dds[6*j+3+k]<0)
                    {
                        dec[3*j+k]=(-param_dsds[6*j+k+3]*ds*ds-param_const[6*j+k+3]+aLmt)/param_dds[6*j+3+k];
                        acc[3*j+k]=(-param_dsds[6*j+k+3]*ds*ds-param_const[6*j+k+3]-aLmt)/param_dds[6*j+3+k];
                    }
                    /*
                    dec[3*j+k]=param_a2[6*j+3+k]*ds*ds+param_a0L[6*j+3+k];
                    acc[3*j+k]=param_a2[6*j+3+k]*ds*ds+param_a0H[6*j+3+k];
                    */
                }
            }
            max_dec=*std::max_element(dec,dec+9);
            min_acc=*std::min_element(acc,acc+9);

            k_st++;
            if(k_st==kstCount)
            {
                dsBoundFlag_st=true;
                printf("WARNING!!! kstCount=%d is too small!!!\n",kstCount);
            }
            else
            {
                if(min_acc<max_dec)
                {
                    for(int k_st2=0;k_st2<10000;k_st2++)
                    {
                        ds=ds_upBound_aLmt[i][0]+0.001*0.0001*k_st2;
                        for (int j=0;j<3;j++)
                        {
                            for (int k=0;k<3;k++)
                            {
                                if(param_dds[6*j+3+k]>0)
                                {
                                    dec[3*j+k]=(-param_dsds[6*j+k+3]*ds*ds-param_const[6*j+k+3]-aLmt)/param_dds[6*j+3+k];
                                    acc[3*j+k]=(-param_dsds[6*j+k+3]*ds*ds-param_const[6*j+k+3]+aLmt)/param_dds[6*j+3+k];
                                }
                                else if(param_dds[6*j+3+k]<0)
                                {
                                    dec[3*j+k]=(-param_dsds[6*j+k+3]*ds*ds-param_const[6*j+k+3]+aLmt)/param_dds[6*j+3+k];
                                    acc[3*j+k]=(-param_dsds[6*j+k+3]*ds*ds-param_const[6*j+k+3]-aLmt)/param_dds[6*j+3+k];
                                }
                                /*
                                dec[3*j+k]=param_a2[6*j+3+k]*ds*ds+param_a0L[6*j+3+k];
                                acc[3*j+k]=param_a2[6*j+3+k]*ds*ds+param_a0H[6*j+3+k];
                                */
                            }
                        }
                        max_dec=*std::max_element(dec,dec+9);
                        min_acc=*std::min_element(acc,acc+9);

                        if(min_acc<max_dec)
                        {
                            //printf("stance ds_upBound_aLmt=%.6f\t",ds_upBound_aLmt[i][0]);
                            dsBoundFlag_st=true;
                            break;
                        }
                        else
                        {
                            dds_lowBound[i][0]=max_dec;
                            dds_upBound[i][0]=min_acc;
                            ds_upBound_aLmt[i][0]=ds;
                        }
                    }
                }
                else
                {
                    ds_upBound_aLmt[i][0]=ds;
                }
            }
        }

        ds_upBound_vLmt[i][0]=vLmt/(*std::max_element(abs_param_dds,abs_param_dds+9));
        ds_upBound[i][0]=std::min(ds_upBound_aLmt[i][0],ds_upBound_vLmt[i][0]);
        ds_lowBound[i][0]=ds_lowBound_aLmt[i][0];

        for (int j=0;j<3;j++)
        {
            memcpy(*output_dsds1+18*i+6*j+3, param_dsds1+6*j+3, 3*sizeof(double));
            memcpy(*output_dsds2+18*i+6*j+3, param_dsds2+6*j+3, 3*sizeof(double));
            memcpy( *output_dsds+18*i+6*j+3,  param_dsds+6*j+3, 3*sizeof(double));
            memcpy(  *output_dds+18*i+6*j+3,   param_dds+6*j+3, 3*sizeof(double));
            memcpy(   *output_a2+18*i+6*j+3,    param_a2+6*j+3, 3*sizeof(double));
            memcpy(  *output_a0L+18*i+6*j+3,   param_a0L+6*j+3, 3*sizeof(double));
            memcpy(  *output_a0H+18*i+6*j+3,   param_a0H+6*j+3, 3*sizeof(double));
        }
    }

    /********** StanceLeg : find switch point - param_dds=0 point & tangent point **********/
    //initialize
    for(int i=0;i<sTotalCount;i++)
    {
        tangentPoint[i][0]=-1;
        switchPoint[i][0]=-1;
        for(int j=0;j<3;j++)
        {
            for(int k=0;k<3;k++)
            {
                paramdds0Point[i][6*j+k+3]=-1;
            }
        }
    }

    //calculate the slope of ds_upBound
    slopedsBound[0][0]=(ds_upBound[1][0]-ds_upBound[0][0])/delta_s;
    for(int i=1;i<sTotalCount-1;i++)
    {
        slopedsBound[i][0]=(ds_upBound[i+1][0]-ds_upBound[i-1][0])/2/delta_s;
    }
    slopedsBound[sTotalCount-1][0]=(ds_upBound[sTotalCount-1][0]-ds_upBound[sTotalCount-2][0])/delta_s;

    for(int i=0;i<sTotalCount;i++)
    {
        slopeDelta[i][0]=slopedsBound[i][0]-(dds_upBound[i][0]+dds_lowBound[i][0])/2;
    }

    for(int i=0;i<sTotalCount-1;i++)
    {
        for(int j=0;j<3;j++)
        {
            for(int k=0;k<3;k++)
            {
                if(output_dds[i+1][6*j+k+3]*output_dds[i][6*j+k+3]<0 || output_dds[i][6*j+k+3]==0)
                {
                    paramdds0Point[paramdds0Count[6*j+k+3]][6*j+k+3]=i;
                    paramdds0Count[6*j+k+3]++;
                }
            }
        }

        if(slopeDelta[i+1][0]*slopeDelta[i][0]<0 || slopeDelta[i][0]==0)
        {
            tangentPoint[tangentCount[0]][0]=i;
            tangentCount[0]++;
        }
    }

    //merge tangentPoint & paramdds0Point into switchPoint
    for(int i=0;i<tangentCount[0];i++)
    {
        switchPoint[i][0]=tangentPoint[i][0];
    }
    switchCount[0]=tangentCount[0];
    for(int j=0;j<3;j++)
    {
        for(int k=0;k<3;k++)
        {
            for(int i=0;i<paramdds0Count[6*j+k+3];i++)
            {
                switchPoint[i+switchCount[0]][0]=paramdds0Point[i][6*j+k+3];
            }
            switchCount[0]+=paramdds0Count[6*j+k+3];
        }
    }

    //filtering the same point & sorting by the value
    for(int i=0;i<switchCount[0];i++)
    {
        for(int j=i+1;j<switchCount[0];j++)
        {
            if(switchPoint[j][0]<switchPoint[i][0])
            {
                int tmp;
                tmp=switchPoint[i][0];
                switchPoint[i][0]=switchPoint[j][0];
                switchPoint[j][0]=tmp;
            }
            else if(switchPoint[j][0]==switchPoint[i][0])
            {
                switchPoint[j][0]=switchPoint[switchCount[0]-1][0];
                switchPoint[switchCount[0]-1][0]=-1;
                j--;
                switchCount[0]--;
            }
        }
    }

    printf("\nStanceLeg Switch Point:");
    for(int i=0;i<switchCount[0]+1;i++)
    {
        printf("%d,",switchPoint[i][0]);
    }
    printf("\n");


    /********** StanceLeg : numerical integration to calculate ds **********/
    for(int i=0;i<sTotalCount;i++)
    {
        real_ds[i][0]=ds_upBound[i][0];
        real_dds[i][0]=dds_upBound[i][0];
    }
    for(int m=0;m<switchCount[0];m++)
    {
        int k_st {switchPoint[m][0]};
        stopFlag=false;
        ds_backward[k_st][0]=ds_upBound[k_st][0];
        while(stopFlag==false)
        {
            double dec[9] {0};
            for (int j=0;j<3;j++)
            {
                for (int k=0;k<3;k++)
                {
                    if(output_dds[k_st][6*j+3+k]>0)
                    {
                        dec[3*j+k]=(-output_dsds[k_st][6*j+k+3]*ds_backward[k_st][0]*ds_backward[k_st][0]-output_const[k_st][6*j+k+3]-aLmt)/output_dds[k_st][6*j+3+k];
                    }
                    else if(output_dds[k_st][6*j+3+k]<0)
                    {
                        dec[3*j+k]=(-output_dsds[k_st][6*j+k+3]*ds_backward[k_st][0]*ds_backward[k_st][0]-output_const[k_st][6*j+k+3]+aLmt)/output_dds[k_st][6*j+3+k];
                    }
                    //dec[3*j+k]=output_a2[k_st][6*j+3+k]*ds_backward[k_st][0]*ds_backward[k_st][0]+output_a0L[k_st][6*j+3+k];
                }
            }
            dds_backward[k_st][0]=*std::max_element(dec,dec+9);
            ds_backward[k_st-1][0]=sqrt(ds_backward[k_st][0]*ds_backward[k_st][0]-2*dds_backward[k_st][0]*delta_s);

            if(ds_backward[k_st-1][0]>ds_upBound[k_st-1][0])
            {
                stopFlag=true;
                printf("StanceLeg backward touching upBound at %d, quit switchPoint %d\n",k_st-1,switchPoint[m][0]);
            }
            else if(k_st==1)
            {
                for (int j=0;j<3;j++)
                {
                    for (int k=0;k<3;k++)
                    {
                        if(output_dds[k_st][6*j+3+k]>0)
                        {
                            dec[3*j+k]=(-output_dsds[k_st][6*j+k+3]*ds_backward[k_st][0]*ds_backward[k_st][0]-output_const[k_st][6*j+k+3]-aLmt)/output_dds[k_st][6*j+3+k];
                        }
                        else if(output_dds[k_st][6*j+3+k]<0)
                        {
                            dec[3*j+k]=(-output_dsds[k_st][6*j+k+3]*ds_backward[k_st][0]*ds_backward[k_st][0]-output_const[k_st][6*j+k+3]+aLmt)/output_dds[k_st][6*j+3+k];
                        }
                        //dec[3*j+k]=output_a2[k_st-1][6*j+3+k]*ds_backward[k_st-1][0]*ds_backward[k_st-1][0]+output_a0L[k_st-1][6*j+3+k];
                    }
                }
                dds_backward[k_st-1][0]=*std::max_element(dec,dec+9);
                for(int i=k_st-1;i<switchPoint[m][0]+1;i++)
                {
                    real_ds[i][0]=ds_backward[i][0];
                    real_dds[i][0]=dds_backward[i][0];
                }
                stopFlag=true;
                printf("StanceLeg backward touching 0, from switchPoint %d\n",switchPoint[m][0]);
            }
            else if(ds_backward[k_st-1][0]>=real_ds[k_st-1][0])
            {
                for (int j=0;j<3;j++)
                {
                    for (int k=0;k<3;k++)
                    {
                        if(output_dds[k_st][6*j+3+k]>0)
                        {
                            dec[3*j+k]=(-output_dsds[k_st][6*j+k+3]*ds_backward[k_st][0]*ds_backward[k_st][0]-output_const[k_st][6*j+k+3]-aLmt)/output_dds[k_st][6*j+3+k];
                        }
                        else if(output_dds[k_st][6*j+3+k]<0)
                        {
                            dec[3*j+k]=(-output_dsds[k_st][6*j+k+3]*ds_backward[k_st][0]*ds_backward[k_st][0]-output_const[k_st][6*j+k+3]+aLmt)/output_dds[k_st][6*j+3+k];
                        }
                        //dec[3*j+k]=output_a2[k_st-1][6*j+3+k]*ds_backward[k_st-1][0]*ds_backward[k_st-1][0]+output_a0L[k_st-1][6*j+3+k];
                    }
                }
                dds_backward[k_st-1][0]=*std::max_element(dec,dec+9);
                for(int i=k_st-1;i<switchPoint[m][0]+1;i++)
                {
                    real_ds[i][0]=ds_backward[i][0];
                    real_dds[i][0]=dds_backward[i][0];
                }
                stopFlag=true;
                printf("StanceLeg backward touching last curve at %d, from switchPoint %d\n",k_st-1,switchPoint[m][0]);
            }
            else
            {
                k_st--;
            }
        }
        if(ds_backward[k_st-1][0]>ds_upBound[k_st-1][0])
        {
            continue;
        }

        bool isEqual0Point {false};
        for(int j=0;j<3;j++)
        {
            for(int k=0;k<3;k++)
            {
                for(int i=0;i<paramdds0Count[6*j+k+3];i++)
                {
                    if(output_dds[switchPoint[i][0]][6*j+k+3]==0 || slopeDelta[switchPoint[i][0]][0]==0)
                    {
                        isEqual0Point=true;
                    }
                }
            }
        }
        if(isEqual0Point==false)
        {
            k_st=switchPoint[m][0]+1;
        }
        else
        {
            k_st=switchPoint[m][0];
        }
        stopFlag=false;
        ds_forward[k_st][0]=ds_upBound[k_st][0];
        while(stopFlag==false)
        {
            double acc[9] {0};
            for (int j=0;j<3;j++)
            {
                for (int k=0;k<3;k++)
                {
                    if(output_dds[k_st][6*j+3+k]>0)
                    {
                        acc[3*j+k]=(-output_dsds[k_st][6*j+k+3]*ds_forward[k_st][0]*ds_forward[k_st][0]-output_const[k_st][6*j+k+3]+aLmt)/output_dds[k_st][6*j+3+k];
                    }
                    else if(output_dds[k_st][6*j+3+k]<0)
                    {
                        acc[3*j+k]=(-output_dsds[k_st][6*j+k+3]*ds_forward[k_st][0]*ds_forward[k_st][0]-output_const[k_st][6*j+k+3]-aLmt)/output_dds[k_st][6*j+3+k];
                    }
                    //acc[3*j+k]=output_a2[k_st][6*j+3+k]*ds_forward[k_st][0]*ds_forward[k_st][0]+output_a0H[k_st][6*j+3+k];
                }
            }
            dds_forward[k_st][0]=*std::min_element(acc,acc+9);
            ds_forward[k_st+1][0]=sqrt(ds_forward[k_st][0]*ds_forward[k_st][0]+2*dds_forward[k_st][0]*delta_s);

            if(ds_forward[k_st+1][0]>ds_upBound[k_st+1][0] || k_st==sTotalCount-2)
            {
                for(int i=switchPoint[m][0]+1;i<k_st+1;i++)
                {
                    real_ds[i][0]=ds_forward[i][0];
                    real_dds[i][0]=dds_forward[i][0];
                }
                stopFlag=true;
                printf("StanceLeg forward touching upBound at %d, from switchPoint %d\n",k_st,switchPoint[m][0]);
            }
            else
            {
                k_st++;
            }
        }
    }

    /*//backward integration
    stopFlag=false;
    ki_back=sTotalCount-1;
    ds_backward[ki_back][0]=ds_upBound[ki_back][0];
    while (stopFlag==false && ki_back>=0)
    {
        double dec[9] {0};
        for (int j=0;j<3;j++)
        {
            for (int k=0;k<3;k++)
            {
                dec[3*j+k]=output_a2[ki_back][6*j+3+k]*ds_backward[ki_back][0]*ds_backward[ki_back][0]
                          +output_a0L[ki_back][6*j+3+k];
            }
        }
        dds_backward[ki_back][0]=*std::max_element(dec,dec+9);
        ds_backward[ki_back-1][0]=sqrt(ds_backward[ki_back][0]*ds_backward[ki_back][0]-2*dds_backward[ki_back][0]*delta_s);

        if (ds_backward[ki_back-1][0]>ds_upBound[ki_back-1][0])
        {
            stopFlag=true;
            stop_back=ki_back;
            printf("StanceLeg Backward Integration ends at k=%d, ds_backward:%.4f; ds_backward[ki_back-1][0]=%.4f > ds_upBound[ki_back-1][0]=%.4f\n",
                   ki_back,ds_backward[ki_back][0],ds_backward[ki_back-1][0],ds_upBound[ki_back-1][0]);
        }
        else
        {
            ki_back--;
        }
    }

    //forward integration
    stopFlag=false;
    accFlag=true;
    ki_for=0;
    cycleCount=0;
    ds_forward[ki_for][0]=ds_upBound[ki_for][0];
    std::fill_n(min_dist,sTotalCount,1);
    while (stopFlag==false)
    {
        if (accFlag==true)
        {
            double acc[9] {0};
            for (int j=0;j<3;j++)
            {
                for (int k=0;k<3;k++)
                {
                    acc[3*j+k]=output_a2[ki_for][6*j+3+k]*ds_forward[ki_for][0]*ds_forward[ki_for][0]
                              +output_a0H[ki_for][6*j+3+k];
                }
            }
            dds_forward[ki_for][0]=*std::min_element(acc,acc+9);
            ds_forward[ki_for+1][0]=sqrt(ds_forward[ki_for][0]*ds_forward[ki_for][0]+2*dds_forward[ki_for][0]*delta_s);
            //printf("acc,sqrt:%.4f\n",ds_forward[ki_for][0]*ds_forward[ki_for][0]+2*dds_forward[ki_for][0]*delta_s);

            if (ds_forward[ki_for+1][0]>ds_upBound[ki_for+1][0])
            {
                accFlag=false;
                dec_start=ki_for;
                printf("StanceLeg acc reach bound at k=%d, ds_forward=%.4f; ",ki_for,ds_forward[ki_for][0]);
            }
            else
            {
                ki_for++;
            }
        }
        else
        {
            double dec[9] {0};
            for (int j=0;j<3;j++)
            {
                for (int k=0;k<3;k++)
                {
                    dec[3*j+k]=output_a2[ki_for][6*j+3+k]*ds_forward[ki_for][0]*ds_forward[ki_for][0]
                              +output_a0L[ki_for][6*j+3+k];
                }
            }
            dds_forward[ki_for][0]=*std::max_element(dec,dec+9);
            ds_forward[ki_for+1][0]=sqrt(ds_forward[ki_for][0]*ds_forward[ki_for][0]+2*dds_forward[ki_for][0]*delta_s);

            if (ds_forward[ki_for+1][0]>ds_upBound[ki_for+1][0])
            {
                //printf("dec trying\n");
                if(dec_start>0)
                {
                    dec_start--;
                    ki_for=dec_start;
                }
                else
                {
                    printf("dec_start=%d\n",dec_start);
                    ki_for=dec_start;
                    ds_forward[0][0]-=ds_upBound[0][0]/1000;
                }
            }
            else
            {
                if ((ds_forward[ki_for][0]*ds_forward[ki_for][0]+2*dds_forward[ki_for][0]*delta_s)<=(ds_lowBound[ki_for+1][0]*ds_lowBound[ki_for+1][0]) || ki_for==(sTotalCount-2))
                {

                    accFlag=true;
                    for(int k=dec_start;k<(ki_for+2);k++)
                    {
                        min_dist[k]=ds_upBound[k][0]-ds_forward[k][0];
                    }
                    dec_end=std::min_element(min_dist+dec_start+1,min_dist+ki_for+2)-min_dist;
                    //dec_start must be ignored, if dec_start is the min_dist, the calculation will cycle between dec_start & dec_start+1
                    ki_for=dec_end-1;
                    printf("dec finished, start at k=%d, end at k=%d, ds_forward=%.4f\n",dec_start,dec_end,ds_forward[ki_for+1][0]);
                }
                ki_for++;
            }
        }

        //loop end condition
        if(ki_for==sTotalCount-1)
        {
            stopFlag=true;
            for(int i=0;i<sTotalCount;i++)
            {
                real_ds[i][0]=ds_forward[i][0];
                real_dds[i][0]=dds_forward[i][0];
                //printf("i=%d, ds_forward:%.4f, dds_forward:%.4f\n",i,ds_forward[i][0],dds_forward[i][0]);
            }
            printf("StanceLeg forward reach the end, and never encounter with the backward, cycleCount=%u < %u\n",cycleCount,0x0FFFFFFF);
        }
        if(ki_for>=stop_back && ds_forward[ki_for][0]>=ds_backward[ki_for][0])
        {
            stopFlag=true;
            for(int i=0;i<ki_for;i++)
            {
                real_ds[i][0]=ds_forward[i][0];
                real_dds[i][0]=dds_forward[i][0];
            }
            for(int i=ki_for;i<sTotalCount;i++)
            {
                real_ds[i][0]=ds_backward[i][0];
                real_dds[i][0]=dds_backward[i][0];
            }
            printf("StanceLeg forward & backward encounters at k=%d\n",ki_for);
        }
        cycleCount++;
        if(cycleCount==0x0FFFFFFF)
        {
            stopFlag=true;
            printf("WARNING!!! StanceLeg integration takes too long, force stop!!! ki=%d\n",cycleCount);
        }
    }*/


    //pb_sw, vb_sw, ab_sw initialized here, and need to be updated during iteration
    totalTime[0]=0;
    for (int i=0;i<sTotalCount;i++)
    {
        if(i!=0)
        {
            totalTime[0]+=2*delta_s/(real_ds[i-1][0]+real_ds[i][0]);
            timeArray[i][0]=totalTime[0];
        }
        pb_sw[i]=b_st[i];
        vb_sw[i]=db_st[i]*real_ds[i][0];
        ab_sw[i]=ddb_st[i]*real_ds[i][0]*real_ds[i][0]+db_st[i]*real_dds[i][0];
    }


    /******************************* Iteration ************************************/
    int iterCount {0};
    while(maxTotalCount!=maxTotalCount_last)
    {
        iterCount++;
        maxTotalCount_last=maxTotalCount;
        /********** SwingLeg : generate the traj & calculate the bound of ds **********/
        for (int i=0;i<sTotalCount;i++)
        {
            s_w=i*PI/(sTotalCount-1);
            f_sw[0]=0;
            f_sw[1]=stepH*sin(s_w);
            f_sw[2]=stepD/2*cos(s_w);
            df_sw[0]=0;
            df_sw[1]=stepH*cos(s_w);
            df_sw[2]=-stepD/2*sin(s_w);
            ddf_sw[0]=0;
            ddf_sw[1]=-stepH*sin(s_w);
            ddf_sw[2]=-stepD/2*cos(s_w);

            pEB[2]=initPeb[2]+pb_sw[i];
            rbt.SetPeb(pEB);
            for (int j=0;j<3;j++)//pEE of swing leg
            {
                pEE[6*j]=initPee[6*j]+f_sw[0];
                pEE[6*j+1]=initPee[6*j+1]+f_sw[1];
                pEE[6*j+2]=initPee[6*j+2]+f_sw[2];
                rbt.pLegs[j]->SetPee(pEE+6*j);
            }

            //calculate param of swing leg
            for (int j=0;j<3;j++)
            {
                rbt.pLegs[2*j]->GetJvi(Jvi,rbt.body());
                rbt.pLegs[2*j]->GetdJacOverPee(dJvi_x,dJvi_y,dJvi_z,"B");
                rbt.pLegs[2*j]->GetPee(*output_PeeB+18*i+6*j,rbt.body());

                double vb_sw_tmp[3] {0};
                vb_sw_tmp[2]=vb_sw[i];
                double ab_sw_tmp[3] {0};
                ab_sw_tmp[2]=ab_sw[i];

                std::fill_n(dJvi_dot_f,9,0);
                aris::dynamic::s_daxpy(9,df_sw[0],dJvi_x,1,dJvi_dot_f,1);
                aris::dynamic::s_daxpy(9,df_sw[1],dJvi_y,1,dJvi_dot_f,1);
                aris::dynamic::s_daxpy(9,df_sw[2],dJvi_z,1,dJvi_dot_f,1);

                std::fill_n(dJvi_dot_vb,9,0);
                aris::dynamic::s_daxpy(9,vb_sw_tmp[0],dJvi_x,1,dJvi_dot_vb,1);
                aris::dynamic::s_daxpy(9,vb_sw_tmp[1],dJvi_y,1,dJvi_dot_vb,1);
                aris::dynamic::s_daxpy(9,vb_sw_tmp[2],dJvi_z,1,dJvi_dot_vb,1);

                std::fill_n(param_dds+6*j,3,0);
                aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,df_sw,1,1,param_dds+6*j,1);

                std::fill_n(param_dsds1+6*j,3,0);
                aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,ddf_sw,1,1,param_dsds1+6*j,1);
                std::fill_n(param_dsds2+6*j,3,0);
                aris::dynamic::s_dgemm(3,1,3,1,dJvi_dot_f,3,df_sw,1,1,param_dsds2+6*j,1);

                std::fill_n(param_ds1+6*j,3,0);
                aris::dynamic::s_dgemm(3,1,3,1,dJvi_dot_vb,3,df_sw,1,1,param_ds1+6*j,1);
                std::fill_n(param_ds2+6*j,3,0);
                aris::dynamic::s_dgemm(3,1,3,1,dJvi_dot_f,3,vb_sw_tmp,1,1,param_ds2+6*j,1);

                std::fill_n(param_const1+6*j,3,0);
                aris::dynamic::s_dgemm(3,1,3,1,dJvi_dot_vb,3,vb_sw_tmp,1,1,param_const1+6*j,1);
                std::fill_n(param_const2+6*j,3,0);
                aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,ab_sw_tmp,1,1,param_const2+6*j,1);

                std::fill_n(param_dsds+6*j,3,0);
                std::fill_n(param_ds+6*j,3,0);
                for (int k=0;k<3;k++)
                {
                    param_dsds[6*j+k]=param_dsds1[6*j+k]+param_dsds2[6*j+k];
                    param_ds[6*j+k]=-param_ds1[6*j+k]-param_ds2[6*j+k];
                    param_const[6*j+k]=param_const1[6*j+k]-param_const2[6*j+k];

                    if(param_dds[6*j+k]>0)
                    {
                        param_a2[6*j+k]=-param_dsds[6*j+k]/param_dds[6*j+k];
                        param_a1[6*j+k]=-param_ds[6*j+k]/param_dds[6*j+k];
                        param_a0L[6*j+k]=(-aLmt-param_const[6*j+k])/param_dds[6*j+k];
                        param_a0H[6*j+k]=(aLmt-param_const[6*j+k])/param_dds[6*j+k];
                    }
                    else if(param_dds[6*j+k]<0)
                    {
                        param_a2[6*j+k]=-param_dsds[6*j+k]/param_dds[6*j+k];
                        param_a1[6*j+k]=-param_ds[6*j+k]/param_dds[6*j+k];
                        param_a0L[6*j+k]=(aLmt-param_const[6*j+k])/param_dds[6*j+k];
                        param_a0H[6*j+k]=(-aLmt-param_const[6*j+k])/param_dds[6*j+k];
                    }
                    else
                    {
                        printf("WARNING!!! param_dds equals zero!!! SwingLeg : i=%d \n",i);
                    }
                }
            }

            /********** calculate ds bound of swingLegs 1 by 1 **********/
            for (int j=0;j<3;j++)
            {
                bool ds_lowBoundFlag_sw {false};
                bool ds_upBoundFlag_sw {false};
                int k_sw {0};
                const int kswCount {30000};
                while (ds_upBoundFlag_sw==false)
                {
                    double ds=0.001*k_sw;
                    double dec[3] {0};
                    double acc[3] {0};
                    double max_dec {0};
                    double min_acc {0};
                    for (int k=0;k<3;k++)
                    {
                        if(param_dds[6*j+k]>0)
                        {
                            dec[k]=(-param_dsds[6*j+k]*ds*ds-param_ds[6*j+k]*ds-param_const[6*j+k]-aLmt)/param_dds[6*j+k];
                            acc[k]=(-param_dsds[6*j+k]*ds*ds-param_ds[6*j+k]*ds-param_const[6*j+k]+aLmt)/param_dds[6*j+k];
                        }
                        else if(param_dds[6*j+k]<0)
                        {
                            dec[k]=(-param_dsds[6*j+k]*ds*ds-param_ds[6*j+k]*ds-param_const[6*j+k]+aLmt)/param_dds[6*j+k];
                            acc[k]=(-param_dsds[6*j+k]*ds*ds-param_ds[6*j+k]*ds-param_const[6*j+k]-aLmt)/param_dds[6*j+k];
                        }
                        //dec[k]=param_a2[6*j+k]*ds*ds+param_a1[6*j+k]*ds+param_a0L[6*j+k];
                        //acc[k]=param_a2[6*j+k]*ds*ds+param_a1[6*j+k]*ds+param_a0H[6*j+k];
                    }
                    max_dec=*std::max_element(dec,dec+3);
                    min_acc=*std::min_element(acc,acc+3);

                    k_sw++;
                    if(k_sw==kswCount)
                    {
                        ds_upBoundFlag_sw=true;
                        printf("WARNING!!! kswCount=%d is too small!!! Leg:%d, i=%d\n",kswCount,2*j,i);
                    }
                    else
                    {
                        if(ds_lowBoundFlag_sw==false && ds_upBoundFlag_sw==false && min_acc>max_dec)
                        {
                            ds_lowBound_aLmt[i][j+1]=ds;

                            if(k_sw==1)
                            {
                                //printf("swing ds_lowBound_aLmt=%.6f\n",ds_lowBound_aLmt[i][j+1]);
                                ds_lowBoundFlag_sw=true;
                            }
                            else
                            {
                                for(int k_sw2=0;k_sw2<10000;k_sw2++)
                                {
                                    ds=ds_lowBound_aLmt[i][j+1]-0.001*0.0001*k_sw2;
                                    for (int k=0;k<3;k++)
                                    {
                                        if(param_dds[6*j+k]>0)
                                        {
                                            dec[k]=(-param_dsds[6*j+k]*ds*ds-param_ds[6*j+k]*ds-param_const[6*j+k]-aLmt)/param_dds[6*j+k];
                                            acc[k]=(-param_dsds[6*j+k]*ds*ds-param_ds[6*j+k]*ds-param_const[6*j+k]+aLmt)/param_dds[6*j+k];
                                        }
                                        else if(param_dds[6*j+k]<0)
                                        {
                                            dec[k]=(-param_dsds[6*j+k]*ds*ds-param_ds[6*j+k]*ds-param_const[6*j+k]+aLmt)/param_dds[6*j+k];
                                            acc[k]=(-param_dsds[6*j+k]*ds*ds-param_ds[6*j+k]*ds-param_const[6*j+k]-aLmt)/param_dds[6*j+k];
                                        }
                                        //dec[k]=param_a2[6*j+k]*ds*ds+param_a1[6*j+k]*ds+param_a0L[6*j+k];
                                        //acc[k]=param_a2[6*j+k]*ds*ds+param_a1[6*j+k]*ds+param_a0H[6*j+k];
                                    }
                                    max_dec=*std::max_element(dec,dec+3);
                                    min_acc=*std::min_element(acc,acc+3);

                                    if(min_acc<max_dec)
                                    {
                                        //printf("swing ds_lowBound_aLmt=%.6f\n",ds_lowBound_aLmt[i][j+1]);
                                        ds_lowBoundFlag_sw=true;
                                        break;
                                    }
                                    else
                                    {
                                        ds_lowBound_aLmt[i][j+1]=ds;
                                    }
                                }
                            }
                        }
                        else if(ds_lowBoundFlag_sw==true && ds_upBoundFlag_sw==false && min_acc>=max_dec)
                        {
                            dds_lowBound[i][j+1]=max_dec;
                            dds_upBound[i][j+1]=min_acc;
                            ds_upBound_aLmt[i][j+1]=ds;
                        }
                        else if(ds_lowBoundFlag_sw==true && ds_upBoundFlag_sw==false && min_acc<max_dec)
                        {
                            for(int k_sw2=0;k_sw2<10000;k_sw2++)
                            {
                                ds=ds_upBound_aLmt[i][j+1]+0.001*0.0001*k_sw2;
                                for (int k=0;k<3;k++)
                                {
                                    if(param_dds[6*j+k]>0)
                                    {
                                        dec[k]=(-param_dsds[6*j+k]*ds*ds-param_ds[6*j+k]*ds-param_const[6*j+k]-aLmt)/param_dds[6*j+k];
                                        acc[k]=(-param_dsds[6*j+k]*ds*ds-param_ds[6*j+k]*ds-param_const[6*j+k]+aLmt)/param_dds[6*j+k];
                                    }
                                    else if(param_dds[6*j+k]<0)
                                    {
                                        dec[k]=(-param_dsds[6*j+k]*ds*ds-param_ds[6*j+k]*ds-param_const[6*j+k]+aLmt)/param_dds[6*j+k];
                                        acc[k]=(-param_dsds[6*j+k]*ds*ds-param_ds[6*j+k]*ds-param_const[6*j+k]-aLmt)/param_dds[6*j+k];
                                    }
                                    //dec[k]=param_a2[6*j+k]*ds*ds+param_a1[6*j+k]*ds+param_a0L[6*j+k];
                                    //acc[k]=param_a2[6*j+k]*ds*ds+param_a1[6*j+k]*ds+param_a0H[6*j+k];
                                }
                                max_dec=*std::max_element(dec,dec+3);
                                min_acc=*std::min_element(acc,acc+3);

                                if(min_acc<max_dec)
                                {
                                    //printf("swing ds_upBound_aLmt=%.6f,i=%d,j=%d\n",ds_upBound_aLmt[i][j+1],i,j);
                                    ds_upBoundFlag_sw=true;
                                    break;
                                }
                                else
                                {
                                    dds_lowBound[i][j+1]=max_dec;
                                    dds_upBound[i][j+1]=min_acc;
                                    ds_upBound_aLmt[i][j+1]=ds;
                                }
                            }
                        }
                    }
                }

                double vLmt_value[3] {0};
                double Jvi_dot_vb[3] {0};
                double vb_sw_tmp[3] {0};
                vb_sw_tmp[2]=vb_sw[i];
                aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,vb_sw_tmp,1,1,Jvi_dot_vb,1);
                for (int k=0;k<3;k++)
                {
                    if(param_dds[6*j+k]>0)
                    {
                        vLmt_value[k]=(Jvi_dot_vb[k]+vLmt)/param_dds[6*j+k];
                    }
                    else if(param_dds[6*j+k]<0)
                    {
                        vLmt_value[k]=(Jvi_dot_vb[k]-vLmt)/param_dds[6*j+k];
                    }
                    else
                    {
                        printf("WARNING!!! param_dds equals zero!!! SwingLeg : i=%d \n",i);
                    }

                }
                ds_upBound_vLmt[i][j+1]=*std::max_element(vLmt_value,vLmt_value+3);
                ds_upBound[i][j+1]=std::min(ds_upBound_aLmt[i][j+1],ds_upBound_vLmt[i][j+1]);
                ds_lowBound[i][j+1]=ds_lowBound_aLmt[i][j+1];

                memcpy(*output_dsds1+18*i+6*j, param_dsds1+6*j, 3*sizeof(double));
                memcpy(*output_dsds2+18*i+6*j, param_dsds2+6*j, 3*sizeof(double));
                memcpy(*output_dsds+18*i+6*j,  param_dsds+6*j,  3*sizeof(double));
                memcpy(*output_ds1+18*i+6*j, param_ds1+6*j, 3*sizeof(double));
                memcpy(*output_ds2+18*i+6*j, param_ds2+6*j, 3*sizeof(double));
                memcpy(*output_ds+18*i+6*j,  param_ds+6*j,  3*sizeof(double));
                memcpy(*output_const1+18*i+6*j, param_const1+6*j, 3*sizeof(double));
                memcpy(*output_const2+18*i+6*j, param_const2+6*j, 3*sizeof(double));
                memcpy(*output_const+18*i+6*j,  param_const+6*j,  3*sizeof(double));
                memcpy(*output_dds+18*i+6*j, param_dds+6*j, 3*sizeof(double));
                memcpy(*output_a2+18*i+6*j,  param_a2+6*j,  3*sizeof(double));
                memcpy(*output_a1+18*i+6*j,  param_a1+6*j,  3*sizeof(double));
                memcpy(*output_a0L+18*i+6*j, param_a0L+6*j, 3*sizeof(double));
                memcpy(*output_a0H+18*i+6*j, param_a0H+6*j, 3*sizeof(double));
            }
        }

        /********** SwingLeg : find switch point - param_dds=0 point & tangent point **********/
        for(int j=0;j<3;j++)
        {
            //initialize
            for(int k=0;k<3;k++)
            {
                paramdds0Count[6*j+k]=0;
            }
            tangentCount[j+1]=0;
            switchCount[j+1]=0;
            for(int i=0;i<sTotalCount;i++)
            {
                tangentPoint[i][j+1]=-1;
                switchPoint[i][j+1]=-1;
                for(int k=0;k<3;k++)
                {
                    paramdds0Point[i][6*j+k]=-1;
                }
            }

            slopedsBound[0][j+1]=(ds_upBound[1][j+1]-ds_upBound[0][j+1])/delta_s;
            for(int i=1;i<sTotalCount-1;i++)
            {
                slopedsBound[i][j+1]=(ds_upBound[i+1][j+1]-ds_upBound[i-1][j+1])/2/delta_s;
            }
            slopedsBound[sTotalCount-1][j+1]=(ds_upBound[sTotalCount-1][j+1]-ds_upBound[sTotalCount-2][j+1])/delta_s;

            for(int i=0;i<sTotalCount;i++)
            {
                slopeDelta[i][j+1]=slopedsBound[i][j+1]-(dds_upBound[i][j+1]+dds_lowBound[i][j+1])/2;
            }

            for(int i=0;i<sTotalCount-1;i++)
            {
                for(int k=0;k<3;k++)
                {
                    if(output_dds[i+1][6*j+k]*output_dds[i][6*j+k]<0 || output_dds[i][6*j+k]==0)
                    {
                        paramdds0Point[paramdds0Count[6*j+k]][6*j+k]=i;
                        paramdds0Count[6*j+k]++;
                    }
                }

                if(slopeDelta[i+1][j+1]*slopeDelta[i][j+1]<0 || slopeDelta[i][j+1]==0)
                {
                    tangentPoint[tangentCount[j+1]][j+1]=i;
                    tangentCount[j+1]++;
                }
            }

            //merge tangentPoint & paramdds0Point into switchPoint
            for(int i=0;i<tangentCount[j+1];i++)
            {
                switchPoint[i][j+1]=tangentPoint[i][j+1];
            }
            switchCount[j+1]=tangentCount[j+1];
            for(int k=0;k<3;k++)
            {
                for(int i=0;i<paramdds0Count[6*j+k];i++)
                {
                    switchPoint[i+switchCount[j+1]][j+1]=paramdds0Point[i][6*j+k];
                }
                switchCount[j+1]+=paramdds0Count[6*j+k];
            }

            //filtering the same point & sorting by the value
            for(int i=0;i<switchCount[j+1];i++)
            {
                for(int k=i+1;k<switchCount[j+1];k++)
                {
                    if(switchPoint[k][j+1]<switchPoint[i][j+1])
                    {
                        int tmp;
                        tmp=switchPoint[i][j+1];
                        switchPoint[i][j+1]=switchPoint[k][j+1];
                        switchPoint[k][j+1]=tmp;
                    }
                    else if(switchPoint[k][j+1]==switchPoint[i][j+1])
                    {
                        switchPoint[k][j+1]=switchPoint[switchCount[j+1]-1][j+1];
                        switchPoint[switchCount[j+1]-1][j+1]=-1;
                        k--;
                        switchCount[j+1]--;
                    }
                }
            }

            printf("\nSwingLeg Switch Point:");
            for(int i=0;i<switchCount[j+1]+1;i++)
            {
                printf("%d,",switchPoint[i][j+1]);
            }
            printf("\n");
        }

        /********** SwingLeg : numerical integration to calculate ds**********/

        for(int j=0;j<3;j++)
        {
            for(int i=0;i<sTotalCount;i++)
            {
                real_ds[i][j+1]=ds_upBound[i][j+1];
                real_dds[i][j+1]=dds_upBound[i][j+1];
            }
            for(int m=0;m<switchCount[j+1];m++)
            {
                int k_sw {switchPoint[m][j+1]};
                stopFlag=false;
                ds_backward[k_sw][j+1]=ds_upBound[k_sw][j+1];
                while(stopFlag==false)
                {
                    double dec[3] {0};
                    for (int k=0;k<3;k++)
                    {
                        if(output_dds[k_sw][6*j+k]>0)
                        {
                            dec[k]=(-output_dsds[k_sw][6*j+k]*ds_backward[k_sw][j+1]*ds_backward[k_sw][j+1]
                                    -output_ds[k_sw][6*j+k]*ds_backward[k_sw][j+1]
                                    -output_const[k_sw][6*j+k]-aLmt)/output_dds[k_sw][6*j+k];
                        }
                        else if(output_dds[k_sw][6*j+k]<0)
                        {
                            dec[k]=(-output_dsds[k_sw][6*j+k]*ds_backward[k_sw][j+1]*ds_backward[k_sw][j+1]
                                    -output_ds[k_sw][6*j+k]*ds_backward[k_sw][j+1]
                                    -output_const[k_sw][6*j+k]+aLmt)/output_dds[k_sw][6*j+k];
                        }
                        //dec[k]=output_a2[k_sw][6*j+k]*ds_backward[k_sw][j+1]*ds_backward[k_sw][j+1]
                         //     +output_a1[k_sw][6*j+k]*ds_backward[k_sw][j+1]
                           //   +output_a0L[k_sw][6*j+k];
                    }
                    dds_backward[k_sw][j+1]=*std::max_element(dec,dec+3);
                    ds_backward[k_sw-1][j+1]=sqrt(ds_backward[k_sw][j+1]*ds_backward[k_sw][j+1]-2*dds_backward[k_sw][j+1]*delta_s);

                    if(ds_backward[k_sw-1][j+1]>ds_upBound[k_sw-1][j+1])
                    {
                        stopFlag=true;
                        printf("SwingLeg backward touching upBound at %d, quit switchPoint %d\n",k_sw-1,switchPoint[m][j+1]);
                    }
                    else if(k_sw==1)
                    {
                        for (int k=0;k<3;k++)
                        {
                            if(output_dds[k_sw][6*j+k]>0)
                            {
                                dec[k]=(-output_dsds[k_sw][6*j+k]*ds_backward[k_sw][j+1]*ds_backward[k_sw][j+1]
                                        -output_ds[k_sw][6*j+k]*ds_backward[k_sw][j+1]
                                        -output_const[k_sw][6*j+k]-aLmt)/output_dds[k_sw][6*j+k];
                            }
                            else if(output_dds[k_sw][6*j+k]<0)
                            {
                                dec[k]=(-output_dsds[k_sw][6*j+k]*ds_backward[k_sw][j+1]*ds_backward[k_sw][j+1]
                                        -output_ds[k_sw][6*j+k]*ds_backward[k_sw][j+1]
                                        -output_const[k_sw][6*j+k]+aLmt)/output_dds[k_sw][6*j+k];
                            }
                            //dec[k]=output_a2[k_sw][6*j+k]*ds_backward[k_sw][j+1]*ds_backward[k_sw][j+1]
                              //    +output_a1[k_sw][6*j+k]*ds_backward[k_sw][j+1]
                                //  +output_a0L[k_sw][6*j+k];
                        }
                        dds_backward[k_sw-1][j+1]=*std::max_element(dec,dec+3);
                        for(int i=k_sw-1;i<switchPoint[m][j+1]+1;i++)
                        {
                            real_ds[i][j+1]=ds_backward[i][j+1];
                            real_dds[i][j+1]=dds_backward[i][j+1];
                        }
                        stopFlag=true;
                        printf("SwingLeg backward touching 0, from switchPoint %d\n",switchPoint[m][j+1]);
                    }
                    else if(ds_backward[k_sw-1][j+1]>=real_ds[k_sw-1][j+1])
                    {
                        for (int k=0;k<3;k++)
                        {
                            if(output_dds[k_sw][6*j+k]>0)
                            {
                                dec[k]=(-output_dsds[k_sw][6*j+k]*ds_backward[k_sw][j+1]*ds_backward[k_sw][j+1]
                                        -output_ds[k_sw][6*j+k]*ds_backward[k_sw][j+1]
                                        -output_const[k_sw][6*j+k]-aLmt)/output_dds[k_sw][6*j+k];
                            }
                            else if(output_dds[k_sw][6*j+k]<0)
                            {
                                dec[k]=(-output_dsds[k_sw][6*j+k]*ds_backward[k_sw][j+1]*ds_backward[k_sw][j+1]
                                        -output_ds[k_sw][6*j+k]*ds_backward[k_sw][j+1]
                                        -output_const[k_sw][6*j+k]+aLmt)/output_dds[k_sw][6*j+k];
                            }
                            //dec[k]=output_a2[k_sw-1][6*j+k]*ds_backward[k_sw-1][j+1]*ds_backward[k_sw-1][j+1]
                              //    +output_a1[k_sw-1][6*j+k]*ds_backward[k_sw-1][j+1]
                                //  +output_a0L[k_sw-1][6*j+k];
                        }
                        dds_backward[k_sw-1][j+1]=*std::max_element(dec,dec+3);
                        for(int i=k_sw-1;i<switchPoint[m][j+1]+1;i++)
                        {
                            real_ds[i][j+1]=ds_backward[i][j+1];
                            real_dds[i][j+1]=dds_backward[i][j+1];
                        }
                        stopFlag=true;
                        printf("SwingLeg backward touching last curve at %d, from switchPoint %d\n",k_sw-1,switchPoint[m][j+1]);
                    }
                    else
                    {
                        k_sw--;
                    }
                }
                if(ds_backward[k_sw-1][j+1]>ds_upBound[k_sw-1][j+1])
                {
                    continue;
                }

                bool isEqual0Point {false};
                for(int k=0;k<3;k++)
                {
                    for(int i=0;i<paramdds0Count[6*j+k];i++)
                    {
                        if(output_dds[switchPoint[i][j+1]][6*j+k]==0 || slopeDelta[switchPoint[i][j+1]][j+1]==0)
                        {
                            isEqual0Point=true;
                        }
                    }
                }
                if(isEqual0Point==false)
                {
                    k_sw=switchPoint[m][j+1]+1;
                }
                else
                {
                    k_sw=switchPoint[m][j+1];
                }
                stopFlag=false;
                ds_forward[k_sw][j+1]=ds_upBound[k_sw][j+1];
                while(stopFlag==false)
                {
                    double acc[3] {0};
                    for (int k=0;k<3;k++)
                    {
                        if(output_dds[k_sw][6*j+k]>0)
                        {
                            acc[k]=(-output_dsds[k_sw][6*j+k]*ds_forward[k_sw][j+1]*ds_forward[k_sw][j+1]
                                    -output_ds[k_sw][6*j+k]*ds_forward[k_sw][j+1]
                                    -output_const[k_sw][6*j+k]+aLmt)/output_dds[k_sw][6*j+k];
                        }
                        else if(output_dds[k_sw][6*j+k]<0)
                        {
                            acc[k]=(-output_dsds[k_sw][6*j+k]*ds_forward[k_sw][j+1]*ds_forward[k_sw][j+1]
                                    -output_ds[k_sw][6*j+k]*ds_forward[k_sw][j+1]
                                    -output_const[k_sw][6*j+k]-aLmt)/output_dds[k_sw][6*j+k];
                        }
                        //acc[k]=output_a2[k_sw][6*j+k]*ds_forward[k_sw][j+1]*ds_forward[k_sw][j+1]
                          //    +output_a1[k_sw][6*j+k]*ds_forward[k_sw][j+1]
                            //  +output_a0H[k_sw][6*j+k];
                    }
                    dds_forward[k_sw][j+1]=*std::min_element(acc,acc+3);
                    ds_forward[k_sw+1][j+1]=sqrt(ds_forward[k_sw][j+1]*ds_forward[k_sw][j+1]+2*dds_forward[k_sw][j+1]*delta_s);

                    if(ds_forward[k_sw+1][j+1]>ds_upBound[k_sw+1][j+1] || k_sw==sTotalCount-2)
                    {
                        for(int i=switchPoint[m][j+1]+1;i<k_sw+1;i++)
                        {
                            real_ds[i][j+1]=ds_forward[i][j+1];
                            real_dds[i][j+1]=dds_forward[i][j+1];
                        }
                        stopFlag=true;
                        printf("SwingLeg forward touching upBound at %d, from switchPoint %d\n",k_sw,switchPoint[m][j+1]);
                    }
                    else
                    {
                        k_sw++;
                    }
                }
            }
            /*
            //backward integration
            stopFlag=false;
            ki_back=sTotalCount-1;
            ds_backward[ki_back][j+1]=0;
            while (stopFlag==false && ki_back>=0)
            {
                double dec[3] {0};
                for (int k=0;k<3;k++)
                {
                    dec[k]=output_a2[ki_back][6*j+k]*ds_backward[ki_back][j+1]*ds_backward[ki_back][j+1]
                          +output_a1[ki_back][6*j+k]*ds_backward[ki_back][j+1]
                          +output_a0L[ki_back][6*j+k];
                }
                dds_backward[ki_back][j+1]=*std::max_element(dec,dec+3);
                ds_backward[ki_back-1][j+1]=sqrt(ds_backward[ki_back][j+1]*ds_backward[ki_back][j+1]-2*dds_backward[ki_back][j+1]*delta_s);

                if (ds_backward[ki_back-1][j+1]>ds_upBound[ki_back-1][j+1])
                {
                    stopFlag=true;
                    stop_back=ki_back;
                    //printf("SwingLeg Backward Iteration ends at k=%d, ds_backward:%.4f\n",ki_back,ds_backward[ki_back][j+1]);
                    //printf("ds_backward[ki_back-1]:%.4f; ds_bound[ki_back-1]:%.4f\n",ds_backward[ki_back-1][j+1],ds_upBound[ki_back-1][j+1]);
                }
                else
                {
                    ki_back--;
                }
            }

            //forward integration
            stopFlag=false;
            accFlag=true;
            ki_for=0;
            cycleCount=0;
            ds_forward[ki_for][j+1]=0;
            std::fill_n(min_dist,sTotalCount,1);
            while (stopFlag==false)
            {
                if (accFlag==true)
                {
                    double acc[3] {0};
                    for (int k=0;k<3;k++)
                    {
                        acc[k]=output_a2[ki_for][6*j+k]*ds_forward[ki_for][j+1]*ds_forward[ki_for][j+1]
                              +output_a1[ki_for][6*j+k]*ds_forward[ki_for][j+1]
                              +output_a0H[ki_for][6*j+k];
                    }
                    dds_forward[ki_for][j+1]=*std::min_element(acc,acc+3);
                    ds_forward[ki_for+1][j+1]=sqrt(ds_forward[ki_for][j+1]*ds_forward[ki_for][j+1]+2*dds_forward[ki_for][j+1]*delta_s);

                    if (ds_forward[ki_for+1][j+1]>ds_upBound[ki_for+1][j+1])
                    {
                        accFlag=false;
                        dec_start=ki_for;
                        //printf("SwingLeg acc touching at k=%d, ds_forward=%.4f; ",ki_for,ds_forward[ki_for][j+1]);
                    }
                    else
                    {
                        ki_for++;
                    }
                }
                else
                {
                    double dec[3] {0};
                    for (int k=0;k<3;k++)
                    {
                        dec[k]=output_a2[ki_for][6*j+k]*ds_forward[ki_for][j+1]*ds_forward[ki_for][j+1]
                              +output_a1[ki_for][6*j+k]*ds_forward[ki_for][j+1]
                              +output_a0L[ki_for][6*j+k];
                    }
                    dds_forward[ki_for][j+1]=*std::max_element(dec,dec+3);
                    ds_forward[ki_for+1][j+1]=sqrt(ds_forward[ki_for][j+1]*ds_forward[ki_for][j+1]+2*dds_forward[ki_for][j+1]*delta_s);

                    if (ds_forward[ki_for+1][j+1]>ds_upBound[ki_for+1][j+1])
                    {
                        //printf("dec trying\n");
                        if(dec_start>0)
                        {
                            dec_start--;
                            ki_for=dec_start;
                        }
                        else
                        {
                            //printf("dec_start=%d, ki_for=%d, ds_forward=%.4f, ds_bound=%.4f\n",dec_start,ki_for+1,ds_forward[ki_for+1][j+1],ds_upBound[ki_for+1][j+1]);
                            ki_for=dec_start;
                            ds_forward[0][j+1]-=ds_upBound[0][j+1]/1000;
                        }
                    }
                    else
                    {
                        //printf("dec ending,ki_for=%d, ds_forward=%.4f\n",ki_for,ds_forward[ki_for+1][j+1]);
                        if ((ds_forward[ki_for][j+1]*ds_forward[ki_for][j+1]+2*dds_forward[ki_for][j+1]*delta_s)<=(ds_lowBound[ki_for+1][j+1]*ds_lowBound[ki_for+1][j+1]) || ki_for==(sTotalCount-2))
                        {
                            accFlag=true;
                            for(int k=dec_start;k<(ki_for+2);k++)
                            {
                                min_dist[k]=ds_upBound[k][j+1]-ds_forward[k][j+1];
                            }
                            dec_end=std::min_element(min_dist+dec_start+1,min_dist+ki_for+2)-min_dist;
                            //dec_start must be ignored, if dec_start is the min_dist, the calculation will cycle between dec_start & dec_start+1
                            ki_for=dec_end-1;
                            //printf("dec finished, start at k=%d, end at k=%d, ds_forward=%.4f\n",dec_start,dec_end,ds_forward[ki_for+1][j+1]);
                        }
                        ki_for++;
                    }
                }

                if(ki_for==sTotalCount-1)
                {
                    stopFlag=true;
                    for(int i=0;i<sTotalCount;i++)
                    {
                        real_ds[i][j+1]=ds_forward[i][j+1];
                        real_dds[i][j+1]=dds_forward[i][j+1];
                    }
                    //printf("SwingLeg forward reach the end, and never encounter with the backward, ki=%u < %u\n",cycleCount,0x0FFFFFFF);
                }
                if(ki_for>=stop_back && ds_forward[ki_for][j+1]>=ds_backward[ki_for][j+1])
                {
                    stopFlag=true;
                    for(int i=0;i<ki_for;i++)
                    {
                        real_ds[i][j+1]=ds_forward[i][j+1];
                        real_dds[i][j+1]=dds_forward[i][j+1];
                    }
                    for(int i=ki_for;i<sTotalCount;i++)
                    {
                        real_ds[i][j+1]=ds_backward[i][j+1];
                        real_dds[i][j+1]=dds_backward[i][j+1];
                    }
                    //printf("SwingLeg forward & backward encounters at k=%d\n",ki_for);
                }
                cycleCount++;
                if(cycleCount==0x0FFFFFFF)
                {
                    stopFlag=true;
                    printf("WARNING!!! SwingLeg integration takes too long, force stop!!! ki=%d\n",cycleCount);
                }
            }*/

            totalTime[j+1]=0;
            for (int i=1;i<sTotalCount;i++)
            {
                totalTime[j+1]+=2*delta_s/(real_ds[i-1][j+1]+real_ds[i][j+1]);
                timeArray[i][j+1]=totalTime[j+1];
            }
        }

        maxTime=*std::max_element(totalTime,totalTime+4);
        maxTotalCount=(int)(maxTime*1000)+1;
        maxTime=0.001*maxTotalCount;
        maxTimeID=std::max_element(totalTime,totalTime+4)-totalTime;
        printf("totalTime: %.4f, %.4f, %.4f, %.4f; maxTime:%.4f\n",totalTime[0],totalTime[1],totalTime[2],totalTime[3],maxTime);

        //update pb_sw, vb_sw, ab_sw here
        //double timeArray_tmp[sTotalCount][4] {0};
        for(int i=0;i<sTotalCount;i++)
        {
            memcpy(*timeArray_tmp+4*i,*timeArray+4*i,4*sizeof(double));
            for(int j=0;j<4;j++)
            {
                timeArray_tmp[i][j]*=maxTime/totalTime[j];
            }
        }
        int j_start {0};
        double pb_sw_tmp[sTotalCount] {0};
        double vb_sw_tmp[sTotalCount] {0};
        double ab_sw_tmp[sTotalCount] {0};
        for(int i=0;i<(sTotalCount-1);i++)
        {
            //if(i%100==99)
                //printf("j_start=%d, ",j_start);
            for(int j=j_start;j<(sTotalCount-1);j++)
            {
                if(timeArray_tmp[i][maxTimeID]>=timeArray_tmp[j][0] && timeArray_tmp[i][maxTimeID]<timeArray_tmp[j+1][0])
                {
                    j_start=j;
                    pb_sw_tmp[i]=(pb_sw[j+1]-pb_sw[j])/(timeArray_tmp[j+1][0]-timeArray_tmp[j][0])*(timeArray_tmp[i][maxTimeID]-timeArray_tmp[j][0])+pb_sw[j];
                    vb_sw_tmp[i]=(vb_sw[j+1]-vb_sw[j])/(timeArray_tmp[j+1][0]-timeArray_tmp[j][0])*(timeArray_tmp[i][maxTimeID]-timeArray_tmp[j][0])+vb_sw[j];
                    vb_sw_tmp[i]*=totalTime[0]/maxTime;
                    ab_sw_tmp[i]=(ab_sw[j+1]-ab_sw[j])/(timeArray_tmp[j+1][0]-timeArray_tmp[j][0])*(timeArray_tmp[i][maxTimeID]-timeArray_tmp[j][0])+ab_sw[j];
                    ab_sw_tmp[i]*=totalTime[0]/maxTime*totalTime[0]/maxTime;
                    break;
                }
            }
            //if(i%100==99)
                //printf("j_end=%d\n",j_start);
        }
        pb_sw_tmp[sTotalCount-1]=pb_sw[sTotalCount-1];
        vb_sw_tmp[sTotalCount-1]=vb_sw[sTotalCount-1]*totalTime[0]/maxTime;
        ab_sw_tmp[sTotalCount-1]=ab_sw[sTotalCount-1]*totalTime[0]/maxTime*totalTime[0]/maxTime;

        memcpy(pb_sw,pb_sw_tmp,sTotalCount*sizeof(double));
        memcpy(vb_sw,vb_sw_tmp,sTotalCount*sizeof(double));
        memcpy(ab_sw,ab_sw_tmp,sTotalCount*sizeof(double));
    }

    printf("Iteration finished, iterCount=%d, maxTime=%.4f, timeArray:%.4f,%.4f\n\n",iterCount,maxTime,timeArray[0][0],timeArray[1][0]);

    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_PeeB.txt",*output_PeeB,sTotalCount,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_dsds1.txt",*output_dsds1,sTotalCount,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_dsds2.txt",*output_dsds2,sTotalCount,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_dsds.txt", *output_dsds, sTotalCount,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_ds1.txt",*output_ds1,sTotalCount,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_ds2.txt",*output_ds2,sTotalCount,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_ds.txt", *output_ds, sTotalCount,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_const1.txt",*output_const1,sTotalCount,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_const2.txt",*output_const2,sTotalCount,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_const.txt", *output_const, sTotalCount,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_dds.txt",*output_dds,sTotalCount,18);

    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a2.txt",*output_a2,sTotalCount,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a1.txt",*output_a1,sTotalCount,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a0L.txt",*output_a0L,sTotalCount,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a0H.txt",*output_a0H,sTotalCount,18);

    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/max_ValueL.txt",*dds_lowBound,sTotalCount,4);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/min_ValueH.txt",*dds_upBound,sTotalCount,4);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_upBound_aLmt.txt",*ds_upBound_aLmt,sTotalCount,4);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_lowBound_aLmt.txt",*ds_lowBound_aLmt,sTotalCount,4);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_bound_vLmt.txt",*ds_upBound_vLmt,sTotalCount,4);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_lowBound.txt",*ds_lowBound,sTotalCount,4);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_upBound.txt",*ds_upBound,sTotalCount,4);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_max.txt",*real_ddsMax,sTotalCount,4);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_min.txt",*real_ddsMin,sTotalCount,4);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_forward.txt",*ds_forward,sTotalCount,4);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_backward.txt",*ds_backward,sTotalCount,4);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_forward.txt",*dds_forward,sTotalCount,4);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_backward.txt",*dds_backward,sTotalCount,4);

    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_ds.txt",*real_ds,sTotalCount,4);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_dds.txt",*real_dds,sTotalCount,4);

    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/timeArray.txt",*timeArray_tmp,sTotalCount,4);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/pb_sw.txt",pb_sw,sTotalCount,1);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/vb_sw.txt",vb_sw,sTotalCount,1);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ab_sw.txt",ab_sw,sTotalCount,1);

    //aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Pee.txt",*output_Pee,sTotalCount,9);
    //aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Pin.txt",*output_Pin,sTotalCount,9);
    //aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Vin.txt",*output_Vin,sTotalCount,9);
    //aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Ain.txt",*output_Ain,sTotalCount,9);


    gettimeofday(&tpend,NULL);
    tused=tpend.tv_sec-tpstart.tv_sec+(double)(tpend.tv_usec-tpstart.tv_usec)/1000000;
    printf("UsedTime:%f\n",tused);
}
