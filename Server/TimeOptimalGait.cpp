#include "TimeOptimalGait.h"

TimeOptimalGait::TimeOptimalGait(){}
TimeOptimalGait::~TimeOptimalGait(){}

void TimeOptimalGait::GetStanceLegParam(int count, int legID, double s)
{
    double pEB[6] {0};
    double pEE[18] {0};

    double dJvi_x[9] {0};
    double dJvi_y[9] {0};
    double dJvi_z[9] {0};
    double dJvi[9] {0};//used in stanceLeg

    double param_dsds1[18] {0};
    double param_dsds2[18] {0};
    double param_dsds[18] {0};
    double param_const[18] {0};//for aLmt of swingLeg

    double b_sb_tmp=-stepD*s;
    double db_sb_tmp=-stepD;
    double ddb_sb_tmp=0;

    pEB[2]=initPeb[2]+b_sb_tmp;
    if(legID%2==1)//135
    {
        pEE[3*legID]=initPee[3*legID];
        pEE[3*legID+1]=initPee[3*legID+1];
        pEE[3*legID+2]=initPee[3*legID+2]-stepD/4;
    }
    else//024
    {
        pEE[3*legID]=initPee[3*legID];
        pEE[3*legID+1]=initPee[3*legID+1];
        pEE[3*legID+2]=initPee[3*legID+2]+stepD/4-stepD;
    }

    rbt.SetPeb(pEB);
    rbt.pLegs[legID]->SetPee(pEE+3*legID);

    rbt.pLegs[legID]->GetJvi(Jvi,rbt.body());
    rbt.pLegs[legID]->GetdJacOverPee(dJvi_x,dJvi_y,dJvi_z,"B");

    double db_sb3[3] {0};//vel and acc in Body CS
    double ddb_sb3[3] {0};
    db_sb3[2]=-db_sb_tmp;
    ddb_sb3[2]=-ddb_sb_tmp;

    std::fill_n(dJvi,9,0);
    aris::dynamic::s_daxpy(9,db_sb3[0],dJvi_x,1,dJvi,1);//for s_b
    aris::dynamic::s_daxpy(9,db_sb3[1],dJvi_y,1,dJvi,1);
    aris::dynamic::s_daxpy(9,db_sb3[2],dJvi_z,1,dJvi,1);

    std::fill_n(param_dds+3*legID,3,0);
    aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,db_sb3,1,1,param_dds+3*legID,1);

    std::fill_n(param_dsds1+3*legID,3,0);
    aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,ddb_sb3,1,1,param_dsds1+3*legID,1);
    std::fill_n(param_dsds2+3*legID,3,0);
    aris::dynamic::s_dgemm(3,1,3,1,dJvi,3,db_sb3,1,1,param_dsds2+3*legID,1);
    std::fill_n(param_dsds+3*legID,3,0);
    for (int k=0;k<3;k++)
    {
        param_dsds[3*legID+k]=param_dsds1[3*legID+k]+param_dsds2[3*legID+k];
        abs_param_dds[3*legID+k]=fabs(param_dds[3*legID+k]);

        if(param_dds[3*legID+k]>0)
        {
            param_a2[3*legID+k]=-param_dsds[3*legID+k]/param_dds[3*legID+k];
            param_a0L[3*legID+k]=(-param_const[3*legID+k]-aLmt)/param_dds[3*legID+k];
            param_a0H[3*legID+k]=(-param_const[3*legID+k]+aLmt)/param_dds[3*legID+k];
        }
        else if(param_dds[3*legID+k]<0)
        {
            param_a2[3*legID+k]=-param_dsds[3*legID+k]/param_dds[3*legID+k];
            param_a0L[3*legID+k]=(-param_const[3*legID+k]+aLmt)/param_dds[3*legID+k];
            param_a0H[3*legID+k]=(-param_const[3*legID+k]-aLmt)/param_dds[3*legID+k];
        }
        else
        {
            if(count!=-1)
            {
//                switchPoint_body[switchCount_body]=count;
//                switchScrewID_body[switchCount_body]=3*legID+k;
//                switchCount_body++;

                isParamddsExact0_body[count][3*legID+k]=1;
                printf("WARNING!!! param_dds equals zero!!! StanceLeg : %d \n",count);
            }
        }
    }

    if(count!=-1)
    {
        b_sb[count]=b_sb_tmp;
        db_sb[count]=db_sb_tmp;
        ddb_sb[count]=ddb_sb_tmp;
        rbt.pLegs[legID]->GetPee(*output_PeeB+18*count+3*legID,rbt.body());
        memcpy( *output_dsds+18*count+3*legID,  param_dsds+3*legID, 3*sizeof(double));
        memcpy(  *output_dds+18*count+3*legID,   param_dds+3*legID, 3*sizeof(double));
        memcpy(   *output_a2+18*count+3*legID,    param_a2+3*legID, 3*sizeof(double));
        memcpy(  *output_a0L+18*count+3*legID,   param_a0L+3*legID, 3*sizeof(double));
        memcpy(  *output_a0H+18*count+3*legID,   param_a0H+3*legID, 3*sizeof(double));
    }
}

double TimeOptimalGait::GetStanceMaxDec(int count, double ds)
{
    double dec[18] {0};
    std::fill_n(dec,18,-1e6);

    for (int j=0;j<6;j++)
    {
        for (int k=0;k<3;k++)
        {
            if(isParamddsExact0_body[count][3*j+k]==-1)
            {
                if(count<901)
                {
                    if(j%2==1)//135
                    {
                        dec[3*j+k]=output_a2[count][3*j+k]*ds*ds+output_a0L[count][3*j+k];
                    }
                }
                else if(count<1351)
                {
                    dec[3*j+k]=output_a2[count][3*j+k]*ds*ds+output_a0L[count][3*j+k];
                }
                else
                {
                    if(j%2==0)//024
                    {
                        dec[3*j+k]=output_a2[count][3*j+k]*ds*ds+output_a0L[count][3*j+k];
                    }
                }
            }
        }
    }

    return *std::max_element(dec,dec+18);
}

double TimeOptimalGait::GetStanceMinAcc(int count, double ds)
{
    double acc[18] {0};
    std::fill_n(acc,18,1e6);

    for (int j=0;j<6;j++)
    {
        for (int k=0;k<3;k++)
        {
            if(isParamddsExact0_body[count][3*j+k]==-1)
            {
                if(count<901)
                {
                    if(j%2==1)
                    {
                        acc[3*j+k]=output_a2[count][3*j+k]*ds*ds+output_a0H[count][3*j+k];
                    }
                }
                else if(count<1351)
                {
                    acc[3*j+k]=output_a2[count][3*j+k]*ds*ds+output_a0H[count][3*j+k];
                }
                else
                {
                    if(j%2==0)
                    {
                        acc[3*j+k]=output_a2[count][3*j+k]*ds*ds+output_a0H[count][3*j+k];
                    }
                }
            }
        }
    }

    return *std::min_element(acc,acc+18);
}

void TimeOptimalGait::GetStanceDsBound(int count)
{
    int k_st {0};
    bool dsBoundFlag_st {false};
    const int kstCount {15000};
    while (dsBoundFlag_st==false)
    {
        double ds=0.001*k_st;
        double max_dec=GetStanceMaxDec(count,ds);
        double min_acc=GetStanceMinAcc(count,ds);

        k_st++;
        if(k_st==kstCount)
        {
            dsBoundFlag_st=true;
            printf("WARNING!!! kstCount=%d is too small!!!\n",kstCount);
        }
        else
        {
            if(min_acc>max_dec)
            {
                ds_upBound_aLmt_body[count]=ds;
            }
            else
            {
                for(int k_st2=0;k_st2<1000;k_st2++)
                {
                    ds=ds_upBound_aLmt_body[count]+0.001*0.001*k_st2;
                    max_dec=GetStanceMaxDec(count,ds);
                    min_acc=GetStanceMinAcc(count,ds);
                    if(min_acc>=max_dec)
                    {
                        dds_lowBound_body[count]=max_dec;
                        dds_upBound_body[count]=min_acc;
                        ds_upBound_aLmt_body[count]=ds;
                    }
                    else
                    {
                        //printf("stance ds_upBound_aLmt=%.6f\t",ds_upBound_aLmt[i][0]);
                        dsBoundFlag_st=true;
                        break;
                    }
                }
            }
        }
    }

    ds_upBound_vLmt_body[count]=vLmt/(*std::max_element(abs_param_dds,abs_param_dds+18));
    ds_upBound_body[count]=std::min(ds_upBound_aLmt_body[count],ds_upBound_vLmt_body[count]);
}

void TimeOptimalGait::GetStanceSwitchPoint()
{
    double slopedsBound_body[2251] {0};
    double paramdds0Point_body[2251] {0};
    int paramdds0Count_body {0};

    double tangentPoint_body[2251] {0};
    int tangentCount_body {0};
    int switchScrewID_body_tmp[2251] {0};

    //initialize
    for(int i=0;i<2251;i++)
    {
        tangentPoint_body[i]=-1;
        paramdds0Point_body[i]=-1;
        switchScrewID_body_tmp[i]=-1;
    }
    tangentCount_body=0;
    paramdds0Count_body=0;

    //calculate the slope of ds_upBound
    slopedsBound_body[0]=(ds_upBound_body[1]-ds_upBound_body[0])/(s_b[1]-s_b[0]);
    for(int i=1;i<2250;i++)
    {
        slopedsBound_body[i]=( (ds_upBound_body[i+1]-ds_upBound_body[i])/(s_b[i+1]-s_b[i])
                              +(ds_upBound_body[i]-ds_upBound_body[i-1])/(s_b[i]-s_b[i-1]) )/2;
    }
    slopedsBound_body[2250]=(ds_upBound_body[2250]-ds_upBound_body[2249])/(s_b[2250]-s_b[2249]);
    for(int i=0;i<2251;i++)
    {
        slopeDelta_body[i]=slopedsBound_body[i]-(dds_upBound_body[i]+dds_lowBound_body[i])/2;
    }

    for(int i=0;i<2250;i++)
    {
        bool isParamdds0 {false};
        if(ds_upBound_vLmt_body[i]>=ds_upBound_aLmt_body[i])
        {
            for(int j=0;j<3;j++)
            {
                int k_start=i<901 ? 3 : 0;
                int k_end=i<1351 ? 6 : 3;
                for(int k=k_start;k<k_end;k++)
                {
                    if(output_dds[i+1][6*j+k]*output_dds[i][6*j+k]<0 || output_dds[i][6*j+k]==0)
                    {
                        paramdds0Point_body[paramdds0Count_body]=i
                                +fabs(output_dds[i][6*j+k])/(fabs(output_dds[i][6*j+k])+fabs(output_dds[i+1][6*j+k]));
                        switchScrewID_body_tmp[paramdds0Count_body]=6*j+k;
                        paramdds0Count_body++;
                        isParamdds0=true;
                    }
                }
            }
        }

        if((slopeDelta_body[i+1]*slopeDelta_body[i]<0 || slopeDelta_body[i]==0) && isParamdds0==false)
        {
            tangentPoint_body[tangentCount_body]=i
                    +fabs(slopeDelta_body[i])/(fabs(slopeDelta_body[i])+fabs(slopeDelta_body[i+1]));
            tangentCount_body++;
        }
    }

    printf("StanceLeg Tangent Switch Point:");
    for(int i=0;i<tangentCount_body+1;i++)
    {
        printf("%.2f,",tangentPoint_body[i]);
    }
    printf("\n");

    printf("StanceLeg ZeroInertia Switch Point:");
    for(int i=0;i<paramdds0Count_body+1;i++)
    {
        printf("%.2f,",paramdds0Point_body[i]);
    }
    printf("\n");

    //merge tangentPoint & paramdds0Point into switchPoint
    switchPoint_body[switchCount_body]=0;
    switchCount_body++;
    switchPoint_body[switchCount_body]=2250;
    switchCount_body++;
    for(int i=0;i<tangentCount_body;i++)
    {
        switchPoint_body[i+switchCount_body]=tangentPoint_body[i];
        //switchScrewID_body[i+switchCount_body]=switchScrewID_body_tmp[i];
    }
    switchCount_body+=tangentCount_body;
    for(int i=0;i<paramdds0Count_body;i++)
    {
        switchPoint_body[i+switchCount_body]=paramdds0Point_body[i];
        switchScrewID_body[i+switchCount_body]=switchScrewID_body_tmp[i];
    }
    switchCount_body+=paramdds0Count_body;

    //filtering the same point & sorting by the value
    for(int i=0;i<switchCount_body;i++)
    {
        for(int j=i+1;j<switchCount_body;j++)
        {
            if(switchPoint_body[j]<switchPoint_body[i])
            {
                auto tmp1=switchPoint_body[i];
                switchPoint_body[i]=switchPoint_body[j];
                switchPoint_body[j]=tmp1;

                auto tmp2=switchScrewID_body[i];
                switchScrewID_body[i]=switchScrewID_body[j];
                switchScrewID_body[j]=tmp2;
            }
//            else if((int)switchPoint_body[j]==(int)switchPoint_body[i])
//            {
//                switchPoint_body[j]=switchPoint_body[switchCount_body-1];
//                switchPoint_body[switchCount_body-1]=-1;
//                j--;
//                switchCount_body--;
//            }
        }
    }

    printf("StanceLeg Switch Point:");
    for(int i=0;i<switchCount_body+1;i++)
    {
        printf("%.2f,",switchPoint_body[i]);
    }
    printf("\n");
}

double TimeOptimalGait::GetStanceSwitchMaxDec(int switchID, double ds)
{
    double dec[18] {0};
    std::fill_n(dec,18,-1e6);

    for (int j=0;j<6;j++)
    {
        for (int k=0;k<3;k++)
        {
            if(switchPoint_body[switchID]<901)
            {
                if(j%2==1)//135
                {
                    dec[3*j+k]=param_a2[3*j+k]*ds*ds+param_a0L[3*j+k];
                }
            }
            else if(switchPoint_body[switchID]<1351)
            {
                dec[3*j+k]=param_a2[3*j+k]*ds*ds+param_a0L[3*j+k];
            }
            else
            {
                if(j%2==0)//024
                {
                    dec[3*j+k]=param_a2[3*j+k]*ds*ds+param_a0L[3*j+k];
                }
            }
        }
    }
    if(switchScrewID_body[switchID]>0)
    {
        dec[switchScrewID_body[switchID]]=-1e6;
    }

    return *std::max_element(dec,dec+18);
}

double TimeOptimalGait::GetStanceSwitchMinAcc(int switchID, double ds)
{
    double acc[18] {0};
    std::fill_n(acc,18,1e6);

    for (int j=0;j<6;j++)
    {
        for (int k=0;k<3;k++)
        {
            if(switchPoint_body[switchID]<901)
            {
                if(j%2==1)
                {
                    acc[3*j+k]=param_a2[3*j+k]*ds*ds+param_a0H[3*j+k];
                }
            }
            else if(switchPoint_body[switchID]<1351)
            {
                acc[3*j+k]=param_a2[3*j+k]*ds*ds+param_a0H[3*j+k];
            }
            else
            {
                if(j%2==0)
                {
                    acc[3*j+k]=param_a2[3*j+k]*ds*ds+param_a0H[3*j+k];
                }
            }
        }
    }
    if(switchScrewID_body[switchID]>0)
    {
        acc[switchScrewID_body[switchID]]=1e6;
    }

    return *std::min_element(acc,acc+18);
}

double TimeOptimalGait::GetStanceSwitchDsBound(int switchID)
{
    double ds_a;
    int k_st {0};
    bool dsBoundFlag_st {false};
    const int kstCount {22500};
    while (dsBoundFlag_st==false)
    {
        ds_a=0.001*k_st;
        double max_dec=GetStanceSwitchMaxDec(switchID,ds_a);
        double min_acc=GetStanceSwitchMinAcc(switchID,ds_a);

        if(min_acc<max_dec)
        {
            k_st--;
            for(int k_st2=0;k_st2<1000;k_st2++)
            {
                ds_a=0.001*k_st+0.001*0.001*k_st2;
                max_dec=GetStanceSwitchMaxDec(switchID,ds_a);
                min_acc=GetStanceSwitchMinAcc(switchID,ds_a);
                if(min_acc<max_dec)
                {
                    ds_a=0.001*k_st+0.001*0.001*(k_st2-1);
                    dsBoundFlag_st=true;
                    break;
                }
            }
        }

        if(k_st==kstCount)
        {
            dsBoundFlag_st=true;
            printf("WARNING!!! kstCount=%d is too small!!!\n",kstCount);
        }
        k_st++;
    }

    double ds_v=vLmt/(*std::max_element(abs_param_dds,abs_param_dds+18));
    return std::min(ds_a,ds_v);
}

void TimeOptimalGait::GetStanceTwoPointAtSwitch(double *lowPoint, double *upPoint)
{
    lowPoint[0]=-1;
    lowPoint[switchCount_body-1]=ds_upBound_body[2250];
    upPoint[0]=ds_upBound_body[0];
    upPoint[switchCount_body-1]=-1;

    for(int i=1;i<switchCount_body-1;i++)
    {
//        double s;
//        if(switchPoint_body[i]<901)
//        {
//            s=s_b1/900*switchPoint_body[i];
//            for(int j=0;j<3;j++)
//            {
//                GetStanceLegParam(-1,2*j+1,s);//135
//            }
//        }
//        else if(switchPoint_body[i]<1351)
//        {
//            s=s_b1+(s_b2-s_b1)/450*(switchPoint_body[i]-900);
//            for(int j=0;j<6;j++)
//            {
//                GetStanceLegParam(-1,j,s);
//            }
//        }
//        else
//        {
//            s=s_b2+(1-s_b2)/900*(switchPoint_body[i]-1350);
//            for(int j=0;j<3;j++)
//            {
//                GetStanceLegParam(-1,2*j,s);//024
//            }
//        }
//        double ds=GetStanceSwitchDsBound(i);
//        double dds_back=GetStanceSwitchMaxDec(i,ds);
//        double dds_for=GetStanceSwitchMinAcc(i,ds);
        int num=round(switchPoint_body[i]);
//        int num=(int)switchPoint_body[i];
        if(slopeDelta_body[num]==0)
        {
            printf("Amazing!!! double equals 0\n");
            upPoint[i]=lowPoint[i]=ds_upBound_body[num];
        }
        else
        {
//            lowPoint[i]=std::min(ds_upBound_body[num-1],sqrt(ds*ds-2*dds_back*(switchPoint_body[i]+1-num)*(s_b[num]-s_b[num-1])));
//            upPoint[i]=std::min(ds_upBound_body[num+1],sqrt(ds*ds+2*dds_for*(num+1-switchPoint_body[i])*(s_b[num+1]-s_b[num])));
            auto tmp=std::min(ds_upBound_body[num-1],ds_upBound_body[num]);
            upPoint[i]=lowPoint[i]=std::min(tmp,ds_upBound_body[num+1]);//std::min(ds_upBound_body[num],ds_upBound_body[num+1]);
        }
    }
}

void TimeOptimalGait::GetStanceOptimalDsBySwitchPoint()
{
    bool stopFlag {false};
    bool quitSwitchPoint {false};
    double *lowPoint=new double [switchCount_body];
    double *upPoint=new double [switchCount_body];
    GetStanceTwoPointAtSwitch(lowPoint,upPoint);

    printf("lowPoint:");
    for(int i=0;i<switchCount_body;i++)
    {
        printf("%.2f,",lowPoint[i]);
    }
    printf("\n");
    printf("upPoint:");
    for(int i=0;i<switchCount_body;i++)
    {
        printf("%.2f,",upPoint[i]);
    }
    printf("\n");

    for(int i=0;i<2251;i++)
    {
        real_ds_body[i]=ds_upBound_body[i];
        real_dds_body[i]=dds_upBound_body[i];
    }
    for(int m=0;m<switchCount_body;m++)
    {
        int k_st=round(switchPoint_body[m]);
        real_dds_body[k_st]=0;

        //start of backward
        if(ds_upBound_body[k_st]>real_ds_body[k_st] && k_st!=2250)
        {
            printf("StanceLeg backward start at a passed point, quit switchPoint %.1f\n",switchPoint_body[m]);
            continue;
        }
        quitSwitchPoint=false;
        stopFlag=false;

        if(k_st!=0 && k_st!=2250)
        {
            k_st--;
        }
        if(k_st!=0)
        {
            ds_backward_body[k_st]=lowPoint[m];
            while(stopFlag==false)
            {
                dds_backward_body[k_st]=GetStanceMaxDec(k_st,ds_backward_body[k_st]);
                ds_backward_body[k_st-1]=sqrt(ds_backward_body[k_st]*ds_backward_body[k_st]-2*dds_backward_body[k_st]*(s_b[k_st]-s_b[k_st-1]));

                if(ds_backward_body[k_st-1]>ds_upBound_body[k_st-1])
                {
                    stopFlag=true;
                    quitSwitchPoint=true;
                    printf("StanceLeg backward touching upBound at %d, quit switchPoint %.1f\n",k_st-1,switchPoint_body[m]);
                }
                else if(k_st==1)
                {
                    dds_backward_body[k_st-1]=GetStanceMaxDec(k_st-1,ds_backward_body[k_st-1]);
                    for(int i=k_st-1;i<round(switchPoint_body[m]);i++)
                    {
                        real_ds_body[i]=ds_backward_body[i];
                        real_dds_body[i]=dds_backward_body[i];
                    }
                    stopFlag=true;
                    printf("StanceLeg backward touching 0, from switchPoint %.1f\n",switchPoint_body[m]);
                }
                else if(ds_backward_body[k_st-1]>=real_ds_body[k_st-1])
                {
                    real_dds_body[k_st-1]=(real_ds_body[k_st-1]-ds_backward_body[k_st])*(real_ds_body[k_st-1]+ds_backward_body[k_st])/2/(s_b[k_st]-s_b[k_st-1]);
                    for(int i=k_st;i<round(switchPoint_body[m]);i++)
                    {
                        real_ds_body[i]=ds_backward_body[i];
                        real_dds_body[i]=dds_backward_body[i];
                    }
                    stopFlag=true;
                    printf("StanceLeg backward touching last curve at %d, from switchPoint %.1f\n",k_st-1,switchPoint_body[m]);
                }
                else
                {
                    k_st--;
                }
            }
            if(ds_backward_body[k_st-1]>ds_upBound_body[k_st-1] || round(switchPoint_body[m])==2250)//end of backward
            {
                continue;
            }
        }

        //start of forward
        k_st=round(switchPoint_body[m]);
        if(k_st!=0)
        {
            k_st++;
        }
        stopFlag=false;
        ds_forward_body[k_st]=upPoint[m];
        while(stopFlag==false)
        {
            dds_forward_body[k_st]=GetStanceMinAcc(k_st,ds_forward_body[k_st]);
            ds_forward_body[k_st+1]=sqrt(ds_forward_body[k_st]*ds_forward_body[k_st]+2*dds_forward_body[k_st]*(s_b[k_st+1]-s_b[k_st]));

            if(ds_forward_body[k_st+1]>ds_upBound_body[k_st+1] || k_st==2249)
            {
                for(int i=round(switchPoint_body[m])+1;i<k_st+1;i++)
                {
                    real_ds_body[i]=ds_forward_body[i];
                    real_dds_body[i]=dds_forward_body[i];
                }
//                if(round(switchPoint_body[m])==0)
//                {
//                    real_ds_body[0]=ds_forward_body[0];
//                    real_dds_body[0]=dds_forward_body[0];
//                }
                if(k_st==2249)
                {
                    if(ds_forward_body[k_st+1]>ds_upBound_body[k_st+1])
                    {
                        real_ds_body[2250]=ds_upBound_body[2250];
                        real_dds_body[2250]=dds_upBound_body[2250];
                    }
                    else
                    {
                        real_ds_body[2250]=ds_forward_body[2250];
                        real_dds_body[2250]=GetStanceMinAcc(2250,ds_forward_body[2250]);
                    }
                }
                stopFlag=true;
                printf("StanceLeg forward touching upBound at %d, from switchPoint %.4f\n",k_st,switchPoint_body[m]);
            }
            else
            {
                k_st++;
            }
        }
    }

    delete [] lowPoint;
    delete [] upPoint;
}

void TimeOptimalGait::GetStanceOptimalDsByDirectNI()
{
    bool stopFlag {false};
    bool accFlag {true};
    int dec_start {0};
    int dec_end {0};
    unsigned int cycleCount {0};
    int stop_back {0};
    int ki_back {0};
    int ki_for {0};

    //backward integration
    ki_back=2250;
    ds_backward_body[ki_back]=ds_upBound_body[ki_back];
    while (stopFlag==false && ki_back>=0)
    {
        dds_backward_body[ki_back]=GetStanceMaxDec(ki_back,ds_backward_body[ki_back]);
        ds_backward_body[ki_back-1]=sqrt(ds_backward_body[ki_back]*ds_backward_body[ki_back]-2*dds_backward_body[ki_back]*(s_b[ki_back]-s_b[ki_back-1]));

        if (ds_backward_body[ki_back-1]>ds_upBound_body[ki_back-1])
        {
            stopFlag=true;
            stop_back=ki_back;
            printf("StanceLeg Backward Integration ends at k=%d, ds_backward:%.4f; ds_backward[ki_back-1][0]=%.4f > ds_upBound[ki_back-1][0]=%.4f\n",
                   ki_back,ds_backward_body[ki_back],ds_backward_body[ki_back-1],ds_upBound_body[ki_back-1]);
        }
        else
        {
            ki_back--;
        }
    }

    //forward integration
    stopFlag=false;
    ki_for=0;
    ds_forward_body[ki_for]=ds_upBound_body[ki_for];
    double min_dist[2251] {0};
    std::fill_n(min_dist,2251,1);
    while (stopFlag==false)
    {
        if (accFlag==true)
        {
            dds_forward_body[ki_for]=GetStanceMinAcc(ki_for,ds_forward_body[ki_for]);
            ds_forward_body[ki_for+1]=sqrt(ds_forward_body[ki_for]*ds_forward_body[ki_for]+2*dds_forward_body[ki_for]*(s_b[ki_for+1]-s_b[ki_for]));
            //printf("acc,sqrt:%.4f\n",ds_forward[ki_for][0]*ds_forward[ki_for][0]+2*dds_forward[ki_for][0]*(s_b[ki_for+1]-s_b[ki_for]));

            if (ds_forward_body[ki_for+1]>ds_upBound_body[ki_for+1])
            {
                accFlag=false;
                dec_start=ki_for;
                printf("StanceLeg acc reach bound at k=%d, ds_forward=%.4f; ",ki_for,ds_forward_body[ki_for]);
            }
            else
            {
                ki_for++;
            }
        }
        else
        {
            dds_forward_body[ki_for]=GetStanceMaxDec(ki_for,ds_forward_body[ki_for]);
            ds_forward_body[ki_for+1]=sqrt(ds_forward_body[ki_for]*ds_forward_body[ki_for]+2*dds_forward_body[ki_for]*(s_b[ki_for+1]-s_b[ki_for]));

            if (ds_forward_body[ki_for+1]>ds_upBound_body[ki_for+1])
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
                    ds_forward_body[0]-=ds_upBound_body[0]/1000;
                }
            }
            else
            {
                if ((ds_forward_body[ki_for]*ds_forward_body[ki_for]+2*dds_forward_body[ki_for]*(s_b[ki_for+1]-s_b[ki_for]))
                        <=0 || ki_for==2249)
                {
                    accFlag=true;
                    for(int k=dec_start;k<(ki_for+2);k++)
                    {
                        min_dist[k]=ds_upBound_body[k]-ds_forward_body[k];
                    }
                    dec_end=std::min_element(min_dist+dec_start+1,min_dist+ki_for+2)-min_dist;
                    //dec_start must be ignored, if dec_start is the min_dist, the calculation will cycle between dec_start & dec_start+1
                    ki_for=dec_end-1;
                    printf("dec finished, start at k=%d, end at k=%d, ds_forward=%.4f\n",dec_start,dec_end,ds_forward_body[ki_for+1]);
                }
                ki_for++;
            }
        }

        //loop end condition
        if(ki_for==2250)
        {
            stopFlag=true;
            for(int i=0;i<2251;i++)
            {
                real_ds_body[i]=ds_forward_body[i];
                real_dds_body[i]=dds_forward_body[i];
                //printf("i=%d, ds_forward:%.4f, dds_forward:%.4f\n",i,ds_forward[i][0],dds_forward[i][0]);
            }
            printf("StanceLeg forward reach the end, and never encounter with the backward, cycleCount=%u < %u\n",cycleCount,0x0FFFFFFF);
        }
        if(ki_for>=stop_back && ds_forward_body[ki_for]>=ds_backward_body[ki_for])
        {
            stopFlag=true;
            for(int i=0;i<ki_for;i++)
            {
                real_ds_body[i]=ds_forward_body[i];
                real_dds_body[i]=dds_forward_body[i];
            }
            for(int i=ki_for;i<2251;i++)
            {
                real_ds_body[i]=ds_backward_body[i];
                real_dds_body[i]=dds_backward_body[i];
            }
            printf("StanceLeg forward & backward encounters at k=%d\n",ki_for);
        }
        cycleCount++;
        if(cycleCount==0x0FFFFFFF)
        {
            stopFlag=true;
            printf("WARNING!!! StanceLeg integration takes too long, force stop!!! ki=%d\n",cycleCount);
        }
    }
}

void TimeOptimalGait::GetStanceOptimalDsByMinorIteration()
{
    double Tsb1 {0};
    double Tsb2 {0};
    s_b1=1-dutyCycle;
    s_b2=dutyCycle;

    int k=0;

    //while((fabs(Tsb1-(1-dutyCycle)*Tstep)>1e-4 || fabs(Tsb2-dutyCycle*Tstep)>1e-4) && k<1000)
    {
        for(int i=0;i<2250;i++)
        {
            if(timeArray_body[i]<=(1-dutyCycle)*Tstep && timeArray_body[i+1]>(1-dutyCycle)*Tstep)
            {
                s_b1=s_b[i]+(s_b[i+1]-s_b[i])/(timeArray_body[i+1]-timeArray_body[i])*((1-dutyCycle)*Tstep-timeArray_body[i]);
            }
            if(timeArray_body[i]<=dutyCycle*Tstep && timeArray_body[i+1]>dutyCycle*Tstep)
            {
                s_b2=s_b[i]+(s_b[i+1]-s_b[i])/(timeArray_body[i+1]-timeArray_body[i])*(dutyCycle*Tstep-timeArray_body[i]);
            }
        }

        Tstep=0;
        switchCount_body=0;
        std::fill_n(switchPoint_body,2251,-1);
        std::fill_n(switchScrewID_body,2251,-1);
        std::fill_n(*isParamddsExact0_body,2251*18,-1);

        for (int i=0;i<2251;i++)
        {
            std::fill_n(abs_param_dds,18,0);
            if(i<901)
            {
                s_b[i]=s_b1/900*i;
                for(int j=0;j<3;j++)
                {
                    GetStanceLegParam(i,2*j+1,s_b[i]);//135
                }
            }
            else if(i<1351)
            {
                s_b[i]=s_b1+(s_b2-s_b1)/450*(i-900);
                for(int j=0;j<6;j++)
                {
                    GetStanceLegParam(i,j,s_b[i]);
                }
            }
            else
            {
                s_b[i]=s_b2+(1-s_b2)/900*(i-1350);
                for(int j=0;j<3;j++)
                {
                    GetStanceLegParam(i,2*j,s_b[i]);//024
                }
            }

            GetStanceDsBound(i);

        }

        GetStanceSwitchPoint();
        GetStanceOptimalDsBySwitchPoint();
        //GetStanceOptimalDsByDirectNI();

        for(int i=0;i<2251;i++)
        {
            if(i!=0)
            {
                Tstep+=2*(s_b[i]-s_b[i-1])/(real_ds_body[i-1]+real_ds_body[i]);
                timeArray_body[i]=Tstep;
            }
        }
        Tsb1=timeArray_body[900];
        Tsb2=timeArray_body[1350];
        k++;

        printf("Tsb1=%.4f,Tsb2=%.4f,Tstep=%.4f\n",Tsb1,Tsb2,Tstep);
    }

    printf("Minor iteration count: %d, s_b1=%.4f, s_b2=%.4f\n",k,s_b1,s_b2);

}

void TimeOptimalGait::GetSwingLegParam(int count, int legID, double sw, double *pva_body)
{
    double pEB[6] {0};
    double pEE[18] {0};

    double dJvi_x[9] {0};
    double dJvi_y[9] {0};
    double dJvi_z[9] {0};
    double dJvi_dot_f[9] {0};//used in swingLeg
    double dJvi_dot_vb[9] {0};//used in swingLeg

    double param_dsds1[18] {0};
    double param_dsds2[18] {0};
    double param_ds1[18] {0};
    double param_ds2[18] {0};
    double param_const1[18] {0};
    double param_const2[18] {0};
    double param_dsds[18] {0};
    double param_ds[18] {0};//for aLmt of swingLeg
    double param_const[18] {0};//for aLmt of swingLeg

    f_sw[0]=0;
    f_sw[1]=stepH*sin(sw);
    f_sw[2]=stepD/2*cos(sw);
    df_sw[0]=0;
    df_sw[1]=stepH*cos(sw);
    df_sw[2]=-stepD/2*sin(sw);
    ddf_sw[0]=0;
    ddf_sw[1]=-stepH*sin(sw);
    ddf_sw[2]=-stepD/2*cos(sw);

    pEB[2]=initPeb[2]+pva_body[0];
    if(legID%2==1)//135
    {
        pEE[3*legID]=initPee[3*legID]+f_sw[0];
        pEE[3*legID+1]=initPee[3*legID+1]+f_sw[1];
        pEE[3*legID+2]=initPee[3*legID+2]-stepD/4+f_sw[2]-stepD/2;
    }
    else//024
    {
        pEE[3*legID]=initPee[3*legID]+f_sw[0];
        pEE[3*legID+1]=initPee[3*legID+1]+f_sw[1];
        pEE[3*legID+2]=initPee[3*legID+2]+stepD/4+f_sw[2]-stepD/2;
    }
    rbt.SetPeb(pEB);
    rbt.pLegs[legID]->SetPee(pEE+3*legID);

    rbt.pLegs[legID]->GetJvi(Jvi,rbt.body());
    rbt.pLegs[legID]->GetdJacOverPee(dJvi_x,dJvi_y,dJvi_z,"B");

    double vb_sw3[3] {0};
    double ab_sw3[3] {0};
    vb_sw3[2]=pva_body[1];
    ab_sw3[2]=pva_body[2];

    std::fill_n(dJvi_dot_f,9,0);
    aris::dynamic::s_daxpy(9,df_sw[0],dJvi_x,1,dJvi_dot_f,1);
    aris::dynamic::s_daxpy(9,df_sw[1],dJvi_y,1,dJvi_dot_f,1);
    aris::dynamic::s_daxpy(9,df_sw[2],dJvi_z,1,dJvi_dot_f,1);

    std::fill_n(dJvi_dot_vb,9,0);
    aris::dynamic::s_daxpy(9,vb_sw3[0],dJvi_x,1,dJvi_dot_vb,1);
    aris::dynamic::s_daxpy(9,vb_sw3[1],dJvi_y,1,dJvi_dot_vb,1);
    aris::dynamic::s_daxpy(9,vb_sw3[2],dJvi_z,1,dJvi_dot_vb,1);

    std::fill_n(param_dds+3*legID,3,0);
    aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,df_sw,1,1,param_dds+3*legID,1);

    std::fill_n(param_dsds1+3*legID,3,0);
    aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,ddf_sw,1,1,param_dsds1+3*legID,1);
    std::fill_n(param_dsds2+3*legID,3,0);
    aris::dynamic::s_dgemm(3,1,3,1,dJvi_dot_f,3,df_sw,1,1,param_dsds2+3*legID,1);

    std::fill_n(param_ds1+3*legID,3,0);
    aris::dynamic::s_dgemm(3,1,3,1,dJvi_dot_vb,3,df_sw,1,1,param_ds1+3*legID,1);
    std::fill_n(param_ds2+3*legID,3,0);
    aris::dynamic::s_dgemm(3,1,3,1,dJvi_dot_f,3,vb_sw3,1,1,param_ds2+3*legID,1);

    std::fill_n(param_const1+3*legID,3,0);
    aris::dynamic::s_dgemm(3,1,3,1,dJvi_dot_vb,3,vb_sw3,1,1,param_const1+3*legID,1);
    std::fill_n(param_const2+3*legID,3,0);
    aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,ab_sw3,1,1,param_const2+3*legID,1);

    std::fill_n(param_dsds+3*legID,3,0);
    std::fill_n(param_ds+3*legID,3,0);
    for (int k=0;k<3;k++)
    {
        param_dsds[3*legID+k]=param_dsds1[3*legID+k]+param_dsds2[3*legID+k];
        param_ds[3*legID+k]=-param_ds1[3*legID+k]-param_ds2[3*legID+k];
        param_const[3*legID+k]=param_const1[3*legID+k]-param_const2[3*legID+k];

        if(param_dds[3*legID+k]>0)
        {
            param_a2[3*legID+k]=-param_dsds[3*legID+k]/param_dds[3*legID+k];
            param_a1[3*legID+k]=-param_ds[3*legID+k]/param_dds[3*legID+k];
            param_a0L[3*legID+k]=(-aLmt-param_const[3*legID+k])/param_dds[3*legID+k];
            param_a0H[3*legID+k]=(aLmt-param_const[3*legID+k])/param_dds[3*legID+k];
        }
        else if(param_dds[3*legID+k]<0)
        {
            param_a2[3*legID+k]=-param_dsds[3*legID+k]/param_dds[3*legID+k];
            param_a1[3*legID+k]=-param_ds[3*legID+k]/param_dds[3*legID+k];
            param_a0L[3*legID+k]=(aLmt-param_const[3*legID+k])/param_dds[3*legID+k];
            param_a0H[3*legID+k]=(-aLmt-param_const[3*legID+k])/param_dds[3*legID+k];
        }
        else
        {
            if(count!=-1)
            {
                int swCount = count<901 ? count : (count-1350);
                isParamddsExact0[swCount][legID]=1;
                printf("WARNING!!! param_dds equals zero!!! SwingLeg : %d \n",count);
            }
        }
    }

    if(count!=-1)
    {
        rbt.pLegs[legID]->GetPee(*output_PeeB+18*count+3*legID,rbt.body());
        memcpy(*output_dsds+18*count+3*legID,  param_dsds+3*legID,  3*sizeof(double));
        memcpy(*output_ds+18*count+3*legID,  param_ds+3*legID,  3*sizeof(double));
        memcpy(*output_const+18*count+3*legID,  param_const+3*legID,  3*sizeof(double));
        memcpy(*output_const1+18*count+3*legID,  param_const1+3*legID,  3*sizeof(double));
        memcpy(*output_const2+18*count+3*legID,  param_const2+3*legID,  3*sizeof(double));
        memcpy(*output_dds+18*count+3*legID, param_dds+3*legID, 3*sizeof(double));
        memcpy(*output_a2+18*count+3*legID,  param_a2+3*legID,  3*sizeof(double));
        memcpy(*output_a1+18*count+3*legID,  param_a1+3*legID,  3*sizeof(double));
        memcpy(*output_a0L+18*count+3*legID, param_a0L+3*legID, 3*sizeof(double));
        memcpy(*output_a0H+18*count+3*legID, param_a0H+3*legID, 3*sizeof(double));
    }
}

double TimeOptimalGait::GetSwingMaxDec(int count, double ds, int legID)
{
    int swCount = count<901 ? count : (count-1350);
    double dec[3] {0};
    std::fill_n(dec,3,-1e6);
    for (int k=0;k<3;k++)
    {
        if(isParamddsExact0[swCount][3*legID+k]==-1)
        {
            dec[k]=output_a2[count][3*legID+k]*ds*ds+output_a1[count][3*legID+k]+output_a0L[count][3*legID+k];
        }
    }
    return *std::max_element(dec,dec+3);
}

double TimeOptimalGait::GetSwingMinAcc(int count, double ds, int legID)
{
    int swCount = count<901 ? count : (count-1350);
    double acc[3] {0};
    std::fill_n(acc,3,1e6);
    for (int k=0;k<3;k++)
    {
        if(isParamddsExact0[swCount][3*legID+k]==-1)
        {
            acc[k]=output_a2[count][3*legID+k]*ds*ds+output_a1[count][3*legID+k]+output_a0H[count][3*legID+k];
        }
    }
    return *std::min_element(acc,acc+3);
}

void TimeOptimalGait::GetSwingDsBound(int count, int legID)
{
    bool ds_lowBoundFlag_sw {false};
    bool ds_upBoundFlag_sw {false};
    int k_sw {0};
    int swCount = count<901 ? count : (count-1350);
    const int kswCount {30000};

    while (ds_upBoundFlag_sw==false)
    {
        double ds=0.001*k_sw;
        double max_dec=GetSwingMaxDec(count,ds,legID);
        double min_acc=GetSwingMinAcc(count,ds,legID);

        k_sw++;
        if(k_sw==kswCount)
        {
            ds_upBoundFlag_sw=true;
            printf("WARNING!!! kswCount=%d is too small!!! Leg:%d, count=%d\n",kswCount,legID,count);
        }
        else
        {
            if(ds_lowBoundFlag_sw==false && ds_upBoundFlag_sw==false && min_acc>max_dec)
            {
                ds_lowBound_aLmt[swCount][legID]=ds;

                if(k_sw==1)
                {
                    ds_lowBoundFlag_sw=true;
                }
                else
                {
                    for(int k_sw2=0;k_sw2<1000;k_sw2++)
                    {
                        ds=ds_lowBound_aLmt[swCount][legID]-0.001*0.001*k_sw2;
                        max_dec=GetSwingMaxDec(count,ds,legID);
                        min_acc=GetSwingMinAcc(count,ds,legID);

                        if(min_acc<max_dec)
                        {
                            ds_lowBoundFlag_sw=true;
                            break;
                        }
                        else
                        {
                            ds_lowBound_aLmt[swCount][legID]=ds;
                        }
                    }
                }
            }
            else if(ds_lowBoundFlag_sw==true && ds_upBoundFlag_sw==false && min_acc>=max_dec)
            {
                dds_lowBound[swCount][legID]=max_dec;
                dds_upBound[swCount][legID]=min_acc;
                ds_upBound_aLmt[swCount][legID]=ds;
            }
            else if(ds_lowBoundFlag_sw==true && ds_upBoundFlag_sw==false && min_acc<max_dec)
            {
                for(int k_sw2=0;k_sw2<1000;k_sw2++)
                {
                    ds=ds_upBound_aLmt[swCount][legID]+0.001*0.001*k_sw2;
                    max_dec=GetSwingMaxDec(count,ds,legID);
                    min_acc=GetSwingMinAcc(count,ds,legID);

                    if(min_acc<max_dec)
                    {
                        ds_upBoundFlag_sw=true;
                        break;
                    }
                    else
                    {
                        dds_lowBound[swCount][legID]=max_dec;
                        dds_upBound[swCount][legID]=min_acc;
                        ds_upBound_aLmt[swCount][legID]=ds;
                    }
                }
            }
        }
    }

    double vLmt_value[3] {0};
    double Jvi_dot_vb[3] {0};
    double vb_sw3[3] {0};
    vb_sw3[2]=pva_b[count][1];
    aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,vb_sw3,1,1,Jvi_dot_vb,1);
    for (int k=0;k<3;k++)
    {
        if(output_dds[count][3*legID+k]>0)
        {
            vLmt_value[k]=(Jvi_dot_vb[k]+vLmt)/output_dds[count][3*legID+k];
        }
        else if(output_dds[count][3*legID+k]<0)
        {
            vLmt_value[k]=(Jvi_dot_vb[k]-vLmt)/output_dds[count][3*legID+k];
        }
        else
        {
            printf("WARNING!!! param_dds equals zero!!! SwingLeg : %d \n",count);
        }

    }
    ds_upBound_vLmt[swCount][legID]=*std::min_element(vLmt_value,vLmt_value+3);
    ds_upBound[swCount][legID]=std::min(ds_upBound_aLmt[swCount][legID],ds_upBound_vLmt[swCount][legID]);
    ds_lowBound[swCount][legID]=ds_lowBound_aLmt[swCount][legID];
}

void TimeOptimalGait::GetSwingSwitchPoint(int legID)
{
    double slopedsBound[901] {0};
    double paramdds0Point[901] {0};
    double tangentPoint[901] {0};
    int paramdds0Count {0};
    int tangentCount {0};
    int switchScrewID_tmp[901] {0};

    //initialize
    paramdds0Count=0;
    tangentCount=0;
    for(int i=0;i<901;i++)
    {
        tangentPoint[i]=-1;
        switchScrewID_tmp[i]=-1;
        paramdds0Point[i]=-1;
    }

    slopedsBound[0]=(ds_upBound[1][legID]-ds_upBound[0][legID])/(s_w[1][legID]-s_w[0][legID]);
    for(int i=1;i<900;i++)
    {
        slopedsBound[i]=( (ds_upBound[i+1][legID]-ds_upBound[i][legID])/(s_w[i+1][legID]-s_w[i][legID])
                         +(ds_upBound[i][legID]-ds_upBound[i-1][legID])/(s_w[i][legID]-s_w[i-1][legID]) )/2;
    }
    slopedsBound[900]=(ds_upBound[900][legID]-ds_upBound[899][legID])/(s_w[900][legID]-s_w[899][legID]);

    for(int i=0;i<901;i++)
    {
        slopeDelta[i][legID]=slopedsBound[i]-(dds_upBound[i][legID]+dds_lowBound[i][legID])/2;
    }

    for(int i=0;i<900;i++)
    {
        bool isParamdds0 {false};
        if(ds_upBound_vLmt[i][legID]>ds_upBound_aLmt[i][legID])
        {
            int count=(legID%2==0) ? i : (i+1350);
            for(int k=0;k<3;k++)
            {
                if(output_dds[count+1][3*legID+k]*output_dds[count][3*legID+k]<0 || output_dds[count][3*legID+k]==0)
                {
                    paramdds0Point[paramdds0Count]=i
                            +fabs(output_dds[i][3*legID+k])/(fabs(output_dds[i][3*legID+k])+fabs(output_dds[i+1][3*legID+k]));
                    switchScrewID_tmp[paramdds0Count]=k;
                    paramdds0Count++;
                    isParamdds0=true;
                }
            }
        }

        if((slopeDelta[i+1][legID]*slopeDelta[i][legID]<0 || slopeDelta[i][legID]==0) && isParamdds0==false)
        {
            tangentPoint[tangentCount]=i
                    +fabs(slopeDelta[i][legID])/(fabs(slopeDelta[i][legID])+fabs(slopeDelta[i+1][legID]));
            tangentCount++;
            if((fabs(slopeDelta[i][legID])+fabs(slopeDelta[i+1][legID]))==0)
            {
                printf("WARNING! tangentPoint!!!\n");
            }
        }
    }

    printf("SwingLeg %d Tangent Switch Point:",legID);
    for(int i=0;i<tangentCount+1;i++)
    {
        printf("%.1f,",tangentPoint[i]);
    }
    printf("\n");

    printf("SwingLeg %d ZeroInertia Switch Point:",legID);
    for(int i=0;i<paramdds0Count+1;i++)
    {
        printf("%.1f,",paramdds0Point[i]);
    }
    printf("\n");

    //merge tangentPoint & paramdds0Point into switchPoint
    switchPoint[0][legID]=0;
    switchPoint[1][legID]=900;
    switchCount[legID]=2;
    for(int i=0;i<tangentCount;i++)
    {
        switchPoint[i+switchCount[legID]][legID]=tangentPoint[i];
    }
    switchCount[legID]+=tangentCount;
    for(int i=0;i<paramdds0Count;i++)
    {
        switchPoint[i+switchCount[legID]][legID]=paramdds0Point[i];
        switchScrewID[i+switchCount[legID]][legID]=switchScrewID_tmp[i];
    }
    switchCount[legID]+=paramdds0Count;
//    for(int i=0;i<switchCount_body;i++)
//    {
//        if(switchPoint_body[i]<=900 && legID%2==1)//135
//        {
//            switchPoint[switchCount[legID]][legID]=switchPoint_body[i];
//            switchCount[legID]++;
//        }
//        if(switchPoint_body[i]>1350 && legID%2==0)
//        {
//            switchPoint[switchCount[legID]][legID]=switchPoint_body[i];
//            switchCount[legID]++;
//        }
//    }

    //filtering the same point & sorting by the value
    for(int i=0;i<switchCount[legID];i++)
    {
        for(int k=i+1;k<switchCount[legID];k++)
        {
            if((int)switchPoint[k][legID]<(int)switchPoint[i][legID])
            {
                auto tmp1=switchPoint[i][legID];
                switchPoint[i][legID]=switchPoint[k][legID];
                switchPoint[k][legID]=tmp1;

                auto tmp2=switchScrewID[i][legID];
                switchScrewID[k][legID]=switchScrewID[k][legID];
                switchScrewID[k][legID]=tmp2;
            }
//            else if((int)switchPoint[k][legID]==(int)switchPoint[i][legID])
//            {
//                switchPoint[k][legID]=switchPoint[switchCount[legID]-1][legID];
//                switchPoint[switchCount[legID]-1][legID]=-1;
//                k--;
//                switchCount[legID]--;
//            }
        }
    }

    printf("SwingLeg %d Switch Point:",legID);
    for(int i=0;i<switchCount[legID]+1;i++)
    {
        printf("%.1f,",switchPoint[i][legID]);
    }
    printf("\n");
}

double TimeOptimalGait::GetSwingSwitchMaxDec(int switchID, double ds, int legID)
{
    double dec[3] {0};
    std::fill_n(dec,3,-1e6);
    for (int k=0;k<3;k++)
    {
        if(k!=switchScrewID[switchID][legID])
        {
            dec[k]=param_a2[3*legID+k]*ds*ds+param_a1[3*legID+k]+param_a0L[3*legID+k];
        }
    }
    return *std::max_element(dec,dec+3);
}

double TimeOptimalGait::GetSwingSwitchMinAcc(int switchID, double ds, int legID)
{
    double acc[3] {0};
    std::fill_n(acc,3,1e6);
    for (int k=0;k<3;k++)
    {
        if(k!=switchScrewID[switchID][legID])
        {
            acc[k]=param_a2[3*legID+k]*ds*ds+param_a1[3*legID+k]+param_a0H[3*legID+k];
        }
    }
    return *std::min_element(acc,acc+3);
}

double TimeOptimalGait::GetSwingSwitchDsBound(int switchID, int legID, double *pva_body)
{
    bool ds_lowBoundFlag_sw {false};
    bool ds_upBoundFlag_sw {false};
    int k_sw {0};
    const int kswCount {30000};
    double ds_a;
    while (ds_upBoundFlag_sw==false)
    {
        ds_a=0.001*k_sw;
        double max_dec=GetSwingSwitchMaxDec(switchID,ds_a,legID);
        double min_acc=GetSwingSwitchMinAcc(switchID,ds_a,legID);

        k_sw++;
        if(k_sw==kswCount)
        {
            ds_upBoundFlag_sw=true;
            printf("WARNING!!! kswCount=%d is too small!!! Leg:%d, switchPoint=%.1f, switchScrewID:%d\n",kswCount,legID,switchPoint[switchID][legID],switchScrewID[switchID][legID]);
            for(int k=0;k<3;k++)
            {
            printf("paramdds=%.1f,parama2=%.1f,parama1=%.1f,parama0H=%.1f,parama0L=%.1f\n",
                   param_dds[3*legID+k],param_a2[3*legID+k],param_a1[3*legID+k],param_a0H[3*legID+k],param_a0L[3*legID+k]);
            }
        }
        else
        {
            if(k_sw%1000==0)
            {
                //printf("min_acc=%.1f, max_dec=%.1f\n",min_acc,max_dec);
            }
            if(ds_lowBoundFlag_sw==false && ds_upBoundFlag_sw==false && min_acc>=max_dec)
            {
                ds_lowBoundFlag_sw=true;
            }
            else if(ds_lowBoundFlag_sw==true && ds_upBoundFlag_sw==false && min_acc<max_dec)
            {
                k_sw--;
                for(int k_sw2=0;k_sw2<1000;k_sw2++)
                {
                    ds_a=0.001*k_sw+0.001*0.001*k_sw2;
                    max_dec=GetSwingSwitchMaxDec(switchID,ds_a,legID);
                    min_acc=GetSwingSwitchMinAcc(switchID,ds_a,legID);

                    if(min_acc<max_dec)
                    {
                        ds_a=0.001*k_sw+0.001*0.001*(k_sw2-1);
                        ds_upBoundFlag_sw=true;
                        break;
                    }
                }
            }
        }
    }

    double vLmt_value[3] {0};
    double Jvi_dot_vb[3] {0};
    double vb_sw3[3] {0};
    vb_sw3[2]=pva_body[1];
    aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,vb_sw3,1,1,Jvi_dot_vb,1);
    for (int k=0;k<3;k++)
    {
        if(k!=switchScrewID[switchID][legID])
        {
            if(param_dds[3*legID+k]>0)
            {
                vLmt_value[k]=(Jvi_dot_vb[k]+vLmt)/param_dds[3*legID+k];
            }
            else if(param_dds[3*legID+k]<0)
            {
                vLmt_value[k]=(Jvi_dot_vb[k]-vLmt)/param_dds[3*legID+k];
            }
        }
        vLmt_value[switchScrewID[switchID][legID]]=1e6;
    }
    double ds_v=*std::min_element(vLmt_value,vLmt_value+3);
    return std::min(ds_a,ds_v);
}

void TimeOptimalGait::GetSwingTwoPointAtSwitch(int legID, double *lowPoint, double *upPoint)
{
    lowPoint[0]=-1;
    lowPoint[switchCount[legID]-1]=0;
    upPoint[0]=0;
    upPoint[switchCount[legID]-1]=-1;

    for(int i=1;i<switchCount[legID]-1;i++)
    {
        int num=round(switchPoint[i][legID]);
//        int swCount=(int)switchPoint[i][legID];
//        int count = legID%2==0 ? swCount : (swCount+1350);
//        double alpha=switchPoint[i][legID]-swCount;
//        double sb = s_b[count]+alpha*(s_b[count+1]-s_b[count]);
//        double sw = s_w[swCount][legID]+alpha*(s_w[swCount+1][legID]-s_w[swCount][legID]);
//        double pva_b_switch[3];
//        for(int k=0;k<3;k++)
//        {
//            pva_b_switch[k]=pva_b[count][k]+alpha*(pva_b[count+1][k]-pva_b[count][k]);
//        }

//        GetSwingLegParam(-1,legID,sw,pva_b_switch);

//        double ds1=GetSwingSwitchDsBound(i,legID,pva_b_switch);
//        double ds_tmp=std::min(ds_upBound[num-1][legID],ds_upBound[num][legID]);
//        double ds2=std::min(ds_tmp,ds_upBound[num+1][legID]);

//        double dds_back=GetSwingMaxDec(i,ds2,legID);
//        double dds_for=GetSwingMinAcc(i,ds2,legID);
//        if(dds_for<dds_back)
//        {
//            printf("\nWARNING!!!dds_for < dds_back at switchPoint=%.1f,dds_back=%.1f,dds_for=%.1f,ds1=%.8f,ds2=%.8f\n",switchPoint[i][legID],dds_back,dds_for,ds1,ds2);
//            for(int k=0;k<3;k++)
//            {
//                printf("param_dds=%.4f,param_a2=%.4f,param_a1=%.4f,param_a0L=%.4f,param_a0H=%.4f\n\n",param_dds[3*legID+k],param_a2[3*legID+k],param_a1[3*legID+k],param_a0L[3*legID+k],param_a0H[3*legID+k]);
//            }
//        }

//        int num=(int)switchPoint[i][legID];
        if(slopeDelta[num][legID]==0)
        {
            printf("Amazing!!! double equals 0\n");
            upPoint[i]=lowPoint[i]=ds_upBound[num][legID];
        }
        else
        {
//            lowPoint[i]=std::min(ds_upBound[num-1][legID],sqrt(ds*ds-2*dds_back*(switchPoint[i][legID]+1-num)*(s_w[num][legID]-s_w[num-1][legID])));
//            upPoint[i]=std::min(ds_upBound[num+1][legID],sqrt(ds*ds+2*dds_for*(num+1-switchPoint[i][legID])*(s_w[num+1][legID]-s_w[num][legID])));
            auto tmp=std::min(ds_upBound[num-1][legID],ds_upBound[num][legID]);
            upPoint[i]=lowPoint[i]=std::min(tmp,ds_upBound[num+1][legID]);//std::min(ds_upBound[num][legID],ds_upBound[num+1][legID]);
        }
    }
}

void TimeOptimalGait::GetSwingOptimalDsBySwitchPoint(int legID)
{
    bool stopFlag {false};
    double *lowPoint=new double [switchCount[legID]];
    double *upPoint=new double [switchCount[legID]];
    GetSwingTwoPointAtSwitch(legID,lowPoint,upPoint);

    printf("lowPoint:");
    for(int i=0;i<switchCount[legID];i++)
    {
        printf("%.2f,",lowPoint[i]);
    }
    printf("\n");
    printf("upPoint:");
    for(int i=0;i<switchCount[legID];i++)
    {
        printf("%.2f,",upPoint[i]);
    }
    printf("\n");


    for(int i=0;i<901;i++)
    {
        real_ds[i][legID]=ds_upBound[i][legID];
        real_dds[i][legID]=dds_upBound[i][legID];
    }
    for(int m=0;m<switchCount[legID];m++)
    {
        int k_sw=round(switchPoint[m][legID]);
        real_dds[k_sw][legID]=0;

        if(ds_upBound[k_sw][legID]>real_ds[k_sw][legID] && k_sw!=900)//start of backward
        {
            printf("SwingLeg backward start at a passed point, quit switchPoint %.1f\n",switchPoint[m][legID]);
            continue;
        }
        stopFlag=false;
        if(k_sw!=0 && k_sw!=900)
        {
            k_sw--;
        }
        if(k_sw!=0)
        {
            ds_backward[k_sw][legID]=lowPoint[m];
            while(stopFlag==false)
            {
                int count=(legID%2==0) ? k_sw : (k_sw+1350);
                dds_backward[k_sw][legID]=GetSwingMaxDec(count,ds_backward[k_sw][legID],legID);
                ds_backward[k_sw-1][legID]=sqrt(ds_backward[k_sw][legID]*ds_backward[k_sw][legID]-2*dds_backward[k_sw][legID]*(s_w[k_sw][legID]-s_w[k_sw-1][legID]));

                if(ds_backward[k_sw-1][legID]>ds_upBound[k_sw-1][legID])
                {
                    stopFlag=true;
                    printf("SwingLeg backward touching upBound at %d, quit switchPoint %.4f\n",k_sw-1,switchPoint[m][legID]);
                }
                else if(k_sw==1)
                {
                    dds_backward[k_sw-1][legID]=GetSwingMaxDec(count-1,ds_backward[k_sw-1][legID],legID);
                    for(int i=k_sw-1;i<round(switchPoint[m][legID]);i++)
                    {
                        real_ds[i][legID]=ds_backward[i][legID];
                        real_dds[i][legID]=dds_backward[i][legID];
                    }
                    stopFlag=true;
                    printf("SwingLeg backward touching 0, from switchPoint %.4f\n",switchPoint[m][legID]);
                }
                else if(ds_backward[k_sw-1][legID]>=real_ds[k_sw-1][legID])
                {
                    real_dds[k_sw-1][legID]=(real_ds[k_sw-1][legID]-ds_backward[k_sw][legID])*(real_ds[k_sw-1][legID]+ds_backward[k_sw][legID])/2/(s_w[k_sw][legID]-s_w[k_sw-1][legID]);
                    for(int i=k_sw;i<round(switchPoint[m][legID]);i++)
                    {
                        real_ds[i][legID]=ds_backward[i][legID];
                        real_dds[i][legID]=dds_backward[i][legID];
                    }
                    stopFlag=true;
                    printf("SwingLeg backward touching last curve at %d, from switchPoint %.4f\n",k_sw-1,switchPoint[m][legID]);
                }
                else
                {
                    k_sw--;
                }
            }
            if(ds_backward[k_sw-1][legID]>ds_upBound[k_sw-1][legID] || round(switchPoint[m][legID])==900)
            {
                continue;
            }
        }

        k_sw=round(switchPoint[m][legID]);
        if(k_sw!=0)
        {
            k_sw++;
        }
        stopFlag=false;
        ds_forward[k_sw][legID]=upPoint[m];
        while(stopFlag==false)
        {
            int count=(legID%2==0) ? k_sw : (k_sw+1350);
            dds_forward[k_sw][legID]=GetSwingMinAcc(count,ds_forward[k_sw][legID],legID);
            ds_forward[k_sw+1][legID]=sqrt(ds_forward[k_sw][legID]*ds_forward[k_sw][legID]+2*dds_forward[k_sw][legID]*(s_w[k_sw+1][legID]-s_w[k_sw][legID]));

            if(ds_forward[k_sw+1][legID]>ds_upBound[k_sw+1][legID] || k_sw==899)
            {
                for(int i=round(switchPoint[m][legID])+1;i<k_sw+1;i++)
                {
                    real_ds[i][legID]=ds_forward[i][legID];
                    real_dds[i][legID]=dds_forward[i][legID];
                }
                if(round(switchPoint[m][legID])==0)
                {
                    real_ds[0][legID]=ds_forward[0][legID];
                    real_dds[0][legID]=dds_forward[0][legID];
                }
                if(k_sw==899)
                {
                    if(ds_forward[k_sw+1][legID]>ds_upBound[k_sw+1][legID])
                    {
                        real_ds[900][legID]=ds_upBound[900][legID];
                        real_dds[900][legID]=dds_upBound[900][legID];
                    }
                    else
                    {
                        real_ds[900][legID]=ds_forward[900][legID];
                        real_dds[900][legID]=GetSwingMinAcc(900,ds_forward[900][legID],legID);
                    }
                }
                stopFlag=true;
                printf("SwingLeg forward touching upBound at %d, from switchPoint %.4f\n",k_sw,switchPoint[m][legID]);
            }
            else
            {
                k_sw++;
            }
        }
    }
    delete [] lowPoint;
    delete [] upPoint;
}

void TimeOptimalGait::GetSwingOptimalDsByDirectNI(int legID)
{
    bool stopFlag {false};
    bool accFlag {true};
    int dec_start {0};
    int dec_end {0};
    unsigned int cycleCount {0};
    int stop_back {0};
    int ki_back {0};
    int ki_for {0};

    //backward integration
    ki_back=900;
    ds_backward[ki_back][legID]=0;
    while (stopFlag==false && ki_back>=0)
    {
        int count=(legID%2==0) ? ki_back : (ki_back+1350);
        dds_backward[ki_back][legID]=GetSwingMaxDec(count,ds_backward[ki_back][legID],legID);
        ds_backward[ki_back-1][legID]=sqrt(ds_backward[ki_back][legID]*ds_backward[ki_back][legID]-2*dds_backward[ki_back][legID]*(s_w[ki_back][legID]-s_w[ki_back-1][legID]));

        if (ds_backward[ki_back-1][legID]>ds_upBound[ki_back-1][legID])
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
    ki_for=0;
    ds_forward[ki_for][legID]=0;
    double min_dist[901] {0};
    std::fill_n(min_dist,901,1);
    while (stopFlag==false)
    {
        if (accFlag==true)
        {
            int count=(legID%2==0) ? ki_for : (ki_for+1350);
            dds_forward[ki_for][legID]=GetSwingMinAcc(count,ds_forward[ki_for][legID],legID);
            ds_forward[ki_for+1][legID]=sqrt(ds_forward[ki_for][legID]*ds_forward[ki_for][legID]+2*dds_forward[ki_for][legID]*(s_w[ki_for+1][legID]-s_w[ki_for][legID]));

            if (ds_forward[ki_for+1][legID]>ds_upBound[ki_for+1][legID])
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
            int count=(legID%2==0) ? ki_for : (ki_for+1350);
            dds_forward[ki_for][legID]=GetSwingMaxDec(count,ds_forward[ki_for][legID],legID);
            ds_forward[ki_for+1][legID]=sqrt(ds_forward[ki_for][legID]*ds_forward[ki_for][legID]+2*dds_forward[ki_for][legID]*(s_w[ki_for+1][legID]-s_w[ki_for][legID]));

            if (ds_forward[ki_for+1][legID]>ds_upBound[ki_for+1][legID])
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
                    ds_forward[0][legID]-=ds_upBound[0][legID]/1000;
                }
            }
            else
            {
                //printf("dec ending,ki_for=%d, ds_forward=%.4f\n",ki_for,ds_forward[ki_for+1][j+1]);
                if ((ds_forward[ki_for][legID]*ds_forward[ki_for][legID]+2*dds_forward[ki_for][legID]*(s_w[ki_for+1][legID]-s_w[ki_for][legID]))
                        <=(ds_lowBound[ki_for+1][legID]*ds_lowBound[ki_for+1][legID]) || ki_for==899)
                {
                    accFlag=true;
                    for(int k=dec_start;k<(ki_for+2);k++)
                    {
                        min_dist[k]=ds_upBound[k][legID]-ds_forward[k][legID];
                    }
                    dec_end=std::min_element(min_dist+dec_start+1,min_dist+ki_for+2)-min_dist;
                    //dec_start must be ignored, if dec_start is the min_dist, the calculation will cycle between dec_start & dec_start+1
                    ki_for=dec_end-1;
                    //printf("dec finished, start at k=%d, end at k=%d, ds_forward=%.4f\n",dec_start,dec_end,ds_forward[ki_for+1][j+1]);
                }
                ki_for++;
            }
        }

        if(ki_for==900)
        {
            stopFlag=true;
            for(int i=0;i<901;i++)
            {
                real_ds[i][legID]=ds_forward[i][legID];
                real_dds[i][legID]=dds_forward[i][legID];
            }
            //printf("SwingLeg forward reach the end, and never encounter with the backward, ki=%u < %u\n",cycleCount,0x0FFFFFFF);
        }
        if(ki_for>=stop_back && ds_forward[ki_for][legID]>=ds_backward[ki_for][legID])
        {
            stopFlag=true;
            for(int i=0;i<ki_for;i++)
            {
                real_ds[i][legID]=ds_forward[i][legID];
                real_dds[i][legID]=dds_forward[i][legID];
            }
            for(int i=ki_for;i<901;i++)
            {
                real_ds[i][legID]=ds_backward[i][legID];
                real_dds[i][legID]=dds_backward[i][legID];
            }
            //printf("SwingLeg forward & backward encounters at k=%d\n",ki_for);
        }
        cycleCount++;
        if(cycleCount==0x0FFFFFFF)
        {
            stopFlag=true;
            printf("WARNING!!! SwingLeg integration takes too long, force stop!!! ki=%d\n",cycleCount);
        }
    }
}

//void TimeOptimalGait::GetPvabAndSw(int i)
//{
//    pva_b[i][0]=b_sb[i];
//    pva_b[i][1]=db_sb[i]*real_ds_body[i];
//    pva_b[i][2]=ddb_sb[i]*real_ds_body[i]*real_ds_body[i]+db_sb[i]*real_dds_body[i];
//}

void TimeOptimalGait::GetOptimalDs()
{
    rbt.loadXml("/home/hex/Desktop/mygit/RobotVIII_demo/resource/Robot_VIII.xml");

    timeval tpstart,tpend;
    float tused;
    gettimeofday(&tpstart,NULL);

    GetStanceOptimalDsByMinorIteration();

    //pva_b & s_w initialized here, and need to be updated during major iteration
    for (int i=0;i<2251;i++)
    {
        output_init_pva_b[i][0]=pva_b[i][0]=b_sb[i];
        output_init_pva_b[i][1]=pva_b[i][1]=db_sb[i]*real_ds_body[i];
        output_init_pva_b[i][2]=pva_b[i][2]=ddb_sb[i]*real_ds_body[i]*real_ds_body[i]+db_sb[i]*real_dds_body[i];

        for(int j=0;j<3;j++)
        {
            if(i<901)
            {
                s_w[i][2*j]=PI/s_b1*s_b[i];
            }
            if(i>=1350 && i<2251)
            {
                s_w[i-1350][2*j+1]=PI/(1-s_b2)*(s_b[i]-s_b2);
            }
        }
    }


    /******************************* SwingLeg & Iteration ************************************/
    int iterCount {0};
    bool stopFlag {false};
    double Tstep_tmp {Tstep};
    double totalTime[6] {0};
    double totalTime_last[6] {0};
    while(stopFlag==false && iterCount<=100)
    {
        printf("\n");
        for(int i=0;i<901;i++)
        {
            for(int j=0;j<6;j++)
            {
                ds_backward[i][j]=0;
                ds_forward[i][j]=0;
            }
        }
        iterCount++;

        //generate the traj & calculate the bound of ds
        std::fill_n(totalTime,6,0);
        std::fill_n(switchCount,6,0);
        std::fill_n(*switchPoint,901*6,-1);
        std::fill_n(*switchScrewID,901*6,-1);
        std::fill_n(*isParamddsExact0,901*18,-1);

        for (int j=0;j<6;j++)
        {
            for (int i=0;i<901;i++)//024
            {
                int count = (j%2==0) ? i : (i+1350);
                GetSwingLegParam(count,j,s_w[i][j],*pva_b+3*count);
                GetSwingDsBound(count,j);
            }

            GetSwingSwitchPoint(j);
            GetSwingOptimalDsBySwitchPoint(j);
            //GetSwingOptimalDsByDirectNI(j);

            for (int i=1;i<901;i++)
            {
                totalTime[j]+=2*(s_w[i][j]-s_w[i-1][j])/(real_ds[i-1][j]+real_ds[i][j]);
                timeArray[i][j]=totalTime[j];
            }
        }
        maxTime=*std::max_element(totalTime,totalTime+6);
        maxTotalCount=(int)(maxTime*1000)+1;
        maxTime=0.001*maxTotalCount;
        printf("totalTime: %.4f, %.4f, %.4f, %.4f, %.4f, %.4f; maxTime:%.4f\n",
               totalTime[0],totalTime[1],totalTime[2],totalTime[3],totalTime[4],totalTime[5],maxTime);

        stopFlag=true;
        for(int i=0;i<6;i++)
        {
            //if((int)(totalTime_last[i]*1000)!=(int)(totalTime[i]*1000))
            if(fabs(totalTime_last[i]-totalTime[i])>1e-4)
            {
                stopFlag=false;
            }
        }

        //update pva_b & s_w here
        double timeArray_body_tmp[2251] {0};
        double timeArray_tmp[901][6] {0};
        double pva_b_tmp[2251][3] {0};
        for(int i=0;i<2251;i++)
        {
            timeArray_body_tmp[i]=timeArray_body[i]*maxTime/((1-dutyCycle)*Tstep);
            pva_b_tmp[i][0]=pva_b[i][0];
            pva_b_tmp[i][1]=pva_b[i][1]*(1-dutyCycle)*Tstep/maxTime;
            pva_b_tmp[i][2]=pva_b[i][2]*(1-dutyCycle)*Tstep/maxTime*(1-dutyCycle)*Tstep/maxTime;
        }
        Tstep_tmp=timeArray_body_tmp[2250];

        for(int i=0;i<901;i++)
        {
            for(int j=0;j<6;j++)
            {
                if(totalTime[j]!=0)
                {
                    timeArray_tmp[i][j]=timeArray[i][j]*maxTime/totalTime[j];
                }
            }
        }

        int k_start {0};
        double s_w_tmp[901][6] {0};
        for(int j=0;j<3;j++)
        {
            k_start=0;
            for(int i=0;i<900;i++)
            {
                for(int k=k_start;k<900;k++)
                {
                    if(timeArray_tmp[k][2*j]<=timeArray_body_tmp[i] && timeArray_tmp[k+1][2*j]>timeArray_body_tmp[i])
                    {
                        k_start=k;
                        s_w_tmp[i][2*j]=s_w[k][2*j]+(s_w[k+1][2*j]-s_w[k][2*j])
                                *(timeArray_body_tmp[i]-timeArray_tmp[k][2*j])/(timeArray_tmp[k+1][2*j]-timeArray_tmp[k][2*j]);
                        break;
                    }
                }
            }
        }
        for(int j=0;j<3;j++)
        {
            k_start=0;
            for(int i=1350;i<2250;i++)
            {
                for(int k=k_start;k<900;k++)
                {
                    if(timeArray_tmp[k][2*j+1]<=(timeArray_body_tmp[i]-dutyCycle*Tstep_tmp)
                            && timeArray_tmp[k+1][2*j+1]>(timeArray_body_tmp[i]-dutyCycle*Tstep_tmp))
                    {
                        k_start=k;
                        s_w_tmp[i-1350][2*j+1]=s_w[k][2*j+1]+(s_w[k+1][2*j+1]-s_w[k][2*j+1])
                                *(timeArray_body_tmp[i]-dutyCycle*Tstep_tmp-timeArray_tmp[k][2*j+1])/(timeArray_tmp[k+1][2*j+1]-timeArray_tmp[k][2*j+1]);
                        break;
                    }
                }
            }
        }
        for(int i=0;i<6;i++)
        {
            s_w_tmp[900][i]=PI;
        }

        memcpy(*s_w,*s_w_tmp,901*6*sizeof(double));
        memcpy(*pva_b,*pva_b_tmp,2251*3*sizeof(double));
        memcpy(totalTime_last,totalTime,6*sizeof(double));
    }

    printf("Iteration finished, iterCount=%d, maxTime=%.4f, Tstep_tmp:%.4f\n\n",iterCount,maxTime,Tstep_tmp);

    gettimeofday(&tpend,NULL);
    tused=tpend.tv_sec-tpstart.tv_sec+(double)(tpend.tv_usec-tpstart.tv_usec)/1000000;
    printf("UsedTime:%f\n",tused);
}

void TimeOptimalGait::GetOptimalGait2s()
{
//    double pEB[6] {0};
//    double pEE[18] {0};
//    double vEB[6] {0};
//    double vEE[18] {0};
//    double aEB[6] {0};
//    double aEE[18] {0};

//    for (int i=0;i<901;i++)
//    {
//        s_w=0.1*i*PI/180;//degree to rad

//        f_sw[0]=0;
//        f_sw[1]=stepH*sin(s_w);
//        f_sw[2]=stepD/2*cos(s_w);
//        df_sw[0]=0;
//        df_sw[1]=stepH*cos(s_w);
//        df_sw[2]=-stepD/2*sin(s_w);
//        ddf_sw[0]=0;
//        ddf_sw[1]=-stepH*sin(s_w);
//        ddf_sw[2]=-stepD/2*cos(s_w);

//        pEB[2]=initPeb[2]+pb_sw_tmp[i];
//        vEB[2]=vb_sw_tmp[i];
//        aEB[2]=ab_sw_tmp[i];
//        for(int j=0;j<3;j++)
//        {
//            //stanceLeg
//            pEE[6*j+3]=initPee[6*j+3];
//            pEE[6*j+4]=initPee[6*j+4];
//            pEE[6*j+5]=initPee[6*j+5];

//            //swingLeg
//            pEE[6*j]=initPee[6*j]+f_sw[0];
//            pEE[6*j+1]=initPee[6*j+1]+f_sw[1];
//            pEE[6*j+2]=initPee[6*j+2]+f_sw[2];

//            vEE[6*j]=df_sw[0]*real_ds[i][2*j];
//            vEE[6*j+1]=df_sw[1]*real_ds[i][2*j];
//            vEE[6*j+2]=df_sw[2]*real_ds[i][2*j];

//            aEE[6*j]=ddf_sw[0]*real_ds[i][2*j]*real_ds[i][2*j]+df_sw[0]*real_dds[i][2*j];
//            aEE[6*j+1]=ddf_sw[1]*real_ds[i][2*j]*real_ds[i][2*j]+df_sw[1]*real_dds[i][2*j];
//            aEE[6*j+2]=ddf_sw[2]*real_ds[i][2*j]*real_ds[i][2*j]+df_sw[2]*real_dds[i][2*j];
//        }

//        rbt.SetPeb(pEB);
//        rbt.SetVb(vEB);
//        rbt.SetAb(aEB);
//        rbt.SetPee(pEE);
//        rbt.SetVee(vEE);
//        rbt.SetAee(aEE);

//        rbt.GetPee(*output_Pee+18*i,rbt.body());
//        rbt.GetPin(*output_Pin+18*i);
//        rbt.GetVin(*output_Vin+18*i);
//        rbt.GetAin(*output_Ain+18*i);

//        //printf("output:Pee=%.4f,Pin=%.4f,Vin=%.4f,Ain=%.4f\n",output_Pee[i][0],output_Pin[i][0],output_Vin[i][0],output_Ain[i][0]);
//    }
}

void TimeOptimalGait::GetOptimalGait2t()
{
//    double * real_s=new double [6*maxTotalCount];
//    double * real_Pee=new double [18*maxTotalCount];
//    double * real_Pin=new double [18*maxTotalCount];

//    for(int j=0;j<6;j++)
//    {
//        real_s[j]=0;
//        for(int i=0;i<901;i++)
//        {
//            real_ds_tmp[i][j]=real_ds[i][j]*totalTime[j]/maxTime;
//        }
//    }

//    for (int i=1;i<maxTotalCount;i++)
//    {
//        double ds[6];
//        for(int j=0;j<6;j++)
//        {
//            ds[j]=0.5*(real_ds_tmp[(int)(real_s[6*(i-1)+j]/(PI/900))][j]+real_ds_tmp[(int)(real_s[6*(i-1)+j]/(PI/900))+1][j]);
//            real_s[6*i+j]=real_s[6*(i-1)+j]+ds[j]*0.001;
//        }
//    }

//    for (int i=0;i<maxTotalCount;i++)
//    {
//        double pEB[6];
//        memcpy(pEB,initPeb,6*sizeof(double));
//        double pEE[18];
//        double b_s=stepD/4-stepD/2*(real_s[6*i+stanceLegID[0]]/PI);
//        pEB[2]=initPeb[2]+b_s;
//        for(int j=0;j<3;j++)
//        {
//            //stanceLeg
//            pEE[6*j+3]=initPee[6*j+3];
//            pEE[6*j+4]=initPee[6*j+4];
//            pEE[6*j+5]=initPee[6*j+5];

//            //swingLeg
//            double f_s[3];
//            f_s[0]=0;
//            f_s[1]=stepH*sin(real_s[6*i+2*j]);
//            f_s[2]=stepD/2*cos(real_s[6*i+2*j]);
//            pEE[6*j]=initPee[6*j]+f_s[0];
//            pEE[6*j+1]=initPee[6*j+1]+f_s[1];
//            pEE[6*j+2]=initPee[6*j+2]+f_s[2];
//        }
//        rbt.SetPeb(pEB);
//        rbt.SetPee(pEE);

//        rbt.GetPee(real_Pee+18*i,rbt.body());
//        rbt.GetPin(real_Pin+18*i);
//    }

//    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_s.txt",real_s,maxTotalCount,6);
//    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_Pee.txt",real_Pee,maxTotalCount,18);
//    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_Pin.txt",real_Pin,maxTotalCount,18);
//    delete [] real_s;
//    delete [] real_Pee;
//    delete [] real_Pin;
}

void TimeOptimalGait::OutputData()
{
    printf("Start output data...\n");
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_upBound_aLmt_body.txt",ds_upBound_aLmt_body,2251,1);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_upBound_vLmt_body.txt",ds_upBound_vLmt_body,2251,1);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_upBound_body.txt",dds_upBound_body,2251,1);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_upBound_body.txt",dds_lowBound_body,2251,1);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_forward_body.txt",ds_forward_body,2251,1);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_backward_body.txt",ds_backward_body,2251,1);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_forward_body.txt",dds_forward_body,2251,1);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_backward_body.txt",dds_backward_body,2251,1);

    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_PeeB.txt",*output_PeeB,2251,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_dsds.txt", *output_dsds, 2251,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_ds.txt", *output_ds, 2251,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_const.txt", *output_const, 2251,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_const1.txt", *output_const1, 2251,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_const2.txt", *output_const2, 2251,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_dds.txt",*output_dds,2251,18);

    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a2.txt",*output_a2,2251,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a1.txt",*output_a1,2251,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a0L.txt",*output_a0L,2251,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a0H.txt",*output_a0H,2251,18);

    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/max_ValueL.txt",*dds_lowBound,901,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/min_ValueH.txt",*dds_upBound,901,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_upBound_aLmt.txt",*ds_upBound_aLmt,901,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_lowBound_aLmt.txt",*ds_lowBound_aLmt,901,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_upBound_vLmt.txt",*ds_upBound_vLmt,901,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_lowBound.txt",*ds_lowBound,901,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_upBound.txt",*ds_upBound,901,6);
    //aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_max.txt",*real_ddsMax,901,6);
    //aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_min.txt",*real_ddsMin,901,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_forward.txt",*ds_forward,901,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_backward.txt",*ds_backward,901,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_forward.txt",*dds_forward,901,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_backward.txt",*dds_backward,901,6);

    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_ds_body.txt",real_ds_body,2251,1);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_dds_body.txt",real_dds_body,2251,1);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_ds.txt",*real_ds,901,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_dds.txt",*real_dds,901,6);

    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/s_w.txt",*s_w,901,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/pva_b.txt",*pva_b,2251,3);

//    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/timeArray.txt",*timeArray_tmp,901,6);
//    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/pb_sw.txt",pb_sw_tmp,901,1);
//    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/vb_sw.txt",vb_sw_tmp,901,1);
//    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ab_sw.txt",ab_sw_tmp,901,1);

//    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Pee.txt",*output_Pee,901,18);
//    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Pin.txt",*output_Pin,901,18);
//    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Vin.txt",*output_Vin,901,18);
//    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Ain.txt",*output_Ain,901,18);
    printf("Finish output data\n");
}
