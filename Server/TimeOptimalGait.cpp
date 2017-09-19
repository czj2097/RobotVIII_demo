#include "TimeOptimalGait.h"

TimeOptimalGait::TimeOptimalGait(){}
TimeOptimalGait::~TimeOptimalGait(){}

void TimeOptimalGait::GetStanceLegParam(int count, int legID)
{
    double pEB[6] {0};
    double pEE[18] {0};

    b_sb[count]=stepD/4-stepD/2*(s_b/PI);
    db_sb[count]=-stepD/2/PI;
    ddb_sb[count]=0;

    pEB[2]=initPeb[2]+b_sb[count];
    pEE[3*legID]=initPee[3*legID];
    pEE[3*legID+1]=initPee[3*legID+1];
    pEE[3*legID+2]=initPee[3*legID+2];
    rbt.SetPeb(pEB);
    rbt.pLegs[legID]->SetPee(pEE+3*legID);

    rbt.pLegs[legID]->GetJvi(Jvi,rbt.body());
    rbt.pLegs[legID]->GetdJacOverPee(dJvi_x,dJvi_y,dJvi_z,"B");
    rbt.pLegs[legID]->GetPee(*output_PeeB+18*count+3*legID,rbt.body());

    double db_st3[3] {0};
    db_st3[2]=-db_sb[count];
    double ddb_st3[3] {0};
    ddb_st3[2]=-ddb_sb[count];

    std::fill_n(dJvi,9,0);
    aris::dynamic::s_daxpy(9,db_st3[0],dJvi_x,1,dJvi,1);//for s_t
    aris::dynamic::s_daxpy(9,db_st3[1],dJvi_y,1,dJvi,1);
    aris::dynamic::s_daxpy(9,db_st3[2],dJvi_z,1,dJvi,1);

    std::fill_n(param_dds+3*legID,3,0);
    aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,db_st3,1,1,param_dds+3*legID,1);

    std::fill_n(param_dsds1+3*legID,3,0);
    aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,ddb_st3,1,1,param_dsds1+3*legID,1);
    std::fill_n(param_dsds2+3*legID,3,0);
    aris::dynamic::s_dgemm(3,1,3,1,dJvi,3,db_st3,1,1,param_dsds2+3*legID,1);
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
            printf("WARNING!!! param_dds equals zero!!! StanceLeg : %d \n",count);
        }
    }
}

void TimeOptimalGait::GetSwingLegParam(int count, int legID)
{
    double pEB[6] {0};
    double pEE[18] {0};

    f_sw[0]=0;
    f_sw[1]=stepH*sin(s_w);
    f_sw[2]=stepD/2*cos(s_w);
    df_sw[0]=0;
    df_sw[1]=stepH*cos(s_w);
    df_sw[2]=-stepD/2*sin(s_w);
    ddf_sw[0]=0;
    ddf_sw[1]=-stepH*sin(s_w);
    ddf_sw[2]=-stepD/2*cos(s_w);

    pEB[2]=initPeb[2]+pb_sw_tmp[count];
    pEE[3*legID]=initPee[3*legID]+f_sw[0];
    pEE[3*legID+1]=initPee[3*legID+1]+f_sw[1];
    pEE[3*legID+2]=initPee[3*legID+2]+f_sw[2];
    rbt.SetPeb(pEB);
    rbt.pLegs[legID]->SetPee(pEE+3*legID);

    rbt.pLegs[legID]->GetJvi(Jvi,rbt.body());
    rbt.pLegs[legID]->GetdJacOverPee(dJvi_x,dJvi_y,dJvi_z,"B");
    rbt.pLegs[legID]->GetPee(*output_PeeB+18*count+3*legID,rbt.body());

    double vb_sw3[3] {0};
    vb_sw3[2]=vb_sw_tmp[count];
    double ab_sw3[3] {0};
    ab_sw3[2]=ab_sw_tmp[count];

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
            printf("WARNING!!! param_dds equals zero!!! SwingLeg : %d \n",count);
        }
    }
}

void TimeOptimalGait::GetStanceDsBound(int count)
{
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
                dec[3*j+k]=param_a2[3*stanceLegID[j]+k]*ds*ds+param_a0L[3*stanceLegID[j]+k];
                acc[3*j+k]=param_a2[3*stanceLegID[j]+k]*ds*ds+param_a0H[3*stanceLegID[j]+k];
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
                    ds=ds_upBound_aLmt[count][stanceLegID[0]]+0.001*0.0001*k_st2;
                    for (int j=0;j<3;j++)
                    {
                        for (int k=0;k<3;k++)
                        {
                            dec[3*j+k]=param_a2[3*stanceLegID[j]+k]*ds*ds+param_a0L[3*stanceLegID[j]+k];
                            acc[3*j+k]=param_a2[3*stanceLegID[j]+k]*ds*ds+param_a0H[3*stanceLegID[j]+k];
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
                        dds_lowBound[count][stanceLegID[0]]=max_dec;
                        dds_upBound[count][stanceLegID[0]]=min_acc;
                        ds_upBound_aLmt[count][stanceLegID[0]]=ds;
                    }
                }
            }
            else
            {
                ds_upBound_aLmt[count][stanceLegID[0]]=ds;
            }
        }
    }

    ds_upBound_vLmt[count][stanceLegID[0]]=vLmt/(*std::max_element(abs_param_dds,abs_param_dds+18));
    ds_upBound[count][stanceLegID[0]]=std::min(ds_upBound_aLmt[count][stanceLegID[0]],ds_upBound_vLmt[count][stanceLegID[0]]);
    ds_lowBound[count][stanceLegID[0]]=ds_lowBound_aLmt[count][stanceLegID[0]];
}

void TimeOptimalGait::GetSwingDsBound(int count, int legID)
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
            dec[k]=param_a2[3*legID+k]*ds*ds+param_a1[3*legID+k]*ds+param_a0L[3*legID+k];
            acc[k]=param_a2[3*legID+k]*ds*ds+param_a1[3*legID+k]*ds+param_a0H[3*legID+k];
        }
        max_dec=*std::max_element(dec,dec+3);
        min_acc=*std::min_element(acc,acc+3);

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
                ds_lowBound_aLmt[count][legID]=ds;

                if(k_sw==1)
                {
                    //printf("swing ds_lowBound_aLmt=%.6f\n",ds_lowBound_aLmt[i][j+1]);
                    ds_lowBoundFlag_sw=true;
                }
                else
                {
                    for(int k_sw2=0;k_sw2<10000;k_sw2++)
                    {
                        ds=ds_lowBound_aLmt[count][legID]-0.001*0.0001*k_sw2;
                        for (int k=0;k<3;k++)
                        {
                            dec[k]=param_a2[3*legID+k]*ds*ds+param_a1[3*legID+k]*ds+param_a0L[3*legID+k];
                            acc[k]=param_a2[3*legID+k]*ds*ds+param_a1[3*legID+k]*ds+param_a0H[3*legID+k];
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
                            ds_lowBound_aLmt[count][legID]=ds;
                        }
                    }
                }
            }
            else if(ds_lowBoundFlag_sw==true && ds_upBoundFlag_sw==false && min_acc>=max_dec)
            {
                dds_lowBound[count][legID]=max_dec;
                dds_upBound[count][legID]=min_acc;
                ds_upBound_aLmt[count][legID]=ds;
            }
            else if(ds_lowBoundFlag_sw==true && ds_upBoundFlag_sw==false && min_acc<max_dec)
            {
                for(int k_sw2=0;k_sw2<10000;k_sw2++)
                {
                    ds=ds_upBound_aLmt[count][legID]+0.001*0.0001*k_sw2;
                    for (int k=0;k<3;k++)
                    {
                        dec[k]=param_a2[3*legID+k]*ds*ds+param_a1[3*legID+k]*ds+param_a0L[3*legID+k];
                        acc[k]=param_a2[3*legID+k]*ds*ds+param_a1[3*legID+k]*ds+param_a0H[3*legID+k];
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
                        dds_lowBound[count][legID]=max_dec;
                        dds_upBound[count][legID]=min_acc;
                        ds_upBound_aLmt[count][legID]=ds;
                    }
                }
            }
        }
    }

    double vLmt_value[3] {0};
    double Jvi_dot_vb[3] {0};
    double vb_sw3[3] {0};
    vb_sw3[2]=vb_sw_tmp[count];
    aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,vb_sw3,1,1,Jvi_dot_vb,1);
    for (int k=0;k<3;k++)
    {
        if(param_dds[3*legID+k]>0)
        {
            vLmt_value[k]=(Jvi_dot_vb[k]+vLmt)/param_dds[3*legID+k];
        }
        else if(param_dds[3*legID+k]<0)
        {
            vLmt_value[k]=(Jvi_dot_vb[k]-vLmt)/param_dds[3*legID+k];
        }
        else
        {
            printf("WARNING!!! param_dds equals zero!!! SwingLeg : %d \n",count);
        }

    }
    ds_upBound_vLmt[count][legID]=*std::max_element(vLmt_value,vLmt_value+3);
    ds_upBound[count][legID]=std::min(ds_upBound_aLmt[count][legID],ds_upBound_vLmt[count][legID]);
    ds_lowBound[count][legID]=ds_lowBound_aLmt[count][legID];
}

void TimeOptimalGait::GetStanceSwitchPoint()
{
    //initialize
    for(int i=0;i<1801;i++)
    {
        tangentPoint[i][stanceLegID[0]]=-1;
        switchPoint[i][stanceLegID[0]]=-1;
        vaCrossPoint[i][stanceLegID[0]]=-1;
        for(int j=0;j<3;j++)
        {
            for(int k=0;k<3;k++)
            {
                paramdds0Point[i][3*stanceLegID[j]+k]=-1;
            }
        }
    }

    //calculate the slope of ds_upBound
    slopedsBound[0][stanceLegID[0]]=(ds_upBound[1][stanceLegID[0]]-ds_upBound[0][stanceLegID[0]])/delta_s;
    for(int i=1;i<1800;i++)
    {
        slopedsBound[i][stanceLegID[0]]=(ds_upBound[i+1][stanceLegID[0]]-ds_upBound[i-1][stanceLegID[0]])/2/delta_s;
    }
    slopedsBound[1800][stanceLegID[0]]=(ds_upBound[1800][stanceLegID[0]]-ds_upBound[1799][stanceLegID[0]])/delta_s;
    for(int i=0;i<1801;i++)
    {
        slopeDelta[i][stanceLegID[0]]=slopedsBound[i][stanceLegID[0]]-(dds_upBound[i][stanceLegID[0]]+dds_lowBound[i][stanceLegID[0]])/2;
        vaDelta[i][stanceLegID[0]]=ds_upBound_vLmt[i][stanceLegID[0]]-ds_upBound_aLmt[i][stanceLegID[0]];
    }

    for(int i=0;i<1800;i++)
    {
        if(ds_upBound_vLmt[i][stanceLegID[0]]>ds_upBound_aLmt[i][stanceLegID[0]])
        {
            for(int j=0;j<3;j++)
            {
                for(int k=0;k<3;k++)
                {
                    if(output_dds[i+1][3*stanceLegID[j]+k]*output_dds[i][3*stanceLegID[j]+k]<0 || output_dds[i][3*stanceLegID[j]+k]==0)
                    {
                        paramdds0Point[paramdds0Count[3*stanceLegID[j]+k]][3*stanceLegID[j]+k]=i
                                +fabs(output_dds[i][3*stanceLegID[j]+k])/(fabs(output_dds[i][3*stanceLegID[j]+k])+fabs(output_dds[i+1][3*stanceLegID[j]+k]));
                        paramdds0Count[6*j+k+3]++;
                    }
                }
            }
        }

        if(slopeDelta[i+1][stanceLegID[0]]*slopeDelta[i][stanceLegID[0]]<0 || slopeDelta[i][stanceLegID[0]]==0)
        {
            tangentPoint[tangentCount[stanceLegID[0]]][stanceLegID[0]]=i
                    +fabs(slopeDelta[i][stanceLegID[0]])/(fabs(slopeDelta[i][stanceLegID[0]])+fabs(slopeDelta[i+1][stanceLegID[0]]));
            tangentCount[0]++;
        }

        if(vaDelta[i+1][stanceLegID[0]]*vaDelta[i][stanceLegID[0]]<0 || vaDelta[i][stanceLegID[0]]==0)
        {
            vaCrossPoint[vaCrossCount[stanceLegID[0]]][stanceLegID[0]]=i
                    +fabs(vaDelta[i][stanceLegID[0]])/(fabs(vaDelta[i][stanceLegID[0]])+fabs(vaDelta[i+1][stanceLegID[0]]));
            vaCrossCount[0]++;
        }
    }

    //merge tangentPoint & paramdds0Point into switchPoint
    switchPoint[0][stanceLegID[0]]=0;
    switchPoint[1][stanceLegID[0]]=1800;
    switchCount[stanceLegID[0]]=2;
    for(int i=0;i<tangentCount[stanceLegID[0]];i++)
    {
        switchPoint[i+switchCount[stanceLegID[0]]][stanceLegID[0]]=tangentPoint[i][stanceLegID[0]];
    }
    switchCount[stanceLegID[0]]+=tangentCount[stanceLegID[0]];
    for(int i=0;i<vaCrossCount[stanceLegID[0]];i++)
    {
        switchPoint[i+switchCount[stanceLegID[0]]][stanceLegID[0]]=vaCrossPoint[i][stanceLegID[0]];
    }
    switchCount[stanceLegID[0]]+=vaCrossCount[stanceLegID[0]];
    for(int j=0;j<3;j++)
    {
        for(int k=0;k<3;k++)
        {
            for(int i=0;i<paramdds0Count[3*stanceLegID[j]+k];i++)
            {
                switchPoint[i+switchCount[stanceLegID[0]]][stanceLegID[0]]=paramdds0Point[i][3*stanceLegID[j]+k];
            }
            switchCount[stanceLegID[0]]+=paramdds0Count[3*stanceLegID[j]+k];
        }
    }

    //filtering the same point & sorting by the value
    for(int i=0;i<switchCount[stanceLegID[0]];i++)
    {
        for(int j=i+1;j<switchCount[stanceLegID[0]];j++)
        {
            if((int)switchPoint[j][stanceLegID[0]]<(int)switchPoint[i][stanceLegID[0]])
            {
                auto tmp=switchPoint[i][stanceLegID[0]];
                switchPoint[i][stanceLegID[0]]=switchPoint[j][stanceLegID[0]];
                switchPoint[j][stanceLegID[0]]=tmp;
            }
            else if((int)switchPoint[j][stanceLegID[0]]==(int)switchPoint[i][stanceLegID[0]])
            {
                switchPoint[j][stanceLegID[0]]=switchPoint[switchCount[stanceLegID[0]]-1][stanceLegID[0]];
                switchPoint[switchCount[stanceLegID[0]]-1][stanceLegID[0]]=-1;
                j--;
                switchCount[stanceLegID[0]]--;
            }
        }
    }

    printf("StanceLeg Switch Point:");
    for(int i=0;i<switchCount[stanceLegID[0]]+1;i++)
    {
        printf("%.1f,",switchPoint[i][stanceLegID[0]]);
    }
    printf("\n");
}

void TimeOptimalGait::GetSwingSwitchPoint(int legID)
{
    //initialize
    for(int k=0;k<3;k++)
    {
        paramdds0Count[3*legID+k]=0;
    }
    tangentCount[legID]=0;
    switchCount[legID]=0;
    vaCrossCount[legID]=0;
    for(int i=0;i<1801;i++)
    {
        vaCrossPoint[i][legID]=-1;
        tangentPoint[i][legID]=-1;
        switchPoint[i][legID]=-1;
        for(int k=0;k<3;k++)
        {
            paramdds0Point[i][3*legID+k]=-1;
        }
    }

    slopedsBound[0][legID]=(ds_upBound[1][legID]-ds_upBound[0][legID])/delta_s;
    for(int i=1;i<1800;i++)
    {
        slopedsBound[i][legID]=(ds_upBound[i+1][legID]-ds_upBound[i-1][legID])/2/delta_s;
    }
    slopedsBound[1800][legID]=(ds_upBound[1800][legID]-ds_upBound[1799][legID])/delta_s;

    for(int i=0;i<1801;i++)
    {
        slopeDelta[i][legID]=slopedsBound[i][legID]-(dds_upBound[i][legID]+dds_lowBound[i][legID])/2;
        vaDelta[i][legID]=ds_upBound_vLmt[i][legID]-ds_upBound_aLmt[i][legID];
    }

    for(int i=0;i<1800;i++)
    {
        if(ds_upBound_vLmt[i][legID]>ds_upBound_aLmt[i][legID])
        {
            for(int k=0;k<3;k++)
            {
                if(output_dds[i+1][3*legID+k]*output_dds[i][3*legID+k]<0 || output_dds[i][3*legID+k]==0)
                {
                    paramdds0Point[paramdds0Count[3*legID+k]][3*legID+k]=i+fabs(output_dds[i][3*legID+k])/(fabs(output_dds[i][3*legID+k])+fabs(output_dds[i+1][3*legID+k]));;
                    paramdds0Count[3*legID+k]++;
                }
            }
        }

        if(slopeDelta[i+1][legID]*slopeDelta[i][legID]<0 || slopeDelta[i][legID]==0)
        {
            tangentPoint[tangentCount[legID]][legID]=i+fabs(slopeDelta[i][legID])/(fabs(slopeDelta[i][legID])+fabs(slopeDelta[i+1][legID]));
            tangentCount[legID]++;
        }

        if(vaDelta[i+1][legID]*vaDelta[i][legID]<0 || vaDelta[i][legID]==0)
        {
            vaCrossPoint[vaCrossCount[legID]][legID]=i+fabs(vaDelta[i][legID])/(fabs(vaDelta[i][legID])+fabs(vaDelta[i+1][legID]));
            vaCrossCount[0]++;
        }
    }

    //merge tangentPoint & paramdds0Point into switchPoint
    switchPoint[0][legID]=0;
    switchPoint[1][legID]=1800;
    switchCount[legID]=2;
    for(int i=0;i<tangentCount[legID];i++)
    {
        switchPoint[i+switchCount[legID]][legID]=tangentPoint[i][legID];
    }
    switchCount[legID]+=tangentCount[legID];
    for(int i=0;i<vaCrossCount[legID];i++)
    {
        switchPoint[i+switchCount[legID]][legID]=vaCrossPoint[i][legID];
    }
    switchCount[legID]+=vaCrossCount[legID];
    for(int k=0;k<3;k++)
    {
        for(int i=0;i<paramdds0Count[3*legID+k];i++)
        {
            switchPoint[i+switchCount[legID]][legID]=paramdds0Point[i][3*legID+k];
        }
        switchCount[legID]+=paramdds0Count[3*legID+k];
    }

    //filtering the same point & sorting by the value
    for(int i=0;i<switchCount[legID];i++)
    {
        for(int k=i+1;k<switchCount[legID];k++)
        {
            if((int)switchPoint[k][legID]<(int)switchPoint[i][legID])
            {
                auto tmp=switchPoint[i][legID];
                switchPoint[i][legID]=switchPoint[k][legID];
                switchPoint[k][legID]=tmp;
            }
            else if((int)switchPoint[k][legID]==(int)switchPoint[i][legID])
            {
                switchPoint[k][legID]=switchPoint[switchCount[legID]-1][legID];
                switchPoint[switchCount[legID]-1][legID]=-1;
                k--;
                switchCount[legID]--;
            }
        }
    }

    printf("SwingLeg Switch Point:");
    for(int i=0;i<switchCount[legID]+1;i++)
    {
        printf("%.1f,",switchPoint[i][legID]);
    }
    printf("\n");
}

void TimeOptimalGait::GetStanceOptimalDsBySwitchPoint()
{
    for(int i=0;i<1801;i++)
    {
        real_ds[i][stanceLegID[0]]=ds_upBound[i][stanceLegID[0]];
        real_dds[i][stanceLegID[0]]=dds_upBound[i][stanceLegID[0]];
    }
    for(int m=0;m<switchCount[stanceLegID[0]];m++)
    {
        int k_st {(int)switchPoint[m][stanceLegID[0]]};
        if(ds_upBound[k_st][stanceLegID[0]]>real_ds[k_st][stanceLegID[0]] && (int)switchPoint[m][stanceLegID[0]]!=1800)//start of backward
        {
            printf("StanceLeg backward start at a passed point, quit switchPoint %.1f\n",switchPoint[m][stanceLegID[0]]);
            continue;
        }
        stopFlag=false;
        quitSwitchPoint=false;
        ds_backward[k_st][stanceLegID[0]]=ds_upBound[k_st][stanceLegID[0]];

        while(stopFlag==false && (int)switchPoint[m][stanceLegID[0]]!=0)
        {
            double dec[9] {0};
            for (int j=0;j<3;j++)
            {
                for (int k=0;k<3;k++)
                {
                    dec[3*j+k]=output_a2[k_st][3*stanceLegID[j]+k]*ds_backward[k_st][stanceLegID[0]]*ds_backward[k_st][stanceLegID[0]]+output_a0L[k_st][3*stanceLegID[j]+k];
                }
            }
            dds_backward[k_st][stanceLegID[0]]=*std::max_element(dec,dec+9);
            ds_backward[k_st-1][stanceLegID[0]]=sqrt(ds_backward[k_st][stanceLegID[0]]*ds_backward[k_st][stanceLegID[0]]-2*dds_backward[k_st][stanceLegID[0]]*delta_s);

            if(ds_backward[k_st-1][stanceLegID[0]]>ds_upBound[k_st-1][stanceLegID[0]])
            {
                stopFlag=true;
                quitSwitchPoint=true;
                printf("StanceLeg backward touching upBound at %d, quit switchPoint %.4f\n",k_st-1,switchPoint[m][stanceLegID[0]]);
            }
            else if(k_st==1)
            {
                for (int j=0;j<3;j++)
                {
                    for (int k=0;k<3;k++)
                    {
                        dec[3*j+k]=output_a2[k_st-1][3*stanceLegID[j]+k]*ds_backward[k_st-1][stanceLegID[0]]*ds_backward[k_st-1][stanceLegID[0]]+output_a0L[k_st-1][3*stanceLegID[j]+k];
                    }
                }
                dds_backward[k_st-1][stanceLegID[0]]=*std::max_element(dec,dec+9);
                for(int i=k_st-1;i<switchPoint[m][stanceLegID[0]]+1;i++)
                {
                    real_ds[i][stanceLegID[0]]=ds_backward[i][stanceLegID[0]];
                    real_dds[i][stanceLegID[0]]=dds_backward[i][stanceLegID[0]];
                }
                stopFlag=true;
                printf("StanceLeg backward touching 0, from switchPoint %.4f\n",switchPoint[m][stanceLegID[0]]);
            }
            else if(ds_backward[k_st-1][stanceLegID[0]]>=real_ds[k_st-1][stanceLegID[0]])
            {
                for (int j=0;j<3;j++)
                {
                    for (int k=0;k<3;k++)
                    {
                        dec[3*j+k]=output_a2[k_st-1][3*stanceLegID[j]+k]*ds_backward[k_st-1][stanceLegID[0]]*ds_backward[k_st-1][stanceLegID[0]]+output_a0L[k_st-1][3*stanceLegID[j]+k];
                    }
                }
                dds_backward[k_st-1][stanceLegID[0]]=*std::max_element(dec,dec+9);
                for(int i=k_st-1;i<switchPoint[m][stanceLegID[0]]+1;i++)
                {
                    real_ds[i][stanceLegID[0]]=ds_backward[i][stanceLegID[0]];
                    real_dds[i][stanceLegID[0]]=dds_backward[i][stanceLegID[0]];
                }
                stopFlag=true;
                printf("StanceLeg backward touching last curve at %d, from switchPoint %.4f\n",k_st-1,switchPoint[m][stanceLegID[0]]);
            }
            else
            {
                k_st--;
            }
        }
        if(quitSwitchPoint==true || (int)switchPoint[m][stanceLegID[0]]==1800)//end of backward
        {
            continue;
        }

        bool isEqual0Point {false};
        for(int j=0;j<3;j++)
        {
            for(int k=0;k<3;k++)
            {
                if(output_dds[(int)switchPoint[m][stanceLegID[0]]][3*stanceLegID[j]+k]==0
                        || slopeDelta[(int)switchPoint[m][stanceLegID[0]]][stanceLegID[0]]==0
                        || vaDelta[(int)switchPoint[m][stanceLegID[0]]][stanceLegID[0]]==0)
                {
                    isEqual0Point=true;
                }
            }
        }
        if(isEqual0Point==false  && switchPoint[m][stanceLegID[0]]!=0)
        {
            k_st=switchPoint[m][stanceLegID[0]]+1;
        }
        else
        {
            k_st=switchPoint[m][stanceLegID[0]];
        }
        stopFlag=false;
        ds_forward[k_st][stanceLegID[0]]=ds_upBound[k_st][stanceLegID[0]];
        while(stopFlag==false)
        {
            double acc[9] {0};
            for (int j=0;j<3;j++)
            {
                for (int k=0;k<3;k++)
                {
                    acc[3*j+k]=output_a2[k_st][3*stanceLegID[j]+k]*ds_forward[k_st][stanceLegID[0]]*ds_forward[k_st][stanceLegID[0]]+output_a0H[k_st][3*stanceLegID[j]+k];
                }
            }
            dds_forward[k_st][stanceLegID[0]]=*std::min_element(acc,acc+9);
            ds_forward[k_st+1][stanceLegID[0]]=sqrt(ds_forward[k_st][stanceLegID[0]]*ds_forward[k_st][stanceLegID[0]]+2*dds_forward[k_st][stanceLegID[0]]*delta_s);

            if(ds_forward[k_st+1][stanceLegID[0]]>ds_upBound[k_st+1][stanceLegID[0]] || k_st==1799)
            {
                for(int i=switchPoint[m][stanceLegID[0]]+1;i<k_st+1;i++)
                {
                    real_ds[i][stanceLegID[0]]=ds_forward[i][stanceLegID[0]];
                    real_dds[i][stanceLegID[0]]=dds_forward[i][stanceLegID[0]];
                }
                if(switchPoint[m][stanceLegID[0]]==0)
                {
                    real_ds[0][stanceLegID[0]]=ds_forward[0][stanceLegID[0]];
                    real_dds[0][stanceLegID[0]]=dds_forward[0][stanceLegID[0]];
                }
                if(k_st==1799)
                {
                    real_ds[1800][stanceLegID[0]]=std::min(ds_upBound[1800][stanceLegID[0]],ds_forward[1800][stanceLegID[0]]);
                    real_dds[1800][stanceLegID[0]]=std::min(dds_upBound[1800][stanceLegID[0]],dds_forward[1800][stanceLegID[0]]);
                }
                stopFlag=true;
                printf("StanceLeg forward touching upBound at %d, from switchPoint %.4f\n",k_st,switchPoint[m][stanceLegID[0]]);
            }
            else
            {
                k_st++;
            }
        }
    }
}

void TimeOptimalGait::GetStanceOptimalDsByIteration()
{
    //backward integration
        stopFlag=false;
        ki_back=1800;
        ds_backward[ki_back][stanceLegID[0]]=ds_upBound[ki_back][stanceLegID[0]];
        while (stopFlag==false && ki_back>=0)
        {
            double dec[9] {0};
            for (int j=0;j<3;j++)
            {
                for (int k=0;k<3;k++)
                {
                    dec[3*j+k]=output_a2[ki_back][3*stanceLegID[j]+k]*ds_backward[ki_back][stanceLegID[0]]*ds_backward[ki_back][stanceLegID[0]]
                              +output_a0L[ki_back][3*stanceLegID[j]+k];
                }
            }
            dds_backward[ki_back][stanceLegID[0]]=*std::max_element(dec,dec+9);
            ds_backward[ki_back-1][stanceLegID[0]]=sqrt(ds_backward[ki_back][stanceLegID[0]]*ds_backward[ki_back][stanceLegID[0]]-2*dds_backward[ki_back][stanceLegID[0]]*delta_s);

            if (ds_backward[ki_back-1][stanceLegID[0]]>ds_upBound[ki_back-1][stanceLegID[0]])
            {
                stopFlag=true;
                stop_back=ki_back;
                printf("StanceLeg Backward Integration ends at k=%d, ds_backward:%.4f; ds_backward[ki_back-1][0]=%.4f > ds_upBound[ki_back-1][0]=%.4f\n",
                       ki_back,ds_backward[ki_back][stanceLegID[0]],ds_backward[ki_back-1][stanceLegID[0]],ds_upBound[ki_back-1][stanceLegID[0]]);
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
        ds_forward[ki_for][stanceLegID[0]]=ds_upBound[ki_for][stanceLegID[0]];
        std::fill_n(min_dist,1801,1);
        while (stopFlag==false)
        {
            if (accFlag==true)
            {
                double acc[9] {0};
                for (int j=0;j<3;j++)
                {
                    for (int k=0;k<3;k++)
                    {
                        acc[3*j+k]=output_a2[ki_for][3*stanceLegID[j]+k]*ds_forward[ki_for][stanceLegID[0]]*ds_forward[ki_for][stanceLegID[0]]
                                  +output_a0H[ki_for][3*stanceLegID[j]+k];
                    }
                }
                dds_forward[ki_for][stanceLegID[0]]=*std::min_element(acc,acc+9);
                ds_forward[ki_for+1][stanceLegID[0]]=sqrt(ds_forward[ki_for][stanceLegID[0]]*ds_forward[ki_for][stanceLegID[0]]+2*dds_forward[ki_for][stanceLegID[0]]*delta_s);
                //printf("acc,sqrt:%.4f\n",ds_forward[ki_for][0]*ds_forward[ki_for][0]+2*dds_forward[ki_for][0]*delta_s);

                if (ds_forward[ki_for+1][stanceLegID[0]]>ds_upBound[ki_for+1][stanceLegID[0]])
                {
                    accFlag=false;
                    dec_start=ki_for;
                    printf("StanceLeg acc reach bound at k=%d, ds_forward=%.4f; ",ki_for,ds_forward[ki_for][stanceLegID[0]]);
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
                        dec[3*j+k]=output_a2[ki_for][3*stanceLegID[j]+k]*ds_forward[ki_for][stanceLegID[0]]*ds_forward[ki_for][stanceLegID[0]]
                                  +output_a0L[ki_for][3*stanceLegID[j]+k];
                    }
                }
                dds_forward[ki_for][stanceLegID[0]]=*std::max_element(dec,dec+9);
                ds_forward[ki_for+1][stanceLegID[0]]=sqrt(ds_forward[ki_for][stanceLegID[0]]*ds_forward[ki_for][stanceLegID[0]]+2*dds_forward[ki_for][stanceLegID[0]]*delta_s);

                if (ds_forward[ki_for+1][stanceLegID[0]]>ds_upBound[ki_for+1][stanceLegID[0]])
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
                        ds_forward[0][stanceLegID[0]]-=ds_upBound[0][stanceLegID[0]]/1000;
                    }
                }
                else
                {
                    if ((ds_forward[ki_for][stanceLegID[0]]*ds_forward[ki_for][stanceLegID[0]]+2*dds_forward[ki_for][stanceLegID[0]]*delta_s)
                            <=(ds_lowBound[ki_for+1][stanceLegID[0]]*ds_lowBound[ki_for+1][stanceLegID[0]]) || ki_for==(1799))
                    {
                        accFlag=true;
                        for(int k=dec_start;k<(ki_for+2);k++)
                        {
                            min_dist[k]=ds_upBound[k][stanceLegID[0]]-ds_forward[k][stanceLegID[0]];
                        }
                        dec_end=std::min_element(min_dist+dec_start+1,min_dist+ki_for+2)-min_dist;
                        //dec_start must be ignored, if dec_start is the min_dist, the calculation will cycle between dec_start & dec_start+1
                        ki_for=dec_end-1;
                        printf("dec finished, start at k=%d, end at k=%d, ds_forward=%.4f\n",dec_start,dec_end,ds_forward[ki_for+1][stanceLegID[0]]);
                    }
                    ki_for++;
                }
            }

            //loop end condition
            if(ki_for==1800)
            {
                stopFlag=true;
                for(int i=0;i<1801;i++)
                {
                    real_ds[i][stanceLegID[0]]=ds_forward[i][stanceLegID[0]];
                    real_dds[i][stanceLegID[0]]=dds_forward[i][stanceLegID[0]];
                    //printf("i=%d, ds_forward:%.4f, dds_forward:%.4f\n",i,ds_forward[i][0],dds_forward[i][0]);
                }
                printf("StanceLeg forward reach the end, and never encounter with the backward, cycleCount=%u < %u\n",cycleCount,0x0FFFFFFF);
            }
            if(ki_for>=stop_back && ds_forward[ki_for][stanceLegID[0]]>=ds_backward[ki_for][stanceLegID[0]])
            {
                stopFlag=true;
                for(int i=0;i<ki_for;i++)
                {
                    real_ds[i][stanceLegID[0]]=ds_forward[i][stanceLegID[0]];
                    real_dds[i][stanceLegID[0]]=dds_forward[i][stanceLegID[0]];
                }
                for(int i=ki_for;i<1801;i++)
                {
                    real_ds[i][stanceLegID[0]]=ds_backward[i][stanceLegID[0]];
                    real_dds[i][stanceLegID[0]]=dds_backward[i][stanceLegID[0]];
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

void TimeOptimalGait::GetSwingOptimalDsBySwitchPoint(int legID)
{
    for(int i=0;i<1801;i++)
    {
        real_ds[i][legID]=ds_upBound[i][legID];
        real_dds[i][legID]=dds_upBound[i][legID];
    }
    for(int m=0;m<switchCount[legID];m++)
    {
        int k_sw {(int)switchPoint[m][legID]};
        if(ds_upBound[k_sw][legID]>real_ds[k_sw][legID] && (int)switchPoint[m][legID]!=1800)//start of backward
        {
            printf("SwingLeg backward start at a passed point, quit switchPoint %.1f\n",switchPoint[m][legID]);
            continue;
        }
        stopFlag=false;
        if((int)switchPoint[m][legID]==1800)
        {
            ds_backward[k_sw][legID]=0;
        }
        else
        {
            ds_backward[k_sw][legID]=ds_upBound[k_sw][legID];
        }

        while(stopFlag==false && (int)switchPoint[m][legID]!=0)
        {
            double dec[3] {0};
            for (int k=0;k<3;k++)
            {
                dec[k]=output_a2[k_sw][3*legID+k]*ds_backward[k_sw][legID]*ds_backward[k_sw][legID]
                      +output_a1[k_sw][3*legID+k]*ds_backward[k_sw][legID]
                      +output_a0L[k_sw][3*legID+k];
            }
            dds_backward[k_sw][legID]=*std::max_element(dec,dec+3);
            ds_backward[k_sw-1][legID]=sqrt(ds_backward[k_sw][legID]*ds_backward[k_sw][legID]-2*dds_backward[k_sw][legID]*delta_s);

            if(ds_backward[k_sw-1][legID]>ds_upBound[k_sw-1][legID])
            {
                stopFlag=true;
                printf("SwingLeg backward touching upBound at %d, quit switchPoint %.4f\n",k_sw-1,switchPoint[m][legID]);
            }
            else if(k_sw==1)
            {
                for (int k=0;k<3;k++)
                {
                    dec[k]=output_a2[k_sw][3*legID+k]*ds_backward[k_sw][legID]*ds_backward[k_sw][legID]
                          +output_a1[k_sw][3*legID+k]*ds_backward[k_sw][legID]
                          +output_a0L[k_sw][3*legID+k];
                }
                dds_backward[k_sw-1][legID]=*std::max_element(dec,dec+3);
                for(int i=k_sw-1;i<switchPoint[m][legID]+1;i++)
                {
                    real_ds[i][legID]=ds_backward[i][legID];
                    real_dds[i][legID]=dds_backward[i][legID];
                }
                stopFlag=true;
                printf("SwingLeg backward touching 0, from switchPoint %.4f\n",switchPoint[m][legID]);
            }
            else if(ds_backward[k_sw-1][legID]>=real_ds[k_sw-1][legID])
            {
                for (int k=0;k<3;k++)
                {
                    dec[k]=output_a2[k_sw-1][3*legID+k]*ds_backward[k_sw-1][legID]*ds_backward[k_sw-1][legID]
                          +output_a1[k_sw-1][3*legID+k]*ds_backward[k_sw-1][legID]
                          +output_a0L[k_sw-1][3*legID+k];
                }
                dds_backward[k_sw-1][legID]=*std::max_element(dec,dec+3);
                for(int i=k_sw-1;i<switchPoint[m][legID]+1;i++)
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
        if(ds_backward[k_sw-1][legID]>ds_upBound[k_sw-1][legID] || (int)switchPoint[m][legID]==1800)
        {
            continue;
        }

        bool isEqual0Point {false};
        for(int k=0;k<3;k++)
        {
                if(output_dds[(int)switchPoint[m][legID]][3*legID+k]==0
                        || slopeDelta[(int)switchPoint[m][legID]][legID]==0
                        || vaDelta[(int)switchPoint[m][legID]][legID]==0)
                {
                    printf("switchPoint integral\n\n");
                    isEqual0Point=true;
                }
        }
        if(isEqual0Point==false && switchPoint[m][legID]!=0)
        {
            k_sw=switchPoint[m][legID]+1;
        }
        else
        {
            k_sw=switchPoint[m][legID];
        }
        stopFlag=false;
        if((int)switchPoint[m][legID]==0)
        {
            ds_forward[k_sw][legID]=0;
        }
        else
        {
            ds_forward[k_sw][legID]=ds_upBound[k_sw][legID];
        }
        while(stopFlag==false)
        {
            double acc[3] {0};
            for (int k=0;k<3;k++)
            {
                acc[k]=output_a2[k_sw][3*legID+k]*ds_forward[k_sw][legID]*ds_forward[k_sw][legID]
                      +output_a1[k_sw][3*legID+k]*ds_forward[k_sw][legID]
                      +output_a0H[k_sw][3*legID+k];
            }
            dds_forward[k_sw][legID]=*std::min_element(acc,acc+3);
            ds_forward[k_sw+1][legID]=sqrt(ds_forward[k_sw][legID]*ds_forward[k_sw][legID]+2*dds_forward[k_sw][legID]*delta_s);

            if(ds_forward[k_sw+1][legID]>ds_upBound[k_sw+1][legID] || k_sw==1799)
            {
                for(int i=switchPoint[m][legID]+1;i<k_sw+1;i++)
                {
                    real_ds[i][legID]=ds_forward[i][legID];
                    real_dds[i][legID]=dds_forward[i][legID];
                }
                if(switchPoint[m][legID]==0)
                {
                    real_ds[0][legID]=ds_forward[0][legID];
                    real_dds[0][legID]=dds_forward[0][legID];
                }
                if(k_sw==1799)
                {
                    real_ds[1800][legID]=std::min(ds_upBound[1800][legID],ds_forward[1800][legID]);
                    real_dds[1800][legID]=std::min(dds_upBound[1800][legID],dds_forward[1800][legID]);
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
}

void TimeOptimalGait::GetSwingOptimalDsByIteration(int legID)
{
    //backward integration
    stopFlag=false;
    ki_back=1800;
    ds_backward[ki_back][legID]=0;
    while (stopFlag==false && ki_back>=0)
    {
        double dec[3] {0};
        for (int k=0;k<3;k++)
        {
            dec[k]=output_a2[ki_back][3*legID+k]*ds_backward[ki_back][legID]*ds_backward[ki_back][legID]
                  +output_a1[ki_back][3*legID+k]*ds_backward[ki_back][legID]
                  +output_a0L[ki_back][3*legID+k];
        }
        dds_backward[ki_back][legID]=*std::max_element(dec,dec+3);
        ds_backward[ki_back-1][legID]=sqrt(ds_backward[ki_back][legID]*ds_backward[ki_back][legID]-2*dds_backward[ki_back][legID]*delta_s);

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
    accFlag=true;
    ki_for=0;
    cycleCount=0;
    ds_forward[ki_for][legID]=0;
    std::fill_n(min_dist,1801,1);
    while (stopFlag==false)
    {
        if (accFlag==true)
        {
            double acc[3] {0};
            for (int k=0;k<3;k++)
            {
                acc[k]=output_a2[ki_for][3*legID+k]*ds_forward[ki_for][legID]*ds_forward[ki_for][legID]
                      +output_a1[ki_for][3*legID+k]*ds_forward[ki_for][legID]
                      +output_a0H[ki_for][3*legID+k];
            }
            dds_forward[ki_for][legID]=*std::min_element(acc,acc+3);
            ds_forward[ki_for+1][legID]=sqrt(ds_forward[ki_for][legID]*ds_forward[ki_for][legID]+2*dds_forward[ki_for][legID]*delta_s);

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
            double dec[3] {0};
            for (int k=0;k<3;k++)
            {
                dec[k]=output_a2[ki_for][3*legID+k]*ds_forward[ki_for][legID]*ds_forward[ki_for][legID]
                      +output_a1[ki_for][3*legID+k]*ds_forward[ki_for][legID]
                      +output_a0L[ki_for][3*legID+k];
            }
            dds_forward[ki_for][legID]=*std::max_element(dec,dec+3);
            ds_forward[ki_for+1][legID]=sqrt(ds_forward[ki_for][legID]*ds_forward[ki_for][legID]+2*dds_forward[ki_for][legID]*delta_s);

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
                if ((ds_forward[ki_for][legID]*ds_forward[ki_for][legID]+2*dds_forward[ki_for][legID]*delta_s)
                        <=(ds_lowBound[ki_for+1][legID]*ds_lowBound[ki_for+1][legID]) || ki_for==(1799))
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

        if(ki_for==1800)
        {
            stopFlag=true;
            for(int i=0;i<1801;i++)
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
            for(int i=ki_for;i<1801;i++)
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


void TimeOptimalGait::GetOptimalDs()
{
    rbt.loadXml("/home/hex/Desktop/mygit/RobotVIII_demo/resource/Robot_VIII.xml");

    timeval tpstart,tpend;
    float tused;
    gettimeofday(&tpstart,NULL);

    /*********************************** StanceLeg ***************************************/
    //generate the traj & calculate the bound of ds
    for (int i=0;i<1801;i++)
    {
        s_b=i*PI/1800;//degree to rad

        for(int j=0;j<3;j++)
        {
            GetStanceLegParam(i,2*j+1);
        }
        GetStanceDsBound(i);

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

    GetStanceSwitchPoint();
    GetStanceOptimalDsBySwitchPoint();
    //GetStanceOptimalDsByIteration();

    //pb_sw, vb_sw, ab_sw initialized here, and need to be updated during iteration
    totalTime[stanceLegID[0]]=0;
    for (int i=0;i<1801;i++)
    {
        if(i!=0)
        {
            totalTime[stanceLegID[0]]+=2*delta_s/(real_ds[i-1][stanceLegID[0]]+real_ds[i][stanceLegID[0]]);
            timeArray[i][stanceLegID[0]]=totalTime[stanceLegID[0]];
        }
        pb_sw_tmp[i]=pb_sw[i]=b_sb[i];
        vb_sw_tmp[i]=vb_sw[i]=db_sb[i]*real_ds[i][stanceLegID[0]];
        ab_sw_tmp[i]=ab_sw[i]=ddb_sb[i]*real_ds[i][stanceLegID[0]]*real_ds[i][stanceLegID[0]]+db_sb[i]*real_dds[i][stanceLegID[0]];
    }


    /******************************* SwingLeg & Iteration ************************************/
    int iterCount {0};
    while(maxTotalCount!=maxTotalCount_last)
    {
        printf("\n");
        for(int i=0;i<1801;i++)
        {
            for(int j=0;j<3;j++)
            {
                ds_backward[i][2*j]=0;
                ds_forward[i][2*j]=0;
            }
        }
        iterCount++;
        maxTotalCount_last=maxTotalCount;
        //generate the traj & calculate the bound of ds
        for (int i=0;i<1801;i++)
        {
            s_w=i*PI/1800;
            for (int j=0;j<3;j++)
            {
                GetSwingLegParam(i,2*j);
            }

            for (int j=0;j<3;j++)
            {
                GetSwingDsBound(i,2*j);

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

        for(int j=0;j<3;j++)
        {
            GetSwingSwitchPoint(2*j);
        }

        for(int j=0;j<3;j++)
        {
            GetSwingOptimalDsBySwitchPoint(2*j);
            //GetSwingOptimalDsByIteration(2*j);

            totalTime[2*j]=0;
            for (int i=1;i<1801;i++)
            {
                totalTime[2*j]+=2*delta_s/(real_ds[i-1][2*j]+real_ds[i][2*j]);
                timeArray[i][2*j]=totalTime[2*j];
            }
        }
        maxTime=*std::max_element(totalTime,totalTime+6);
        maxTotalCount=(int)(maxTime*1000)+1;
        maxTime=0.001*maxTotalCount;
        maxTimeID=std::max_element(totalTime,totalTime+6)-totalTime;
        printf("totalTime: %.4f, %.4f, %.4f, %.4f, %.4f, %.4f; maxTime:%.4f\n",totalTime[0],totalTime[1],totalTime[2],totalTime[3],totalTime[4],totalTime[5],maxTime);

        //update pb_sw, vb_sw, ab_sw here
        for(int i=0;i<1801;i++)
        {
            memcpy(*timeArray_tmp+6*i,*timeArray+6*i,6*sizeof(double));
            for(int j=0;j<6;j++)
            {
                if(totalTime[j]!=0)
                {
                    timeArray_tmp[i][j]*=maxTime/totalTime[j];
                }
            }
        }
        int j_start {0};
        for(int i=0;i<1800;i++)
        {
            //if(i%100==99)
                //printf("j_start=%d, ",j_start);
            for(int j=j_start;j<1800;j++)
            {
                if(timeArray_tmp[i][maxTimeID]>=timeArray_tmp[j][stanceLegID[0]] && timeArray_tmp[i][maxTimeID]<timeArray_tmp[j+1][stanceLegID[0]])
                {
                    j_start=j;
                    pb_sw_tmp[i]=(pb_sw[j+1]-pb_sw[j])/(timeArray_tmp[j+1][stanceLegID[0]]-timeArray_tmp[j][stanceLegID[0]])*(timeArray_tmp[i][maxTimeID]-timeArray_tmp[j][stanceLegID[0]])+pb_sw[j];
                    vb_sw_tmp[i]=(vb_sw[j+1]-vb_sw[j])/(timeArray_tmp[j+1][stanceLegID[0]]-timeArray_tmp[j][stanceLegID[0]])*(timeArray_tmp[i][maxTimeID]-timeArray_tmp[j][stanceLegID[0]])+vb_sw[j];
                    vb_sw_tmp[i]*=totalTime[stanceLegID[0]]/maxTime;
                    ab_sw_tmp[i]=(ab_sw[j+1]-ab_sw[j])/(timeArray_tmp[j+1][stanceLegID[0]]-timeArray_tmp[j][stanceLegID[0]])*(timeArray_tmp[i][maxTimeID]-timeArray_tmp[j][stanceLegID[0]])+ab_sw[j];
                    ab_sw_tmp[i]*=totalTime[stanceLegID[0]]/maxTime*totalTime[stanceLegID[0]]/maxTime;
                    break;
                }
            }
            //if(i%100==99)
                //printf("j_end=%d\n",j_start);
        }
        pb_sw_tmp[1800]=pb_sw[1800];
        vb_sw_tmp[1800]=vb_sw[1800]*totalTime[stanceLegID[0]]/maxTime;
        ab_sw_tmp[1800]=ab_sw[1800]*totalTime[stanceLegID[0]]/maxTime*totalTime[stanceLegID[0]]/maxTime;

        //memcpy(pb_sw,pb_sw_tmp,(1801)*sizeof(double));
        //memcpy(vb_sw,vb_sw_tmp,(1801)*sizeof(double));
        //memcpy(ab_sw,ab_sw_tmp,(1801)*sizeof(double));
    }

    printf("Iteration finished, iterCount=%d, maxTime=%.4f, timeArray:%.4f,%.4f\n\n",iterCount,maxTime,timeArray[0][stanceLegID[0]],timeArray[1][stanceLegID[0]]);

    gettimeofday(&tpend,NULL);
    tused=tpend.tv_sec-tpstart.tv_sec+(double)(tpend.tv_usec-tpstart.tv_usec)/1000000;
    printf("UsedTime:%f\n",tused);
}

void TimeOptimalGait::GetOptimalGait2s()
{
    double pEB[6] {0};
    double pEE[18] {0};
    double vEB[6] {0};
    double vEE[18] {0};
    double aEB[6] {0};
    double aEE[18] {0};

    for (int i=0;i<1801;i++)
    {
        s_w=0.1*i*PI/180;//degree to rad

        f_sw[0]=0;
        f_sw[1]=stepH*sin(s_w);
        f_sw[2]=stepD/2*cos(s_w);
        df_sw[0]=0;
        df_sw[1]=stepH*cos(s_w);
        df_sw[2]=-stepD/2*sin(s_w);
        ddf_sw[0]=0;
        ddf_sw[1]=-stepH*sin(s_w);
        ddf_sw[2]=-stepD/2*cos(s_w);

        pEB[2]=initPeb[2]+pb_sw_tmp[i];
        vEB[2]=vb_sw_tmp[i];
        aEB[2]=ab_sw_tmp[i];
        for(int j=0;j<3;j++)
        {
            //stanceLeg
            pEE[6*j+3]=initPee[6*j+3];
            pEE[6*j+4]=initPee[6*j+4];
            pEE[6*j+5]=initPee[6*j+5];

            //swingLeg
            pEE[6*j]=initPee[6*j]+f_sw[0];
            pEE[6*j+1]=initPee[6*j+1]+f_sw[1];
            pEE[6*j+2]=initPee[6*j+2]+f_sw[2];

            vEE[6*j]=df_sw[0]*real_ds[i][2*j];
            vEE[6*j+1]=df_sw[1]*real_ds[i][2*j];
            vEE[6*j+2]=df_sw[2]*real_ds[i][2*j];

            aEE[6*j]=ddf_sw[0]*real_ds[i][2*j]*real_ds[i][2*j]+df_sw[0]*real_dds[i][2*j];
            aEE[6*j+1]=ddf_sw[1]*real_ds[i][2*j]*real_ds[i][2*j]+df_sw[1]*real_dds[i][2*j];
            aEE[6*j+2]=ddf_sw[2]*real_ds[i][2*j]*real_ds[i][2*j]+df_sw[2]*real_dds[i][2*j];
        }

        rbt.SetPeb(pEB);
        rbt.SetVb(vEB);
        rbt.SetAb(aEB);
        rbt.SetPee(pEE);
        rbt.SetVee(vEE);
        rbt.SetAee(aEE);

        rbt.GetPee(*output_Pee+18*i,rbt.body());
        rbt.GetPin(*output_Pin+18*i);
        rbt.GetVin(*output_Vin+18*i);
        rbt.GetAin(*output_Ain+18*i);

        //printf("output:Pee=%.4f,Pin=%.4f,Vin=%.4f,Ain=%.4f\n",output_Pee[i][0],output_Pin[i][0],output_Vin[i][0],output_Ain[i][0]);
}

}

void TimeOptimalGait::GetOptimalGait2t()
{
    double * real_s=new double [6*maxTotalCount];
    double * real_Pee=new double [18*maxTotalCount];
    double * real_Pin=new double [18*maxTotalCount];

    for(int j=0;j<6;j++)
    {
        real_s[j]=0;
        for(int i=0;i<1801;i++)
        {
            real_ds_tmp[i][j]=real_ds[i][j]*totalTime[j]/maxTime;
        }
    }

    for (int i=1;i<maxTotalCount;i++)
    {
        double ds[6];
        for(int j=0;j<6;j++)
        {
            ds[j]=0.5*(real_ds_tmp[(int)(real_s[6*(i-1)+j]/(PI/1800))][j]+real_ds_tmp[(int)(real_s[6*(i-1)+j]/(PI/1800))+1][j]);
            real_s[6*i+j]=real_s[6*(i-1)+j]+ds[j]*0.001;
        }
    }

    for (int i=0;i<maxTotalCount;i++)
    {
        double pEB[6];
        memcpy(pEB,initPeb,6*sizeof(double));
        double pEE[18];
        double b_s=stepD/4-stepD/2*(real_s[6*i+stanceLegID[0]]/PI);
        pEB[2]=initPeb[2]+b_s;
        for(int j=0;j<3;j++)
        {
            //stanceLeg
            pEE[6*j+3]=initPee[6*j+3];
            pEE[6*j+4]=initPee[6*j+4];
            pEE[6*j+5]=initPee[6*j+5];

            //swingLeg
            double f_s[3];
            f_s[0]=0;
            f_s[1]=stepH*sin(real_s[6*i+2*j]);
            f_s[2]=stepD/2*cos(real_s[6*i+2*j]);
            pEE[6*j]=initPee[6*j]+f_s[0];
            pEE[6*j+1]=initPee[6*j+1]+f_s[1];
            pEE[6*j+2]=initPee[6*j+2]+f_s[2];
        }
        rbt.SetPeb(pEB);
        rbt.SetPee(pEE);

        rbt.GetPee(real_Pee+18*i,rbt.body());
        rbt.GetPin(real_Pin+18*i);
    }

    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_s.txt",real_s,maxTotalCount,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_Pee.txt",real_Pee,maxTotalCount,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_Pin.txt",real_Pin,maxTotalCount,18);
    delete [] real_s;
    delete [] real_Pee;
    delete [] real_Pin;
}

void TimeOptimalGait::OutputData()
{
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_PeeB.txt",*output_PeeB,1801,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_dsds1.txt",*output_dsds1,1801,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_dsds2.txt",*output_dsds2,1801,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_dsds.txt", *output_dsds, 1801,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_ds1.txt",*output_ds1,1801,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_ds2.txt",*output_ds2,1801,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_ds.txt", *output_ds, 1801,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_const1.txt",*output_const1,1801,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_const2.txt",*output_const2,1801,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_const.txt", *output_const, 1801,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_dds.txt",*output_dds,1801,18);

    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a2.txt",*output_a2,1801,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a1.txt",*output_a1,1801,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a0L.txt",*output_a0L,1801,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a0H.txt",*output_a0H,1801,18);

    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/max_ValueL.txt",*dds_lowBound,1801,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/min_ValueH.txt",*dds_upBound,1801,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_upBound_aLmt.txt",*ds_upBound_aLmt,1801,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_lowBound_aLmt.txt",*ds_lowBound_aLmt,1801,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_bound_vLmt.txt",*ds_upBound_vLmt,1801,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_lowBound.txt",*ds_lowBound,1801,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_upBound.txt",*ds_upBound,1801,6);
    //aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_max.txt",*real_ddsMax,1801,6);
    //aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_min.txt",*real_ddsMin,1801,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_forward.txt",*ds_forward,1801,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_backward.txt",*ds_backward,1801,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_forward.txt",*dds_forward,1801,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_backward.txt",*dds_backward,1801,6);

    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_ds.txt",*real_ds,1801,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_dds.txt",*real_dds,1801,6);

    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/timeArray.txt",*timeArray_tmp,1801,6);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/pb_sw.txt",pb_sw_tmp,1801,1);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/vb_sw.txt",vb_sw_tmp,1801,1);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ab_sw.txt",ab_sw_tmp,1801,1);

    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Pee.txt",*output_Pee,1801,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Pin.txt",*output_Pin,1801,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Vin.txt",*output_Vin,1801,18);
    aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Ain.txt",*output_Ain,1801,18);
}
