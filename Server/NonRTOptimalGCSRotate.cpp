#include "NonRTOptimalGCSRotate.h"

namespace TimeOptimal
{
    void NonRTOptimalGCSRotate::GetTimeOptimalGait(double step_length, double step_height, double step_alpha, double step_beta,
                                             double &duty_cycle, double acc_limit, double vel_limit, double *init_tippos,
                                             double *out_tippos, double &out_bodyvel, double &out_period)
    {
        stepD=step_length;
        stepH=step_height;
        stepAlpha=step_alpha;
        stepBeta=step_beta;
        dutyCycle=duty_cycle;
        aLmt=acc_limit;
        vLmt=vel_limit;
        memcpy(initPee,init_tippos,18*sizeof(double));

        GetOptimalDsByMajorIteration();
        GetOptimalGait2t(out_tippos,out_bodyvel,out_period);
        //OutputData();
        duty_cycle=dutyCycle;
    }

	/*获取约束方程中dds,ds^2,ds以及常数项的系数，针对身体运动受到的站立腿约束*/
    void NonRTOptimalGCSRotate::GetStanceLegParam(int count, int legID, double s)
    {
        double pEB[6] {0};
        double pEE[18] {0};

        double dJvi_x[9] {0};
        double dJvi_y[9] {0};
        double dJvi_z[9] {0};
        double dJvi_dot_C1[9] {0};//used in stanceLeg

        double param_dsds1[18] {0};
        double param_dsds2[18] {0};
        double param_dsds[18] {0};
        double param_const[18] {0};//for aLmt of swingLeg

        double b_sb_tmp[3];//{gama,z,x}
        double db_sb_tmp[3];
        double ddb_sb_tmp[3];

        b_sb_tmp[0]=stepBeta*s;
        b_sb_tmp[1]=0;
        b_sb_tmp[2]=0;

        db_sb_tmp[0]=stepBeta;
        db_sb_tmp[1]=0;
        db_sb_tmp[2]=0;

        ddb_sb_tmp[0]=0;
        ddb_sb_tmp[1]=0;
        ddb_sb_tmp[2]=0;

        pEB[0]=initPeb[0]+b_sb_tmp[2];
        pEB[2]=initPeb[2]+b_sb_tmp[1];
        pEB[3]=initPeb[3]+b_sb_tmp[0];
        if(legID%2==1)//135
        {
            if(count<count3)
            {
                pEE[3*legID]=initPee[3*legID]*cos(stepBeta/4)+initPee[3*legID+2]*sin(stepBeta/4);
                pEE[3*legID+1]=initPee[3*legID+1];
                pEE[3*legID+2]=-initPee[3*legID]*sin(stepBeta/4)+initPee[3*legID+2]*cos(stepBeta/4);
            }
            else
            {
                pEE[3*legID]=initPee[3*legID]*cos(5*stepBeta/4)+initPee[3*legID+2]*sin(5*stepBeta/4);
                pEE[3*legID+1]=initPee[3*legID+1];
                pEE[3*legID+2]=-initPee[3*legID]*sin(5*stepBeta/4)+initPee[3*legID+2]*cos(5*stepBeta/4);
            }
        }
        else//024
        {
            if(count<count1)
            {
                pEE[3*legID]=initPee[3*legID]*cos(-stepBeta/4)+initPee[3*legID+2]*sin(-stepBeta/4);
                pEE[3*legID+1]=initPee[3*legID+1];
                pEE[3*legID+2]=-initPee[3*legID]*sin(-stepBeta/4)+initPee[3*legID+2]*cos(-stepBeta/4);
            }
            else
            {
                pEE[3*legID]=initPee[3*legID]*cos(3*stepBeta/4)+initPee[3*legID+2]*sin(3*stepBeta/4);
                pEE[3*legID+1]=initPee[3*legID+1];
                pEE[3*legID+2]=-initPee[3*legID]*sin(3*stepBeta/4)+initPee[3*legID+2]*cos(3*stepBeta/4);
            }
        }

        rbt.SetPeb(pEB,"213");
        rbt.pLegs[legID]->SetPee(pEE+3*legID);

        rbt.pLegs[legID]->GetJvi(Jvi,rbt.body());
        rbt.pLegs[legID]->GetdJacOverPee(dJvi_x,dJvi_y,dJvi_z,"B");

        double C1_tmp[3] {0};//C1 in Thesis
        double C2_tmp[3] {0};//C2 in Thesis
        C1_tmp[0]=-db_sb_tmp[2]-db_sb_tmp[0]*(pEE[3*legID+2]-b_sb_tmp[1]);
        C1_tmp[2]=-db_sb_tmp[1]+db_sb_tmp[0]*(pEE[3*legID]  -b_sb_tmp[2]);
        C2_tmp[0]=-ddb_sb_tmp[0]*(pEE[3*legID+2]-b_sb_tmp[1])+db_sb_tmp[0]*db_sb_tmp[0]*(pEE[3*legID]  -b_sb_tmp[2])-2*db_sb_tmp[0]*C1_tmp[2]-ddb_sb_tmp[2];
        C2_tmp[2]= ddb_sb_tmp[0]*(pEE[3*legID]  -b_sb_tmp[2])+db_sb_tmp[0]*db_sb_tmp[0]*(pEE[3*legID+2]-b_sb_tmp[1])+2*db_sb_tmp[0]*C1_tmp[0]-ddb_sb_tmp[1];

        double C1[3] {0};
        double C2[3] {0};
        C1[0]=C1_tmp[0]*cos(b_sb_tmp[0])-C1_tmp[2]*sin(b_sb_tmp[0]);
        C1[2]=C1_tmp[0]*sin(b_sb_tmp[0])+C1_tmp[2]*cos(b_sb_tmp[0]);
        C2[0]=C2_tmp[0]*cos(b_sb_tmp[0])-C2_tmp[2]*sin(b_sb_tmp[0]);
        C2[2]=C2_tmp[0]*sin(b_sb_tmp[0])+C2_tmp[2]*cos(b_sb_tmp[0]);

        memcpy(*output_C+36*count+3*legID   ,C1,3*sizeof(double));
        memcpy(*output_C+36*count+18+3*legID,C2,3*sizeof(double));

        std::fill_n(dJvi_dot_C1,9,0);
        aris::dynamic::s_daxpy(9,C1[0],dJvi_x,1,dJvi_dot_C1,1);//for s_b
        aris::dynamic::s_daxpy(9,C1[1],dJvi_y,1,dJvi_dot_C1,1);
        aris::dynamic::s_daxpy(9,C1[2],dJvi_z,1,dJvi_dot_C1,1);

        std::fill_n(param_dds+3*legID,3,0);
        aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,C1,1,1,param_dds+3*legID,1);

        std::fill_n(param_dsds1+3*legID,3,0);
        aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,C2,1,1,param_dsds1+3*legID,1);//G*C2
        std::fill_n(param_dsds2+3*legID,3,0);
        aris::dynamic::s_dgemm(3,1,3,1,dJvi_dot_C1,3,C1,1,1,param_dsds2+3*legID,1);//H*C1
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
                    isParamddsExact0_body[count][3*legID+k]=1;
                    printf("WARNING!!! param_dds equals zero!!! StanceLeg : %d \n",count);
                }
            }
        }

        if(count!=-1)
        {
            for(int i=0;i<3;i++)
            {
                b_sb[count][i]=b_sb_tmp[i];
                db_sb[count][i]=db_sb_tmp[i];
                ddb_sb[count][i]=ddb_sb_tmp[i];
            }
            rbt.pLegs[legID]->GetPee(*output_PeeB+18*count+3*legID,rbt.body());
            memcpy( *output_dsds+18*count+3*legID,  param_dsds+3*legID, 3*sizeof(double));
            memcpy(  *output_dds+18*count+3*legID,   param_dds+3*legID, 3*sizeof(double));
            memcpy(   *output_a2+18*count+3*legID,    param_a2+3*legID, 3*sizeof(double));
            memcpy(  *output_a0L+18*count+3*legID,   param_a0L+3*legID, 3*sizeof(double));
            memcpy(  *output_a0H+18*count+3*legID,   param_a0H+3*legID, 3*sizeof(double));
        }
    }

	/*获取加速度约束下限值的最大值*/
    double NonRTOptimalGCSRotate::GetStanceMaxDec(int count, double ds)
    {
        double dec[18] {0};
        std::fill_n(dec,18,-1e6);

        for (int j=0;j<6;j++)
        {
            for (int k=0;k<3;k++)
            {
                if(isParamddsExact0_body[count][3*j+k]==-1)
                {
                    if(count>=count1 && count<count2)//stance 135
                    {
                        if(j%2==1)
                        dec[3*j+k]=output_a2[count][3*j+k]*ds*ds+output_a0L[count][3*j+k];
                    }
                    else if(count>=count3 && count<count4)//stance 024
                    {
                        if(j%2==0)
                        dec[3*j+k]=output_a2[count][3*j+k]*ds*ds+output_a0L[count][3*j+k];
                    }
                    else
                    {
                        dec[3*j+k]=output_a2[count][3*j+k]*ds*ds+output_a0L[count][3*j+k];
                    }
                }
            }
        }

        return *std::max_element(dec,dec+18);
    }

	/*获取加速度约束上限值的最小值*/
    double NonRTOptimalGCSRotate::GetStanceMinAcc(int count, double ds)
    {
        double acc[18] {0};
        std::fill_n(acc,18,1e6);

        for (int j=0;j<6;j++)
        {
            for (int k=0;k<3;k++)
            {
                if(isParamddsExact0_body[count][3*j+k]==-1)
                {
                    if(count>=count1 && count<count2)//stance 135
                    {
                        if(j%2==1)
                        acc[3*j+k]=output_a2[count][3*j+k]*ds*ds+output_a0H[count][3*j+k];
                    }
                    else if(count>=count3 && count<count4)//stance 024
                    {
                        if(j%2==0)
                        acc[3*j+k]=output_a2[count][3*j+k]*ds*ds+output_a0H[count][3*j+k];
                    }
                    else
                    {
                        acc[3*j+k]=output_a2[count][3*j+k]*ds*ds+output_a0H[count][3*j+k];
                    }
                }
            }
        }

        return *std::min_element(acc,acc+18);
    }

	/*遍历法，令加速度约束上限的最小值等于下限的最大值求得速度上限值，表现求为多条二次函数曲线的最小交点*/
    void NonRTOptimalGCSRotate::GetStanceDsBound(int count)
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
                            //printf("stance ds_upBound_aLmt=%.6f\n",ds_upBound_aLmt_body[count]);
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

	/*牛顿迭代法，令加速度约束上限的最小值等于下限的最大值求得速度上限值，表现求为多条二次函数曲线的最小交点*/
    void NonRTOptimalGCSRotate::GetStanceDsBoundByNewton(int count)
    {
        double ds0 {1e6};
        double value0;
        double ds1;
        double value_up;
        double value_low;
        double slope;

        int k {0};
        value0=GetStanceMinAcc(count,ds0)-GetStanceMaxDec(count,ds0);
        while (value0<0 && k<10)
        {
            k++;
            ds0/=10;
            value0=GetStanceMinAcc(count,ds0)-GetStanceMaxDec(count,ds0);
        }
        if(k==0)
        {
            printf("WARNING!!! ds0=%f too small, please set it larger and try again!",ds0);
            std::abort();
        }
        else
        {
            ds1=ds0*10;
        }
        //printf("ds0=%f,ds1=%f\n",ds0,ds1);

        k=0;
        value_low=-1;
        while(value_low<=0 && k<50)
        {
            k++;
            ds0=ds1;
            value_low=GetStanceMinAcc(count,ds0-1e-3)-GetStanceMaxDec(count,ds0-1e-3);
            value_up=GetStanceMinAcc(count,ds0+1e-3)-GetStanceMaxDec(count,ds0+1e-3);
    //        value0=GetStanceMinAcc(count,ds0)-GetStanceMaxDec(count,ds0);
    //        slope=(value_up-value_low)/2e-3;
    //        ds1=ds0-value0/slope;
            ds1=(value_up*(ds0-1e-3)-value_low*(ds0+1e-3))/(value_up-value_low);
        }
        //printf("First Newton Iteration of StanceLeg stops at k=%d, value_low=%f, value_up=%f, ds0=%f, ds1=%f\n",k,value_low,value_up,ds0,ds1);

        ds1=ds0;
        k=0;
        value_low=-1;
        while(value_low<=0 && k<50)
        {
            k++;
            ds0=ds1;
            value_low=GetStanceMinAcc(count,ds0-1e-6)-GetStanceMaxDec(count,ds0-1e-6);
            value_up=GetStanceMinAcc(count,ds0+1e-6)-GetStanceMaxDec(count,ds0+1e-6);
    //        value0=GetStanceMinAcc(count,ds0)-GetStanceMaxDec(count,ds0);
    //        slope=(value_up-value_low)/2e-6;
    //        ds1=ds0-value0/slope;
            ds1=(value_up*(ds0-1e-6)-value_low*(ds0+1e-6))/(value_up-value_low);
        }

        ds_upBound_aLmt_body[count]=ds1;
        ds_upBound_vLmt_body[count]=vLmt/(*std::max_element(abs_param_dds,abs_param_dds+18));
        ds_upBound_body[count]=std::min(ds_upBound_aLmt_body[count],ds_upBound_vLmt_body[count]);

        dds_lowBound_body[count]=GetStanceMaxDec(count,ds_upBound_body[count]);
        dds_upBound_body[count]=GetStanceMinAcc(count,ds_upBound_body[count]);

        //if(count>1000 && count<=1100)
        //printf("%d,Second Newton Iteration of StanceLeg stops at k=%d, value_low=%f,value_up=%f, ds0=%f, ds1=%f\n",count,k,value_low,value_up,ds0,ds1);
    }

	/*求二次函数曲线的所有交点，取其最小值，其中某些曲线可能没交点*/
    void NonRTOptimalGCSRotate::GetStanceDsBoundByFunc(int count)
    {
        double a2[18] {0};
        double a0H[18] {0};
        double a0L[18] {0};
        int num {0};
        double FuncA2;
        double FuncA0;
        double ds;
        double minDs {1e6};

        for (int j=0;j<6;j++)
        {
            for (int k=0;k<3;k++)
            {
                if(isParamddsExact0_body[count][3*j+k]==-1)
                {
                    if(count>=count1 && count<count2)//stance 135
                    {
                        if(j%2==1)
                        {
                            a2[num]=output_a2[count][3*j+k];
                            a0H[num]=output_a0H[count][3*j+k];
                            a0L[num]=output_a0L[count][3*j+k];
                            num++;
                        }
                    }
                    else if(count>=count3 && count<count4)//stance 024
                    {
                        if(j%2==0)
                        {
                            a2[num]=output_a2[count][3*j+k];
                            a0H[num]=output_a0H[count][3*j+k];
                            a0L[num]=output_a0L[count][3*j+k];
                            num++;
                        }
                    }
                    else
                    {
                        a2[num]=output_a2[count][3*j+k];
                        a0H[num]=output_a0H[count][3*j+k];
                        a0L[num]=output_a0L[count][3*j+k];
                        num++;
                    }
                }
            }
        }
        //printf("num:%d\n",num);

        for(int i=0;i<num;i++)
        {
            for(int j=0;j<num;j++)
            {
                if(i!=j)
                {
                    FuncA2=a2[i]-a2[j];
                    FuncA0=a0H[i]-a0L[j];
                    if(FuncA2==0 || FuncA0*FuncA2>0)
                    {
                        ds=-1;
                    }
                    else
                    {
                        ds=sqrt(-FuncA0/FuncA2);
                    }
                    if(minDs>ds && ds>0)
                    {
                        minDs=ds;
                    }
                }
            }
        }

        if(minDs==1e6)
        {
            printf("WARNING! Stance minDs too small! count=%d\n",count);
        }

        ds_upBound_aLmt_body[count]=minDs;
        ds_upBound_vLmt_body[count]=vLmt/(*std::max_element(abs_param_dds,abs_param_dds+18));
        ds_upBound_body[count]=std::min(ds_upBound_aLmt_body[count],ds_upBound_vLmt_body[count]);

        dds_lowBound_body[count]=GetStanceMaxDec(count,ds_upBound_body[count]);
        dds_upBound_body[count]=GetStanceMinAcc(count,ds_upBound_body[count]);
    }

	/*针对身体运动受到的站立腿约束，获取三类切换点*/
    void NonRTOptimalGCSRotate::GetStanceSwitchPoint()
    {
        double slopedsBound_for_body[2201] {0};
        double slopedsBound_back_body[2201] {0};
        double paramdds0Point_body[2201] {0};
        int paramdds0Count_body {0};

        double tangentPoint_body[2201] {0};
        int tangentCount_body {0};
        int switchScrewID_body_tmp[2201] {0};

        //initialize
        for(int i=0;i<count5+1;i++)
        {
            tangentPoint_body[i]=-1;
            paramdds0Point_body[i]=-1;
            switchScrewID_body_tmp[i]=-1;
        }
        tangentCount_body=0;
        paramdds0Count_body=0;

        //known switch point
        switchPoint_body[switchCount_body]=0;
        switchType_body[switchCount_body]='b';
        switchCount_body++;
    //    switchPoint_body[switchCount_body]=count5/2;
    //    switchType_body[switchCount_body]='b';
    //    switchCount_body++;
        switchPoint_body[switchCount_body]=count5;
        switchType_body[switchCount_body]='b';
        switchCount_body++;

        switchPoint_body[switchCount_body] = count1;//discontinuous form count1-1 to count1
        switchType_body[switchCount_body]='d';
        switchCount_body++;
        switchPoint_body[switchCount_body] = count2;
        switchType_body[switchCount_body]='d';
        switchCount_body++;
        switchPoint_body[switchCount_body] = count3;
        switchType_body[switchCount_body]='d';
        switchCount_body++;
        switchPoint_body[switchCount_body] = count4;
        switchType_body[switchCount_body]='d';
        switchCount_body++;

        //calculate the slope of ds_upBound, devided by t, not s
        slopedsBound_for_body[0]=slopedsBound_back_body[0]=(ds_upBound_body[1]*ds_upBound_body[1]-ds_upBound_body[0]*ds_upBound_body[0])/2/(s_b[1]-s_b[0]);
        for(int i=1;i<count5;i++)
        {
            slopedsBound_for_body[i]=(ds_upBound_body[i+1]*ds_upBound_body[i+1]-ds_upBound_body[i]*ds_upBound_body[i])/2/(s_b[i+1]-s_b[i]);
            slopedsBound_back_body[i]=(ds_upBound_body[i]*ds_upBound_body[i]-ds_upBound_body[i-1]*ds_upBound_body[i-1])/2/(s_b[i]-s_b[i-1]);
        }
        slopedsBound_for_body[count5]=slopedsBound_back_body[count5]=(ds_upBound_body[count5]*ds_upBound_body[count5]-ds_upBound_body[count5-1]*ds_upBound_body[count5-1])/2/(s_b[count5]-s_b[count5-1]);
    //    for(int i=0;i<count5+1;i++)
    //    {
    //        slopeDelta_body[i]=slopedsBound_for_body[i]-dds_upBound_body[i];
    //    }

        for(int i=0;i<count5;i++)
        {
            bool isParamdds0 {false};
            if(ds_upBound_vLmt_body[i]>=ds_upBound_aLmt_body[i])
            {
                int k_start;
                int k_end;

                if(i>=count1 && i<count2)//stance 135
                {
                    k_start=3;
                    k_end=6;
                }
                else if(i>=count3 && i<count4)//stance 024
                {
                    k_start=0;
                    k_end=3;
                }
                else// stance all
                {
                    k_start=0;
                    k_end=6;
                }

                for(int j=0;j<3;j++)
                {
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

            if(slopedsBound_for_body[i]<=dds_upBound_body[i] && slopedsBound_for_body[i+1]>dds_upBound_body[i+1]
                    && isParamdds0==false && i!=count1-1 && i!=count2-1  && i!=count3-1  && i!=count4-1)
            {
                tangentPoint_body[tangentCount_body]=i
                        +fabs(slopedsBound_for_body[i]-dds_upBound_body[i])/(fabs(slopedsBound_for_body[i]-dds_upBound_body[i])+fabs(slopedsBound_for_body[i+1]-dds_upBound_body[i+1]));
                tangentCount_body++;
            }
            //special switch point on vel bound curve
            if(slopedsBound_back_body[i]<=dds_lowBound_body[i] && slopedsBound_back_body[i+1]>dds_lowBound_body[i+1]
                    && isParamdds0==false && i!=count1-1 && i!=count2-1 && i!=count3-1 && i!=count4-1)
            {
                tangentPoint_body[tangentCount_body]=i
                        +fabs(slopedsBound_back_body[i]-dds_lowBound_body[i])/(fabs(slopedsBound_back_body[i]-dds_lowBound_body[i])+fabs(slopedsBound_back_body[i+1]-dds_lowBound_body[i+1]));
                tangentCount_body++;
            }
        }

//        printf("StanceLeg Tangent Switch Point:");
//        for(int i=0;i<tangentCount_body+1;i++)
//        {
//            printf("%.2f,",tangentPoint_body[i]);
//        }
//        printf("\n");

//        printf("StanceLeg ZeroInertia Switch Point:");
//        for(int i=0;i<paramdds0Count_body+1;i++)
//        {
//            printf("%.2f,",paramdds0Point_body[i]);
//        }
//        printf("\n");

        //merge tangentPoint & paramdds0Point into switchPoint
        for(int i=0;i<tangentCount_body;i++)
        {
            switchPoint_body[i+switchCount_body]=tangentPoint_body[i];
            switchType_body[i+switchCount_body]='t';
            //switchScrewID_body[i+switchCount_body]=switchScrewID_body_tmp[i];
        }
        switchCount_body+=tangentCount_body;
        for(int i=0;i<paramdds0Count_body;i++)
        {
            switchPoint_body[i+switchCount_body]=paramdds0Point_body[i];
            switchType_body[i+switchCount_body]='z';
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

                    auto tmp3=switchType_body[i];
                    switchType_body[i]=switchType_body[j];
                    switchType_body[j]=tmp3;
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
            printf("%.4f,",s_b[(int)switchPoint_body[i]]+(switchPoint_body[i]-(int)switchPoint_body[i])*(s_b[(int)switchPoint_body[i]+1]-s_b[(int)switchPoint_body[i]]));
        }
        printf("\n");
        printf("StanceLeg Switch Type:");
        for(int i=0;i<switchCount_body+1;i++)
        {
            printf("%c,",switchType_body[i]);
        }
        printf("\n");
    }

	/*获取切换点前后向前积分和向后积分的起始点，直接从切换点积分在某些时候会因离散积分误差而无法通过*/
    void NonRTOptimalGCSRotate::GetStanceTwoPointAtSwitch(double *lowPoint, double *upPoint)
    {
        lowPoint[0]=-1;
        lowPoint[switchCount_body-1]=ds_upBound_body[count5];
        upPoint[0]=ds_upBound_body[0];
        upPoint[switchCount_body-1]=-1;

        for(int i=1;i<switchCount_body-1;i++)
        {
            if(switchType_body[i]=='t')
            {
                int num=(int)switchPoint_body[i];
                upPoint[i]=ds_upBound_body[num+1];
                lowPoint[i]=ds_upBound_body[num];
            }
            else if(switchType_body[i]=='d')
            {
                int num=(int)switchPoint_body[i];
                upPoint[i]=lowPoint[i]=std::min(ds_upBound_body[num],ds_upBound_body[num-1]);
            }
            else if(switchType_body[i]=='z')
            {
                int num=round(switchPoint_body[i]);
                auto tmp=std::min(ds_upBound_body[num-1],ds_upBound_body[num]);
                upPoint[i]=lowPoint[i]=std::min(tmp,ds_upBound_body[num+1]);
            }
        }
    }

	/*从切换点开始前向后向积分获取最佳ds曲线*/
    void NonRTOptimalGCSRotate::GetStanceOptimalDsBySwitchPoint()
    {
        bool stopFlag {false};
        int forwardEnd_s {-1};
        double forwardEnd_ds {0};
        double *lowPoint=new double [switchCount_body];
        double *upPoint=new double [switchCount_body];
        GetStanceTwoPointAtSwitch(lowPoint,upPoint);

//        printf("lowPoint:");
//        for(int i=0;i<switchCount_body;i++)
//        {
//            printf("%.2f,",lowPoint[i]);
//        }
//        printf("\n");
//        printf("upPoint:");
//        for(int i=0;i<switchCount_body;i++)
//        {
//            printf("%.2f,",upPoint[i]);
//        }
//        printf("\n");

        for(int i=0;i<count5+1;i++)
        {
            real_ds_body[i]=ds_upBound_body[i];
            real_dds_body[i]=dds_upBound_body[i];
        }
        for(int m=0;m<switchCount_body;m++)
        {
            //start of backward
            bool ignoreBackward {false};
            int k_st=(int)switchPoint_body[m];
            int k_st_start=k_st;
            if(switchType_body[m]=='z')
            {
                k_st_start=k_st=round(switchPoint_body[m])-1;
            }
            else if(switchType_body[m]=='d')
            {
                if(ds_upBound_body[k_st_start]>ds_upBound_body[k_st_start-1])
                {
                    k_st_start=k_st=(int)switchPoint_body[m]-1;
                }
            }

            if(k_st_start==0)
            {
                ignoreBackward=true;
            }
            else if(k_st_start<forwardEnd_s)
            {
                ignoreBackward=true;
                //printf("StanceLeg backward start at a passed point, quit switchPoint %.1f\n",switchPoint_body[m]);
            }
            else if(k_st_start==forwardEnd_s && lowPoint[m]>forwardEnd_ds)
            {
                ignoreBackward=true;
                if(switchType_body[m]=='z')
                {
                    real_dds_body[k_st_start+1]=0;
                    real_ds_body[k_st_start+1]=lowPoint[m];
                }
            }

            if(ignoreBackward==false)
            {
                if(switchType_body[m]=='z')
                {
                    real_dds_body[k_st_start+1]=0;
                    real_ds_body[k_st_start+1]=lowPoint[m];
                }
                stopFlag=false;
                ds_backward_body[k_st]=lowPoint[m];
                while(stopFlag==false)
                {
                    dds_backward_body[k_st]=GetStanceMaxDec(k_st,ds_backward_body[k_st]);
                    ds_backward_body[k_st-1]=sqrt(ds_backward_body[k_st]*ds_backward_body[k_st]-2*dds_backward_body[k_st]*(s_b[k_st]-s_b[k_st-1]));

                    if(ds_backward_body[k_st-1]>ds_upBound_body[k_st-1])
                    {
                        real_dds_body[k_st-1]=(real_ds_body[k_st-1]-ds_backward_body[k_st])*(real_ds_body[k_st-1]+ds_backward_body[k_st])/2/(s_b[k_st]-s_b[k_st-1]);
                        for(int i=k_st;i<k_st_start+1;i++)
                        {
                            real_ds_body[i]=ds_backward_body[i];
                            real_dds_body[i]=dds_backward_body[i];
                        }
                        stopFlag=true;
                        printf("StanceLeg backward touching upBound at %d, from switchPoint %.1f\n",k_st-1,switchPoint_body[m]);
                    }
                    else if(k_st==1)
                    {
                        dds_backward_body[k_st-1]=GetStanceMaxDec(k_st-1,ds_backward_body[k_st-1]);
                        for(int i=k_st-1;i<k_st_start+1;i++)
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
                        for(int i=k_st;i<k_st_start+1;i++)
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
    //            //if(ds_backward_body[k_st-1]>ds_upBound_body[k_st-1] || round(switchPoint_body[m])==2200)//end of backward
    //            if(k_st_start==2200)
    //            {
    //                continue;
    //            }
            }

            //start of forward
            bool ignoreForward {false};
            k_st_start=k_st=(int)switchPoint_body[m];
            if(switchType_body[m]=='t')
            {
                k_st_start=k_st=(int)switchPoint_body[m]+1;
            }
            else if(switchType_body[m]=='z')
            {
                k_st_start=k_st=round(switchPoint_body[m])+1;
            }
            else if(switchType_body[m]=='d')
            {
                if(ds_upBound_body[k_st_start]>ds_upBound_body[k_st_start-1])
                {
                    k_st_start=k_st=(int)switchPoint_body[m]-1;
                }
            }

            if(k_st_start==count5+1 || k_st_start==count5)
            {
                ignoreForward=true;
            }
            else if(k_st_start<forwardEnd_s)
            {
                ignoreForward=true;
                //printf("StanceLeg forward start at a passed point, quit switchPoint %.1f\n",switchPoint_body[m]);
            }
            else if(k_st_start==forwardEnd_s)
            {
                if(upPoint[m]>forwardEnd_ds)
                {
                    printf("How possible! StanceLeg forward curve should not stop here!\n");
                }
            }

            if(ignoreForward==false)
            {
                stopFlag=false;
                ds_forward_body[k_st]=upPoint[m];
                while(stopFlag==false)
                {
                    dds_forward_body[k_st]=GetStanceMinAcc(k_st,ds_forward_body[k_st]);
                    ds_forward_body[k_st+1]=sqrt(ds_forward_body[k_st]*ds_forward_body[k_st]+2*dds_forward_body[k_st]*(s_b[k_st+1]-s_b[k_st]));

                    if(ds_forward_body[k_st+1]>ds_upBound_body[k_st+1] || k_st==count5-1)
                    {
                        forwardEnd_s=k_st;
                        forwardEnd_ds=ds_forward_body[k_st];
                        for(int i=k_st_start;i<k_st+1;i++)
                        {
                            real_ds_body[i]=ds_forward_body[i];
                            real_dds_body[i]=dds_forward_body[i];
                        }
                        if(k_st==count5-1)
                        {
                            if(ds_forward_body[k_st+1]>ds_upBound_body[k_st+1])
                            {
                                real_ds_body[count5]=ds_upBound_body[count5];
                                real_dds_body[count5]=dds_upBound_body[count5];
                            }
                            else
                            {
                                real_ds_body[count5]=ds_forward_body[count5];
                                real_dds_body[count5]=GetStanceMinAcc(count5,ds_forward_body[count5]);
                                forwardEnd_s=count5;
                                forwardEnd_ds=ds_forward_body[count5];
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
        }

        for(int i=0;i<count5+1;i++)
        {
            real_ddsMax_body[i]=GetStanceMinAcc(i,real_ds_body[i]);
            real_ddsMin_body[i]=GetStanceMaxDec(i,real_ds_body[i]);
        }

        delete [] lowPoint;
        delete [] upPoint;
    }

	/*因离散积分误差存在，最优曲线的起始点和终止点值有差异，但行走过程要求这两个值相等，因此取其较大侧修正，并积分*/
    void NonRTOptimalGCSRotate::ApplyStanceExtraItegration()
    {
        bool stopFlag {false};
        int k_st {0};
        int k_st_start {0};
        double real_ds_body_tmp {0};
        if(real_ds_body[0]>real_ds_body[count5])
        {
            //forward
            real_ds_body[0]=real_ds_body[count5];
            k_st_start=k_st=0;
            stopFlag=false;
            while(stopFlag==false)
            {
                real_dds_body[k_st]=GetStanceMinAcc(k_st,real_ds_body[k_st]);
                real_ds_body_tmp=sqrt(real_ds_body[k_st]*real_ds_body[k_st]+2*real_dds_body[k_st]*(s_b[k_st+1]-s_b[k_st]));

                if(real_ds_body_tmp>real_ds_body[k_st+1])
                {
                    stopFlag=true;
                }
                else
                {
                    real_ds_body[k_st+1]=real_ds_body_tmp;
                    k_st++;
                }
            }
        }
        else if(real_ds_body[0]<real_ds_body[count5])
        {
            //backward
            real_ds_body[count5]=real_ds_body[0];
            k_st_start=k_st=count5;
            stopFlag=false;
            while(stopFlag==false)
            {
                real_dds_body[k_st]=GetStanceMaxDec(k_st,real_ds_body[k_st]);
                real_ds_body_tmp=sqrt(real_ds_body[k_st]*real_ds_body[k_st]-2*real_dds_body[k_st]*(s_b[k_st]-s_b[k_st-1]));

                if(real_ds_body_tmp>real_ds_body[k_st-1])
                {
                    stopFlag=true;
                }
                else
                {
                    real_ds_body[k_st-1]=real_ds_body_tmp;
                    k_st--;
                }
            }
        }
        else
        {
            printf("Amazing!!! Ds[start] equal Ds[end].No need to apply extra itegration.\n");
        }

    }

	/*不使用切换点，从起点和终点开始尝试积分，超出ds上限后退回再尝试，如此实现求解。该方法开发于不了解切换点时，计算耗时大于切换点且计算误差大，不建议使用*/
    void NonRTOptimalGCSRotate::GetStanceOptimalDsByDirectNI()
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
        ki_back=count5;
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
        double min_dist[2201] {0};
        std::fill_n(min_dist,count5+1,1);
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
                            <=0 || ki_for==count5-1)
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
            if(ki_for==count5)
            {
                stopFlag=true;
                for(int i=0;i<count5+1;i++)
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
                for(int i=ki_for;i<count5+1;i++)
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

	/*调用前面的函数，计算给定s_b1,s_b2,s_b3,s_b4时的最优曲线*/
    void NonRTOptimalGCSRotate::GetStanceOptimalDsAtSb(double s_b1, double s_b2, double s_b3, double s_b4)
    {
        std::fill_n(ds_backward_body,count5+1,0);
        std::fill_n(ds_forward_body,count5+1,0);
        Tstep=0;
        switchCount_body=0;
        std::fill_n(switchPoint_body,count5+1,-1);
        std::fill_n(switchType_body,count5+1,'0');
        std::fill_n(switchScrewID_body,count5+1,-1);
        std::fill_n(*isParamddsExact0_body,(count5+1)*18,-1);

        timeval tpstart,tpend;
        float tused;
        gettimeofday(&tpstart,NULL);

        for (int i=0;i<count5+1;i++)
        {
            std::fill_n(abs_param_dds,18,0);
            if(i<count1)
            {
                s_b[i]=s_b1/count1*i;
                for(int j=0;j<6;j++)
                {
                    GetStanceLegParam(i,j,s_b[i]);
                }
            }
            else if(i<count2)
            {
                s_b[i]=s_b1+(s_b2-s_b1)/(count2-count1)*(i-count1);
                for(int j=0;j<3;j++)
                {
                    GetStanceLegParam(i,2*j+1,s_b[i]);//135
                }
            }
            else if(i<count5/2)
            {
                s_b[i]=s_b2+(0.5-s_b2)/(count5/2-count2)*(i-count2);
                for(int j=0;j<6;j++)
                {
                    GetStanceLegParam(i,j,s_b[i]);
                }
            }
            else if(i<count3)
            {
                s_b[i]=0.5+(s_b3-0.5)/(count3-count5/2)*(i-count5/2);
                for(int j=0;j<6;j++)
                {
                    GetStanceLegParam(i,j,s_b[i]);
                }
            }
            else if(i<count4)
            {
                s_b[i]=s_b3+(s_b4-s_b3)/(count4-count3)*(i-count3);
                for(int j=0;j<3;j++)
                {
                    GetStanceLegParam(i,2*j,s_b[i]);//024
                }
            }
            else
            {
                s_b[i]=s_b4+(1-s_b4)/(count5-count4)*(i-count4);
                for(int j=0;j<6;j++)
                {
                    GetStanceLegParam(i,j,s_b[i]);
                }
            }

            //GetStanceDsBound(i);
            //GetStanceDsBoundByNewton(i);
            GetStanceDsBoundByFunc(i);
        }
        //output data to analyse when error
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_PeeB.txt",*output_PeeB,2201,18);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_dsds.txt", *output_dsds, 2201,18);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_ds.txt", *output_ds, 2201,18);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_const.txt", *output_const, 2201,18);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_dds.txt",*output_dds,2201,18);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a2.txt",*output_a2,2201,18);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a1.txt",*output_a1,2201,18);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a0L.txt",*output_a0L,2201,18);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a0H.txt",*output_a0H,2201,18);

        gettimeofday(&tpend,NULL);
        tused=tpend.tv_sec-tpstart.tv_sec+(double)(tpend.tv_usec-tpstart.tv_usec)/1000000;
        printf("UsedTime param:%f\n",tused);

        GetStanceSwitchPoint();
        GetStanceOptimalDsBySwitchPoint();
        ApplyStanceExtraItegration();
        //GetStanceOptimalDsByDirectNI();

        for(int i=0;i<count5+1;i++)
        {
            if(i!=0)
            {
                Tstep+=2*(s_b[i]-s_b[i-1])/(real_ds_body[i-1]+real_ds_body[i]);
                timeArray_body[i]=Tstep;
            }
        }

    }

	/*因s_b1,s_b2,s_b3,s_b4未知，需要迭代求解修正。该函数使用自己的算法，不能保证收敛性。实测收敛*/
    void NonRTOptimalGCSRotate::GetStanceOptimalDsByMinorIteration()
    {
        double t1 {0};
        double t2 {0};
        double t3 {0};
        double t4 {0};
        int k {0};
    //    s_b1=0.5*dutyCycle-0.25;/*(dutyCycle-0.5)/2*/
    //    s_b2=0.75-0.5*dutyCycle;/*(dutyCycle-0.5)/2+(1-dutyCycle)*/
    //    s_b3=0.5*dutyCycle+0.25;/*(dutyCycle-0.5)/2+0.5*/
    //    s_b4=1.25-0.5*dutyCycle;/*1-(dutyCycle-0.5)/2*/
        s_b1=0;
        s_b2=0.5;
        s_b3=0.5;
        s_b4=1;

        while((fabs(t1-(0.5*dutyCycle-0.25)*Tstep)>1e-6 || fabs(t2-(0.75-0.5*dutyCycle)*Tstep)>1e-6 ||
               fabs(t3-(0.5*dutyCycle+0.25)*Tstep)>1e-6 || fabs(t4-(1.25-0.5*dutyCycle)*Tstep)>1e-6) && k<30)
        {
            if(k!=0)
            {
                for(int i=0;i<count5;i++)
                {
                    if(timeArray_body[i]<=(0.5*dutyCycle-0.25)*Tstep && timeArray_body[i+1]>(0.5*dutyCycle-0.25)*Tstep)
                    {
                        s_b1=s_b[i]+(s_b[i+1]-s_b[i])/(timeArray_body[i+1]-timeArray_body[i])*((0.5*dutyCycle-0.25)*Tstep-timeArray_body[i]);
                    }
                    if(timeArray_body[i]<=(0.75-0.5*dutyCycle)*Tstep && timeArray_body[i+1]>(0.75-0.5*dutyCycle)*Tstep)
                    {
                        s_b2=s_b[i]+(s_b[i+1]-s_b[i])/(timeArray_body[i+1]-timeArray_body[i])*((0.75-0.5*dutyCycle)*Tstep-timeArray_body[i]);
                    }
                    if(timeArray_body[i]<=(0.5*dutyCycle+0.25)*Tstep && timeArray_body[i+1]>(0.5*dutyCycle+0.25)*Tstep)
                    {
                        s_b3=s_b[i]+(s_b[i+1]-s_b[i])/(timeArray_body[i+1]-timeArray_body[i])*((0.5*dutyCycle+0.25)*Tstep-timeArray_body[i]);
                    }
                    if(timeArray_body[i]<=(1.25-0.5*dutyCycle)*Tstep && timeArray_body[i+1]>(1.25-0.5*dutyCycle)*Tstep)
                    {
                        s_b4=s_b[i]+(s_b[i+1]-s_b[i])/(timeArray_body[i+1]-timeArray_body[i])*((1.25-0.5*dutyCycle)*Tstep-timeArray_body[i]);
                    }
                }
            }

            GetStanceOptimalDsAtSb(s_b1, s_b2, s_b3, s_b4);

            t1=timeArray_body[count1];
            t2=timeArray_body[count2];
            t3=timeArray_body[count3];
            t4=timeArray_body[count4];
            k++;

            printf("Tsb1=%.6f,Tsb2=%.6f,Tsb3=%.6f,Tsb4=%.6f,Tstep=%.6f\n",t1,t2,t3,t4,Tstep);
        }

        printf("Minor iteration count: %d, s_b1=%.4f, s_b2=%.4f, s_b3=%.4f, s_b4=%.4f\n",k,s_b1,s_b2,s_b3,s_b4);
        //aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/s_body.txt",s_b,count5+1,1);
    }

	/*因s_b1,s_b2,s_b3,s_b4未知，需要迭代求解修正。该函数使用拟牛顿法，能保证收敛性。*/
    void NonRTOptimalGCSRotate::GetStanceOptimalDsByNewton()
    {
        double t1 {0};
        double t2 {0};
        double t3 {0};
        double t4 {0};
        double D1 {1};
        double D2 {1};
        double D3 {1};
        double D4 {1};

        int k {0};
        s_b1=0.5*dutyCycle-0.25;/*(dutyCycle-0.5)/2*/
        s_b2=0.75-0.5*dutyCycle;/*(dutyCycle-0.5)/2+(1-dutyCycle)*/
        s_b3=0.5*dutyCycle+0.25;/*(dutyCycle-0.5)/2+0.5*/
        s_b4=1.25-0.5*dutyCycle;/*1-(dutyCycle-0.5)/2*/
//        s_b1=0;
//        s_b2=0.5;
//        s_b3=0.5;
//        s_b4=1;
        GetStanceOptimalDsAtSb(s_b1, s_b2, s_b3, s_b4);
        t1=timeArray_body[count1]-(0.5*dutyCycle-0.25)*Tstep;
        t2=timeArray_body[count2]-(0.75-0.5*dutyCycle)*Tstep;
        t3=timeArray_body[count3]-(0.5*dutyCycle+0.25)*Tstep;
        t4=timeArray_body[count4]-(1.25-0.5*dutyCycle)*Tstep;

        while((fabs(t1)>1e-6 || fabs(t2)>1e-6 || fabs(t3)>1e-6 || fabs(t4)>1e-6) && k<30)
        {
            timeval tpstart,tpend;
            float tused;
            gettimeofday(&tpstart,NULL);

            double d1=-D1*t1;
            double d2=-D2*t2;
            double d3=-D3*t3;
            double d4=-D4*t4;
            s_b1+=d1;
            s_b2+=d2;
            s_b3+=d3;
            s_b4+=d4;
            printf("\n\nsb:%f,%f,%f,%f\n",s_b1,s_b2,s_b3,s_b4);

            GetStanceOptimalDsAtSb(s_b1, s_b2, s_b3, s_b4);
            double t1_tmp=timeArray_body[count1]-(0.5*dutyCycle-0.25)*Tstep;
            double t2_tmp=timeArray_body[count2]-(0.75-0.5*dutyCycle)*Tstep;
            double t3_tmp=timeArray_body[count3]-(0.5*dutyCycle+0.25)*Tstep;
            double t4_tmp=timeArray_body[count4]-(1.25-0.5*dutyCycle)*Tstep;

            D1=d1/(t1_tmp-t1);
            D2=d2/(t2_tmp-t2);
            D3=d3/(t3_tmp-t3);
            D4=d4/(t4_tmp-t4);
            t1=t1_tmp;
            t2=t2_tmp;
            t3=t3_tmp;
            t4=t4_tmp;
            k++;

            printf("Tsb1=%.7f,Tsb2=%.7f,Tsb3=%.7f,Tsb4=%.7f,Tstep=%.7f\n",t1,t2,t3,t4,Tstep);

            gettimeofday(&tpend,NULL);
            tused=tpend.tv_sec-tpstart.tv_sec+(double)(tpend.tv_usec-tpstart.tv_usec)/1000000;
            printf("UsedTime body:%f\n",tused);
        }

        printf("Quasi-Newton Method count: %d, s_b1=%.4f, s_b2=%.4f, s_b3=%.4f, s_b4=%.4f\n",k,s_b1,s_b2,s_b3,s_b4);
        //output data to analyse when error
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_upBound_aLmt_body.txt",ds_upBound_aLmt_body,2201,1);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_upBound_vLmt_body.txt",ds_upBound_vLmt_body,2201,1);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_upBound_body.txt",dds_upBound_body,2201,1);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_upBound_body.txt",dds_lowBound_body,2201,1);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_forward_body.txt",ds_forward_body,2201,1);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_backward_body.txt",ds_backward_body,2201,1);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_forward_body.txt",dds_forward_body,2201,1);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_backward_body.txt",dds_backward_body,2201,1);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_ds_body.txt",real_ds_body,2201,1);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_dds_body.txt",real_dds_body,2201,1);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_ddsMax_body.txt",real_ddsMax_body,2201,1);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_ddsMin_body.txt",real_ddsMin_body,2201,1);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/timeArray_body_tmp.txt",timeArray_body_tmp,2201,1);
    }

	/*获取约束方程中dds,ds^2,ds以及常数项的系数，针对已求得身体运动后摆动腿受到的约束*/
    void NonRTOptimalGCSRotate::GetSwingLegParam(int count, int legID, double sw, double *sb_mtx)
    {
        double pEB[6] {0};
        double pEE[18] {0};

        double dJvi_x[9] {0};
        double dJvi_y[9] {0};
        double dJvi_z[9] {0};
        double dJvi_dot_df[9] {0};//used in swingLeg
        double dJvi_dot_C1[9] {0};//used in swingLeg

        double param_dsds1[18] {0};
        double param_dsds2[18] {0};
        double param_ds1[18] {0};
        double param_ds2[18] {0};
        double param_const1[18] {0};
        double param_const2[18] {0};
        double param_const3[18] {0};
        double param_dsds[18] {0};
        double param_ds[18] {0};//for aLmt of swingLeg
        double param_const[18] {0};//for aLmt of swingLeg

        double ratio_f[3] {0};
        if(legID%2==1)
        {
            ratio_f[0]=initPee[3*legID]*(cos(5*stepBeta/4)-cos(stepBeta/4))+initPee[3*legID+2]*(sin(5*stepBeta/4)-sin(stepBeta/4));
            ratio_f[2]=-initPee[3*legID]*(sin(5*stepBeta/4)-sin(stepBeta/4))+initPee[3*legID+2]*(cos(5*stepBeta/4)-cos(stepBeta/4));
        }
        else
        {
            ratio_f[0]=initPee[3*legID]*(cos(3*stepBeta/4)-cos(-stepBeta/4))+initPee[3*legID+2]*(sin(3*stepBeta/4)-sin(-stepBeta/4));
            ratio_f[2]=-initPee[3*legID]*(sin(3*stepBeta/4)-sin(-stepBeta/4))+initPee[3*legID+2]*(cos(3*stepBeta/4)-cos(-stepBeta/4));
        }

        f_sw[0]=-ratio_f[0]/2*(cos(sw)-1);
        f_sw[1]=stepH*sin(sw);
        f_sw[2]=-ratio_f[2]/2*(cos(sw)-1);

        df_sw[0]=ratio_f[0]/2*sin(sw);
        df_sw[1]=stepH*cos(sw);
        df_sw[2]=ratio_f[2]/2*sin(sw);

        ddf_sw[0]=ratio_f[0]/2*cos(sw);
        ddf_sw[1]=-stepH*sin(sw);
        ddf_sw[2]=ratio_f[2]/2*cos(sw);

        //b_sb:{gama,z,x}
        pEB[0]=initPeb[0]+b_sb[count][2];
        pEB[2]=initPeb[2]+b_sb[count][1];
        pEB[3]=initPeb[3]+b_sb[count][0];
        if(legID%2==1)//135
        {
            pEE[3*legID]=initPee[3*legID]*cos(stepBeta/4)+initPee[3*legID+2]*sin(stepBeta/4)+f_sw[0];
            pEE[3*legID+1]=initPee[3*legID+1]+f_sw[1];
            pEE[3*legID+2]=-initPee[3*legID]*sin(stepBeta/4)+initPee[3*legID+2]*cos(stepBeta/4)+f_sw[2];
        }
        else//024
        {
            pEE[3*legID]=initPee[3*legID]*cos(-stepBeta/4)+initPee[3*legID+2]*sin(-stepBeta/4)+f_sw[0];
            pEE[3*legID+1]=initPee[3*legID+1]+f_sw[1];
            pEE[3*legID+2]=-initPee[3*legID]*sin(-stepBeta/4)+initPee[3*legID+2]*cos(-stepBeta/4)+f_sw[2];
        }
        rbt.SetPeb(pEB,"213");
        rbt.pLegs[legID]->SetPee(pEE+3*legID);

        memcpy(*output_PebPee+24*count,pEB,6*sizeof(double));
        memcpy(*output_PebPee+24*count+6+3*legID,pEE+3*legID,3*sizeof(double));

        rbt.pLegs[legID]->GetJvi(Jvi,rbt.body());
        rbt.pLegs[legID]->GetdJacOverPee(dJvi_x,dJvi_y,dJvi_z,"B");

        //b_sb,db_sb,ddb_sb:{gama,z,x}
        double C1_tmp[3] {0};//C1 in Thesis
        double C2_tmp[3] {0};//C2 in Thesis
        C1_tmp[0]=-db_sb[count][2]-db_sb[count][0]*(pEE[3*legID+2]-b_sb[count][1]);
        C1_tmp[2]=-db_sb[count][1]+db_sb[count][0]*(pEE[3*legID]-b_sb[count][2]);
        C2_tmp[0]=-ddb_sb[count][0]*(pEE[3*legID+2]-b_sb[count][1])+db_sb[count][0]*db_sb[count][0]*(pEE[3*legID]  -b_sb[count][2])-2*db_sb[count][0]*C1_tmp[2]-ddb_sb[count][2];
        C2_tmp[2]= ddb_sb[count][0]*(pEE[3*legID]  -b_sb[count][2])+db_sb[count][0]*db_sb[count][0]*(pEE[3*legID+2]-b_sb[count][1])+2*db_sb[count][0]*C1_tmp[0]-ddb_sb[count][1];

        double C1[3] {0};
        double C2[3] {0};
        C1[0]=C1_tmp[0]*cos(b_sb[count][0])-C1_tmp[2]*sin(b_sb[count][0]);
        C1[2]=C1_tmp[0]*sin(b_sb[count][0])+C1_tmp[2]*cos(b_sb[count][0]);
        C2[0]=C2_tmp[0]*cos(b_sb[count][0])-C2_tmp[2]*sin(b_sb[count][0]);
        C2[2]=C2_tmp[0]*sin(b_sb[count][0])+C2_tmp[2]*cos(b_sb[count][0]);

        memcpy(*output_C+36*count+3*legID   ,C1,3*sizeof(double));
        memcpy(*output_C+36*count+18+3*legID,C2,3*sizeof(double));

        double dgama_cross_df_tmp[3] {db_sb[count][0]*df_sw[2],0,-db_sb[count][0]*df_sw[0]};
        double dgama_cross_df[3] {0};
        dgama_cross_df[0]=dgama_cross_df_tmp[0]*cos(b_sb[count][0])-dgama_cross_df_tmp[2]*sin(b_sb[count][0]);
        dgama_cross_df[2]=dgama_cross_df_tmp[0]*sin(b_sb[count][0])+dgama_cross_df_tmp[2]*cos(b_sb[count][0]);

        double df_sw_inB[3] {0};
        df_sw_inB[0]=df_sw[0]*cos(b_sb[count][0])-df_sw[2]*sin(b_sb[count][0]);
        df_sw_inB[1]=df_sw[1];
        df_sw_inB[2]=df_sw[0]*sin(b_sb[count][0])+df_sw[2]*cos(b_sb[count][0]);

        double ddf_sw_inB[3] {0};
        ddf_sw_inB[0]=ddf_sw[0]*cos(b_sb[count][0])-ddf_sw[2]*sin(b_sb[count][0]);
        ddf_sw_inB[1]=ddf_sw[1];
        ddf_sw_inB[2]=ddf_sw[0]*sin(b_sb[count][0])+ddf_sw[2]*cos(b_sb[count][0]);

        std::fill_n(dJvi_dot_df,9,0);
        aris::dynamic::s_daxpy(9,df_sw_inB[0],dJvi_x,1,dJvi_dot_df,1);
        aris::dynamic::s_daxpy(9,df_sw_inB[1],dJvi_y,1,dJvi_dot_df,1);
        aris::dynamic::s_daxpy(9,df_sw_inB[2],dJvi_z,1,dJvi_dot_df,1);

        std::fill_n(dJvi_dot_C1,9,0);
        aris::dynamic::s_daxpy(9,C1[0],dJvi_x,1,dJvi_dot_C1,1);
        aris::dynamic::s_daxpy(9,C1[1],dJvi_y,1,dJvi_dot_C1,1);
        aris::dynamic::s_daxpy(9,C1[2],dJvi_z,1,dJvi_dot_C1,1);

        std::fill_n(param_dds+3*legID,3,0);
        aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,df_sw_inB,1,1,param_dds+3*legID,1);

        std::fill_n(param_dsds1+3*legID,3,0);
        aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,ddf_sw_inB,1,1,param_dsds1+3*legID,1);
        std::fill_n(param_dsds2+3*legID,3,0);
        aris::dynamic::s_dgemm(3,1,3,1,dJvi_dot_df,3,df_sw_inB,1,1,param_dsds2+3*legID,1);

        std::fill_n(param_ds1+3*legID,3,0);
        aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,dgama_cross_df,1,1,param_ds1+3*legID,1);
        std::fill_n(param_ds2+3*legID,3,0);
        aris::dynamic::s_dgemm(3,1,3,1,dJvi_dot_df,3,C1,1,1,param_ds2+3*legID,1);

        std::fill_n(param_const1+3*legID,3,0);
        aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,C1,1,1,param_const1+3*legID,1);
        std::fill_n(param_const2+3*legID,3,0);
        aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,C2,1,1,param_const2+3*legID,1);//G*C2
        std::fill_n(param_const3+3*legID,3,0);
        aris::dynamic::s_dgemm(3,1,3,1,dJvi_dot_C1,3,C1,1,1,param_const3+3*legID,1);//H*C1

        std::fill_n(param_dsds+3*legID,3,0);
        std::fill_n(param_ds+3*legID,3,0);
        for (int k=0;k<3;k++)
        {
            param_dsds[3*legID+k]=param_dsds1[3*legID+k]+param_dsds2[3*legID+k];
            param_ds[3*legID+k]=2*(-param_ds1[3*legID+k]+param_ds2[3*legID+k])*sb_mtx[1];
            param_const[3*legID+k]=param_const1[3*legID+k]*sb_mtx[2]+(param_const2[3*legID+k]+param_const3[3*legID+k])*sb_mtx[1]*sb_mtx[1];

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
                    int swCount = count<count3 ? (count-count1) : (count-count3);
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
//            memcpy(*output_ds1+18*count+3*legID, param_ds1+3*legID, 3*sizeof(double));
//            memcpy(*output_ds2+18*count+3*legID, param_ds2+3*legID, 3*sizeof(double));
            memcpy(*output_const+18*count+3*legID,  param_const+3*legID,  3*sizeof(double));
//            memcpy(*output_const1+18*count+3*legID,  param_const1+3*legID,  3*sizeof(double));
//            memcpy(*output_const2+18*count+3*legID,  param_const2+3*legID,  3*sizeof(double));
//            memcpy(*output_const3+18*count+3*legID,  param_const3+3*legID,  3*sizeof(double));
            memcpy(*output_dds+18*count+3*legID, param_dds+3*legID, 3*sizeof(double));
            memcpy(*output_a2+18*count+3*legID,  param_a2+3*legID,  3*sizeof(double));
            memcpy(*output_a1+18*count+3*legID,  param_a1+3*legID,  3*sizeof(double));
            memcpy(*output_a0L+18*count+3*legID, param_a0L+3*legID, 3*sizeof(double));
            memcpy(*output_a0H+18*count+3*legID, param_a0H+3*legID, 3*sizeof(double));
        }
    }

	/*获取摆动腿加速度下限的最大值*/
    double NonRTOptimalGCSRotate::GetSwingMaxDec(int count, double ds, int legID)
    {
        int swCount = count<count3 ? (count-count1) : (count-count3);
        double dec[3] {0};
        std::fill_n(dec,3,-1e6);
        for (int k=0;k<3;k++)
        {
            if(isParamddsExact0[swCount][3*legID+k]==-1)
            {
                dec[k]=output_a2[count][3*legID+k]*ds*ds+output_a1[count][3*legID+k]*ds+output_a0L[count][3*legID+k];
            }
        }
        return *std::max_element(dec,dec+3);
    }

	/*获取摆动腿加速度上限的最小值*/
    double NonRTOptimalGCSRotate::GetSwingMinAcc(int count, double ds, int legID)
    {
        int swCount = count<count3 ? (count-count1) : (count-count3);
        double acc[3] {0};
        std::fill_n(acc,3,1e6);
        for (int k=0;k<3;k++)
        {
            if(isParamddsExact0[swCount][3*legID+k]==-1)
            {
                acc[k]=output_a2[count][3*legID+k]*ds*ds+output_a1[count][3*legID+k]*ds+output_a0H[count][3*legID+k];
            }
        }
        return *std::min_element(acc,acc+3);
    }

	/*遍历法，获取摆动腿速度上限值*/
    void NonRTOptimalGCSRotate::GetSwingDsBound(int count, int legID)
    {
        bool ds_lowBoundFlag_sw {false};
        bool ds_upBoundFlag_sw {false};
        int k_sw {0};
        int swCount = count<count3 ? (count-count1) : (count-count3);
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
        double Jvi_dot_C1[3] {0};
        double pEE[18];
        rbt.pLegs[legID]->GetPee(pEE+3*legID);

        double C1_tmp[3] {0};//C1 in Thesis
        C1_tmp[0]=-db_sb[count][2]-db_sb[count][0]*(pEE[3*legID+2]-b_sb[count][1]);
        C1_tmp[2]=-db_sb[count][1]+db_sb[count][0]*(pEE[3*legID]-b_sb[count][2]);
        double C1[3] {0};
        C1[0]=C1_tmp[0]*cos(b_sb[count][0])-C1_tmp[2]*sin(b_sb[count][2]);
        C1[2]=C1_tmp[0]*sin(b_sb[count][0])+C1_tmp[2]*cos(b_sb[count][2]);

        aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,C1,1,1,Jvi_dot_C1,1);
        for (int k=0;k<3;k++)
        {
            if(output_dds[count][3*legID+k]>0)
            {
                vLmt_value[k]=(Jvi_dot_C1[k]+vLmt)/output_dds[count][3*legID+k];
            }
            else if(output_dds[count][3*legID+k]<0)
            {
                vLmt_value[k]=(Jvi_dot_C1[k]-vLmt)/output_dds[count][3*legID+k];
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

	/*牛顿迭代法，获取摆动腿速度上限值*/
    void NonRTOptimalGCSRotate::GetSwingDsBoundByNewton(int count, int legID)
    {
        int swCount = count<count3 ? (count-count1) : (count-count3);

        //lowBound
        double ds0 {0};
        double value0;
        double ds1;
        double value_up;
        double value_low;
        double slope;

        int k;
        value0=GetSwingMinAcc(count,ds0,legID)-GetSwingMaxDec(count,ds0,legID);
        if(value0>=0)
        {
            ds_lowBound[swCount][legID]=ds_lowBound_aLmt[swCount][legID]=0;
        }
        else
        {
            k=0;
            value_up=-1;
            ds1=1e-3;
            while(value_up<=0 && k<50)
            {
                k++;
                ds0=ds1;
                value_low=GetSwingMinAcc(count,ds0-1e-3,legID)-GetSwingMaxDec(count,ds0-1e-3,legID);
                value_up=GetSwingMinAcc(count,ds0+1e-3,legID)-GetSwingMaxDec(count,ds0+1e-3,legID);
    //            value0=GetSwingMinAcc(count,ds0,legID)-GetSwingMaxDec(count,ds0,legID);
    //            slope=(value_up-value_low)/2e-3;
    //            ds1=ds0-value0/slope;
                ds1=(value_up*(ds0-1e-3)-value_low*(ds0+1e-3))/(value_up-value_low);
            }
            //printf("First Newton Iteration of lowBound for SwingLeg stops at k=%d, value_low=%f, value_up=%f, ds0=%f, ds1=%f\n",k,value_low,value_up,ds0,ds1);

            ds1=ds0;
            k=0;
            while(value_up<=0 && k<50)
            {
                k++;
                ds0=ds1;
                value_low=GetSwingMinAcc(count,ds0-1e-6,legID)-GetSwingMaxDec(count,ds0-1e-6,legID);
                value_up=GetSwingMinAcc(count,ds0+1e-6,legID)-GetSwingMaxDec(count,ds0+1e-6,legID);
    //            value0=GetSwingMinAcc(count,ds0,legID)-GetSwingMaxDec(count,ds0,legID);
    //            slope=(value_up-value_low)/2e-6;
    //            ds1=ds0-value0/slope;
                ds1=(value_up*(ds0-1e-6)-value_low*(ds0+1e-6))/(value_up-value_low);
            }
            //printf("Second Newton Iteration of lowBound for SwingLeg stops at k=%d, value_low=%f,value_up=%f, ds0=%f, ds1=%f\n",k,value_low,value_up,ds0,ds1);
            ds_lowBound[swCount][legID]=ds_lowBound_aLmt[swCount][legID]=ds1;
        }

        //upBound
        ds0=ds_lowBound[swCount][legID]+1e-3;//when ds_lowBound=0, ds0=ds1=0 all the time. use +1e-3 to avoid this.
        value0=GetSwingMinAcc(count,ds0,legID)-GetSwingMaxDec(count,ds0,legID);
        k=0;
        while (value0>=0 && k<10)
        {
            k++;
            ds0*=10;
            value0=GetSwingMinAcc(count,ds0,legID)-GetSwingMaxDec(count,ds0,legID);
        }
        if(k==0)
        {
            printf("Error, ds_lowBound=%f not right, please check!",ds0);
            std::abort();
        }
        else
        {
            ds1=ds0;
        }
        //printf("ds0=%f,ds1=%f,value0=%f\n",ds0,ds1,value0);

        k=0;
        value_low=-1;
        while(value_low<=0 && k<50)
        {
            k++;
            ds0=ds1;
            value_low=GetSwingMinAcc(count,ds0-1e-3,legID)-GetSwingMaxDec(count,ds0-1e-3,legID);
            value_up=GetSwingMinAcc(count,ds0+1e-3,legID)-GetSwingMaxDec(count,ds0+1e-3,legID);
    //        value0=GetSwingMinAcc(count,ds0,legID)-GetSwingMaxDec(count,ds0,legID);
    //        slope=(value_up-value_low)/2e-3;
    //        ds1=ds0-value0/slope;
            ds1=(value_up*(ds0-1e-3)-value_low*(ds0+1e-3))/(value_up-value_low);
        }
        //printf("First Newton iteration of upBound for SwingLeg stops at k=%d, value_low=%f, value_up=%f, ds0=%f, ds1=%f\n",k,value_low,value_up,ds0,ds1);

        ds1=ds0;
        k=0;
        value_low=-1;
        while(value_low<=0 && k<50)
        {
            k++;
            ds0=ds1;
            value_low=GetSwingMinAcc(count,ds0-1e-6,legID)-GetSwingMaxDec(count,ds0-1e-6,legID);
            value_up=GetSwingMinAcc(count,ds0+1e-6,legID)-GetSwingMaxDec(count,ds0+1e-6,legID);
    //        value0=GetSwingMinAcc(count,ds0,legID)-GetSwingMaxDec(count,ds0,legID);
    //        slope=(value_up-value_low)/2e-6;
    //        ds1=ds0-value0/slope;
            ds1=(value_up*(ds0-1e-6)-value_low*(ds0+1e-6))/(value_up-value_low);
        }
        //printf("Second Newton iteration of upBound for SwingLeg stops at k=%d, value_low=%f,value_up=%f, ds0=%f, ds1=%f\n",k,value_low,value_up,ds0,ds1);
        ds_upBound_aLmt[swCount][legID]=ds1;

        double vLmt_value[3] {0};
        double Jvi_dot_C1[3] {0};
        double pEE[18];
        rbt.pLegs[legID]->GetPee(pEE+3*legID);

        double C1_tmp[3] {0};//C1 in Thesis
        C1_tmp[0]=-db_sb[count][2]-db_sb[count][0]*(pEE[3*legID+2]-b_sb[count][1]);
        C1_tmp[2]=-db_sb[count][1]+db_sb[count][0]*(pEE[3*legID]-b_sb[count][2]);
        double C1[3] {0};
        C1[0]=C1_tmp[0]*cos(b_sb[count][0])-C1_tmp[2]*sin(b_sb[count][2]);
        C1[2]=C1_tmp[0]*sin(b_sb[count][0])+C1_tmp[2]*cos(b_sb[count][2]);
        aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,C1,1,1,Jvi_dot_C1,1);
        for (int k=0;k<3;k++)
        {
            if(output_dds[count][3*legID+k]>0)
            {
                vLmt_value[k]=(Jvi_dot_C1[k]+vLmt)/output_dds[count][3*legID+k];
            }
            else if(output_dds[count][3*legID+k]<0)
            {
                vLmt_value[k]=(Jvi_dot_C1[k]-vLmt)/output_dds[count][3*legID+k];
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

	/*显式求解二次函数曲线交点，取最小值，某些曲线可能无交点*/
    void NonRTOptimalGCSRotate::GetSwingDsBoundByFunc(int count, int legID)
    {
        int swCount = count<count3 ? (count-count1) : (count-count3);

        double a2[3] {0};
        double a1[3] {0};
        double a0H[3] {0};
        double a0L[3] {0};
        int num {0};
        double FuncA2;
        double FuncA1;
        double FuncA0;
        double ds;
        double minDs {1e6};

        for (int k=0;k<3;k++)
        {
            if(isParamddsExact0[swCount][3*legID+k]==-1)
            {
                a2[num]=output_a2[count][3*legID+k];
                a1[num]=output_a1[count][3*legID+k];
                a0H[num]=output_a0H[count][3*legID+k];
                a0L[num]=output_a0L[count][3*legID+k];
                num++;
            }
        }

        for(int i=0;i<num;i++)
        {
            for(int j=0;j<num;j++)
            {
                if(i!=j)
                {
                    FuncA2=a2[i]-a2[j];
                    FuncA1=a1[i]-a1[j];
                    FuncA0=a0H[i]-a0L[j];
                    if(FuncA2==0 || FuncA1*FuncA1-4*FuncA0*FuncA2<0)
                    {
                        ds=-1;
                    }
                    else
                    {
                        double ds1=(sqrt(FuncA1*FuncA1-4*FuncA0*FuncA2)-FuncA1)/(2*FuncA2);
                        double ds2=(-sqrt(FuncA1*FuncA1-4*FuncA0*FuncA2)-FuncA1)/(2*FuncA2);
                        ds=std::max(ds1,ds2);
                    }
                    if(minDs>ds && ds>0 && GetSwingMaxDec(count,ds-1e-3,legID)<GetSwingMinAcc(count,ds-1e-3,legID))
                    {
                        minDs=ds;
                    }
                }
            }
        }

        if(minDs==1e6)
        {
            printf("WARNING!!! Swing minDs too small or minDs donot exist at count %d!\n",count);
        }

        ds_upBound_aLmt[swCount][legID]=minDs;
        //printf("ds_upBound_aLmt finished!\n");

        double vLmt_value[3] {0};
        double Jvi_dot_C1[3] {0};
        double pEE[18];
        rbt.pLegs[legID]->GetPee(pEE+3*legID);

        double C1_tmp[3] {0};//C1 in Thesis
        C1_tmp[0]=-db_sb[count][2]-db_sb[count][0]*(pEE[3*legID+2]-b_sb[count][1]);
        C1_tmp[2]=-db_sb[count][1]+db_sb[count][0]*(pEE[3*legID]-b_sb[count][2]);
        double C1[3] {0};
        C1[0]=C1_tmp[0]*cos(b_sb[count][0])-C1_tmp[2]*sin(b_sb[count][2]);
        C1[2]=C1_tmp[0]*sin(b_sb[count][0])+C1_tmp[2]*cos(b_sb[count][2]);

        aris::dynamic::s_dgemm(3,1,3,1,Jvi,3,C1,1,1,Jvi_dot_C1,1);
        for (int k=0;k<3;k++)
        {
            if(output_dds[count][3*legID+k]>0)
            {
                vLmt_value[k]=(Jvi_dot_C1[k]+vLmt)/output_dds[count][3*legID+k];
            }
            else if(output_dds[count][3*legID+k]<0)
            {
                vLmt_value[k]=(Jvi_dot_C1[k]-vLmt)/output_dds[count][3*legID+k];
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

	/*获取摆动腿三类切换点*/
    void NonRTOptimalGCSRotate::GetSwingSwitchPoint(int legID)
    {
        double slopedsBound_for[901] {0};
        double slopedsBound_back[901] {0};
        double paramdds0Point[901] {0};
        double tangentPoint[901] {0};
        int paramdds0Count {0};
        int tangentCount {0};
        int switchScrewID_tmp[901] {0};

        //initialize
        paramdds0Count=0;
        tangentCount=0;
        for(int i=0;i<swingSCount+1;i++)
        {
            tangentPoint[i]=-1;
            switchScrewID_tmp[i]=-1;
            paramdds0Point[i]=-1;
        }

        slopedsBound_for[0]=slopedsBound_back[0]=(ds_upBound[1][legID]*ds_upBound[1][legID]-ds_upBound[0][legID]*ds_upBound[0][legID])/2/(s_w[1][legID]-s_w[0][legID]);
        for(int i=1;i<swingSCount;i++)
        {
            slopedsBound_for[i]=(ds_upBound[i+1][legID]*ds_upBound[i+1][legID]-ds_upBound[i][legID]*ds_upBound[i][legID])/2/(s_w[i+1][legID]-s_w[i][legID]);
            slopedsBound_back[i]=(ds_upBound[i][legID]*ds_upBound[i][legID]-ds_upBound[i-1][legID]*ds_upBound[i-1][legID])/2/(s_w[i][legID]-s_w[i-1][legID]);
        }
        slopedsBound_for[swingSCount]=slopedsBound_back[swingSCount]=(ds_upBound[swingSCount][legID]*ds_upBound[swingSCount][legID]-ds_upBound[swingSCount-1][legID]*ds_upBound[swingSCount-1][legID])/2/(s_w[swingSCount][legID]-s_w[swingSCount-1][legID]);

    //    for(int i=0;i<swingCount+1;i++)
    //    {
    //        slopeDelta[i][legID]=slopedsBound_for[i]-dds_upBound[i][legID];
    //    }

        for(int i=0;i<swingSCount;i++)
        {
            bool isParamdds0 {false};
            if(ds_upBound_vLmt[i][legID]>ds_upBound_aLmt[i][legID])
            {
                int count=(legID%2==0) ? (i+count1) : (i+count3);
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

            //if((slopeDelta[i+1][legID]*slopeDelta[i][legID]<0 || slopeDelta[i][legID]==0) && isParamdds0==false)
            if(slopedsBound_for[i]<=dds_upBound[i][legID] && slopedsBound_for[i+1]>dds_upBound[i+1][legID] && isParamdds0==false)
            {
                tangentPoint[tangentCount]=i
                        +fabs(slopedsBound_for[i]-dds_upBound[i][legID])/(fabs(slopedsBound_for[i]-dds_upBound[i][legID])+fabs(slopedsBound_for[i+1]-dds_upBound[i+1][legID]));
                tangentCount++;
                if(fabs(slopedsBound_for[i]-dds_upBound[i][legID])+fabs(slopedsBound_for[i+1]-dds_upBound[i+1][legID])==0)
                {
                    printf("WARNING! tangentPoint!!!\n");
                }
            }
            if(slopedsBound_back[i]<=dds_lowBound[i][legID] && slopedsBound_back[i+1]>dds_lowBound[i+1][legID]  && isParamdds0==false)
            {
                tangentPoint[tangentCount]=i
                        +fabs(slopedsBound_back[i]-dds_lowBound[i][legID])/(fabs(slopedsBound_back[i]-dds_lowBound[i][legID])+fabs(slopedsBound_back[i+1]-dds_lowBound[i+1][legID]));
                tangentCount++;
            }
        }

//        printf("SwingLeg %d Tangent Switch Point:",legID);
//        for(int i=0;i<tangentCount+1;i++)
//        {
//            printf("%.1f,",tangentPoint[i]);
//        }
//        printf("\n");

//        printf("SwingLeg %d ZeroInertia Switch Point:",legID);
//        for(int i=0;i<paramdds0Count+1;i++)
//        {
//            printf("%.1f,",paramdds0Point[i]);
//        }
//        printf("\n");

        //merge tangentPoint & paramdds0Point into switchPoint
        switchPoint[0][legID]=0;
        switchPoint[1][legID]=swingSCount;
        switchCount[legID]=2;
        for(int i=0;i<tangentCount;i++)
        {
            switchPoint[i+switchCount[legID]][legID]=tangentPoint[i];
            switchType[i+switchCount[legID]][legID]='t';
        }
        switchCount[legID]+=tangentCount;
        for(int i=0;i<paramdds0Count;i++)
        {
            switchPoint[i+switchCount[legID]][legID]=paramdds0Point[i];
            switchType[i+switchCount[legID]][legID]='z';
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
                    switchScrewID[i][legID]=switchScrewID[k][legID];
                    switchScrewID[k][legID]=tmp2;

                    auto tmp3=switchType[i][legID];
                    switchType[i][legID]=switchType[k][legID];
                    switchType[k][legID]=tmp3;
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
            printf("%.4f,",s_w[(int)switchPoint[i][legID]][legID]+(switchPoint[i][legID]-(int)switchPoint[i][legID])*(s_w[(int)switchPoint[i][legID]+1][legID]-s_w[(int)switchPoint[i][legID]][legID]));
        }
        printf("\n");
        printf("SwingLeg %d Switch Type:",legID);
        for(int i=0;i<switchCount[legID]+1;i++)
        {
            printf("%c,",switchType[i][legID]);
        }
        printf("\n");
    }

	/*获取摆动腿前向后向积分的起点*/
    void NonRTOptimalGCSRotate::GetSwingTwoPointAtSwitch(int legID, double *lowPoint, double *upPoint)
    {
        lowPoint[0]=-1;
        //lowPoint[switchCount[legID]-1]=0;
        lowPoint[switchCount[legID]-1]=ds_upBound[count2-count1][legID];
        //upPoint[0]=0;
        upPoint[0]=ds_upBound[0][legID];
        upPoint[switchCount[legID]-1]=-1;

        for(int i=1;i<switchCount[legID]-1;i++)
        {

            if(switchType[i][legID]=='t')
            {
                int num=(int)switchPoint[i][legID];
                lowPoint[i]=ds_upBound[num][legID];
                upPoint[i]=ds_upBound[num+1][legID];
            }
            else if(switchType[i][legID]=='z')
            {

                int num=round(switchPoint[i][legID]);
                auto tmp=std::min(ds_upBound[num-1][legID],ds_upBound[num][legID]);
                upPoint[i]=lowPoint[i]=std::min(tmp,ds_upBound[num+1][legID]);//std::min(ds_upBound[num][legID],ds_upBound[num+1][legID]);
            }
        }
    }

	/*切换点法，获取最优速度曲线*/
    void NonRTOptimalGCSRotate::GetSwingOptimalDsBySwitchPoint(int legID)
    {
        bool stopFlag {false};
        int forwardEnd_s {-1};
        double forwardEnd_ds {0};
        double *lowPoint=new double [switchCount[legID]];
        double *upPoint=new double [switchCount[legID]];
        GetSwingTwoPointAtSwitch(legID,lowPoint,upPoint);

//        printf("lowPoint:");
//        for(int i=0;i<switchCount[legID];i++)
//        {
//            printf("%.2f,",lowPoint[i]);
//        }
//        printf("\n");
//        printf("upPoint:");
//        for(int i=0;i<switchCount[legID];i++)
//        {
//            printf("%.2f,",upPoint[i]);
//        }
//        printf("\n");


        for(int i=0;i<swingSCount+1;i++)
        {
            real_ds[i][legID]=ds_upBound[i][legID];
            real_dds[i][legID]=dds_upBound[i][legID];
        }
        for(int m=0;m<switchCount[legID];m++)
        {
            bool ignoreBackward {false};
            int k_sw=(int)switchPoint[m][legID];
            int k_sw_start=k_sw;
            if(switchType[m][legID]=='z')
            {
                k_sw_start=k_sw=round(switchPoint[m][legID])-1;
            }

            if(k_sw_start==0)
            {
                ignoreBackward=true;
            }
            else if(k_sw_start<forwardEnd_s)
            {
                ignoreBackward=true;
                //printf("SwingLeg backward start at a passed point, quit switchPoint %.1f\n",switchPoint[m][legID]);
            }
            else if(k_sw_start==forwardEnd_s && lowPoint[m]>forwardEnd_ds)
            {
                ignoreBackward=true;
                if(switchType[m][legID]=='z')
                {
                    real_dds[k_sw_start+1][legID]=0;
                    real_ds[k_sw_start+1][legID]=lowPoint[m];
                }
            }

            if(ignoreBackward==false)
            {
                if(switchType[m][legID]=='z')
                {
                    real_dds[k_sw_start+1][legID]=0;
                    real_ds[k_sw_start+1][legID]=lowPoint[m];
                }
                stopFlag=false;
                ds_backward[k_sw][legID]=lowPoint[m];
                while(stopFlag==false)
                {
                    int count=(legID%2==0) ? (k_sw+count1) : (k_sw+count3);
                    dds_backward[k_sw][legID]=GetSwingMaxDec(count,ds_backward[k_sw][legID],legID);
                    ds_backward[k_sw-1][legID]=sqrt(ds_backward[k_sw][legID]*ds_backward[k_sw][legID]-2*dds_backward[k_sw][legID]*(s_w[k_sw][legID]-s_w[k_sw-1][legID]));

                    if(ds_backward[k_sw-1][legID]>ds_upBound[k_sw-1][legID])
                    {
                        real_dds[k_sw-1][legID]=(real_ds[k_sw-1][legID]-ds_backward[k_sw][legID])*(real_ds[k_sw-1][legID]+ds_backward[k_sw][legID])/2/(s_w[k_sw][legID]-s_w[k_sw-1][legID]);
                        for(int i=k_sw;i<k_sw_start+1;i++)
                        {
                            real_ds[i][legID]=ds_backward[i][legID];
                            real_dds[i][legID]=dds_backward[i][legID];
                        }
                        stopFlag=true;
                        printf("SwingLeg backward touching upBound at %d, from switchPoint %.4f\n",k_sw-1,switchPoint[m][legID]);
                    }
                    else if(k_sw==1)
                    {
                        dds_backward[k_sw-1][legID]=GetSwingMaxDec(count-1,ds_backward[k_sw-1][legID],legID);
                        for(int i=k_sw-1;i<k_sw_start+1;i++)
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
                        for(int i=k_sw;i<k_sw_start+1;i++)
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
    //            //if(ds_backward[k_sw-1][legID]>ds_upBound[k_sw-1][legID] || round(switchPoint[m][legID])==900)
    //            if(k_sw_start==900)
    //            {
    //                continue;
    //            }
            }

            bool ignoreForward {false};
            k_sw_start=k_sw=(int)switchPoint[m][legID];
            if(switchType[m][legID]=='t')
            {
                k_sw_start=k_sw=(int)switchPoint[m][legID]+1;
            }
            else if(switchType[m][legID]=='z')
            {
                k_sw_start=k_sw=round(switchPoint[m][legID])+1;
            }

            if(k_sw_start==swingSCount+1 || k_sw_start==swingSCount)
            {
                ignoreForward=true;
            }
            else if(k_sw_start<forwardEnd_s)
            {
                ignoreForward=true;
                //printf("SwingLeg forward start at a passed point, quit switchPoint %.1f\n",switchPoint[m][legID]);
            }
            else if(k_sw_start==forwardEnd_s && upPoint[m]>forwardEnd_ds)
            {
                printf("How possible! SwingLeg forward curve should not stop here!\n");
            }

            if(ignoreForward==false)
            {
                stopFlag=false;
                ds_forward[k_sw][legID]=upPoint[m];
                while(stopFlag==false)
                {
                    int count=(legID%2==0) ? (k_sw+count1) : (k_sw+count3);
                    dds_forward[k_sw][legID]=GetSwingMinAcc(count,ds_forward[k_sw][legID],legID);
                    ds_forward[k_sw+1][legID]=sqrt(ds_forward[k_sw][legID]*ds_forward[k_sw][legID]+2*dds_forward[k_sw][legID]*(s_w[k_sw+1][legID]-s_w[k_sw][legID]));

                    if(ds_forward[k_sw+1][legID]>ds_upBound[k_sw+1][legID] || k_sw==swingSCount-1)
                    {
                        for(int i=k_sw_start;i<k_sw+1;i++)
                        {
                            real_ds[i][legID]=ds_forward[i][legID];
                            real_dds[i][legID]=dds_forward[i][legID];
                        }
        //                if(round(switchPoint[m][legID])==0)
        //                {
        //                    real_ds[0][legID]=ds_forward[0][legID];
        //                    real_dds[0][legID]=dds_forward[0][legID];
        //                }
                        if(k_sw==swingSCount-1)
                        {
                            if(ds_forward[k_sw+1][legID]>ds_upBound[k_sw+1][legID])
                            {
                                real_ds[swingSCount][legID]=ds_upBound[swingSCount][legID];
                                real_dds[swingSCount][legID]=dds_upBound[swingSCount][legID];
                            }
                            else
                            {
                                real_ds[swingSCount][legID]=ds_forward[swingSCount][legID];
                                real_dds[swingSCount][legID]=GetSwingMinAcc(count+1,ds_forward[swingSCount][legID],legID);
                            }
                        }
                        stopFlag=true;
                        forwardEnd_s=k_sw;
                        forwardEnd_ds=ds_forward[k_sw][legID];
                        printf("SwingLeg forward touching upBound at %d, from switchPoint %.4f\n",k_sw,switchPoint[m][legID]);
                    }
                    else
                    {
                        k_sw++;
                    }
                }
            }
        }

        for(int i=0;i<swingSCount+1;i++)
        {
            int count=(legID%2==0) ? (i+count1) : (i+count3);
            real_ddsMax[i][legID]=GetSwingMinAcc(count,real_ds[i][legID],legID);
            real_ddsMin[i][legID]=GetSwingMaxDec(count,real_ds[i][legID],legID);
        }

        delete [] lowPoint;
        delete [] upPoint;
    }

	/*摆动腿边界条件要求速度为零，因此需要从零分别再进行积分*/
    void NonRTOptimalGCSRotate::ApplySwingExtraIntegration(int legID)
    {
        bool stopFlag {false};
        int k_sw {0};
        int k_sw_start {0};
        double real_ds_tmp {0};

        //forward
        real_ds[0][legID]=0;
        k_sw_start=k_sw=0;
        stopFlag=false;
        while(stopFlag==false && k_sw<count2-count1)
        {
            int count=(legID%2==0) ? (k_sw+count1) : (k_sw+count3);
            real_dds[k_sw][legID]=GetSwingMinAcc(count,real_ds[k_sw][legID],legID);
            real_ds_tmp=sqrt(real_ds[k_sw][legID]*real_ds[k_sw][legID]+2*real_dds[k_sw][legID]*(s_w[k_sw+1][legID]-s_w[k_sw][legID]));
            if(real_ds_tmp>real_ds[k_sw+1][legID])
            {
                stopFlag=true;
            }
            else
            {
                real_ds[k_sw+1][legID]=real_ds_tmp;
                k_sw++;
            }
        }
        printf("ExtraIntegration(forward) stops at k_sw=%d\n",k_sw);

        //backward
        real_ds[count2-count1][legID]=0;
        k_sw_start=k_sw=count2-count1;
        stopFlag=false;
        while(stopFlag==false  && k_sw>0)
        {
            int count=(legID%2==0) ? (k_sw+count1) : (k_sw+count3);
            real_dds[k_sw][legID]=GetSwingMaxDec(count,real_ds[k_sw][legID],legID);
            real_ds_tmp=sqrt(real_ds[k_sw][legID]*real_ds[k_sw][legID]-2*real_dds[k_sw][legID]*(s_w[k_sw][legID]-s_w[k_sw-1][legID]));
            if(real_ds_tmp>real_ds[k_sw-1][legID])
            {
                stopFlag=true;
            }
            else
            {
                real_ds[k_sw-1][legID]=real_ds_tmp;
                k_sw--;
            }
        }
        printf("ExtraIntegration(backward) stops at k_sw=%d\n",k_sw);
    }

	/*不使用切换点，从起点和终点开始尝试积分，超出ds上限后退回再尝试，如此实现求解。该方法开发于不了解切换点时，计算耗时大于切换点且计算误差大，不建议使用*/
    void NonRTOptimalGCSRotate::GetSwingOptimalDsByDirectNI(int legID)
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
        ki_back=swingSCount;
        ds_backward[ki_back][legID]=0;
        while (stopFlag==false && ki_back>=0)
        {
            int count=(legID%2==0) ? (ki_back+count1) : (ki_back+count3);
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
        std::fill_n(min_dist,swingSCount+1,1);
        while (stopFlag==false)
        {
            if (accFlag==true)
            {
                int count=(legID%2==0) ? (ki_for+count1) : (ki_for+count3);
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
                int count=(legID%2==0) ? (ki_for+count1) : (ki_for+count3);
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
                            <=(ds_lowBound[ki_for+1][legID]*ds_lowBound[ki_for+1][legID]) || ki_for==swingSCount-1)
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

            if(ki_for==swingSCount)
            {
                stopFlag=true;
                for(int i=0;i<swingSCount+1;i++)
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
                for(int i=ki_for;i<swingSCount+1;i++)
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

	/*身体和摆动腿最优速度耗时不同，需拉伸同步，拉伸后s_b改变，因此需要迭代求解*/
    void NonRTOptimalGCSRotate::GetOptimalDsByMajorIteration()
    {
        rbt.loadXml("../../../RobotVIII_demo/resource/RobotEDU2.xml");

        timeval tpstart,tpend;
        float tused;
        gettimeofday(&tpstart,NULL);

        //AnalyseStanceOptimalTime();
        //GetStanceOptimalDsByMinorIteration();
        GetStanceOptimalDsByNewton();

        //pva_b & s_w initialized here, and need to be updated during major iteration
        for (int i=0;i<count5+1;i++)
        {
            s_b_mtrx_scale[i][0]=s_b[i];
            s_b_mtrx_scale[i][1]=real_ds_body[i];
            s_b_mtrx_scale[i][2]=real_dds_body[i];
            for(int j=0;j<3;j++)
            {
                output_init_pva_b[i][j]=pva_b[i][j]=b_sb[i][j];
                output_init_pva_b[i][j+3]=pva_b[i][j+3]=db_sb[i][j]*real_ds_body[i];
                output_init_pva_b[i][j+6]=pva_b[i][j+6]=ddb_sb[i][j]*real_ds_body[i]*real_ds_body[i]+db_sb[i][j]*real_dds_body[i];
            }

            for(int j=0;j<3;j++)
            {
                if(i>=count1 && i<count2+1)
                {
                    s_w[i-count1][2*j]=PI/(s_b2-s_b1)*(s_b[i]-s_b1);
                }
                if(i>=count3 && i<count4+1)
                {
                    s_w[i-count3][2*j+1]=PI/(s_b4-s_b3)*(s_b[i]-s_b3);
                }
            }
        }

        /******************************* SwingLeg & Iteration ************************************/
        int iterCount {0};
        bool stopFlag {false};
        double Tstep_tmp {Tstep};
        double totalTime[6] {0};
        double maxTime_last {0};
        //double totalTime_last[6] {0};
        double real_ds_scale_tmp[901][6] {0};
        while(stopFlag==false && iterCount<=30)
        {
            timeval tpstart1,tpend1;
            float tused1;
            gettimeofday(&tpstart1,NULL);

            printf("\n");
            iterCount++;

            //generate the traj & calculate the bound of ds
            std::fill_n(*ds_backward,(swingSCount+1)*6,0);
            std::fill_n(*ds_forward,(swingSCount+1)*6,0);
            std::fill_n(totalTime,6,0);
            std::fill_n(switchCount,6,0);
            std::fill_n(*switchPoint,(swingSCount+1)*6,-1);
            std::fill_n(*switchScrewID,(swingSCount+1)*6,-1);
            std::fill_n(*switchType,(swingSCount+1)*6,'0');
            std::fill_n(*isParamddsExact0,(swingSCount+1)*18,-1);

            for (int j=0;j<6;j++)
            {
                for (int i=0;i<swingSCount+1;i++)//024
                {
                    int count = (j%2==0) ? (i+count1) : (i+count3);
                    GetSwingLegParam(count,j,s_w[i][j],*s_b_mtrx_scale+3*count);
                    //GetSwingDsBound(count,j);
                    //GetSwingDsBoundByNewton(count,j);
                    GetSwingDsBoundByFunc(count,j);
                }

                GetSwingSwitchPoint(j);
                GetSwingOptimalDsBySwitchPoint(j);
                ApplySwingExtraIntegration(j);
                //GetSwingOptimalDsByDirectNI(j);

                for (int i=1;i<swingSCount+1;i++)
                {
                    totalTime[j]+=2*(s_w[i][j]-s_w[i-1][j])/(real_ds[i-1][j]+real_ds[i][j]);
                    timeArray[i][j]=totalTime[j];
                }

            }
            maxSwingTime=*std::max_element(totalTime,totalTime+6);
            printf("totalTime: %.4f, %.4f, %.4f, %.4f, %.4f, %.4f; maxTime:%.4f\n",
                   totalTime[0],totalTime[1],totalTime[2],totalTime[3],totalTime[4],totalTime[5],maxSwingTime);

            stopFlag=true;
    //        for(int i=0;i<6;i++)
    //        {
    //            //if((int)(totalTime_last[i]*1000)!=(int)(totalTime[i]*1000))
    //            if(fabs(totalTime_last[i]-totalTime[i])>1e-4)
    //            {
    //                stopFlag=false;
    //                //printf("totalTime_last=%.5f,totalTime=%.5f\n",totalTime_last[i],totalTime[i]);
    //            }
    //        }
            if(fabs(maxTime_last-maxSwingTime)>1e-4)
            {
                stopFlag=false;
                //printf("totalTime_last=%.5f,totalTime=%.5f\n",totalTime_last[i],totalTime[i]);
            }
            if(stopFlag==true)
            {
                stepTimeCount=(int)(maxSwingTime/(1-dutyCycle)*1000)+1;
                if(stepTimeCount%2==1)
                {
                    stepTimeCount++;
                }
                Tstep_tmp=0.001*stepTimeCount;
                swingTimeCount=(int)(maxSwingTime*1000)+1;
                maxSwingTime=0.001*swingTimeCount;
                dutyCycle=1-maxSwingTime/Tstep_tmp;
            }

            //update s_b & s_w here
            double s_b_mtrx_tmp[2201][3] {0};
            double pva_b_tmp[2201][9] {0};
            for(int i=0;i<count5+1;i++)
            {
                timeArray_body_tmp[i]=timeArray_body[i]*maxSwingTime/((1-dutyCycle)*Tstep);
                s_b_mtrx_tmp[i][0]=s_b[i];
                s_b_mtrx_tmp[i][1]=real_ds_body[i]*(1-dutyCycle)*Tstep/maxSwingTime;
                s_b_mtrx_tmp[i][2]=real_dds_body[i]*(1-dutyCycle)*Tstep/maxSwingTime*(1-dutyCycle)*Tstep/maxSwingTime;

                for(int j=0;j<3;j++)
                {
                    pva_b_tmp[i][j]  =output_init_pva_b[i][j];
                    pva_b_tmp[i][j+3]=output_init_pva_b[i][j+3]*(1-dutyCycle)*Tstep/maxSwingTime;
                    pva_b_tmp[i][j+6]=output_init_pva_b[i][j+6]*(1-dutyCycle)*Tstep/maxSwingTime*(1-dutyCycle)*Tstep/maxSwingTime;
                }
            }
            Tstep_tmp=timeArray_body_tmp[count5];

            for(int i=0;i<swingSCount+1;i++)
            {
                for(int j=0;j<6;j++)
                {
                    if(totalTime[j]!=0)
                    {
                        timeArray_tmp[i][j]=timeArray[i][j]*maxSwingTime/totalTime[j];
                    }
                }
            }

            int k_start {0};
            double s_w_tmp[901][6] {0};
            for(int j=0;j<3;j++)
            {
                k_start=0;
                for(int i=count1;i<count2+1;i++)
                {
                    for(int k=k_start;k<swingSCount;k++)
                    {
                        if(timeArray_tmp[k][2*j]<=(timeArray_body_tmp[i]-(0.5*dutyCycle-0.25)*Tstep_tmp)
                                && timeArray_tmp[k+1][2*j]>(timeArray_body_tmp[i]-(0.5*dutyCycle-0.25)*Tstep_tmp))
                        {
                            k_start=k;
                            s_w_tmp[i-count1][2*j]=s_w[k][2*j]+(s_w[k+1][2*j]-s_w[k][2*j])*(timeArray_body_tmp[i]-(0.5*dutyCycle-0.25)*Tstep_tmp-timeArray_tmp[k][2*j])/(timeArray_tmp[k+1][2*j]-timeArray_tmp[k][2*j]);
                            real_ds_scale[i-count1][2*j]=(real_ds[k][2*j]+(real_ds[k+1][2*j]-real_ds[k][2*j])*(timeArray_body_tmp[i]-(0.5*dutyCycle-0.25)*Tstep_tmp-timeArray_tmp[k][2*j])/(timeArray_tmp[k+1][2*j]-timeArray_tmp[k][2*j]))
                                    *totalTime[2*j]/maxSwingTime;
                            real_dds_scale[i-count1][2*j]=(real_dds[k][2*j]+(real_dds[k+1][2*j]-real_dds[k][2*j])*(timeArray_body_tmp[i]-(0.5*dutyCycle-0.25)*Tstep_tmp-timeArray_tmp[k][2*j])/(timeArray_tmp[k+1][2*j]-timeArray_tmp[k][2*j]))
                                    *totalTime[2*j]/maxSwingTime*totalTime[2*j]/maxSwingTime;
                            break;
                        }
                    }
                }
            }
            for(int j=0;j<3;j++)
            {
                k_start=0;
                for(int i=count3;i<count4+1;i++)
                {
                    for(int k=k_start;k<swingSCount;k++)
                    {
                        if(timeArray_tmp[k][2*j+1]<=(timeArray_body_tmp[i]-(0.5*dutyCycle+0.25)*Tstep_tmp)
                                && timeArray_tmp[k+1][2*j+1]>(timeArray_body_tmp[i]-(0.5*dutyCycle+0.25)*Tstep_tmp))
                        {
                            k_start=k;
                            s_w_tmp[i-count3][2*j+1]=s_w[k][2*j+1]+(s_w[k+1][2*j+1]-s_w[k][2*j+1])*(timeArray_body_tmp[i]-(0.5*dutyCycle+0.25)*Tstep_tmp-timeArray_tmp[k][2*j+1])/(timeArray_tmp[k+1][2*j+1]-timeArray_tmp[k][2*j+1]);
                            real_ds_scale[i-count3][2*j+1]=(real_ds[k][2*j+1]+(real_ds[k+1][2*j+1]-real_ds[k][2*j+1])*(timeArray_body_tmp[i]-(0.5*dutyCycle+0.25)*Tstep_tmp-timeArray_tmp[k][2*j+1])/(timeArray_tmp[k+1][2*j+1]-timeArray_tmp[k][2*j+1]))
                                    *totalTime[2*j+1]/maxSwingTime;
                            real_dds_scale[i-count3][2*j+1]=(real_dds[k][2*j+1]+(real_dds[k+1][2*j+1]-real_dds[k][2*j+1])*(timeArray_body_tmp[i]-(0.5*dutyCycle+0.25)*Tstep_tmp-timeArray_tmp[k][2*j+1])/(timeArray_tmp[k+1][2*j+1]-timeArray_tmp[k][2*j+1]))
                                    *totalTime[2*j+1]/maxSwingTime*totalTime[2*j+1]/maxSwingTime;
                            break;
                        }
                    }
                }
            }
            for(int i=0;i<6;i++)
            {
                s_w_tmp[swingSCount][i]=PI;
                real_ds_scale[swingSCount][i]=real_ds[swingSCount][i]*totalTime[i]/maxSwingTime;
                real_dds_scale[swingSCount][i]=real_dds[swingSCount][i]*totalTime[i]/maxSwingTime*totalTime[i]/maxSwingTime;
            }

            memcpy(*real_ds_scale_tmp,*real_ds_scale,(swingSCount+1)*6*sizeof(double));
            memcpy(*s_w,*s_w_tmp,(swingSCount+1)*6*sizeof(double));
            memcpy(*pva_b,*pva_b_tmp,(count5+1)*9*sizeof(double));
            memcpy(*s_b_mtrx_scale,*s_b_mtrx_tmp,(count5+1)*3*sizeof(double));
            maxTime_last=maxSwingTime;
            //memcpy(totalTime_last,totalTime,6*sizeof(double));

            gettimeofday(&tpend1,NULL);
            tused1=tpend1.tv_sec-tpstart1.tv_sec+(double)(tpend1.tv_usec-tpstart1.tv_usec)/1000000;
            printf("UsedTime leg:%f\n",tused1);
        }

        printf("Iteration finished, iterCount=%d, maxTime=%.4f, Tstep_tmp:%.4f, dutyCycle:%4f\n\n",iterCount,maxSwingTime,Tstep_tmp,dutyCycle);

        gettimeofday(&tpend,NULL);
        tused=tpend.tv_sec-tpstart.tv_sec+(double)(tpend.tv_usec-tpstart.tv_usec)/1000000;
        printf("UsedTime:%f\n",tused);
    }

	/*获取位置-s轨迹*/
    void NonRTOptimalGCSRotate::GetOptimalGait2s()
    {
        printf("start GetOptimalGait2s\n");
        double pEB[6] {0};
        double pEE[18] {0};
        double vEB[6] {0};
        double vEE[18] {0};
        double aEB[6] {0};
        double aEE[18] {0};
        double ratio_f[3] {0};

        std::fill_n(vEE,18,0);
        std::fill_n(aEE,18,0);
        for (int i=0;i<count1;i++)
        {
            pEB[0]=initPeb[0]+pva_b[i][2];
            pEB[2]=initPeb[2]+pva_b[i][1];
            pEB[3]=initPeb[3]+pva_b[i][0];
            vEB[0]=pva_b[i][5];
            vEB[2]=pva_b[i][4];
            vEB[3]=pva_b[i][3];
            aEB[0]=pva_b[i][8];
            aEB[2]=pva_b[i][7];
            aEB[3]=pva_b[i][6];
            for(int j=0;j<3;j++)
            {
                pEE[6*j]=initPee[6*j]*cos(-stepBeta/4)+initPee[6*j+2]*sin(-stepBeta/4);
                pEE[6*j+1]=initPee[6*j+1];
                pEE[6*j+2]=-initPee[6*j]*sin(-stepBeta/4)+initPee[6*j+2]*cos(-stepBeta/4);

                pEE[6*j+3]=initPee[6*j+3]*cos(stepBeta/4)+initPee[6*j+5]*sin(stepBeta/4);
                pEE[6*j+4]=initPee[6*j+4];
                pEE[6*j+5]=-initPee[6*j+3]*sin(stepBeta/4)+initPee[6*j+5]*cos(stepBeta/4);
            }
            rbt.SetPeb(pEB,"213");
            rbt.SetVb(vEB);
            rbt.SetAb(aEB);
            rbt.SetPee(pEE);
            rbt.SetVee(vEE);
            rbt.SetAee(aEE);

            rbt.GetPee(*Pee_s+18*i,rbt.body());
            rbt.GetPin(*Pin_s+18*i);
            rbt.GetVin(*Vin_s+18*i);
            rbt.GetAin(*Ain_s+18*i);
        }

        std::fill_n(vEE,18,0);
        std::fill_n(aEE,18,0);
        for (int i=count1;i<count2+1;i++)
        {
            pEB[0]=initPeb[0]+pva_b[i][2];
            pEB[2]=initPeb[2]+pva_b[i][1];
            pEB[3]=initPeb[3]+pva_b[i][0];
            vEB[0]=pva_b[i][5];
            vEB[2]=pva_b[i][4];
            vEB[3]=pva_b[i][3];
            aEB[0]=pva_b[i][8];
            aEB[2]=pva_b[i][7];
            aEB[3]=pva_b[i][6];
            for(int j=0;j<3;j++)
            {
                ratio_f[0]=initPee[6*j]*(cos(3*stepBeta/4)-cos(-stepBeta/4))+initPee[6*j+2]*(sin(3*stepBeta/4)-sin(-stepBeta/4));
                ratio_f[2]=-initPee[6*j]*(sin(3*stepBeta/4)-sin(-stepBeta/4))+initPee[6*j+2]*(cos(3*stepBeta/4)-cos(-stepBeta/4));

                f_sw[0]=-ratio_f[0]/2*(cos(s_w[i-count1][2*j])-1);
                f_sw[1]=stepH*sin(s_w[i-count1][2*j]);
                f_sw[2]=-ratio_f[2]/2*(cos(s_w[i-count1][2*j])-1);

                df_sw[0]=ratio_f[0]/2*sin(s_w[i-count1][2*j]);
                df_sw[1]=stepH*cos(s_w[i-count1][2*j]);
                df_sw[2]=ratio_f[2]/2*sin(s_w[i-count1][2*j]);

                ddf_sw[0]=ratio_f[0]/2*cos(s_w[i-count1][2*j]);
                ddf_sw[1]=-stepH*sin(s_w[i-count1][2*j]);
                ddf_sw[2]=ratio_f[2]/2*cos(s_w[i-count1][2*j]);

                pEE[6*j]=initPee[6*j]*cos(-stepBeta/4)+initPee[6*j+2]*sin(-stepBeta/4)+f_sw[0];
                pEE[6*j+1]=initPee[6*j+1]+f_sw[1];
                pEE[6*j+2]=-initPee[6*j]*sin(-stepBeta/4)+initPee[6*j+2]*cos(-stepBeta/4)+f_sw[2];

                pEE[6*j+3]=initPee[6*j+3]*cos(stepBeta/4)+initPee[6*j+5]*sin(stepBeta/4);
                pEE[6*j+4]=initPee[6*j+4];
                pEE[6*j+5]=-initPee[6*j+3]*sin(stepBeta/4)+initPee[6*j+5]*cos(stepBeta/4);

                vEE[6*j]=df_sw[0]*real_ds_scale[i-count1][2*j];
                vEE[6*j+1]=df_sw[1]*real_ds_scale[i-count1][2*j];
                vEE[6*j+2]=df_sw[2]*real_ds_scale[i-count1][2*j];

                aEE[6*j]=ddf_sw[0]*real_ds_scale[i-count1][2*j]*real_ds_scale[i-count1][2*j]+df_sw[0]*real_dds_scale[i-count1][2*j];
                aEE[6*j+1]=ddf_sw[1]*real_ds_scale[i-count1][2*j]*real_ds_scale[i-count1][2*j]+df_sw[1]*real_dds_scale[i-count1][2*j];
                aEE[6*j+2]=ddf_sw[2]*real_ds_scale[i-count1][2*j]*real_ds_scale[i-count1][2*j]+df_sw[2]*real_dds_scale[i-count1][2*j];
            }

            rbt.SetPeb(pEB,"213");
            rbt.SetVb(vEB);
            rbt.SetAb(aEB);
            rbt.SetPee(pEE);
            rbt.SetVee(vEE);
            rbt.SetAee(aEE);

            rbt.GetPee(*Pee_s+18*i,rbt.body());
            rbt.GetPin(*Pin_s+18*i);
            rbt.GetVin(*Vin_s+18*i);
            rbt.GetAin(*Ain_s+18*i);
        }

        std::fill_n(vEE,18,0);
        std::fill_n(aEE,18,0);
        for(int i=count2+1;i<count3;i++)
        {
            pEB[0]=initPeb[0]+pva_b[i][2];
            pEB[2]=initPeb[2]+pva_b[i][1];
            pEB[3]=initPeb[3]+pva_b[i][0];
            vEB[0]=pva_b[i][5];
            vEB[2]=pva_b[i][4];
            vEB[3]=pva_b[i][3];
            aEB[0]=pva_b[i][8];
            aEB[2]=pva_b[i][7];
            aEB[3]=pva_b[i][6];
            for(int j=0;j<3;j++)
            {
                pEE[6*j]=initPee[6*j]*cos(3*stepBeta/4)+initPee[6*j+2]*sin(3*stepBeta/4);
                pEE[6*j+1]=initPee[6*j+1];
                pEE[6*j+2]=-initPee[6*j]*sin(3*stepBeta/4)+initPee[6*j+2]*cos(3*stepBeta/4);

                pEE[6*j+3]=initPee[6*j+3]*cos(stepBeta/4)+initPee[6*j+5]*sin(stepBeta/4);
                pEE[6*j+4]=initPee[6*j+4];
                pEE[6*j+5]=-initPee[6*j+3]*sin(stepBeta/4)+initPee[6*j+5]*cos(stepBeta/4);
            }
            rbt.SetPeb(pEB,"213");
            rbt.SetVb(vEB);
            rbt.SetAb(aEB);
            rbt.SetPee(pEE);
            rbt.SetVee(vEE);
            rbt.SetAee(aEE);

            rbt.GetPee(*Pee_s+18*i,rbt.body());
            rbt.GetPin(*Pin_s+18*i);
            rbt.GetVin(*Vin_s+18*i);
            rbt.GetAin(*Ain_s+18*i);
        }

        std::fill_n(vEE,18,0);
        std::fill_n(aEE,18,0);
        for(int i=count3;i<count4+1;i++)
        {
            pEB[0]=initPeb[0]+pva_b[i][2];
            pEB[2]=initPeb[2]+pva_b[i][1];
            pEB[3]=initPeb[3]+pva_b[i][0];
            vEB[0]=pva_b[i][5];
            vEB[2]=pva_b[i][4];
            vEB[3]=pva_b[i][3];
            aEB[0]=pva_b[i][8];
            aEB[2]=pva_b[i][7];
            aEB[3]=pva_b[i][6];
            for(int j=0;j<3;j++)
            {
                ratio_f[0]=initPee[6*j+3]*(cos(5*stepBeta/4)-cos(stepBeta/4))+initPee[6*j+5]*(sin(5*stepBeta/4)-sin(stepBeta/4));
                ratio_f[2]=-initPee[6*j+3]*(sin(5*stepBeta/4)-sin(stepBeta/4))+initPee[6*j+5]*(cos(5*stepBeta/4)-cos(stepBeta/4));

                f_sw[0]=-ratio_f[0]/2*(cos(s_w[i-count1][2*j])-1);
                f_sw[1]=stepH*sin(s_w[i-count1][2*j]);
                f_sw[2]=-ratio_f[2]/2*(cos(s_w[i-count1][2*j])-1);

                df_sw[0]=ratio_f[0]/2*sin(s_w[i-count1][2*j]);
                df_sw[1]=stepH*cos(s_w[i-count1][2*j]);
                df_sw[2]=ratio_f[2]/2*sin(s_w[i-count1][2*j]);

                ddf_sw[0]=ratio_f[0]/2*cos(s_w[i-count1][2*j]);
                ddf_sw[1]=-stepH*sin(s_w[i-count1][2*j]);
                ddf_sw[2]=ratio_f[2]/2*cos(s_w[i-count1][2*j]);

                pEE[6*j]=initPee[6*j]*cos(3*stepBeta/4)+initPee[6*j+2]*sin(3*stepBeta/4);
                pEE[6*j+1]=initPee[6*j+1];
                pEE[6*j+2]=-initPee[6*j]*sin(3*stepBeta/4)+initPee[6*j+2]*cos(3*stepBeta/4);

                pEE[6*j+3]=initPee[6*j+3]*cos(stepBeta/4)+initPee[6*j+5]*sin(stepBeta/4)+f_sw[0];
                pEE[6*j+4]=initPee[6*j+4]+f_sw[1];
                pEE[6*j+5]=-initPee[6*j+3]*sin(stepBeta/4)+initPee[6*j+5]*cos(stepBeta/4)+f_sw[2];

                vEE[6*j+3]=df_sw[0]*real_ds_scale[i-count3][2*j+1];
                vEE[6*j+4]=df_sw[1]*real_ds_scale[i-count3][2*j+1];
                vEE[6*j+5]=df_sw[2]*real_ds_scale[i-count3][2*j+1];

                aEE[6*j+3]=ddf_sw[0]*real_ds_scale[i-count3][2*j+1]*real_ds_scale[i-count3][2*j+1]+df_sw[0]*real_dds_scale[i-count3][2*j+1];
                aEE[6*j+4]=ddf_sw[1]*real_ds_scale[i-count3][2*j+1]*real_ds_scale[i-count3][2*j+1]+df_sw[1]*real_dds_scale[i-count3][2*j+1];
                aEE[6*j+5]=ddf_sw[2]*real_ds_scale[i-count3][2*j+1]*real_ds_scale[i-count3][2*j+1]+df_sw[2]*real_dds_scale[i-count3][2*j+1];
            }
            rbt.SetPeb(pEB,"213");
            rbt.SetVb(vEB);
            rbt.SetAb(aEB);
            rbt.SetPee(pEE);
            rbt.SetVee(vEE);
            rbt.SetAee(aEE);

            rbt.GetPee(*Pee_s+18*i,rbt.body());
            rbt.GetPin(*Pin_s+18*i);
            rbt.GetVin(*Vin_s+18*i);
            rbt.GetAin(*Ain_s+18*i);
        }

        std::fill_n(vEE,18,0);
        std::fill_n(aEE,18,0);
        for(int i=count4+1;i<count5+1;i++)
        {
            pEB[0]=initPeb[0]+pva_b[i][2];
            pEB[2]=initPeb[2]+pva_b[i][1];
            pEB[3]=initPeb[3]+pva_b[i][0];
            vEB[0]=pva_b[i][5];
            vEB[2]=pva_b[i][4];
            vEB[3]=pva_b[i][3];
            aEB[0]=pva_b[i][8];
            aEB[2]=pva_b[i][7];
            aEB[3]=pva_b[i][6];
            for(int j=0;j<3;j++)
            {

                pEE[6*j]=initPee[6*j]*cos(3*stepBeta/4)+initPee[6*j+2]*sin(3*stepBeta/4);
                pEE[6*j+1]=initPee[6*j+1];
                pEE[6*j+2]=-initPee[6*j]*sin(3*stepBeta/4)+initPee[6*j+2]*cos(3*stepBeta/4);

                pEE[6*j+3]=initPee[6*j+3]*cos(5*stepBeta/4)+initPee[6*j+5]*sin(5*stepBeta/4);
                pEE[6*j+4]=initPee[6*j+4];
                pEE[6*j+5]=-initPee[6*j+3]*sin(5*stepBeta/4)+initPee[6*j+5]*cos(5*stepBeta/4);
            }
            rbt.SetPeb(pEB,"213");
            rbt.SetVb(vEB);
            rbt.SetAb(aEB);
            rbt.SetPee(pEE);
            rbt.SetVee(vEE);
            rbt.SetAee(aEE);

            rbt.GetPee(*Pee_s+18*i,rbt.body());
            rbt.GetPin(*Pin_s+18*i);
            rbt.GetVin(*Vin_s+18*i);
            rbt.GetAin(*Ain_s+18*i);
        }
        printf("finish GetOptimalGait2s\n");
    }

	/*获取位置-时间轨迹*/
    void NonRTOptimalGCSRotate::GetOptimalGait2t()
    {
        printf("start GetOptimalGait2t\n");
        double * pva_b_t=new double [9*stepTimeCount];
        double * s_w_t=new double [6*swingTimeCount];
        double * ds_w_t=new double [6*swingTimeCount];
        double * dds_w_t=new double [6*swingTimeCount];
        double * Pee_t=new double [18*stepTimeCount];
        double * Pin_t=new double [18*stepTimeCount];
        double * Vin_t=new double [18*stepTimeCount];
        double * Ain_t=new double [18*stepTimeCount];

        double * VeeMinus_t=new double [18*stepTimeCount];
        double * VeeB_t=new double [18*stepTimeCount];
        double * AeeMinus_t=new double [18*stepTimeCount];
        double * AeeB_t=new double [18*stepTimeCount];

        int k_start {0};
        for(int i=0;i<stepTimeCount;i++)
        {
            for(int k=k_start;k<count5+1;k++)
            {
                if(0.001*i>=timeArray_body_tmp[k] && 0.001*i<timeArray_body_tmp[k+1])
                {
                    k_start=k;
                    for(int j=0;j<9;j++)
                    {
                        *(pva_b_t+9*i+j)=pva_b[k][j]+(pva_b[k+1][j]-pva_b[k][j])
                                *(0.001*i-timeArray_body_tmp[k])/(timeArray_body_tmp[k+1]-timeArray_body_tmp[k]);
                    }
                    break;
                }
            }
        }
        for(int j=0;j<6;j++)
        {
            k_start=0;
            for(int i=0;i<swingTimeCount;i++)
            {
                for(int k=k_start;k<swingSCount+1;k++)
                {
                    if(0.001*i>=timeArray_tmp[k][j] && 0.001*i<timeArray_tmp[k+1][j])
                    {
                        k_start=k;
                        *(s_w_t+6*i+j)=s_w[k][j]+(s_w[k+1][j]-s_w[k][j])
                                *(0.001*i-timeArray_tmp[k][j])/(timeArray_tmp[k+1][j]-timeArray_tmp[k][j]);
                        *(ds_w_t+6*i+j)=real_ds_scale[k][j]+(real_ds_scale[k+1][j]-real_ds_scale[k][j])
                                *(0.001*i-timeArray_tmp[k][j])/(timeArray_tmp[k+1][j]-timeArray_tmp[k][j]);
                        *(dds_w_t+6*i+j)=real_dds_scale[k][j]+(real_dds_scale[k+1][j]-real_dds_scale[k][j])
                                *(0.001*i-timeArray_tmp[k][j])/(timeArray_tmp[k+1][j]-timeArray_tmp[k][j]);
                        break;
                    }
                }
            }
        }

        double pEB[6] {0};
        double vEB[6] {0};
        double aEB[6] {0};
        double pEE[18] {0};
        double vEE[18] {0};
        double aEE[18] {0};
        double ratio_f[3] {0};

        std::fill_n(vEE,18,0);
        std::fill_n(aEE,18,0);
        int iStart;
        int iEnd;
        iEnd=round(swingTimeCount/(1-dutyCycle)*(2*dutyCycle-1)/4);
        for (int i=0;i<iEnd;i++)
        {
            pEB[0]=initPeb[0]+*(pva_b_t+9*i+2);
            pEB[2]=initPeb[2]+*(pva_b_t+9*i+1);
            pEB[3]=initPeb[3]+*(pva_b_t+9*i+0);
            vEB[0]=*(pva_b_t+9*i+5);
            vEB[2]=*(pva_b_t+9*i+4);
            vEB[3]=*(pva_b_t+9*i+3);
            aEB[0]=*(pva_b_t+9*i+8);
            aEB[2]=*(pva_b_t+9*i+7);
            aEB[3]=*(pva_b_t+9*i+6);
            for(int j=0;j<3;j++)
            {
                pEE[6*j]=initPee[6*j]*cos(-stepBeta/4)+initPee[6*j+2]*sin(-stepBeta/4);
                pEE[6*j+1]=initPee[6*j+1];
                pEE[6*j+2]=-initPee[6*j]*sin(-stepBeta/4)+initPee[6*j+2]*cos(-stepBeta/4);

                pEE[6*j+3]=initPee[6*j+3]*cos(stepBeta/4)+initPee[6*j+5]*sin(stepBeta/4);
                pEE[6*j+4]=initPee[6*j+4];
                pEE[6*j+5]=-initPee[6*j+3]*sin(stepBeta/4)+initPee[6*j+5]*cos(stepBeta/4);
            }
            rbt.SetPeb(pEB,"213");
            rbt.SetVb(vEB);
            rbt.SetAb(aEB);
            rbt.SetPee(pEE);
            rbt.SetVee(vEE);
            rbt.SetAee(aEE);

            for(int j=0;j<6;j++)
            {
                *(VeeMinus_t+18*i+3*j)=vEE[3*j];
                *(VeeMinus_t+18*i+3*j+1)=vEE[3*j+1];
                *(VeeMinus_t+18*i+3*j+2)=vEE[3*j+2]-vEB[2];
                *(AeeMinus_t+18*i+3*j)=aEE[3*j];
                *(AeeMinus_t+18*i+3*j+1)=aEE[3*j+1];
                *(AeeMinus_t+18*i+3*j+2)=aEE[3*j+2]-aEB[2];
            }
            rbt.GetVee(VeeB_t+18*i,rbt.body());
            rbt.GetAee(AeeB_t+18*i,rbt.body());

            rbt.GetPee(Pee_t+18*i,rbt.body());
            rbt.GetPin(Pin_t+18*i);
            rbt.GetVin(Vin_t+18*i);
            rbt.GetAin(Ain_t+18*i);
        }

        std::fill_n(vEE,18,0);
        std::fill_n(aEE,18,0);
        iStart=iEnd;
        iEnd=round(swingTimeCount/(1-dutyCycle)*(2*dutyCycle-1)/4+swingTimeCount);
        if(iEnd-iStart>swingTimeCount)
        {
            printf("WARNING!!! s_w short than swingTimeCount!\n");
            iEnd=iStart+swingTimeCount;
        }
        for (int i=iStart;i<iEnd;i++)
        {
            pEB[0]=initPeb[0]+*(pva_b_t+9*i+2);
            pEB[2]=initPeb[2]+*(pva_b_t+9*i+1);
            pEB[3]=initPeb[3]+*(pva_b_t+9*i+0);
            vEB[0]=*(pva_b_t+9*i+5);
            vEB[2]=*(pva_b_t+9*i+4);
            vEB[3]=*(pva_b_t+9*i+3);
            aEB[0]=*(pva_b_t+9*i+8);
            aEB[2]=*(pva_b_t+9*i+7);
            aEB[3]=*(pva_b_t+9*i+6);
            for(int j=0;j<3;j++)
            {
                ratio_f[0]=initPee[6*j]*(cos(3*stepBeta/4)-cos(-stepBeta/4))+initPee[6*j+2]*(sin(3*stepBeta/4)-sin(-stepBeta/4));
                ratio_f[2]=-initPee[6*j]*(sin(3*stepBeta/4)-sin(-stepBeta/4))+initPee[6*j+2]*(cos(3*stepBeta/4)-cos(-stepBeta/4));

                f_sw[0]=-ratio_f[0]/2*(cos(*(s_w_t+6*(i-iStart)+2*j))-1);
                f_sw[1]=stepH*sin(*(s_w_t+6*(i-iStart)+2*j));
                f_sw[2]=-ratio_f[2]/2*(cos(*(s_w_t+6*(i-iStart)+2*j))-1);

                df_sw[0]=ratio_f[0]/2*sin(*(s_w_t+6*(i-iStart)+2*j));
                df_sw[1]=stepH*cos(*(s_w_t+6*(i-iStart)+2*j));
                df_sw[2]=ratio_f[2]/2*sin(*(s_w_t+6*(i-iStart)+2*j));

                ddf_sw[0]=ratio_f[0]/2*cos(*(s_w_t+6*(i-iStart)+2*j));
                ddf_sw[1]=-stepH*sin(*(s_w_t+6*(i-iStart)+2*j));
                ddf_sw[2]=ratio_f[2]/2*cos(*(s_w_t+6*(i-iStart)+2*j));

                pEE[6*j]=initPee[6*j]*cos(-stepBeta/4)+initPee[6*j+2]*sin(-stepBeta/4)+f_sw[0];
                pEE[6*j+1]=initPee[6*j+1]+f_sw[1];
                pEE[6*j+2]=-initPee[6*j]*sin(-stepBeta/4)+initPee[6*j+2]*cos(-stepBeta/4)+f_sw[2];

                pEE[6*j+3]=initPee[6*j+3]*cos(stepBeta/4)+initPee[6*j+5]*sin(stepBeta/4);
                pEE[6*j+4]=initPee[6*j+4];
                pEE[6*j+5]=-initPee[6*j+3]*sin(stepBeta/4)+initPee[6*j+5]*cos(stepBeta/4);

                vEE[6*j]=df_sw[0]*(*(ds_w_t+6*(i-iStart)+2*j));
                vEE[6*j+1]=df_sw[1]*(*(ds_w_t+6*(i-iStart)+2*j));
                vEE[6*j+2]=df_sw[2]*(*(ds_w_t+6*(i-iStart)+2*j));

                aEE[6*j]=ddf_sw[0]*(*(ds_w_t+6*(i-iStart)+2*j))*(*(ds_w_t+6*(i-iStart)+2*j))+df_sw[0]*(*(dds_w_t+6*(i-iStart)+2*j));
                aEE[6*j+1]=ddf_sw[1]*(*(ds_w_t+6*(i-iStart)+2*j))*(*(ds_w_t+6*(i-iStart)+2*j))+df_sw[1]*(*(dds_w_t+6*(i-iStart)+2*j));
                aEE[6*j+2]=ddf_sw[2]*(*(ds_w_t+6*(i-iStart)+2*j))*(*(ds_w_t+6*(i-iStart)+2*j))+df_sw[2]*(*(dds_w_t+6*(i-iStart)+2*j));
            }
            rbt.SetPeb(pEB,"213");
            rbt.SetVb(vEB);
            rbt.SetAb(aEB);
            rbt.SetPee(pEE);
            rbt.SetVee(vEE);
            rbt.SetAee(aEE);

            for(int j=0;j<6;j++)
            {
                *(VeeMinus_t+18*i+3*j)=vEE[3*j];
                *(VeeMinus_t+18*i+3*j+1)=vEE[3*j+1];
                *(VeeMinus_t+18*i+3*j+2)=vEE[3*j+2]-vEB[2];
                *(AeeMinus_t+18*i+3*j)=aEE[3*j];
                *(AeeMinus_t+18*i+3*j+1)=aEE[3*j+1];
                *(AeeMinus_t+18*i+3*j+2)=aEE[3*j+2]-aEB[2];
            }
            rbt.GetVee(VeeB_t+18*i,rbt.body());
            rbt.GetAee(AeeB_t+18*i,rbt.body());

            rbt.GetPee(Pee_t+18*i,rbt.body());
            rbt.GetPin(Pin_t+18*i);
            rbt.GetVin(Vin_t+18*i);
            rbt.GetAin(Ain_t+18*i);
        }

        std::fill_n(vEE,18,0);
        std::fill_n(aEE,18,0);
        iStart=iEnd;
        iEnd=round(3*swingTimeCount/(1-dutyCycle)*(2*dutyCycle-1)/4+swingTimeCount);
        for(int i=iStart;i<iEnd;i++)
        {
            pEB[0]=initPeb[0]+*(pva_b_t+9*i+2);
            pEB[2]=initPeb[2]+*(pva_b_t+9*i+1);
            pEB[3]=initPeb[3]+*(pva_b_t+9*i+0);
            vEB[0]=*(pva_b_t+9*i+5);
            vEB[2]=*(pva_b_t+9*i+4);
            vEB[3]=*(pva_b_t+9*i+3);
            aEB[0]=*(pva_b_t+9*i+8);
            aEB[2]=*(pva_b_t+9*i+7);
            aEB[3]=*(pva_b_t+9*i+6);
            for(int j=0;j<3;j++)
            {
                pEE[6*j]=initPee[6*j]*cos(3*stepBeta/4)+initPee[6*j+2]*sin(3*stepBeta/4);
                pEE[6*j+1]=initPee[6*j+1];
                pEE[6*j+2]=-initPee[6*j]*sin(3*stepBeta/4)+initPee[6*j+2]*cos(3*stepBeta/4);

                pEE[6*j+3]=initPee[6*j+3]*cos(stepBeta/4)+initPee[6*j+5]*sin(stepBeta/4);
                pEE[6*j+4]=initPee[6*j+4];
                pEE[6*j+5]=-initPee[6*j+3]*sin(stepBeta/4)+initPee[6*j+5]*cos(stepBeta/4);
            }
            rbt.SetPeb(pEB,"213");
            rbt.SetVb(vEB);
            rbt.SetAb(aEB);
            rbt.SetPee(pEE);
            rbt.SetVee(vEE);
            rbt.SetAee(aEE);

            for(int j=0;j<6;j++)
            {
                *(VeeMinus_t+18*i+3*j)=vEE[3*j];
                *(VeeMinus_t+18*i+3*j+1)=vEE[3*j+1];
                *(VeeMinus_t+18*i+3*j+2)=vEE[3*j+2]-vEB[2];
                *(AeeMinus_t+18*i+3*j)=aEE[3*j];
                *(AeeMinus_t+18*i+3*j+1)=aEE[3*j+1];
                *(AeeMinus_t+18*i+3*j+2)=aEE[3*j+2]-aEB[2];
            }
            rbt.GetVee(VeeB_t+18*i,rbt.body());
            rbt.GetAee(AeeB_t+18*i,rbt.body());

            rbt.GetPee(Pee_t+18*i,rbt.body());
            rbt.GetPin(Pin_t+18*i);
            rbt.GetVin(Vin_t+18*i);
            rbt.GetAin(Ain_t+18*i);
        }

        std::fill_n(vEE,18,0);
        std::fill_n(aEE,18,0);
        iStart=iEnd;
        iEnd=round(3*swingTimeCount/(1-dutyCycle)*(2*dutyCycle-1)/4+2*swingTimeCount);
        if(iEnd-iStart>swingTimeCount)
        {
            printf("WARNING!!! s_w short than swingTimeCount!\n");
            iEnd=iStart+swingTimeCount;
        }
        for(int i=iStart;i<iEnd;i++)
        {
            pEB[0]=initPeb[0]+*(pva_b_t+9*i+2);
            pEB[2]=initPeb[2]+*(pva_b_t+9*i+1);
            pEB[3]=initPeb[3]+*(pva_b_t+9*i+0);
            vEB[0]=*(pva_b_t+9*i+5);
            vEB[2]=*(pva_b_t+9*i+4);
            vEB[3]=*(pva_b_t+9*i+3);
            aEB[0]=*(pva_b_t+9*i+8);
            aEB[2]=*(pva_b_t+9*i+7);
            aEB[3]=*(pva_b_t+9*i+6);
            for(int j=0;j<3;j++)
            {
                ratio_f[0]=initPee[6*j+3]*(cos(5*stepBeta/4)-cos(stepBeta/4))+initPee[6*j+5]*(sin(5*stepBeta/4)-sin(stepBeta/4));
                ratio_f[2]=-initPee[6*j+3]*(sin(5*stepBeta/4)-sin(stepBeta/4))+initPee[6*j+5]*(cos(5*stepBeta/4)-cos(stepBeta/4));

                f_sw[0]=-ratio_f[0]/2*(cos(*(s_w_t+6*(i-iStart)+2*j+1))-1);
                f_sw[1]=stepH*sin(*(s_w_t+6*(i-iStart)+2*j+1));
                f_sw[2]=-ratio_f[2]/2*(cos(*(s_w_t+6*(i-iStart)+2*j+1))-1);

                df_sw[0]=ratio_f[0]/2*sin(*(s_w_t+6*(i-iStart)+2*j+1));
                df_sw[1]=stepH*cos(*(s_w_t+6*(i-iStart)+2*j+1));
                df_sw[2]=ratio_f[2]/2*sin(*(s_w_t+6*(i-iStart)+2*j+1));

                ddf_sw[0]=ratio_f[0]/2*cos(*(s_w_t+6*(i-iStart)+2*j+1));
                ddf_sw[1]=-stepH*sin(*(s_w_t+6*(i-iStart)+2*j+1));
                ddf_sw[2]=ratio_f[2]/2*cos(*(s_w_t+6*(i-iStart)+2*j+1));

                pEE[6*j]=initPee[6*j]*cos(3*stepBeta/4)+initPee[6*j+2]*sin(3*stepBeta/4);
                pEE[6*j+1]=initPee[6*j+1];
                pEE[6*j+2]=-initPee[6*j]*sin(3*stepBeta/4)+initPee[6*j+2]*cos(3*stepBeta/4);

                pEE[6*j+3]=initPee[6*j+3]*cos(stepBeta/4)+initPee[6*j+5]*sin(stepBeta/4)+f_sw[0];
                pEE[6*j+4]=initPee[6*j+4]+f_sw[1];
                pEE[6*j+5]=-initPee[6*j+3]*sin(stepBeta/4)+initPee[6*j+5]*cos(stepBeta/4)+f_sw[2];

                vEE[6*j+3]=df_sw[0]*(*(ds_w_t+6*(i-iStart)+2*j+1));
                vEE[6*j+4]=df_sw[1]*(*(ds_w_t+6*(i-iStart)+2*j+1));
                vEE[6*j+5]=df_sw[2]*(*(ds_w_t+6*(i-iStart)+2*j+1));

                aEE[6*j+3]=ddf_sw[0]*(*(ds_w_t+6*(i-iStart)+2*j+1))*(*(ds_w_t+6*(i-iStart)+2*j+1))+df_sw[0]*(*(dds_w_t+6*(i-iStart)+2*j+1));
                aEE[6*j+4]=ddf_sw[1]*(*(ds_w_t+6*(i-iStart)+2*j+1))*(*(ds_w_t+6*(i-iStart)+2*j+1))+df_sw[1]*(*(dds_w_t+6*(i-iStart)+2*j+1));
                aEE[6*j+5]=ddf_sw[2]*(*(ds_w_t+6*(i-iStart)+2*j+1))*(*(ds_w_t+6*(i-iStart)+2*j+1))+df_sw[2]*(*(dds_w_t+6*(i-iStart)+2*j+1));
            }
            rbt.SetPeb(pEB,"213");
            rbt.SetVb(vEB);
            rbt.SetAb(aEB);
            rbt.SetPee(pEE);
            rbt.SetVee(vEE);
            rbt.SetAee(aEE);

            for(int j=0;j<6;j++)
            {
                *(VeeMinus_t+18*i+3*j)=vEE[3*j];
                *(VeeMinus_t+18*i+3*j+1)=vEE[3*j+1];
                *(VeeMinus_t+18*i+3*j+2)=vEE[3*j+2]-vEB[2];
                *(AeeMinus_t+18*i+3*j)=aEE[3*j];
                *(AeeMinus_t+18*i+3*j+1)=aEE[3*j+1];
                *(AeeMinus_t+18*i+3*j+2)=aEE[3*j+2]-aEB[2];
            }
            rbt.GetVee(VeeB_t+18*i,rbt.body());
            rbt.GetAee(AeeB_t+18*i,rbt.body());

            rbt.GetPee(Pee_t+18*i,rbt.body());
            rbt.GetPin(Pin_t+18*i);
            rbt.GetVin(Vin_t+18*i);
            rbt.GetAin(Ain_t+18*i);
        }

        std::fill_n(vEE,18,0);
        std::fill_n(aEE,18,0);
        iStart=iEnd;
        iEnd=stepTimeCount;
        for(int i=iStart;i<iEnd;i++)
        {
            pEB[0]=initPeb[0]+*(pva_b_t+9*i+2);
            pEB[2]=initPeb[2]+*(pva_b_t+9*i+1);
            pEB[3]=initPeb[3]+*(pva_b_t+9*i+0);
            vEB[0]=*(pva_b_t+9*i+5);
            vEB[2]=*(pva_b_t+9*i+4);
            vEB[3]=*(pva_b_t+9*i+3);
            aEB[0]=*(pva_b_t+9*i+8);
            aEB[2]=*(pva_b_t+9*i+7);
            aEB[3]=*(pva_b_t+9*i+6);
            for(int j=0;j<3;j++)
            {
                pEE[6*j]=initPee[6*j]*cos(3*stepBeta/4)+initPee[6*j+2]*sin(3*stepBeta/4);
                pEE[6*j+1]=initPee[6*j+1];
                pEE[6*j+2]=-initPee[6*j]*sin(3*stepBeta/4)+initPee[6*j+2]*cos(3*stepBeta/4);

                pEE[6*j+3]=initPee[6*j+3]*cos(5*stepBeta/4)+initPee[6*j+5]*sin(5*stepBeta/4);
                pEE[6*j+4]=initPee[6*j+4];
                pEE[6*j+5]=-initPee[6*j+3]*sin(5*stepBeta/4)+initPee[6*j+5]*cos(5*stepBeta/4);
            }
            rbt.SetPeb(pEB,"213");
            rbt.SetVb(vEB);
            rbt.SetAb(aEB);
            rbt.SetPee(pEE);
            rbt.SetVee(vEE);
            rbt.SetAee(aEE);

            for(int j=0;j<6;j++)
            {
                *(VeeMinus_t+18*i+3*j)=vEE[3*j];
                *(VeeMinus_t+18*i+3*j+1)=vEE[3*j+1];
                *(VeeMinus_t+18*i+3*j+2)=vEE[3*j+2]-vEB[2];
                *(AeeMinus_t+18*i+3*j)=aEE[3*j];
                *(AeeMinus_t+18*i+3*j+1)=aEE[3*j+1];
                *(AeeMinus_t+18*i+3*j+2)=aEE[3*j+2]-aEB[2];
            }
            rbt.GetVee(VeeB_t+18*i,rbt.body());
            rbt.GetAee(AeeB_t+18*i,rbt.body());

            rbt.GetPee(Pee_t+18*i,rbt.body());
            rbt.GetPin(Pin_t+18*i);
            rbt.GetVin(Vin_t+18*i);
            rbt.GetAin(Ain_t+18*i);
        }

        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/pva_b_t.txt",pva_b_t,stepTimeCount,9);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/s_w_t.txt",s_w_t,swingTimeCount,6);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Pee_t.txt",Pee_t,stepTimeCount,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Pin_t.txt",Pin_t,stepTimeCount,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Vin_t.txt",Vin_t,stepTimeCount,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Ain_t.txt",Ain_t,stepTimeCount,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/VeeB_t.txt",VeeB_t,stepTimeCount,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/VeeMinus_t.txt",VeeMinus_t,stepTimeCount,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/AeeB_t.txt",AeeB_t,stepTimeCount,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/AeeMinus_t.txt",AeeMinus_t,stepTimeCount,18);
        delete [] pva_b_t;
        delete [] s_w_t;
        delete [] Pee_t;
        delete [] Pin_t;
        delete [] Vin_t;
        delete [] Ain_t;
        delete [] AeeB_t;
        delete [] AeeMinus_t;
        delete [] VeeB_t;
        delete [] VeeMinus_t;
        printf("finish GetOptimalGait2t\n");
    }

	/*获取位置-时间轨迹，重构，可以输出参数*/
    void NonRTOptimalGCSRotate::GetOptimalGait2t(double *out_tippos, double &out_bodyvel, double &out_period)
    {
        printf("start GetOptimalGait2t\n");
        double * pva_b_t=new double [9*stepTimeCount];
        double * s_w_t=new double [6*swingTimeCount];
        double * ds_w_t=new double [6*swingTimeCount];
        double * dds_w_t=new double [6*swingTimeCount];
        double * Pee_t=new double [18*stepTimeCount];
        double * Pin_t=new double [18*stepTimeCount];
        double * Vin_t=new double [18*stepTimeCount];
        double * Ain_t=new double [18*stepTimeCount];

        double * VeeMinus_t=new double [18*stepTimeCount];
        double * VeeB_t=new double [18*stepTimeCount];
        double * AeeMinus_t=new double [18*stepTimeCount];
        double * AeeB_t=new double [18*stepTimeCount];

        int k_start {0};
        for(int i=0;i<stepTimeCount;i++)
        {
            for(int k=k_start;k<count5+1;k++)
            {
                if(0.001*i>=timeArray_body_tmp[k] && 0.001*i<timeArray_body_tmp[k+1])
                {
                    k_start=k;
                    for(int j=0;j<9;j++)
                    {
                        *(pva_b_t+9*i+j)=pva_b[k][j]+(pva_b[k+1][j]-pva_b[k][j])
                                *(0.001*i-timeArray_body_tmp[k])/(timeArray_body_tmp[k+1]-timeArray_body_tmp[k]);
                    }
                    break;
                }
            }
        }
        for(int j=0;j<6;j++)
        {
            k_start=0;
            for(int i=0;i<swingTimeCount;i++)
            {
                for(int k=k_start;k<swingSCount+1;k++)
                {
                    if(0.001*i>=timeArray_tmp[k][j] && 0.001*i<timeArray_tmp[k+1][j])
                    {
                        k_start=k;
                        *(s_w_t+6*i+j)=s_w[k][j]+(s_w[k+1][j]-s_w[k][j])
                                *(0.001*i-timeArray_tmp[k][j])/(timeArray_tmp[k+1][j]-timeArray_tmp[k][j]);
                        *(ds_w_t+6*i+j)=real_ds_scale[k][j]+(real_ds_scale[k+1][j]-real_ds_scale[k][j])
                                *(0.001*i-timeArray_tmp[k][j])/(timeArray_tmp[k+1][j]-timeArray_tmp[k][j]);
                        *(dds_w_t+6*i+j)=real_dds_scale[k][j]+(real_dds_scale[k+1][j]-real_dds_scale[k][j])
                                *(0.001*i-timeArray_tmp[k][j])/(timeArray_tmp[k+1][j]-timeArray_tmp[k][j]);
                        break;
                    }
                }
            }
        }

        double pEB[6] {0};
        double vEB[6] {0};
        double aEB[6] {0};
        double pEE[18] {0};
        double vEE[18] {0};
        double aEE[18] {0};
        double ratio_f[3] {0};

        std::fill_n(vEE,18,0);
        std::fill_n(aEE,18,0);
        int iStart;
        int iEnd;
        iEnd=round(swingTimeCount/(1-dutyCycle)*(2*dutyCycle-1)/4);
        for (int i=0;i<iEnd;i++)
        {
            pEB[0]=initPeb[0]+*(pva_b_t+9*i+2);
            pEB[2]=initPeb[2]+*(pva_b_t+9*i+1);
            pEB[3]=initPeb[3]+*(pva_b_t+9*i+0);
            vEB[0]=*(pva_b_t+9*i+5);
            vEB[2]=*(pva_b_t+9*i+4);
            vEB[3]=*(pva_b_t+9*i+3);
            aEB[0]=*(pva_b_t+9*i+8);
            aEB[2]=*(pva_b_t+9*i+7);
            aEB[3]=*(pva_b_t+9*i+6);
            for(int j=0;j<3;j++)
            {
                pEE[6*j]=initPee[6*j]*cos(-stepBeta/4)+initPee[6*j+2]*sin(-stepBeta/4);
                pEE[6*j+1]=initPee[6*j+1];
                pEE[6*j+2]=-initPee[6*j]*sin(-stepBeta/4)+initPee[6*j+2]*cos(-stepBeta/4);

                pEE[6*j+3]=initPee[6*j+3]*cos(stepBeta/4)+initPee[6*j+5]*sin(stepBeta/4);
                pEE[6*j+4]=initPee[6*j+4];
                pEE[6*j+5]=-initPee[6*j+3]*sin(stepBeta/4)+initPee[6*j+5]*cos(stepBeta/4);
            }
            rbt.SetPeb(pEB,"213");
            rbt.SetVb(vEB);
            rbt.SetAb(aEB);
            rbt.SetPee(pEE);
            rbt.SetVee(vEE);
            rbt.SetAee(aEE);

            for(int j=0;j<6;j++)
            {
                *(VeeMinus_t+18*i+3*j)=vEE[3*j];
                *(VeeMinus_t+18*i+3*j+1)=vEE[3*j+1];
                *(VeeMinus_t+18*i+3*j+2)=vEE[3*j+2]-vEB[2];
                *(AeeMinus_t+18*i+3*j)=aEE[3*j];
                *(AeeMinus_t+18*i+3*j+1)=aEE[3*j+1];
                *(AeeMinus_t+18*i+3*j+2)=aEE[3*j+2]-aEB[2];
            }
            rbt.GetVee(VeeB_t+18*i,rbt.body());
            rbt.GetAee(AeeB_t+18*i,rbt.body());

            rbt.GetPee(Pee_t+18*i,rbt.body());
            rbt.GetPin(Pin_t+18*i);
            rbt.GetVin(Vin_t+18*i);
            rbt.GetAin(Ain_t+18*i);
        }

        std::fill_n(vEE,18,0);
        std::fill_n(aEE,18,0);
        iStart=iEnd;
        iEnd=round(swingTimeCount/(1-dutyCycle)*(2*dutyCycle-1)/4+swingTimeCount);
        if(iEnd-iStart>swingTimeCount)
        {
            printf("WARNING!!! s_w short than swingTimeCount!\n");
            iEnd=iStart+swingTimeCount;
        }
        for (int i=iStart;i<iEnd;i++)
        {
            pEB[0]=initPeb[0]+*(pva_b_t+9*i+2);
            pEB[2]=initPeb[2]+*(pva_b_t+9*i+1);
            pEB[3]=initPeb[3]+*(pva_b_t+9*i+0);
            vEB[0]=*(pva_b_t+9*i+5);
            vEB[2]=*(pva_b_t+9*i+4);
            vEB[3]=*(pva_b_t+9*i+3);
            aEB[0]=*(pva_b_t+9*i+8);
            aEB[2]=*(pva_b_t+9*i+7);
            aEB[3]=*(pva_b_t+9*i+6);
            for(int j=0;j<3;j++)
            {
                ratio_f[0]=initPee[6*j]*(cos(3*stepBeta/4)-cos(-stepBeta/4))+initPee[6*j+2]*(sin(3*stepBeta/4)-sin(-stepBeta/4));
                ratio_f[2]=-initPee[6*j]*(sin(3*stepBeta/4)-sin(-stepBeta/4))+initPee[6*j+2]*(cos(3*stepBeta/4)-cos(-stepBeta/4));

                f_sw[0]=-ratio_f[0]/2*(cos(*(s_w_t+6*(i-iStart)+2*j))-1);
                f_sw[1]=stepH*sin(*(s_w_t+6*(i-iStart)+2*j));
                f_sw[2]=-ratio_f[2]/2*(cos(*(s_w_t+6*(i-iStart)+2*j))-1);

                df_sw[0]=ratio_f[0]/2*sin(*(s_w_t+6*(i-iStart)+2*j));
                df_sw[1]=stepH*cos(*(s_w_t+6*(i-iStart)+2*j));
                df_sw[2]=ratio_f[2]/2*sin(*(s_w_t+6*(i-iStart)+2*j));

                ddf_sw[0]=ratio_f[0]/2*cos(*(s_w_t+6*(i-iStart)+2*j));
                ddf_sw[1]=-stepH*sin(*(s_w_t+6*(i-iStart)+2*j));
                ddf_sw[2]=ratio_f[2]/2*cos(*(s_w_t+6*(i-iStart)+2*j));

                pEE[6*j]=initPee[6*j]*cos(-stepBeta/4)+initPee[6*j+2]*sin(-stepBeta/4)+f_sw[0];
                pEE[6*j+1]=initPee[6*j+1]+f_sw[1];
                pEE[6*j+2]=-initPee[6*j]*sin(-stepBeta/4)+initPee[6*j+2]*cos(-stepBeta/4)+f_sw[2];

                pEE[6*j+3]=initPee[6*j+3]*cos(stepBeta/4)+initPee[6*j+5]*sin(stepBeta/4);
                pEE[6*j+4]=initPee[6*j+4];
                pEE[6*j+5]=-initPee[6*j+3]*sin(stepBeta/4)+initPee[6*j+5]*cos(stepBeta/4);

                vEE[6*j]=df_sw[0]*(*(ds_w_t+6*(i-iStart)+2*j));
                vEE[6*j+1]=df_sw[1]*(*(ds_w_t+6*(i-iStart)+2*j));
                vEE[6*j+2]=df_sw[2]*(*(ds_w_t+6*(i-iStart)+2*j));

                aEE[6*j]=ddf_sw[0]*(*(ds_w_t+6*(i-iStart)+2*j))*(*(ds_w_t+6*(i-iStart)+2*j))+df_sw[0]*(*(dds_w_t+6*(i-iStart)+2*j));
                aEE[6*j+1]=ddf_sw[1]*(*(ds_w_t+6*(i-iStart)+2*j))*(*(ds_w_t+6*(i-iStart)+2*j))+df_sw[1]*(*(dds_w_t+6*(i-iStart)+2*j));
                aEE[6*j+2]=ddf_sw[2]*(*(ds_w_t+6*(i-iStart)+2*j))*(*(ds_w_t+6*(i-iStart)+2*j))+df_sw[2]*(*(dds_w_t+6*(i-iStart)+2*j));
            }
            rbt.SetPeb(pEB,"213");
            rbt.SetVb(vEB);
            rbt.SetAb(aEB);
            rbt.SetPee(pEE);
            rbt.SetVee(vEE);
            rbt.SetAee(aEE);

            for(int j=0;j<6;j++)
            {
                *(VeeMinus_t+18*i+3*j)=vEE[3*j];
                *(VeeMinus_t+18*i+3*j+1)=vEE[3*j+1];
                *(VeeMinus_t+18*i+3*j+2)=vEE[3*j+2]-vEB[2];
                *(AeeMinus_t+18*i+3*j)=aEE[3*j];
                *(AeeMinus_t+18*i+3*j+1)=aEE[3*j+1];
                *(AeeMinus_t+18*i+3*j+2)=aEE[3*j+2]-aEB[2];
            }
            rbt.GetVee(VeeB_t+18*i,rbt.body());
            rbt.GetAee(AeeB_t+18*i,rbt.body());

            rbt.GetPee(Pee_t+18*i,rbt.body());
            rbt.GetPin(Pin_t+18*i);
            rbt.GetVin(Vin_t+18*i);
            rbt.GetAin(Ain_t+18*i);
        }

        std::fill_n(vEE,18,0);
        std::fill_n(aEE,18,0);
        iStart=iEnd;
        iEnd=round(3*swingTimeCount/(1-dutyCycle)*(2*dutyCycle-1)/4+swingTimeCount);
        for(int i=iStart;i<iEnd;i++)
        {
            pEB[0]=initPeb[0]+*(pva_b_t+9*i+2);
            pEB[2]=initPeb[2]+*(pva_b_t+9*i+1);
            pEB[3]=initPeb[3]+*(pva_b_t+9*i+0);
            vEB[0]=*(pva_b_t+9*i+5);
            vEB[2]=*(pva_b_t+9*i+4);
            vEB[3]=*(pva_b_t+9*i+3);
            aEB[0]=*(pva_b_t+9*i+8);
            aEB[2]=*(pva_b_t+9*i+7);
            aEB[3]=*(pva_b_t+9*i+6);
            for(int j=0;j<3;j++)
            {
                pEE[6*j]=initPee[6*j]*cos(3*stepBeta/4)+initPee[6*j+2]*sin(3*stepBeta/4);
                pEE[6*j+1]=initPee[6*j+1];
                pEE[6*j+2]=-initPee[6*j]*sin(3*stepBeta/4)+initPee[6*j+2]*cos(3*stepBeta/4);

                pEE[6*j+3]=initPee[6*j+3]*cos(stepBeta/4)+initPee[6*j+5]*sin(stepBeta/4);
                pEE[6*j+4]=initPee[6*j+4];
                pEE[6*j+5]=-initPee[6*j+3]*sin(stepBeta/4)+initPee[6*j+5]*cos(stepBeta/4);
            }
            rbt.SetPeb(pEB,"213");
            rbt.SetVb(vEB);
            rbt.SetAb(aEB);
            rbt.SetPee(pEE);
            rbt.SetVee(vEE);
            rbt.SetAee(aEE);

            for(int j=0;j<6;j++)
            {
                *(VeeMinus_t+18*i+3*j)=vEE[3*j];
                *(VeeMinus_t+18*i+3*j+1)=vEE[3*j+1];
                *(VeeMinus_t+18*i+3*j+2)=vEE[3*j+2]-vEB[2];
                *(AeeMinus_t+18*i+3*j)=aEE[3*j];
                *(AeeMinus_t+18*i+3*j+1)=aEE[3*j+1];
                *(AeeMinus_t+18*i+3*j+2)=aEE[3*j+2]-aEB[2];
            }
            rbt.GetVee(VeeB_t+18*i,rbt.body());
            rbt.GetAee(AeeB_t+18*i,rbt.body());

            rbt.GetPee(Pee_t+18*i,rbt.body());
            rbt.GetPin(Pin_t+18*i);
            rbt.GetVin(Vin_t+18*i);
            rbt.GetAin(Ain_t+18*i);
        }

        std::fill_n(vEE,18,0);
        std::fill_n(aEE,18,0);
        iStart=iEnd;
        iEnd=round(3*swingTimeCount/(1-dutyCycle)*(2*dutyCycle-1)/4+2*swingTimeCount);
        if(iEnd-iStart>swingTimeCount)
        {
            printf("WARNING!!! s_w short than swingTimeCount!\n");
            iEnd=iStart+swingTimeCount;
        }
        for(int i=iStart;i<iEnd;i++)
        {
            pEB[0]=initPeb[0]+*(pva_b_t+9*i+2);
            pEB[2]=initPeb[2]+*(pva_b_t+9*i+1);
            pEB[3]=initPeb[3]+*(pva_b_t+9*i+0);
            vEB[0]=*(pva_b_t+9*i+5);
            vEB[2]=*(pva_b_t+9*i+4);
            vEB[3]=*(pva_b_t+9*i+3);
            aEB[0]=*(pva_b_t+9*i+8);
            aEB[2]=*(pva_b_t+9*i+7);
            aEB[3]=*(pva_b_t+9*i+6);
            for(int j=0;j<3;j++)
            {
                ratio_f[0]=initPee[6*j+3]*(cos(5*stepBeta/4)-cos(stepBeta/4))+initPee[6*j+5]*(sin(5*stepBeta/4)-sin(stepBeta/4));
                ratio_f[2]=-initPee[6*j+3]*(sin(5*stepBeta/4)-sin(stepBeta/4))+initPee[6*j+5]*(cos(5*stepBeta/4)-cos(stepBeta/4));

                f_sw[0]=-ratio_f[0]/2*(cos(*(s_w_t+6*(i-iStart)+2*j+1))-1);
                f_sw[1]=stepH*sin(*(s_w_t+6*(i-iStart)+2*j+1));
                f_sw[2]=-ratio_f[2]/2*(cos(*(s_w_t+6*(i-iStart)+2*j+1))-1);

                df_sw[0]=ratio_f[0]/2*sin(*(s_w_t+6*(i-iStart)+2*j+1));
                df_sw[1]=stepH*cos(*(s_w_t+6*(i-iStart)+2*j+1));
                df_sw[2]=ratio_f[2]/2*sin(*(s_w_t+6*(i-iStart)+2*j+1));

                ddf_sw[0]=ratio_f[0]/2*cos(*(s_w_t+6*(i-iStart)+2*j+1));
                ddf_sw[1]=-stepH*sin(*(s_w_t+6*(i-iStart)+2*j+1));
                ddf_sw[2]=ratio_f[2]/2*cos(*(s_w_t+6*(i-iStart)+2*j+1));

                pEE[6*j]=initPee[6*j]*cos(3*stepBeta/4)+initPee[6*j+2]*sin(3*stepBeta/4);
                pEE[6*j+1]=initPee[6*j+1];
                pEE[6*j+2]=-initPee[6*j]*sin(3*stepBeta/4)+initPee[6*j+2]*cos(3*stepBeta/4);

                pEE[6*j+3]=initPee[6*j+3]*cos(stepBeta/4)+initPee[6*j+5]*sin(stepBeta/4)+f_sw[0];
                pEE[6*j+4]=initPee[6*j+4]+f_sw[1];
                pEE[6*j+5]=-initPee[6*j+3]*sin(stepBeta/4)+initPee[6*j+5]*cos(stepBeta/4)+f_sw[2];

                vEE[6*j+3]=df_sw[0]*(*(ds_w_t+6*(i-iStart)+2*j+1));
                vEE[6*j+4]=df_sw[1]*(*(ds_w_t+6*(i-iStart)+2*j+1));
                vEE[6*j+5]=df_sw[2]*(*(ds_w_t+6*(i-iStart)+2*j+1));

                aEE[6*j+3]=ddf_sw[0]*(*(ds_w_t+6*(i-iStart)+2*j+1))*(*(ds_w_t+6*(i-iStart)+2*j+1))+df_sw[0]*(*(dds_w_t+6*(i-iStart)+2*j+1));
                aEE[6*j+4]=ddf_sw[1]*(*(ds_w_t+6*(i-iStart)+2*j+1))*(*(ds_w_t+6*(i-iStart)+2*j+1))+df_sw[1]*(*(dds_w_t+6*(i-iStart)+2*j+1));
                aEE[6*j+5]=ddf_sw[2]*(*(ds_w_t+6*(i-iStart)+2*j+1))*(*(ds_w_t+6*(i-iStart)+2*j+1))+df_sw[2]*(*(dds_w_t+6*(i-iStart)+2*j+1));
            }
            rbt.SetPeb(pEB,"213");
            rbt.SetVb(vEB);
            rbt.SetAb(aEB);
            rbt.SetPee(pEE);
            rbt.SetVee(vEE);
            rbt.SetAee(aEE);

            for(int j=0;j<6;j++)
            {
                *(VeeMinus_t+18*i+3*j)=vEE[3*j];
                *(VeeMinus_t+18*i+3*j+1)=vEE[3*j+1];
                *(VeeMinus_t+18*i+3*j+2)=vEE[3*j+2]-vEB[2];
                *(AeeMinus_t+18*i+3*j)=aEE[3*j];
                *(AeeMinus_t+18*i+3*j+1)=aEE[3*j+1];
                *(AeeMinus_t+18*i+3*j+2)=aEE[3*j+2]-aEB[2];
            }
            rbt.GetVee(VeeB_t+18*i,rbt.body());
            rbt.GetAee(AeeB_t+18*i,rbt.body());

            rbt.GetPee(Pee_t+18*i,rbt.body());
            rbt.GetPin(Pin_t+18*i);
            rbt.GetVin(Vin_t+18*i);
            rbt.GetAin(Ain_t+18*i);
        }

        std::fill_n(vEE,18,0);
        std::fill_n(aEE,18,0);
        iStart=iEnd;
        iEnd=stepTimeCount;
        for(int i=iStart;i<iEnd;i++)
        {
            pEB[0]=initPeb[0]+*(pva_b_t+9*i+2);
            pEB[2]=initPeb[2]+*(pva_b_t+9*i+1);
            pEB[3]=initPeb[3]+*(pva_b_t+9*i+0);
            vEB[0]=*(pva_b_t+9*i+5);
            vEB[2]=*(pva_b_t+9*i+4);
            vEB[3]=*(pva_b_t+9*i+3);
            aEB[0]=*(pva_b_t+9*i+8);
            aEB[2]=*(pva_b_t+9*i+7);
            aEB[3]=*(pva_b_t+9*i+6);
            for(int j=0;j<3;j++)
            {
                pEE[6*j]=initPee[6*j]*cos(3*stepBeta/4)+initPee[6*j+2]*sin(3*stepBeta/4);
                pEE[6*j+1]=initPee[6*j+1];
                pEE[6*j+2]=-initPee[6*j]*sin(3*stepBeta/4)+initPee[6*j+2]*cos(3*stepBeta/4);

                pEE[6*j+3]=initPee[6*j+3]*cos(5*stepBeta/4)+initPee[6*j+5]*sin(5*stepBeta/4);
                pEE[6*j+4]=initPee[6*j+4];
                pEE[6*j+5]=-initPee[6*j+3]*sin(5*stepBeta/4)+initPee[6*j+5]*cos(5*stepBeta/4);
            }
            rbt.SetPeb(pEB,"213");
            rbt.SetVb(vEB);
            rbt.SetAb(aEB);
            rbt.SetPee(pEE);
            rbt.SetVee(vEE);
            rbt.SetAee(aEE);

            for(int j=0;j<6;j++)
            {
                *(VeeMinus_t+18*i+3*j)=vEE[3*j];
                *(VeeMinus_t+18*i+3*j+1)=vEE[3*j+1];
                *(VeeMinus_t+18*i+3*j+2)=vEE[3*j+2]-vEB[2];
                *(AeeMinus_t+18*i+3*j)=aEE[3*j];
                *(AeeMinus_t+18*i+3*j+1)=aEE[3*j+1];
                *(AeeMinus_t+18*i+3*j+2)=aEE[3*j+2]-aEB[2];
            }
            rbt.GetVee(VeeB_t+18*i,rbt.body());
            rbt.GetAee(AeeB_t+18*i,rbt.body());

            rbt.GetPee(Pee_t+18*i,rbt.body());
            rbt.GetPin(Pin_t+18*i);
            rbt.GetVin(Vin_t+18*i);
            rbt.GetAin(Ain_t+18*i);
        }

        memcpy(out_tippos,Pee_t,18*stepTimeCount*sizeof(double));
        out_bodyvel=*(pva_b_t+1);
        out_period=0.001*stepTimeCount;

        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/pva_b_t.txt",pva_b_t,stepTimeCount,9);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/s_w_t.txt",s_w_t,swingTimeCount,6);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Pee_t.txt",Pee_t,stepTimeCount,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Pin_t.txt",Pin_t,stepTimeCount,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Vin_t.txt",Vin_t,stepTimeCount,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Ain_t.txt",Ain_t,stepTimeCount,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/VeeB_t.txt",VeeB_t,stepTimeCount,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/VeeMinus_t.txt",VeeMinus_t,stepTimeCount,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/AeeB_t.txt",AeeB_t,stepTimeCount,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/AeeMinus_t.txt",AeeMinus_t,stepTimeCount,18);
        delete [] pva_b_t;
        delete [] s_w_t;
        delete [] Pee_t;
        delete [] Pin_t;
        delete [] Vin_t;
        delete [] Ain_t;
        delete [] AeeB_t;
        delete [] AeeMinus_t;
        delete [] VeeB_t;
        delete [] VeeMinus_t;
        printf("finish GetOptimalGait2t\n");
    }

    void NonRTOptimalGCSRotate::OutputData()
    {
        printf("Start output data...\n");
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_upBound_aLmt_body.txt",ds_upBound_aLmt_body,2201,1);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_upBound_vLmt_body.txt",ds_upBound_vLmt_body,2201,1);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_upBound_body.txt",dds_upBound_body,2201,1);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_upBound_body.txt",dds_lowBound_body,2201,1);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_forward_body.txt",ds_forward_body,2201,1);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_backward_body.txt",ds_backward_body,2201,1);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_forward_body.txt",dds_forward_body,2201,1);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_backward_body.txt",dds_backward_body,2201,1);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_ds_body.txt",real_ds_body,2201,1);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_dds_body.txt",real_dds_body,2201,1);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_ddsMax_body.txt",real_ddsMax_body,2201,1);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_ddsMin_body.txt",real_ddsMin_body,2201,1);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/timeArray_body_tmp.txt",timeArray_body_tmp,2201,1);

        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_C.txt",*output_C,2201,36);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_PebPee.txt",*output_PebPee,2201,24);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_PeeB.txt",*output_PeeB,2201,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_dsds.txt", *output_dsds, 2201,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_ds.txt", *output_ds, 2201,18);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_ds1.txt",*output_ds1,2201,18);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_ds2.txt",*output_ds2,2201,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_const.txt", *output_const, 2201,18);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_const1.txt", *output_const1, 2201,18);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_const2.txt", *output_const2, 2201,18);
//        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_const3.txt", *output_const3, 2201,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_dds.txt",*output_dds,2201,18);

        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a2.txt",*output_a2,2201,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a1.txt",*output_a1,2201,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a0L.txt",*output_a0L,2201,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/param_a0H.txt",*output_a0H,2201,18);

        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/max_ValueL.txt",*dds_lowBound,901,6);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/min_ValueH.txt",*dds_upBound,901,6);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_upBound_aLmt.txt",*ds_upBound_aLmt,901,6);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_lowBound_aLmt.txt",*ds_lowBound_aLmt,901,6);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_upBound_vLmt.txt",*ds_upBound_vLmt,901,6);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_lowBound.txt",*ds_lowBound,901,6);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_upBound.txt",*ds_upBound,901,6);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_forward.txt",*ds_forward,901,6);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/ds_backward.txt",*ds_backward,901,6);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_forward.txt",*dds_forward,901,6);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/dds_backward.txt",*dds_backward,901,6);


        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_ds.txt",*real_ds,901,6);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_dds.txt",*real_dds,901,6);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_ddsMax.txt",*real_ddsMax,901,6);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/real_ddsMin.txt",*real_ddsMin,901,6);

        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/s_w.txt",*s_w,901,6);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/pva_b.txt",*pva_b,2201,9);

        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/timeArray_tmp.txt",*timeArray_tmp,901,6);

        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Pee_s.txt",*Pee_s,2201,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Pin_s.txt",*Pin_s,2201,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Vin_s.txt",*Vin_s,2201,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Ain_s.txt",*Ain_s,2201,18);
        printf("Finish output data\n");
    }

	/*利用以上求解拼接出完整的匀速段轨迹，并与未优化轨迹对比*/
    void NonRTOptimalGCSRotate::GetNormalGait(double &out_period)
    {
        printf("Start GetNormalGait\n");
        double pEB[6] {0};
        double pEE[18] {0};
        int TotalCount=stepTimeCount;
        int countLen {100};
        double maxAin {0};
        double minAin {0};
        double maxAbsAin {0};
        double maxVin {0};
        double minVin {0};
        double maxAbsVin {0};
        int k {0};
        bool stopFlag {false};
        double ratio_f[2] {0};

        while(stopFlag==false && k<=100)
        {
            double * normalPee = new double [(TotalCount+1)*18];
            double * normalPin = new double [(TotalCount+1)*18];
            double * normalVin = new double [TotalCount*18];
            double * normalAin = new double [(TotalCount-1)*18];

            int iStart;
            int iEnd;
            iEnd=round(TotalCount*(2*dutyCycle-1)/4);
            for(int i=0;i<iEnd;i++)
            {
                pEB[0]=initPeb[0];
                pEB[2]=initPeb[2];
                pEB[3]=initPeb[3]+stepBeta*i/TotalCount;
                for(int j=0;j<3;j++)
                {
                    pEE[6*j]=initPee[6*j]*cos(-stepBeta/4)+initPee[6*j+2]*sin(-stepBeta/4);
                    pEE[6*j+1]=initPee[6*j+1];
                    pEE[6*j+2]=-initPee[6*j]*sin(-stepBeta/4)+initPee[6*j+2]*cos(-stepBeta/4);

                    pEE[6*j+3]=initPee[6*j+3]*cos(stepBeta/4)+initPee[6*j+5]*sin(stepBeta/4);
                    pEE[6*j+4]=initPee[6*j+4];
                    pEE[6*j+5]=-initPee[6*j+3]*sin(stepBeta/4)+initPee[6*j+5]*cos(stepBeta/4);
                }
                rbt.SetPeb(pEB,"213");
                rbt.SetPee(pEE);

                rbt.GetPee(normalPee+18*i,rbt.body());
                rbt.GetPin(normalPin+18*i);
            }

            iStart=iEnd;
            iEnd=round(TotalCount*((2*dutyCycle-1)/4+(1-dutyCycle)));
            for (int i=iStart;i<iEnd;i++)
            {
                pEB[0]=initPeb[0];
                pEB[2]=initPeb[2];
                pEB[3]=initPeb[3]+stepBeta*i/TotalCount;
                for(int j=0;j<3;j++)
                {
                    ratio_f[0]=initPee[6*j]*(cos(3*stepBeta/4)-cos(-stepBeta/4))+initPee[6*j+2]*(sin(3*stepBeta/4)-sin(-stepBeta/4));
                    ratio_f[2]=-initPee[6*j]*(sin(3*stepBeta/4)-sin(-stepBeta/4))+initPee[6*j+2]*(cos(3*stepBeta/4)-cos(-stepBeta/4));

                    pEE[6*j]=initPee[6*j]*cos(-stepBeta/4)+initPee[6*j+2]*sin(-stepBeta/4)-ratio_f[0]/2*(cos(PI*(i-iStart)/(iEnd-iStart))-1);
                    pEE[6*j+1]=initPee[6*j+1]-stepH/2*(cos(2*PI*(i-iStart)/(iEnd-iStart))-1);
                    pEE[6*j+2]=-initPee[6*j]*sin(-stepBeta/4)+initPee[6*j+2]*cos(-stepBeta/4)+-ratio_f[2]/2*(cos(PI*(i-iStart)/(iEnd-iStart))-1);

                    pEE[6*j+3]=initPee[6*j+3]*cos(stepBeta/4)+initPee[6*j+5]*sin(stepBeta/4);
                    pEE[6*j+4]=initPee[6*j+4];
                    pEE[6*j+5]=-initPee[6*j+3]*sin(stepBeta/4)+initPee[6*j+5]*cos(stepBeta/4);
                }
                rbt.SetPeb(pEB,"213");
                rbt.SetPee(pEE);

                rbt.GetPee(normalPee+18*i,rbt.body());
                rbt.GetPin(normalPin+18*i);
            }

            iStart=iEnd;
            iEnd=round(TotalCount*(3*(2*dutyCycle-1)/4+(1-dutyCycle)));
            for(int i=iStart;i<iEnd;i++)
            {
                pEB[0]=initPeb[0];
                pEB[2]=initPeb[2];
                pEB[3]=initPeb[3]+stepBeta*i/TotalCount;
                for(int j=0;j<3;j++)
                {
                    pEE[6*j]=initPee[6*j]*cos(3*stepBeta/4)+initPee[6*j+2]*sin(3*stepBeta/4);
                    pEE[6*j+1]=initPee[6*j+1];
                    pEE[6*j+2]=-initPee[6*j]*sin(3*stepBeta/4)+initPee[6*j+2]*cos(3*stepBeta/4);

                    pEE[6*j+3]=initPee[6*j+3]*cos(stepBeta/4)+initPee[6*j+5]*sin(stepBeta/4);
                    pEE[6*j+4]=initPee[6*j+4];
                    pEE[6*j+5]=-initPee[6*j+3]*sin(stepBeta/4)+initPee[6*j+5]*cos(stepBeta/4);
                }
                rbt.SetPeb(pEB,"213");
                rbt.SetPee(pEE);

                rbt.GetPee(normalPee+18*i,rbt.body());
                rbt.GetPin(normalPin+18*i);
            }

            iStart=iEnd;
            iEnd=round(TotalCount*(3*(2*dutyCycle-1)/4+2*(1-dutyCycle)));
            for(int i=iStart;i<iEnd;i++)
            {
                pEB[0]=initPeb[0];
                pEB[2]=initPeb[2];
                pEB[3]=initPeb[3]+stepBeta*i/TotalCount;
                for(int j=0;j<3;j++)
                {
                    ratio_f[0]=initPee[6*j+3]*(cos(5*stepBeta/4)-cos(stepBeta/4))+initPee[6*j+5]*(sin(5*stepBeta/4)-sin(stepBeta/4));
                    ratio_f[2]=-initPee[6*j+3]*(sin(5*stepBeta/4)-sin(stepBeta/4))+initPee[6*j+5]*(cos(5*stepBeta/4)-cos(stepBeta/4));

                    pEE[6*j]=initPee[6*j]*cos(3*stepBeta/4)+initPee[6*j+2]*sin(3*stepBeta/4);
                    pEE[6*j+1]=initPee[6*j+1];
                    pEE[6*j+2]=-initPee[6*j]*sin(3*stepBeta/4)+initPee[6*j+2]*cos(3*stepBeta/4);

                    pEE[6*j+3]=initPee[6*j+3]*cos(stepBeta/4)+initPee[6*j+5]*sin(stepBeta/4)-ratio_f[0]/2*(cos(PI*(i-iStart)/(iEnd-iStart))-1);
                    pEE[6*j+4]=initPee[6*j+4]-stepH/2*(cos(2*PI*(i-iStart)/(iEnd-iStart))-1);
                    pEE[6*j+5]=-initPee[6*j+3]*sin(stepBeta/4)+initPee[6*j+5]*cos(stepBeta/4)-ratio_f[2]/2*(cos(PI*(i-iStart)/(iEnd-iStart))-1);
                }
                rbt.SetPeb(pEB,"213");
                rbt.SetPee(pEE);

                rbt.GetPee(normalPee+18*i,rbt.body());
                rbt.GetPin(normalPin+18*i);
            }

            iStart=iEnd;
            iEnd=TotalCount+1;
            for(int i=iStart;i<iEnd;i++)
            {
                pEB[0]=initPeb[0];
                pEB[2]=initPeb[2];
                pEB[3]=initPeb[3]+stepBeta*i/TotalCount;
                for(int j=0;j<3;j++)
                {
                    pEE[6*j]=initPee[6*j]*cos(3*stepBeta/4)+initPee[6*j+2]*sin(3*stepBeta/4);
                    pEE[6*j+1]=initPee[6*j+1];
                    pEE[6*j+2]=-initPee[6*j]*sin(3*stepBeta/4)+initPee[6*j+2]*cos(3*stepBeta/4);

                    pEE[6*j+3]=initPee[6*j+3]*cos(5*stepBeta/4)+initPee[6*j+5]*sin(5*stepBeta/4);
                    pEE[6*j+4]=initPee[6*j+4];
                    pEE[6*j+5]=-initPee[6*j+3]*sin(5*stepBeta/4)+initPee[6*j+5]*cos(5*stepBeta/4);
                }
                rbt.SetPeb(pEB,"213");
                rbt.SetPee(pEE);

                rbt.GetPee(normalPee+18*i,rbt.body());
                rbt.GetPin(normalPin+18*i);
            }

            for(int i=0;i<TotalCount;i++)
            {
                for(int j=0;j<18;j++)
                {
                    *(normalVin+18*i+j)=(*(normalPin+18*i+18+j)-*(normalPin+18*i+j))*1000;
                }
            }
            for(int i=0;i<TotalCount-1;i++)
            {
                for(int j=0;j<18;j++)
                {
                    *(normalAin+18*i+j)=(*(normalVin+18*i+18+j)-*(normalVin+18*i+j))*1000;
                }
            }

            maxVin=*std::max_element(normalVin,normalVin+18*TotalCount);
            minVin=*std::min_element(normalVin,normalVin+18*TotalCount);
            maxAbsVin=std::max(fabs(maxVin),fabs(minVin));

            maxAin=*std::max_element(normalAin,normalAin+18*(TotalCount-1));
            minAin=*std::min_element(normalAin,normalAin+18*(TotalCount-1));
            maxAbsAin=std::max(fabs(maxAin),fabs(minAin));

            if(maxAbsAin<aLmt && maxAbsVin<vLmt)
            {
                if(countLen==1)
                {
                    if(k==0)
                    {
                        printf("How impossible! NormalGait is faster than OptimalGait!\n");
                    }
                    else
                    {
                        out_period=TotalCount*0.001;
                        stopFlag=true;
                        printf("Ain=%.4f or Vin=%.4f reach the maximum(A:%.4f,V:%.4f) at TotalCount=%d(from %d), iteration time:%d\n",maxAin,maxVin,aLmt,vLmt,TotalCount,stepTimeCount,k);
                        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/normalPee.txt",normalPee,TotalCount+1,18);
                        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/normalPin.txt",normalPin,TotalCount+1,18);
                        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/normalVin.txt",normalVin,TotalCount,18);
                        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/normalAin.txt",normalAin,TotalCount-1,18);
                    }
                }
                else
                {
                    TotalCount-=countLen;
                    countLen/=10;
                    //printf("countLen change to %d\n",countLen);
                }
            }
            else
            {
                TotalCount+=countLen;
                k++;
            }

            if(k==101)
            {
                //output data to analyse when error
                aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/normalPee.txt",normalPee,TotalCount+1,18);
                aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/normalPin.txt",normalPin,TotalCount+1,18);
                aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/normalVin.txt",normalVin,TotalCount,18);
                aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/normalAin.txt",normalAin,TotalCount-1,18);
            }
            delete [] normalPee;
            delete [] normalPin;
            delete [] normalVin;
            delete [] normalAin;
        }


        printf("Finish GetNormalGait,k=%d,Ain=%.4f or Vin=%.4f, TotalCount=%d\n",k,maxAbsAin,maxAbsVin,TotalCount);
    }

	/*包括加减速段*/
    /******only for alpha=0 and beta=0**********/
    void NonRTOptimalGCSRotate::GetEntireGait()
    {
        printf("Start GetEntireGait\n");
        double c3;
        double c2;
        double c1;
        double c0;
        double pEB[6];
        double pEE[18];
        memcpy(pEB,initPeb,6*sizeof(double));

        double * entirePeb=new double [6*2*stepTimeCount];
        double * entirePee=new double [18*2*stepTimeCount];
        double * entirePin=new double [18*2*stepTimeCount];
        double * Pin_const=new double [18*stepTimeCount];
        double * pva_b_const=new double [3*stepTimeCount];

        //acc for stepTimeCount/2
        c3=(pva_b[0][1]+stepD/2/(0.5*stepTimeCount*0.001))/(0.5*stepTimeCount*0.001)/(0.5*stepTimeCount*0.001);
        c2=(pva_b[0][1]-3*(pva_b[0][1]+stepD/2/(0.5*stepTimeCount*0.001)))/(stepTimeCount*0.001);
        int iStart;
        int iEnd;
        iEnd=(int)(stepTimeCount*(2*dutyCycle-1)/4)+1;
        for(int i=0;i<iEnd;i++)
        {
            pEB[2]=c3*1e-9*i*i*i+c2*1e-6*i*i;
            memcpy(pEE,initPee,18*sizeof(double));

            rbt.SetPeb(pEB);
            rbt.SetPee(pEE);

            memcpy(entirePeb+6*i,pEB,6*sizeof(double));
            memcpy(entirePee+18*i,pEE,18*sizeof(double));
            rbt.GetPin(entirePin+18*i);
        }
        iStart=iEnd;
        iEnd=(int)(stepTimeCount*((2*dutyCycle-1)/4+(1-dutyCycle)))+1;
        for(int i=iStart;i<iEnd;i++)
        {
            pEB[2]=c3*1e-9*i*i*i+c2*1e-6*i*i;
            memcpy(pEE,initPee,18*sizeof(double));
            for(int j=0;j<3;j++)
            {
                pEE[6*j+4]=initPee[6*j+4]-stepH/2*(cos(2*PI*(i-iStart)/(iEnd-iStart))-1);
                pEE[6*j+5]=initPee[6*j+5]+stepD/4*(cos(PI*(i-iStart)/(iEnd-iStart))-1);
            }

            rbt.SetPeb(pEB);
            rbt.SetPee(pEE);

            memcpy(entirePeb+6*i,pEB,6*sizeof(double));
            memcpy(entirePee+18*i,pEE,18*sizeof(double));
            rbt.GetPin(entirePin+18*i);
        }
        iStart=iEnd;
        iEnd=stepTimeCount/2;
        for(int i=iStart;i<iEnd;i++)
        {
            pEB[2]=c3*1e-9*i*i*i+c2*1e-6*i*i;
            memcpy(pEE,initPee,18*sizeof(double));
            for(int j=0;j<3;j++)
            {
                pEE[6*j+5]=initPee[6*j+5]-stepD/2;
            }

            rbt.SetPeb(pEB);
            rbt.SetPee(pEE);

            memcpy(entirePeb+6*i,pEB,6*sizeof(double));
            memcpy(entirePee+18*i,pEE,18*sizeof(double));
            rbt.GetPin(entirePin+18*i);
        }

        //const for stepTimeCount
        aris::dynamic::dlmread("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/Pin_t.txt",Pin_const);
        aris::dynamic::dlmread("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/pva_b_t.txt",pva_b_const);
        for(int i=0;i<stepTimeCount;i++)
        {
            pEB[2]=-stepD/4+*(pva_b_const+3*i);
            rbt.SetPeb(pEB);
            rbt.SetPin(Pin_const+18*i);

            memcpy(entirePeb+6*(i+stepTimeCount/2),pEB,6*sizeof(double));
            rbt.GetPee(entirePee+18*(i+stepTimeCount/2));
            memcpy(entirePin+18*(i+stepTimeCount/2),Pin_const+18*i,18*sizeof(double));
        }

        //dec for stepTimeCount/2
        c0=-stepD-stepD/4;
        c1=pva_b[count5][1];
        c2=(-3*stepD/4-2*c1*(0.5*stepTimeCount*0.001))/(0.5*stepTimeCount*0.001)/(0.5*stepTimeCount*0.001);
        c3=(stepD/2+c1*(0.5*stepTimeCount*0.001))/(0.5*stepTimeCount*0.001)/(0.5*stepTimeCount*0.001)/(0.5*stepTimeCount*0.001);

        iEnd=(int)(stepTimeCount*(2*dutyCycle-1)/4)+1;
        for(int i=0;i<iEnd;i++)
        {
            pEB[2]=c3*1e-9*i*i*i+c2*1e-6*i*i+c1*1e-3*i+c0;
            memcpy(pEE,initPee,18*sizeof(double));
            for(int j=0;j<3;j++)
            {
                pEE[6*j+2]=initPee[6*j+2]-stepD;
                pEE[6*j+5]=initPee[6*j+5]-stepD/2-stepD;
            }

            rbt.SetPeb(pEB);
            rbt.SetPee(pEE);

            memcpy(entirePeb+6*(i+3*stepTimeCount/2),pEB,6*sizeof(double));
            memcpy(entirePee+18*(i+3*stepTimeCount/2),pEE,18*sizeof(double));
            rbt.GetPin(entirePin+18*(i+3*stepTimeCount/2));
        }
        iStart=iEnd;
        iEnd=(int)(stepTimeCount*((2*dutyCycle-1)/4+(1-dutyCycle)))+1;
        for(int i=iStart;i<iEnd;i++)
        {
            pEB[2]=c3*1e-9*i*i*i+c2*1e-6*i*i+c1*1e-3*i+c0;
            memcpy(pEE,initPee,18*sizeof(double));
            for(int j=0;j<3;j++)
            {
                pEE[6*j+1]=initPee[6*j+1]-stepH/2*(cos(2*PI*(i-iStart)/(iEnd-iStart))-1);
                pEE[6*j+2]=initPee[6*j+2]-stepD+stepD/4*(cos(PI*(i-iStart)/(iEnd-iStart))-1);

                pEE[6*j+5]=initPee[6*j+5]-stepD/2-stepD;
            }
            rbt.SetPeb(pEB);
            rbt.SetPee(pEE);

            memcpy(entirePeb+6*(i+3*stepTimeCount/2),pEB,6*sizeof(double));
            memcpy(entirePee+18*(i+3*stepTimeCount/2),pEE,18*sizeof(double));
            rbt.GetPin(entirePin+18*(i+3*stepTimeCount/2));
        }
        iStart=iEnd;
        iEnd=stepTimeCount/2;
        for(int i=iStart;i<iEnd;i++)
        {
            pEB[2]=c3*1e-9*i*i*i+c2*1e-6*i*i+c1*1e-3*i+c0;
            memcpy(pEE,initPee,18*sizeof(double));
            for(int j=0;j<3;j++)
            {
                pEE[6*j+2]=initPee[6*j+2]-stepD/2-stepD;
                pEE[6*j+5]=initPee[6*j+5]-stepD/2-stepD;
            }

            rbt.SetPeb(pEB);
            rbt.SetPee(pEE);

            memcpy(entirePeb+6*(i+3*stepTimeCount/2),pEB,6*sizeof(double));
            memcpy(entirePee+18*(i+3*stepTimeCount/2),pEE,18*sizeof(double));
            rbt.GetPin(entirePin+18*(i+3*stepTimeCount/2));
        }

        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/entirePeb.txt",entirePeb,2*stepTimeCount,6);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/entirePee.txt",entirePee,2*stepTimeCount,18);
        aris::dynamic::dlmwrite("/home/hex/Desktop/mygit/RobotVIII_demo/build/bin/entirePin.txt",entirePin,2*stepTimeCount,18);

        delete [] entirePeb;
        delete [] entirePee;
        delete [] entirePin;
        delete [] Pin_const;
        delete [] pva_b_const;
        printf("Finish GetEntireGait\n");
    }


    /*
    <fwg default="fwg_param">
        <fwg_param type="group">
            <totalCount abbreviation="t" type="int" default="2000"/>
            <n abbreviation="n" type="int" default="2"/>
            <distance abbreviation="d" type="double" default="0.5"/>
            <height abbreviation="h" type="double" default="0.05"/>
        </fwg_param>
    </fwg>

    rs.addCmd("fwg",TimeOptimal::FastWalkGCS::parseFastWalk,TimeOptimal::FastWalkGCS::fastWalk);
    */
/*
    double FastWalkGCS::initBodyPos[6];
    double FastWalkGCS::initPee[18];
    double FastWalkGCS::swingPee_scale[3001][18];
    double FastWalkGCS::bodyConstVel;
    double FastWalkGCS::dutyCycle;
    FastWalkGCSParam FastWalkGCS::fwgParam;
    Pipe<FastWalkGCSParam> FastWalkGCS::fwgPipe(true);
    std::thread FastWalkGCS::fwgThread;

    void FastWalkGCS::parseFastWalk(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
    {
        Robots::WalkParam param;

        for (auto i : params)
        {
            if (i.first == "distance")
            {
                param.d = std::stod(i.second);
            }
            else if (i.first == "height")
            {
                param.h = std::stod(i.second);
            }
            else if (i.first == "n")
            {
                param.n = std::stoi(i.second);
            }
            else if (i.first == "totalCount")
            {
                param.totalCount = std::stoi(i.second);
            }
        }

        std::fill_n(initBodyPos,6,0);
//        double initPee[18] { -0.3, -0.85, -0.65,
//                         -0.45, -0.85,     0,
//                          -0.3, -0.85,  0.65,
//                           0.3, -0.85, -0.65,
//                          0.45, -0.85,     0,
//                           0.3, -0.85,  0.65 };
        double initTipPos[18] { -0.30, -0.58, -0.52,
                             -0.60, -0.58,  0,
                             -0.30, -0.58,  0.52,
                              0.30, -0.58, -0.52,
                              0.60, -0.58,  0,
                              0.30, -0.58,  0.52 };
        memcpy(initPee,initTipPos,18*sizeof(double));

        GetScalingPath(param);
        msg.copyStruct(param);

    }

    int  FastWalkGCS::fastWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
    {
        auto &robot = static_cast<Robots::RobotBase &>(model);
        auto &param = static_cast<const Robots::WalkParam &>(param_in);

        double pEB[6];
        double pEE[18];
        memcpy(pEB,initBodyPos,6*sizeof(double));
        memcpy(pEE,initPee,18*sizeof(double));

        double c3;
        double c2;
        double c1;
        double c0;
        int period_count=param.count%(param.totalCount/2);
        int iStart;
        int iEnd;
        iStart=(int)(param.totalCount*(2*dutyCycle-1)/4)+1;
        iEnd=(int)(param.totalCount*((2*dutyCycle-1)/4+(1-dutyCycle)))+1;
 
        if(param.count/(param.totalCount/2)==0)//acc 024 swing
        {
            //acc for param.totalCount/2
            c3=(bodyConstVel+param.d/2/(0.5*param.totalCount*0.001))/(0.5*param.totalCount*0.001)/(0.5*param.totalCount*0.001);
            c2=(bodyConstVel-3*(bodyConstVel+param.d/2/(0.5*param.totalCount*0.001)))/(param.totalCount*0.001);
            pEB[2]=c3*1e-9*period_count*period_count*period_count+c2*1e-6*period_count*period_count;

            if(period_count>=0 && period_count<(int)(param.totalCount*(2*dutyCycle-1)/4)+1)
            {
                pEB[2]=c3*1e-9*period_count*period_count*period_count+c2*1e-6*period_count*period_count;

                robot.SetPeb(pEB);
                robot.SetPee(pEE);
            }
            else if(period_count<(int)(param.totalCount*((2*dutyCycle-1)/4+(1-dutyCycle)))+1)
            {
               pEB[2]=c3*1e-9*period_count*period_count*period_count+c2*1e-6*period_count*period_count;
                for(int j=0;j<3;j++)
                {
                    pEE[6*j+4]=initPee[6*j+4]-param.h/2*(cos(2*PI*(period_count-iStart)/(iEnd-iStart))-1);
                    pEE[6*j+5]=initPee[6*j+5]+param.d/4*(cos(PI*(period_count-iStart)/(iEnd-iStart))-1);
                }

                robot.SetPeb(pEB);
                robot.SetPee(pEE);
            }
            else
            {
                pEB[2]=c3*1e-9*period_count*period_count*period_count+c2*1e-6*period_count*period_count;
                for(int j=0;j<3;j++)
                {
                    pEE[6*j+5]=initPee[6*j+5]-param.d/2;
                }

                robot.SetPeb(pEB);
                robot.SetPee(pEE);
            }
        }
        else if(param.count/(param.totalCount/2)==2*param.n-1)//dec 135 swing
        {
            c0=-param.d-param.d/4;
            c1=bodyConstVel;
            c2=(-3*param.d/4-2*c1*(0.5*param.totalCount*0.001))/(0.5*param.totalCount*0.001)/(0.5*param.totalCount*0.001);
            c3=(param.d/2+c1*(0.5*param.totalCount*0.001))/(0.5*param.totalCount*0.001)/(0.5*param.totalCount*0.001)/(0.5*param.totalCount*0.001);

            if(period_count>=0 && period_count<(int)(param.totalCount*(2*dutyCycle-1)/4)+1)
            {
                pEB[2]=c3*1e-9*period_count*period_count*period_count+c2*1e-6*period_count*period_count+c1*1e-3*period_count+c0;
                for(int j=0;j<3;j++)
                {
                    pEE[6*j+2]=initPee[6*j+2]-param.d;
                    pEE[6*j+5]=initPee[6*j+5]-param.d/2-param.d;
                }

                robot.SetPeb(pEB);
                robot.SetPee(pEE);
            }
            else if(period_count<(int)(param.totalCount*((2*dutyCycle-1)/4+(1-dutyCycle)))+1)
            {
                pEB[2]=c3*1e-9*period_count*period_count*period_count+c2*1e-6*period_count*period_count+c1*1e-3*period_count+c0;
                for(int j=0;j<3;j++)
                {
                    pEE[6*j+1]=initPee[6*j+1]-param.h/2*(cos(2*PI*(period_count-iStart)/(iEnd-iStart))-1);
                    pEE[6*j+2]=initPee[6*j+2]-param.d+param.d/4*(cos(PI*(period_count-iStart)/(iEnd-iStart))-1);

                    pEE[6*j+5]=initPee[6*j+5]-param.d/2-param.d;
                }
                robot.SetPeb(pEB);
                robot.SetPee(pEE);
            }
            else
            {
                pEB[2]=c3*1e-9*period_count*period_count*period_count+c2*1e-6*period_count*period_count+c1*1e-3*period_count+c0;
                for(int j=0;j<3;j++)
                {
                    pEE[6*j+2]=initPee[6*j+2]-param.d/2-param.d;
                    pEE[6*j+5]=initPee[6*j+5]-param.d/2-param.d;
                }

                robot.SetPeb(pEB);
                robot.SetPee(pEE);
            }
        }
        else//const
        {
            period_count=(param.count-param.totalCount/2)%param.totalCount;
            pEB[2]=-param.d/4;

            robot.SetPeb(pEB);
            robot.SetPee(*swingPee_scale+18*period_count,robot.body());
        }
        //rt_printf("count=%d\n",param.count);

        fwgParam.Peb=pEB[2];
        robot.GetPee(fwgParam.Pee,robot.body());
        robot.GetPin(fwgParam.Pin);
        fwgPipe.sendToNrt(fwgParam);

        return param.n*param.totalCount-param.count-1;
    }

    void FastWalkGCS::GetScalingPath(Robots::WalkParam &param)
    {
        double aLmt {5.4};
        double vLmt {1.5};
        double out_TipPos[3000][18] {0};
        double out_period {0};

        dutyCycle=0.5;
        TimeOptimal::NonRTOptimalGCSRotate planner;
        planner.GetTimeOptimalGait(param.d,param.h,dutyCycle,aLmt,vLmt,initPee,*out_TipPos,bodyConstVel,out_period);
        printf("out_period is %.3f\n",out_period);

        int out_period_count=round(out_period*1000);
        int totalCount = param.totalCount;
        int k {0};
        double ratio = (double)totalCount/out_period_count;//>=1
        if(ratio<1)
        {
            printf("Warning! TotalCount(%d) smaller than the minimum count(%d)! We have made TotalCount=minimum count!\n",totalCount,out_period_count);
            //totalCount=out_period_count;
        }

        bodyConstVel/=ratio;
        for(int i=0; i<totalCount-1; i++)
        {
            for(int j=k; j<out_period_count; j++)
            {
                if(i >= ratio*j && i < ratio*(j+1))
                {
                    if(i==totalCount-2 || i==totalCount-1)
                    {
                        printf("j=%d,ratio*j=%.2f\n",j,ratio*j);
                    }
                    for(int leg_id=0;leg_id<6;leg_id++)
                    {
                        swingPee_scale[i][3*leg_id]=out_TipPos[j][3*leg_id]+(out_TipPos[j+1][3*leg_id]-out_TipPos[j][3*leg_id])*(i-ratio*j)/ratio;
                        swingPee_scale[i][3*leg_id+1]=out_TipPos[j][3*leg_id+1]+(out_TipPos[j+1][3*leg_id+1]-out_TipPos[j][3*leg_id+1])*(i-ratio*j)/ratio;
                        swingPee_scale[i][3*leg_id+2]=out_TipPos[j][3*leg_id+2]+(out_TipPos[j+1][3*leg_id+2]-out_TipPos[j][3*leg_id+2])*(i-ratio*j)/ratio;
                    }
                    k=j;
                    break;
                }
            }
        }
        for(int leg_id=0;leg_id<6;leg_id++)
        {
            int i=totalCount-1;
            int j=out_period_count-1;
            swingPee_scale[i][3*leg_id]=out_TipPos[j][3*leg_id]+(out_TipPos[0][3*leg_id]-out_TipPos[j][3*leg_id])*(i-ratio*j)/ratio;
            swingPee_scale[i][3*leg_id+1]=out_TipPos[j][3*leg_id+1]+(out_TipPos[0][3*leg_id+1]-out_TipPos[j][3*leg_id+1])*(i-ratio*j)/ratio;
            swingPee_scale[i][3*leg_id+2]=out_TipPos[j][3*leg_id+2]+(out_TipPos[0][3*leg_id+2]-out_TipPos[j][3*leg_id+2])*(i-ratio*j)/ratio;
        }
        aris::dynamic::dlmwrite("/var/log/aris/swingPee_scale.txt",*swingPee_scale,3000,18);
        aris::dynamic::dlmwrite("/var/log/aris/out_TipPos.txt",*out_TipPos,3000,18);

    }

    void FastWalkGCS::recordData()
    {
        fwgThread = std::thread([&]()
        {
            struct FastWalkGCSParam param;
            static std::fstream fileGait;
            std::string name = aris::core::logFileName();
            name.replace(name.rfind("log.txt"), std::strlen("FastWalkGCSData.txt"), "FastWalkGCSData.txt");
            fileGait.open(name.c_str(), std::ios::out | std::ios::trunc);

            long long count = -1;
            while (1)
            {
                fwgPipe.recvInNrt(param);

                fileGait << ++count << " ";
                fileGait << param.Peb << " ";
                for (int i=0;i<18;i++)
                {
                    fileGait << param.Pee[i] << "  ";
                }
                for (int i=0;i<18;i++)
                {
                    fileGait << param.Pin[i] << "  ";
                }
                fileGait << std::endl;
            }

            fileGait.close();
        });
    }
*/
}
