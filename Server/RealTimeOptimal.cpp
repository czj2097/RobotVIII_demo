#include "plan.h"
#include <cmath>
#include <cstdlib>
#include <algorithm>

void dlmwrite(const char *filename, const double *mtx, const int m, const int n)
{
    std::ofstream file;

    file.open(filename);

    file << std::setprecision(15);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            file << mtx[n*i + j] << "   ";
        }
        file << std::endl;
    }
}

RTOptimal::RTOptimal()
{
    rbt.loadXml("/usr/Robots/resource/Robot_Type_I/Robot_III.xml");
    rbt.SetPeb(initPeb);
    rbt.SetVb(initVeb);
    isCurPntPassed=true;
    lstPntTime=0;
}
RTOptimal::~RTOptimal(){}

/*find gait with maxVel in joint space by iteration*/

//pnts means the knots of the leg, which is a pntsNum*3 matrix
bool RTOptimal::GetTrajOneLeg(double *pnts, int pntsNum, double *startV, double *endV, int legID, int count, double *pIn, double *pEE)
{
    //GetParamInFromParamEE(pnts,pntsNum,startV,endV,legID);

    if(isCurPntPassed==true)
    {
        isCurPntPassed=false;
        lstPntTime+=nxtPntMinTime;

        for(int i=0;i<3;i++)
        {
            JudgeLineType(lineParam[3*legID+i]);
            switch(lineParam[3*legID+i].lineType)
            {
            case PolylineType::I:
                GetTwoTimeTypeI(lineParam[3*legID+i],curP2PParam[3*legID+i]);
                break;
            case PolylineType::II:
                GetTwoTimeTypeII(lineParam[3*legID+i],curP2PParam[3*legID+i]);
                break;
            case PolylineType::III:
                GetTwoTimeTypeIII(lineParam[3*legID+i],curP2PParam[3*legID+i]);
                break;
            default:
                break;
            }
        }
        DecideNxtPntTime(lineParam+3*legID);

        for(int i=0;i<3;i++)
        {
            GetNxtPntOptVel(lineParam[3*legID+i],nxtPntParam[3*legID+i]);
            GetNxtPntReachableVel(lineParam[3*legID+i],nxtPntParam[3*legID+i],i);

            RecordNxtPntTime(nxtPntParam[3*legID+i],count,i);
        }

        for(int i=0;i<3;i++)
        {
            printf("pnts of screw %d:",i);
            for(int k=0;k<lineParam[3*legID+i].pntsNum;k++)
            {
                printf("%.4f,",lineParam[3*legID+i].pnts[k]);
            }
            printf("\n");
        }
        printf("\n");
    }

    for(int i=0;i<3;i++)
    {
        pIn[3*legID+i]=GetNxtPin(lineParam[3*legID+i],nxtPntParam[3*legID+i],count);
    }
    rbt.pLegs[legID]->SetPin(pIn+3*legID);
    rbt.pLegs[legID]->GetPee(pEE+3*legID);

    if(lineParam[3*legID].pntsNum==1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void RTOptimal::GetParamInFromParamEE(double *pnts, int pntsNum, double *startV, double *endV, int legID)
{
    for(int i=0;i<3;i++)
    {
        lineParam[3*legID+i].pntsNum=pntsNum;
    }
    double startVin_tmp[3];
    double endVin_tmp[3];
    for(int i=0;i<pntsNum;i++)
    {
        double pIn_tmp[3];

        rbt.pLegs[legID]->SetPee(pnts+3*i,rbt.body());
        rbt.pLegs[legID]->GetPin(pIn_tmp);
        for(int k=0;k<3;k++)
        {
            lineParam[3*legID+k].pnts[i]=pIn_tmp[k];
        }

        if(i==0)
        {
            rbt.pLegs[legID]->SetVee(startV,rbt.body());
            rbt.pLegs[legID]->GetVin(startVin_tmp);
            for(int k=0;k<3;k++)
            {
                lineParam[3*legID+k].startV=startVin_tmp[k];
            }
        }
        if(i==pntsNum-1)
        {
            rbt.pLegs[legID]->SetVee(endV,rbt.body());
            rbt.pLegs[legID]->GetVin(endVin_tmp);
            for(int k=0;k<3;k++)
            {
                lineParam[3*legID+k].endV=endVin_tmp[k];
            }
        }
        printf("pIn of leg %d: %.4f, %.4f, %.4f\n",legID,lineParam[3*legID].pnts[i],lineParam[3*legID+1].pnts[i],lineParam[3*legID+2].pnts[i]);
    }
    printf("startVin: %.4f,%.4f,%.4f; endVin:%.4f,%.4f,%.4f\n",
           lineParam[3*legID].startV,lineParam[3*legID+1].startV,lineParam[3*legID+2].startV,
           lineParam[3*legID].endV,lineParam[3*legID+1].endV,lineParam[3*legID+2].endV);
    printf("startVin_tmp: %.4f,%.4f,%.4f; endVin_tmp:%.4f,%.4f,%.4f\n",
           startVin_tmp[0],startVin_tmp[1],startVin_tmp[2],
           endVin_tmp[0],endVin_tmp[1],endVin_tmp[2]);
    printf("startVee: %.4f,%.4f,%.4f; endVee:%.4f,%.4f,%.4f\n",
           startV[0],startV[1],startV[2],endV[0],endV[1],endV[2]);
}

void RTOptimal::JudgeLineType(PolylineParam &lnParam)
{
    lnParam.fstBrkPntNum=0;
    lnParam.lstBrkPntNum=lnParam.pntsNum-1;

    if(lnParam.pntsNum>2)
    {
        for(int i=1;i<lnParam.pntsNum-1;i++)
        {
            if((lnParam.pnts[i]-lnParam.pnts[i-1])*(lnParam.pnts[i]-lnParam.pnts[i+1])>0)
            {
                lnParam.fstBrkPntNum=i;
                break;
            }
        }
        for(int i=lnParam.pntsNum-2;i>0;i--)
        {
            if((lnParam.pnts[i]-lnParam.pnts[i-1])*(lnParam.pnts[i]-lnParam.pnts[i+1])>0)
            {
                lnParam.lstBrkPntNum=i;
                break;
            }
        }

        if(lnParam.fstBrkPntNum==0 && lnParam.lstBrkPntNum==lnParam.pntsNum-1)
        {
            lnParam.lineType=PolylineType::I;
        }
        else if(lnParam.fstBrkPntNum==lnParam.lstBrkPntNum)
        {
            lnParam.lineType=PolylineType::II;
        }
        else
        {
            lnParam.lineType=PolylineType::III;
        }
    }
    else if(lnParam.pntsNum==2)
    {
        lnParam.lineType=PolylineType::I;
    }
    else
    {
        printf("Error! Trajectory cannot be generated by only one knot!\n");
    }
}

void RTOptimal::GetTwoTimeTypeI(PolylineParam &lnParam, P2PMotionParam &p2pParam)
{
    lnParam.fstBrkPntNum=lnParam.pntsNum-1;
    lnParam.fstBrkPntVel=lnParam.endV;

    GetOptimalP2PMotionAcc(lnParam.pnts[0],lnParam.pnts[lnParam.fstBrkPntNum],lnParam.startV,lnParam.fstBrkPntVel,p2pParam);
    lnParam.totalTime=p2pParam.minTime;
    GetThreeNxtPntTime(lnParam,p2pParam);

    printf("Type I p2pParam time: %.4f,%.4f,%.4f, acc: %.4f,%.4f,%.4f, pos: %.4f,%.4f,%.4f\n",
           p2pParam.trajTime[0],p2pParam.trajTime[1],p2pParam.trajTime[2],
           p2pParam.trajAcc[0],p2pParam.trajAcc[1],p2pParam.trajAcc[2],
           p2pParam.trajPos[0],p2pParam.trajPos[1],p2pParam.trajPos[2]);
}

void RTOptimal::GetTwoTimeTypeII(PolylineParam &lnParam, P2PMotionParam &p2pParam)
{
    P2PMotionParam p2pParam_tmp;
    if(lnParam.pnts[lnParam.fstBrkPntNum]>lnParam.pnts[0])
    {
        double judgeSs=lnParam.pnts[0]+lnParam.startV*lnParam.startV/(2*aLmt);
        double judgeSe=lnParam.pnts[lnParam.pntsNum-1]+lnParam.endV*lnParam.endV/(2*aLmt);
        if(judgeSs<=lnParam.pnts[lnParam.fstBrkPntNum] && judgeSe<=lnParam.pnts[lnParam.fstBrkPntNum])
        {
            lnParam.fstBrkPntVel=lnParam.lstBrkPntVel=0;
        }
        else//When break point is passed twice, get the further time as the reach time
        {
            if(judgeSs>=judgeSe)
            {
                lnParam.lstBrkPntVel=lnParam.fstBrkPntVel
                        =GetBrkPntSndPassVel(judgeSs,lnParam.pnts[lnParam.pntsNum-1],0,lnParam.endV,lnParam.pnts[lnParam.fstBrkPntNum]);
            }
            else
            {
                lnParam.fstBrkPntVel=lnParam.lstBrkPntVel
                        =-sqrt(lnParam.endV*lnParam.endV-2*aLmt*(lnParam.pnts[lnParam.lstBrkPntNum]-lnParam.pnts[lnParam.pntsNum-1]));
            }
        }
    }
    else
    {
        double judgeSs=lnParam.pnts[0]-lnParam.startV*lnParam.startV/(2*aLmt);
        double judgeSe=lnParam.pnts[lnParam.pntsNum-1]-lnParam.endV*lnParam.endV/(2*aLmt);
        if(judgeSs>=lnParam.pnts[lnParam.fstBrkPntNum] && judgeSe>=lnParam.pnts[lnParam.fstBrkPntNum])
        {
            lnParam.fstBrkPntVel=lnParam.lstBrkPntVel=0;
        }
        else//When break point is passed twice, get the further time as the reach time
        {
            if(judgeSs<=judgeSe)
            {
                lnParam.lstBrkPntVel=lnParam.fstBrkPntVel
                        =GetBrkPntSndPassVel(judgeSs,lnParam.pnts[lnParam.pntsNum-1],0,lnParam.endV,lnParam.pnts[lnParam.fstBrkPntNum]);
            }
            else
            {
                lnParam.fstBrkPntVel=lnParam.lstBrkPntVel
                        =sqrt(lnParam.endV*lnParam.endV+2*aLmt*(lnParam.pnts[lnParam.lstBrkPntNum]-lnParam.pnts[lnParam.pntsNum-1]));
            }
        }
    }

    GetOptimalP2PMotionAcc(lnParam.pnts[0],lnParam.pnts[lnParam.fstBrkPntNum],lnParam.startV,lnParam.fstBrkPntVel,p2pParam);
    GetOptimalP2PMotionAcc(lnParam.pnts[lnParam.fstBrkPntNum],lnParam.pnts[lnParam.pntsNum-1],lnParam.fstBrkPntVel,lnParam.endV,p2pParam_tmp);
    lnParam.totalTime=p2pParam.minTime+p2pParam_tmp.minTime;
    GetThreeNxtPntTime(lnParam,p2pParam);

    printf("Type II p2pParam time: %.4f,%.4f,%.4f, acc: %.4f,%.4f,%.4f, pos: %.4f,%.4f,%.4f\n",
           p2pParam.trajTime[0],p2pParam.trajTime[1],p2pParam.trajTime[2],
           p2pParam.trajAcc[0],p2pParam.trajAcc[1],p2pParam.trajAcc[2],
           p2pParam.trajPos[0],p2pParam.trajPos[1],p2pParam.trajPos[2]);
}

void RTOptimal::GetTwoTimeTypeIII(PolylineParam &lnParam, P2PMotionParam &p2pParam)
{
    int sndBrkPntNum;
    for(int i=lnParam.fstBrkPntNum+1;i<=lnParam.lstBrkPntNum;i++)
    {
        if((lnParam.pnts[i]-lnParam.pnts[i-1])*(lnParam.pnts[i]-lnParam.pnts[i+1])>0)
        {
            sndBrkPntNum=i;
        }
    }

    P2PMotionParam p2pParam_tmp;
    if(lnParam.pnts[lnParam.fstBrkPntNum]>lnParam.pnts[0])
    {
        double judgeSs=lnParam.pnts[0]+lnParam.startV*lnParam.startV/(2*aLmt);
        if(judgeSs<=lnParam.pnts[lnParam.fstBrkPntNum])
        {
            lnParam.fstBrkPntVel=0;
        }
        else
        {
            if(sndBrkPntNum==lnParam.lstBrkPntNum)//pnts[lstBrkPntNum]<pnts[pntsNum-1]
            {
                double judgeSe=lnParam.pnts[lnParam.pntsNum-1]-lnParam.endV*lnParam.endV/(2*aLmt);
                if(judgeSe>=lnParam.pnts[lnParam.lstBrkPntNum])
                {
                    lnParam.lstBrkPntVel=0;
                }
                else
                {
                    lnParam.lstBrkPntVel=sqrt(lnParam.endV*lnParam.endV+2*aLmt*(lnParam.pnts[lnParam.lstBrkPntNum]-lnParam.pnts[lnParam.pntsNum-1]));
                }
                lnParam.fstBrkPntVel=GetBrkPntSndPassVel(judgeSs,lnParam.pnts[lnParam.lstBrkPntNum],0,lnParam.lstBrkPntVel,lnParam.pnts[lnParam.fstBrkPntNum]);
            }
            else
            {
                lnParam.fstBrkPntVel=GetBrkPntSndPassVel(judgeSs,lnParam.pnts[sndBrkPntNum],0,0,lnParam.pnts[lnParam.fstBrkPntNum]);
            }
        }
    }
    else
    {
        double judgeSs=lnParam.pnts[0]-lnParam.startV*lnParam.startV/(2*aLmt);
        if(judgeSs>=lnParam.pnts[lnParam.fstBrkPntNum])
        {
            lnParam.fstBrkPntVel=0;
        }
        else
        {
            if(sndBrkPntNum==lnParam.lstBrkPntNum)//pnts[lstBrkPntNum]>pnts[pntsNum-1]
            {
                double judgeSe=lnParam.pnts[lnParam.pntsNum-1]+lnParam.endV*lnParam.endV/(2*aLmt);
                if(judgeSe<=lnParam.pnts[lnParam.lstBrkPntNum])
                {
                    lnParam.lstBrkPntVel=0;
                }
                else
                {
                    lnParam.lstBrkPntVel=-sqrt(lnParam.endV*lnParam.endV-2*aLmt*(lnParam.pnts[lnParam.lstBrkPntNum]-lnParam.pnts[lnParam.pntsNum-1]));
                }
                lnParam.fstBrkPntVel=GetBrkPntSndPassVel(judgeSs,lnParam.pnts[lnParam.lstBrkPntNum],0,lnParam.lstBrkPntVel,lnParam.pnts[lnParam.fstBrkPntNum]);
            }
            else
            {
                lnParam.fstBrkPntVel=GetBrkPntSndPassVel(judgeSs,lnParam.pnts[sndBrkPntNum],0,0,lnParam.pnts[lnParam.fstBrkPntNum]);
            }
        }
    }

    if(lnParam.pnts[lnParam.lstBrkPntNum]>lnParam.pnts[lnParam.pntsNum-1])
    {
        double judgeSe=lnParam.pnts[lnParam.pntsNum-1]+lnParam.endV*lnParam.endV/(2*aLmt);
        if(judgeSe<=lnParam.pnts[lnParam.lstBrkPntNum])
        {
            lnParam.lstBrkPntVel=0;
        }
        else
        {
            lnParam.lstBrkPntVel=-sqrt(lnParam.endV*lnParam.endV-2*aLmt*(lnParam.pnts[lnParam.lstBrkPntNum]-lnParam.pnts[lnParam.pntsNum-1]));
        }
    }
    else
    {
        double judgeSe=lnParam.pnts[lnParam.pntsNum-1]-lnParam.endV*lnParam.endV/(2*aLmt);
        if(judgeSe>=lnParam.pnts[lnParam.lstBrkPntNum])
        {
            lnParam.lstBrkPntVel=0;
        }
        else
        {
            lnParam.lstBrkPntVel=sqrt(lnParam.endV*lnParam.endV+2*aLmt*(lnParam.pnts[lnParam.lstBrkPntNum]-lnParam.pnts[lnParam.pntsNum-1]));
        }
    }

    GetOptimalP2PMotionAcc(lnParam.pnts[0],lnParam.pnts[lnParam.fstBrkPntNum],lnParam.startV,lnParam.fstBrkPntVel,p2pParam);
    lnParam.totalTime=p2pParam.minTime;
    double beginVel=lnParam.fstBrkPntVel;
    double beginPos=lnParam.pnts[lnParam.fstBrkPntNum];
    for(int i=lnParam.fstBrkPntNum+1;i<=lnParam.lstBrkPntNum;i++)
    {
        if((lnParam.pnts[i]-lnParam.pnts[i-1])*(lnParam.pnts[i]-lnParam.pnts[i+1])>0)
        {
            GetOptimalP2PMotionAcc(beginPos,lnParam.pnts[i],beginVel,0,p2pParam_tmp);
            lnParam.totalTime+=p2pParam_tmp.minTime;
            beginPos=lnParam.pnts[i];
            beginVel=0;
        }
        if(i==lnParam.lstBrkPntNum-1)
        {
            GetOptimalP2PMotionAcc(beginPos,lnParam.pnts[lnParam.lstBrkPntNum],beginVel,lnParam.lstBrkPntVel,p2pParam_tmp);
            lnParam.totalTime+=p2pParam_tmp.minTime;
            GetOptimalP2PMotionAcc(lnParam.pnts[lnParam.lstBrkPntNum],lnParam.pnts[lnParam.pntsNum-1],lnParam.lstBrkPntVel,lnParam.endV,p2pParam_tmp);
            lnParam.totalTime+=p2pParam_tmp.minTime;
            break;
        }
    }
    GetThreeNxtPntTime(lnParam,p2pParam);

    printf("Type III p2pParam time: %.4f,%.4f,%.4f, acc: %.4f,%.4f,%.4f, pos: %.4f,%.4f,%.4f\n",
           p2pParam.trajTime[0],p2pParam.trajTime[1],p2pParam.trajTime[2],
           p2pParam.trajAcc[0],p2pParam.trajAcc[1],p2pParam.trajAcc[2],
           p2pParam.trajPos[0],p2pParam.trajPos[1],p2pParam.trajPos[2]);
}

double RTOptimal::GetBrkPntSndPassVel(double startP, double endP, double startV, double endV, double pnt)
{
    double passVel;
    P2PMotionParam p2pParam_tmp;
    GetOptimalP2PMotionAcc(startP,endP,startV,endV,p2pParam_tmp);

    double lmtTimeInterval_tmp[2] {0,p2pParam_tmp.minTime};
    double reachTime_tmp {0};
    reachTime_tmp=LocatePntInLmtInterval(p2pParam_tmp,pnt,lmtTimeInterval_tmp);

    if(reachTime_tmp<=p2pParam_tmp.trajTime[0])
    {
        passVel=p2pParam_tmp.trajAcc[0]*reachTime_tmp;
    }
    else if(reachTime_tmp<p2pParam_tmp.trajTime[0]+p2pParam_tmp.trajTime[1])
    {
        passVel=p2pParam_tmp.maxVel;
    }
    else
    {
        passVel=p2pParam_tmp.maxVel+p2pParam_tmp.trajAcc[2]*(reachTime_tmp-p2pParam_tmp.trajTime[0]-p2pParam_tmp.trajTime[1]);
    }

    return passVel;
}

void RTOptimal::GetThreeNxtPntTime(PolylineParam &lnParam, P2PMotionParam &p2pParam)
{
    //deal with the situation when next point is passed for multi times, most 3 times.
    //three interval: [pnts[0],zeroVelPos[0]],[zeroVelPos[0],zeroVelPos[1]],[zeroVelPos[1],pnts[fstBrkPntNum]]
    int zeroVelPntNum {0};
    double zeroVelPos[2] {0};
    double zeroVelTime[2] {0};
    if((lnParam.pnts[0]+p2pParam.trajPos[0])*lnParam.pnts[0]<0)
    {
        zeroVelPntNum=1;
        zeroVelTime[0]=(0-lnParam.startV)/p2pParam.trajAcc[0];
        zeroVelPos[0]=lnParam.startV*zeroVelTime[0]+0.5*p2pParam.trajAcc[0]*zeroVelTime[0]*zeroVelTime[0];
        if((lnParam.pnts[lnParam.fstBrkPntNum]-p2pParam.trajPos[2])*lnParam.pnts[lnParam.fstBrkPntNum]<0)
        {
            zeroVelPntNum=2;
            zeroVelTime[1]=(lnParam.fstBrkPntVel-0)/p2pParam.trajAcc[2];
            zeroVelPos[1]=lnParam.pnts[lnParam.fstBrkPntNum]-0.5*p2pParam.trajAcc[2]*zeroVelTime[1]*zeroVelTime[1];
        }
    }
    else
    {
        if((lnParam.pnts[lnParam.fstBrkPntNum]-p2pParam.trajPos[2])*lnParam.pnts[lnParam.fstBrkPntNum]<0)
        {
            zeroVelPntNum=1;
            zeroVelTime[0]=(lnParam.fstBrkPntVel-0)/p2pParam.trajAcc[2];
            zeroVelPos[0]=lnParam.pnts[lnParam.fstBrkPntNum]-lnParam.pnts[0]-0.5*p2pParam.trajAcc[2]*zeroVelTime[0]*zeroVelTime[0];
        }
    }

    lnParam.passNxtPntNum=0;
    std::fill_n(lnParam.passNxtPntTime,3,0);
    double d=lnParam.pnts[1]-lnParam.pnts[0];
    double dist=lnParam.pnts[lnParam.fstBrkPntNum]-lnParam.pnts[0];
    double lmtTimeInterval[2];
    if(d*(d-zeroVelPos[0])<=0 && d-zeroVelPos[0]!=0)
    {
        lmtTimeInterval[0]=0;
        lmtTimeInterval[1]=zeroVelTime[0];
        lnParam.passNxtPntTime[lnParam.passNxtPntNum]=LocatePntInLmtInterval(p2pParam,lnParam.pnts[1],lmtTimeInterval);
        lnParam.passNxtPntNum++;
    }
    if((d-zeroVelPos[0])*(d-zeroVelPos[1])<=0 && d-zeroVelPos[1]!=0)
    {
        lmtTimeInterval[0]=zeroVelTime[0];
        lmtTimeInterval[1]=zeroVelTime[1];
        lnParam.passNxtPntTime[lnParam.passNxtPntNum]=LocatePntInLmtInterval(p2pParam,lnParam.pnts[1],lmtTimeInterval);
        lnParam.passNxtPntNum++;
    }
    if((d-zeroVelPos[1])*(d-dist)<=0 && d-dist!=0)
    {
        lmtTimeInterval[0]=zeroVelTime[1];
        lmtTimeInterval[1]=p2pParam.minTime;
        lnParam.passNxtPntTime[lnParam.passNxtPntNum]=LocatePntInLmtInterval(p2pParam,lnParam.pnts[1],lmtTimeInterval);
        lnParam.passNxtPntNum++;
    }
    if(d-dist==0)
    {
        lnParam.passNxtPntTime[lnParam.passNxtPntNum]=p2pParam.minTime;
        lnParam.passNxtPntNum++;
    }
    if(lnParam.passNxtPntNum==4)
    {
        printf("Error! Next point is passed for 4 times. But it should be mostly 3 theoretically.\n");
    }
}

double RTOptimal::LocatePntInLmtInterval(P2PMotionParam &p2pParam, double pnt, double *limitedTimeInterval)
{
    double reachT[5];
    if((pnt-p2pParam.startPos)*(pnt-p2pParam.startPos-p2pParam.trajPos[0])<=0)
    {
        double d=pnt-p2pParam.startPos;
        double v_square=p2pParam.startVel*p2pParam.startVel+2*p2pParam.trajAcc[0]*d;
        if(v_square<0)
        {
            if(v_square>-1e-10)
            {
                v_square=0;
            }
            else
            {
                printf("Error! sqrt applied to negative!\n");
            }
        }
        reachT[0]=(-p2pParam.startVel+sqrt(v_square))/p2pParam.trajAcc[0];
        reachT[1]=(-p2pParam.startVel-sqrt(v_square))/p2pParam.trajAcc[0];
    }
    if((pnt-p2pParam.startPos-p2pParam.trajPos[0])*(pnt-p2pParam.startPos-p2pParam.trajPos[0]-p2pParam.trajPos[1])<=0)
    {
        reachT[2]=(pnt-p2pParam.startPos-p2pParam.trajPos[0])/p2pParam.maxVel;
    }
    if((pnt-p2pParam.startPos-p2pParam.trajPos[0]-p2pParam.trajPos[1])*(pnt-p2pParam.endPos)<=0)
    {
        double d=pnt-p2pParam.startPos-p2pParam.trajPos[0]-p2pParam.trajPos[1];
        double v_square=p2pParam.maxVel*p2pParam.maxVel+2*p2pParam.trajAcc[2]*d<0;
        if(v_square<0)
        {
            if(v_square>-1e-10)
            {
                v_square=0;
            }
            else
            {
                printf("Error! sqrt applied to negative!\n");
            }
        }
        reachT[3]=(-p2pParam.maxVel+sqrt(v_square))/p2pParam.trajAcc[2];
        reachT[4]=(-p2pParam.maxVel-sqrt(v_square))/p2pParam.trajAcc[2];
    }

    bool isFound {false};
    double reachTime {-1};
    for(int i=0;i<5;i++)
    {
        if(reachT[i]>=limitedTimeInterval[0] && reachT[i]<limitedTimeInterval[1])
        {
            isFound=true;
            reachTime=reachT[i];
            break;
        }
    }
    if(isFound==true)
    {
        return reachTime;
    }
    else
    {
        printf("Error! Conflict! Next point time cannot be calaulated, but it is judged to be here!\n");
        return 0;
    }
}

void RTOptimal::DecideNxtPntTime(PolylineParam *lnParam)
{
    //find optimal method to traverse next point in mostly 3*3*3=27 cases.
    double maxTime {0};
    int nums[3] {0};
    for(int i=0;i<lnParam[0].passNxtPntNum;i++)
    {
        for(int j=0;j<lnParam[1].passNxtPntNum;j++)
        {
            for(int k=0;k<lnParam[2].passNxtPntNum;k++)
            {
                double nxtPntTime_tmp[3] {lnParam[0].passNxtPntTime[i],lnParam[1].passNxtPntTime[j],lnParam[2].passNxtPntTime[k]};
                double remainTime_tmp[3] {lnParam[0].totalTime-lnParam[0].passNxtPntTime[i],lnParam[1].totalTime-lnParam[1].passNxtPntTime[j],lnParam[2].totalTime-lnParam[2].passNxtPntTime[k]};
                double maxNxtPntTime=*std::max_element(nxtPntTime_tmp,nxtPntTime_tmp+3);
                double maxRemainTime=*std::max_element(remainTime_tmp,remainTime_tmp+3);
                double maxTime_tmp=maxNxtPntTime+maxRemainTime;
                if(maxTime_tmp>maxTime)
                {
                    maxTime=maxTime_tmp;
                    nums[0]=i;
                    nums[1]=j;
                    nums[2]=k;
                }
            }
        }
    }
    for(int i=0;i<3;i++)
    {
        lnParam[i].nxtPntTime=lnParam[i].passNxtPntTime[nums[i]];
    }
    double nxtPntTime_tmp[3] {lnParam[0].nxtPntTime,lnParam[1].nxtPntTime,lnParam[2].nxtPntTime};
    actScrewID=std::max_element(nxtPntTime_tmp,nxtPntTime_tmp+3)-nxtPntTime_tmp;
    nxtPntMinTime=nxtPntTime_tmp[actScrewID];
    printf("nxtPnttime:%.4f,%.4f,%.4f; nxtPntMinTime:%.4f\n",nxtPntTime_tmp[0],nxtPntTime_tmp[1],nxtPntTime_tmp[2],nxtPntMinTime);
}

double RTOptimal::GetNxtPin(PolylineParam &lnParam, NextPointParam &nxtParam, int count)
{
    double pIn;
    double curT=0.001*count-lstPntTime;
    if(nxtParam.optVel<=nxtParam.reachableVel[0])
    {
        if(curT<=nxtParam.min_t1)
        {
            pIn=lnParam.pnts[0]+lnParam.startV*curT+0.5*aLmt*curT*curT;
            if(curT<=0.001)
            {
                nxtParam.record_t1=nxtParam.min_t1;
                nxtParam.record_t2=nxtParam.min_t2;
                printf("Unreachable NxtPnt time:%.4f,%.4f, acc:%.4f,%.4f\n",nxtParam.min_t1,nxtParam.min_t2,aLmt,-aLmt);
            }
        }
        else if(curT<=nxtParam.min_t1+nxtParam.min_t2)
        {
            pIn=lnParam.pnts[0]+lnParam.startV*nxtParam.min_t1+0.5*aLmt*nxtParam.min_t1*nxtParam.min_t1
                    +nxtParam.min_vm*(curT-nxtParam.min_t1)-0.5*aLmt*(curT-nxtParam.min_t1)*(curT-nxtParam.min_t1);

            if(curT+0.001>nxtParam.min_t1+nxtParam.min_t2)
            {
                printf("Current target point passed! Go on to next one!\n");
                for(int j=0;j<lnParam.pntsNum+1;j++)
                {
                    lnParam.pnts[j]=lnParam.pnts[j+1];
                }
                lnParam.pntsNum--;
                lnParam.startV=nxtParam.reachableVel[0];
                isCurPntPassed=true;
            }
        }
        else
        {
            printf("Error! Impossible to enter here! currentT=%.4f, min_T=%.4f,%.4f\n",curT,nxtParam.min_t1,nxtParam.min_t2);
        }
    }
    else if(nxtParam.optVel<=nxtParam.reachableVel[1])
    {
        if(nxtParam.optVel==lnParam.startV)
        {
            nxtParam.rch_t1=nxtPntMinTime/2;
        }
        else
        {
            double s=lnParam.pnts[1]-lnParam.pnts[0];
            double delta=(lnParam.startV*nxtPntMinTime-s)*(lnParam.startV*nxtPntMinTime-s)/2
                        +(nxtParam.optVel*nxtPntMinTime-s)*(nxtParam.optVel*nxtPntMinTime-s)/2;
            double t1=(-(nxtParam.optVel*nxtPntMinTime-s)+sqrt(delta))/(lnParam.startV-nxtParam.optVel);
            double t2=(-(nxtParam.optVel*nxtPntMinTime-s)-sqrt(delta))/(lnParam.startV-nxtParam.optVel);
            if(t1>=0 && t1<=nxtPntMinTime)
            {
                nxtParam.rch_t1=t1;
            }
            else if(t2>=0 && t2<=nxtPntMinTime)
            {
                nxtParam.rch_t1=t2;
            }
            else
            {
                printf("Error! optVel cannot be reached! t1=%.4f & t2=%.4f are negative!\n",t1,t2);
            }
        }

        double acc=(nxtParam.optVel-lnParam.startV)/(2*nxtParam.rch_t1-nxtPntMinTime);
        nxtParam.rch_t2=nxtPntMinTime-nxtParam.rch_t1;
        nxtParam.rch_vm=lnParam.startV+acc*nxtParam.rch_t1;

        if(curT<=nxtParam.rch_t1)
        {
            pIn=lnParam.pnts[0]+lnParam.startV*curT+0.5*acc*curT*curT;
            if(curT<=0.001)
            {
                nxtParam.record_t1=nxtParam.rch_t1;
                nxtParam.record_t2=nxtParam.rch_t2;
                printf("Reachable NxtPnt time:%.4f,%.4f, acc:%.4f,%.4f\n",nxtParam.rch_t1,nxtParam.rch_t2,acc,-acc);
            }
        }
        else if(curT<=nxtParam.rch_t1+nxtParam.rch_t2)
        {
            pIn=lnParam.pnts[0]+lnParam.startV*nxtParam.rch_t1+0.5*acc*nxtParam.rch_t1*nxtParam.rch_t1
                    +nxtParam.rch_vm*(curT-nxtParam.rch_t1)-0.5*acc*(curT-nxtParam.rch_t1)*(curT-nxtParam.rch_t1);

            if(curT+0.001>nxtParam.rch_t1+nxtParam.rch_t2)
            {
                printf("Current target point passed! Go on to next one!\n");
                for(int j=0;j<lnParam.pntsNum+1;j++)
                {
                    lnParam.pnts[j]=lnParam.pnts[j+1];
                }
                lnParam.pntsNum--;
                lnParam.startV=nxtParam.optVel;
                isCurPntPassed=true;
            }
        }
        else
        {
            printf("Error! Impossible to enter here! currentT=%.4f, rch_T=%.4f,%.4f\n",curT,nxtParam.rch_t1,nxtParam.rch_t2);
        }
    }
    else
    {
        if(curT<=nxtParam.max_t1)
        {
            pIn=lnParam.pnts[0]+lnParam.startV*curT-0.5*aLmt*curT*curT;
            if(curT<=0.001)
            {
                nxtParam.record_t1=nxtParam.max_t1;
                nxtParam.record_t2=nxtParam.max_t2;
                printf("UnReachable NxtPnt time:%.4f,%.4f, acc:%.4f,%.4f\n",nxtParam.max_t1,nxtParam.max_t2,-aLmt,aLmt);
            }
        }
        else if(curT<nxtParam.max_t1+nxtParam.max_t2)
        {
            pIn=lnParam.pnts[0]+lnParam.startV*nxtParam.max_t1-0.5*aLmt*nxtParam.max_t1*nxtParam.max_t1
                    +nxtParam.max_vm*(curT-nxtParam.max_t1)+0.5*aLmt*(curT-nxtParam.max_t1)*(curT-nxtParam.max_t1);

            if(curT+0.001>nxtParam.max_t1+nxtParam.max_t2)
            {
                printf("Current target point passed! Go on to next one!\n");
                for(int j=0;j<lnParam.pntsNum+1;j++)
                {
                    lnParam.pnts[j]=lnParam.pnts[j+1];
                }
                lnParam.pntsNum--;
                lnParam.startV=nxtParam.reachableVel[1];
                isCurPntPassed=true;
            }
        }
        else
        {
            printf("Error! Impossible to enter here! currentT=%.4f, max_T=%.4f,%.4f\n",curT,nxtParam.max_t1,nxtParam.max_t2);
        }
    }

    return pIn;
}

void RTOptimal::GetNxtPntOptVel(PolylineParam &lnParam, NextPointParam &nxtParam)
{
    double nxtBrkPnt;
    double nxtBrkPntVel;
    switch(lnParam.lineType)
    {
    case PolylineType::I:
        nxtBrkPnt=lnParam.pnts[lnParam.pntsNum-1];
        nxtBrkPntVel=lnParam.endV;
        break;
    case PolylineType::II:
    case PolylineType::III:
        nxtBrkPnt=lnParam.pnts[lnParam.fstBrkPntNum];
        nxtBrkPntVel=lnParam.fstBrkPntVel;
        break;
    default:
        break;
    }
    if(lnParam.pntsNum==2)
    {
        nxtParam.optVel=lnParam.endV;
    }
    else
    {
        if(nxtBrkPnt>=lnParam.pnts[0])
        {
            nxtParam.optVel=sqrt(nxtBrkPntVel*nxtBrkPntVel+2*aLmt*(nxtBrkPnt-lnParam.pnts[1]));
            if(nxtParam.optVel>vLmt)
            {
                nxtParam.optVel=vLmt;
            }
        }
        else
        {
            nxtParam.optVel=-sqrt(nxtBrkPntVel*nxtBrkPntVel-2*aLmt*(nxtBrkPnt-lnParam.pnts[1]));
            if(nxtParam.optVel<-vLmt)
            {
                nxtParam.optVel=-vLmt;
            }
        }
    }
    printf("OptVel:%.4f\n",nxtParam.optVel);
}

bool RTOptimal::GetNxtPntReachableVel(PolylineParam &lnParam, NextPointParam &nxtParam, int screwID)
{
    double s=lnParam.pnts[1]-lnParam.pnts[0];
    double max_s=lnParam.startV*nxtPntMinTime+0.5*aLmt*nxtPntMinTime*nxtPntMinTime;
    double min_s=lnParam.startV*nxtPntMinTime-0.5*aLmt*nxtPntMinTime*nxtPntMinTime;

    //nxtPntMinTime in inoperative interval, new time need to be figured out
    if(s>max_s)//startV<0,s<0,a>0
    {
        printf("Error! Screw %d cannot reach next point within the min time %.4f.\n",actScrewID,nxtPntMinTime);
        //nxtPntMinTime=(-lnParam.startV-sqrt(lnParam.startV*lnParam.startV+2*aLmt*s))/aLmt;
        //actScrewID=screwID;
        return false;
    }
    else if(s<min_s)//startV>0,s>0,a<0
    {
        printf("Error! Screw %d cannot reach next point within the min time %.4f.\n",actScrewID,nxtPntMinTime);
        //nxtPntMinTime=(-lnParam.startV+sqrt(lnParam.startV*lnParam.startV-2*aLmt*s))/(-aLmt);
        //actScrewID=screwID;
        return false;
    }
    else
    {
        nxtParam.min_t2=sqrt(aLmt*lnParam.startV*nxtPntMinTime+0.5*aLmt*aLmt*nxtPntMinTime*nxtPntMinTime-aLmt*s)/aLmt;
        nxtParam.min_t1=nxtPntMinTime-nxtParam.min_t2;
        nxtParam.min_vm=lnParam.startV+aLmt*nxtParam.min_t1;
        nxtParam.reachableVel[0]=lnParam.startV+aLmt*nxtParam.min_t1-aLmt*nxtParam.min_t2;

        nxtParam.max_t2=sqrt(-aLmt*lnParam.startV*nxtPntMinTime+0.5*aLmt*aLmt*nxtPntMinTime*nxtPntMinTime+aLmt*s)/aLmt;
        nxtParam.max_t1=nxtPntMinTime-nxtParam.max_t2;
        nxtParam.max_vm=lnParam.startV-aLmt*nxtParam.max_t1;
        nxtParam.reachableVel[1]=lnParam.startV-aLmt*nxtParam.max_t1+aLmt*nxtParam.max_t2;

        printf("reachableVel:%.4f,%.4f\n",nxtParam.reachableVel[0],nxtParam.reachableVel[1]);
        printf("min_vm=%.4f, max_vm=%.4f\n",nxtParam.min_vm,nxtParam.max_vm);
        return true;
    }
}

void RTOptimal::GetOptimalP2PMotionAcc(double startP, double endP, double startV, double endV, P2PMotionParam &p2pParam)
{
    if(startP<=endP)
    {
        double midV {0};
        double judgeDist1 {0};
        double judgeDist2 {0};

        p2pParam.startPos=startP;
        p2pParam.startVel=startV;
        p2pParam.endPos=endP;
        p2pParam.endVel=endV;
        p2pParam.inopInter[0]=-1;
        p2pParam.inopInter[1]=-1;

        //type judgement of inoperative time interval
        if(startV<=0)
        {
            if(endV<=0)
            {
                //Type I -- acc-pos-dec
                p2pParam.inopInterNum=0;
                p2pParam.velType=VelType::AccDec;
                midV=sqrt(startV*startV/2+endV*endV/2+aLmt*(endP-startP));
                GetP2PMotionParam(p2pParam,midV);
                p2pParam.minTime=p2pParam.trajTime[0]+p2pParam.trajTime[1]+p2pParam.trajTime[2];
            }
            else
            {
                judgeDist1=(endV*endV-startV*startV)/(2*aLmt);
                if(endP-startP<judgeDist1)//Type III -- dec-neg-acc
                {
                    p2pParam.inopInterNum=0;
                    p2pParam.velType=VelType::DecAcc;
                    midV=-sqrt(startV*startV/2+endV*endV/2-aLmt*(endP-startP));
                    GetP2PMotionParam(p2pParam,midV);
                    p2pParam.minTime=p2pParam.trajTime[0]+p2pParam.trajTime[1]+p2pParam.trajTime[2];
                }
                else//Type I -- acc-pos-dec
                {
                    p2pParam.inopInterNum=0;
                    p2pParam.velType=VelType::AccDec;
                    midV=sqrt(startV*startV/2+endV*endV/2+aLmt*(endP-startP));
                    GetP2PMotionParam(p2pParam,midV);
                    p2pParam.minTime=p2pParam.trajTime[0]+p2pParam.trajTime[1]+p2pParam.trajTime[2];
                }
            }
        }
        else
        {
            if(endV<=0)
            {
                judgeDist1=(startV*startV-endV*endV)/(2*aLmt);
                if(endP-startP<judgeDist1)//Type III -- dec-neg-acc
                {
                    p2pParam.inopInterNum=0;
                    p2pParam.velType=VelType::DecAcc;
                    midV=-sqrt(startV*startV/2+endV*endV/2-aLmt*(endP-startP));
                    GetP2PMotionParam(p2pParam,midV);
                    p2pParam.minTime=p2pParam.trajTime[0]+p2pParam.trajTime[1]+p2pParam.trajTime[2];
                }
                else//Type I -- acc-pos-dec
                {
                    p2pParam.inopInterNum=0;
                    p2pParam.velType=VelType::AccDec;
                    midV=sqrt(startV*startV/2+endV*endV/2+aLmt*(endP-startP));
                    GetP2PMotionParam(p2pParam,midV);
                    p2pParam.minTime=p2pParam.trajTime[0]+p2pParam.trajTime[1]+p2pParam.trajTime[2];
                }
            }
            else
            {
                judgeDist1=fabs(startV*startV-endV*endV)/(2*aLmt);
                judgeDist2=(startV*startV+endV*endV)/(2*aLmt);
                if(endP-startP<=judgeDist1)//Type III -- dec-neg-acc
                {
                    p2pParam.inopInterNum=0;
                    p2pParam.velType=VelType::DecAcc;
                    midV=-sqrt(startV*startV/2+endV*endV/2-aLmt*(endP-startP));
                    GetP2PMotionParam(p2pParam,midV);
                    p2pParam.minTime=p2pParam.trajTime[0]+p2pParam.trajTime[1]+p2pParam.trajTime[2];
                }
                else if(endP-startP>=judgeDist2)//Type I -- acc-pos-dec
                {
                    p2pParam.inopInterNum=0;
                    p2pParam.velType=VelType::AccDec;
                    midV=sqrt(startV*startV/2+endV*endV/2+aLmt*(endP-startP));
                    GetP2PMotionParam(p2pParam,midV);
                    p2pParam.minTime=p2pParam.trajTime[0]+p2pParam.trajTime[1]+p2pParam.trajTime[2];
                }
                else//All type
                {
                    p2pParam.inopInterNum=1;

                    //Type II -- dec-acc
                    p2pParam.velType=VelType::DecAcc;
                    midV=sqrt(startV*startV/2+endV*endV/2-aLmt*(endP-startP));
                    GetP2PMotionParam(p2pParam,midV);
                    p2pParam.inopInter[0]=p2pParam.trajTime[0]+p2pParam.trajTime[1]+p2pParam.trajTime[2];

                    //Type III -- dec-neg-acc
                    p2pParam.velType=VelType::DecAcc;
                    midV=-sqrt(startV*startV/2+endV*endV/2-aLmt*(endP-startP));
                    GetP2PMotionParam(p2pParam,midV);
                    p2pParam.inopInter[1]=p2pParam.trajTime[0]+p2pParam.trajTime[1]+p2pParam.trajTime[2];

                    //Type I -- acc-pos-dec
                    p2pParam.velType=VelType::AccDec;
                    midV=sqrt(startV*startV/2+endV*endV/2+aLmt*(endP-startP));
                    GetP2PMotionParam(p2pParam,midV);
                    p2pParam.minTime=p2pParam.trajTime[0]+p2pParam.trajTime[1]+p2pParam.trajTime[2];
                }
            }
        }
    }
    else
    {
        GetOptimalP2PMotionAcc(-startP,-endP,-startV,-endV,p2pParam);
        for(int i=0;i<3;i++)
        {
            p2pParam.trajAcc[i]=-p2pParam.trajAcc[i];
            p2pParam.trajPos[i]=-p2pParam.trajPos[i];
        }
        p2pParam.startPos=-p2pParam.startPos;
        p2pParam.startVel=-p2pParam.startVel;
        p2pParam.endPos=-p2pParam.endPos;
        p2pParam.endVel=-p2pParam.endVel;
        p2pParam.maxVel=-p2pParam.maxVel;
        if(p2pParam.velType==VelType::AccDec)
        {
            p2pParam.velType=VelType::DecAcc;
        }
        else if(p2pParam.velType==VelType::DecAcc)
        {
            p2pParam.velType=VelType::AccDec;
        }
        else if(p2pParam.velType==VelType::AccConstDec)
        {
            p2pParam.velType=VelType::DecConstAcc;
        }
        else if(p2pParam.velType==VelType::DecConstAcc)
        {
            p2pParam.velType=VelType::AccConstDec;
        }
        else
        {
            printf("Error! Unknown VelType!\n");
        }
    }
}

void RTOptimal::GetP2PMotionParam(P2PMotionParam &p2pParam, double midV)
{
    switch(p2pParam.velType)
    {
    case VelType::DecAcc:
        if(midV<-vLmt)//constant phase
        {
            p2pParam.velType=VelType::DecConstAcc;
            p2pParam.maxVel=-vLmt;
            p2pParam.trajTime[0]=(p2pParam.startVel+vLmt)/aLmt;
            p2pParam.trajAcc[0]=-aLmt;
            p2pParam.trajTime[2]=(p2pParam.endVel+vLmt)/aLmt;
            p2pParam.trajAcc[2]=aLmt;

            p2pParam.trajPos[0]=(vLmt*vLmt-p2pParam.startVel*p2pParam.startVel)/(-2*aLmt);
            p2pParam.trajPos[2]=(p2pParam.endVel*p2pParam.endVel-vLmt*vLmt)/(2*aLmt);
            p2pParam.trajPos[1]=p2pParam.endPos-p2pParam.startPos-p2pParam.trajPos[0]-p2pParam.trajPos[2];
            if(p2pParam.trajPos[1]>0)
            {
                printf("Error! Remain s should be negative. remainS=%.2f.\n",p2pParam.trajPos[1]);
                std::abort();
            }
            else
            {
                p2pParam.trajTime[1]=p2pParam.trajPos[1]/(-vLmt);
                p2pParam.trajAcc[1]=0;

            }
        }
        else//none constant phase
        {
            p2pParam.maxVel=midV;
            p2pParam.trajTime[0]=(p2pParam.startVel-midV)/aLmt;
            p2pParam.trajAcc[0]=-aLmt;
            p2pParam.trajPos[0]=(midV*midV-p2pParam.startVel*p2pParam.startVel)/(-2*aLmt);
            p2pParam.trajTime[1]=0;
            p2pParam.trajAcc[1]=0;
            p2pParam.trajPos[1]=0;
            p2pParam.trajTime[2]=(p2pParam.endVel-midV)/aLmt;
            p2pParam.trajAcc[2]=aLmt;
            p2pParam.trajPos[2]=(p2pParam.endVel*p2pParam.endVel-midV*midV)/(2*aLmt);
        }

        break;
    case VelType::AccDec:
        if(midV>vLmt)//constant phase
        {
            p2pParam.velType=VelType::AccConstDec;
            p2pParam.maxVel=vLmt;
            p2pParam.trajTime[0]=(vLmt-p2pParam.startVel)/aLmt;
            p2pParam.trajAcc[0]=aLmt;
            p2pParam.trajTime[2]=(vLmt-p2pParam.endVel)/aLmt;
            p2pParam.trajAcc[2]=-aLmt;

            p2pParam.trajPos[0]=(vLmt*vLmt-p2pParam.startVel*p2pParam.startVel)/(2*aLmt);
            p2pParam.trajPos[2]=(vLmt*vLmt-p2pParam.endVel*p2pParam.endVel)/(2*aLmt);
            p2pParam.trajPos[1]=p2pParam.endPos-p2pParam.startPos-p2pParam.trajPos[0]-p2pParam.trajPos[2];
            if(p2pParam.trajPos[1]<0)
            {
                printf("Error! Remain s should be positive. remainS=%.2f.\n",p2pParam.trajPos[1]);
                std::abort();
            }
            else
            {
                p2pParam.trajTime[1]=p2pParam.trajPos[1]/vLmt;
                p2pParam.trajAcc[1]=0;
            }
        }
        else//none constant phase
        {
            p2pParam.maxVel=midV;
            p2pParam.trajTime[0]=(midV-p2pParam.startVel)/aLmt;
            p2pParam.trajAcc[0]=aLmt;
            p2pParam.trajPos[0]=(midV*midV-p2pParam.startVel*p2pParam.startVel)/(2*aLmt);
            p2pParam.trajTime[1]=0;
            p2pParam.trajAcc[1]=0;
            p2pParam.trajPos[1]=0;
            p2pParam.trajTime[2]=(midV-p2pParam.endVel)/aLmt;
            p2pParam.trajAcc[2]=-aLmt;
            p2pParam.trajPos[2]=(p2pParam.endVel*p2pParam.endVel-midV*midV)/(-2*aLmt);
        }

        break;
    default:
        break;
    }

    if(p2pParam.trajTime[0]<0 || p2pParam.trajTime[1]<0 || p2pParam.trajTime[2]<0)
    {
        printf("Error! Time impossible to be negative. T=(%.2f,%.2f,%.2f).\n",
              p2pParam.trajTime[0],p2pParam.trajTime[1],p2pParam.trajTime[2]);
        std::abort();
    }
}

void RTOptimal::GetOptimalP2PMotionJerk(double startP, double endP, double startV, double endV)
{

}

void RTOptimal::RecordNxtPntTime(NextPointParam &nxtParam, int count, int screwID)
{
    nxtPntTime[count][3*screwID]=nxtParam.record_t1;
    nxtPntTime[count][3*screwID+1]=nxtParam.record_t2;
    nxtPntTime[count][3*screwID+2]=nxtParam.record_t1+nxtParam.record_t2;
}
