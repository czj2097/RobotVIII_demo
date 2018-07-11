#include "RTOptimal.h"
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
bool RTOptimal::GetTrajOneLeg(int legID, int count, double *pIn, double *pEE)
{
    if(isCurPntPassed==true)
    {
        printf("\n\n\n\n");
        GetAllParamForNxtPin(legID);
    }

    bool isNxtPinGot {false};
    while(isNxtPinGot==false)
    {
        isNxtPinGot=GetNxtPin(lineParam[3*legID],nxtPntParam[3*legID],count,pIn[3*legID]);
        if(isNxtPinGot==true)
        {
            GetNxtPin(lineParam[3*legID+1],nxtPntParam[3*legID+1],count,pIn[3*legID+1]);
            GetNxtPin(lineParam[3*legID+2],nxtPntParam[3*legID+2],count,pIn[3*legID+2]);
        }
        else
        {
            GetAllParamForNxtPin(legID);
        }
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

void RTOptimal::GetAllParamForNxtPin(int legID)
{
//    for(int i=0;i<3;i++)
//    {
//        printf("pnts of screw %d:",i);
//        for(int k=0;k<lineParam[3*legID+i].pntsNum;k++)
//        {
//            printf("%.4f,",lineParam[3*legID+i].pnts[k]);
//        }
//        printf("\n");
//    }
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

//        if(lineParam[3*legID].lineType==PolylineType::I && lineParam[3*legID+1].lineType==PolylineType::I && lineParam[3*legID+2].lineType==PolylineType::I)
//        {
//            double p2pMotTime_tmp[3] {curP2PParam[3*legID].minTime,curP2PParam[3*legID+1].minTime,curP2PParam[3*legID+2].minTime};
//            int p2pMaxID=std::max_element(p2pMotTime_tmp,p2pMotTime_tmp+3)-p2pMotTime_tmp;
//            double p2pMotMaxTime=p2pMotTime_tmp[p2pMaxID];
//            double p2pScaleTime=0.001*((int)((p2pMotMaxTime+lstPntTime)*1000)+1)-lstPntTime;
//            ScaleP2PMotionAcc(curP2PParam[3*legID+p2pMaxID],p2pScaleTime);
//        }

    GetFinalNxtPntTime(lineParam+3*legID);
    GetNxtPntReachableVel(lineParam+3*legID,nxtPntParam+3*legID);
    for(int i=0;i<3;i++)
    {
        GetNxtPntOptVel(lineParam[3*legID+i],nxtPntParam[3*legID+i]);
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
    }

//    for(int i=0;i<pntsNum;i++)
//    {
//        printf("knot %d:",i);
//        for(int k=0;k<3;k++)
//        {
//            printf("%.4f,",pnts[3*i+k]);
//        }
//        printf("\n");
//    }
    printf("startVin: %.4f,%.4f,%.4f; endVin:%.4f,%.4f,%.4f\n",
           lineParam[3*legID].startV,lineParam[3*legID+1].startV,lineParam[3*legID+2].startV,
           lineParam[3*legID].endV,lineParam[3*legID+1].endV,lineParam[3*legID+2].endV);
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
    if(zeroVelPntNum==0)
    {
        printf("num=0\n");
        //printf("num=0, d-dist=%.10f\n",d-dist);
        if(d*(d-dist)<=0 && d-dist!=0)
        {
            printf("passNxtPntTime=0\n");
            lmtTimeInterval[0]=0;
            lmtTimeInterval[1]=p2pParam.minTime;
            lnParam.passNxtPntTime[lnParam.passNxtPntNum]=LocatePntInLmtInterval(p2pParam,lnParam.pnts[1],lmtTimeInterval);
            lnParam.passNxtPntNum++;
        }
        if(d-dist==0)
        {
            printf("passNxtPntTime=1\n");
            lnParam.passNxtPntTime[lnParam.passNxtPntNum]=p2pParam.minTime;
            lnParam.passNxtPntNum++;
        }
    }
    else if(zeroVelPntNum==1)
    {
        printf("num=1\n");
        if(d*(d-zeroVelPos[0])<=0 && d-zeroVelPos[0]!=0)
        {
            printf("passNxtPntTime=0\n");
            lmtTimeInterval[0]=0;
            lmtTimeInterval[1]=zeroVelTime[0];
            lnParam.passNxtPntTime[lnParam.passNxtPntNum]=LocatePntInLmtInterval(p2pParam,lnParam.pnts[1],lmtTimeInterval);
            lnParam.passNxtPntNum++;
        }
        if((d-zeroVelPos[0])*(d-dist)<=0 && d-dist!=0)
        {
            printf("passNxtPntTime=1\n");
            lmtTimeInterval[0]=zeroVelTime[0];
            lmtTimeInterval[1]=zeroVelTime[1];
            lnParam.passNxtPntTime[lnParam.passNxtPntNum]=LocatePntInLmtInterval(p2pParam,lnParam.pnts[1],lmtTimeInterval);
            lnParam.passNxtPntNum++;
        }
        if(d-dist==0)
        {
            printf("passNxtPntTime=2\n");
            lnParam.passNxtPntTime[lnParam.passNxtPntNum]=p2pParam.minTime;
            lnParam.passNxtPntNum++;
        }
    }
    else
    {
        printf("num=2\n");
        if(d*(d-zeroVelPos[0])<=0 && d-zeroVelPos[0]!=0)
        {
            printf("passNxtPntTime=0\n");
            lmtTimeInterval[0]=0;
            lmtTimeInterval[1]=zeroVelTime[0];
            lnParam.passNxtPntTime[lnParam.passNxtPntNum]=LocatePntInLmtInterval(p2pParam,lnParam.pnts[1],lmtTimeInterval);
            lnParam.passNxtPntNum++;
        }
        if((d-zeroVelPos[0])*(d-zeroVelPos[1])<=0 && d-zeroVelPos[1]!=0)
        {
            printf("passNxtPntTime=1\n");
            lmtTimeInterval[0]=zeroVelTime[0];
            lmtTimeInterval[1]=zeroVelTime[1];
            lnParam.passNxtPntTime[lnParam.passNxtPntNum]=LocatePntInLmtInterval(p2pParam,lnParam.pnts[1],lmtTimeInterval);
            lnParam.passNxtPntNum++;
        }
        if((d-zeroVelPos[1])*(d-dist)<=0 && d-dist!=0)
        {
            printf("passNxtPntTime=2\n");
            lmtTimeInterval[0]=zeroVelTime[1];
            lmtTimeInterval[1]=p2pParam.minTime;
            lnParam.passNxtPntTime[lnParam.passNxtPntNum]=LocatePntInLmtInterval(p2pParam,lnParam.pnts[1],lmtTimeInterval);
            lnParam.passNxtPntNum++;
        }
        if(d-dist==0)
        {
            printf("passNxtPntTime=3\n");
            lnParam.passNxtPntTime[lnParam.passNxtPntNum]=p2pParam.minTime;
            lnParam.passNxtPntNum++;
        }
    }

    if(lnParam.passNxtPntNum==4)
    {
        printf("Error! Next point is passed for 4 times. But it should be mostly 3 theoretically.\n");
    }
    printf("passNxtPntTime:%f,%f,%f,%f\n",lnParam.passNxtPntTime[0],lnParam.passNxtPntTime[1],lnParam.passNxtPntTime[2],lnParam.passNxtPntTime[3]);
}

double RTOptimal::LocatePntInLmtInterval(P2PMotionParam &p2pParam, double pnt, double *limitedTimeInterval)
{
    double reachT[5];
    std::fill_n(reachT,5,-1);
    if((pnt-p2pParam.startPos)*(pnt-p2pParam.startPos-p2pParam.trajPos[0])<=0)
    {
        printf("locate=0\n");
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
        printf("startVel:%.10f,acc:%.10f,reachT:%.10f,%.10f\n",p2pParam.startVel,p2pParam.trajAcc[0],reachT[0],reachT[1]);
    }
    if((pnt-p2pParam.startPos-p2pParam.trajPos[0])*(pnt-p2pParam.startPos-p2pParam.trajPos[0]-p2pParam.trajPos[1])<=0)
    {
        printf("locate=1\n");
        reachT[2]=p2pParam.trajTime[0]+(pnt-p2pParam.startPos-p2pParam.trajPos[0])/p2pParam.maxVel;
    }
    if((pnt-p2pParam.startPos-p2pParam.trajPos[0]-p2pParam.trajPos[1])*(pnt-p2pParam.endPos)<=0)
    {
        printf("locate=2\n");
        double d=pnt-p2pParam.startPos-p2pParam.trajPos[0]-p2pParam.trajPos[1];
        double v_square=p2pParam.maxVel*p2pParam.maxVel+2*p2pParam.trajAcc[2]*d;
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
        reachT[3]=p2pParam.trajTime[0]+p2pParam.trajTime[1]+(-p2pParam.maxVel+sqrt(v_square))/p2pParam.trajAcc[2];
        reachT[4]=p2pParam.trajTime[0]+p2pParam.trajTime[1]+(-p2pParam.maxVel-sqrt(v_square))/p2pParam.trajAcc[2];
    }

    bool isFound {false};
    double reachTime {-1};
    for(int i=0;i<5;i++)
    {
        if(reachT[i]>=limitedTimeInterval[0] && reachT[i]<limitedTimeInterval[1])//Time interval:[limitedTimeInterval[0],limitedTimeInterval[1])
        {
            isFound=true;
            reachTime=reachT[i];
            break;
        }
    }
    if(isFound==true)
    {
        printf("Return reachTime found: %f in interval [%f,%f]\n",reachTime,limitedTimeInterval[0],limitedTimeInterval[1]);
        return reachTime;
    }
    else
    {
        printf("Error! Conflict! Next point time cannot be calaulated, but it is judged to be here!\n");
        return 0;
    }
}

void RTOptimal::GetFinalNxtPntTime(PolylineParam *lnParam)
{
    //find optimal method to traverse next point in mostly 3*3*3=27 cases.
    printf("Next point is passed for %d,%d,%d times.\n",lnParam[0].passNxtPntNum,lnParam[1].passNxtPntNum,lnParam[2].passNxtPntNum);
    double maxTime {1e10};
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
                if(maxTime_tmp<maxTime)
                {
                    maxTime=maxTime_tmp;
                    nums[0]=i;
                    nums[1]=j;
                    nums[2]=k;
                }
            }
        }
    }

    printf("nums:%d,%d,%d\n",nums[0],nums[1],nums[2]);
    for(int i=0;i<3;i++)
    {
        lnParam[i].nxtPntTime=lnParam[i].passNxtPntTime[nums[i]];
    }
    double nxtPntTime_tmp[3] {lnParam[0].nxtPntTime,lnParam[1].nxtPntTime,lnParam[2].nxtPntTime};
    actScrewID=std::max_element(nxtPntTime_tmp,nxtPntTime_tmp+3)-nxtPntTime_tmp;
    nxtPntMinTime=nxtPntTime_tmp[actScrewID];
    printf("nxtPnttime:%.4f,%.4f,%.4f; nxtPntMinTime:%.4f\n",nxtPntTime_tmp[0],nxtPntTime_tmp[1],nxtPntTime_tmp[2],nxtPntMinTime);
}

bool RTOptimal::GetNxtPin(PolylineParam &lnParam, NextPointParam &nxtParam, int count, double &pIn)
{
    double curT=0.001*count-lstPntTime;
    if(nxtParam.optVel<=nxtParam.reachableVel[0])
    {
        if(curT<=nxtParam.minT[0])
        {
            pIn=lnParam.pnts[0]+lnParam.startV*curT+0.5*aLmt*curT*curT;
        }
        else if(curT<=nxtParam.minT[0]+nxtParam.minT[1])
        {
            pIn=lnParam.pnts[0]+lnParam.startV*nxtParam.minT[0]+0.5*aLmt*nxtParam.minT[0]*nxtParam.minT[0]
                    +nxtParam.min_vm*(curT-nxtParam.minT[0]);
        }
        else if(curT<=nxtParam.minT[0]+nxtParam.minT[1]+nxtParam.minT[2])
        {
            pIn=lnParam.pnts[0]+lnParam.startV*nxtParam.minT[0]+0.5*aLmt*nxtParam.minT[0]*nxtParam.minT[0]
                    +nxtParam.min_vm*(curT-nxtParam.minT[0])-0.5*aLmt*(curT-nxtParam.minT[0]-nxtParam.minT[1])*(curT-nxtParam.minT[0]-nxtParam.minT[1]);
        }
        else
        {
            printf("Current target point passed in less than 0.001s! Ignore the point and go on to next one!\n");
            for(int j=0;j<lnParam.pntsNum+1;j++)
            {
                lnParam.pnts[j]=lnParam.pnts[j+1];
            }
            lnParam.pntsNum--;
            lnParam.startV=nxtParam.reachableVel[0];
            isCurPntPassed=true;
            return false;
        }

        if(curT<=0.001)
        {
            nxtParam.record_t1=nxtParam.minT[0];
            nxtParam.record_t2=nxtParam.minT[2];
            printf("Unreachable NxtPnt time:%.4f,%.4f, acc:%.4f,%.4f\n",nxtParam.minT[0],nxtParam.minT[2],aLmt,-aLmt);
        }
        if(curT+0.001>nxtParam.minT[0]+nxtParam.minT[1]+nxtParam.minT[2])
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
    else if(nxtParam.optVel<=nxtParam.reachableVel[1])
    {
        double s=lnParam.pnts[1]-lnParam.pnts[0];
        double acc;
        if(nxtParam.optVel==lnParam.startV)
        {
            nxtParam.rchT[2]=nxtParam.rchT[0]=nxtPntMinTime/2;
            acc=2*(s/2-lnParam.startV*nxtParam.rchT[0])/(nxtParam.rchT[0]*nxtParam.rchT[0]);
            nxtParam.rch_vm=lnParam.startV+acc*nxtParam.rchT[0];
        }
        else
        {
            double delta=(lnParam.startV*nxtPntMinTime-s)*(lnParam.startV*nxtPntMinTime-s)/2
                        +(nxtParam.optVel*nxtPntMinTime-s)*(nxtParam.optVel*nxtPntMinTime-s)/2;
            double t1=(-(nxtParam.optVel*nxtPntMinTime-s)+sqrt(delta))/(lnParam.startV-nxtParam.optVel);
            double t2=(-(nxtParam.optVel*nxtPntMinTime-s)-sqrt(delta))/(lnParam.startV-nxtParam.optVel);
            if(t1>=0 && t1<=nxtPntMinTime)
            {
                nxtParam.rchT[0]=t1;
            }
            else if(t2>=0 && t2<=nxtPntMinTime)
            {
                nxtParam.rchT[0]=t2;
            }
            else
            {
                printf("Error! optVel cannot be reached! t1=%.4f & t2=%.4f are negative!\n",t1,t2);
            }
            acc=(nxtParam.optVel-lnParam.startV)/(2*nxtParam.rchT[0]-nxtPntMinTime);
            nxtParam.rchT[2]=nxtPntMinTime-nxtParam.rchT[0];
            nxtParam.rch_vm=lnParam.startV+acc*nxtParam.rchT[0];
        }

        if(curT<=nxtParam.rchT[0])
        {
            pIn=lnParam.pnts[0]+lnParam.startV*curT+0.5*acc*curT*curT;
        }
        else if(curT<=nxtParam.rchT[0]+nxtParam.rchT[2])
        {
            pIn=lnParam.pnts[0]+lnParam.startV*nxtParam.rchT[0]+0.5*acc*nxtParam.rchT[0]*nxtParam.rchT[0]
                    +nxtParam.rch_vm*(curT-nxtParam.rchT[0])-0.5*acc*(curT-nxtParam.rchT[0])*(curT-nxtParam.rchT[0]);
        }
        else
        {
            printf("Current target point passed in less than 0.001s! Ignore the point and go on to next one!\n");
            for(int j=0;j<lnParam.pntsNum+1;j++)
            {
                lnParam.pnts[j]=lnParam.pnts[j+1];
            }
            lnParam.pntsNum--;
            lnParam.startV=nxtParam.reachableVel[0];
            isCurPntPassed=true;
            return false;
        }

        if(curT<=0.001)
        {
            nxtParam.record_t1=nxtParam.rchT[0];
            nxtParam.record_t2=nxtParam.rchT[2];
            printf("Reachable NxtPnt time:%.4f,%.4f, acc:%.4f,%.4f\n",nxtParam.rchT[0],nxtParam.rchT[2],acc,-acc);
        }
        if(curT+0.001>nxtParam.rchT[0]+nxtParam.rchT[2])
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
        if(curT<=nxtParam.maxT[0])
        {
            pIn=lnParam.pnts[0]+lnParam.startV*curT-0.5*aLmt*curT*curT;
        }
        else if(curT<nxtParam.maxT[0]+nxtParam.maxT[1])
        {
            pIn=lnParam.pnts[0]+lnParam.startV*nxtParam.maxT[0]-0.5*aLmt*nxtParam.maxT[0]*nxtParam.maxT[0]
                    +nxtParam.max_vm*(curT-nxtParam.maxT[0]);
        }
        else if(curT<nxtParam.maxT[0]+nxtParam.maxT[1]+nxtParam.maxT[2])
        {
            pIn=lnParam.pnts[0]+lnParam.startV*nxtParam.maxT[0]-0.5*aLmt*nxtParam.maxT[0]*nxtParam.maxT[0]
                    +nxtParam.max_vm*(curT-nxtParam.maxT[0])+0.5*aLmt*(curT-nxtParam.maxT[0]-nxtParam.maxT[1])*(curT-nxtParam.maxT[0]-nxtParam.maxT[1]);
        }
        else
        {
            printf("Current target point passed in less than 0.001s! Ignore the point and go on to next one!\n");
            for(int j=0;j<lnParam.pntsNum+1;j++)
            {
                lnParam.pnts[j]=lnParam.pnts[j+1];
            }
            lnParam.pntsNum--;
            lnParam.startV=nxtParam.reachableVel[0];
            isCurPntPassed=true;
            return false;
        }

        if(curT<=0.001)
        {
            nxtParam.record_t1=nxtParam.maxT[0];
            nxtParam.record_t2=nxtParam.maxT[2];
            printf("UnReachable NxtPnt time:%.4f,%.4f, acc:%.4f,%.4f\n",nxtParam.maxT[0],nxtParam.maxT[2],-aLmt,aLmt);
        }
        if(curT+0.001>nxtParam.maxT[0]+nxtParam.maxT[2])
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

    return true;
}

void RTOptimal::GetNxtPntOptVel(PolylineParam &lnParam, NextPointParam &nxtParam)
{
    if(lnParam.pntsNum==2)
    {
        nxtParam.optVel=lnParam.endV;
    }
    else
    {
        if(lnParam.pnts[lnParam.fstBrkPntNum]>=lnParam.pnts[0])
        {
            nxtParam.optVel=sqrt(lnParam.fstBrkPntVel*lnParam.fstBrkPntVel+2*aLmt*(lnParam.pnts[lnParam.fstBrkPntNum]-lnParam.pnts[1]));
            if(nxtParam.optVel>vLmt)
            {
                nxtParam.optVel=vLmt;
            }
        }
        else
        {
            nxtParam.optVel=-sqrt(lnParam.fstBrkPntVel*lnParam.fstBrkPntVel-2*aLmt*(lnParam.pnts[lnParam.fstBrkPntNum]-lnParam.pnts[1]));
            if(nxtParam.optVel<-vLmt)
            {
                nxtParam.optVel=-vLmt;
            }
        }
    }
    printf("OptVel:%.4f\n",nxtParam.optVel);
}

void RTOptimal::GetNxtPntReachableVel(PolylineParam *lnParam, NextPointParam *nxtParam)
{
    int k=0;
    bool isNxtVelGot[3];
    while(k<5)
    {
        isNxtVelGot[0]=GetNxtPntParam(lnParam[0],nxtParam[0],0);
        if(isNxtVelGot[0]==false)
        {
            printf("MinTime for screw 0 modified! New nxtPntMinTime=%f.\n",nxtPntMinTime);
            k++;
            continue;
        }
        else
        {
            isNxtVelGot[1]=GetNxtPntParam(lnParam[1],nxtParam[1],1);
            if(isNxtVelGot[1]==false)
            {
                printf("MinTime for screw 1 modified! New nxtPntMinTime=%f.\n",nxtPntMinTime);
                k++;
                continue;
            }
            else
            {
                isNxtVelGot[2]=GetNxtPntParam(lnParam[2],nxtParam[2],2);
                if(isNxtVelGot[1]==false)
                {
                    printf("MinTime for screw 2 modified! New nxtPntMinTime=%f.\n",nxtPntMinTime);
                    k++;
                    continue;
                }
                else
                {
                    break;
                }
            }
        }
    }
    if(isNxtVelGot[0]==false || isNxtVelGot[1]==false || isNxtVelGot[2]==false || k==5)
    {
        printf("Error! Something wrong when calculating next point reachable velocity.k=%d\n",k);
    }
}

bool RTOptimal::GetNxtPntParam(PolylineParam &lnParam, NextPointParam &nxtParam, int screwID)
{
    double s=lnParam.pnts[1]-lnParam.pnts[0];
    double max_s=lnParam.startV*nxtPntMinTime+0.5*aLmt*nxtPntMinTime*nxtPntMinTime;
    double min_s=lnParam.startV*nxtPntMinTime-0.5*aLmt*nxtPntMinTime*nxtPntMinTime;

    //nxtPntMinTime in inoperative interval, new time need to be figured out
    if(s<0 && s>max_s)//startV<0,s<0,a>0
    {
        printf("startV=%f,s=%f,max_s=%.10f,min_s=%.10f\n",lnParam.startV,s,max_s,min_s);
        printf("Error! Screw %d cannot reach d=%.10f within the min time %.10f. max_d=%.10f\n",screwID,s,nxtPntMinTime,max_s);
        nxtPntMinTime=(-lnParam.startV-sqrt(lnParam.startV*lnParam.startV+2*aLmt*s))/aLmt;
        actScrewID=screwID;
        printf("Change new nxtPntMinTime to %f\n",nxtPntMinTime);
        return false;
    }
    else if(s>0 && s<min_s)//startV>0,s>0,a<0
    {
        printf("startV=%f,s=%f,max_s=%.10f,min_s=%.10f\n",lnParam.startV,s,max_s,min_s);
        printf("Error! Screw %d cannot reach d=%.10f within the min time %.10f. min_s=%.10f\n",screwID,s,nxtPntMinTime,min_s);
        nxtPntMinTime=(-lnParam.startV+sqrt(lnParam.startV*lnParam.startV-2*aLmt*s))/(-aLmt);
        printf("Change new nxtPntMinTime to %f\n",nxtPntMinTime);
        actScrewID=screwID;
        return false;
    }
    else
    {
        double v_square=aLmt*lnParam.startV*nxtPntMinTime+0.5*aLmt*aLmt*nxtPntMinTime*nxtPntMinTime-aLmt*s;
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
        nxtParam.minT[2]=sqrt(v_square)/aLmt;
        nxtParam.minT[1]=0;
        nxtParam.minT[0]=nxtPntMinTime-nxtParam.minT[2];
        nxtParam.min_vm=lnParam.startV+aLmt*nxtParam.minT[0];
        if(nxtParam.min_vm>vLmt)
        {
            nxtParam.min_vm=vLmt;
            nxtParam.minT[0]=(vLmt-lnParam.startV)/aLmt;
            double s_remain=s-(vLmt*vLmt-lnParam.startV*lnParam.startV)/(2*aLmt);
            nxtParam.minT[2]=2*(vLmt*(nxtPntMinTime-nxtParam.minT[0]-s_remain))/aLmt;
            nxtParam.minT[1]=nxtPntMinTime-nxtParam.minT[0]-nxtParam.minT[2];
        }
        nxtParam.reachableVel[0]=nxtParam.min_vm-aLmt*nxtParam.minT[2];

        v_square=-aLmt*lnParam.startV*nxtPntMinTime+0.5*aLmt*aLmt*nxtPntMinTime*nxtPntMinTime+aLmt*s;
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
        nxtParam.maxT[2]=sqrt(v_square)/aLmt;
        nxtParam.maxT[1]=0;
        nxtParam.maxT[0]=nxtPntMinTime-nxtParam.maxT[2];
        nxtParam.max_vm=lnParam.startV-aLmt*nxtParam.maxT[0];
        if(nxtParam.max_vm<-vLmt)
        {
            nxtParam.max_vm=-vLmt;
            nxtParam.maxT[0]=(-vLmt-lnParam.startV)/(-aLmt);
            double s_remain=s-(vLmt*vLmt-lnParam.startV*lnParam.startV)/(-2*aLmt);
            nxtParam.maxT[2]=2*(-vLmt*(nxtPntMinTime-nxtParam.maxT[0]-s_remain))/(-aLmt);
            nxtParam.maxT[1]=nxtPntMinTime-nxtParam.maxT[0]-nxtParam.maxT[2];
        }
        nxtParam.reachableVel[1]=nxtParam.max_vm+aLmt*nxtParam.maxT[2];

        printf("minT=%.4f,%.4f,%.4f, maxT=%.4f,%.4f,%.4f\n",
               nxtParam.minT[0],nxtParam.minT[1],nxtParam.minT[2],nxtParam.maxT[0],nxtParam.maxT[1],nxtParam.maxT[2]);
        printf("min_vm=%.4f, max_vm=%.4f\n",nxtParam.min_vm,nxtParam.max_vm);
        printf("reachableVel:%.4f,%.4f\n",nxtParam.reachableVel[0],nxtParam.reachableVel[1]);

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
                GetP2PMotionParamAcc(p2pParam,midV);
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
                    GetP2PMotionParamAcc(p2pParam,midV);
                    p2pParam.minTime=p2pParam.trajTime[0]+p2pParam.trajTime[1]+p2pParam.trajTime[2];
                }
                else//Type I -- acc-pos-dec
                {
                    p2pParam.inopInterNum=0;
                    p2pParam.velType=VelType::AccDec;
                    midV=sqrt(startV*startV/2+endV*endV/2+aLmt*(endP-startP));
                    GetP2PMotionParamAcc(p2pParam,midV);
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
                    GetP2PMotionParamAcc(p2pParam,midV);
                    p2pParam.minTime=p2pParam.trajTime[0]+p2pParam.trajTime[1]+p2pParam.trajTime[2];
                }
                else//Type I -- acc-pos-dec
                {
                    p2pParam.inopInterNum=0;
                    p2pParam.velType=VelType::AccDec;
                    midV=sqrt(startV*startV/2+endV*endV/2+aLmt*(endP-startP));
                    GetP2PMotionParamAcc(p2pParam,midV);
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
                    GetP2PMotionParamAcc(p2pParam,midV);
                    p2pParam.minTime=p2pParam.trajTime[0]+p2pParam.trajTime[1]+p2pParam.trajTime[2];
                }
                else if(endP-startP>=judgeDist2)//Type I -- acc-pos-dec
                {
                    p2pParam.inopInterNum=0;
                    p2pParam.velType=VelType::AccDec;
                    midV=sqrt(startV*startV/2+endV*endV/2+aLmt*(endP-startP));
                    GetP2PMotionParamAcc(p2pParam,midV);
                    p2pParam.minTime=p2pParam.trajTime[0]+p2pParam.trajTime[1]+p2pParam.trajTime[2];
                }
                else//All type
                {
                    p2pParam.inopInterNum=1;

                    //Type II -- dec-acc
                    p2pParam.velType=VelType::DecAcc;
                    midV=sqrt(startV*startV/2+endV*endV/2-aLmt*(endP-startP));
                    GetP2PMotionParamAcc(p2pParam,midV);
                    p2pParam.inopInter[0]=p2pParam.trajTime[0]+p2pParam.trajTime[1]+p2pParam.trajTime[2];

                    //Type III -- dec-neg-acc
                    p2pParam.velType=VelType::DecAcc;
                    midV=-sqrt(startV*startV/2+endV*endV/2-aLmt*(endP-startP));
                    GetP2PMotionParamAcc(p2pParam,midV);
                    p2pParam.inopInter[1]=p2pParam.trajTime[0]+p2pParam.trajTime[1]+p2pParam.trajTime[2];

                    //Type I -- acc-pos-dec
                    p2pParam.velType=VelType::AccDec;
                    midV=sqrt(startV*startV/2+endV*endV/2+aLmt*(endP-startP));
                    GetP2PMotionParamAcc(p2pParam,midV);
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

void RTOptimal::GetP2PMotionParamAcc(P2PMotionParam &p2pParam, double midV)
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

void RTOptimal::ScaleP2PMotionAcc(P2PMotionParam &p2pParam, double scaleTime)
{

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
