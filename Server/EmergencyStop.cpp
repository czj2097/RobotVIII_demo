#include "EmergencyStop.h"

bool EmergencyStop::needRecord {true};
double EmergencyStop::sndLastPin[18] {0};
double EmergencyStop::lastPin[18];
double EmergencyStop::lastVin[18];
int EmergencyStop::startCount {0};
int EmergencyStop::totalCount {0};

EmergencyStop::EmergencyStop(){}
EmergencyStop::~EmergencyStop(){}

void EmergencyStop::doRecord(Robots::RobotBase &robot)
{
    if(needRecord==true)
    {
        robot.GetPin(lastPin);
        for(int i=0;i<18;i++)
        {
            lastVin[i]=(lastPin[i]-sndLastPin[i])*1000;
        }
        memcpy(sndLastPin,lastPin,18*sizeof(double));
    }
}

int EmergencyStop::stop(Robots::RobotBase &robot, int currentCount, double maxAin)
{
    double nextPin[18] {0};
    double Ain[18] {0};

    if(needRecord==true)
    {
        double absLastVin[18] {0};
        double maxAbsLastVin {0};
        for(int i=0;i<18;i++)
        {
            absLastVin[i]=fabs(lastVin[i]);
        }
        maxAbsLastVin=*std::max_element(absLastVin,absLastVin+18);
        totalCount=(int)(maxAbsLastVin/fabs(maxAin)*1000)+1;
        startCount=currentCount-1;
    }
    needRecord=false;

    for(int i=0;i<18;i++)
    {
        Ain[i]=-lastVin[i]/(0.001*totalCount);
        nextPin[i]=lastPin[i]+lastVin[i]*(currentCount-startCount)*1e-3+0.5*Ain[i]*(currentCount-startCount)*(currentCount-startCount)*1e-6;
    }
    robot.SetPin(nextPin);

    return totalCount - (currentCount-startCount) - 1;
}
