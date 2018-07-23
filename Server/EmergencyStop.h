#ifndef EMERGENCYSTOP_H
#define EMERGENCYSTOP_H

#include <Robot_Type_I.h>

class EmergencyStop
{
public:
    EmergencyStop();
    ~EmergencyStop();
    void doRecord(Robots::RobotBase &robot);
    int stop(Robots::RobotBase &robot, int currentCount, double maxAin);

private:
    static bool needRecord;
    static double sndLastPin[18];
    static double lastPin[18];
    static double lastVin[18];
    static int startCount;
    static int totalCount;
};

/*
//init
EmergencyStop stoper;
stoper.doRecord(robot);

//use
int ret=stoper.stop(robot,param.count,3.2);
return ret;
*/

#endif
