#include "NormalGait.h"

/*
<mwr default="mwr_param">
    <mwr_param type="group">
        <totalCount abbreviation="t" type="int" default="1000"/>
        <u abbreviation="u" type="double" default="0"/>
        <v abbreviation="v" type="double" default="0"/>
        <w abbreviation="w" type="double" default="0"/>
        <pitch abbreviation="p" type="double" default="0"/>
        <yaw abbreviation="y" type="double" default="0"/>
        <roll abbreviation="r" type="double" default="0"/>
    </mwr_param>
</mwr>
<arc default="arc_param">
    <arc_param type="group">
        <leg_param type="unique" default="all">
           <all abbreviation="a"/>
           <leg abbreviation="l" type="int" default="0"/>
        </leg_param>
        <distance abbreviation="d" type="double" default="0"/>
        <totalCount abbreviation="t" type="int" default="1000"/>
    </arc_param>
</arc>
<cwk default="cwk_param">
    <cwk_param type="group">
        <startangle abbreviation="a" type="double" default="0"/>
        <radius abbreviation="r" type="double" default="1"/>
        <totalCount abbreviation="t" type="int" default="1000"/>
        <n abbreviation="n" type="int" default="1"/>
        <h abbreviation="h" type="double" default="0.05"/>
        <beta abbreviation="b" type="double" default="0.05"/>
        <direction abbreviation="d" type="int" default="1"/>
    </cwk_param>
</cwk>
<ppw default="ppw_param">
    <ppw_param type="group">
        <totalCount abbreviation="t" type="int" default="1000"/>
        <height abbreviation="h" type="double" default="0.05"/>
        <x abbreviation="x" type="double" default="0"/>
        <y abbreviation="y" type="double" default="0"/>
        <z abbreviation="z" type="double" default="0"/>
    </ppw_param>
</ppw>


rs.addCmd("mwr",NormalGait::parseMoveWithRotate,NormalGait::moveWithRotate);
rs.addCmd("arc",NormalGait::parseAdjustRc,NormalGait::adjustRc);
rs.addCmd("cwk",NormalGait::parseCircleWalk,NormalGait::circleWalk);
rs.addCmd("ppw",NormalGait::parseP2PWalk,NormalGait::p2pWalk);
*/

namespace NormalGait
{
    void parseMoveWithRotate(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
    {
        MoveRotateParam param;

        for(auto &i:params)
        {
            if(i.first=="u")
            {
                param.targetBodyPE213[0]=stod(i.second);
            }
            else if(i.first=="v")
            {
                param.targetBodyPE213[1]=stod(i.second);
            }
            else if(i.first=="w")
            {
                param.targetBodyPE213[2]=stod(i.second);
            }
            else if(i.first=="yaw")
            {
                param.targetBodyPE213[3]=stod(i.second)*PI/180;
            }
            else if(i.first=="pitch")
            {
                param.targetBodyPE213[4]=stod(i.second)*PI/180;
            }
            else if(i.first=="roll")
            {
                param.targetBodyPE213[5]=stod(i.second)*PI/180;
            }
            else if(i.first=="totalCount")
            {
                param.totalCount=stoi(i.second);
            }
            else
            {
                std::cout<<"parse failed"<<std::endl;
            }
        }

        msg.copyStruct(param);

        std::cout<<"finished parse"<<std::endl;
    }

    int moveWithRotate(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
    {
        auto &robot = static_cast<Robots::RobotBase &>(model);
        auto &param = static_cast<const MoveRotateParam &>(param_in);

        /*****test IMU*****/
        double peImuGrnd2BodyGrnd[6] {0,0,0,PI/6,-PI/2,0};
        double pmImuGrnd2BodyGrnd[16];
        double pmBody2Imu[16];
        double peImu2ImuGrnd[6] {0};
        peImu2ImuGrnd[3]=param.imu_data->at(0).euler[2];
        peImu2ImuGrnd[4]=param.imu_data->at(0).euler[1];
        peImu2ImuGrnd[5]=param.imu_data->at(0).euler[0];
        double pmImu2ImuGrnd[16];
        double bodyPe[6];//213
        double bodyPm[16];
        aris::dynamic::s_pe2pm(peImuGrnd2BodyGrnd,pmImuGrnd2BodyGrnd,"213");
        aris::dynamic::s_pe2pm(peImu2ImuGrnd,pmImu2ImuGrnd,"321");
        aris::dynamic::s_inv_pm(pmImuGrnd2BodyGrnd,pmBody2Imu);
        aris::dynamic::s_pm_dot_pm(pmImuGrnd2BodyGrnd,pmImu2ImuGrnd,pmBody2Imu,bodyPm);
        aris::dynamic::s_pm2pe(bodyPm,bodyPe,"213");
        //printf("bodyPe:%f,%f,%f\n",bodyPe[3],bodyPe[4],bodyPe[5]);
        /*******************/

        static double beginBodyPE213[6];
        static double pEE[18];
        if(param.count==0)
        {
            robot.GetPeb(beginBodyPE213,"213");
            robot.GetPee(pEE);
        }

        double realBodyPE213[6];
        for(int i=0;i<6;i++)
        {
            double s = -(param.targetBodyPE213[i] / 2)*cos(PI * (param.count + 1) / param.totalCount ) + param.targetBodyPE213[i] / 2;
            realBodyPE213[i]=beginBodyPE213[i]+s; //target of current ms
        }

        robot.SetPeb(realBodyPE213,"213");
        robot.SetPee(pEE);

        return param.totalCount - param.count - 1;
    }

    void parseAdjustRc(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
    {
        AdjustRcParam param;

        for(auto &i:params)
        {
            if(i.first=="distance")
            {
                param.distance=std::stod(i.second);
            }
            else if(i.first=="all")
            {
                param.isAll=true;
            }
            else if(i.first=="leg")
            {
                param.legID=std::stoi(i.second);
		std::fill_n(param.actLegs,6,false);
		param.actLegs[param.legID]=true;
            }
            else if(i.first=="totalCount")
            {
                param.totalCount=std::stoi(i.second);
            }
            else
            {
                std::cout<<"parse failed"<<std::endl;
            }
        }

        msg.copyStruct(param);

        std::cout<<"finished parse"<<std::endl;
    }

    int adjustRc(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
    {
        auto &robot = static_cast<Robots::RobotBase &>(model);
        auto &param = static_cast<const AdjustRcParam &>(param_in);

        static double StartPeb[6];
        static double StartPee[18];
        double currentPee[18];
        if(param.count==0)
        {
            robot.GetPeb(StartPeb);
            robot.GetPee(StartPee);
        }

        if(param.isAll==true)
        {
            for(int i=0;i<3;i++)//0 2 4
            {
                if(param.count<param.totalCount)
                {
                    currentPee[6*i]=StartPee[6*i]+param.distance/2*sin(PI/6+2*i*PI/3)*(1-cos(param.count*PI/param.totalCount));
                    currentPee[6*i+1]=StartPee[6*i+1]+0.025*(1-cos(param.count*2*PI/param.totalCount));
                    currentPee[6*i+2]=StartPee[6*i+2]+param.distance/2*cos(PI/6+2*i*PI/3)*(1-cos(param.count*PI/param.totalCount));

                    currentPee[6*i+3]=StartPee[6*i+3];
                    currentPee[6*i+4]=StartPee[6*i+4];
                    currentPee[6*i+5]=StartPee[6*i+5];
                }
                else
                {
                    currentPee[6*i]=StartPee[6*i]+param.distance*sin(PI/6+2*i*PI/3);
                    currentPee[6*i+1]=StartPee[6*i+1];
                    currentPee[6*i+2]=StartPee[6*i+2]+param.distance*cos(PI/6+2*i*PI/3);

                    currentPee[6*i+3]=StartPee[6*i+3]+param.distance/2*sin(PI/2-2*i*PI/3)*(1-cos((param.count-param.totalCount)*PI/param.totalCount));
                    currentPee[6*i+4]=StartPee[6*i+4]+0.025*(1-cos((param.count-param.totalCount)*2*PI/param.totalCount));
                    currentPee[6*i+5]=StartPee[6*i+5]+param.distance/2*cos(PI/2-2*i*PI/3)*(1-cos((param.count-param.totalCount)*PI/param.totalCount));
                }
            }
        }
        else
        {
            double k[6] {PI/6,PI/2,5*PI/6,-PI/6,3*PI/2,-5*PI/6};
            for(int i=0;i<6;i++)
            {
                if(param.actLegs[i]==true)
                {
                    currentPee[3*i]=StartPee[3*i]+param.distance/2*sin(k[i])*(1-cos(param.count*PI/param.totalCount));
                    currentPee[3*i+1]=StartPee[3*i+1]+0.025*(1-cos(param.count*2*PI/param.totalCount));
                    currentPee[3*i+2]=StartPee[3*i+2]+param.distance/2*cos(k[i])*(1-cos(param.count*PI/param.totalCount));
                }
                else
                {
                    currentPee[3*i]=StartPee[3*i];
                    currentPee[3*i+1]=StartPee[3*i+1];
                    currentPee[3*i+2]=StartPee[3*i+2];
                }
            }
        }

        robot.SetPeb(StartPeb);
        robot.SetPee(currentPee);
        if(param.isAll==true)
        {
            return 2*param.totalCount - param.count - 1;
        }
        else
        {
            return param.totalCount - param.count - 1;
        }
    }

    void parseP2PWalk(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
    {
        P2PWalkParam param;

        for(auto &i:params)
        {
            if(i.first=="x")
            {
                param.targetPointInB[0]=std::stod(i.second);
            }
            else if(i.first=="y")
            {
                param.targetPointInB[1]=std::stod(i.second);
            }
            else if(i.first=="z")
            {
                param.targetPointInB[2]=std::stod(i.second);
            }
            else if(i.first=="totalCount")
            {
                param.totalCount=std::stoi(i.second);
            }
            else if(i.first=="height")
            {
                param.height=std::stod(i.second);
            }
            else
            {
                std::cout<<"parse failed"<<std::endl;
            }
        }

        msg.copyStruct(param);

        std::cout<<"finished parse"<<std::endl;
    }

    int p2pWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
    {
        auto &robot = static_cast<Robots::RobotBase &>(model);
        auto &param = static_cast<const P2PWalkParam &>(param_in);

        static double beginPeb[6] {0};
        if(param.count==0)
        {
            robot.GetPeb(beginPeb);
        }
        double Peb[6] {0};
        double dist = sqrt(param.targetPointInB[0]*param.targetPointInB[0]+param.targetPointInB[2]*param.targetPointInB[2]);

        Robots::WalkParam walkParam;
        walkParam.alpha=atan2(param.targetPointInB[0],param.targetPointInB[2])+PI;
        walkParam.beta=0;
        walkParam.h=param.height;
        if(dist<=0.2)
        {
            walkParam.n=1;
        }
        else
        {
            walkParam.n=(int)((dist-0.2)/0.4)+2;
        }
        walkParam.d=dist/(walkParam.n-0.5);
        walkParam.count=param.count;
        walkParam.totalCount=param.totalCount;

        Robots::walkGait(robot,walkParam);
        robot.GetPeb(Peb);
        Peb[1]=beginPeb[1]+param.targetPointInB[1]*(-cos(param.count/(2*walkParam.n*param.totalCount)*PI)+1)/2;

        robot.SetPeb(Peb);

        return 2 * walkParam.n * param.totalCount - param.count - 1;
    }

    void parseCircleWalk(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg)
    {
        CircleWalkParam param;

        for(auto &i:params)
        {
            if(i.first=="startangle")
            {
                param.startAngle=std::stod(i.second);
            }
            else if(i.first=="radius")
            {
                param.radius=std::stod(i.second);
            }
            else if(i.first=="totalCount")
            {
                param.totalCount=std::stoi(i.second);
            }
            else if(i.first=="n")
            {
                param.n=std::stoi(i.second);
            }
            else if(i.first=="h")
            {
                param.h=std::stod(i.second);
            }
            else if(i.first=="beta")
            {
                param.beta=std::stod(i.second);
            }
            else if(i.first=="direction")
            {
                param.direction=std::stoi(i.second);
            }
            else
            {
                std::cout<<"parse failed"<<std::endl;
            }
        }

        msg.copyStruct(param);

        std::cout<<"finished parse"<<std::endl;
    }

    int circleWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in)
    {
        auto &robot = static_cast<Robots::RobotBase &>(model);
        auto &param = static_cast<const CircleWalkParam &>(param_in);

        static aris::dynamic::FloatMarker beginMak{ robot.ground() };
        static double beginPee[18];

        if (param.count%param.totalCount == 0)
        {
            beginMak.setPrtPm(*robot.body().pm());
            beginMak.update();
            robot.GetPee(beginPee, beginMak);
        }

        double Peb[6] {0};
        double Pee[18] {0};
        memcpy(Pee,beginPee,18*sizeof(double));

        int period_count = param.count%param.totalCount;
        int half_period = param.count/param.totalCount;
        int period = param.count/(2*param.totalCount);
        double bodyStartAngle=param.startAngle+(0.5*half_period-0.25)*param.beta*param.direction;
        double legStartAngle0=param.startAngle+(period-0.5)*param.beta*param.direction;
        double legStartAngle1=param.startAngle+period*param.beta*param.direction;

        if(param.count<param.totalCount)//acc
        {
            Peb[0]=param.radius*(sin(param.startAngle)-sin(param.startAngle+0.25*param.direction*param.beta*param.count*param.count/param.totalCount/param.totalCount));
            Peb[2]=param.radius*(cos(param.startAngle)-cos(param.startAngle+0.25*param.direction*param.beta*param.count*param.count/param.totalCount/param.totalCount));

            for(int i=0;i<3;i++)
            {
                if(param.direction==-1)
                {
                    Pee[6*i]=beginPee[6*i]+param.radius*(sin(param.startAngle)-sin(param.startAngle+0.5*param.direction*param.beta))/2*(1-cos(param.count*PI/param.totalCount));
                    Pee[6*i+1]=beginPee[6*i+1]+0.05*(1-cos(param.count*2*PI/param.totalCount));
                    Pee[6*i+2]=beginPee[6*i+2]+param.radius*(cos(param.startAngle)-cos(param.startAngle+0.5*param.direction*param.beta))/2*(1-cos(param.count*PI/param.totalCount));
                }
                else
                {
                    Pee[6*i+3]=beginPee[6*i+3]+param.radius*(sin(param.startAngle)-sin(param.startAngle+0.5*param.direction*param.beta))/2*(1-cos(param.count*PI/param.totalCount));
                    Pee[6*i+4]=beginPee[6*i+4]+0.05*(1-cos(param.count*2*PI/param.totalCount));
                    Pee[6*i+5]=beginPee[6*i+5]+param.radius*(cos(param.startAngle)-cos(param.startAngle+0.5*param.direction*param.beta))/2*(1-cos(param.count*PI/param.totalCount));
                }
            }
        }
        else if(param.count<(2*param.n-1)*param.totalCount)//const
        {
            Peb[0]=param.radius*(sin(bodyStartAngle)-sin(bodyStartAngle+0.5*param.direction*param.beta*period_count/param.totalCount));
            Peb[2]=param.radius*(cos(bodyStartAngle)-cos(bodyStartAngle+0.5*param.direction*param.beta*period_count/param.totalCount));

            for(int i=0;i<3;i++)
            {
                if(half_period%2==1)
                {
                    if(param.direction==-1)
                    {
                        Pee[6*i+3]=beginPee[6*i+3]+param.radius*(sin(legStartAngle1)-sin(legStartAngle1+param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                        Pee[6*i+4]=beginPee[6*i+4]+0.05*(1-cos(period_count*2*PI/param.totalCount));
                        Pee[6*i+5]=beginPee[6*i+5]+param.radius*(cos(legStartAngle1)-cos(legStartAngle1+param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                    }
                    else
                    {
                        Pee[6*i]=beginPee[6*i]+param.radius*(sin(legStartAngle1)-sin(legStartAngle1+param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                        Pee[6*i+1]=beginPee[6*i+1]+0.05*(1-cos(period_count*2*PI/param.totalCount));
                        Pee[6*i+2]=beginPee[6*i+2]+param.radius*(cos(legStartAngle1)-cos(legStartAngle1+param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                    }
                }
                else
                {
                    if(param.direction==-1)
                    {
                        Pee[6*i]=beginPee[6*i]+param.radius*(sin(legStartAngle0)-sin(legStartAngle0+param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                        Pee[6*i+1]=beginPee[6*i+1]+0.05*(1-cos(period_count*2*PI/param.totalCount));
                        Pee[6*i+2]=beginPee[6*i+2]+param.radius*(cos(legStartAngle0)-cos(legStartAngle0+param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                    }
                    else
                    {
                        Pee[6*i+3]=beginPee[6*i+3]+param.radius*(sin(legStartAngle0)-sin(legStartAngle0+param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                        Pee[6*i+4]=beginPee[6*i+4]+0.05*(1-cos(period_count*2*PI/param.totalCount));
                        Pee[6*i+5]=beginPee[6*i+5]+param.radius*(cos(legStartAngle0)-cos(legStartAngle0+param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                    }
                }
            }
        }
        else//dec
        {
            Peb[0]=param.radius*(sin(bodyStartAngle)-sin(bodyStartAngle+0.5*param.direction*param.beta*period_count/param.totalCount-0.25*param.direction*param.beta*param.count*param.count/param.totalCount/param.totalCount));
            Peb[2]=param.radius*(cos(bodyStartAngle)-cos(bodyStartAngle+0.5*param.direction*param.beta*period_count/param.totalCount-0.25*param.direction*param.beta*param.count*param.count/param.totalCount/param.totalCount));

            for(int i=0;i<3;i++)
            {
                if(param.direction==-1)
                {
                    Pee[6*i+3]=beginPee[6*i+3]+param.radius*(sin(legStartAngle1)-sin(legStartAngle1+0.5*param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                    Pee[6*i+4]=beginPee[6*i+4]+0.05*(1-cos(period_count*2*PI/param.totalCount));
                    Pee[6*i+5]=beginPee[6*i+5]+param.radius*(cos(legStartAngle1)-cos(legStartAngle1+0.5*param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                }
                else
                {
                    Pee[6*i]=beginPee[6*i]+param.radius*(sin(legStartAngle1)-sin(legStartAngle1+0.5*param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                    Pee[6*i+1]=beginPee[6*i+1]+0.05*(1-cos(period_count*2*PI/param.totalCount));
                    Pee[6*i+2]=beginPee[6*i+2]+param.radius*(cos(legStartAngle1)-cos(legStartAngle1+0.5*param.direction*param.beta))/2*(1-cos(period_count*PI/param.totalCount));
                }
            }
        }
        return param.count-2*param.n*param.totalCount-1;
    }
}
