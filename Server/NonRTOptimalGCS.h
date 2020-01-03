#ifndef NONRT_OPTIMAL_GCS_H
#define NONRT_OPTIMAL_GCS_H
#include <sys/time.h>

#include <aris_dynamic.h>
#include <Robot_Type_I.h>
#include <Robot_Gait.h>

using namespace aris::control;

namespace TimeOptimal
{
    class NonRTOptimalGCS
    {
    public:
        void GetTimeOptimalGait(double step_length, double step_height, double &duty_cycle, double acc_limit, double vel_limit, double *init_tippos,
                                double *out_tippos, double &out_bodyvel, double &out_period);

        void GetStanceLegParam(int count, int legID, double s);
        void GetStanceDsBound(int count);
        void GetStanceDsBoundByNewton(int count);
        void GetStanceDsBoundByFunc(int count);
        void GetStanceSwitchPoint();
        void GetStanceOptimalDsBySwitchPoint();
        void GetStanceOptimalDsByDirectNI();
        void GetStanceOptimalDsAtSb(double s_b1, double s_b2, double s_b3, double s_b4);
        void GetStanceOptimalDsByMinorIteration();
        void GetStanceOptimalDsByNewton();
        void AnalyseStanceOptimalTime();

        void GetSwingLegParam(int count, int legID, double sw, double *pva_body);
        void GetSwingDsBound(int count, int legID);
        void GetSwingDsBoundByNewton(int count, int legID);
        void GetSwingDsBoundByFunc(int count, int legID);
        void GetSwingSwitchPoint(int legID);
        void GetSwingOptimalDsBySwitchPoint(int legID);
        void GetSwingOptimalDsByDirectNI(int legID);

        void GetOptimalDsByMajorIteration();
        void GetOptimalGait2s();
        void GetOptimalGait2t();
        void GetOptimalGait2t(double *out_tippos, double &out_bodyvel, double &out_period);

        void OutputData();

        void GetNormalGait();
        void GetEntireGait();

    private:
        double GetStanceSwitchMaxDec(int switchID, double ds);
        double GetStanceSwitchMinAcc(int switchID, double ds);
        double GetStanceSwitchDsBound(int switchID);
        void GetStanceTwoPointAtSwitch(double *lowPoint, double *upPoint);
        double GetStanceMaxDec(int count, double ds);
        double GetStanceMinAcc(int count, double ds);
        void ApplyStanceExtraItegration();

        double GetSwingSwitchMaxDec(int switchID, double ds, int legID);
        double GetSwingSwitchMinAcc(int switchID, double ds, int legID);
        double GetSwingSwitchDsBound(int switchID, int legID, double *pva_body);
        void GetSwingTwoPointAtSwitch(int legID, double *lowPoint, double *upPoint);
        double GetSwingMaxDec(int count, double ds, int legID);
        double GetSwingMinAcc(int count, double ds, int legID);
        void ApplySwingExtraIntegration(int legID);

        Robots::RobotTypeI rbt;

        double vLmt {1.0};
        double aLmt {5.0};
        double stepH {0.05};
        double stepD {0.4};
        double dutyCycle {0.6};
        double initPeb[6] {0};
    //    double initPee[18] {  -0.3, -0.85, -0.65,
    //                         -0.45, -0.85,     0,
    //                          -0.3, -0.85,  0.65,
    //                           0.3, -0.85, -0.65,
    //                          0.45, -0.85,     0,
    //                           0.3, -0.85,  0.65 };//024 swing, 135 stance
        double initPee[18] { -0.30, -0.58, -0.52,
                             -0.60, -0.58,  0,
                             -0.30, -0.58,  0.52,
                              0.30, -0.58, -0.52,
                              0.60, -0.58,  0,
                              0.30, -0.58,  0.52 };

        const int swingSCount {900};
        const int count1 {100};
        const int count2 {1000};
        const int count3 {1200};
        const int count4 {2100};
        const int count5 {2200};


        double s_b[2201] {0};
        double b_sb[2201] {0};
        double db_sb[2201] {0};
        double ddb_sb[2201] {0};

        double Jvi[9] {0};
        double abs_param_dds[18] {0};//for vLmt of stanceLeg
        double param_dds[18] {0};
        double param_a2[18] {0};//coefficient of ds*ds, param_dsds divided by param_dds
        double param_a1[18] {0};//coefficient of ds, param_ds divided by param_dds
        double param_a0L[18] {0};
        double param_a0H[18] {0};

        double output_PeeB[2201][18] {{0}};
        double output_dsds[2201][18] {{0}};
        double output_ds[2201][18] {{0}};
        double output_const[2201][18] {{0}};
        double output_const1[2201][18] {{0}};
        double output_const2[2201][18] {{0}};
        double output_dds[2201][18] {{0}};
        double output_a2[2201][18] {{0}};
        double output_a1[2201][18] {{0}};
        double output_a0L[2201][18] {{0}};
        double output_a0H[2201][18] {{0}};

        double ds_upBound_aLmt_body[2201] {0};
        double ds_upBound_vLmt_body[2201] {0};
        double dds_upBound_body[2201] {0};
        double dds_lowBound_body[2201] {0};
        double ds_upBound_body[2201] {0};

        double slopeDelta_body[2201] {0};
        int isParamddsExact0_body[2201][18] {{0}};
        double switchPoint_body[2201] {0};
        char switchType_body[2201] {0};
        int switchScrewID_body[2201] {0};
        int switchCount_body {0};

        double real_ds_body[2201] {0};
        double real_dds_body[2201] {0};
        double real_ddsMax_body[2201] {0};
        double real_ddsMin_body[2201] {0};
        double ds_forward_body[2201] {0};
        double ds_backward_body[2201] {0};
        double dds_forward_body[2201] {0};
        double dds_backward_body[2201] {0};

        double timeArray_body[2201] {0};
        double Tstep {1};
        double s_b1 {0};
        double s_b2 {0};
        double s_b3 {0};
        double s_b4 {0};

        double output_init_pva_b[2201][3] {{0}};
        double pva_b[2201][3] {{0}};

        int swCount {0};
        double s_w[901][6] {{0}};
        double f_sw[3] {0};
        double df_sw[3] {0};
        double ddf_sw[3] {0};

        double ds_upBound_aLmt[901][6] {{0}};
        double ds_lowBound_aLmt[901][6] {{0}};
        double ds_upBound_vLmt[901][6] {{0}};
        double ds_upBound[901][6] {{0}};
        double ds_lowBound[901][6] {{0}};
        double dds_lowBound[901][6] {{0}};
        double dds_upBound[901][6] {{0}};

        double slopeDelta[901][6] {{0}};
        int isParamddsExact0[901][18] {{0}};
        double switchPoint[901][6] {{0}};
        char switchType[901][6] {{0}};
        int switchScrewID[901][6] {{0}};
        int switchCount[6] {0};

        double ds_forward[901][6] {{0}};
        double ds_backward[901][6] {{0}};
        double dds_forward[901][6] {{0}};
        double dds_backward[901][6] {{0}};
        double real_ds[901][6] {{0}};
        double real_dds[901][6] {{0}};
        double real_ddsMax[901][6] {{0}};
        double real_ddsMin[901][6] {{0}};
        double real_ds_scale[901][6] {{0}};
        double real_dds_scale[901][6] {{0}};

        double timeArray[901][6] {{0}};
        double maxSwingTime {0};
        int swingTimeCount {0};
        int stepTimeCount {0};
        double timeArray_tmp[901][6] {{0}};
        double timeArray_body_tmp[2201] {0};

        double Pee_s[2201][18] {{0}};
        double Pin_s[2201][18] {{0}};
        double Vin_s[2201][18] {{0}};
        double Ain_s[2201][18] {{0}};
    };

    struct FastWalkGCSParam
    {
        double Peb;
        double Pee[18];
        double Pin[18];
    };
    class FastWalkGCS
    {
    public:
        static int  fastWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);
        static void parseFastWalk(const std::string &cmd, const std::map<std::string, std::string> &params, aris::core::Msg &msg);
        static void GetScalingPath(Robots::WalkParam &param);
        static void recordData();

    private:
        static double initBodyPos[6];
        static double initPee[18];
        static double swingPee_scale[3001][18];
        static double bodyConstVel;
        static double dutyCycle;

        static FastWalkGCSParam fwgParam;
        static Pipe<FastWalkGCSParam> fwgPipe;
        static std::thread fwgThread;
    };
}

#endif
