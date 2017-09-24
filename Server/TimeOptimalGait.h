#ifndef TIMEOPTIMALGAIT_H
#define TIMEOPTIMALGAIT_H

#include <aris.h>
#include <Robot_Type_I.h>
#include <sys/time.h>

class TimeOptimalGait
{
public:
    TimeOptimalGait();
    ~TimeOptimalGait();

    void GetStanceLegParam(int count, int legID, double s);
    void GetStanceDsBound(int count);
    void GetStanceSwitchPoint();
    void GetStanceOptimalDsBySwitchPoint();
    void GetStanceOptimalDsByDirectNI();
    void GetStanceOptimalDsByMinorIteration();

    void GetSwingLegParam(int count, int legID, double s);
    void GetSwingDsBound(int count, int legID);
    void GetSwingSwitchPoint(int legID);
    void GetSwingOptimalDsBySwitchPoint(int legID);
    void GetSwingOptimalDsByIteration(int legID);

    void GetOptimalDs();
    void GetOptimalGait2s();
    void GetOptimalGait2t();

    void OutputData();

private:
    double GetStanceSwitchMaxDec(int switchID, double ds);
    double GetStanceSwitchMinAcc(int switchID, double ds);
    double GetStanceSwitchDsBound(int switchID);
    void GetStanceTwoPointAtSwitch(double *lowPoint, double *upPoint);

    double GetStanceMaxDec(int count, double ds);
    double GetStanceMinAcc(int count, double ds);


    double GetSwingMaxDec(int count, double ds, int legID);
    double GetSwingMinAcc(int count, double ds, int legID);

    Robots::RobotTypeI rbt;

    const double vLmt {0.9};
    const double aLmt {3.2};
    const double delta_s {PI/900};
    double stepH {0.1};
    double stepD {0.8};
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
    const int stanceLegID[3] {1,3,5};


    double s_b[2251] {0};
    double b_sb[2251] {0};
    double db_sb[2251] {0};
    double ddb_sb[2251] {0};

    double Jvi[9] {0};
    double abs_param_dds[18] {0};//for vLmt of stanceLeg
    double param_a2[18] {0};//coefficient of ds*ds, param_dsds divided by param_dds
    double param_a1[18] {0};//coefficient of ds, param_ds divided by param_dds
    double param_a0L[18] {0};
    double param_a0H[18] {0};

    double output_PeeB[2251][18] {{0}};
    double output_dsds[2251][18] {{0}};
    double output_ds[2251][18] {{0}};
    double output_const[2251][18] {{0}};
    double output_const1[2251][18] {{0}};
    double output_const2[2251][18] {{0}};
    double output_dds[2251][18] {{0}};
    double output_a2[2251][18] {{0}};
    double output_a1[2251][18] {{0}};
    double output_a0L[2251][18] {{0}};
    double output_a0H[2251][18] {{0}};

    double ds_upBound_aLmt_body[2251] {0};
    double ds_upBound_vLmt_body[2251] {0};
    double dds_upBound_body[2251] {0};
    double dds_lowBound_body[2251] {0};
    double ds_upBound_body[2251] {0};

//    double slopedsBound_body[2251] {0};
//    double paramdds0Point_body[2251][18] {{0}};
//    double tangentPoint_body[2251] {0};
//    int paramdds0Count_body[18] {0};
//    int tangentCount_body {0};

    double slopeDelta_body[2251] {0};
    int isParamddsExact0_body[2251][18] {{0}};
    double switchPoint_body[2251] {0};
    char switchPointType_body[2251] {0};
    int switchScrewID_body[2251] {0};
    int switchCount_body {0};

    double real_ds_body[2251] {0};
    double real_dds_body[2251] {0};
    double ds_forward_body[2251] {0};
    double ds_backward_body[2251] {0};
    double dds_forward_body[2251] {0};
    double dds_backward_body[2251] {0};

    double timeArray_body[2251] {0};
    double Tstep {1};
    double s_b1 {0};
    double s_b2 {0};
    double pb_sw[2251] {0};
    double vb_sw[2251] {0};
    double ab_sw[2251] {0};
    double pb_sw_tmp[2251] {0};
    double vb_sw_tmp[2251] {0};
    double ab_sw_tmp[2251] {0};
    double pva_b[2251][3] {{0}};

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

//    double slopedsBound[901][6] {{0}};
//    double paramdds0Point[901][18] {{0}};
//    double tangentPoint[901][6] {{0}};
//    int paramdds0Count[18] {0};
//    int tangentCount[6] {0};

    double slopeDelta[901][6] {{0}};
    int isParamddsExact0[901][18] {{0}};
    double switchPoint[901][6] {{0}};
    int switchCount[6] {0};

    bool quitSwitchPoint {false};
    bool stopFlag {false};
    bool accFlag {true};//true acc, false dec
    int dec_start {0};
    int dec_end {0};
    unsigned int cycleCount {0};//used to limit the cycleCount of forward integration
    int stop_back {0};
    int ki_back {0};
    int ki_for {0};

    double ds_forward[901][6] {{0}};
    double ds_backward[901][6] {{0}};
    double dds_forward[901][6] {{0}};
    double dds_backward[901][6] {{0}};
    double real_ds[901][6] {{0}};
    double real_dds[901][6] {{0}};
//    double real_ddsMax[901][6] {{0}};
//    double real_ddsMin[901][6] {{0}};
    double real_ds_tmp[901][6] {{0}};

    double timeArray[901][6] {{0}};
    double timeArray_tmp[901][6] {{0}};
    double totalTime[6] {0};
    double maxTime {0};
    int maxTotalCount {0};
    int maxTotalCount_last {1};
    int maxTimeID {0};

//    double output_Pee[901][18] {{0}};
//    double output_Pin[901][18] {{0}};
//    double output_Vin[901][18] {{0}};
//    double output_Ain[901][18] {{0}};
};

#endif