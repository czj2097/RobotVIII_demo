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

    void GetStanceLegParam(int count, int legID);
    void GetSwingLegParam(int count, int legID);
    void GetStanceDsBound(int count);
    void GetSwingDsBound(int count, int legID);
    void GetStanceSwitchPoint();
    void GetSwingSwitchPoint(int legID);
    void GetStanceOptimalDsBySwitchPoint();
    void GetSwingOptimalDsBySwitchPoint(int legID);
    void GetStanceOptimalDsByIteration();
    void GetSwingOptimalDsByIteration(int legID);

    void GetOptimalDs();
    void GetOptimalGait2s();
    void GetOptimalGait2t();

    void OutputData();

private:
    Robots::RobotTypeI rbt;

    const int totalCount {1800};
    const double vLmt {0.9};
    const double aLmt {3.2};
    const double delta_s {PI/1800};
    const int stanceLegID[3] {1,3,5};

    double stepH {0.1};
    double stepD {0.8};
    double initPeb[6] {0};
    double initPee[18] {  -0.3, -0.85, -0.65,
                         -0.45, -0.85,     0,
                          -0.3, -0.85,  0.65,
                           0.3, -0.85, -0.65,
                          0.45, -0.85,     0,
                           0.3, -0.85,  0.65 };//024 swing, 135 stance

    double s_t {0};
    double s_w {0};
    double b_st[1801] {0};
    double db_st[1801] {0};
    double ddb_st[1801] {0};
    double pb_sw[1801] {0};
    double vb_sw[1801] {0};
    double ab_sw[1801] {0};
    double f_sw[3] {0};
    double df_sw[3] {0};
    double ddf_sw[3] {0};

    double Jvi[9] {0};
    double dJvi_x[9] {0};
    double dJvi_y[9] {0};
    double dJvi_z[9] {0};
    double dJvi[9] {0};//used in stanceLeg
    double dJvi_dot_f[9] {0};//used in swingLeg
    double dJvi_dot_vb[9] {0};//used in swingLeg
    double param_dsds[18] {0};
    double param_dsds1[18] {0};
    double param_dsds2[18] {0};
    double param_dds[18] {0};
    double abs_param_dds[18] {0};//for vLmt of stanceLeg
    double param_ds[18] {0};//for aLmt of swingLeg
    double param_ds1[18] {0};
    double param_ds2[18] {0};
    double param_const[18] {0};//for aLmt of swingLeg
    double param_const1[18] {0};
    double param_const2[18] {0};
    double param_a2[18] {0};//coefficient of ds*ds, param_dsds divided by param_dds
    double param_a1[18] {0};//coefficient of ds, param_ds divided by param_dds
    double param_a0L[18] {0};
    double param_a0H[18] {0};

    double output_PeeB[1801][18] {{0}};
    double output_dsds[1801][18] {{0}};
    double output_dsds1[1801][18] {{0}};
    double output_dsds2[1801][18] {{0}};
    double output_ds[1801][18] {{0}};
    double output_ds1[1801][18] {{0}};
    double output_ds2[1801][18] {{0}};
    double output_const[1801][18] {{0}};
    double output_const1[1801][18] {{0}};
    double output_const2[1801][18] {{0}};
    double output_dds[1801][18] {{0}};
    double output_a2[1801][18] {{0}};
    double output_a1[1801][18] {{0}};
    double output_a0L[1801][18] {{0}};
    double output_a0H[1801][18] {{0}};

    double ds_upBound_aLmt[1801][6] {{0}};
    double ds_lowBound_aLmt[1801][6] {{0}};
    double ds_upBound_vLmt[1801][6] {{0}};
    double ds_upBound[1801][6] {{0}};
    double ds_lowBound[1801][6] {{0}};
    double dds_lowBound[1801][6] {{0}};
    double dds_upBound[1801][6] {{0}};

    double slopedsBound[1801][6] {{0}};
    double slopeDelta[1801][6] {{0}};
    double vaDelta[1801][6] {{0}};
    double paramdds0Point[1801][18] {{0}};
    double vaCrossPoint[1801][6] {{0}};
    double tangentPoint[1801][6] {{0}};
    double switchPoint[1801][6] {{0}};
    int paramdds0Count[18] {{0}};
    int vaCrossCount[6] {{0}};
    int tangentCount[6] {{0}};
    int switchCount[6] {{0}};
    bool quitSwitchPoint {false};
    double ds_forward[1801][6] {{0}};
    double ds_backward[1801][6] {{0}};
    double dds_forward[1801][6] {{0}};
    double dds_backward[1801][6] {{0}};
    bool stopFlag {false};
    bool accFlag {true};//true acc, false dec
    int dec_start {0};
    int dec_end {0};
    double min_dist[1801];
    unsigned int cycleCount {0};//used to limit the cycleCount of forward integration
    int stop_back {0};
    int ki_back {0};
    int ki_for {0};

    double real_ds[1801][6] {{0}};
    double real_dds[1801][6] {{0}};
    double real_ddsMax[1801][6] {{0}};
    double real_ddsMin[1801][6] {{0}};

    double timeArray[1801][6] {{0}};
    double timeArray_tmp[1801][6] {{0}};
    double totalTime[6] {0};
    double maxTime {0};
    int maxTotalCount {0};
    int maxTotalCount_last {1};
    int maxTimeID {0};

    double output_Pee[1800][9] {{0}};
    double output_Pin[1800][9] {{0}};
    double output_Vin[1800][9] {{0}};
    double output_Ain[1800][9] {{0}};
};

#endif
