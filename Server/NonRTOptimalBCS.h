#ifndef NONRT_OPTIMAL_BCS_H
#define NONRT_OPTIMAL_BCS_H
#include <stdio.h>
#include <stdexcept>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/time.h>
#include <algorithm>

#include <aris.h>
#include <Robot_Type_I.h>
#include <Robot_Gait.h>

namespace time_optimal
{
    void matrix_dot_matrix(double *matrix1, int matrix1_row, int matrix1_col, double *matrix2, int matrix2_col, double *matrix_out);
    void dlmwrite(const char *filename, const double *mtx, const int m, const int n);
    void dlmread(const char *FileName, double *pMatrix);

    class NonRTOptimalBCS
    {
    public:
        void GetTimeOptimalGait(double step_length, double step_height, double acc_limit, double vel_limit, int leg_id, double *init_tippos, double *out_tippos, double &out_period);

    private:
        void Initialize(double step_length, double step_height, double acc_limit, double vel_limit, int leg_id, double *init_tippos);
        void GetParam();
        void GetDsBound(int count);
        void GetSwitchPoint();
        void GetOptimalDsBySwitchPoint();
        void ExtraItegrationToBoundaryPoint(int s_count, double ds);
        void ExtraItegrationToMakeBoundaryDsEqual();
        void GetConstVelocityGait();
        void GetOptimalGait2t(double *out_tippos, double &out_period);
        void outputData();
        void GetNormalGait();
        //void GetEntireGait();

        double GetMaxDec(int count, double ds);
        double GetMinAcc(int count, double ds);
        void GetTwoPointAtSwitch(double *lowPoint, double *upPoint);

        Robots::RobotTypeI rbt;
        double initPeb[6];
        int legID;

        double s[901];
        double stepD;
        double stepH;
        double aLmt;
        double vLmt;
        double initTipPos[3];
        double TipPos[901][3];

        double param_dds[901][3];
        double abs_param_dds[901][3];
        double param_dsds[901][3];
        double param_a2[901][3];
        double param_a0L[901][3];
        double param_a0H[901][3];
        int isParamddsExact0[901][3];

        double ds_upBound[901];
        double ds_upBound_aLmt[901];
        double ds_upBound_vLmt[901];
        double dds_upBound[901];
        double dds_lowBound[901];

        int dis_count1;
        int dis_count2;
        int switchCount;
        double switchPoint[901];
        char switchType[901];
        double slopeDelta[901];
        int switchScrewID[901];

        double real_ds[901];
        double real_dds[901];
        double real_ddsMax[901];
        double real_ddsMin[901];
        double ds_forward[901];
        double ds_backward[901];
        double dds_forward[901];
        double dds_backward[901];

        int totalCount;
        double v0;
        double vt;
    };

    class FastWalk
    {
    public:
        static int  fastWalk(aris::dynamic::Model &model, const aris::dynamic::PlanParamBase &param_in);
        static void parseFastWalk(const std::string &cmd, const map<std::string, std::string> &params, aris::core::Msg &msg);
        static bool GetScalingPath(Robots::WalkParam &param, int leg_id);

    private:
        static double initBodyPos[6];
        static double initTipPos[18];
        static double swingPee_scale[3001][18];
    };
}


#endif
