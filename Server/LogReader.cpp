#include <fstream>
#include <iostream>
#include "aris.h"

int main(int argc, char* argv[])
{
    using namespace std;
    aris::control::data_emitter::Data data;
    
    if (argc < 2)
    {
        cout << "Usage: LogReader <path-to-your-log-file>" << endl;
        exit(0);
    }

    ifstream fin(argv[1], ios::in | ios::binary);
    ofstream fout("ParsedFile.txt");

    if (fin.fail())
    {
        cout << "Open file error: " << strerror(errno) << endl;
        return errno;
    }
    while(fin.read((char *)&data, sizeof(data)))
    {
        fout << data.timecount << "  ";

        // ForceData
        for(int j = 0; j < data.FOR_NUM; j++)
        {
            fout << data.force_data[j].Fx << "  ";
            fout << data.force_data[j].Fy << "  ";
            fout << data.force_data[j].Fz << "  ";
            fout << data.force_data[j].Mx << "  ";
            fout << data.force_data[j].My << "  ";
            fout << data.force_data[j].Mz << "  ";
        }

        fout << data.imu_data.roll << "  ";
        fout << data.imu_data.pitch<< "  ";
        fout << data.imu_data.yaw  << "  ";


        for(int i = 0; i < data.MOT_NUM; i++)
        {
            fout << data.motor_data[i].target_pos << "  ";
            fout << data.motor_data[i].feedback_pos<< "  ";
            fout << data.motor_data[i].target_vel << "  ";
            fout << data.motor_data[i].feedback_vel<< "  ";
            fout << data.motor_data[i].target_cur << "  ";
            fout << data.motor_data[i].feedback_cur<< "  ";
        }
        fout << endl;
    }

    fin.close();
    fout.flush();
    fout.close();

    cout << "Log file reading finished" << endl;
    return 0;
}
