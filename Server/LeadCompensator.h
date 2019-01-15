#ifndef Lead_COMPENSATOR_H
#define Lead_COMPENSATOR_H

/*
static LeadCompensator<6> compensator;
compensator.Initialize();
compensator.SetLeadParam(1.1,1.6,0.4,1000);
compensator.Compensate(raw_data,compensate_data);
*/

template <size_t WIDTH>
class LeadCompensator
{
    public:
        LeadCompensator();
        ~LeadCompensator();
        void SetLeadParam(double Kc, double T, double alpha, double sampleRate);
        void Initialize();
        void Compensate(const double (&signalIn)[WIDTH], double (&signalOut)[WIDTH]);

    private:
        double Kc;
        double T;
        double alpha;
        double dt;
        double lstIn[WIDTH];
        double inte[WIDTH];
        double diff[WIDTH];
};

template <size_t WIDTH>
LeadCompensator<WIDTH>::LeadCompensator()
{
}

template <size_t WIDTH>
LeadCompensator<WIDTH>::~LeadCompensator()
{
}

template <size_t WIDTH>
void LeadCompensator<WIDTH>::SetLeadParam(double setKc, double setT, double setAlpha, double sampleRate)
{
    Kc=setKc;
    T=setT;
    alpha=setAlpha;
    dt=1.0/sampleRate;
}

template <size_t WIDTH>
void LeadCompensator<WIDTH>::Compensate(const double (&signalIn)[WIDTH], double (&signalOut)[WIDTH])
{
    for (int i = 0; i < WIDTH; ++i)
    {
        diff[i]=(signalIn[i]-lstIn[i])/dt+1/T*signalIn[i];
        inte[i]=inte[i]+(diff[i]-inte[i]/(alpha*T))*dt;

        lstIn[i]=signalIn[i];

        signalOut[i]=Kc*inte[i];
    }
}

template <size_t WIDTH>
void LeadCompensator<WIDTH>::Initialize()
{
    for (int i = 0; i < WIDTH; ++i)
    {
        lstIn[i] = 0;
        inte[i]  = 0;
        diff[i]  = 0;
    }
}
#endif
