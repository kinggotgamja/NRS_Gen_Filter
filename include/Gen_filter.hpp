#include <iostream>
#include <vector>
#include <numeric>


/*** Kalman Filter Class ***/
class NRS_KalmanFilter
{
    public:

        /***  General Kalman Filter Parameters ***/

        /* Transition matrix: 2x2 */
        float Phi_matrix[4];
        /* Q covariance plant noise matrix: 2x2 */
        float Q_matrix[4];
        /* Sensitivity matrix: 1X2 */
        float H_matrix[2];
        /* Observation noise: R covariance matrix 1x1 */
        float R_matrix;
        /* P plus current covariance matrix 2x2: estimate error */
        float P_plus[4];
        /* x plus current state vector 2x1: value, speed */
        float x_plus[2];

        /***  1D Kalman filter parameters ***/
        float x_pre, p_pre; 
        float Q, R;

        /***  General Kalman Filter Functions ***/

        NRS_KalmanFilter() {}
        ~NRS_KalmanFilter() {}

        float KalmanFilter(float input); 
        float KalmanFilter1D(float input);


    private:

};

/*** Moving Everage Filter Class ***/
class NRS_MovFilter
{
    public: 
        /***  Moving Everage Filter Parameters ***/

        int mv_num;
        int counter = 0;
        std::vector<float> saved_data;

        /***  Moving Everage Filter Functions ***/

        NRS_MovFilter(int mv_num_);
        NRS_MovFilter() : NRS_MovFilter(0) {}  // 기본 생성자 (다른 class 에서 array로 instance 사용시 필수)
        ~NRS_MovFilter() {}

        float MovFilter(float input);

    private:
        float output = 0;

};


/*** Frequency Filter Class ***/
class NRS_FreqFilter
{
    public:
    /*** Frequency Filter Parameters ***/
    /* Common parameters */
    double Ts; // sampling time(s)

    /* High Pass Filter Parameters */
    double HFP_timeZone[3] = {0,};  // Do not change (defualt:0)

	double HPF_cutF = 5; // cut-off frequency(Hz)
	double HPF_zeta = 0.7; // damping ratio
	
    /* Low Pass Filter Parameters */
    double PastInput = 0; // Do not change (defualt:0)
    double PastOutput = 0; // Do not change (defualt:0)

	double LPF_cutF = 5; // cut-off frequency(Hz)

    /* Band Stop Filter Parameters */
    double BSF_timeZone[3] = {0,}; // Do not change (defualt:0)

	double BSF_cutF = 5; // cut-off frequency, stop frequency(Hz)
	double BSF_BW = 5; // stop frequency width(Hz)

    /*** Frequency Filter Functions ***/
    NRS_FreqFilter(double Ts_);
    NRS_FreqFilter() : NRS_FreqFilter(0.0) {}  // 기본 생성자 (다른 class 에서 array로 instance 사용시 필수)

    ~NRS_FreqFilter() {}

    float HPF(float input);
    float LPF(float input);
    float BSF(float input);
};
