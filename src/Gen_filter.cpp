#include "Gen_filter.hpp"

/*** Kalman Filter Functions ***/

float NRS_KalmanFilter::KalmanFilter(float input)
{
    float P_minus[4]; /* matrix 2x2 */
    float x_minus[2]; /* vector 2x1 */
    float K_gain[2];  /* matrix 2x1 */
    float temp_help;
    
    /* Prediction Step */
    x_minus[0] = Phi_matrix[0]*x_plus[0] + Phi_matrix[1]*x_plus[1];
    x_minus[1] = Phi_matrix[2]*x_plus[0] + Phi_matrix[3]*x_plus[1];
    P_minus[0] = (Phi_matrix[0]*P_plus[0] + Phi_matrix[1]*P_plus[2])*Phi_matrix[0];
    P_minus[0] += (Phi_matrix[0]*P_plus[1] + Phi_matrix[1]*P_plus[3])*Phi_matrix[1];
    P_minus[0] += Q_matrix[0];
    P_minus[1] = (Phi_matrix[0]*P_plus[0] + Phi_matrix[1]*P_plus[2])*Phi_matrix[2];
    P_minus[1] += (Phi_matrix[0]*P_plus[1] + Phi_matrix[1]*P_plus[3])*Phi_matrix[3];
    P_minus[1] += Q_matrix[1];
    P_minus[2] = (Phi_matrix[2]*P_plus[0] + Phi_matrix[3]*P_plus[2])*Phi_matrix[0];
    P_minus[2] += (Phi_matrix[2]*P_plus[1] + Phi_matrix[3]*P_plus[3])*Phi_matrix[1];
    P_minus[2] += Q_matrix[2];
    P_minus[3] = (Phi_matrix[2]*P_plus[0] + Phi_matrix[3]*P_plus[2])*Phi_matrix[2];
    P_minus[3] += (Phi_matrix[2]*P_plus[1] + Phi_matrix[3]*P_plus[3])*Phi_matrix[3];
    P_minus[3] += Q_matrix[3];
    /* Kalman Gain */
    temp_help = (H_matrix[0]*P_minus[0] + H_matrix[1]*P_minus[2])*H_matrix[0];
    temp_help += (H_matrix[0]*P_minus[1] + H_matrix[1]*P_minus[3])*H_matrix[1];
    temp_help += R_matrix;
    K_gain[0] = (H_matrix[0]*P_minus[0] + H_matrix[1]*P_minus[1])/temp_help; /* temp_help shall be !=0 */
    K_gain[1] = (H_matrix[0]*P_minus[2] + H_matrix[1]*P_minus[3])/temp_help;
    /* Correction Step */
    P_plus[0] = (1.0 - K_gain[0]*H_matrix[0])*P_minus[0] - K_gain[0]*H_matrix[1]*P_minus[2];
    P_plus[1] = (1.0 - K_gain[0]*H_matrix[0])*P_minus[1] - K_gain[0]*H_matrix[1]*P_minus[3];
    P_plus[2] = -K_gain[1]*H_matrix[0]*P_minus[0] + (1.0 - K_gain[1]*H_matrix[1])*P_minus[2];
    P_plus[3] = -K_gain[1]*H_matrix[0]*P_minus[1] + (1.0 - K_gain[1]*H_matrix[1])*P_minus[3];
    x_plus[0] = x_minus[0] + K_gain[0]*(input - x_minus[0]);
    x_plus[1] = x_minus[1] + K_gain[1]*(input - x_minus[0]);
    
    return x_plus[0];
}

float NRS_KalmanFilter::KalmanFilter1D(float input)
{
    float x_mi, p_mi, K, x, p;

    x_mi = x_pre;
    p_mi = p_pre + Q;
    
    K = p_mi/(p_mi+R);
    x = x_mi + K*(input-x_mi);
    p = (1-K)*p_mi;
    
    x_pre = x;
    p_pre = p;

    return x;
}

/*** Moving Everage Filter Functions ***/

NRS_MovFilter::NRS_MovFilter(int mv_num_)
: mv_num(mv_num_) {}

float NRS_MovFilter::MovFilter(float input)
{    
    if(mv_num > counter)
    {   
        saved_data.push_back(input);
        counter++;
    }
    else
    {
        output = accumulate(saved_data.begin(), saved_data.end(), 0.0); // begin, end, initial sum value
        output /= (float)mv_num;

        /* Initialization */
        fill(saved_data.begin(),saved_data.end(),0.0);
        counter = 0;
    }

    return output;
}

NRS_FreqFilter::NRS_FreqFilter(double Ts_)
:Ts(Ts_) {}

float NRS_FreqFilter::HPF(float input)
{
    double HPF_out = 0;

    double w0, T, Q;
    double a0_, a1_, a2_, b0_, b1_, b2_;
    double a1, a2, b0, b1, b2;
    double H0 = 1;
    double sum0;
    
    w0 = 2*(3.141592)*HPF_cutF;
    T = Ts; // sampling time
    Q = 1/(2*HPF_zeta);
    
    a0_ = 4/(T*T) + 2*w0/(Q*T) + (w0*w0);
    a1_ = -8/(T*T) + 2*(w0*w0);
    a2_ = 4/(T*T) - 2*w0/(Q*T) + (w0*w0);
    
    b0_ = 4*H0/(T*T);
    b1_ = -8*H0/(T*T);
    b2_ = 4*H0/(T*T);
    
    a1 = a1_/a0_;
    a2 = a2_/a0_;
    b0 = b0_/a0_;
    b1 = b1_/a0_;
    b2 = b2_/a0_;
    
    
    sum0 = -a1*HFP_timeZone[1] - a2*HFP_timeZone[0];
    HFP_timeZone[2] = input + sum0;
    HPF_out = b0*HFP_timeZone[2] + b1*HFP_timeZone[1] + b2*HFP_timeZone[0];

    HFP_timeZone[0] = HFP_timeZone[1];
    HFP_timeZone[1] = HFP_timeZone[2];

    return HPF_out;
}

float NRS_FreqFilter::LPF(float input)
{
    double a1,b0,b1,w0;
    double LPF_out;
    double SamplingFrequency = (1/Ts);
    w0 = 2*3.14*LPF_cutF;
    a1 = (w0 - 2*SamplingFrequency)/(2*SamplingFrequency + w0);
    b0 = w0/(2*SamplingFrequency + w0);
    b1 = b0;

    LPF_out = b0*(input) + b1*(PastInput) - a1*(PastOutput);
    PastOutput = LPF_out;
    PastInput = input;

    return LPF_out;
}

float NRS_FreqFilter::BSF(float input)
{
    double BSF_out = 0;

    double w0_peak, Q;
    double a0_, a1_, a2_, b0_, b1_, b2_;
    double a1, a2, b0, b1, b2;
    double H0 = 1;
    double sum0;

    w0_peak = 2*(3.141592)*BSF_cutF;
    Q = BSF_cutF/BSF_BW;
    
    b0_ = H0*4/(Ts*Ts) + H0*(w0_peak*w0_peak);
    b1_ = -2*H0*4/(Ts*Ts) + 2*(w0_peak*w0_peak);
    b2_ = H0*4/(Ts*Ts) + H0*(w0_peak*w0_peak);
    
    a0_ = 4/(Ts*Ts)+2*w0_peak/(Q*Ts)+(w0_peak*w0_peak);
    a1_ = -8/(Ts*Ts)+2*(w0_peak*w0_peak);
    a2_ = 4/(Ts*Ts)-2*w0_peak/(Q*Ts)+(w0_peak*w0_peak);
    
    a1 = a1_/a0_;
    a2 = a2_/a0_;
    b0 = b0_/a0_;
    b1 = b1_/a0_;
    b2 = b2_/a0_;
    
    
    sum0 = -a1*BSF_timeZone[1] - a2*BSF_timeZone[0];
    BSF_timeZone[2] = input + sum0;
    BSF_out = b0*BSF_timeZone[2] + b1*BSF_timeZone[1] + b2*BSF_timeZone[0];

    BSF_timeZone[0] = BSF_timeZone[1];
    BSF_timeZone[1] = BSF_timeZone[2];

    return BSF_out;
}