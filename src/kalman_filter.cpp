#include <iostream>
#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  
  // Predict state x and P

  x_ = F_ * x_;   // x = Fx + u, where here u = 0;
  P_ = F_ * P_ * F_.transpose() + Q_;
  

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

  // Laser uses normal Update()
  VectorXd z_pred = H_ * x_;    // Note: z is from measurement pack
  VectorXd y = z - z_pred;   

  // R_ is set to R_laser_ in FusionEKF::ProcessMeasurement() 
  MatrixXd S = H_ * P_ * H_.transpose() + R_; 
  MatrixXd K =  P_ * H_.transpose() * S.inverse();

  // Update state x and new P with a new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */


  // For EKF:
  // y = z - h(x) instead of y = z - H*x
  // Also: For radar (EKF), use Hj to calc S, K, P; for lidar (reg KF), use H to calculate y, S, K, P

  // Define variables in h(x) to map state to polar coordinates
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  // Check division by zero
  if (fabs(px*px + py*py) < 0.0001) {
    std::cout << "UpdateEKF(): CalculateJacobian() - Error - Division by Zero" << std::endl;
    return;
  }


  VectorXd z_pred = VectorXd(3);    // z_pred is another name for h(x) for radar case
  float rho = sqrt(px*px + py*py); 
  float phi = atan2(py,px);         // note atan2() returns degrees in range [-pi,pi]
  float rho_dot = (px*vx + py*vy)/rho;  
  z_pred << rho, phi, rho_dot;

  VectorXd y = z - z_pred;

  // Check to make sure phi in vector y is in the range [-pi, pi]...
  // ... otherwise +/- 2pi until it is
  // Note: y(1) is phi
  
  // normalize the angle = atan2(sin(angle),cos(angle));
  y(1) = atan2(sin(y(1)), cos(y(1)));
  const double pi = 3.14159265358979323846;
  if (y(1) < -pi) {
    y(1) = y(1) + 2*pi;
  }
  else if (y(1) > pi) {
    y(1) = y(1) - 2*pi;
  }



  // Use Hj (Jacobian) to calculate S, K, P
  // H_ is set to Hj_ in FusionEKF::ProcessMeasurement() before calling this function
  // R_ is set to R_radar_ in FusionEKF::ProcessMeasurement() before calling this function
  MatrixXd S = H_ * P_ * H_.transpose() + R_; 
  MatrixXd K = P_ * H_.transpose() * S.inverse();

  // Update state x and new P with a new estimate
  x_ = x_ + (K * y);

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}