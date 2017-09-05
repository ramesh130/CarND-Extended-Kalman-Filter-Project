#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  // Measurement matrix H. 
  // H_laser is  the H matrix from the lesson
  // Note that Hj is calculated below in ProcessMeasurement() 
  H_laser_ << 1, 0, 0, 0,   // From lesson 5 section 10
              0, 1, 0, 0;

  // State covariance P
  // initialize P_ - I used values of 1 and 1000 from lesson
  ekf_.P_ = MatrixXd(4,4);
  ekf_.P_ << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1000, 0,
              0, 0, 0, 1000;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;    // Last 2 values can supposedly affect RMSE, need to experiment with them.

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      // Load measurement rho, phi, rho_dot into px and py.
      // Use rho and phi to obtain px and py...
      // ... but can't calculate initial vx, vy, so set them = 0.

      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];

      float px = rho * cos(phi);
      float py = rho * sin(phi);


      ekf_.x_(0) = px;    // Initialized vx, vy with 0
      ekf_.x_(1) = py;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */

      // MY CODE
      // Set initial state and zero velocity directly from sensor, since it's laser
      // Load just px, py, and vx=vy=0 since it's a laser 
      
      ekf_.x_(0) = measurement_pack.raw_measurements_[0];
      ekf_.x_(1) = measurement_pack.raw_measurements_[1];
    }
    
    ekf_.F_ = MatrixXd::Identity(4,4);

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  // MY CODE below

  // compute elapsed time in seconds between previous and current measurement
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;   // dt - expressed in seconds
      
  previous_timestamp_ = measurement_pack.timestamp_;

  // Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;   // See Lesson 8 Section 5
  ekf_.F_(1, 3) = dt;

  // Update process covariance Q with time dt and noise_ax & noise_ay
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  float noise_ax = 9;   // Noise is defined as sigma squared, so don't square it again in the matrix!
  float noise_ay = 9;

  ekf_.Q_ = MatrixXd(4,4); 
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
         0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
         dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
         0, dt_3/2*noise_ay, 0, dt_2*noise_ay;


  // Update predicted state x and P matrix
  ekf_.Predict();
  

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
   
    cout << "Enter Measurement Update RADAR" << endl;
    Hj_ = tools.CalculateJacobian(ekf_.x_);  
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);   // Update for radar


  } else {
    // Laser updates

    cout << "Enter Measurement Update LASER" << endl;
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);   // Update just for laser
  

  } 
}