#include "FusionEKF.h"
#include <iostream>
#include <math.h>
#include "Eigen/Dense"
#include "tools.h"

#define EPS 0.0001 


using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
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
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */


}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}


float check_angle_range(float angle)
{
	  float _correctedAngle =angle;
	  if (angle > M_PI)
	    _correctedAngle-=2.0*M_PI;
	  if (angle < -M_PI)
	     _correctedAngle+=2.0*M_PI;
	  return _correctedAngle;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) 
  {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
	{
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
	  float rho = measurement_pack.raw_measurements_[0]; // range
	  float phi = measurement_pack.raw_measurements_[1]; // bearing
	  float rho_dot = measurement_pack.raw_measurements_[2]; // velocity of rho
	  
	  phi = check_angle_range(phi);
	  std::cout<<"phi verified[-pi to +pi]: " << phi << std::endl;

	  // Coordinates convertion from polar to cartesian
	  float x = rho * cos(phi); 
	  float y = rho * sin(phi);

	  float vx = rho_dot * cos(phi);
	  float vy = rho_dot * sin(phi);

	  ekf_.x_ << x, y, vx , vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) 
	{
      // TODO: Initialize state.
      // We don't know velocities from the first measurement of the LIDAR, so, we use zeros instead
	  float x = measurement_pack.raw_measurements_[0];
	  float y = measurement_pack.raw_measurements_[1];
      ekf_.x_ << x, y, 0, 0;
    }
	
	if (fabs(ekf_.x_(0)) < EPS and fabs(ekf_.x_(1)) < EPS)
	{
		ekf_.x_(0) = EPS;
		ekf_.x_(1) = EPS;
	}

    // Initial covariance matrix
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
			   0, 1, 0, 0,
			   0, 0, 1000, 0,
			   0, 0, 0, 1000;

    // Print the initialization results
    cout << "EKF init: " << ekf_.x_ << endl;

    // Save the initiall timestamp for dt calculation
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  
  

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
   // Calculate the timestep between measurements in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_);
  dt /= 1000000.0; // convert micros to s
  previous_timestamp_ = measurement_pack.timestamp_;
  std::cout<<std::endl;
  std::cout<<"DeltaTime: "<<std::endl<<dt<<std::endl;
  

  // State transition matrix update
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, dt, 0,
			0, 1, 0, dt,
			0, 0, 1, 0,
			0, 0, 0, 1;

  // Noise covariance matrix computation	

  // Noise values from the task
  float noise_ax = 9.0;
  float noise_ay = 9.0;

  // Precompute some usefull values to speed up calculations of Q
  float dt_2 = dt * dt; //dt^2
  float dt_3 = dt_2 * dt; //dt^3
  float dt_4 = dt_3 * dt; //dt^4
  float dt_4_4 = dt_4 / 4; //dt^4/4
  float dt_3_2 = dt_3 / 2; //dt^3/2

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4_4 * noise_ax, 0, dt_3_2 * noise_ax, 0,
	         0, dt_4_4 * noise_ay, 0, dt_3_2 * noise_ay,
	         dt_3_2 * noise_ax, 0, dt_2 * noise_ax, 0,
 	         0, dt_3_2 * noise_ay, 0, dt_2 * noise_ay;
  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
  {
    // TODO: Radar updates
	// Use Jacobian instead of H
	ekf_.H_ = tools.CalculateJacobian(ekf_.x_);  // the H(x_) instead of the normal H because the non-linear behaviour of the radar sensor.
	ekf_.R_ = R_radar_;
	ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    std::cout<<"Radar"<<std::endl;
  } 
  else 
  {
    // Laser updates
	ekf_.H_ = H_laser_;
	ekf_.R_ = R_laser_;
	ekf_.Update(measurement_pack.raw_measurements_);
	std::cout<<"Laser"<<std::endl;
  }

  // print the output
  cout << "H_ ="<<std::endl << ekf_.H_<<std::endl;
  //cout << "R_ ="<<std::endl << ekf_.R_<<std::endl;
  cout << "x_ = " <<std::endl << ekf_.x_ << endl;
  cout << "P_ = " <<std::endl << ekf_.P_ << endl;
}
