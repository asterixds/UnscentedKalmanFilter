#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  // Initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // Predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  // Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);

  // The current NIS (consistency check) for laser
  NIS_laser_ = 0.0;
  
  // The current NIS (consistency check) for radar
  NIS_radar_ = 0.0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

    double rho = meas_package.raw_measurements_.coeff(0); //rho
    double phi = meas_package.raw_measurements_.coeff(1); //phi
    x_ << rho * cos(phi), rho * sin(phi), 0, 0, 0;

  }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    x_ << meas_package.raw_measurements_.coeff(0), meas_package.raw_measurements_.coeff(1), 0, 0, 0;
  }

   //Initialize state covariance matrix
   P_ << 1, 0, 0, 0, 0, 
          0, 1, 0, 0, 0, 
          0, 0, 1, 0, 0,  
          0, 0, 0, 1, 0,  
          0, 0, 0, 0, 1;  

   time_us_ = meas_package.timestamp_;
   is_initialized_ = true;
   return;
  }
  // Calculate new elasped time and Time is measured in seconds
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;    //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;
  Prediction(dt);
  
  // Update step based on Lidar or Radar measurement
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } 
  else {
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */

void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0);
  P_aug.topLeftCorner( 5,5 ) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;
  
  //create square root matrix
  MatrixXd L =P_aug.llt().matrixL();
  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  
  //predict sigma points
  float EPS = 0.001;
  double delta_t2 = delta_t*delta_t;
  for ( int i = 0; i < 2* n_aug_ +1; i++){
      double p_x = Xsig_aug(0,i);
      double p_y = Xsig_aug(1,i);
      double v = Xsig_aug(2,i);
      double yaw = Xsig_aug(3,i);
      double yaw_rate = Xsig_aug(4,i);;
      double nu_a = Xsig_aug(5,i);
      double nu_dot = Xsig_aug(6,i);
  
      double p_xp, p_yp;
      //avoid division by zero
      if (fabs(yaw_rate) > EPS){
        p_xp = p_x + v*(sin(yaw + yaw_rate*delta_t) - sin(yaw)) /yaw_rate;
        p_yp = p_y + v*(-cos(yaw + yaw_rate*delta_t) + cos(yaw)) /yaw_rate;
      }
      else {
        p_xp = p_x + v*cos(yaw)*delta_t;
        p_yp = p_y + v*sin(yaw)*delta_t;
      }
      
      p_xp += 0.5*delta_t2*cos(yaw)*nu_a;
      p_yp += 0.5*delta_t2*sin(yaw)*nu_a;
      
      double vp = v + delta_t*nu_a;
      double yawp = yaw + yaw_rate*delta_t + 0.5*delta_t2*nu_dot; 
      double yaw_ratep = yaw_rate + delta_t*nu_dot;
      //write predicted sigma points into right column
      Xsig_pred_(0,i) = p_xp;
      Xsig_pred_(1,i) = p_yp;
      Xsig_pred_(2,i) = vp;
      Xsig_pred_(3,i) = yawp;
      Xsig_pred_(4,i) = yaw_ratep;
  }

  //set weights_
  double w0 = lambda_/(lambda_ + n_aug_);
  weights_(0) = w0; 
  for (int i = 1; i < 2*n_aug_ +1; i++){
      weights_(i) = 0.5/(lambda_+ n_aug_);
  }
  //predict state mean
  x_.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++){
    x_ =   x_ + weights_(i)*Xsig_pred_.col(i);
  }

  //predict state covariance matrix
  P_.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++){
    VectorXd diff = Xsig_pred_.col(i) - x_;
    while (diff(3)> M_PI) diff(3)-=2.*M_PI;
    while (diff(3)<-M_PI) diff(3)+=2.*M_PI;
    P_ =   P_+ weights_(i)*diff*diff.transpose();
  }  
  
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    Zsig(0,i) = Xsig_pred_(0,i);
    Zsig(1,i) = Xsig_pred_(1,i);
  }
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i)*Zsig.col(i);
  }
  //calculate measurement covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  
  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_laspx_ * std_laspx_, 0,
        0, std_laspy_ * std_laspy_;
  S += R;
  
  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(this->n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) { 

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();
  //update state mean and covariance matrix
  // calculate radar NIS
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;
  
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  NIS_laser_ = (z_diff).transpose() * S.inverse() * (z_diff);

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double yaw_rate = Xsig_pred_(4,i);
    
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);
    Zsig(1,i) =  atan2(p_y,p_x);
    Zsig(2,i) = (p_x*cos(yaw) +p_y*sin(yaw))*v/Zsig(0,i);
  }
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i)*Zsig.col(i);
  }
  //calculate measurement covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
          
  S += R;
  
  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(this->n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) { 

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();
  //update state mean and covariance matrix
  // calculate radar NIS
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;
  
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  NIS_radar_ = (z_diff).transpose() * S.inverse() * (z_diff);
  
}