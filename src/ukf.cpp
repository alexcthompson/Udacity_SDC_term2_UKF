#include "ukf.h"
#include "Eigen/Dense"
#include "tools.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // dims and spreading lambda
  n_x_ = 5;
  n_aug_ = 7;
  n_radar_ = 3;
  n_lidar_ = 2;
  lambda_ = 3 - n_aug_;

  // state vectors
  x_ = VectorXd::Zero(n_x_);
  z_pred_radar_ = VectorXd::Zero(n_radar_);
  z_pred_lidar_ = VectorXd::Zero(n_lidar_);
  z_meas_radar_ = VectorXd::Zero(n_radar_);
  z_meas_lidar_ = VectorXd::Zero(n_lidar_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // measurement covariance matrices
  S_radar_ = MatrixXd::Zero(n_radar_, n_radar_);
  S_lidar_ = MatrixXd::Zero(n_lidar_, n_lidar_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  // EKF initialization
  is_initialized_ = false;

  // sigma point matrices
  Xsig_aug_ = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);
  Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);

  // predicted measurement matrices
  Zsig_radar_ = MatrixXd::Zero(n_radar_, 2 * n_aug_ + 1);
  Zsig_lidar_ = MatrixXd::Zero(n_lidar_, 2 * n_aug_ + 1);

  // set up the weights, since they remain constant throughout
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  weights_.tail(2 * n_aug_) = VectorXd::Constant(2 * n_aug_, 1 / (2 * (lambda_ + n_aug_)));
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

  // INITIALIZATION
  if (!is_initialized_) {
    cout << "==== Initializing EKF on " << meas_package.sensor_type_ << " ====" << endl;
    // velocity, yaw, or yaw rate cannot be fully inferred regardless of measurement source
    x_(2) = 0.0;
    x_(3) = 0.0;
    x_(4) = 0.0;

    time_us_ = meas_package.timestamp_;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      float ro = meas_package.raw_measurements_(0);
      float theta = meas_package.raw_measurements_(1);

      x_(0) = cos(theta) * ro;
      x_(1) = sin(theta) * ro;

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    }

    // set velocity, yaw, yaw_rate high to allow early measurements to update fast
    P_ << 0.01,    0,    0,    0,    0,
         0, 0.01,    0,    0,    0,
         0,    0,   30,    0,    0,
         0,    0,    0,    9,    0,
         0,    0,    0,    0,    2;

    is_initialized_ = true;
    return;
  }

  // Prediction
  float delta_t = (meas_package.timestamp_ - time_us_) / 1.0e6;
  Prediction(delta_t);

  // Update
  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
  else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  GenAugmentedSigmaPoints();
  GenMatrixPrediction(delta_t);
  UpdateMeanCovInStateSpace();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  PredictLidarMeasurement();
  PredictMeanCovInLidarSpace();
  time_us_ = meas_package.timestamp_;
  z_meas_lidar_ = meas_package.raw_measurements_;
  UpdateStateWithLidar();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  // if bicycle is too close to radar, skip radar update = too much variation in yaw
  if (abs(meas_package.raw_measurements_[0]) < 0.1) { return;}

  PredictRadarMeasurement();
  PredictMeanCovInRadarSpace();
  time_us_ = meas_package.timestamp_;
  z_meas_radar_ = meas_package.raw_measurements_;
  UpdateStateWithRadar();
}

// generates augmented sigma points
void UKF::GenAugmentedSigmaPoints() {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  MatrixXd P_aug = MatrixXd::Zero(7, 7);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  // square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0) = x_aug;
  
  for(int i = 0; i < n_aug_; i++){
    Xsig_aug_.col(i + 1)          = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
    Xsig_aug_.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }
}

// applies the process model to a point
VectorXd UKF::ProcessModel(VectorXd point, double delta_t) {
  double v     = point(2);
  double psi   = point(3);
  double psid  = point(4);
  double nua   = point(5);
  double nupsi = point(6);

  VectorXd process = VectorXd(5);
  VectorXd process_noise =VectorXd(5);

  if(psid == 0.0)
  {
    process(0) = v * cos(psi) * delta_t;
    process(1) = v * sin(psi) * delta_t;
    process(2) = 0.0;
    process(3) = psid * delta_t;
    process(4) = 0.0;

    process_noise(0) = 0.5 * delta_t * delta_t * cos(psi) * nua;
    process_noise(1) = 0.5 * delta_t * delta_t * sin(psi) * nua;
    process_noise(2) = delta_t * nua;
    process_noise(3) = 0.5 * delta_t * delta_t * nupsi;
    process_noise(4) = delta_t * nupsi;
  }
  else
  {
    process(0) = (v / psid) * (  sin(psi + psid * delta_t) - sin(psi));
    process(1) = (v / psid) * (- cos(psi + psid * delta_t) + cos(psi));
    process(2) = 0.0;
    process(3) = psid * delta_t;
    process(4) = 0.0;

    process_noise(0) = 0.5 * delta_t * delta_t * cos(psi) * nua;
    process_noise(1) = 0.5 * delta_t * delta_t * sin(psi) * nua;
    process_noise(2) = delta_t * nua;
    process_noise(3) = 0.5 * delta_t * delta_t * nupsi;
    process_noise(4) = delta_t * nupsi;
  }
  
  return point.head(5) + process + process_noise;
}

// applies the process model to columns of a matrix
void UKF::GenMatrixPrediction(double delta_t) {
  for(int i = 0; i < Xsig_aug_.cols(); i++){
    VectorXd temp = Xsig_aug_.col(i);
    Xsig_pred_.col(i) = ProcessModel(temp, delta_t);
  }
}

void UKF::UpdateMeanCovInStateSpace() {
  // predict state mean
  x_ = VectorXd::Zero(n_x_);

  // build up x_
  for(int i = 0; i < 2 * n_aug_ + 1; i++){
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //predict state covariance matrix
  P_ = MatrixXd::Zero(n_x_, n_x_);

  // for angle wrapping
  Tools tools;

  for(int i = 0; i < 2 * n_aug_ + 1; i++){
    VectorXd dev = Xsig_pred_.col(i) - x_;
    dev[3] = tools.WrapAngle(dev[3]);
    P_ += weights_(i) * dev * dev.transpose();
  }

  // if angle has wrapped for x yaw, fix it
  x_(3) = tools.WrapAngle(x_(3));
}

// RADAR UPDATE HELPERS

VectorXd UKF::MapToRadarSpace(VectorXd point) {
  VectorXd res = VectorXd::Zero(3);
  
  double px = point[0];
  double py = point[1];
  double v = point[2];
  double yaw = point[3];
  
  res[0] = sqrt(px * px + py * py);
  res[1] = atan2(py, px);
  res[2] = (px * cos(yaw) * v + py * sin(yaw) * v) / res[0];

  return res;
}

void UKF::PredictRadarMeasurement() {
  for(int i=0; i < 2 * n_aug_ + 1; i++) {
    Zsig_radar_.col(i) = MapToRadarSpace(Xsig_pred_.col(i));
  }
}

void UKF::PredictMeanCovInRadarSpace() {
  // mean predicted measurement
  z_pred_radar_ = VectorXd::Zero(n_radar_);
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred_radar_ += weights_(i) * Zsig_radar_.col(i);
  }

  //calculate covariance matrix S
  S_radar_ = MatrixXd::Zero(n_radar_, n_radar_);

  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd zdiff = Zsig_radar_.col(i) - z_pred_radar_;
    S_radar_ += weights_(i) * zdiff * zdiff.transpose();
  }


  S_radar_(0, 0) += std_radr_ * std_radr_;
  S_radar_(1, 1) += std_radphi_ * std_radphi_;
  S_radar_(2, 2) += std_radrd_ * std_radrd_;
}

void UKF::UpdateStateWithRadar() {
  //create matrix for cross correlation Tc
  MatrixXd T = MatrixXd::Zero(n_x_, n_radar_);

  //calculate cross correlation matrix
  for(int i = 0; i < 2 * n_aug_ + 1; i++){
    T += weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig_radar_.col(i) - z_pred_radar_).transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = T * S_radar_.inverse();

  //update state mean and covariance matrix
  x_ = x_ + K * (z_meas_radar_ - z_pred_radar_);
  P_ = P_ - K * S_radar_ * K.transpose();
}

// LIDAR UPDATE HELPERS

VectorXd UKF::MapToLidarSpace(VectorXd point) {
  VectorXd res = VectorXd::Zero(2);
  
  res[0] = point[0];
  res[1] = point[1];

  return res;
}

void UKF::PredictLidarMeasurement() {
  for(int i=0; i < 2 * n_aug_ + 1; i++) {
    Zsig_lidar_.col(i) = MapToLidarSpace(Xsig_pred_.col(i));
  }
}

void UKF::PredictMeanCovInLidarSpace() {
  // mean predicted measurement
  z_pred_lidar_ = VectorXd::Zero(n_lidar_);

  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred_lidar_ += weights_(i) * Zsig_lidar_.col(i);
  }

  //calculate covariance matrix S
  S_lidar_ = MatrixXd::Zero(n_lidar_, n_lidar_);

  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd zdiff = Zsig_lidar_.col(i) - z_pred_lidar_;
    S_lidar_ += weights_(i) * zdiff * zdiff.transpose();
  }

  S_lidar_(0, 0) += std_laspx_ * std_laspx_;
  S_lidar_(1, 1) += std_laspy_ * std_laspy_;
}

void UKF::UpdateStateWithLidar() {
  //create matrix for cross correlation Tc
  MatrixXd T = MatrixXd::Zero(n_x_, n_lidar_);

  //calculate cross correlation matrix
  for(int i = 0; i < 2 * n_aug_ + 1; i++){
    T += weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig_lidar_.col(i) - z_pred_lidar_).transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = T * S_lidar_.inverse();

  //update state mean and covariance matrix
  x_ = x_ + K * (z_meas_lidar_ - z_pred_lidar_);
  P_ = P_ - K * S_lidar_ * K.transpose();
}
