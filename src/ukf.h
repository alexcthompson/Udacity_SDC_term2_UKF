#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* predicted measurements
  VectorXd z_pred_radar_;
  VectorXd z_pred_lidar_;

  ///* actual measurements
  VectorXd z_meas_radar_;
  VectorXd z_meas_lidar_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* measurement covariance matrices
  MatrixXd S_radar_;
  MatrixXd S_lidar_;

  ///* augmented sigma points
  MatrixXd Xsig_aug_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* predicted measurement sigma points, radar
  MatrixXd Zsig_radar_;

  ///* predicted measurement sigma points, lidar
  MatrixXd Zsig_lidar_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Radar space dimension
  int n_radar_;

  ///* Lidar space dimension
  int n_lidar_;

  ///* Sigma point spreading parameter
  double lambda_;


  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  // Creates sigma points for the prediction update
  void GenAugmentedSigmaPoints();

  // Applies process model to a given augmented vector
  VectorXd ProcessModel(VectorXd point, double delta_t);

  // Applies process model to a matrix of points
  void GenMatrixPrediction(double delta_t);

  // Finds mean and covariance of predicted points in state space
  void UpdateMeanCovInStateSpace();

  // Maps a point to radar measurement space
  VectorXd MapToRadarSpace(VectorXd point);

  // Maps a matrix of columns to radar measurement space
  void PredictRadarMeasurement();

  // Finds mean and covariance of predicted sigma points in radar space
  void PredictMeanCovInRadarSpace();

  // Finishes the radar measurement update once sigma points and measurement are in place
  void UpdateStateWithRadar();

  // Maps a point to lidear measurement space
  VectorXd MapToLidarSpace(VectorXd point);

  // Maps a matrix of columns to lidar measurement space
  void PredictLidarMeasurement();

  // Finds mean and covariance of predicted sigma points in lidar space
  void PredictMeanCovInLidarSpace();

  // Finishes the lidar measurement update once sigma points and measurement are in place
  void UpdateStateWithLidar();

};

#endif /* UKF_H */
