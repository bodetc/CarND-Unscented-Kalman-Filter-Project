#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "measurement_update.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_ = false;

  ///* if this is false, laser measurements will be ignored (except for init)
  const bool use_laser_ = true;

  ///* if this is false, radar measurements will be ignored (except for init)
  const bool use_radar_ = true;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  const double std_a_ = 0.2;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  const double std_yawdd_ = 0.2;

  ///* Laser measurement noise standard deviation position1 in m
  const double std_laspx_ = 0.15;

  ///* Laser measurement noise standard deviation position2 in m
  const double std_laspy_ = 0.15;

  ///* Radar measurement noise standard deviation radius in m
  const double std_radr_ = 0.3;

  ///* Radar measurement noise standard deviation angle in rad
  const double std_radphi_ = 0.03;

  ///* Radar measurement noise standard deviation radius change in m/s
  const double std_radrd_ = 0.3;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  const int n_x_ = 5;

  ///* Augmented state dimension
  const int n_aug_ = n_x_ + 2;

  ///* Sigma point spreading parameter
  const double lambda_ = 3. - n_aug_;

  ///* set measurement dimension, radar can measure r, phi, and r_dot
  const int n_z_radar_ = 3;

  ///* set measurement dimension, lidar measures p_x and p_y
  const int n_z_lidar_ = 2;

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
   * @param meas_package The latest measurement data of either radar or laser
   */
  void Prediction(MeasurementPackage meas_package);

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

private:
  /**
   * Initialize
   * @param meas_package The latest measurement data of either radar or laser
   */
  void Initialize(MeasurementPackage meas_package);

  MatrixXd CreateAugmentedSigmaPoints();

  MatrixXd GenerateSigmaPoints(const VectorXd& x, const MatrixXd& P, int n);

  void PredictSigmaPoint(const MatrixXd& Xsig_aug, double delta_t);

  void PredictMeanAndCovariance();

  MatrixXd TransformSigmaPointsToRadarSpace();
  MatrixXd TransformSigmaPointsToLidarSpace();

  void UpdateState(const MatrixXd& Zsig, const VectorXd& z);
};

#endif /* UKF_H */
