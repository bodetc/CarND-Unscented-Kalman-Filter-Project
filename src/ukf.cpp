#include "ukf.h"
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

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Initial Sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

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

  // set weights for the calculation of the mean and covariance in the prediction step
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/(lambda_+n_aug_);
  const double weight = 0.5/(n_aug_+lambda_);
  for (int i=1; i<2*n_aug_+1; i++) {
    weights_(i) = weight;
  }
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    Initialize(meas_package);
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  //compute the time elapsed between the current and previous measurements
  const double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.;
  time_us_ = meas_package.timestamp_;

  //Call the Kalman Filter predict() function
  Prediction(delta_t);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  switch (meas_package.sensor_type_) {
    case MeasurementPackage::LASER:
      if (use_laser_) {
        UpdateLidar(meas_package);
      }
      break;
    case MeasurementPackage::RADAR:
      if (use_radar_) {
        UpdateRadar(meas_package);
      }
      break;
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}

/**
 * Initializes the
 * @param {MeasurementPackage} meas_package
 */
void UKF::Initialize(MeasurementPackage meas_package) {
  // first measurement
  cout << "EKF: " << endl;

  double px = 0.;
  double py = 0.;
  double v = 0.;
  double v_phi = 0.;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Convert radar from polar to cartesian coordinates and initialize state.
    double rho = meas_package.raw_measurements_(0);
    double phi = meas_package.raw_measurements_(1);

    px = rho * cos(phi);
    py = rho * sin(phi);

    v = meas_package.raw_measurements_(2);
    v_phi = phi; // Assumes the car drives exactly away from us

  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    px = meas_package.raw_measurements_(0);
    py = meas_package.raw_measurements_(1);
  }

  time_us_ = meas_package.timestamp_;

  // Initialize state.
  x_ << px, py, v, v_phi, 0.;
  cout << "Intitial x_ = " << x_ << endl;

  // done initializing, no need to predict or update
  is_initialized_ = true;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  MatrixXd Xsig_aug = CreateAugmentedSigmaPoints();
  PredictSigmaPoint(Xsig_aug, delta_t);
  PredictMeanAndCovariance();
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
}

MatrixXd UKF::CreateAugmentedSigmaPoints() {
  //create augmented mean state
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;

  //create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(n_x_,n_x_) = std_a_*std_a_;
  P_aug(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;

  return GenerateSigmaPoints(x_aug, P_aug, n_aug_);
}

MatrixXd UKF::GenerateSigmaPoints(VectorXd x, MatrixXd P, int n) {
  MatrixXd Xsig = MatrixXd(n, 2 * n + 1);
  Xsig.col(0)=x;

  //calculate square root of P
  MatrixXd A = P.llt().matrixL();

  for(int i = 0; i<n; i++) {
    MatrixXd Delta = sqrt(lambda_+n)*A.col(i);

    Xsig.col(1+i)=x+Delta;
    Xsig.col(1+n+i)=x-Delta;
  }

  return Xsig;
}

void UKF::PredictSigmaPoint(MatrixXd Xsig_aug, double delta_t) {
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin(yaw + yawd*delta_t) - sin(yaw)) ;
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}

void UKF::PredictMeanAndCovariance() {

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3)=remainder(x_diff(3), 2*M_PI);

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}
