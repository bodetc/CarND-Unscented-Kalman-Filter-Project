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
  P_ = MatrixXd(5, 5);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 10, 0, 0,
        0, 0, 0, 10, 0,
        0, 0, 0, 0, 10;

  // Initial Sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

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
   *  Prediction & Update
   ****************************************************************************/

  switch (meas_package.sensor_type_) {
    case MeasurementPackage::LASER:
      if (use_laser_) {
        Prediction(meas_package);
        UpdateLidar(meas_package);
      }
      break;
    case MeasurementPackage::RADAR:
      if (use_radar_) {
        Prediction(meas_package);
        UpdateRadar(meas_package);
      }
      break;
  }
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

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Convert radar from polar to cartesian coordinates and initialize state.
    double rho = meas_package.raw_measurements_(0);
    double phi = meas_package.raw_measurements_(1);

    px = rho * cos(phi);
    py = rho * sin(phi);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    px = meas_package.raw_measurements_(0);
    py = meas_package.raw_measurements_(1);
  }

  time_us_ = meas_package.timestamp_;

  // Initialize state.
  x_ << px, py, 0, 0, 0.;
  cout << "Intitial x_ = " << x_ << endl;

  // done initializing, no need to predict or update
  is_initialized_ = true;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(MeasurementPackage meas_package) {
  //compute the time elapsed between the current and previous measurements
  const double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.;
  time_us_ = meas_package.timestamp_;

  // use the Kalman filter prediction
  MatrixXd Xsig_aug = CreateAugmentedSigmaPoints();
  PredictSigmaPoint(Xsig_aug, delta_t);
  PredictMeanAndCovariance();
}

/**
 * Creates augmented Sigma points for the predicion steps
 */
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

/**
 * Predicts the position of the augmented sigma point after a timestep delta_t
 */
MatrixXd UKF::GenerateSigmaPoints(const VectorXd& x, const MatrixXd& P, int n) {
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

/**
 * Generates a Sigma point matrix
 * @param x state vector
 * @param P state covariance matrix
 * @param n number of dimensions of the state vector and covariance matrix
 * @return the generated Sigma points
 */
void UKF::PredictSigmaPoint(const MatrixXd& Xsig_aug, double delta_t) {
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

/**
 * Predicts the new mean vector and covariance matrix from the augmented Sigma points
 */
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

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  MatrixXd Zsig = TransformSigmaPointsToRadarSpace();

  //extract measurement
  VectorXd z = meas_package.raw_measurements_;

  UpdateState(Zsig, z);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  MatrixXd Zsig = TransformSigmaPointsToLidarSpace();

  //extract measurement
  VectorXd z = meas_package.raw_measurements_;

  UpdateState(Zsig, z);
}

/**
 * Transforms the Sigma points to radar space
 */
MatrixXd UKF::TransformSigmaPointsToRadarSpace() {
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  return Zsig;
}

/**
 * Transforms the Sigma points to lidar space
 */
MatrixXd UKF::TransformSigmaPointsToLidarSpace() {
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_lidar_, 2 * n_aug_ + 1);

  //simply take the first two rows (p_x and p_y)
  Zsig = Xsig_pred_.topRows(2);

  return Zsig;
}

/**
 * Update the state after a measurement and calculate the NIS
 * @param Zsig the Sigma points in measurement space
 * @param z the actual measurement
 */
void UKF::UpdateState(const MatrixXd& Zsig, const VectorXd& z) {

  // Extract sigma space dimensions
  const int n_z = Zsig.rows();
  const bool isRadar = n_z == n_z_radar_;

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    if(isRadar) {
      z_diff(1)=remainder(z_diff(1), 2*M_PI);
    }

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  if(isRadar) {
    R <<    std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0,std_radrd_*std_radrd_;
  } else {
    R <<  std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;
  }
  S = S + R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    if(isRadar) {
      z_diff(1)=remainder(z_diff(1), 2*M_PI);
    }

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3)=remainder(x_diff(3), 2*M_PI);

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  if(isRadar) {
    z_diff(1)=remainder(z_diff(1), 2*M_PI);
  }

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  // calculate NIS
  double nis = z_diff.transpose() * S.inverse() * z_diff;

  std::cout << time_us_ << ", " << nis << std::endl;
}