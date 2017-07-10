//
// Created by Cédric Bodet on 10.07.17.
//

#ifndef UNSCENTEDKF_MEASUREMENT_UPDATE_H
#define UNSCENTEDKF_MEASUREMENT_UPDATE_H

#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class MeasurementUpdate {
public:
  const int n_z_;
  const MatrixXd Zsig_;
  const VectorXd z_pred_;
  const MatrixXd S_;

  /**
   * Constructor
   */
  MeasurementUpdate(int n_z, const MatrixXd &Zsig, const VectorXd& z_pred, const MatrixXd& S)
      : n_z_(n_z), Zsig_(Zsig), z_pred_(z_pred), S_(S) {}

  /**
   * Destructor
   */
  virtual ~MeasurementUpdate() {}
};


#endif //UNSCENTEDKF_MEASUREMENT_UPDATE_H
