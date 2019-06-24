#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;



/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in, 
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &R_ra, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Ra_ = R_ra;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  
  VectorXd z_ = H_ * x_;
  VectorXd y_ = z - z_ ; 
  MatrixXd Ht = H_.transpose();
  MatrixXd S_ = H_ * P_ * Ht + R_;
  MatrixXd Si = S_.inverse();
  MatrixXd K_ = (P_ * Ht) * Si;
  
  x_ = x_ + (K_ * y_);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K_ * H_) * P_ ;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);
  
  Hj_ = tools_.CalculateJacobian( x_);
  
  VectorXd h(3);
  
  double rho = sqrt((px*px + py*py));
  double theta = atan2(py,px);
  double rate = (px*vx + py*vy)/rho;
  
  h << rho, theta,rate;
  
  VectorXd y = z - h;
  
  if(y(1) > M_PI){
    y(1) -= M_PI;
  }
  
  if(y(1) < M_PI){
    y(1) += M_PI;
  }
  
  MatrixXd Hj_tr = Hj_.transpose();
  MatrixXd S_ = Hj_ * P_ * Hj_tr + Ra_;
  MatrixXd Si = S_.inverse();
  MatrixXd K_ = (P_ * Hj_tr) * Si;
  
  x_ = x_ + (K_ * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K_ * Hj_) * P_ ;
   
}