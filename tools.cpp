#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  //RMSE Calculation  
  
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  if((estimations.size() != ground_truth.size()) || estimations.size() == 0 ){
    std::cout << " Invalid estimation " << std::endl;
    return rmse;
  }
  
  for(unsigned int i = 0; i < estimations.size() ; ++i){
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
    }
   
  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();
  
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  // Jacobian Calculation.
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);
  
  double a1 = 1/sqrt((px*px + py*py));
  double a2 = 1/(px*px + py*py);
  double a3 = 1/sqrt((px*px + py*py) * (px*px + py*py) * (px*px + py*py));
  
  if (fabs(1/a2) < 0.0001){
    std::cout << "Error - Division by zero " << std::endl;
  }
  
  MatrixXd Hj(3,4);
  
  Hj << (px*a1) , (py * a1), 0 , 0,
        -(py * a2), (px * a2),0, 0,
        py*(vx*py - vy*px) *a3, px*(vy*px - vx*py)*a3, px*a1 , py * a1;
  
  return Hj;
        
}
