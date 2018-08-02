#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse = VectorXd::Zero(4);
  if (estimations.size() == 0 || estimations.size() != ground_truth.size()){
    return rmse;
  }

  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
  VectorXd res = estimations[i] - ground_truth[i];
    rmse = rmse.array() +  res.array() * res.array();
  }

  rmse = rmse.array() / estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;
}

float Tools::WrapAngle(float theta) {
  float new_theta = theta;
  while(new_theta > M_PI) {
    new_theta -= 2 * M_PI;
  }

  while(new_theta < - M_PI){
    new_theta += 2 * M_PI;
  }

  return new_theta;
}