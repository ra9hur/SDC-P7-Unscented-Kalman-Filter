#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  // ... your code here

  if (estimations.size() == 0) {
	  cout << "Size of estimation vector is zero" << endl;
  } 
  else if (estimations.size() != ground_truth.size()) {
	  cout << "Estimation vector size is not equal to ground truth vector size" << endl;
  }	
  else {

	  VectorXd error(4), error_sq(4), mean(4);

	  size_t n = estimations.size();

	  for(int i=0; (i < n); ++i){
      // residual error
		  error = estimations[i]-ground_truth[i];
		  error = error.array().abs();
      error_sq = error.array()*error.array();

		  //accumulate squared residuals
		  mean += error_sq;
	  }
	  //calculate the mean
	  mean = (1./n)*mean;

	  //calculate the squared root
    mean = abs(mean.array());
	  rmse = mean.array().sqrt();
  }

  //return the result
  return rmse;  
}


VectorXd Tools::AngleNormalize(VectorXd diff, int i) {
  /**
  TODO:
    * Calculations to normalize given angle here.
  */
  while (diff(i)> M_PI) diff(i)-=2.*M_PI;
  while (diff(i)<-M_PI) diff(i)+=2.*M_PI;

  return diff;
}