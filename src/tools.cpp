#include <iostream>
#include "tools.h"

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

	if (estimations.size() == 0){
	    cout << "error happened size of measurement is 0!!!!" << endl;
	    return rmse;

	}

    if (estimations.size() != ground_truth.size()){
        cout << "Error the dimensions don't fit!!!"<< endl;
        return rmse;

    }

	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
        VectorXd temp(4);
		temp = ((estimations[i] - ground_truth[i]).array() )* ((estimations[i] - ground_truth[i]).array());
	    rmse += temp;
	}

	//calculate the mean
    rmse = rmse / estimations.size();
	//calculate the squared root
    rmse = rmse.array().sqrt();

    return rmse;
}
