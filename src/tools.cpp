#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

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

	// check validity of inputs
	if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
		cout << "Invalid estimation or ground truth data." << endl;
		return rmse;
	}

	// Accumulate squared residuals
	for (unsigned int i=0; i < estimations.size(); i++) {

		VectorXd residual = estimations[i] - ground_truth[i];

		// coefficient-wise multiplication
		residual = residual.array() * residual.array();
		rmse += residual;
	}

	// Calculate the mean
	rmse = rmse/estimations.size();

	// Sqrt
	rmse = rmse.array().sqrt();

	// Done
	if (rmse(0) > 0.11 || rmse(1) > 0.11) {
		cout << "x/y RMSE > 0.11: " << rmse(0) << ", " << rmse(1) << endl;
	}
	if (rmse(2) > 0.52 || rmse(3) > 0.52) {
		cout << "vx/vy RMSE > 0.52: " << rmse(2) << ", " << rmse(3) << endl;
	}

	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */

	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	float c1 = px*px + py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero
	if (fabs(c1) < 0.0001) {
		cout << "Error: CalculateJacobian() - Division by Zero or Close to Zero" << endl;
		return Hj;
	}

	// Lesson version
	Hj << (px/c2), (py/c2), 0, 0,
		  -(py/c1), (px/c1), 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	return Hj;
}