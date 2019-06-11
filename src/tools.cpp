#include "tools.h"
#include <iostream>

#define EPS 0.3  // 10 zu schlecht, 7 war gut
#define EPS2 0.00001 //was 0.0000001

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
  /**
   * TODO: Calculate the RMSE here.
   */
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;
	// Check the validity of the following inputs:
   	if ((estimations.size() == 0.0) || (estimations.size() != ground_truth.size()))
	{
		std::cout << "calculateRMSE() - Error: - estimation vector is size zero or not equal to ground truth!" << std::endl;
		return rmse;
	}

	for (int i = 0; i < estimations.size(); ++i)
	{
		VectorXd residual = estimations[i] - ground_truth[i];
		residual = residual.array()*residual.array();
		rmse += residual;
	}
	rmse = rmse.array() / estimations.size();
	rmse = rmse.array().sqrt();
	
    std::cout << "RMSE: " <<std::endl<< rmse << std::endl;
	return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) 
{
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
	MatrixXd Hj(3, 4);
	// recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	// TODO: YOUR CODE HERE 
	//pre-compute a set of terms to avoid repeated calculation
	// Code from lectures quizes
	
   //check of px and py for not beeing to small    
   if (fabs(px) < EPS) 
   {
	  px = EPS;
   }
   if (fabs(py) < EPS)
   {   
	  py = EPS;
   }

   // Pre-compute a set of terms to avoid repeated calculation
   float c1 = px*px+py*py;
  
   /* // Check division by zero
   if(fabs(c1) < EPS2)
   {
	c1 = EPS2;
   } */
   float c2 = sqrt(c1);
   float c3 = (c1*c2);


    //compute the Jacobian matrix
    Hj<<(px/c2),                   (py/c2),                0,      0,
       -(py/c1),                   (px/c1),                0,      0, 
       py*(vx*py - vy*px)/c3,  px*(px*vy - py*vx)/c3,  px/c2,  py/c2;
    return Hj;
}
