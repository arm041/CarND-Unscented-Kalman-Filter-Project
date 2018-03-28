#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  //set weights

  is_initialized_ = false;
  time_us_ = 0;
  n_x_ = 5;
  n_aug_ = 7;

  Xsig_pred_ = MatrixXd (n_x_, 2 * n_aug_ + 1);

  weights_ = VectorXd(2*n_aug_+1);
  lambda_ = 3 - n_aug_;
  weights_ (0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1 ; i < 2 * n_aug_ + 1 ; i++)

    weights_ (i) = 1 / (2 * (lambda_ + n_aug_));

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */



	  /*****************************************************************************
	   *  Initialization
	   ****************************************************************************/
	  if (!is_initialized_) {
	    // first measurement
	    cout << "UKF: " << endl;

	    x_ << 1, 1, 1, 1, 1;
	    P_ << 1, 0, 0, 0, 0,
	    	  0, 1, 0, 0, 0,
	    	  0, 0, 1, 0, 0,
	    	  0, 0, 0, 1, 0,
	    	  0, 0, 0, 0, 1;

	    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
	      /**
	      Convert radar from polar to cartesian coordinates and initialize state.
	      */
	    	cout << "In Radar initializing";
	    	float rho = meas_package.raw_measurements_(0);
	    	float phi = meas_package.raw_measurements_(1);
	    	float rhodot = meas_package.raw_measurements_(2);
	    	x_ << rho * sin(phi), rho * cos (phi), 0, 0, 0;//rhodot * cos(phi), rhodot * sin (phi);
	    }
	    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
	      /**
	      Initialize state.
	      */
	    	cout << "In lidar initializing!!!" << endl;;
	    	float px = meas_package.raw_measurements_(0);
	    	float py = meas_package.raw_measurements_(1);
	    	x_ << px, py, 0, 0, 0;


	    }


	    time_us_ = meas_package.timestamp_;
	    // done initializing, no need to predict or update
	    is_initialized_ = true;
	    return;
	  }

	  /*****************************************************************************
	   *  Prediction
	   ****************************************************************************/

		double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
		time_us_ = meas_package.timestamp_;

		Prediction (dt);


	  /*****************************************************************************
	   *  Update
	   ****************************************************************************/


	  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
	    // Radar updates
		  UpdateRadar(meas_package);
	  } else {
	    // Laser updates
		  UpdateLidar(meas_package);
	  }

	  // print the output
	  cout << "x_ = " << x_ << endl;
	  cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

	//First generating the sigma points in augmented mode
	MatrixXd Xsigma_aug = MatrixXd (n_aug_, 2 * n_aug_ + 1);
	AugmentedSigmaPoints(&Xsigma_aug);

	//Use prediction for the sigma points
	SigmaPointPrediction(&Xsigma_aug, delta_t);


	  //create vector for predicted state
	  VectorXd x = VectorXd(n_x_);

	  //create covariance matrix for prediction
	  MatrixXd P = MatrixXd(n_x_, n_x_);


	  //predict state mean
	  //predict state covariance matrix
	  //predicted state mean
	   x.fill(0.0);
	   for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
	     x = x+ weights_(i) * Xsig_pred_.col(i);
	   }

	   //predicted state covariance matrix
	   P.fill(0.0);
	   for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

	      // state difference
	      VectorXd x_diff = Xsig_pred_.col(i) - x;
	      //angle normalization
	      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
	      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

	      P = P + weights_(i) * x_diff * x_diff.transpose() ;
	    }
	  x_ = x;
	  P_ = P;

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


	 VectorXd x = VectorXd(n_x_);
	  x = x_;

	  //create example matrix for predicted state covariance
	  MatrixXd P = MatrixXd(n_x_,n_x_);
	  P = P_;
	  //set measurement dimension, radar can measure r, phi, and r_dot
	  int n_z = 2;

	  VectorXd z = VectorXd(n_z);
	  z = meas_package.raw_measurements_;

	  //create matrix for sigma points in measurement space
	  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	  //transform sigma points into measurement space
	  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

	    // extract values for better readibility
	    double p_x = Xsig_pred_(0,i);
	    double p_y = Xsig_pred_(1,i);

	    // measurement model
	    Zsig(0,i) = p_x;
	    Zsig(1,i) = p_y;
	  }

	  //mean predicted measurement
	  VectorXd z_pred = VectorXd(n_z);
	  z_pred.fill(0.0);
	  for (int i=0; i < 2*n_aug_+1; i++) {
	      z_pred = z_pred + weights_(i) * Zsig.col(i);
	  }

	  //innovation covariance matrix S
	  MatrixXd S = MatrixXd(n_z,n_z);
	  S.fill(0.0);
	  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
	    //residual
	    VectorXd z_diff = Zsig.col(i) - z_pred;

	    //angle normalization
	    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
	    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

	    S = S + weights_(i) * z_diff * z_diff.transpose();
	  }

	  //add measurement noise covariance matrix
	  MatrixXd R = MatrixXd(n_z,n_z);
	  R <<    std_laspx_ * std_laspx_, 0,
	          0,std_laspy_ * std_laspy_;
	  S = S + R;

	  //create matrix for cross correlation Tc
	  MatrixXd Tc = MatrixXd(n_x_, n_z);

	   //calculate cross correlation matrix
	   Tc.fill(0.0);
	   for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

	     //residual
	     VectorXd z_diff = Zsig.col(i) - z_pred;
	     //angle normalization
	     while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
	     while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

	     // state difference
	     VectorXd x_diff = Xsig_pred_.col(i) - x_;
	     //angle normalization
	     while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
	     while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

	     Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	   }

	   //Kalman gain K;
	   MatrixXd K = Tc * S.inverse();

	   //residual
	   VectorXd z_diff = z - z_pred;

	   //angle normalization
	   while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
	   while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

	   //update state mean and covariance matrix
	   x = x + K * z_diff;
	   P = P - K*S*K.transpose();

	   x_ = x;
	   P_ = P;

	   double lidar_NIS = z_diff.transpose() * S.inverse() * z_diff;
	   cout << "Lidar NIS is: " << lidar_NIS << endl;

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


	 VectorXd x = VectorXd(n_x_);
	  x = x_;

	  //create example matrix for predicted state covariance
	  MatrixXd P = MatrixXd(n_x_,n_x_);
	  P = P_;
	  //set measurement dimension, radar can measure r, phi, and r_dot
	  int n_z = 3;

	  VectorXd z = VectorXd(n_z);
	  z = meas_package.raw_measurements_;

	  //create matrix for sigma points in measurement space
	  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	  //transform sigma points into measurement space
	  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

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

	  //mean predicted measurement
	  VectorXd z_pred = VectorXd(n_z);
	  z_pred.fill(0.0);
	  for (int i=0; i < 2*n_aug_+1; i++) {
	      z_pred = z_pred + weights_(i) * Zsig.col(i);
	  }

	  //innovation covariance matrix S
	  MatrixXd S = MatrixXd(n_z,n_z);
	  S.fill(0.0);
	  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
	    //residual
	    VectorXd z_diff = Zsig.col(i) - z_pred;

	    //angle normalization
	    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
	    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

	    S = S + weights_(i) * z_diff * z_diff.transpose();
	  }

	  //add measurement noise covariance matrix
	  MatrixXd R = MatrixXd(n_z,n_z);
	  R <<    std_radr_*std_radr_, 0, 0,
	          0, std_radphi_*std_radphi_, 0,
	          0, 0,std_radrd_*std_radrd_;
	  S = S + R;

	  //create matrix for cross correlation Tc
	  MatrixXd Tc = MatrixXd(n_x_, n_z);

	   //calculate cross correlation matrix
	   Tc.fill(0.0);
	   for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

	     //residual
	     VectorXd z_diff = Zsig.col(i) - z_pred;
	     //angle normalization
	     while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
	     while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

	     // state difference
	     VectorXd x_diff = Xsig_pred_.col(i) - x_;
	     //angle normalization
	     while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
	     while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

	     Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	   }

	   //Kalman gain K;
	   MatrixXd K = Tc * S.inverse();

	   //residual
	   VectorXd z_diff = z - z_pred;

	   //angle normalization
	   while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
	   while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

	   //update state mean and covariance matrix
	   x = x + K * z_diff;
	   P = P - K*S*K.transpose();

	   x_ = x;
	   P_ = P;

	   double radar_NIS = z_diff.transpose() * S.inverse() * z_diff;
	   cout << "Radar NIS is: " << radar_NIS << endl;


}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {

	//create augmented mean vector
	VectorXd x_aug = VectorXd(7);

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(7, 7);

	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	/*******************************************************************************
	 * Student part begin
	 ******************************************************************************/

	MatrixXd Q = MatrixXd (2, 2);
	Q << std_a_ * std_a_, 0,
			0, std_yawdd_ * std_yawdd_;

	//create augmented mean state
	x_aug.head (n_x_) = x_;
	x_aug (n_aug_ - 2) = 0;
	x_aug (n_aug_ - 1) = 0;

	//create augmented covariance matrix
	P_aug.topLeftCorner (P_.rows(),P_.cols()) = P_;
	P_aug.topRightCorner (7 - Q.rows(), 7 - P_.cols()).setZero();
	P_aug.bottomRightCorner (Q.rows(), Q.cols()) = Q;
	P_aug.bottomLeftCorner (7 - P_.rows(), 7 - Q.cols()).setZero();

	//create square root matrix
	MatrixXd A = MatrixXd (7,7);
	A = P_aug.llt().matrixL();

	//create augmented sigma points
	Xsig_aug.col (0) = x_aug;
	for (int i = 0; i < n_aug_ ; i++)
	{
		Xsig_aug.col (i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
		Xsig_aug.col (n_aug_ + i + 1) = x_aug - sqrt (lambda_ + n_aug_) * A.col(i);

	}

	//write result
	*Xsig_out = Xsig_aug;

}



void UKF::SigmaPointPrediction(MatrixXd *Xsig_aug, double delta_t) {

  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column


  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  VectorXd Noise = VectorXd (n_x_);
  double delta_t_2 = delta_t * delta_t;

  VectorXd Process = VectorXd (n_x_);

  for (int i = 0 ; i < 2 * n_aug_ + 1 ; i++)
  {
      Noise << delta_t_2 * cos (Xsig_aug->col(i)(3)) * Xsig_aug->col(i) (5) / 2,
           delta_t_2 * sin (Xsig_aug->col(i)(3)) * Xsig_aug->col(i) (5) / 2,
           delta_t * Xsig_aug->col(i) (5),
           delta_t_2 * Xsig_aug->col(i) (6) / 2,
           delta_t * Xsig_aug->col(i) (6);
        Process << Xsig_aug->col(i) (2) * cos(Xsig_aug->col(i)(3)) * delta_t,
        		Xsig_aug->col(i) (2) * sin(Xsig_aug->col(i)(3)) * delta_t,
             0,
             delta_t * Xsig_aug->col(i)(4),
             0;
      if (Xsig_aug->col(i) (4) > 0.001)
      {
          Process (0) = Xsig_aug->col(i)(2) / Xsig_aug->col(i)(4) * (sin(Xsig_aug->col(i)(3) + delta_t * Xsig_aug->col(i)(4))-sin(Xsig_aug->col(i)(3)));
            Process (1) = Xsig_aug->col(i)(2) / Xsig_aug->col(i)(4) * (cos(Xsig_aug->col(i)(3))-cos(Xsig_aug->col(i)(3) + delta_t * Xsig_aug->col(i)(4)));
      }


        Xsig_pred.col (i) = Xsig_aug->col(i).head (n_x_) + Process + Noise;




  }


  //write result
  Xsig_pred_ = Xsig_pred;

}
