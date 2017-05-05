#include "ukf.h"
#include "tools.h"
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
  
  // This initialization missed in the initial submit
  // Had no issues with Ubuntu, so overseen this completely
  is_initialized_ = false;
  
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

	//state covariance matrix P
  // px,py are available initially and hence less uncertainty when 
  // compared v, yaw, yaw_rate
	P_ << 1., 0., 0., 0., 0.,
			  0., 1., 0., 0., 0.,
			  0., 0., 1., 0., 0.,
			  0., 0., 0., 1., 0.,
        0., 0., 0., 0., 1.;

  //set state dimension
  n_x_ = 5;

  //set state dimension
  n_aug_ = 7;

  //set spreading parameter
  lambda_ = 3. - n_x_;

  //set vector for weights
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; (i<(2*n_aug_+1)); i++) {  
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }  

  time_us_ = 0.;

  // As suggested by reviewer, this line is redundant
  //Tools tools;

  // the current NIS for radar
  NIS_radar_ = 0.0;

  // the current NIS for laser
  NIS_laser_ = 0.0;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  float rho, phi, rho_dot;
  float px, py; //, vx, vy, v; //, yaw, yaw_dot;


/*****************************************************************************
   *  Initialize
****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      rho = meas_package.raw_measurements_[0];
      phi = meas_package.raw_measurements_[1];
      rho_dot = meas_package.raw_measurements_[2];

      px = rho * cos(phi);
      py = rho * sin(phi);

      x_ << px,py, 0., 0., 0.;

      time_us_ = meas_package.timestamp_;
      
      // Check if initial measurements are zeros
      if ((rho==0) && (phi==0) && (rho_dot==0)) {
        is_initialized_ = false;
      } 
      // done initializing, no need to predict or update
      else {
        is_initialized_ = true;
      }

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      px = meas_package.raw_measurements_[0];
      py = meas_package.raw_measurements_[1];

      x_ << px, py, 0., 0., 0.;
 
      time_us_ = meas_package.timestamp_;
      
      // Check if initial measurements are zeros
      if ((px==0) && (py==0)) {
        is_initialized_ = false;
      } 
      // done initializing, no need to predict or update
      else {
        is_initialized_ = true;
      }
    }

    // done initializing, no need to predict or update
    return;  
  }

/*****************************************************************************
   *  Predict
****************************************************************************/
  float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  /*
   * As suggested by reviewer, this check is redundant
   * Required only for "sample-laser-radar-measurement-data-2.txt"
  while (delta_t > 0.09)
  {
    const double dt = 0.05;
    Prediction(dt);
    delta_t -= dt;
  }
  */
  Prediction(delta_t); 


/*****************************************************************************
   *  Update
****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates

    UpdateRadar(meas_package);

  } 
  
  else {
    // Laser updates

    UpdateLidar(meas_package);

  }  

}


/**
* AugmentedSigmaPoints: Generates augmented sigma points for a
* given state and co-variance matrix
* @param Xsig_aug matrix to store generated sigma points
*/
// function signature changed as per review comments
MatrixXd UKF::AugmentedSigmaPoints(MatrixXd &Xsig_aug) {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create augmented mean state
  x_aug << x_, 0., 0.;
  //x_aug.head(5) = x;
  //x_aug(5) = 0;
  //x_aug(6) = 0;

  //create augmented covariance matrix
  MatrixXd Q = MatrixXd(2,2);
  Q << std_a_*std_a_, 0, 0, std_yawdd_*std_yawdd_;
  
  MatrixXd P_ak = MatrixXd::Zero(7,7);
  P_ak << P_,
         MatrixXd::Zero(5,2),
         MatrixXd::Zero(2,5),
         Q;
  //P_aug.fill(0.0);
  //P_aug.topLeftCorner(5,5) = P;
  //P_aug(5,5) = std_a*std_a;
  //P_aug(6,6) = std_yawdd*std_yawdd;

  //create square root matrix
  MatrixXd A = P_ak.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  float sqrt_lambda = sqrt(lambda_+n_aug_);

  for (int i=0;i<n_aug_;i++) {
    Xsig_aug.col(i+1) = x_aug + sqrt_lambda*A.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt_lambda*A.col(i);
  }
  //print result
  //std::cout << "Xsig_aug = " << std::endl << Xsig_aug.col(0) << std::endl;

  return Xsig_aug;
}


/**
 * SigmaPointPrediction: Predict sigma points using process model
 * given augmented sigma points
 * @param Xsig_aug matrix to store generated sigma points
 */
// void SigmaPointPrediction(MatrixXd* Xsig_aug, MatrixXd* Xsig_pred_, double dt);
// Review comments - removing "MatrixXd* Xsig_pred_"
// Redundant - anyway available inside the function
void UKF::SigmaPointPrediction(const MatrixXd &Xsig_aug, double dt) {

  //predict sigma points
  MatrixXd X = MatrixXd(n_x_,2 * n_aug_ + 1);
  VectorXd Nu_a = VectorXd(2 * n_aug_ + 1);
  VectorXd Nu_yaw = VectorXd(2 * n_aug_ + 1);

  //X = Xsig_aug.topLeftCorner(5,15);
  X = Xsig_aug.topRows(5);
  Nu_a = Xsig_aug.row(5);
  Nu_yaw = Xsig_aug.row(6);

  VectorXd dX = VectorXd(n_x_);
  VectorXd noise = VectorXd(n_x_);
  //avoid division by zero
  //write predicted sigma points into right column
  for (int i=0; i<(2 * n_aug_ + 1); i++) {
    
    float v = X(2,i);
    float yaw = X(3,i);
    float yaw_dot = X(4,i);
    //std::cout << "Calculating  noise" << std::endl;
    noise << 0.5*dt*dt*cos(yaw)*Nu_a(i),
             0.5*dt*dt*sin(yaw)*Nu_a(i),
             dt*Nu_a(i),
             0.5*dt*dt*Nu_yaw(i),
             dt*Nu_yaw(i);

    if ((fabs(yaw_dot) < 0.001)) {
      //std::cout << "yaw_dot==0 dX" << std::endl;
      dX << v*cos(yaw)*dt, v*sin(yaw)*dt, 0., yaw_dot*dt, 0.;
    } else {
      //std::cout << "yaw_dot!=0 dX" << std::endl;
      dX << (v/yaw_dot)*(sin(yaw+yaw_dot*dt)-sin(yaw)),
            (v/yaw_dot)*(-cos(yaw+yaw_dot*dt)+cos(yaw)),
            0., yaw_dot*dt, 0.;
    }
    //std::cout << "Xsig_pred" << std::endl;
    Xsig_pred_.col(i) << X.col(i) + dX + noise;
    //std::cout << Xsig_pred << std::endl << std::endl;
  }
  //print result
  //std::cout << "Xsig_pred = " << std::endl << Xsig_pred.col(4) << std::endl;
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  /*****************************************************************************
  *  Generate and Augment Sigma Points
  ****************************************************************************/
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug = AugmentedSigmaPoints(Xsig_aug);


  /*****************************************************************************
  *  Predict Sigma Points
  ****************************************************************************/
  SigmaPointPrediction(Xsig_aug, dt);


  /*****************************************************************************
  *  Convert Predicted Sigma Points to Mean/Covariance
  ****************************************************************************/

  //predict state mean
  x_ = Xsig_pred_.col(0) * weights_(0);

  for (int i=1; i<(2 * n_aug_ + 1); i++) {
    x_ += Xsig_pred_.col(i) * weights_(i);
  }

  //predict state covariance matrix
  P_ = weights_(0) * (Xsig_pred_.col(0) - x_) * ((Xsig_pred_.col(0) - x_).transpose());

  for (int i=1; i<(2 * n_aug_ + 1); i++) {
    
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    x_diff = tools.AngleNormalize(x_diff, 3);

    P_ += weights_(i) * x_diff * x_diff.transpose() ;
  }

  //print result
  //std::cout << "Predicted state" << std::endl;
  //std::cout << x_ << std::endl;
  //std::cout << "Predicted covariance matrix" << std::endl;
  //std::cout << P_ << std::endl;

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

  //measurement matrix

	MatrixXd H_ = MatrixXd(2, 5);
	H_ << 1, 0, 0, 0, 0,
		    0, 1, 0, 0, 0;

  MatrixXd R_ = MatrixXd(2, 2);
  R_ << 0.0225, 0,
        0, 0.0225;

  VectorXd z = VectorXd(2);
  z << meas_package.raw_measurements_[0], // px in m
       meas_package.raw_measurements_[1]; // py in rad

  float px = x_(0);
  float py = x_(1);
  float sq = (px*px + py*py);

  if (fabs(sq) < 0.0001) {
		cout << "Update() - Error - Division by Zero" << endl;
  } else {
  	VectorXd z_pred = H_ * x_;
		VectorXd y = z - z_pred;
		MatrixXd Ht = H_.transpose();
    MatrixXd PHt = P_ * Ht;
		MatrixXd S = H_ * PHt + R_;
		MatrixXd Si = S.inverse();
		MatrixXd K = PHt * Si;

    // Normalised Innovation Squared (NIS)
    NIS_laser_ = y.transpose() * Si * y;

		//new estimate
		x_ = x_ + (K * y);
		long x_size = x_.size();
		MatrixXd I = MatrixXd::Identity(x_size, x_size);
		P_ = (I - K * H_) * P_;
	}

  //print result
  //std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  //std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

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

  /*****************************************************************************
  *  Predict Radar Measurement
  ****************************************************************************/
  //set measurement dimension, radar can measure r, phi, and r_dot
  const int n_z = 3;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  // Incoming radar measurement values
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], // rho in m
       meas_package.raw_measurements_[1], // phi in rad
       meas_package.raw_measurements_[2]; // rho_dot in m/s

  //measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z,n_z);

  double px, py, v, yaw, sq_1_2, rho, phi, rho_dot;

  //transform sigma points into measurement space
  for (int i=0;i<(2*n_aug_+1);i++) {
    px = Xsig_pred_(0,i);
    py = Xsig_pred_(1,i);
    v = Xsig_pred_(2,i);
    yaw = Xsig_pred_(3,i);

    sq_1_2 = sqrt(px*px + py*py);
    if (sq_1_2>0.01){
      rho = sq_1_2;
      phi = atan2(py,px);
      rho_dot = (px*cos(yaw)*v + py*sin(yaw)*v)/rho;
    } else {
      rho = sqrt(0.01);
      phi = 0.;
      rho_dot = (px*cos(yaw)*v + py*sin(yaw)*v)/rho;
    }
    Zsig.col(i) << rho, phi, rho_dot;
  }

  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0;i<(2*n_aug_+1);i++) {
    z_pred += Zsig.col(i) * weights_(i);
  }

  //calculate measurement covariance matrix S
  MatrixXd R = MatrixXd::Zero(n_z,n_z);
  R(0,0) = std_radr_*std_radr_;
  R(1,1) = std_radphi_*std_radphi_;
  R(2,2) = std_radrd_*std_radrd_;

  for (int i=0;i<(2*n_aug_+1);i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    z_diff = tools.AngleNormalize(z_diff, 1);
    S += weights_(i) * z_diff * z_diff.transpose() ;
  }

  // Add measurment noise co-variance matrix
  S += R;

  //print result
  //std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  //std::cout << "S: " << std::endl << S << std::endl;


  /*****************************************************************************
  *  UKF Update for Radar
  ****************************************************************************/
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);

  //calculate cross correlation matrix
  for (int i=0;i<(2*n_aug_+1);i++) {
    // residual
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    x_diff = tools.AngleNormalize(x_diff, 3);
    z_diff = tools.AngleNormalize(z_diff, 1);

    Tc += weights_(i) * x_diff * z_diff.transpose() ;
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  z_diff = tools.AngleNormalize(z_diff, 1);

  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();


  //print result
  //std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  //std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

}
