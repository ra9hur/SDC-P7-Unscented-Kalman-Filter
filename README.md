# Unscented Kalman Filter Project

In this project, Unscented Kalman filter (UKF) is implemented using the Constant Turn-Rate and velocity magnitude (CTRV) motion model.

All Kalman filters have the same three steps:

1. Initialization
2. Prediction
3. Update

A standard Kalman filter can only handle linear equations. Both the extended Kalman filter (EKF) and the unscented Kalman filter (UKF) allow you to use non-linear equations; the difference between EKF and UKF is how they handle non-linear equations. But the basics are the same: initialize, predict, update.

EKF uses the method called first order Taylor expansion to obtain linear approximation of the non-linear measurements. In  highly non-linear systems, covariance after update step, may not be gaussian. Applying EKF may cause significant errors because of the propagation of uncertainty through the nonlinear system.

The idea with UKF is to produce several sampling points (Sigma points) around the current state estimate based on its covariance. Then, propagating these points through the nonlinear map to get more accurate estimation of the mean and covariance of the mapping results. In this way, it avoids the need to calculate the Jacobian, hence incurs only the similar computation load as the EKF.

---

## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt`


