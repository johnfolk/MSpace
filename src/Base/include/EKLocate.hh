// = AUTHOR(S)
//    John Folkesson
//    
//    April, 2005
 //
//    Copyright (c) 2005 John Folkesson
//    

#ifndef CURE_EKLOCATE_HH
#define CURE_EKLOCATE_HH




#include <iostream>
#include <stdlib.h> 
#include "PosedFeature.hh"
#include "Measurement.hh"
#include "Match.hh"
#include "MapFeature.hh"
#include "MapBank.hh"
#include "Matrix.hh" 
#include "Pose3D.hh"


namespace Cure{

  /**
   * A General EKF Localizer
   *
   * Initialization info:
   * 1. call initRobotPose 
   * 2. call initSensorOffset
   * 
   * @author John Folkesson 
   */
  class EKLocate
  {
  public:
    /** Pointer to the bank that manages all map features */
    MapBank *Bank;

    /** @todo comments this variable */
    int Check;


    /** This will multiply all odometry errors*/
    double Gain;

    /** The covariance matrix */
    Cure::Matrix C;

    /** 
     * @todo comment this variable. Maybe it should be called Dxp
     * since we only use P-space coordinates?
     *
     * The state vector. When using the P-space this vector contains
     * the P-space coordinates.
     */
    Cure::Matrix Dx;

    /**
     * The first row for the parameter dimensions
     *
     * @todo given that RobotRows gives the number of rows with
     * features, maybe ParameterStartRow is a better name?
     */
    int ParameterRows;

    /**
     * An array of MapFeature Keys.
     */
    //     LongList *Keys;

    /**
     * Number of sensor offset poses being used max 5.
     */
    int Sensors;

    /**
     * The rows of C cooresponding to the various offsets poses.  
     */
    int SensorRows[5];

    /**
     * The current robot pose values;
     */
    Cure::Pose3D Robot;

    /**
     * The current sensor offset pose values;
     */
    Cure::Pose3D Offsets[5];

  protected:
    /** Type of covariance for robot state variables */
    unsigned short RobotCovType;

    /** @todo Comment difference between SensorCovTypes and OffsetCovTypes */
    unsigned short SensorCovTypes[5];
    unsigned short OffsetCovTypes[5];

    /** Number of rows (state variables) for robot pose */
    int RobotRows;

    /** @todo is it even used anymore? */
    //double *AllocatedC;

  public:
    EKLocate(MapBank *b,int dim);
    ~EKLocate();
    void initRobotPose(Cure::Pose3D &pose);
    /**
     * get the coorest covarainace type for the pose of resetRobotPose.
     * @return The RobotCovType
     */
    unsigned short getRobotCovType(){return RobotCovType;}
    /**
     * Reset the pose and its uncertainty.
     * @see EKLocate::getRobotCovType.
     *
     * @param pose a the pose to reset to with the right Covariance
     * size and type.
     * @return 0 if ok else 1 
     */
    int resetRobotPose(const Cure::Pose3D &pose);

    /**
     * Here you can define offsets, for up to 5 sensors, from the robot frame
     * and have the offset updated by the EKF.
     * So the sensor pose = robotPose + offsetPose.
     * The robot pose must be i initialized first. 
     * The returned index can then be used when updating with measurments
     * from this sensor to pass which sensor info to the EKF.
     * @param offsetPose The initial value w/uncertianty of the offset.
     * @param sensroCovType this is the covariance type of the sensor pose,
     *                      sensorpose=robotpose+offset
     * @return the sensorIndex for this sensor or -1 on fail.
     */
    int initSensorOffset(Cure::Pose3D & offsetPose, int sensorCovType);
    void reallocate(int dim);


    MapFeature * getFeature(int i);
    /**
     * does pose=pose+inc and updates C.
     * @param pose the Pose state before the update in and the updated pose out
     * @param inc the incremental change in Pose with uncertainty.
     */ 
    void incrementalPredict(Cure::Pose3D & pose, Cure::Pose3D & inc); 
    /**
     * does offset[sensorIndex]=offsetpose and updates C.
     * This assumes that the new offset is independent of the old one.
     * That is ok sometimes like a typical pan/tilt uncertainty with 
     * known static offset.
     * @param offsetpose The new sensor offset with its uncertainty 
     * @param sensorIndex the index into Offsets, returned by initSensorOffset.
     */ 
    void offsetPredict(Cure::Pose3D & offsetpose, int sensorIndex); 
    /**
     * This updates the map with the measurments and adds dense 
     * info in a given frame
     * @param matches and array of length n of measurments<->PosedFeatures
     * @param n number of matches.
     * @param psensor the pose state of the measurement sensor.
     * @param map2info the transformation to take the dense info
     * from the map to the info frame.
     * @param mahalanobisTest If the M-distance is greater than this the match
     *        will be ignored.  If this is 0 the test is skipped.
     *        If the measurement fails this test no dense info is used.
     * 
     */
    void getSensorPose(Cure::Pose3D &psensor, int sensorIndex)
    {
      psensor.setCovType(SensorCovTypes[sensorIndex]);
      psensor.add_(Robot,Offsets[sensorIndex]);
    }


    /**
     * This updates the map with the measurments and adds dense 
     * int the map frame and extends the features if possible.
     *
     * @param matches and array of length n of measurments<->PosedFeatures
     * @param n number of matches.
     * @param psensor the pose state of the measurement sensor.
     *        from the map to the info frame.
     * @param map2info transformation from map coordinates to the frame to
     *        store dense information in
     * @param mahalanobisTest If the M-distance is greater than this the match
     *        will be ignored.  If this is 0 the test is skipped.
     *        If the measurement fails this test no dense info is used.
      * 
     */
    void update(Match *matches, int n, Cure::Pose3D &psensor,
		Cure::Transformation3D & map2info, 
		double mahalanobisTest);

    /**
     * This updates the map with the measurments and adds dense 
     * int the map frame.
     * @param matches and array of length n of measurments<->PosedFeatures
     * @param n number of matches.
     * @param psensor the pose state of the measurement sensor.
     * from the map to the info frame.
     * @param mahalanobisTest If the M-distance is greater than this the match
     *        will be ignored.  If this is 0 the test is skipped.
     *        If the measurement fails this test no dense info is used.
      * 
     */
    void update(Match *matches, int n, Cure::Pose3D &psensor,
		double mahalanobisTest);


    /**
     * This updates the map with the single measurment.
     * @param matches contains measurments<->PosedFeatures
     * @param robotpose the pose state of the robot.
     * @param mahalanobisTest If the M-distance is greater than this the match
     *        will be ignored.  If this is 0 the test is skipped.
     * @param sensorIndex the index into Offsets, returned by initSensorOffset.
     * @return 0 if succeed 1, if mahalanobis fails and 
     *         -1 if feature is not found.
     */
    //    int updateSensor(Match &mat, double mahalanobisTest,
    //	      int sensorIndex);

    void  updateSensor(Match *matches, int n, Pose3D &probot,
		       Transformation3D & map2info, 
		       double manhanalobisTest,int sensorIndex);
    
    //int extend(PosedFeature &pf );
    // int extend(PosedFeature &pf,int sensorIndex);
    //int extend(Pose3D &pose, PosedFeature &pf, int sensorIndex);    
    /*
      int extend(Cure::Pose3D &pose, 
	       PosedFeature &pf,int poserow=0);


      */


    //    int extend(Match &mat, int sensorIndex);
    /**
     *	do S= A*(jac1,jac2)^T + cov
     *	(xs,xp)->(xs,xp)+W*de
     *	C->C-W*A =(I-WJ)C=(I-A^TS^-1J)C
     *	A=(0,...0,jac1,0,...0,jac2,0...0)*C
     *	W=A^T*S^-1
     *
     *
     * @param A =jacobians*C
     * @param Sinv =Inverse(A*jacobians^T+cov)
     * @param cov The covariance of the measurments
     * @param jacobians[i][m] i runs 1 to numblocks m and runs from 1 to numJac
     *        for each i all the jacobians have the same number of rows while the
     *        columns are indexed int C by index. 
     * @param index  The index[i][m][j] gives the C row of the jth 
     *        column of jacobian[i][m].
     * @param numblocks This is the number of measurement blocks,
     * @param numJac  This is the number of jacobians in each block.
     */
    void getMetric(Cure::Matrix & A, Cure::Matrix & Sinv, 
			Cure::Matrix & cov, Cure::Matrix *** jacobians,
			int *** index,  int numblocks,int numJac);

    /**
     * This function goes through all the features and propagates
     * changes to the P-space coordinates to changes in the feature
     * coordinates. It also removes space in the Dx and C
     * corresponding to
     */
    //    void updateP();

    int checkEigen();
 
    /**
     * Propagates chnages in P-space to changes in the feature
     * coordinates of a certain feature.
     */
    int updateP(PosedFeature & pf);
    /**
     * This calculates the sum of log(P(v|map,probot)for the measurments.
     * Actually only upto a normalization term.
     * We have the probability as a Gaussian prob density function.
     * In order to have a probability (dimensionless) 
     * we must integrate over a small
     * interval which we take as a constant very small percentage of the 
     * measurement standard deviation.  To this should be added the 
     * normalization term.  All this doesn't ussualy make any difference
     * since it is differences of these numbers that matter.
     *
     * The arguments are similar to the updateSensor Method but of course
     * do not change the Filter estimate of the Robot pose. 
     *
     * The nomatchlikelihood is pretty important but hard to say what
     * it should be.  If comparing two robot poses where the number of 
     * matched measurements is the same then it does not matter what
     * this is.  If the number of matches is different then this better be
     * <0 or one gets no matches at all as best situation (try something 
     * like -32 my guess).
     *
     * @param matches
     * @param matches and array of length n of measurments<->PosedFeatures
     * @param n number of matches.
     * @param probot the pose state of the robot, 
     *               the sensor offset will be added.
     * @param sensorIndex the index into Offsets, returned by initSensorOffset.
     * @param nomatchlikelihood will be added into result for unmatched 
     *                          measurments.
     * @param matchtest the minimum value for an indidual measurements loglik.
     *                  below which the nomatchlikelihood is substituted.
     * @return the sum of the Gaussian exponents fro the n matches
     *              (ie -innov*(cov)^-1*innov/2).   
     */
    double getloglikelihood(Match *matches, int n,
			    Cure::Transformation3D &probot,
			    int sensorIndex,
			    double nomatchlikelihood=0,
			    double matchtest=-1E100);
    /**
     * This helps to do the above method.
     * @return the loglikelihood of the match given the sensor pose.
     */
    double getloglikelihood(Cure::Measurement *m, PosedFeature *pf,
			    Cure::Transformation3D &sensorpose);
    double getloglikelihood(Match &mat,
			    Cure::Transformation3D &sensorpose);

  protected:
    /**
     * This updates the map with the single measurment.
     * @param mat contains measurments<->PosedFeatures
     * @param pose the pose state of the measurement sensor.
     * @param mahalanobisTest If the M-distance is greater than this the match
     *        will be ignored.  If this is 0 the test is skipped.
        * @return 0 if succeed 1 if aborted.
     */
    int update(Match &mat, double mahalanobisTest,Cure::Pose3D & pose);
    int updateSensor(Match &mat, double mahalanobisTest, int sensorIndex);
    void updatePoses();
    void updatePose();

    void setRobotCov();
 };

} // namespace Cure
#endif
