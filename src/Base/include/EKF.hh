// = AUTHOR(S)
//    John Folkesson
//    
//    August, 2004
 //
//    Copyright (c) 2004 John Folkesson
//    

#ifndef CURE_EKF_HH
#define CURE_EKF_HH




#include <iostream>
#include <stdlib.h> 
#include "PosedFeature.hh"
#include "Measurement.hh"
#include "Match.hh"
#include "MapFeature.hh"
#include "MapBank.hh"
#include "Matrix.hh" 
#include "Pose3D.hh"
#include "FeatureCovarianceInterface.hh"

namespace Cure{

  /**
   * A General SLAM EKF
   *
   * Initialization info:
   * 1. call initRobotPose 
   * 2. call initSensorOffset
   * 
   * @author John Folkesson 
   */
  class EKF : public FeatureCovarianceInterface
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
     * The first row for the MapFeature dimensions
     *
     * @todo given that RobotRows gives the number of rows with
     * features, maybe FeatureStartRow is a better name?
     */
    int FeatureRows;

    /**
     * An array of MapFeature Keys.
     */
     LongList *Keys;

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
    /** @todo Comment difference between SensorCovTypes and OffsetCovTypes */
    unsigned short SensorCovTypes[5];
    unsigned short OffsetCovTypes[5];

    /** Type of covariance for robot state variables */
    unsigned short RobotCovType;

    /** Number of rows (state variables) for robot pose */
    int RobotRows;

    /** @todo is it even used anymore? */
    double *AllocatedC;

  public:
    EKF(MapBank *b,int dim);
    ~EKF();
    void initRobotPose(Cure::Pose3D &pose);
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
    int initVelocity(unsigned short veltype);
    void reallocate(int dim);
    /**
     * Removes dim rows from Keys and C begining with row=start.
     */
    int removeRows(int start,int dim);
    /**
     * Checks that the MapFeatures have not been deleted for all 
     * rows including and above 'start'. 
     */
    void checkKeys(int start);
    /**
     * This checks if some of the features at Keys[i] (for i= ind, ind+1, ind+2)
     * really have Index[k]=i for some k<P-dim.  
     * If not (the point has been merged) the rows of C are deleted.
     * @param ind the row to start checking from.  
     */  
    void checkIndex(long ind);
    MapFeature * getFeature(int i);
    /**
     * does pose=pose+inc and updates C.
     * @param pose the Pose state before the update in and the updated pose out
     * @param inc the incremental change in Pose with uncertainty.
     */ 
    void incrementalPredict(Cure::Pose3D & pose, Cure::Pose3D & inc); 
    /**
     * This does offset[sensorIndex]=offsetpose and updates C.
     * This assumes that the new offset is independent of the old one.
     * That is ok sometimes like a typical pan/tilt uncertainty with 
     * known static offset.
     *
     * Call this if you have a new better estimate of the offset between
     * the robot and sensor frames.  If you have no new info do no call this
     * and the old offset will be used which has been adjusted by the 
     * kalman earlier updates.
     *
     * @param offsetpose The new sensor offset with its uncertainty 
     * @param sensorIndex the index into Offsets, returned by initSensorOffset.
     */ 
    void offsetPredict(Cure::Pose3D & offsetpose, int sensorIndex); 

    /**
     *  This will do a more general predict using a second order
     *  model.  
     *
     *  x(k)=x(k-1)+v(k-1)dt+a(k)(dt*dt/2).
     *  Euler(k)=euler(k-1)+w(k)dt.
     *  v(k)=v(k-1)+a(k)dt.
     *
     *  here the state shoule be the x(k)=(x,y,z), Euler(k)=(theta,phi,psi)
     *  and v(k)=(vx, vy, vz).
     *
     *  The linear acceleration a(k) is tranformed from the robot
     *  frame measurement in accPart to earth frame using the state at k-1.
     *
     *  The angular velocity w(k) is transformed from the the robot
     *  frame measurement in incPart's Euler/dt using the state at k-1.
     *
     *  So the incPart is the incremental Pose over the interval with
     *  uncertianty.  It should have CovType 48.  It is then used to 
     *  infer the angulare velocities and uncertianties by division
     *  with dt.  
     * 
     * @param pose the updated pose out
     * @param inc the incremental change in Pose with uncertainty.
     * @param accPart The first row of Data is the acceleration while
     * the rows beneth the first give the covaraince in the
     * accelleration.
     */
    void incrementalPredict(Pose3D & pose, Pose3D & incPart, 
			    GenericData &accPart);
    /**
     * This returns the current estimate of the sensor pose.
     * @psensor the sensorpose is returned here
     * @sensorIndex This selects which sensor you want.
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
     * @param mat contains measurments<->PosedFeatures
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
                   long *** index,  int numblocks,int numJac);

    /**
     * This function goes through all the features and propagates
     * changes to the P-space coordinates to changes in the feature
     * coordinates. It also removes space in the Dx and C
     * corresponding to
     */
    void updateP();

    int checkEigen();
 
    /**
     * constraint is a*(pl[0].P,pl[1].P)=b 
     * dX=-K(aP-b)=-Ke;
     * e=a*(pl[0].P,pl[1].P)-b 
     * C'=(I-Ka)C
     * K=Ca^T(aCa^T)^-1
     * 
     */
    int merge(PosedFeature *pl[2],Cure::Matrix &a,Cure::Matrix &e,
	      Cure::Matrix & cov, double mahalanobisTest,
	      Cure::Pose3D & pose);	
    void permuteRows(Cure::Matrix & permute,PosedFeature **pl,long number);
    void permuteRows(Cure::Matrix & permute,long *index);


    /**
     * Propagates chnages in P-space to changes in the feature
     * coordinates of a certain feature.
     */
    int updateP(PosedFeature & pf);

    /**
     * Checks if enough dense information has been collected to grow
     * the P-space to more dimensions. For example to add an endpoint
     * to a wall line feature.
     */
    int extend(Cure::Pose3D &pose, Cure::PosedFeature &pf);

    int extend(Pose3D &pose, PosedFeature &pf, int sensorIndex){ 
      int r= extend(pf,sensorIndex);
      pose=Robot;
      return r;
    }

    /**
     * Get the covariance of a feature's M-Space.
     * @param returns witht the Covariance of just this feature's
     *        M-Space.
     * @param mf the MapFeature to find the covaraince of
     * @return M-Space dimension of mf 
     */
    int getFeatureCov(Matrix &fCov,MapFeature *mf);


    /**
     * Get the correlation between a feature and the robot
     * @param corr correlation matrix
     * @param mf the MapFeature for which to find correlation to robot for
     * @return M-space dimension of mf
     */
    int getFeatureRobotCorr(Matrix &corr, MapFeature *mf);

  protected:
    /**
     * This updates the map with the measurments and adds dense 
     * info in a given frame
     * @param mat contains measurments<->PosedFeatures
     * @param pose the pose state of the measurement sensor.
     * @param map2info the transformation to take the dense info
     * from the map to the info frame.
     */
    int addinfo(Match &mat,
		Cure::Transformation3D & map2info){
      if (mat.MatchedFeature)
	{
	  updateP(*mat.MatchedFeature);
	  return  mat.MatchedFeature->addInfo(mat,Robot,map2info, 
					      mat.Measure->MeasurementType);
	}
      return 1;
    }
 
    int addinfo(Match &mat,
		unsigned short covType,
		Cure::Transformation3D & map2info){
      if (mat.MatchedFeature)
	{
	  updateP(*mat.MatchedFeature);
	  return  mat.MatchedFeature->addInfo(mat,Robot,covType,map2info, 
					      mat.Measure->MeasurementType);
	}
      return 1;
    }
    int addinfo(Match &mat,Cure::Transformation3D & sensorpose,
		unsigned short covType,
		Cure::Transformation3D & map2info){
      if (mat.MatchedFeature)
	{
	  updateP(*mat.MatchedFeature);
	  return  mat.MatchedFeature->addInfo(mat,sensorpose,covType,map2info, 
					      mat.Measure->MeasurementType);
	}
      return 1;
    }
  
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

    /**
     * Adds blank rows to C and Dx.  This moves a number of rows equal to
     * currentDim and indexed by index to the end of C.
     * Then it adds addedDim rows to it with C=10000.
     * 
     */
    void insertRows(int currentDim,int addedDim, long *index,long featureKey);
    int extend(Match &mat);
    int extend(PosedFeature &pf);
    int extend(PosedFeature &pf,int sensorIndex);
    int extend(Match &mat, int sensorIndex)
    { 
      PosedFeature &pf=*mat.MatchedFeature;
      int dim=extend(pf, sensorIndex);
      return dim;
    }

    void setRobotCov();
 };

} // namespace Cure
#endif
