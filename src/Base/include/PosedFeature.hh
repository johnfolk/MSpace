// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_POSEDFEATURE_H
#define CURE_POSEDFEATURE_H

#include "MapBank.hh"
#include "MapFeatureList.hh"
#include "MapFeature.hh" 
#include "Transformation3D.hh"
#include "ScanFragment.hh"
#include "Measurement.hh"

namespace Cure{

  // Forward declaration
  class Match;

  /**
   *  This class is used to represent the MapFeature in the sensor
   *  frame.  It facilitates matching, updates, initailization and
   *  M-space extentsion of the MapFeatures.  To do this it has
   *  members for the innovation of the feature measuremts and the
   *  jacobians of these innovations with respect to the feature
   *  coordinates and the sensor pose.  
   *
   * Summary of the notation behind the member names
   *
   * This is a MapFeature in the relative frame of a sensor.
   * Define the frames:
   *
   *  Sensor = S
   *  Robot  = R
   *  Map    = M
   * 
   * Define coordinates:
   * 
   *  Innovation  = E (eta)  
   *  Measurement = Z (gamma,rho,pixels, x,y ...)Generic for matching 
   *  Feature     = F (3D,2D,Scalars) in the map frame
   *  Sensor      = S (X of Tms, the transformation from M to S frames)
   *  Observed    = O (Xo=Tms(Xf), the feature coordinates in sensor frame
   *  P-space     = P (Xp=B*Xf, where B is the Projection matrix for the M-space)
   *
   *    dXf=Bdual*dXp  Projects changes in Xp to changes in Xf.
   *    dEta=Jeo*(Jos*(Jsr,Jss),Jof*Bdual) *(dXmr, dXrs, dXp)^T
   *
   *  where Jss is short for Jms rs (the jacobian of Xms wrt Xrs).
   *
   * we can write this a number of ways:
 
   *   dEta= (Jes,Jep)* (dXs, dXp)^T
   *   dEta= (Jeo*Jos,Jep)* (dXs, dXp)^T
   *
   * 1. So we get a measurement calculate z for it then do matching by calling 
   *     transform, roughtMatch and or fineMatch.  We can then do manhalobis 
   *     matching too if needed by calling linearizeMeasurement. 
   * 2. After doing matching we need to determing the measurement type 
   *    which depends on the measurment data and the map feature. 
   * 3. Then we can updateInfo and extend the feature.
   * 4. Then we call linearizeMeasurement to get the jacobians and 
   *    Eta we need for a P update.
   * 5. calculate the P update and call updateP on the map feature.  
   */
  class PosedFeature
  {
  public:
    MapBank *Bank;
    long FeatureKey;
    long ID;
    /**
     * (###########,Je,Eta,Jo,Xo)
     * Transform causes Xo and Z to be calculated.
     * CalcJo causes Jos, Jof and Jop to be calcualted.  
     * CalcEta normally causes Jeo and Jev to calculate  
     * CalcJe causes Jep and Jes to be calcualted.
     *
     * Xo --> Eta --> Jeo --> Jep.
     * Xo -> Jof --> Jop.
     * Xo -> Jos.
     * Jos,Jeo --> Jes.
     */
    unsigned short CalcState;
    /**
     * Transform(pose) causes Xo and Z to be calculated, CalcState=1.
     * One can then do rough and fine matching to the feature.
     */
    Cure::Matrix Xo;
    Cure::Matrix Z;//Used for rough and fine matching
    /**
     * calcJo(dR) causes Jos,Jof and Jop to be calculated, CalcState|=2..
     */
    Cure::Matrix Jos;
    Cure::Matrix Jof;
    Cure::Matrix Jop;

    /**
     * calcEta(v,mtype) causes Eta,Jeo and Jev to be calculated, CalcState|=4..
     */
    Cure::Matrix Eta;
    Cure::Matrix Jeo;
    Cure::Matrix Jev;

    /**
     * calcJe() causes Jep and Jes to be calculated, CalcState|=8..
     */
    Cure::Matrix Jep;
    Cure::Matrix Jes;

    Cure::Point3D **Points3D;
    Cure::Point2D **Points2D;

    /** @todo Comment */
    unsigned short MeasurementType;

    /** @todo Comment */
    unsigned short PoseType;

    /** 
     * Number of 3D and 2D points, number of scalars and the total full
     * dimension of this feature 
     */
    int Number3D,Number2D,NumberScalars,FullDim;
  public:
    PosedFeature(MapBank *b=0, long fKey=-1, long id=-1)
    {
      MeasurementType=0;
      PoseType=11;  // bitmask for x,y and theta (1+2+8)
      CalcState=0;
      ID=id;
      Bank=b;
      FeatureKey=fKey;
    }
    virtual ~PosedFeature(){}
    void equal(PosedFeature &pw);
    virtual PosedFeature * copy(){
      PosedFeature *pf=new PosedFeature();
      pf->equal(*this);
      return pf;
    }
    int getPoints(MapPoint ** points){
      MapFeature *mf =getMapFeature(Bank,FeatureKey);
      if (mf)
	{
	  points=mf->Points;
	  return mf->NumberPoints;
	}
      return 0;
    }
    MapFeature * getFeature(){
      return getMapFeature(Bank,FeatureKey);
    }


    long* getIndex()
    {
      MapFeature *mf =getMapFeature(Bank,FeatureKey);
      if (mf)
	return mf->getIndex();
      return 0;
    }

    int setIndex(long *ind)
    {
      MapFeature *mf =getMapFeature(Bank,FeatureKey);
      if (mf)mf->setIndex(ind);
      else return MAP_OBJECT_INVALID;
      return 0;
    }

    long  getFeatureID()
    {
      MapFeature *mf =getMapFeature(Bank,FeatureKey);
      if (mf)
	return mf->ID;
      return -2;
    }
    int setFeatureID(long id)
    {
      MapFeature *mf =getMapFeature(Bank,FeatureKey);
      if (mf)mf->ID=id;
      else return MAP_OBJECT_INVALID;
      return 0;
    }
    void setMeasurementType(unsigned short type)
    {
      type=getMeasurementType(type);
      if (MeasurementType!=type)
	CalcState=(CalcState&0x3);
      MeasurementType=type;
    }
    /*
     * Finds the best measurment type for this type of measurment data and the 
     * current infomation  collected on the map feature.
     * Meaurement type 0 is no possible linearized measurment.  This is for 
     *  uninitailized features. 
     */
    virtual unsigned short getMeasurementType(unsigned short type)
    {
      MapFeature *mf =getFeature();
      if (mf) return mf->getMeasurementType(type);
      return 0;  
    }
    virtual void narrow(){}
    /**
     * IMPORTANT:
     * This must be called before any other methods each time the
     * sensor frame changes. 
     *
     * @param trans, the map to sensor transformation, Tms. 
     * @return NOT_VISABLE; if not in image plane or behind wall ect. 
     * of  MAP_OBJECT_INVALID if FeatureKey is not good.
     */
    int transform(Cure::Pose3D& pose){
      return transform(pose,pose.getCovType());
    }
    /**
     * IMPORTANT:
     * This must be called before any other methods each time the
     * sensor frame changes. 
     *
     * @param trans, the map to sensor transformation, Tms. 
     * @param covType the Covariance Type of the sensor pose used to
     * determine which coordinates to include in the jacobians wrt the pose.
     * @return NOT_VISABLE; if not in image plane or behind wall ect. 
     * of  MAP_OBJECT_INVALID if FeatureKey is not good.
     */
    virtual int transform(Cure::Transformation3D& pose, unsigned short covType)
    {
      PoseType=covType;
      return calcXo(pose);
    }
    /**
     * An initial elimination of most of the features from the map.
     * This feature is compared to the MeasurmentType specific z using
     * the thresholds to see if this feature could possibly match.
     * @param z The sensor/feature dependent measurment as an array.
     * @param Thresholds Array of max 'errors' in measurements compared to
     * the predicted values. 
     * @return 0 for a match  
     */
    virtual int roughMatch(double  *z,double *thresholds);
    /**
     * Use this to check if the feature is in some rectangle in the
     * sensor frame.
'    * @param rect The corners of the rectangle (minx,miny, maxx,maxy) 
     *  @return 1 if inside
     */
    virtual int inRectangle(double rect[4]){return 0;}
 
    /**
     * A finer test of the features from the map.
     * This feature is compared to the MeasurmentType specific z using
     * the metric sigma and the FineThreshold.
     * @param z The measurement to be matched should be given here, 
     * the predicted measurement is subtracted from this.
     * @param sigma the metric for the error in z
     * @return -1 if fails else the distance estimated using sigma.
     */
    virtual double fineMatch(double * z,Cure::Matrix & sigma);

 
    /**
     * A finer test of the features from the map.
     * This feature is compared to the MeasurmentType specific z using
     * the metric sigma and the FineThreshold.
     * @param m The measurement to be matched should be given here, 
     * the predicted measurement is subtracted from this.
     * @param sigma the metric for the error in z
     * @return -1 if fails else the dist@ance estimated using sigma.
     */
    virtual double fineMatch(Cure::Measurement &m,Cure::Matrix & sigma){
      return fineMatch(m.Z.Element,sigma);
    }

    virtual int consistentPose(Cure::Measurement &m, Cure::Transformation3D &pose){return 1;}
    /**
     * Calculates Eta, Jeo and Jev.  
     *
     * @param v The measurment vector.
     * @param mtype Specifies definiation of Eta
     * @param noJac set to 1 to stop the calc after Eta.
     * @return 1 if fails 0 on success. 
     */
    virtual int calcEta(Cure::Matrix & v,unsigned short mtype, unsigned short noJac=0){
      if (!(CalcState&1))return 1;
      setMeasurementType(mtype);
      if (CalcState&4)return 0;
      return 1;
    } 
    /**
     * calculates Jos, Jof and Jop.
     */
    int calcJo(Cure::Matrix dR[10]);

    /**
     * Calculates Jep and  Jes
     */
    int calcJe();

    /**
     * Transforms the feature and calculates the inovation eta and 
     * all jacobians.
     * @param v The measurement is used to calculate eta and jacobians
     * @param pose The sensor pose.
     * @param mtype Measument type.
     */
    int linearizeMeasurement(Cure::Matrix  & v,
			     Cure::Pose3D & pose,
			     unsigned short mtype=0);
    /**
     * Adds dp to the p part of the MapFeature'scoordinates 
     * for i = 0; i<Bdual.Columns.  
     * It then sets dp=0. It also recenters if needed.
     * @param dp Column matrix of the change in p-space dimension.
     */
    void updateP(Cure::Matrix & dp)
    {
      MapFeature *mf =getMapFeature(Bank,FeatureKey);
      if (mf) mf->updateP(dp);
    }
    /**
     * Adds dxp to the p part of the coordinates and then sets dxp[i] to 0
     * for i = 0; i<B.Columns.
     * @param dxp Array of the change in p-space dimension.
     * @return number of P-space dimensions 
     *         MAP_OBJECT_INVALID for invalid map feature.
     */
    int updateP(double * dxp)
    {
      MapFeature *mf =getMapFeature(Bank,FeatureKey);
      if (mf) return mf->updateP(dxp);
      return MAP_OBJECT_INVALID;
    }
    int updateIndexedP(double *dxp)
    {
      MapFeature *mf =getMapFeature(Bank,FeatureKey);
      if (mf) return mf->updateIndexedP(dxp);
      return MAP_OBJECT_INVALID;
    }
    /**
     * Adds the new collected Info to the cumulation and updates the 
     * LowInfo dimensions and then calls transform(pose).
     * @param mat Contains the Info to add
     * @param pose the pose of the sensor when taking the measurement
     * @param type Meaning can vary depending on subclass.
     * @return 0 if ok
     * -1 if feature is already deleted.
     */
    virtual int  addInfo(Match &mat,
			 Cure::Pose3D & pose,
			 const int type=0){
      Cure::Transformation3D map2Info;
      return addInfo(mat,pose,pose.getCovType(),map2Info,type);
    }
    /**
     * Adds the new collected Info to the cumulation and updates the 
     * LowInfo dimensions and then calls transform(pose).
     * @param mat Contains the Info to add
     * @param pose the current pose of the sensor for this PosedFeature 
     *        to be transformed to aat the end. This is to set up things 
     *        like Eta, Jes and so on.  So it is the CovType and 
     *        Transformation3D information that is used.  The actual Covariance
     *        need not be correct.
     * @param map2info The continuous information frame relative to the
     *                 fixed/global map frame.  These can be the same,
     *                 map2info=I.  Better if this takes one to the dead
     *                 reckoning frame of the robot.
     * @param type Meaning can vary depending on subclass.
     * @return 0 if ok
     * -1 if feature is already deleted.
     */
    virtual int  addInfo(Match &mat,
			 Cure::Pose3D & pose,Cure::Transformation3D & map2info,
			 const int type=0){
      return addInfo(mat,pose,pose.getCovType(),map2info,type);
    }
    /**
     * Adds the new collected Info to the cumulation and updates the 
     * LowInfo dimensions and then calls transform(pose).
     * @param mat Contains the Info to add
     * @param pose the current pose of the sensor for this PosedFeature 
     *        to be transformed to aat the end. This is to set up things 
     *        like Eta, Jes and so on.  
     * @param covType The CovType for Jes Columns.
     * @param map2info The continuous information frame relative to the
     *                 fixed/global map frame.  These can be the same,
     *                 map2info=I.  Better if this takes one to the dead
     *                 reckoning frame of the robot.
     * @param type Meaning can vary depending on subclass.
     * @return 0 if ok
     * -1 if feature is already deleted.
     */
    virtual int  addInfo(Match &mat,
			 Cure::Transformation3D & pose,unsigned short covType,
			 Cure::Transformation3D & map2info,
			 const int type=0){return 0;}

    /*
     * @param t the transform to the sensor frame
     * @param jaces the jacobian of the pseudo measuremnent
     * @param invcov Returns with the (covariance)^-1 of the psuedo measurement.
     * @return number of new dimensions, 
     */
    virtual  int  extend(Cure::Transformation3D & t, 
			 Cure::Matrix & jaces, 
			 Cure::Matrix & invcov);

    /**
     * @param s put Scan index range pairs in s from start to end
     * @return NOT_VISABLE if not scanable (behind wall ect.)
     *     else 0 and 
     */ 
    virtual int scan(ScanFragment& s){
      return NOT_VISABLE;
    }
  
    virtual void print(int level=0);
    int setXo(); 
  protected:
    int calcXo(Cure::Transformation3D & trans); 
  };
}
#endif
