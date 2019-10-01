
// = AUTHOR(S)
//    John Folkesson
//    
//    August 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_POINTFEATUREHELPER_HH
#define CURE_POINTFEATUREHELPER_HH

#include "MapPointFeature.hh"
#include "PosedPointFeature.hh"
#include "MapHelper.hh"
#include "Match.hh"

namespace Cure {

  /** 
   * Helper class for PointFeatures
   *
   * The structure of the Measurements are assumed to be 
   * SensorType=SensorData::SENSORTYPE_CAMERA=2
   * MeasurementType=32
   * V: 3x1 (x,y,f)
   * Z: 2x1 (x,y)
   * W: 1x10 (imwidth, imheight, f, pixwidth, pixheight, 
   *          x_pix, y_pix, xsize_pix, ysize_pix, angle)
   */
  class PointFeatureHelper: public MapHelper 
  {
  public:
    PosedPointFeature * PosedPointFeaturePtr;
    MapPointFeature * MapPointFeaturePtr;
    MapPointFeature TempletFeature;
    Cure::Point3D Center;
    double DistanceGuess;
    /** 
     * The Templet will be given these thresholds which then are used to 
     * initilize the MapFeatures made with this helper. 
     */
    PointFeatureHelper(MapBank *b, double triThreshold=.05,
		       double wgtThreshold=1,
		       double distanceThreshold=20,double trackThreshold=.1,
		       double distanceGuess=2.14):
      MapHelper(b){
      TempletFeature.setBank(b);
      TempletFeature.CastPtr=&MapPointFeaturePtr;
      TempletFeature.TriangleThreshold=triThreshold;
      TempletFeature.WeightThreshold=wgtThreshold;
      TempletFeature.DistanceThreshold=distanceThreshold; 
      TempletFeature.TrackThreshold=trackThreshold; 
      DistanceGuess=distanceGuess;
      FeatureType=MAPPOINTFEATURE_TYPE;
      ObjectSubType=FeatureType;
    }
    ~PointFeatureHelper(){}
    virtual void printConfiguration();

    PosedPointFeature * castPosedPointFeature(PosedFeature *pf){
      PosedPointFeaturePtr=0;
      pf->narrow();
      return PosedPointFeaturePtr;
    }
    MapPointFeature * castMapPointFeature(MapObject *pf){
      MapPointFeaturePtr=0;
      pf->narrow();
      return MapPointFeaturePtr;
    }
    MapPointFeature *getMapPointFeature(long key)
    {
      MapObject *mo=Bank->getMapObject(key);
      if (!mo) return 0;
      return castMapPointFeature(mo);
    } 
    virtual MapFeature *getMapFeature(long key){
      return getMapPointFeature(key);
    }

    virtual MapFeature * makeMapFeature(MapBank *b=0){
      return new MapPointFeature(&TempletFeature,b);
    }

    /**
     * @return A MapPointFeature based on the TempletFeature. 
     */
    MapPointFeature * makeMapPointFeature(){
      return new MapPointFeature(&TempletFeature,Bank);
    }
    /**
     * 
     * @param sensorpose The pose of the sensor with the y axis ahead, 
     *        x is right and z up.
     * @param x The coordinates of the center 
     * in meter units relative to sensor.
     * @return A MapPointFeature based on the TempletFeature and the params. 
     */
    MapPointFeature * makeMapPointFeature(Cure::Pose3D sensorpose,
					  const double imagex[2],
					  double distanceguess=0) 
    {
      if (distanceguess==0)distanceguess=DistanceGuess;
      MapPointFeature *temppt=new MapPointFeature(&TempletFeature,Bank);
      temppt->initializeFromPixels(sensorpose,imagex,   PixelInfo[0], distanceguess);
      return temppt;
    }     

    virtual PosedFeature * makePosedFeature(MapFeature* mf){
      return makePosedPointFeature(PixelInfo[0],mf);}
    virtual MapFeature * makeMapFeature(Cure::Pose3D &sensorpose,
					Cure::Measurement &m,
					bool addToVisable=false);

    //  virtual MapFeature * makeMapFeature(Cure::Pose3D sensorpose,
    //			      double distanceguess=2.0){
    //  return makeMapPointFeature(sensorpose,distanceguess);}
    /**
     * 
     * @param Sensorpose The pose of the sensor with the y axis along the 
     *        camera axis, x is right and z up.
     * @return A MapPointFeature based on the Center, 
     * TempletFeature and the params. 
     */
    MapPointFeature * makeMapPointFeature(Cure::Pose3D & sensorpose,
					  double distanceguess=2.0) 
    {
      MapPointFeature *temppt=new MapPointFeature(&TempletFeature,Bank);
      temppt->initializeFromPixels(sensorpose,Center.X, PixelInfo[0],
				   distanceguess);
      return temppt;
    }     
    PosedPointFeature * makePosedPointFeature(double focal,
					      MapFeature* mf){
      if(mf->Type==MAPPOINTFEATURE_TYPE)
	return new PosedPointFeature(&PosedPointFeaturePtr,focal,
				     castMapPointFeature(mf));
      return 0;
    }
    virtual bool getC(MapFeature  *mf1,MapFeature  *mf2,
		      Cure::Matrix &c1){
      MapPointFeature *mw1=castMapPointFeature(mf1);
      if (!mw1)return false;
      MapPointFeature *mw2=castMapPointFeature(mf2);
      if (!mw2)return false;
      return mw1->getC(mw2,c1);
    }
    virtual bool testMatch(MapFeature  *mf1,MapFeature  *mf2,
			   double tolerance){
      MapPointFeature *mw1=castMapPointFeature(mf1);
      if (!mw1)return false;
      MapPointFeature *mw2=castMapPointFeature(mf2);
      if (!mw2)return false;
      return mw1->testMatch(mw2,tolerance);
    }

    virtual int config(const std::string &arglist);

    virtual bool supportsSubconfig(int sc);
  
    /**
     * @param mat the Match is returned in this.
     * @param v The (1x2) coordinates of the center in the image
     * @param distance The diatance to add info with.
     */
    int addDistance(Match &mat,double distance);

    /*
     * try finding a merge candidates from search list
     * If found remove them from searchlist and return them in 
     * pl[0] and pl[1].  pl[0] is the one to delete.
     * @return  0  for no matches;
     *          else the first 3 bits give pl[0]'s P-dim
     *          0x1=1   pl[0] has 0 p-dim
     *          0x2=2   pl[0] has 1 p-dim
     *          0x4=4   pl[0] has 3 p-dim
     *          and the next 3 bits give pl[1]'s P-dim
     *          0x8=8   pl[1] has 0 p-dim
     *          0x10=16  pl[1] has 1 p-dim
     *          0x20=32  pl[1] has 3 p-dim
     *          0x40=64  if start and endpoints are switched.
     */
    unsigned short findMergeFeatures(PosedFeature * pfp[2],
				     PosedFeatureList * searchList);

    /**
     * constraint is a*(pl[0].P,pl[1].P)=b 
     * dX=-K(ap-b)=Ke;
     * C'=(I-Ka)C
     * K=Ca^T(aCa^T)^-1
     * So the measurement is e with jacobian Jep of a and covaraiance cov. 
     * @param typ the return from findMergeFeatures
     */
    int getMergeContstraint(Cure::Matrix &a, Cure::Matrix &e,
			     Cure::Matrix &cov, double distance,
			     unsigned short typ, PosedFeature * pl[2]);
    virtual int merge(PosedFeature *pf[2], unsigned short typ=1){
      MapPointFeature *mf0=getMapPointFeature(pf[0]->FeatureKey);
      MapPointFeature *mf1=getMapPointFeature(pf[1]->FeatureKey);
      if (!mf1)return -1;
      if (!mf0)return -1;
      int r=mf1->merge(mf0,typ);  
      if (typ&1){
	delete mf0;
      }
      return r;
    }
    static void 
    makeBearings(Cure::Matrix &bearings,bool *frames,
		 Cure::Matrix &v,
		 Cure::Transformation3D  *sensorposes,
		 const double distanceguess=5);



    /**
     * Precalculate the path and deltapath for triangulatePoint. This
     * function is typically used before calling triangulatePoint.
     * Here we are updating the matricies from the triangulation of the 
     * previous frame buffer.  So we remove the oldest (index) and 
     * replace it with the new pose.  This is for a circular buffer
     * that avoid copying.
     *
     * The function will grow the path and deltapath matricies if
     * index is beyont their size and replace if index is within them.
     *
     * @param path the (x,y,z) of camera when bearings taken in 
     *             the info frame (output)
     * @param deltapath the differences in (x,y,z) of camera (output)
     * @param sensorposes the new camera pose
     * @param index the number of the new frame. 
     */
    static void 
    insertPoseTopath(Cure::Matrix &path, 
		     Cure::Matrix &deltapath,
		     const Transformation3D  &sensorpose,
		     int index);
    /**
     * Calculate the predicted image plane coordinates of the meanpoint.
     * @param vhat returns with (x,z,focallength) for all poses that have
     *        vhat(i,2)=focallength>0 when called.  Use the third column
     *        to indicate which poses to use.
     * @param sensorposes An array [n] of the camera poses in the info frame.
     * @param meanpoint the (xyz) of the point in the world frame.
     */
    static void 
    predictPixels(Cure::Matrix &vhat,
		  Cure::Transformation3D  *sensorposes,
		  Cure::Matrix &meanpoint);

    /**
     * This creats a time ordered matrix of innovations with time in 
     * whole frames from oldest frame assuming a cirular buffer.
     * The vhat(i,2)>0 indicates which frames have data.
     * @param innovation (v(i,0)-vhat(i,0), v(i,1)-vhat(i,1),time in frames).
     * @param oldestframe The index of oldest frame 
     */
    static void 
    calcInnovations(Cure::Matrix &innovation,
		    Cure::Matrix &vhat,
		    Cure::Matrix &v, int oldestframe);
    /**
     * @param dx_dt (x(i+1,0)-x(i,0)/x(i+1,2)-x(i,2),
     *               x(i+1,0)-x(i,0)/x(i+1,2)-x(i,2),
     *               x(i+1,2)+x(i,2)/2,
     */
    static void calcTimeDerivative(Cure::Matrix &dx_dt, Cure::Matrix &x);

    /**
     * Triangulates a set of point observations. Use this function to
     * find an estimate of a point feature given a number of
     * bearing-only measurements to it. The measurements are in the form
     * of camera poses, image coordinates and focal length, one for each
     * frame where the feature has been detected. The output is the
     * meanpoint which is the estimate of the point and the return value
     * which gives a measure of how accurate the triangulation.
     *
     * The path and delta can be calculated for a set of frames using calcPath 
     * or growPath.  Then these same matricies can be reused for repeated calls
     * to this function.
     *
     * So if frames[i]==true then bearings(i,k) has a 
     * bearing to the point at path(i,k). 
     *
     * You can pass on pointers to variables to get hold of the value
     * for the different parameters that are checked against the
     * thresholds. This gives you a way of "tuning" by looking at what
     * triangulation points are good and then looking at what values
     * these
     * 
     * @param meanpoint returns the (3,1) column of the x,y,z of the point 
     *                  in the info frame (putput)
     * @param path the matrix (n,3) of camera (x,y,z) positions
     * @param deltapath the (n,3*n) matrix of camera displacements
     * @param bearings the (n,3) matrix of measurements (u,v,focallength)
     * @param weightThreshold used to test sum of sin^2(theta),angle 
     *                       between bearings.
     * @param sqDistanceThreshold  used to test distance of each bearing 
     *                             ray from the meanpoint.
     * @param mseThreshold max mean square error of all the bearings 
     *                     from the meanpooint.
     *
     * @param retWeight pointer to storage variable of sqDist used to
     * compare with the threshold weightThreshold
     * @param retSqDist pointer to storage variable to store max sqDist
     * @param retMse pointer to storage variable of mse used to
     * compare with the threshold mseThreshold
     * @param retSecDer pointer to storage variable of secDer used to
     * compare with the threshold secDerThreshold
     *
     * @return -1 if no baseline, sum sin^2 angle < threshold
     *         -2 if some bearing's best point was further than 
     *            consistencytest from meanpoint
     *         -3 if mse >consistencytest/sqrt(n) mse=mean square error of 
     *            point cloud from meanpoint
     *         else sum(sin^2(theta)) triangulation quality 
     *     
     */
    static double 
    triangulateBuffer(Cure::Matrix &meanpoint,
		      Cure::Matrix &path, 
		      Cure::Matrix &deltapath, 
		      Cure::Matrix &bearings,
		      bool *frames,
		      const double weightThreshold=.1, 
		      const double sqDistanceThreshold=.01,
		      const double mseThreshold=.001,
                      double *retWeight = 0,
                      double *sqDist = 0,
                      double *mse = 0);

    /**
     * Triangulates a set of point observations. Use this function to
     * find an estimate of a point feature given a number of
     * bearing-only measurements to it. The measurements are in the form
     * of camera poses, image coordinates and focal length, one for each
     * frame where the feature has been detected. The output is the
     * meanpoint which is the estimate of the point and the return value
     * which gives a measure of how accurate the triangulation.
     *
     * Here we calculate pairwise and then take a weighted aveerage.
     * The complexity is n^2 as opposed to n^3 for the min distance form.
     *
     * The measurement matrix v can have rows of 0's if so the point is 
     * ignored (v(i,2)==0, is the condition to ignore this sensorpose. 
     *
     * Thus one can have a fixed, (even circular) buffer of poses and 
     * just be sure that the rows of v(i,j) coorespond to the sensorpose[i].
     * Unmatched frames will then have v(i,2)=0; and are ignored.
     *
     * Upon return all frames not used will have v(i,2)<=0. 
     * The orginal empty rows will be unchanged (ie.0)
     * The detected outliers will have v(i,2)==-1.
     *
     * That is to say some of the bearings will not pass the 
     * sqdistanceThreshold, (not point at the meanpoint) and will not
     * be used.  They haave v(i,2) set to -1.  THis information can also 
     * be used by the calling object.
     *
     * You can pass on pointers to variables to get hold of the value
     * for the different parameters that are checked against the
     * thresholds. This gives you a way of "tuning" by looking at what
     * triangulation points are good and then looking at what values
     * these
     * 
     * @param meanpoint returns the (3,1) column of the x,y,z of the point 
     *                  in the info frame (putput)
     * @param path the matrix (n,3) of camera (x,y,z) positions
     * @param deltapath the (n,3*n) matrix of camera displacements  
     * @param v the (n,3) matrix of measurements (u,v,focallength)
     * @param sensorposes An array [n] of the camera poses in the info frame.
     * @param weightThreshold used to test sum of sin^2(theta),angle 
     *                       between bearings.
     * @param sqDistanceThreshold  used to test distance of each bearing 
     *                             ray from the meanpoint.
     * @param mseThreshold max mean square error of all the bearings 
     *                     from the meanpooint.
     * @param secDerThreshold threshold on the second derivateive if <=0 
     *                     this test will be ignored
     * @param distanceguess default 5 (units of v)this is used to get bearings 
     *
     * @param retWeight pointer to storage variable of sqDist used to
     * compare with the threshold weightThreshold
     * @param retSqDist pointer to storage variable to store max sqDist
     * @param retMse pointer to storage variable of mse used to
     * compare with the threshold mseThreshold
     * @param retSecDer pointer to storage variable of secDer used to
     * compare with the threshold secDerThreshold
     * 
     * @return -1 if no baseline, sum sin^2 angle < threshold
     *         -2 if some bearing's best point was further than 
     *            consistencytest from meanpoint
     *         -3 if mse >consistencytest/sqrt(n) mse=mean square error of 
     *            point cloud from meanpoint
     *         -4 if second derivative is two large 
     *         else sum(sin^2(theta)) triangulation quality 
    */
    static double 
      triangulateBuffer(Cure::Matrix &meanpoint,
			Cure::Matrix &path, 
			Cure::Matrix &deltapath, 
			Cure::Matrix &v,
			Cure::Transformation3D  *sensorposes,
			const double weightThreshold=.1, 
			const double sqDistanceThreshold=.01,
			const double mseThreshold=.001,
                        const double secDerThreshold=1e-4,
			const double distanceguess=5,
                        double *retWeight = 0,
                        double *retSqDist = 0,
                        double *retMse = 0,
                        double *retSecDer = 0);
    
    
    /**
     * Same as the above function but requires only the sensorposes
     * array. The path and delatpath are derived automatically from
     * this. This is of course more inefficient that calculating the
     * path and deltapath once and for all.
     *
     * Triangulates a set of point observations. Use this function to
     * find an estimate of a point feature given a number of
     * bearing-only measurements to it. The measurements are in the form
     * of camera poses, image coordinates and focal length, one for each
     * frame where the feature has been detected. The output is the
     * meanpoint which is the estimate of the point and the return value
     * which gives a measure of how accurate the triangulation.
     *
     * The measurement matrix v can have rows of 0's if so the point is 
     * ignored (v(i,2)==0, is the condition to ignore this sensorpose. 
     *
     * Thus one can have a fixed, (even circular) buffer of poses and 
     * just be sure that the rows of v(i,j) coorespond to the sensorpose[i].
     * Unmatched frames will then have v(i,2)<=0; and are ignored.
     *
     * Upon return all frames not used will have v(i,2)=0.
     * That is to say some of the bearings will not pass the 
     * sqdistanceThreshold, (not point at the meanpoint) and will not
     * be used.  They haave v(i,2) set to -1.  THis information can also 
     * be used by the calling object.
     *
     * You can pass on pointers to variables to get hold of the value
     * for the different parameters that are checked against the
     * thresholds. This gives you a way of "tuning" by looking at what
     * triangulation points are good and then looking at what values
     * these
     * 
     * @param meanpoint returns the (3,1) column of the x,y,z of the point 
     *                  in the info frame (putput)
     * @param v the (n,3) matrix of measurements (u,v,focallength)
     * @param sensorposes An array [n] of the camera poses in the info frame.
     * @param weightThreshold used to test sum of sin^2(theta),angle 
     *                       between bearings.
     * @param sqDistanceThreshold  used to test distance of each bearing 
     *                             ray from the meanpoint.
     * @param mseThreshold max mean square error of all the bearings 
     *                     from the meanpooint.
     * @param secDerThreshold threshold on the second derivateive if <=0 
     *                     this test will be ignored
     * @param distanceguess default 5 (units of v)this is used to get bearings 
     * 
     * @return -1 if no baseline, sum sin^2 angle < threshold
     *         -2 if some bearing's best point was further than 
     *            consistencytest from meanpoint
     *         -3 if mse >consistencytest/sqrt(n) mse=mean square error of 
     *            point cloud from meanpoint
     *         -4 if second derivative is two large 
     *         else sum(sin^2(theta)) triangulation quality 
    */
   static double 
   triangulateBuffer(Cure::Matrix &meanpoint,
		     Cure::Matrix &v,
		     Cure::Transformation3D  *sensorposes,
		     const double weightThreshold=.1, 
		     const double sqDistanceThreshold=.01,
		     const double mseThreshold=.001,
                     const double secDerThreshold=1e-4,
		     const double distanceguess=5);

  protected:

    static double makeDerivateCheck(Cure::Matrix &meanpoint,
                                    Cure::Matrix &v,
                                    Cure::Transformation3D *sensorposes);
    
    /**
     * Reads configurations of version 1
     */
    int configVer1(const std::string &arglist);
  };

} // namespace Cure

#endif
