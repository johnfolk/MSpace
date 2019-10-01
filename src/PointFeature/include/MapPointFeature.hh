// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_MAPPOINTFEATURE_HH
#define CURE_MAPPOINTFEATURE_HH


#include "MapFeature.hh"
#include "Base/FeatureFun.hh"
#include "MatrixStuff.hh"
#include "LinkedArray.hh"

namespace Cure {

  /**
   * Class managing point features in the map
   * 
   * @author John Folkesson
   */
  class MapPointFeature: public MapFeature
  {
  public:
    MapPointFeature **CastPtr;
    long AllocatedIndex[3];
    MapPoint * AllocatedPoints[1];
    double TotalWgt;
    /**
     * The 'weights' are sumed for each pair of bearing vectors that 
     * contribute to the weighted triangulation.  These weights are
     * the product of the weight supplied with each bearing vector
     * and the sin^2(angle between the Vectors).
     * This sum is then compared to WeightThreshold and if it is > the 
     * High dim are expanded.
     */
    double WeightThreshold;
    /**
     * The min distance in meters, between 2 camera positions for including in
     * the wieghted sum of triangulation points.
     */
    double TriangleThreshold;
    /**
     * This is a threshold on 1-cos(dtheta) where dtheta is the change in
     * angle of the vector from the camera towards the point.
     */
    double TrackThreshold;

    int VCount;
    /**
     * Vectors.Element's are double arrays with this structure:
     * distance, w, x_i, s_i, s*x_i, r_i, sum_w  (0..10)
     * distance w x_i and s_i are given when calling addInfo.
     * w is a relative weight estimated to approximate the amount of
     * infomation contained in this observation (inverse variance)
     * x_i is the xyz of the camera when the observation was made.
     * s_i is the bearing vector from x_i towards the point (magnitude ignored)
     * s*xi is the dot product of the 2 previous numbers.
     * r_i is the weighted sum of the distances along s_i to the 
     * triangulated intersection with later Vector.Elements.
     * sum_w is the total of the weights in r_i so that
     * the estimate for the point is x=x_i+r_i*s_i/sum_w. 
     */
    Cure::LinkedArray Vectors;

  public:
    MapPointFeature();

    /**
     *  @param wp The templete feature, CastPtr, WeightThreshold, 
     *  TriangleThreshold and DistanceThreshold are copied from the templete.
     */
    MapPointFeature( MapPointFeature *wp,MapBank *b=0);

    virtual void narrow(){*CastPtr=this;}

    void setCenter(MapPoint *val){addPoint(val,0);}

    /**
     * @param c a unit 3-vector pointing at point c, |c|=1,in the map frame
     * @return 0 if this can be reliably matched to based on testValue
     */
    int testVisable(double  c[3]){
      /** @todo explain this! */
      if (Bdual.Columns>0)return 0;

      double testValue=1-(LastBearing[0]*c[0]+
                          LastBearing[1]*c[1]+
                          LastBearing[2]*c[2]);
      if (testValue>TrackThreshold)
	return NOT_VISABLE;
      return 0;
    }

    int inside(double bottomleft[2], double topright[2]){
      return FeatureFcn::pointInRectangle(Points[0]->X,bottomleft, topright);
    }

    unsigned short getMeasurementType(unsigned short type);

    virtual int extend();
    /*
     * Use this to unconditionlly extend to full dimension.
     * The Index[0]=-1 is the flag that this has been done.
     */
    void forceExtend();
    /**
     * Adds the new collected Info to the cumulation and updates the 
     * LowInfo dimensions.
     * @param v a row matrix(1,8)  with this structure:
     * distance, w, x_i, s_i  (0..7)
     * w is a relative weight estimated to approximate the amount of
     * infomation contained in this observation (inverse variance)
     * x_i is the xyz of the camera when the observation was made.
     * s_i is a vector from x_i towards the point (magnitude ignored)
     * @param type 0.  
     * @return 0
     */
    virtual int  addInfo(const Cure::Matrix & v,const int type=0);
    virtual int  addInfo(const Cure::Matrix & v,
			 Cure::Transformation3D & map2info,
			 const int type=0);
   /**
    * clears all the initialization info and adds this instead.
    * @param v a  matrix(r,8)  with this structure for each row:
    * distance, w, x_i, u_i  (0..7)
    * w is a relative weight estimated to approximate the amount of
    * infomation contained in this observation (inverse variance)
    * x_i is the xyz of the camera when the observation was made.
    * u_i is an unit vector from x_i towards the point (magnitude must be one)
    * in the 'info frame' (a global map frame).
    * @return 0 if not enought for extending to 3D point or 
    *         1 if it can now be initialized with extend().
    * @param minangle the min of the max square of the angle in rads between
    * two bearings before initialization is done (extra threshold). 
    */
    int  addFullInfo(const Cure::Matrix & v, double minAngle=0);
    bool testMatch(MapPointFeature *mw, double tolerance);
    bool getC(MapPointFeature *mw,Cure::Matrix &c)const;
   
    virtual void getB(Cure::Matrix &b)const;

    /**
     * @param distanceguess  relative to camera.
     */
    void initializeFromPixels(Cure::Transformation3D &t,
			      const double centerpixels[2],
			      double focallength,
			      double distanceguess);

    int merge(MapPointFeature *mf, unsigned short type=1);
  
    /**
     * Write the MapPoint information to a text file.
     * @param fs the file to write to. It should be opened with std::ios::out
     */
    void write(std::fstream &fs ){
      if (Bdual.Columns<3)return;
      MapFeature::write(fs);
    }
      
    /**
     * This initializes a MapFeture from the information in a text file.
     * @param version 1 only for now.
     * @param fs the file to read from. It should be opened with std::ios::in
     * @param readKey set this to false if key has been read allready.
     * @return 1 if fail else 0
    */
    virtual int read(int version, std::fstream &fs, bool readKey=true){
      if (version!=1)return 1;
      int ret=MapFeature::read(1,fs,readKey);
      if (ret)return ret;
      recenter();
      return 0;
    }
    /**
     * Restores the feature from a GenericData object.  If the
     * MapPoints of gd already exists on the Bank then the data is
     * simply copied to the existing object.  By already exists we
     * mean that a MapPoint is associated with gd's BankID, Key
     * combination.  In other cases the Points on this feature are
     * removed and the new ones added.  If the removed points have no
     * features they are deleted. This doesn't bother with any
     * initialization info since the feature is likely to be initialized.
     * That includes all the thresholds and vectors.
     *
     * @param gd the Data to restore this feature from.
     *
     * @return 0 if ok,1 if ds is formated wrong, MAP_OBJECT_INVALID
     * if the point's keys are associated with the wrong type MapObject.
     * 
     */
    virtual int set(GenericData &gd);

  protected:
    void trackPoint(Cure::Matrix &m);

    /**
     * Removes the old information from Vectors. 
     */
    void prune(double minDistance);

    int eigen_it(Cure::Matrix & cloud, double wsums[4],
		 double ev[9], double lambda[3]);

    double reestimate_it(double oldx[4],double newx[4],double ev[9],
			 double lambda[3]);

    double estimate_it(double wsums[4],double ev[9],double lambda[3]);
  };

} // namespace Cure

#endif
