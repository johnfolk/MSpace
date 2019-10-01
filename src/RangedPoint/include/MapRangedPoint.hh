// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_MAPRANGEDPOINT_HH
#define CURE_MAPRANGEDPOINT_HH


#include "MapFeature.hh"
#include "FeatureFun.hh"
#include "MatrixStuff.hh"
#include "LinkedArray.hh"

namespace Cure {


  /**
   * Class managing point features in the map
   * 
   * @author John Folkesson
   */
  class MapRangedPoint: public MapFeature
  {
  public:
    MapRangedPoint **CastPtr;
    long AllocatedIndex[3];
    MapPoint * AllocatedPoints[1];

  public:
    MapRangedPoint();

    /**
     * 
     *  @param wp The templete feature, CastPtr, 
     *  and DistanceThreshold are copied from the templete.
     */
    MapRangedPoint( MapRangedPoint *wp,MapBank *b=0);
    /**
     * This will make the MapPoint and set the center to x relative to
     * sensorPose.
     * @param sensorPose the transformation to the sensor frame that x is in
     * @param x this point, 'p', will be set to sensorPose.transform(x,p)
     */
    void initialize(Cure::Transformation3D &sensorPose, double x[3]){
      if (!Points[0])
	{
	  MapPoint *p=new MapPoint(Bank);
	  setCenter(p);
	}
      Points[0]->initializeFromRelative(sensorPose,x);
    }
    /**
     * Normally the CastPtr is a pointer on the Helper.  By calling this 
     * the helper can cast this MapFeature to a MapRangedPoint.
     */
    virtual void narrow(){*CastPtr=this;}
    /** 
     * This attaches the MapPoint to the feature.
     */
    void setCenter(MapPoint *val){addPoint(val,0);}

    /**
     * Reamains from vision features.  This might be used eventually here too.
     * @param c a unit 3-vector pointing at point c, |c|=1,in the map frame
     * @return 0 if this can be reliably matched to based on testValue
     */
    int testVisable(double  c[3]){
      return 0;
    }
    /**
     *  Test if point is inside a 2D rectangle in xy.
     */
    int inside(double bottomleft[2], double topright[2]){
      return FeatureFcn::pointInRectangle(Points[0]->X,bottomleft, topright);
    }
    /**
     * This returns 0 unless type is 1 and pDim=3.  Then it returns 1.
    */
    unsigned short getMeasurementType(unsigned short type);
    /**
     * This simple changes the pDim to 3 if it isn't already,
     * @return 3 if changed else 0.
     */
    virtual int extend();
    /*
     * Use this to unconditionlly extend to full dimension.
     * The Index[0]=-1 is the flag that this has been done.
     */
    void forceExtend();
    /**
     * These do nothing right now as we don't do initialization on this yet.
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
    int  addInfo(const Cure::Matrix & v,const int type=0);
    int  addInfo(const Cure::Matrix & v,
			 Cure::Transformation3D & map2info,
			 const int type=0);
    /**
     * @return true if both this and mw are pDim=3.
     */
    bool testMatch(MapRangedPoint *mw, double tolerance);

    /**
     * c is returned as 3X3 I for this.
     * @return true if both this and mw are pDim=3.
     */
    bool getC(MapRangedPoint *mw,Cure::Matrix &c)const;   
    void getB(Cure::Matrix &b)const;

    int merge(MapRangedPoint *mf, unsigned short type=1);
  
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


  protected:
  };

} // namespace Cure

#endif
