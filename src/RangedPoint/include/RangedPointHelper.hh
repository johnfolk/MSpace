
// = AUTHOR(S)
//    John Folkesson
//    
//    June 11, 2006
//
//    Copyright (c) 2006 John Folkesson
//    
#ifndef CURE_RANGEDPOINTHELPER_HH
#define CURE_RANGEDPOINTHELPER_HH

#include "MapRangedPoint.hh"
#include "PosedRangedPoint.hh"
#include "MapHelper.hh"
#include "Match.hh"

namespace Cure {

  /** 
   * Helper class for RangedPoints
   *
   * The structure of the Measurements are assumed to be 
   * SensorType=SensorData::SENSORTYPE_CAMERA=2
   * MeasurementType=32
   * V: 3x1 (x,y,f)
   * Z: 2x1 (x,y)
   * W: 1x10 (imwidth, imheight, f, pixwidth, pixheight, 
   *          x_pix, y_pix, xsize_pix, ysize_pix, angle)
   */
  class RangedPointHelper: public MapHelper 
  {
  public:
    PosedRangedPoint * PosedRangedPointPtr;
    MapRangedPoint * MapRangedPointPtr;
    MapRangedPoint TempletFeature;
    /** 
     * The Templet will be given these thresholds which then are used to 
     * initilize the MapFeatures made with this helper. 
     */
    RangedPointHelper(MapBank *b, 
		       double distanceThreshold=20)
      :MapHelper(b){
      
      TempletFeature.setBank(b);
      TempletFeature.CastPtr=&MapRangedPointPtr;
      TempletFeature.DistanceThreshold=distanceThreshold; 
      FeatureType=MAPRANGEDPOINT_TYPE;
      ObjectSubType=FeatureType;
    }
    ~RangedPointHelper(){}
    virtual void printConfiguration();

    PosedRangedPoint * castPosedRangedPoint(PosedFeature *pf){
      PosedRangedPointPtr=0;
      pf->narrow();
      return PosedRangedPointPtr;
    }
    MapRangedPoint * castMapRangedPoint(MapObject *pf){
      MapRangedPointPtr=0;
      pf->narrow();
      return MapRangedPointPtr;
    }
    MapRangedPoint *getMapRangedPoint(long key)
    {
      MapObject *mo=Bank->getMapObject(key);
      if (!mo) return 0;
      return castMapRangedPoint(mo);
    } 
    virtual MapFeature *getMapFeature(long key){
      return getMapRangedPoint(key);
    }

    virtual MapFeature * makeMapFeature(MapBank *b=0){
      return new MapRangedPoint(&TempletFeature,b);
    }

    /**
     * @return A MapRangedPoint based on the TempletFeature. 
     */
    MapRangedPoint * makeMapRangedPoint(){
      return new MapRangedPoint(&TempletFeature,Bank);
    }
    /**
     * 
     * @param sensorpose The pose of the sensor with the y axis ahead, 
     *        x is right and z up.
     * @param v The xyz coordinates of the MapRangedPoint in relative
     * to the sensor.
     * @return A MapRangedPoint based on the TempletFeature and the params. 
     */
    MapRangedPoint * makeMapRangedPoint(Cure::Transformation3D &sensorpose, 
					double v[3]) 
    {
      MapRangedPoint *temppt=new MapRangedPoint(&TempletFeature,Bank);
      temppt->initialize(sensorpose,v);
      return temppt;
    }     
    virtual PosedFeature * makePosedFeature(MapFeature* mf){
      return makePosedRangedPoint(mf);}
    virtual MapFeature * makeMapFeature(Cure::Pose3D &sensorpose,
					Cure::Measurement &m,
					bool addToVisable=false);
    PosedRangedPoint * makePosedRangedPoint( MapFeature* mf){
      if(mf->Type==MAPRANGEDPOINT_TYPE)
	return new PosedRangedPoint(&PosedRangedPointPtr,
				     castMapRangedPoint(mf));
      return 0;
    }
    virtual bool getC(MapFeature  *mf1,MapFeature  *mf2,
		      Cure::Matrix &c1){
      MapRangedPoint *mw1=castMapRangedPoint(mf1);
      if (!mw1)return false;
      MapRangedPoint *mw2=castMapRangedPoint(mf2);
      if (!mw2)return false;
      return mw1->getC(mw2,c1);
    }
    virtual bool testMatch(MapFeature  *mf1,MapFeature  *mf2,
			   double tolerance){
      MapRangedPoint *mw1=castMapRangedPoint(mf1);
      if (!mw1)return false;
      MapRangedPoint *mw2=castMapRangedPoint(mf2);
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
      MapRangedPoint *mf0=getMapRangedPoint(pf[0]->FeatureKey);
      MapRangedPoint *mf1=getMapRangedPoint(pf[1]->FeatureKey);
      if (!mf1)return -1;
      if (!mf0)return -1;
      int r=mf1->merge(mf0,typ);  
      if (typ&1){
	delete mf0;
      }
      return r;
    }
  protected:

    /**
     * Reads configurations of version 1
     */
    int configVer1(const std::string &arglist);
  };

} // namespace Cure

#endif
