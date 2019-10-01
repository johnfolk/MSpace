// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_HLINEHELPER_HH
#define CURE_HLINEHELPER_HH

#include "MapHLine.hh"
#include "PosedHLine.hh"
#include "MapHelper.hh"
#include "Match.hh"

namespace Cure{

/**
 * Helper class for horizontal lines
 *
 * @author John Folkeeson
 * @see MapHelper
 */
  class HLineHelper: public MapHelper 
  {
  public:
    PosedHLine * PosedHLinePtr;
    MapHLine * MapHLinePtr;
    MapHLine TempletFeature;
    Cure::Line2D Line;
    double DistanceGuess;
    /** 
     * The Templet will be given these thresholds which then are used to 
     * initilize the MapFeatures made with this helper. 
     */
    HLineHelper(MapBank *b, double tangentThreshold=1,double rhoThreshold=1.5,
		double distanceThreshold=20, double trackThreshold=.5)
      :MapHelper(b)
    {
      TempletFeature.setBank(b);
      TempletFeature.CastPtr=&MapHLinePtr;
      TempletFeature.TanThreshold=tangentThreshold;
      TempletFeature.RhoThreshold=rhoThreshold;
      TempletFeature.DistanceThreshold=distanceThreshold;
      TempletFeature.TrackThreshold=trackThreshold;
      SensorType=CAMERA_TYPE;
      DistanceGuess=2;
      FeatureType=MAPHLINE_TYPE;
      ObjectSubType=FeatureType;
    }
    ~HLineHelper(){}
    virtual void printConfiguration();
    PosedHLine * castPosedHLine(PosedFeature *pf){
      PosedHLinePtr=0;
      pf->narrow();
      return PosedHLinePtr;
    }
    MapHLine * castMapHLine(MapObject *pf){
      MapHLinePtr=0;
      pf->narrow();
      return MapHLinePtr;
    }

    MapHLine *getMapHLine(long key)
    {
      MapObject *mo=Bank->getMapObject(key);
      if (!mo) return 0;
      return castMapHLine(mo);
    } 

    virtual MapFeature * makeMapFeature(MapBank *b=0){
      return new MapHLine(&TempletFeature,b);
    }

    /**
     * @return A MapHLine based on the TempletFeature. 
     */
    MapHLine * makeMapHLine(){
      return new MapHLine(&TempletFeature,Bank);
    }

    virtual MapFeature *getMapFeature(long key){
      return getMapHLine(key);
    }
    /**
     * 
     * @param camerapose The pose of the camera with the y axis along the 
     *        camera axis, x is right and z up.
     * @param imagex The coordinates of the start and endpoints in the image
     * in meter units.
     * @param heightguess estimated hieght of the line above the camera pose.
     * @return A MapHLine based on the TempletFeature and the params. 
     */
    MapHLine * makeMapHLine(Cure::Pose3D camerapose,double imagex[4], 
			    double heightguess=2.14) 
    {
      MapHLine *templine=new MapHLine(&TempletFeature,Bank);
      templine->initializeFromPixels(camerapose,imagex,imagex+2,
				     PixelInfo[0],heightguess);
      return templine;
    }     

  


    virtual PosedFeature * makePosedFeature(MapFeature* mf){
      return makePosedHLine(PixelInfo[0],mf);}
    virtual MapFeature * makeMapFeature(Cure::Pose3D & sensorpose,
					Cure::Measurement &m,
					bool addToVisable=false);
    /**
     * 
     * @param camerapose The pose of the camera with the y axis along the 
     *        camera axis, x is right and z up.
     * @param heightguess estimated hieght of the line above the camera pose.
     * @return A MapHLine based on the TempletFeature and the params. 
     */  
    MapHLine * makeMapHLine(Cure::Pose3D & camerapose, 
			    double heightguess=2.0) 
    {
      MapHLine *templine=new MapHLine(&TempletFeature,Bank);
      templine->initializeFromPixels(camerapose,Line.StartPoint.X,
				     Line.EndPoint.X,
				     PixelInfo[0],heightguess);
      return templine;
    }   
  
    PosedHLine * makePosedHLine(double focal,MapFeature* mf){
      if(mf->Type==MAPHLINE_TYPE)
	return new PosedHLine(&PosedHLinePtr,focal,castMapHLine(mf));
      return 0;
    }
    PosedHLine * makePosedHLine(MapFeature* mf){
      return makePosedHLine(PixelInfo[0],mf);
    }

    virtual bool getC(MapFeature  *mf1,MapFeature  *mf2,
		      Cure::Matrix &c1){
      MapHLine *mw1=castMapHLine(mf1);
      if (!mw1)return false;
      MapHLine *mw2=castMapHLine(mf2);
      if (!mw2)return false;
      return mw1->getC(mw2,c1);
    }
    virtual bool testMatch(MapFeature  *mf1,MapFeature  *mf2,
			   double tolerance){
      MapHLine *mw1=castMapHLine(mf1);
      if (!mw1)return false;
      MapHLine *mw2=castMapHLine(mf2);
      if (!mw2)return false;
      return mw1->testMatch(mw2,tolerance);
    }
    virtual void match(Match *matches,int n,
		       double distance,Cure::Pose3D &pcam, LongList &mapFeats);
  
    virtual int config(const std::string &arglist);

    virtual bool supportsSubconfig(int sc);

    int addDistance(Match &mat,double distance);

    /**
     * This will set the SqWidth of the feature based on the two 
     * @param m1 one of the matched Measurments of a MapHline
     * @param m2 The other of the matched Measurments of the MapHline
     */
    void setWidth(Match & m1,Match & m2);
    /**
     * This will set the CovV of the Measurment based on the Width of the 
     * MapHLine and the current CovV. 
     * @param m1 The matched Measurments of a MapHline
     */
    void setCovV(Match & m1);

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
    virtual unsigned short findMergeFeatures(PosedFeature * pfp[2],
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
      MapHLine *mf0=getMapHLine(pf[0]->FeatureKey);
      MapHLine *mf1=getMapHLine(pf[1]->FeatureKey);
      if (!mf1)return -1;
      if (!mf0)return -1;
      int r=mf1->merge(mf0,typ);  
      if (typ&1){
	delete mf0;
      }
      return r;
    }

  protected:
    int configVer1(const std::string &arglist);
  };

} // namespace Cure

#endif
