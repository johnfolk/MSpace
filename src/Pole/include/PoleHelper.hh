// = AUTHOR(S)
//    John Folkesson
//    
//    August 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef POLEHELPER_HH
#define POLEHELPER_HH

#include "MapPole.hh"
#include "PosedPole.hh"
#include "MapHelper.hh"
#include "Measurement.hh"
#include "Point2D.hh"
#include "Trig.hh"
namespace Cure{
  /**
   * THIS HAS NEVER BEEN TESTED!!!
   */
class PoleHelper: public MapHelper 
{
public:
  PosedPole * PosedPolePtr;
  MapPole * MapPolePtr;
  MapPole TempletFeature;
  Cure::Point2D Center;
  Cure::Trig Triger;
   /** 
   * The Templet will be given these thresholds which then are used to 
   * initilize the MapFeatures made with this helper. 
   */
  PoleHelper(MapBank *b, int countThreshold=25, double radiusThreshold=.5,
	     double distanceThreshold=20)
    :MapHelper(b){
    FeatureType=MAPPOLE_TYPE;
    TempletFeature.setBank(b);
    TempletFeature.CastPtr=&MapPolePtr;
    TempletFeature.CountThreshold=countThreshold;
    TempletFeature.RadiusThreshold=radiusThreshold;
    TempletFeature.DistanceThreshold=distanceThreshold; 
    ObjectSubType=FeatureType;
  }
  ~PoleHelper(){}

  PosedPole * castPosedPole(PosedFeature *pf){
    PosedPolePtr=0;
    pf->narrow();
    return PosedPolePtr;
  }
  MapPole * castMapPole(MapObject *pf){
    MapPolePtr=0;
    pf->narrow();
    return MapPolePtr;
  }
  MapPole *getMapPole(long key)
  {
    MapObject *mo=Bank->getMapObject(key);
    if (!mo) return 0;
    return castMapPole(mo);
  } 
  /**
   * @return A MapHLine based on the TempletFeature. 
   */
  MapPole * makeMapPole(){
    return new MapPole(&TempletFeature,Bank);
  }
  /**
   * 
   * @param sensorpose The pose of the sensro with the y axis ahead, 
   *        x is right and z up.
   * @param x The coordinates of the center 
   * in meter units relative to sensor.
   * @return A MapPole based on the TempletFeature and the params. 
  */
  MapPole * makeMapPole(Cure::Pose3D sensorpose,double x[2], double r=0) 
  {
    Center.setXY(x);
    MapPole *temppole=new MapPole(&TempletFeature,Bank);
    temppole->initailizeFromPoint(sensorpose,Center,r);
    return temppole;
}     
  /**
   * 
   * @param Sensorpose The pose of the sensor with the y axis along the 
   *        camera axis, x is right and z up.
    * @return A MapPole based on the TempletFeature and the params. 
  */
  MapPole * makeMapPole(Cure::Pose3D sensorpose) 
  {
    MapPole *temppole=new MapPole(&TempletFeature,Bank);
    temppole->initailizeFromPoint(sensorpose,Center);
    return temppole;
  }     
  PosedPole * makePosedPole(MapFeature* mf){
    if(mf->Type==MAPPOLE_TYPE)
      return new PosedPole(&PosedPolePtr,&Triger, castMapPole(mf));
    return 0;
  }
  /**
   * This takes pixel MapFeatures from the list 'Near' and makes PosedFeatures
   * out of them and adds to 'PixelFeature' list.  
   */
  PosedFeature * makePosedFeature(MapFeature* mf){
    return makePosedPole(mf);
  }
  /**
   * Search  list for closest PoleFeature to a point.
   * @param pp The matched PosedPole is copied to pp.
   * @param pt the point in the sensor frame that is 
   * being matched to.
   * @param sigma the covaraince in gamma/rho
   * @param poles the list to look for matches on, default will check 
   * all Visable. 
   * @return -1 if a match is not found else an estimated distance.
   */
  double matchPole(PosedPole & pp, Cure::Point2D & pt,Cure::Matrix & sigma,
		 PosedFeatureList * poles=0);

 
  /**
   * Search  list for closest PoleFeature to a point.
   * @param pp The matched PosedPole is copied to pp.
   * being matched to.
   * @param sigma the covaraince in gamma/rho
   * @param poles the list to look for matches on, default will check 
   * all Visable. 
   * @return -1 if a match is not found else an estimated distance.
   */
  double matchPole(PosedPole & pp,Cure::Matrix & sigma,
		    PosedFeatureList * poles=0)
  {
    return matchPole(pp,Center,sigma,poles);
  }
 
  
  /**
   * Search poles list for closest PoleFeature to a point.
   * @param matches The matched PosedHLine's are added to matches.
   * @param pr the point in the sensor frame that is being matched to.
   * @param searchlimit the farthest distance away to look for matches,
   * this sets the limits on bounding box, rho and gamma.
   * @param lines the list to look for matches on, default will check 
   * all PixelVisable. 
   */
  void matchPoles(PosedFeatureList & matches, Cure::Point2D & pt, 
		  double searchlimit,
		  PosedFeatureList * poles=0);
  
 
  void matchPoles(PosedFeatureList & matches, 
		   double searchlimit,
		   PosedFeatureList * poles=0){
    matchPoles(matches, Center, searchlimit,poles);
  }

  /**
   * @param x The coordinates of the center and radius in the sensor frame
   * @param distance The distance to add info with.
   * @param m the Mesurement is returned in this.
   * @return
   */
  double makeMeasurement(double x[3],double distance, Cure::Measurement & m);
};
}
#endif
