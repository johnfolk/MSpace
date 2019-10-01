// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef MAPPOLE_H
#define MAPPOLE_H


#include "MapFeature.hh"
#include "FeatureFun.hh"
#include "MatrixStuff.hh"
#include "CircleCloud.hh"

namespace Cure{
  /**
   * THIS HAS NEVER BEEN TESTED!!!
   */
class MapPole: public MapFeature
{
public:
  MapPole **CastPtr;
  double Radius;
  double RadiusThreshold;
  double CountThreshold;
  int InfoInitialized;
  MapPoint * AllocatedPoints[1];
  long AllocatedIndex[3];
  Cure::CircleCloud Cloud; //distance, x,y
public:
  MapPole();
  MapPole(MapPole * wp,MapBank *b=0);
  virtual void narrow(){*CastPtr=this;}
  void setCenter(MapPoint *val);
  void setRadius(double r);
  void initailizeFromPoint(Cure::Transformation3D &t,Cure::Point2D & pt, double r=0);
  int inside(double bottomleft[2], double topright[2]){
    return FeatureFcn::pointInRectangle(Points[0]->X,bottomleft, topright);
  }
  virtual int extend();

  /**
   * Adds the new collected Info to the cumulation and updates the 
   * LowInfo dimensions.
   * Have not implemented type!=0 yet
   * @param v A cloud of (distance,x,y) (N x 3) 
   * used to set center and radius.
   * @param covAdjust for EKF do cov -> covAdjust*cov*covAdjust.
   * @param type 0.  
   */
  int  updateInfo(const Cure::Matrix & v, 
		 Cure::Matrix & covAdjust,const int type=0);
    virtual void getB(Cure::Matrix &b)const;
};
}
#endif
