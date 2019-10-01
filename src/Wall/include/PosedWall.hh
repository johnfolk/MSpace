// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_POSEDWALL_H
#define CURE_POSEDWALL_H

#include "PosedFeature.hh"
#include "MapWall.hh"
#include "Trig.hh"
#include "Line2D.hh"
#include "Match.hh"

namespace Cure {


  /**
   * "Posed" feature of type wall
   *
   * @author John Folkesson
   */
  class PosedWall: public PosedFeature
  {
  public:
    PosedWall **CastPtr;
    Cure::Trig *Triger;
    Cure::Point2D *AllocatedPoints2D[2];
    Cure::Line2D Line;
    double EndPointUpdateThreshold;
  public:
    PosedWall(){}
    PosedWall(PosedWall **wp,Cure::Trig *t,MapBank *b=0, long fKey=-1, long id=-1);
    PosedWall(PosedWall **wp,Cure::Trig *t,MapWall *w,long id=-1);
    PosedWall(PosedWall & pw);
    void operator =(PosedWall & pw);
    PosedFeature * copy(){
      return new PosedWall(*this);
    }
    void narrow(){*CastPtr=this;}
    virtual void init();
    /**
     * This must be called before any other methods each time the
     * sensor frame changes. 
     *
     * @param pose, the map to sensor transformation, Tms.
     * @param covType this sets the columns of Jes. 
     * @return NOT_VISABLE; if not in image plane or behind wall ect. 
     * of  MAP_OBJECT_INVALID if FeatureKey is not good.
     */
    virtual int transform(Cure::Transformation3D& pose, 
			  unsigned short covType);
    int calcZ();
    /**
     * Does a step match to gamma rho and checks for overlap in scan
     * angles.  If Thresholds[2]==0 it skips the overlap part.  That
     * is usful if you want to check if it is say within .01 rads and
     * 5 cm so you could skip/overide mahalanobis check.
     *
     * @param measurement =(Gamma, Rho, Start rads, End rads).
     * @param thresholds =(Max Gamma Error, Max Rho Error, Min overlap in rads)
     * @return 0 for a match  
     */
    virtual int roughMatch(double * z,double *thresholds);
    int inRectangle(double rect[4]){
      return FeatureFcn::lineInRectangle(Line.StartPoint.X,
                                         Line.EndPoint.X, rect,rect+2);
    }
    virtual double fineMatch(Cure::Measurement &m,Cure::Matrix & sigma);
    virtual double fineMatch(double * z,Cure::Matrix & sigma);    
    virtual int consistentPose(Cure::Measurement &m, 
			       Cure::Transformation3D &pose);

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
			 const int type=0);


    /**
     * 
     * Calculates Eta, Jeo and Jev.  
     * @param v (sx,sy, ex, ey) relative to sensor.
     * @param mtype 0 for no dim 1 for 2 dim 3 for start + 2D, 
     *  5 for end + 2D, 7 for 4D  
     * @param noJac set to 1 to stop the calc after Eta.
     * @return 1 if fails 0 on success. 
     */
    virtual int calcEta(Cure::Matrix & v,unsigned short mtype,
			unsigned short noJAc=0);
    int scan(ScanFragment& s);

  };

} // namespace Cure

#endif
