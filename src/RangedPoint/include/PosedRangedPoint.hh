// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_POSEDRANGEDPOINT_HH
#define CURE_POSEDRANGEDPOINT_HH

#include "PosedFeature.hh"
#include "MapRangedPoint.hh"
#include "Point3D.hh"
#include "Match.hh"
namespace Cure {

  /**
   * "Posed" feature of type point
   *
   * @author John Folkesson
   * @see
   */
  class PosedRangedPoint: public PosedFeature
  {
  public:
    /** @todo comment how this works */
    PosedRangedPoint **CastPtr;

    Cure::Point3D Center;
    Cure::Point3D *AllocatedPoints3D[1];
  public:
    PosedRangedPoint(PosedRangedPoint **pcc,MapBank *b=0, 
		      long fKey=-1, long id=-1);
    PosedRangedPoint(PosedRangedPoint **pcc,
		      MapRangedPoint *mcc,long id=-1);
    PosedRangedPoint(PosedRangedPoint & pw);
    void operator =(PosedRangedPoint & pw);
    PosedFeature * copy(){
      return new PosedRangedPoint(*this);
    }

    virtual void narrow(){*CastPtr=this;}
    virtual void init();

    int transform(Cure::Transformation3D& pose,unsigned short covType);
    /**
     * This changes the bearing to the point to the predicted measurement Z.
     * z=phi,theta, r where phi is angle between bearing and sensor y axis
     * theta is angle between x axis and 
     * projection of bearing on x-z plane 
     * (between +x and +z  quad is positive theta.
     */
    int calcZ();

    /**
     * Does a step match .
     * @param measurement =().
     * @param thresholds =()
     * @return 0 for a match  
     */
    virtual int roughMatch(double *z,double *thresholds);

    int inRectangle(double rect[4]){
      return (FeatureFcn::pointInRectangle(Center.X, rect,rect+2));
    }

    virtual int  addInfo(Match &mat,
			 Cure::Transformation3D & pose,unsigned short covType,
			 Cure::Transformation3D & map2info,
			 const int type=0);

    /**
     * Calculates Eta, Jeo and Jev.  
     *
     * @param v The measurment vector.
     * @param mtype Specifies definiation of Eta
     * @param noJac set to 1 to stop the calc after Eta.
     * @return 1 if fails 0 on success. 
     */
    virtual int calcEta(Cure::Matrix & v,unsigned short mtype, 
			unsigned short noJac=0);
 

    int testVisable(double  bottomLeft[2],double topRight[2]);
  };

} // namespace Cure

#endif
