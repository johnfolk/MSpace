// = AUTHOR(S)
//    John Folkesson
//    
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_POSEDHLINE_HH
#define CURE_POSEDHLINE_HH

#include "MapHLine.hh"
#include "PosedFeature.hh"
#include "Line3D.hh"
#include "Match.hh"

namespace Cure{

  /**
   * "Posed" feature of type horizontal line
   *
   * @author John Folkesson
   * @see
   */
  class PosedHLine: public PosedFeature
  {
  public:
    PosedHLine **CastPtr;
    double FocalLength;
    short Reflected;
    Cure::Vector3D Vertical;
    Cure::Line3D Line;
    Cure::Line2D PixelLine;
    Cure::Vector3D Tangent;
    Cure::Point3D *AllocatedPoints3D[2];
  public:
    PosedHLine(){}
    PosedHLine(PosedHLine **pt,MapBank *b=0, long fKey=-1, long id=-1);
    PosedHLine(PosedHLine **pt,double focalLength,
	       MapHLine * mw, long id=-1);
    PosedHLine(PosedHLine & pw);
    void operator =(PosedHLine & pw);
    PosedFeature * copy(){
      return new PosedHLine(*this);
    }
    virtual void  narrow(){*CastPtr=this;}
    virtual void init();

    int transform(Cure::Transformation3D& pose,unsigned short covType);
    int calcZ();
    virtual int roughMatch(double  *z,double *thresholds);  
    /**
     * Check if image in rectangle of image plane
     */
    int inRectangle(double rect[4]){
      return FeatureFcn::lineInRectangle(PixelLine.StartPoint.X,
                                         PixelLine.EndPoint.X, rect,rect+2);
    }
    virtual double fineMatch(double *  z,Cure::Matrix & sigma);
     virtual int addInfo(Match &mat,Cure::Transformation3D & pose, 
			 unsigned short covType,
			Cure::Transformation3D & map2info, 
			const int type=0);

    /**
     * Calculates Eta, Jeo and Jev.  
     *
     * @param v The measurment vector.
     * @param v column vector(5,1) (startPixelX, startPixelY, endPixelX, 
     *                               endPixelY, FocalLength).
     * @param mtype Specifies definiation of Eta
     * @param noJac set to 1 to stop the calc after Eta.
     * @return 1 if fails 0 on success. 
     */
    virtual int calcEta(Cure::Matrix & v,unsigned short mtype, 
			unsigned short noJac=0);

  };

} // namespace Cure

#endif
