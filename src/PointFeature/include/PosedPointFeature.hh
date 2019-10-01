// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_POSEDPOINTFEATURE_HH
#define CURE_POSEDPOINTFEATURE_HH

#include "PosedFeature.hh"
#include "MapPointFeature.hh"
#include "Point3D.hh"
#include "Match.hh"
namespace Cure {

  /**
   * "Posed" feature of type point
   *
   * @author John Folkesson
   * @see
   */
  class PosedPointFeature: public PosedFeature
  {
  public:
    /** @todo comment how this works */
    PosedPointFeature **CastPtr;

    double FocalLength;
    Cure::Point3D Center;
    Cure::Point2D PixelCenter;
    Cure::Point3D *AllocatedPoints3D[1];
  public:
    PosedPointFeature(PosedPointFeature **pcc,MapBank *b=0, 
		      long fKey=-1, long id=-1);
    PosedPointFeature(PosedPointFeature **pcc, double focal,
		      MapPointFeature *mcc,long id=-1);
    PosedPointFeature(PosedPointFeature & pw);
    void operator =(PosedPointFeature & pw);
    PosedFeature * copy(){
      return new PosedPointFeature(*this);
    }

    virtual void narrow(){*CastPtr=this;}
    virtual void init();

    int transform(Cure::Transformation3D& pose,unsigned short covType);
    int calcZ();

    /**
     * Does a step match .
     * @param measurement =().
     * @param thresholds =()
     * @return 0 for a match  
     */
    virtual int roughMatch(double *z,double *thresholds);

    int inRectangle(double rect[4]){
      // On the right side of the camera, i.e. the predicted point
      // pose has y > f
      if (Center(1)<(FocalLength+1E-15))return 0;

      return (FeatureFcn::pointInRectangle(PixelCenter.X, rect,rect+2));
    }

    virtual int  addInfo(Match &mat,
			 Cure::Transformation3D & pose,unsigned short covType,
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
 

    int testVisable(double  bottomLeft[2],double topRight[2]);
    /*
      int pixel(PixelFragment& pix,double  bottomLeft[2],
      double topRight[2],double pixelWidth,
      double pixelHeight);
      void  pixelize(double x[2], int pix[2],double pixelWidth,
      double pixelHeight)
      {
      if (x[0]<0) pix[0]=(int)((x[0]-1)/pixelWidth);
      else pix[0]=(int)(x[0]/pixelWidth);
      if (x[1]<0) pix[1]=(int)((x[1]-1)/pixelHeight);
      else pix[1]=(int)(x[1]/pixelHeight);
      }
    */
  };

} // namespace Cure

#endif
