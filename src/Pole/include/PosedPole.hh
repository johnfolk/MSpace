// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    

#ifndef POSEDTREE_H
#define POSEDTREE_H


#include "MapFeature.hh"
#include "FeatureFun.hh"
#include "PosedFeature.hh"
#include "MapPole.hh"
#include "Trig.hh"
#include "Point2D.hh"

namespace Cure{

/**
 * THIS HAS NEVER BEEN TESTED!!!
 * @author John Folkesson
 */
class PosedPole: public PosedFeature
{
public:
  PosedPole **CastPtr;
  double StartAngle;
  double EndAngle;
  Cure::Trig *Triger;
  double Normal[2];
  double A[4];
  int U[3];
  Cure::Point2D *AllocatedPoints2D[1];
  Cure::Point2D Center;
 public:
  PosedPole(PosedPole **pt,Cure::Trig *t,MapBank *b=0, long fKey=-1, long id=-1);
  PosedPole(PosedPole **pt,Cure::Trig *t,MapPole *tr,long id=-1);
  PosedPole(PosedPole & pt);
  void operator =(PosedPole & pt);
  PosedFeature * copy(){
    return new PosedPole(*this);
  }
  virtual void narrow(){*CastPtr=this;}
  virtual void init();

  /**
   * Does a step match to theta and rho and radius.
   * @param measurement =(theta, rho, radius).
   * @param thresholds =(Max theta Error, Max rho Error, Max radius Error)
   * @return 0 for a match  
   */
  virtual int roughMatch(double *z,double *thresholds);
 
  
  /**
   * using the given sigma and z 
   * @param z theta, rho
   * @param sigma the estimated covariance of the measurement error, (2x2).
   * @return -1 if fails else the mhalanobis distance estimated
   * from the infomation matrix of the MapFeature>=0.
   */
  double fineMatch(double *z,Cure::Matrix  & sigma);
  int transform(Cure::Transformation3D& pose,unsigned short covType);
  /**
   * Calculates Eta, Jeo and Jev.  
   * @param v Theta, rho, radius,
   * @param mtype 0 for just Theta Rho 1 for gamma,rho,radius.
   * @return 0 if ok -1 else.
   */
  virtual int calcEta(Cure::Matrix & v,unsigned short mtype);

  int scan(ScanFragment& s);
  int thetaOfU(double *r,int u);
  int thetaOfU(int u);

protected:
  void recursiv(int u0,int u1,int *a,double *r,int length);
};

} // namespace Cure
 
#endif
