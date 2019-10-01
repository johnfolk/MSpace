// = AUTHOR(S)
//    John Folkesson
//    
//    August 1, 2007
//
//    Copyright (c) 2007 John Folkesson
//    
#ifndef CURE_POSEDRELATIVEPOINTFEATURE_HH
#define CURE_POSEDRELATIVEPOINTFEATURE_HH

#include "PosedFeature.hh"
#include "RelativePointFeature.hh"
#include "Match.hh"
namespace Cure {

  /**
   * This is not exactly like the other PosedFeatures since 
   * it depends on the reference pose as well as the feature 
   * coordinates and sensorpose.  
   *
   * Transform causes Xo and Z to be calced but Xo is the point in the
   * sensor frame while the transformed feature coordinates are
   * scalars.  Z is (unit bearing vector, phi,theta,r,rho).
   *
   * calcEta causes an innovation to be calcualted but no jacobians.
   * The innovation is quite simple to estimate uncertianty for based
   * on for instance pixel resolution or sensor specs.  So for bearing
   * only the Eta measurement uncertianty should be the angular
   * uncertianty.  for range /bearing it is the phi/theta/range
   * uncertinties (polar coordingates), with the possible scaleing for
   * the cubic polynomial in theta.  
   *
   * Eta can then be used to test matches.
   *
   * the further calc is not implemented
   *
   * @author John Folkesson
   * @see
   */
  class PosedRelativePointFeature: public PosedFeature
  {
  public:
    PosedRelativePointFeature **CastPtr;
    Transformation3D SensorPose;
  public:
    PosedRelativePointFeature(PosedRelativePointFeature **pcc,MapBank *b=0, 
			      long fKey=-1, long id=-1);
    PosedRelativePointFeature(PosedRelativePointFeature **pcc,
			      RelativePointFeature *mcc,long id=-1);
    PosedRelativePointFeature(PosedRelativePointFeature & pw);
    void operator =(PosedRelativePointFeature & pw);
    PosedFeature * copy(){
      return new PosedRelativePointFeature(*this);
    }
    RelativePointFeature * getRelativePointFeature(){
      MapFeature *mf =getFeature();
      if (mf)
	if (mf->Type==RELATIVEPOINTFEATURE_TYPE)
	  return (RelativePointFeature*)mf;
      return 0;
    }

    void narrow(){*CastPtr=this;}
    void init();

    /**
     * This will calcuate the Z and Xo needed to do inRectangle and
     * rough match.
     * @param pose the sensor frame pose
     * @param covType the covariance type of the sensor pose.
     */
    int transform(Cure::Transformation3D& pose,unsigned short covType);
    /**
     * Make a first easy calc to cut out the very far measurements.
     * This compares sin^2(angle between predicted and measured
     * bearing)< theshold[0] for bearing only measurments.  For range
     * bearing it compares (phi,theta,range) to predicted and put
     * thresholds[0..2] to the absolute value of the differences.
     *
     * @param z So for Mesurment type 32 z should be a unit bearing vector in
     * the sensor frame and thresholds[0] a thrshold on sin^2(angle).
     * for others z should be (phi,theta,range and thresholds three
     * limits on these errors. 
     *
     * @return 0 for a match else 1;
     */
    int roughMatch(double  *z,double *thresholds);
    /**
     * A finer test of the features from the map.  This feature is
     * compared to the MeasurmentType specific V using the metric
     * sigma and the FineThreshold.  Eta is calcutated for the
     * relative feature type.  
     *
     *@param m The measurement to be matched should be given here, the
     *predicted measurement is subtracted from this.
     *
     * @param sigma the metric for the error in Eta
     * @return -1 if fails else the distance estimated using sigma.
     */
    double fineMatch(Cure::Measurement &m,Cure::Matrix & sigma){
      if (!(CalcState&1))return -1;
      calcEta(m.V,m.MeasurementType,1);
      //Here we need to have Eta not be set by the pdim
      //but by the type
      Matrix t,temp;
      int r=Eta.Rows;
      int sr=sigma.Rows;
      if (r<sr)sigma.offset(0,r);
      else if (r>sr)Eta.offset(0,0,sr,1);
      temp.multiply_(sigma,Eta);
      t.multTranspose_(Eta,temp,1);
      sigma.offset(0,sr);
      Eta.offset(0,0,r,1);
      return t(0,0);
    }
    /**
     * This is not the same method as fineMatch on a Measurement. For
     * MeasurementType 34 and 35 it is a one dimensional check on 
     * the sin^2 of the angle between the predicted and measured bearings.
     * For types 1,3,5,7 it is a 3D  phi, f(theta), range comparision.
     * This is not called by the helper.
     */
    double fineMatch(double * z,Cure::Matrix & sigma);
    /**
     * This calculates the Eta but nothing more.
     * One must call transform before this.
     * No jacobians are computed.
     *
     * For Bearing only.
     *
     * Eta(0,0) is defined as the angle between the measurement
     * bearing and a bearing to the point of closest approach between
     * the features bearing in the reference frame and the measurement
     * bearing.  So the angle to the epipolar plane.  
     *
     * Eta(1,0) is defined as the angle between the feature point and
     * the point of closest approach in the measurement sensor frame.
     * So the two angles are 'perpendicular' for small angles and
     * Eta(1,0) is not defined if the motion was not enough.
     *
     * For range bearing the Eta is the error in (phi,theta,range)
     * or in the case of MeasurmentType 3 and 7 (phi,f(theta),range)
     * where f(thete) is  a cubic polynomial.
     * @v For bearing only the unit bearing vector in the sensor frame.
     * for rangebearing it is (phi,f(theta),range).
     * mtype the type of measurement 32,33 = bearing only 1,3,5,7 range/bearing
    */
    int calcEta(Cure::Matrix & v,
		unsigned short mtype, 
		unsigned short noJac);
     /**
     * Use this to check if the feature is in some rectangle in the
     * sensor frame.  There are 4 cases depending on type of
     * RelativePointFeature:
     *
     * bearingOnly:
     * rect is (low phi, low theta, high phi, high theta)
     * note  that phi is  in (0,pi) and theta is (-pi,pi)
     * so if theta range is around pi then you are in trouble. 
     *
     * phiRange
     * rect is  (low phi, low range, high phi, high range)
     *
     * thetaRange
     * rect is  (low theta, low range, high theta, high range)
     *
     * else
     *
     * rect is (low x, low y, high x, high y)
     *
     * This is in polar coordinates not like PointFeature
     * @param rect The corners of the rectangle (minphi,mintheta, maxphi,maxytheta) 
     *  @return 1 if inside
     */
    int inRectangle(double rect[4]){
      if (!(CalcState&1))return 0;
      RelativePointFeature *mf=getRelativePointFeature();
      if (!mf)return 0;
      if ((mf->bearingOnly())){
	if (rect[0]>Z(3,0))return 0;
	if (rect[1]>Z(4,0))return 0;
	if (rect[2]<Z(3,0))return 0;
	if (rect[3]<Z(4,0))return 0;
	return 1;
      }
      if (mf->phiRange()){
	if (rect[0]>Z(3,0))return 0;
	if (rect[1]>Z(5,0))return 0;
	if (rect[2]<Z(3,0))return 0;
	if (rect[3]<Z(5,0))return 0;
	return 1;
      }
      if (mf->thetaRange()){
	if (rect[0]>Z(4,0))return 0;
	if (rect[1]>Z(5,0))return 0;
	if (rect[2]<Z(4,0))return 0;
	if (rect[3]<Z(5,0))return 0;
	return 1;
      }
      if (Xo(0,0)<rect[0])return 0;
      if (Xo(1,0)<rect[1])return 0;
      if (Xo(0,0)>rect[2])return 0;
      if (Xo(1,0)>rect[1])return 0;
      return 1;
    }
  };

} // namespace Cure

#endif
