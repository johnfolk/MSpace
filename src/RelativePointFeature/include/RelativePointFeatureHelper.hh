// = AUTHOR(S)
//    John Folkesson
//    
//    August 1, 2007
//
//    Copyright (c) 2007 John Folkesson
//    
#ifndef CURE_RELATIVEPOINTFEATUREHELPER_HH
#define CURE_RELATIVEPOINTFEATUREHELPER_HH

#include "RelativePointFeature.hh"
#include "PosedRelativePointFeature.hh"
#include "PointFeatureHelper.hh"
#include "MapHelper.hh"
#include "Match.hh"
#include "VisionData.hh"
namespace Cure {

  /** 
   * Helper class for RelativePointFeatures You ned to use separte ones 
   * for bearing only, phiRange, and thetaRange type relative features.
   *  
   * It can be used to match different sorts of relative features.  It
   * can be used to do Lucas-Kanade tracking between two matched
   * frames So if you have matched two frames and want to interpolate
   * to the frames in between (without needing to compute complex
   * descriptors) this can do it.
   *
   * 'Normal' matching can also be done based on the bearing and if
   * available range to the point.  This does not require a descriptor
   * or a PoseTree but does require you to call match.  One must set
   * the approriate thresholds first.
   *
   * The 'ImagePlane' here must be compatable with the type of
   * RelativePointFeature:
   *
   * bearingOnly:
   * ImagePlane is (low phi, low theta, high phi, high theta)
   * note  that phi is  in (0,pi) and theta is (-pi,pi)
   * so if theta range is around pi then you are in trouble. 
   *
   * phiRange
   * ImagePlane is  (low phi, low range, high phi, high range)
   *
   * thetaRange
   * ImagePlane is  (low theta, low range, high theta, high range)
   *
   * else
   *
   * ImagePlane is (low x, low y, high x, high y)
   *
   * 
   * The match method works but there are some overides of the
   * functions it calls
   *   
   * The match method will do the whole matching job.  If no match is
   * found MatchDistance will be left at -1, if there is a near match
   * it will be set to -2 (This -w is to prevent using this to
   * initilaize a new feature near an existing one.) Otherwise it is
   * set to Eta*metric*Eta So this function first makes a near list
   * then a visable then a rough match then a fine match.
   *
   * PosedFeaatures returned on the matches are created for the
   * match and must be eventually deleted.  The Visable list has its
   * own PosedFeatures which it will delete when match is called
   * again.
   * 
   * The addDistance function is called which will convert the
   * MeasurmentType and add matching threshold information to the
   * matches.  The Near list is formed by a square box in the xy
   * plane 2*NearSearchRange on a side and centered on the sensor.
   *
   * The Visable list is formed by transforming all the near posed
   * features to the sensor frame and calling inRectangle on them
   * using ImagePlane as the rectangle.  That will behave
   * differently for different types of relative features @see
   * RelativePointFeature::inRectangle(...).
   *
   * The rough matching is done usign the bounding box that is
   * created for each match object.  It also works different for
   * different types of features.  The fine matching is done usign
   * the metric that is created for each match object.  It compute
   * an Eta that is the innovation of the measurement
   * @seePosedRelativePointfeature::calcEta.  It also works
   * different for different types of features.
   *
   *
   */
  class RelativePointFeatureHelper: public MapHelper 
  {
  public:
    PointFeatureHelper *m_PointHelper;
    /** Used for testCarteasian */
    double m_CartThreshold;
    /**In Pixels used for bearings only*/
    float m_FocalLength[2];
    /**Pixel xy of image center*/
    float m_BearingZero[2];
    /** this holdes an array of images to be used for matching and tracking.*/
    VisionData m_Images;
    /** 
     * The Current index into m_Images
     * (ie. m_Images.getShort(r,c,index))).
     */
    int m_Current;
    
    PosedRelativePointFeature * PosedRelativePointFeaturePtr;
    RelativePointFeature * RelativePointFeaturePtr;
    /** This can be set to how a new feature should look */
    RelativePointFeature TempletFeature;
    /** The key to the sensor MapPoseTree */
    long m_SensorKey;
    /** The current branch of the tree */
    short m_SensorBranch;
  protected:
    /**  Used to set match's metric for Eta used when calling fineMatch*/
    double m_BearingOnlyMetric[2];
    /**  Used to set match's metric for Eta used when calling fineMatch*/
    double m_PhiRangeMetric[3];
    /**  Used to set match's metric for Eta used when calling fineMatch*/
    double m_ThetaRangeMetric[2];
    /** This is used to set the bounding boxes for the various types.*/
    double m_ZBox[6];
  public:
    /** 
     * The Templet will be given these thresholds which then are used to 
     * initilize the MapFeatures made with this helper. 
     *
     * @param b the mapbank pointer should be provided.
     * @param bearingonly set true if these are bearign only type features.
     *
     * @param nonlineartheta set true if theta is transformed by cubic
     * polynomial.  Ignored if bearing only or thetarange is true.
     *
     * @param thetarange set true if bearingonly is false and the phi
     * mesurements are weak.  If false the theta or range are
     * considered 'weak'.
     *
     */
    RelativePointFeatureHelper(MapBank *b,
			       bool bearingonly=false,
			       bool nonlineartheta=false,
			       bool thetarange=false,	
			       double scale=130,
			       double alpha=1E-10,
			       double distanceThreshold=20)
      :MapHelper(b)
    {
      m_PointHelper=0;
      m_CartThreshold=.2;
      Bank=b;
      TempletFeature.setBank(b);
      TempletFeature.CastPtr=&RelativePointFeaturePtr;
      TempletFeature.DistanceThreshold=distanceThreshold; 
      if (bearingonly)
	TempletFeature.m_Function=(TempletFeature.m_Function|0x101);
      else if (thetarange)
	TempletFeature.m_Function=(TempletFeature.m_Function|0x400);
      else if (nonlineartheta)
	TempletFeature.m_Function=(TempletFeature.m_Function|0x202);
      else TempletFeature.m_Function=(TempletFeature.m_Function|0x200);
      TempletFeature.m_Scale=scale;
      TempletFeature.m_Alpha=alpha;      
      FeatureType=RELATIVEPOINTFEATURE_TYPE;
      ObjectSubType=FeatureType;
      m_FocalLength[0]=0;
      m_FocalLength[1]=0;
      m_BearingZero[0]=0;
      m_BearingZero[1]=0;
      m_Current=0;
      m_SensorKey=-1;
      m_SensorBranch=0;
      m_UseImagePlane=true;
      m_BearingOnlyMetric[0]=10000; 
      m_BearingOnlyMetric[1]=10000; 
      m_PhiRangeMetric[0]=10000;
      m_PhiRangeMetric[1]=25;
      m_PhiRangeMetric[2]=100;
      m_ThetaRangeMetric[0]=100;
      m_ThetaRangeMetric[1]=10000;
      m_ThetaRangeMetric[2]=100;
      m_ZBox[0]=20/sqrt(m_BearingOnlyMetric[0]);
      m_ZBox[1]=20/sqrt(m_BearingOnlyMetric[1]);
      m_ZBox[2]=20/sqrt(m_PhiRangeMetric[0]);
      m_ZBox[3]=20/sqrt(m_PhiRangeMetric[2]);
      m_ZBox[4]=20/sqrt(m_ThetaRangeMetric[1]);
      m_ZBox[5]=20/sqrt(m_ThetaRangeMetric[2]);
    }
    /**
     * This is used for fine matching of bearing only type features.
     * It is used as the diagonal of the Eta information matrix used
     * for fine matching of features.  It also will set the rought
     * matching bounding box to 20 std deviations of this.
     * 
     * @param metric the values to set the metric to.
     */
    void setBearingOnlyMetric(double metric[2]){
      m_BearingOnlyMetric[0]=metric[0];
      m_BearingOnlyMetric[1]=metric[1];
      m_ZBox[0]=20/sqrt(m_BearingOnlyMetric[0]);
      m_ZBox[1]=20/sqrt(m_BearingOnlyMetric[1]);
    }
    /**
     * Get method.
     * @param metric the values to returned in this.
     */
    void getBearingOnlyMetric(double metric[2]){
      metric[0]=m_BearingOnlyMetric[0];
      metric[1]=m_BearingOnlyMetric[1];
    }
    /**
     * This is used for fine matching of phi/range  type features.
     * It is used as the diagonal of the Eta information matrix used
     * for fine matching of features.  It also will set the rough
     * matching bounding box to 20 std deviations of this.  The rough
     * matching is done in the phi and range directions only.
     * 
     * @param metric the values to set the metric to.
     */
    void setPhiRangeMetric(double metric[3]){
      m_PhiRangeMetric[0]=metric[0];
      m_PhiRangeMetric[1]=metric[1];
      m_PhiRangeMetric[2]=metric[2];
      m_ZBox[2]=20/sqrt(m_PhiRangeMetric[0]);
      m_ZBox[3]=20/sqrt(m_PhiRangeMetric[2]);
    }
    /**
     * Get method.
     * @param metric the values to returned in this.
     */
    void getPhiRangeMetric(double metric[3]){
      metric[0]=m_PhiRangeMetric[0];
      metric[1]=m_PhiRangeMetric[1];
      metric[3]=m_PhiRangeMetric[3];
    }
    /**
     * This is used for fine matching of theta/range  type features.
     * It is used as the diagonal of the Eta information matrix used
     * for fine matching of features.  It also will set the rough
     * matching bounding box to 20 std deviations of this.  The rough
     * matching is done in the theta and range directions only.
     * 
     * @param metric the values to set the metric to.
     */
    void setThetaRangeMetric(double metric[3]){
      m_ThetaRangeMetric[0]=metric[0];
      m_ThetaRangeMetric[1]=metric[1];
      m_ThetaRangeMetric[2]=metric[2];
      m_ZBox[4]=20/sqrt(m_ThetaRangeMetric[1]);
      m_ZBox[5]=20/sqrt(m_ThetaRangeMetric[2]);
    }
    /**
     * Get method.
     * @param metric the values to returned in this.
     */
    void getThetaRangeMetric(double metric[3]){
      metric[0]=m_ThetaRangeMetric[0];
      metric[1]=m_ThetaRangeMetric[1];
      metric[3]=m_ThetaRangeMetric[3];
    }
    MapPoseTree * getSensorTree(){
      return getMapPoseTree(Bank,m_SensorKey);
    }
    /**
     * Get method.
     * @return the pointer to the tree sensor leaf's transformation.
     */
    Transformation3D * sensorPose(){
      MapPoseTree *pt=getSensorTree();
      if (pt) return pt->getLeafPose(m_SensorBranch);
      return 0;
    }
    /**
     * Get method.
     * param pose the sensor pose is returned. 
     * @return 0 if ok else 1
     */
   int  getSensorPose(Pose3D &pose){
      MapPoseTree *pt=getSensorTree();
      if (pt){
	unsigned short type=0;
	Transformation3D p;
	pt->getLeafPose(p,type,
			m_SensorBranch);
	pose=p;
	pose.setCovType(type);
	pose.setTime(pt->time());
       	return 0;
      }
      return 1;
    }

    
    /**
     * Get method.
     * @return the pointer to the tree branch's Descriptors.
     */
    FeatureDescriptors *getSensorDescriptors(){
      MapPoseTree *pt=getSensorTree();
      if (pt)
	return pt->getDescriptors(m_SensorBranch);
      return 0;
    }


    ~RelativePointFeatureHelper(){}
    virtual void printConfiguration();
    
    /*
     * @param focal the focal lengths xy for finding bearings
     * @param center the pixel of the camera pose direction.
     */
    void setParameters(float focal[2],float center[2]){
      m_FocalLength[0]=focal[0];
      m_FocalLength[1]=focal[1];
      m_BearingZero[0]=center[0];
      m_BearingZero[1]=center[1];
    }
    PosedRelativePointFeature * castPosedRelativePointFeature(PosedFeature *pf){
      PosedRelativePointFeaturePtr=0;
      pf->narrow();
      return PosedRelativePointFeaturePtr;
    }
    RelativePointFeature * castRelativePointFeature(MapObject *pf){
      RelativePointFeaturePtr=0;
      pf->narrow();
      return RelativePointFeaturePtr;
    }
    RelativePointFeature *getRelativePointFeature(long key)
    {
      MapObject *mo=Bank->getMapObject(key);
      if (!mo) return 0;
      return castRelativePointFeature(mo);
    } 
    virtual MapFeature *getMapFeature(long key){
      return getRelativePointFeature(key);
    }

    virtual MapFeature * makeMapFeature(MapBank *b=0){
      return new RelativePointFeature(&TempletFeature,b);
    }

    /**
     * @return A RelativePointFeature based on the TempletFeature. 
     */
    RelativePointFeature * makeRelativePointFeature(){
      return new RelativePointFeature(&TempletFeature,Bank);
    }

    /**
     * 
     * @param f The FeatureDescriptor in the reference frame
     *  
     * @return A RelativePointFeature based on the TempletFeature and f.
     */
    RelativePointFeature * makeRelativePointFeature
    (Cure::FeatureDescriptor *f) 
    {    
      double dist=f->m_Range;
      if (dist==0)dist=CURE_FEATURE_NOMINAL_RANGE;
      RelativePointFeature *temppt=new RelativePointFeature(&TempletFeature,
							    Bank);
      double xrel[3];
      xrel[0]=f->m_Bearing(0)*dist;
      xrel[1]=f->m_Bearing(1)*dist;
      xrel[2]=f->m_Bearing(2)*dist;
      temppt->setCenter(xrel);
      temppt->addDescriptor(f);      
      return temppt;
    }     

    virtual PosedFeature * makePosedFeature(MapFeature* mf){
      return makePosedRelativePointFeature(mf);}
    virtual MapFeature * makeMapFeature(Cure::Pose3D &sensorpose,
					Cure::Measurement &m,
					bool addToVisable=false);

    //  virtual MapFeature * makeMapFeature(Cure::Pose3D sensorpose,
    //			      double distanceguess=2.0){
    //  return makeRelativePointFeature(sensorpose,distanceguess);}
    PosedRelativePointFeature * makePosedRelativePointFeature(MapFeature* mf){
      if(mf->Type==RELATIVEPOINTFEATURE_TYPE)
	return new PosedRelativePointFeature
	  (&PosedRelativePointFeaturePtr,
	   castRelativePointFeature(mf));
      return 0;
    }
    /**
     * This finds the first good matchs from a FeatureDescriptors
     * object's discriptors to a list of other FeatureDescriptors.
     * The FeatureDescriptors and RelativePointFeatures matched will
     * be set to reference one another.  New features will be made for
     * firstime matches. 
     *
     * This will use both the feature descriptor match and the
     * epipolar match.  After calling this all feature descriptors on
     * currentdes will have there m_FeatureKey set to either a new
     * feature or an existing one.
     *
     * @param currentdes the FeatureDescriptors associtated with the 
     * current pose that you are trying to match to earlier poses.
     *
     * @param des an array of earlier FeatureDescriptors.  The order
     * should be closest one to the current first.
     *
     * @param n the length of des
     *
     * @param errorlimit applied to the return from
     * FeatureDescriptors::match(...)
     * @param radsqerr estimated variance in the bearing angle
     * @param coarseEnergyLimit the threshold on the epipolarconstraint.
     * @return the number of matches found.
     * 
     */
    int findMatches(FeatureDescriptors &currentdes,
		    FeatureDescriptors **des,
		    int n,
		    float errorlimit, 
		    double radsqerr, 
		    double coarseEnergyLimit);
    /**
     * This gets a short value for a pixel of an image.
     * @param r the row index into the images array.
     * @param c the column index into the images array.
     * @param imgindex the index of the image.
     * @return the pixel value as a shsort. 
     */
    virtual short getPixelValue
    (const unsigned short r, 
     const unsigned short c, const unsigned short imgindex=0)
    {
      return m_Images.getShort(r,c,imgindex);
    } 
    virtual bool getC(MapFeature  *mf1,MapFeature  *mf2,
		      Cure::Matrix &c1){
      RelativePointFeature *mw1=castRelativePointFeature(mf1);
      if (!mw1)return false;
      RelativePointFeature *mw2=castRelativePointFeature(mf2);
      if (!mw2)return false;
      return mw1->getC(mw2,c1);
    }
    virtual bool testMatch(MapFeature  *mf1,MapFeature  *mf2,
			   double tolerance){
      RelativePointFeature *mw1=castRelativePointFeature(mf1);
      if (!mw1)return false;
      RelativePointFeature *mw2=castRelativePointFeature(mf2);
      if (!mw2)return false;
      return mw1->testMatch(mw2,tolerance);
    }
    void match(Match *matches,int n,
	       double distance,long sensorKey, short branch, 
	       LongList &mapFeats){
      m_SensorKey=sensorKey;
      m_SensorBranch=branch;
      Pose3D psensor;
      TempletFeature.m_ReferenceKey=sensorKey;
      TempletFeature.m_RefBranch=branch;
      if (getSensorPose(psensor))return;
      MapHelper::match(matches,n,distance,psensor,mapFeats);
    }

    virtual int config(const std::string &arglist);

    virtual bool supportsSubconfig(int sc);
    
    /**
     * Applys tests and if pass makes a world frame MapPointFeature
     * that is equal to the RelativePointFeature.  The test is on the
     * one sigma spread in cartesian space.  if any spread of the xyz
     * in the reference ref-frame exceed m_CartThreshold no point is made.
     * 
     * @jac if replaced the jacobian of the returned feature wrt
     * refernace tree cov coordinates and the rf->Scalars
     * @rf the relative feature to replace
     * @param info the information matix diagonal for the rf Scalars.
     * @return if converged a MapPointFeature that can replace rf.
     */
    MapPointFeature * testCartesian(Matrix &jac,RelativePointFeature *rf, 
				    double info[3])
    {
      if (!rf)return 0;
      if (!m_PointHelper)return 0;
      if (rf->Bdual.Columns<3)return 0;
      double test;
      double jpolar[3];
      int rows=rf->getPolarJac(jpolar);
      double r=rf->getRange();
      if (r<1E-6)return 0;
     
      test=r/sqrt(info[0]);
      if (test>m_CartThreshold)return 0;
      test=r*jpolar[1]*jpolar[1]/info[1];
      if (test>m_CartThreshold)return 0;
      if (rows&4)
	test=r*r*r/info[2];
      else test=sqrt(info[2]);
      if (test>m_CartThreshold)return 0;
      Matrix x;
      if (rf->getCartesian(jac,x))return 0;
      MapPointFeature *mf=m_PointHelper->makeMapPointFeature();
      mf->forceExtend();
      mf->setX(x);
      return mf;
    }

    /**
     * this adds FeatureDescriptor objects for in between poses on a
     * list of FeatureDescriptors by doing lucas-kanade matching on
     * the frames.
     * 
     * @param des an array of FeatureDescriptors in cronological order
     * with the first and last having been matched and the in-between lists
     * being empty. 
     * @param n the length of des
     * @param windowsize size of window is 2*windowsize+1
     * @param offset the number of images the current image is
     * beyond the image of des[n-1].
     * @return the number of added descriptors.
     */
    int lucasKanadeTrack(FeatureDescriptors **des, 
			 int n,
			 unsigned short windowsize=7,
			 int offset=0);
    
    /**
     * @param pix the center of the window
     * @param windowsize size of window is 2*windowsize+1
     * return true if window is within the image
     */
    bool isWithin(short pix[2], unsigned short windowsize){
      int t=pix[0]-windowsize;
      if ((t<0)||(t>=m_Images.getWidth()))return false;
      t=pix[0]+windowsize;
      if ((t<0)||(t>=m_Images.getWidth()))return false;
      t=pix[1]-windowsize;
      if ((t<0)||(t>=m_Images.getHeight()))return false;
      t=pix[1]+windowsize;
      if ((t<0)||(t>=m_Images.getHeight()))return false;
      return true;
    } 
    /**
     * This will do the lucas-kanade tracking when the subclass overrides
     * gettemplet anad getDisplacement..
     *
     * @param originalindex the index of the frame @param trackedindex
     * the index of the frame @param originalpix the pixel of the
     * feature @param trackedpix the pixel estimated of the feature
     * called with a rough idea and returned with an improved one.
     * @param displace returns with the fractional part of the new
     * trackedpixel location.
     * @param windowsize size of window is 2*windowsize+1
     *
     * @return the square of the interpolated pixel 'error' or -1 if fails. 
     */
    float lktrack(int originalindex, int trackedindex, 
		  short originalpix[2],
		  short trackedpix[2],
		  Cure::Matrix & displace,
		  unsigned short windowsize=7);
    /**
     * Gets a templet of the pixels around a center pixel. 
     * @param imageindex the index of the image in some array of images
     * @param templet called with the templet size set to the window size
     * returned with the templet set to the pixels in the window.
     * @param centerx the center of the window 
     * @param centery the center of the window 
     * @return 0 if success else 1.
    */
    void getTemplet(int imageindex,
		   Cure::ShortMatrix &templet,
		   short centerx,short centery){
      int off=(templet.Rows-1)/2;
      int iw=(centerx-off);
      off=(templet.Columns-1)/2;
      for (int ti=0;ti<templet.Rows;ti++,iw++)
	{
	  int j=centery-off;
	  for (int tj=0;tj<templet.Columns;tj++,j++)
	    templet(ti,tj)=getPixelValue(iw,j,imageindex);
	}
    }
    /**
     * gets a short out of the pixel from an imagae of the image array 
     *
     * @param imgindex the index into the image array.
     * @param pixelindex pixel row x m_ImageWidth + pixel cloumn.
     * @return the value as a short. 
     */
    /**
     *This does the lucas-kanade tracking between templet and img2.
     *
     *
     * @param center2
     * @return the square difference in pixel values interpolated to the
     * displace location or -1 if fails.
     */
    float getDisplacement(int imageindex,Cure::Matrix &displace,
			  Cure::ShortMatrix &templet,
			  short center2[2]);
    
    /**
     * This compute Z, V, BoundingBox, Thresholds and Metric  
     *
     *
     * for (TempletFeature.bearingOnly():
     * The structure of the input Measurements are assumed to be 
     * SensorType=SensorData::SENSORTYPE_CAMERA=2
     * MeasurementType=32,33
     * V: 3x1 (x,y,f) // pixels
     * Z: 2x1 (x,y) //pixels
     * W: 1x10 (imwidth, imheight, f, pixwidth, pixheight, 
     *          x_pix, y_pix, xsize_pix, ysize_pix, angle)
     * 
     * or
     *     
     * MeasurementType=34
     * V: 3x1 bx,by,bz unit bearing
     * Z: 3x1 phi,theta,,
     * W: 1x10 (imwidth, imheight, f, pixwidth, pixheight, 
     *          x_pix, y_pix, xsize_pix, ysize_pix, angle)
     *
     *
     * The output has
     * MeasurementType=34,35  (35 is same as 34 except it is forced to extend)
     * This has Z=phi,theta
     * V=unit bearing, 
     *
     *
     * For TempletFeature.phiRange():
     * Z=(phi,theta, range
     * V=(phi,f(theta)range).
     *
     *
     * @param mat the Match is returned in this.
     * @param distance The distance to add info with.
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
      RelativePointFeature *mf0=getRelativePointFeature(pf[0]->FeatureKey);
      RelativePointFeature *mf1=getRelativePointFeature(pf[1]->FeatureKey);
      if (!mf1)return -1;
      if (!mf0)return -1;
      int r=mf1->merge(mf0,typ);  
      if (typ&1){
	delete mf0;
      }
      return r;
    }
    static void 
    makeBearings(Cure::Matrix &bearings,bool *frames,
		 Cure::Matrix &v,
		 Cure::Transformation3D  *sensorposes,
		 const double distanceguess=5);



    /**
     * Precalculate the path and deltapath for triangulatePoint. This
     * function is typically used before calling triangulatePoint.
     * Here we are updating the matricies from the triangulation of the 
     * previous frame buffer.  So we remove the oldest (index) and 
     * replace it with the new pose.  This is for a circular buffer
     * that avoid copying.
     *
     * The function will grow the path and deltapath matricies if
     * index is beyont their size and replace if index is within them.
     *
     * @param path the (x,y,z) of camera when bearings taken in 
     *             the info frame (output)
     * @param deltapath the differences in (x,y,z) of camera (output)
     * @param sensorposes the new camera pose
     * @param index the number of the new frame. 
     */
    static void 
    insertPoseTopath(Cure::Matrix &path, 
		     Cure::Matrix &deltapath,
		     const Transformation3D  &sensorpose,
		     int index);
    /**
     * Calculate the predicted image plane coordinates of the meanpoint.
     * @param vhat returns with (x,z,focallength) for all poses that have
     *        vhat(i,2)=focallength>0 when called.  Use the third column
     *        to indicate which poses to use.
     * @param sensorposes An array [n] of the camera poses in the info frame.
     * @param meanpoint the (xyz) of the point in the world frame.
     */
    static void 
    predictPixels(Cure::Matrix &vhat,
		  Cure::Transformation3D  *sensorposes,
		  Cure::Matrix &meanpoint);

    /**
     * This creats a time ordered matrix of innovations with time in 
     * whole frames from oldest frame assuming a cirular buffer.
     * The vhat(i,2)>0 indicates which frames have data.
     * @param innovation (v(i,0)-vhat(i,0), v(i,1)-vhat(i,1),time in frames).
     * @param oldestframe The index of oldest frame 
     */
    static void 
    calcInnovations(Cure::Matrix &innovation,
		    Cure::Matrix &vhat,
		    Cure::Matrix &v, int oldestframe);
    /**
     * @param dx_dt (x(i+1,0)-x(i,0)/x(i+1,2)-x(i,2),
     *               x(i+1,0)-x(i,0)/x(i+1,2)-x(i,2),
     *               x(i+1,2)+x(i,2)/2,
     */
    static void calcTimeDerivative(Cure::Matrix &dx_dt, Cure::Matrix &x);

    /**
     * Triangulates a set of point observations. Use this function to
     * find an estimate of a point feature given a number of
     * bearing-only measurements to it. The measurements are in the form
     * of camera poses, image coordinates and focal length, one for each
     * frame where the feature has been detected. The output is the
     * meanpoint which is the estimate of the point and the return value
     * which gives a measure of how accurate the triangulation.
     *
     * The path and delta can be calculated for a set of frames using calcPath 
     * or growPath.  Then these same matricies can be reused for repeated calls
     * to this function.
     *
     * So if frames[i]==true then bearings(i,k) has a 
     * bearing to the point at path(i,k). 
     *
     * You can pass on pointers to variables to get hold of the value
     * for the different parameters that are checked against the
     * thresholds. This gives you a way of "tuning" by looking at what
     * triangulation points are good and then looking at what values
     * these
     * 
     * @param meanpoint returns the (3,1) column of the x,y,z of the point 
     *                  in the info frame (putput)
     * @param path the matrix (n,3) of camera (x,y,z) positions
     * @param deltapath the (n,3*n) matrix of camera displacements
     * @param bearings the (n,3) matrix of measurements (u,v,focallength)
     * @param weightThreshold used to test sum of sin^2(theta),angle 
     *                       between bearings.
     * @param sqDistanceThreshold  used to test distance of each bearing 
     *                             ray from the meanpoint.
     * @param mseThreshold max mean square error of all the bearings 
     *                     from the meanpooint.
     *
     * @param retWeight pointer to storage variable of sqDist used to
     * compare with the threshold weightThreshold
     * @param retSqDist pointer to storage variable to store max sqDist
     * @param retMse pointer to storage variable of mse used to
     * compare with the threshold mseThreshold
     * @param retSecDer pointer to storage variable of secDer used to
     * compare with the threshold secDerThreshold
     *
     * @return -1 if no baseline, sum sin^2 angle < threshold
     *         -2 if some bearing's best point was further than 
     *            consistencytest from meanpoint
     *         -3 if mse >consistencytest/sqrt(n) mse=mean square error of 
     *            point cloud from meanpoint
     *         else sum(sin^2(theta)) triangulation quality 
     *     
     */
    static double 
    triangulateBuffer(Cure::Matrix &meanpoint,
		      Cure::Matrix &path, 
		      Cure::Matrix &deltapath, 
		      Cure::Matrix &bearings,
		      bool *frames,
		      const double weightThreshold=.1, 
		      const double sqDistanceThreshold=.01,
		      const double mseThreshold=.001,
                      double *retWeight = 0,
                      double *sqDist = 0,
                      double *mse = 0);

    /**
     * Triangulates a set of point observations. Use this function to
     * find an estimate of a point feature given a number of
     * bearing-only measurements to it. The measurements are in the form
     * of camera poses, image coordinates and focal length, one for each
     * frame where the feature has been detected. The output is the
     * meanpoint which is the estimate of the point and the return value
     * which gives a measure of how accurate the triangulation.
     *
     * Here we calculate pairwise and then take a weighted aveerage.
     * The complexity is n^2 as opposed to n^3 for the min distance form.
     *
     * The measurement matrix v can have rows of 0's if so the point is 
     * ignored (v(i,2)==0, is the condition to ignore this sensorpose. 
     *
     * Thus one can have a fixed, (even circular) buffer of poses and 
     * just be sure that the rows of v(i,j) coorespond to the sensorpose[i].
     * Unmatched frames will then have v(i,2)=0; and are ignored.
     *
     * Upon return all frames not used will have v(i,2)<=0. 
     * The orginal empty rows will be unchanged (ie.0)
     * The detected outliers will have v(i,2)==-1.
     *
     * That is to say some of the bearings will not pass the 
     * sqdistanceThreshold, (not point at the meanpoint) and will not
     * be used.  They haave v(i,2) set to -1.  THis information can also 
     * be used by the calling object.
     *
     * You can pass on pointers to variables to get hold of the value
     * for the different parameters that are checked against the
     * thresholds. This gives you a way of "tuning" by looking at what
     * triangulation points are good and then looking at what values
     * these
     * 
     * @param meanpoint returns the (3,1) column of the x,y,z of the point 
     *                  in the info frame (putput)
     * @param path the matrix (n,3) of camera (x,y,z) positions
     * @param deltapath the (n,3*n) matrix of camera displacements  
     * @param v the (n,3) matrix of measurements (u,v,focallength)
     * @param sensorposes An array [n] of the camera poses in the info frame.
     * @param weightThreshold used to test sum of sin^2(theta),angle 
     *                       between bearings.
     * @param sqDistanceThreshold  used to test distance of each bearing 
     *                             ray from the meanpoint.
     * @param mseThreshold max mean square error of all the bearings 
     *                     from the meanpooint.
     * @param secDerThreshold threshold on the second derivateive if <=0 
     *                     this test will be ignored
     * @param distanceguess default 5 (units of v)this is used to get bearings 
     *
     * @param retWeight pointer to storage variable of sqDist used to
     * compare with the threshold weightThreshold
     * @param retSqDist pointer to storage variable to store max sqDist
     * @param retMse pointer to storage variable of mse used to
     * compare with the threshold mseThreshold
     * @param retSecDer pointer to storage variable of secDer used to
     * compare with the threshold secDerThreshold
     * 
     * @return -1 if no baseline, sum sin^2 angle < threshold
     *         -2 if some bearing's best point was further than 
     *            consistencytest from meanpoint
     *         -3 if mse >consistencytest/sqrt(n) mse=mean square error of 
     *            point cloud from meanpoint
     *         -4 if second derivative is two large 
     *         else sum(sin^2(theta)) triangulation quality 
    */
    static double 
      triangulateBuffer(Cure::Matrix &meanpoint,
			Cure::Matrix &path, 
			Cure::Matrix &deltapath, 
			Cure::Matrix &v,
			Cure::Transformation3D  *sensorposes,
			const double weightThreshold=.1, 
			const double sqDistanceThreshold=.01,
			const double mseThreshold=.001,
                        const double secDerThreshold=1e-4,
			const double distanceguess=5,
                        double *retWeight = 0,
                        double *retSqDist = 0,
                        double *retMse = 0,
                        double *retSecDer = 0);
    
    
    /**
     * Same as the above function but requires only the sensorposes
     * array. The path and delatpath are derived automatically from
     * this. This is of course more inefficient that calculating the
     * path and deltapath once and for all.
     *
     * Triangulates a set of point observations. Use this function to
     * find an estimate of a point feature given a number of
     * bearing-only measurements to it. The measurements are in the form
     * of camera poses, image coordinates and focal length, one for each
     * frame where the feature has been detected. The output is the
     * meanpoint which is the estimate of the point and the return value
     * which gives a measure of how accurate the triangulation.
     *
     * The measurement matrix v can have rows of 0's if so the point is 
     * ignored (v(i,2)==0, is the condition to ignore this sensorpose. 
     *
     * Thus one can have a fixed, (even circular) buffer of poses and 
     * just be sure that the rows of v(i,j) coorespond to the sensorpose[i].
     * Unmatched frames will then have v(i,2)<=0; and are ignored.
     *
     * Upon return all frames not used will have v(i,2)=0.
     * That is to say some of the bearings will not pass the 
     * sqdistanceThreshold, (not point at the meanpoint) and will not
     * be used.  They haave v(i,2) set to -1.  THis information can also 
     * be used by the calling object.
     *
     * You can pass on pointers to variables to get hold of the value
     * for the different parameters that are checked against the
     * thresholds. This gives you a way of "tuning" by looking at what
     * triangulation points are good and then looking at what values
     * these
     * 
     * @param meanpoint returns the (3,1) column of the x,y,z of the point 
     *                  in the info frame (putput)
     * @param v the (n,3) matrix of measurements (u,v,focallength)
     * @param sensorposes An array [n] of the camera poses in the info frame.
     * @param weightThreshold used to test sum of sin^2(theta),angle 
     *                       between bearings.
     * @param sqDistanceThreshold  used to test distance of each bearing 
     *                             ray from the meanpoint.
     * @param mseThreshold max mean square error of all the bearings 
     *                     from the meanpooint.
     * @param secDerThreshold threshold on the second derivateive if <=0 
     *                     this test will be ignored
     * @param distanceguess default 5 (units of v)this is used to get bearings 
     * 
     * @return -1 if no baseline, sum sin^2 angle < threshold
     *         -2 if some bearing's best point was further than 
     *            consistencytest from meanpoint
     *         -3 if mse >consistencytest/sqrt(n) mse=mean square error of 
     *            point cloud from meanpoint
     *         -4 if second derivative is two large 
     *         else sum(sin^2(theta)) triangulation quality 
    */
   static double 
   triangulateBuffer(Cure::Matrix &meanpoint,
		     Cure::Matrix &v,
		     Cure::Transformation3D  *sensorposes,
		     const double weightThreshold=.1, 
		     const double sqDistanceThreshold=.01,
		     const double mseThreshold=.001,
                     const double secDerThreshold=1e-4,
		     const double distanceguess=5);

    

  protected:

    static double makeDerivateCheck(Cure::Matrix &meanpoint,
                                    Cure::Matrix &v,
                                    Cure::Transformation3D *sensorposes);
    
    /**
     * Reads configurations of version 1
     */
    int configVer1(const std::string &arglist);
  };

} // namespace Cure

#endif
