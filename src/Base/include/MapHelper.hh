// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    

#ifndef CURE_MAPHELPER_H
#define CURE_MAPHELPER_H

#include "PosedFeatureList.hh"
#include "MapFeatureList.hh"
#include "FeatureFun.hh"
#include "PosedFeatureList.hh"
#include "Measurement.hh"
#include "Match.hh"

#include <string>

namespace Cure{


  /**
   * Base class for map helper classes that do the real work for each
   * of the feature types in the map.  These work with the
   * PosedFeature and MapFeature Classes to implement data association
   * (matching) and feature creation.
   *
   * @author John Folkesson
   */
  class MapHelper: public MapObjectMaker
  {
  public:
    /** The type of featue this helps */
    short FeatureType;
    /** the list of FeatureKeys of the near list*/
    LongList NearKeys;
    /** A HAsh table for tracked features */
    LongPairList *TrackKeys;
    /** The size of the Hash table*/
    unsigned short TrackHashSize;
    /** 
     * List of map features that are near the predicted pose of the
     * robot. The list holds pointers to these map features.
     */
    MapFeatureList Near;

    /** 
     * List of pointers to PosedFeatures. These are created in the
     * match function in every iteration. You need to be careful to
     * delete the memory. This is done every time in the
     * makePosedFeatures function. You want to explicitly call clear()
     * on this list in any destructor that uses a MapHelper class 
     */
    PosedFeatureList PosedFeatures;

    /** 
     * Fatures that are visible. The list just contains pointers into
     * the PosedFeatures PosedFeatureList 
     */
    PosedFeatureList Visable;

   
    /**
     *  The near list is everything in a square centered on robot and side
     *  2 times NearSearchRange in meters.
     */
    double NearSearchRange;
    /**
     * added to BoundingBox during matchFeatures to do rough match, also
     * used to set the thresholds when addDistance is called with a Match. 
     */
    double RoughSearchRange;
    /**
     * if MatchDistance is between this and MatchThreshold it is set to -2.
     *  This is used to discard measurements near but not matching features.
     */
    double NewFeatureRange;
    /**
     *  This defines a match for matchFeature 
     *  transpose(delZ)*metric*delZ<MatchThreshold
     */
    double MatchThreshold;
    /**
     *  This is the defining rectangle for finding match candidates.
     *  Any pixel prediction with x or y outside this will be considered
     *  out of view.  This should be a little larger than the actual camera
     *  field of view. (minx,miny,maxx,maxy)
     */
    double ImagePlane[4];
    /**
     * PixelInfo[0]=focalLength;
     * PixelInfo[1]=pixwidth; 
     * PixelInfo[2]=pixheight;
     */
    double PixelInfo[3];
    /** 
     * Set this to true to force the imagge plane to be checked when
     * forming the visable list.
     */
    bool m_UseImagePlane;
    //This needs to get into measurement
    int  SensorType;
    MapHelper(MapBank *b=0):MapObjectMaker(b){
      ObjectType=FEATURE_TYPE;
      ObjectSubType=0;
      NearSearchRange=1;
      RoughSearchRange=.5;
      SensorType=0;
      TrackKeys=0;
      TrackHashSize=0;
      FeatureType=0;
      m_UseImagePlane=false;
    }
    virtual ~MapHelper(){
      if (TrackKeys!=0)delete[]TrackKeys;
      TrackHashSize=0;
      TrackKeys=0;
    }
    virtual MapFeature * getMapFeature(long){return 0;}
    virtual PosedFeature * makePosedFeature(MapFeature* ){return 0;}
    virtual MapFeature * makeMapFeature(Cure::Pose3D &,
					Cure::Measurement &, 
					bool addToVisable=false){return 0;}

    virtual MapFeature * makeMapFeature(MapBank *b=0){return 0;}
    MapObject * makeMapObject(){return makeMapFeature();}
    void setTrackHashSize(unsigned short ths)
    {
      if (ths<=0)return;
      if (TrackHashSize!=0)delete[]TrackKeys;
      TrackKeys=new LongPairList[ths];
      TrackHashSize=ths;
    }
    /**
     * The Tracked key allows external match hints to be used by the Helper.
     * @param tkey the tracked key is the Measure->Key if any that was
     * set when this feature was matched last.  It is used for any
     * external matching/tracking of the feaures and allows bypassing
     * the matching in the MapHelper.
     * @param fkey the FeauterKey associated with the tracked key.
     *
     */
    void addTrackedKey(unsigned long tkey, long fkey)
    {
      if (!TrackKeys)setTrackHashSize(200);
      unsigned short n=tkey%TrackHashSize;
      TrackKeys[n].replace(tkey,fkey);
    } 
    /*
     * Finds the associated feature.
     * @param tkey the tracked key is the Measure->Key if any that was
     * set when this feature was matched last.  It is used for any
     * external matching/tracking of the feaures and allows bypassing
     * the matching in the MapHelper.
     * @return the associated FeatureKey or -1 if not found.
     */
    long findTrackedKey(unsigned long tkey)
    {
      if (TrackHashSize==0)return -1;
      unsigned short n=tkey%TrackHashSize;
      return TrackKeys[n].find(tkey);
    } 

    /*
     * Finds the associated feature.
     * @param tkey the tracked key is the Measure->Key if any that was
     * set when this feature was matched last.  It is used for any
     * external matching/tracking of the feaures and allows bypassing
     * the matching in the MapHelper.
     * @return the associated MapFeature or 0 if not found.
     */
    virtual MapFeature * getTrackedFeature(unsigned long tkey){
      if (TrackHashSize==0)return 0;
      unsigned short n=tkey%TrackHashSize;
      long fkey=TrackKeys[n].find(tkey);
      if (fkey!=-1)
	return getMapFeature(fkey);
      return 0;
    }
    virtual void printConfiguration(){}
    /**
     * This increases (or decreases for factor<1) the threshold for 
     * makeing the preliminary match to features..
     * @param factor one is no change, >1 makes it easier to match lines,
     *               <1 makies it harder. 
     */  
    void loosenMatching(double factor)
;
    /**
     * Get the matricies c1 and c2 so that delta=(c1 B1 x1 - c2 B1 x2) is
     * a meaningfull measure of distance between the 2 features.
     * Here B1 is mf1's B matrix and x1 is its state X_f vector.
     * Notice that delta can have less dimensions that B1 has rows. 
     * c2 is rescaled for the difference in lengths for angle dimensions.
     * Also the c2 angle row is multiplied by -1 if the start and end are 
     * defined opposite.
     * c1 is not so c1^T c1 is like a projection to the common subspace.   
     * @return true if any such measure exsists.
     */
    virtual bool getCommon(MapFeature  *mf1,MapFeature  *mf2,
			    Cure::Matrix &c1, Cure::Matrix &c2);
    /**
     * Tests the non-common dimensions for consistancy within
     * tolerance, a distance error in meters;
     */ 
   virtual bool testMatch(MapFeature  *,MapFeature  *,
			    double tolerance){return false;}
    /**
     * Tests the non-common dimensions for consistancy within
     * tolerance, a distance error in meters;
     */ 
   virtual double judgeMatch(MapFeature  *,MapFeature  *,
			    double ){return 0;}
    /*
     * This helps getCommon by calculating the c1 part.
     */
    virtual bool getC(MapFeature  *, MapFeature  *,
			Cure::Matrix &){return false;}

    /**
     * Configure this helper given a list of arguments in a string
     */
    virtual int config(const std::string &arglist) { return 0; }


    virtual bool supportsSubconfig(int subconfig) = 0;

    /**
     *  This will do the whole matching job.  If no match is found
     * MatchDistance will be left at -1, if there is a near match it
     * will be set to -2 (This -w is to prevent using this to
     * initilaize a new feature near an existing one.) Otherwise it is
     * set to dz*metric*dz So this function first makes a near list
     * then a visable then a rough match then a fine match.  The
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
     * differently for different types of features @see
     * MapFeature::inRectangle(...).
     *
     * The rough matching is done usign the bounding box that is
     * created for each match object.  It also works different for
     * different types of features.  The fine matching is done usign
     * the metric that is created for each match object
     * 
     * @param matches the measurements are here when called pointers to 
     *        PosedFeatures are added for all matched measurments.  
     * @param n the numner of measurements. 
     * @param distance the decorrelation distance along the path for addinfo.
     * @param mapFeats as list of the map features to be matched.
     * @param psensor the sensor pose.
     */
    virtual void match(Match *matches,int n,
		       double distance,Cure::Pose3D &psensor, 
                       LongList &mapFeats);

    /**
     * This creates a list 'Near' of MapFeatures within a square out
     * of a list of keys. The list just contains pointers to existing
     * MapFeatures from the map.
     * @param center the x and y for the center of the near square.
     */
    void makeNear(const double center[2],LongList * llist=0)
    {
      if (!llist)llist=&NearKeys;
      MapFeatureList temp;
      FeatureFcn::listFeatures(temp,*Bank,*llist);
      FeatureFcn::near(Near,temp, NearSearchRange, center);
      FeatureFcn::listLongs(Near,NearKeys);
    }

    /**
     * This takes MapFeatures from the list 'Near' and makes
     * PosedFeatures out of them and adds to 'PosedFeature' list. This
     * function thus allocates new memory when constructing the posed
     * features. It is important to delete the allocated objects. 
     */
    void makePosedFeatures()
    {
      PosedFeatures.clear();
      for (MapFeatureList *mlist=&Near; mlist->Next; mlist=mlist->Next)
	{
	  PosedFeature *pf=makePosedFeature(mlist->Element);
	  if (pf)PosedFeatures.add(pf);
	}
    }
    /**
     * This creates a list Visable of transformed PosedFeatures 
     * from PosedFeatures that are in front of sensor. 
     */
    void makeVisable(Cure::Pose3D& pose, int sensortype=0);

    /*
     * Form a List of non-occluded scan-visable features and the 
     * 360 degree predicted scan at one degree resolution.
     *
     * @param visFeatures returns with the scan-visable list.
     * @param scan returns with the predicted range scan.
     */
    void makeFeatureScan(PosedFeatureList &visFeatures,
			 double scan[360]);




    void setImage(double bottomLeft[2],double topRight[2],
		  double focalLength, double pixwidth, double pixheight)
    {
      PixelInfo[0]=focalLength;
      PixelInfo[1]=pixwidth;
      PixelInfo[2]=pixheight;
      ImagePlane[0]=bottomLeft[0];
      ImagePlane[1]=bottomLeft[1];
      ImagePlane[2]=topRight[0];
      ImagePlane[3]=topRight[1];
    }
    /**
     * 
     * @param data =(image width, image height, focal length, pix width,
     *                pix height) The first three are in units of pixels,
     *                The last 2 are the conversion from pixels to meters.
     */
    void setImage(double data[5])
    {
      double tr[2];
      tr[0]=data[0]*data[3]*.65; // 0.65 = 1.3 * 0.5, i.e 30% extra
      tr[1]=data[1]*data[4]*.65;
      double f=data[2]*data[3];
      setImage(tr,f,data[3],data[4]);
    }
    void setImage(double topRight[2],
		  double focalLength, double pixwidth, double pixheight)
    {
      PixelInfo[0]=focalLength;
      PixelInfo[1]=pixwidth;
      PixelInfo[2]=pixheight;
      if (topRight[0]>0)ImagePlane[2]=topRight[0];
      else ImagePlane[2]=-topRight[0];
      if (topRight[1]>0)ImagePlane[3]=topRight[1];
      else ImagePlane[3]=-topRight[1];
      ImagePlane[0]=-ImagePlane[2];
      ImagePlane[1]=-ImagePlane[3];
    }
    void deletePartialFeatures(double distance,LongList * llist=0)
    {
      if (!llist)llist=&NearKeys;
      MapFeatureList temp;
      FeatureFcn::listFeatures(temp,*Bank,*llist);
      for (MapFeatureList *mlist=&temp; mlist->Next; mlist=mlist->Next)
	if (mlist->Element->partial(distance))
	  delete mlist->Element;
    }

    /**
     * This finds the features that match this measurement roughly.
     * Call startMeasument before this.
     * @param results a list of posed features form search list that match m.
     * @param m a Match that has been its Measurement, Thresholds and Metric set.
     * @param searchList a list of visable features of the right type.
     */
    void matchFeatures(PosedFeatureList & results, Match &m,
		       PosedFeatureList * searchList);
    /*
     * Returns with m.MatchFeature set to a new PosedFeature for the
     * closest match and m.MatchDistance set to the distance.
     * @param mat a Matchthat has been 'started' and had Metric set.
     * @param searchList a list of visable features of the right type.
     */
    void matchFeature(Match & mat,PosedFeatureList * searchlist);
    void simpleMatch(Match *matches,int n, Cure::Pose3D &psensor,
		     double threshold=0);
    
    void match(Match & mat,PosedFeatureList *searchlist)
    {
      PosedFeatureList matches;
      matchFeatures(matches,mat,searchlist);
      matchFeature(mat,&matches);
    }
    virtual void getPermutation(Cure::Matrix &,
				PosedFeature *pf[2],
				unsigned short){}


    /**
     * This deletes mf0=pf[0]'s MapFeature if typ&1.
     *
     * The subclass might merge all info from mf0 to pf[1]'s MapFeature.  
     * If merging a point it replaces mf0's point and deletes 
     * the point that is not not left attached to 
     * any feature.  It then change mf0->Index to the 
     * index of the common point.  It also sets of PSharedIndex for both 
     * MapFeatures so that  Points[PSharedIndex[k]] is the point that 
     * has the coordinates equalt to P space dim k: 
     * P(k,0)=Points[PSharedIndex[k]].X(k-n,0) where n=2*PSharedIndex[k].  
     * 
     * So for merging Points it will set for both MapFeatures:  
     *
     * NumberSharedPoints, 
     * PSharedIndex.
     *
     * For mf it will set 
     *   
     * Index,  
     * Points[]
     * And delete the old point.
     *
     * For all merges it makes the merged coordinates equal and combines 
     * the relavant info. 
     *
     * @param pf  These point to the PosedFeatures to be merged or 
     *            that have points to be merged.                  
     * @param typ the return value from findMergeFeatures.
     */
      virtual int  merge(PosedFeature *pf[2],unsigned short typ=1){
      MapFeature *mf=getMapFeature(pf[0]->FeatureKey);
      if (mf)
	if (typ&1){
	  delete mf;
	}
      return 0;
    }

    /** 
     * Typically overloaded by subclasses to do something more than
     * just adding the distance
     */
    virtual int addDistance(Match &mat,double distance){
      mat.PathDistance=distance;
      return 0;
    }    
    /** 
     * Write the TrackKeys to a file.
     * @param fs the file to write to
     */
    int writeTrackKeys(std::fstream &fs );
    /** 
     * Read the TrackKeys from a file.
     * @param fs the file to read from
     */
    void readTrackKeys(std::fstream &fs );
     /**
      * Write the feature for a file
     * @param key the key of the feature to be writen to fs.
     * @param fs the file to write to. It should be opened with std::ios::out
     * @return 0 if ok 1 if fs not good an -1 if key not good.
     */
    virtual int writeFeature(long key,std::fstream &fs );
    /**
     * This makes a MapFeature from the infomation in a text file.
     * @param fs the file to read from. It should be opened with std::ios::in
     * @param version The version shoud have been read in the section header
     *                in the file before the list of Features. 
     * @param type If the next feature is not understood by the helper
     *             (for ex. reading a wall with a pointFeatureHelper) 
     *             the Feature type of the next feature is returned here.    
     * @param readtype if the feature type has already been parsed the 
     *        MapHelper will skip reading it by setting readtype=false.
     * @return 0 if fail else the MapFeature created
     */
    virtual MapFeature * readFeature(std::fstream &fs,int version,
				     int *type,bool readtype=true);
   
    /**
     * Implemented by subclass to find two features froma list that
     * can be merged into on.
     * 
     * @param pfp the features to be merged are returned here,
     * @return meaning is specified in subclass.
     */
    virtual unsigned short findMergeFeatures(PosedFeature * pfp[2],
					     PosedFeatureList *)
    {return 0;}

    /**
     * constraint is a*(pl[0].P,pl[1].P)=b 
     * dX=-K(aP-b)=-Ke;
     * e=a*(pl[0].P,pl[1].P)-b 
     * C'=(I-Ka)C
     * K=Ca^T(aCa^T)^-1
     * 
     * @param typ The Merge type gives information on how to set up a and b
     */
    virtual int getMergeContstraint(Cure::Matrix &a, Cure::Matrix &e,
				     Cure::Matrix &cov,
				     double distance,
				     unsigned short typ, 	
				     PosedFeature * pl[2]){return 0;}

    /**
     * This is implemented by subclasses to force the MapFeatures with
     * points to be merged to both prepare to extend the P space for
     * those points to be 'full rank'.  One must then call extend to
     * actually cause any change to the P space.  
     *
     * @param pf These point to the PosedFeatures that have points to
     * be merged.
     * @param typ the return value from findMergeFeatures.
     */
    virtual void forceExtend(PosedFeature *pf[2], unsigned short typ){}
  };
}
#endif
