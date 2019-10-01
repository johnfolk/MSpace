// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_WALLHELPER_HH
#define CURE_WALLHELPER_HH

#include "MapWall.hh"
#include "PosedWall.hh"
#include "MapHelper.hh"
#include "Match.hh"
#include "Trig.hh"

namespace Cure {

/**
 * Class that manages Wall features in the map
 * 
 * @author John Folkesson
 */
  class WallHelper: public MapHelper 
  {
  public:
    PosedWall * PosedWallPtr;
    MapWall * MapWallPtr;
    MapWall TempletFeature;
    Cure::Line2D Line;
    Cure::Trig Triger;

    /**
     * Flag to use the endpoint detections Measurement type>1 or not.
     * if Not then no corners will be found either.  
     * @default true	    
     */
    bool UseEndpoints;

    /*
     * This will refuse to use endpoint measurements if they are further then
     * this distance in meters from the current estimated endpoint position.
     * @default .01	    
    */
    double EndPointUpdateThreshold;;

    /*
     * This sets how much overlap in meters between two wall in order to
     * Qualify for merging.  This can be negative, then it is max gap.
     * @default .01
     */
    double    MergeOverLap;
 
    /*
     * Maximum difference in the angles of two wall to merge, in radians.
     * @default .03
    */
    double MergeMaxGammaError;
    
    /*
     * Maximum perpendicular distance of two wall to merge, in meters
     * @default .05
     */
    double MergeMaxRhoError;
    
    /*
     * Max separation between two endpoints from two walls for forming 
     * a corner point out of them units square meters.
     * @default .04
     */
    double CornerSeparation;
    
    /** 
     * The Templet will be given these thresholds which then are used to 
     * initilize the MapFeatures made with this helper. 
    */
    WallHelper(MapBank *b, int countThreshold=25,double tightness=0.3,
	       double endThreshold=.1,
	       double lengthThreshold=.5, double varRhoThreshold=1E3,	      
	       double distanceThreshold=20):
    MapHelper(b){
      TempletFeature.setBank(b);
      TempletFeature.CastPtr=&MapWallPtr;
      TempletFeature.CountThreshold=countThreshold;
      TempletFeature.TightnessValue=tightness;
      TempletFeature.EndThreshold=endThreshold;
      TempletFeature.LengthThreshold=lengthThreshold;
      TempletFeature.VarRhoThreshold=varRhoThreshold;
      TempletFeature.DistanceThreshold=distanceThreshold; 
      SensorType=SCAN_TYPE;
      EndPointUpdateThreshold=.01;
      MergeOverLap=.01;
      CornerSeparation=.04;
      MergeMaxGammaError=.03;
      MergeMaxRhoError=.05;
      UseEndpoints=true;
      FeatureType=MAPWALL_TYPE;
      ObjectSubType=FeatureType;
    }
    ~WallHelper(){}

    virtual void printConfiguration();  
  /**
   * This increases (or decreases) the  threshold for 
   * the merging lines.
   * @param factor one is no change, >1 makes it easier to merge lines,
   *               <1 makies it harder. 
   */  
    void loosenMerging(double factor);


    PosedWall * castPosedWall(PosedFeature *pf){
      PosedWallPtr=0;
      pf->narrow();
      return PosedWallPtr;
    }

    MapWall * castMapWall(MapObject *pf){
      MapWallPtr=0;
      pf->narrow();
      return MapWallPtr;
    }

    MapWall *getMapWall(long key){
      MapObject *mo=Bank->getMapObject(key);
      if (!mo) return 0;
      return castMapWall(mo);
    } 

    MapFeature *getMapFeature(long key)
    {
      return getMapWall(key);
    } 
    virtual MapFeature * makeMapFeature(MapBank *b=0){
      return new MapWall(&TempletFeature,b);
     }
    /**
     * @return A MapHLine based on the TempletFeature. 
     */
    MapWall * makeMapWall(){
      return new MapWall(&TempletFeature,Bank);
    }

    /**
     * 
     * @param sensorpose The pose of the sensro with the y axis ahead, 
     *        x is right and z up.
     * @param x The coordinates of the start and endpoints in 
     * in meter units relative to sensor.
     * @return A MapWall based on the TempletFeature and the params. 
     */
    MapWall * makeMapWall(Cure::Pose3D & sickpose,double x[4]) 
    {
      Line.StartPoint.setXY(x);
      Line.EndPoint.setXY(x+2);
      MapWall *tempwall=new MapWall(&TempletFeature,Bank);
      tempwall->initailizeFromLine(sickpose,Line);
      return tempwall;
    }     

    /**
     * This takes a MapFeature and makes a PosedFeature
     * out of it
     */
    virtual PosedFeature * makePosedFeature(MapFeature* mf){
      return makePosedWall(mf);}
    virtual MapFeature * makeMapFeature(Cure::Pose3D &sensorpose,
					Cure::Measurement &m,
					bool addToVisable=false);

    PosedWall * makePosedWall(MapFeature* mf){
      if(mf->Type!=MAPWALL_TYPE)return 0;
      PosedWall *pw= new PosedWall(&PosedWallPtr,&Triger, castMapWall(mf));
      pw->EndPointUpdateThreshold=EndPointUpdateThreshold;
      return pw;
    }

    virtual bool getC(MapFeature  *mf1,MapFeature  *mf2,
		      Cure::Matrix &c1){
	MapWall *mw1=castMapWall(mf1);
	if (!mw1)return false;
	MapWall *mw2=castMapWall(mf2);
	if (!mw2)return false;
	return mw1->getC(mw2,c1);
      }
    virtual bool testMatch(MapFeature  *mf1,MapFeature  *mf2,
			   double tolerance){
	MapWall *mw1=castMapWall(mf1);
	if (!mw1)return false;
	MapWall *mw2=castMapWall(mf2);
	if (!mw2)return false;
	return mw1->testMatch(mw2,tolerance);
      }
    virtual double judgeMatch(MapFeature  *mf1,MapFeature  *mf2,
			      double tolerance){
	MapWall *mw1=castMapWall(mf1);
	if (!mw1)return false;
	MapWall *mw2=castMapWall(mf2);
	if (!mw2)return false;
	return mw1->judgeMatch(mw2,tolerance);
      }
 
    /**
     * Config the helper with a string of parameters where the first
     * parameter should be the version number
     */
    virtual int config(const std::string &arglist);    

    virtual bool supportsSubconfig(int sc); 

    /**
     * Puts the Distance on the cloud and adds a Metric/Thresholds to the
     * Match, based on RoughSearchRange.
     */
    int addDistance(Match & mat,double distance);
  
    /**
     * This searches a list for two PosedFeatures that can be merged.
     * @param pfp The pfp[0] will point to the feature that should be deleted
     *            while pfp[1] will remain.
     * @param searchlist The list to look for features.
     * @return 0 for no merge found.
     *         0x0001 for merge whole features
     *         0x1ba0 for merge Points[a] of pfp[0] to Points[b] of pfp[1] 
     *            
     */
    unsigned short  findMergeFeatures(PosedFeature * pfp[2],
				      PosedFeatureList * searchList);

    /**
     * This calculates the matricies needed to impose the merge constraint on
     * a EKF without making any 'real' changes to the features.
     * constraint is a*(pl[0].Xp,pl[1].Xp)=b 
     * dX=-K(ap-b)=Ke;
     * C'=(I-Ka)C
     * K=Ca^T(aCa^T)^-1
     * So the measurement is e with jacobian Jep of a and covaraiance cov. 
     * @param typ the return from findMergeFeatures
     * @param pl the returned posedFeatures from findMergeFeatures.
     * @param typ the return from findMergeFeatures.
     * @return 1 if match is within thresholds else 0
    */
    int getMergeContstraint(Cure::Matrix &a, Cure::Matrix &e,
			     Cure::Matrix & cov,
			     double distance,unsigned short typ, 
			     PosedFeature * pl[2]);
    /**
     * This makes the P space, LowInfo and shared points of  
     * pf[1]'s MapFeature consistant with pf[0]'s, then deletes pf[0]'s
     * MapFeature. 
     * @param pf pf[0] points to a PosedFeature with the MapFeature to be 
     *            deleted while pf[1]'s will be made consistant with 
     *            the information of pf[0]'s.
     * @typ the return value from findMergeFeatures.
     */  
    virtual int merge(PosedFeature *pf[2], unsigned short typ=1){
      MapWall *mf0=getMapWall(pf[0]->FeatureKey);
      MapWall *mf1=getMapWall(pf[1]->FeatureKey);
      if (!mf1)return -1;
      if (!mf0)return -1;
      int r=mf1->merge(mf0,typ);  
      if (typ&1){
	delete mf0;
      }
      return r;
    }
    /**
     * This force the MapFeatures with points to be merged  to both prepare
     * to extend the  P space for those points to be 'full rank'.  
     * One must then call extend to actually cause any change to the P space.
     * @param pf  These point to the PosedFeatures that have points to
     *            be merged.                  
     * @param typ the return value from findMergeFeatures.
     */
    virtual void forceExtend(PosedFeature *pf[2], unsigned short typ){
      MapWall *mf0=getMapWall(pf[0]->FeatureKey);
      MapWall *mf1=getMapWall(pf[1]->FeatureKey);
      if (typ&0x1000)
	{
	  if (mf0)mf0->forceExtend((0x00F&(typ>>4)));
	  if (mf1)mf1->forceExtend((0x00F&(typ>>8)));
	}
    }
    /**
     * This calculates the matrix P to do P(sigma)P^T to the covariance matix
     * for the EKF.  The columns and rows of permute need to be indexed into
     * the correct columns for the covariance.  These indecies are in the 
     * MapFeature->Index.
     * @param permute the part of P that is not the identity matrix.
     * @param pf  These point to the PosedFeatures that have points to
     *            be merged.                  
     * @typ the return value from findMergeFeatures.
     */ 
    virtual void getPermutation(Cure::Matrix & permute,
				PosedFeature *pf[2],
				unsigned short typ);

    /**
     *  Depreciated use writeFeature
     * @param key the key of the wall to be writen to fs.
     * @param fs the file to write to. It should be opened with std::ios::out
     * @return 0 if ok 1 if fs not good an -1 if key not good.
     */
    int writeWall(long key,std::fstream &fs ){
      return writeFeature(key,fs);
    }

    /**
     *  Depreciated use readFeature
     * This makes a MapWall from the infomation in a text file.
     * @param fs the file to read from. It should be opened with std::ios::in
     * @param version The version shoud have been read in the section header
     *                in the file before the list of walls. 
     * @param type If the next feature is not a wall, 
     *             the Feature type of the next feature is returned here.    
     * @param readtype if the feature type has already been parsed the 
     *        WallHelper will skip reading it by setting readtype=false.
     * @return 0 if fail else the MapWall created
     */
    MapWall * readWall(std::fstream &fs,int version,
		       int *type,bool readtype=true){
      MapFeature *mf=readFeature(fs,version,type,readtype);
      if (mf)
	return castMapWall(mf);
      return 0;
    }

  protected:
    /** 
     * Reads config string for version 1 type configs.
     */
    int configVer1(const std::string &arglist);
  };
}
#endif
