// = AUTHOR(S)
//    John Folkesson
//    
//    June 1, 2006
//
//    Copyright (c) 2006 John Folkesson
//    
#ifndef FEATUREDESCRIPTORS_H
#define FEATUREDESCRIPTORS_H
#include "Matrix.hh"
#include "ShortMatrix.hh"
#include "CureDebug.hh"
#include "Vector3D.hh"
#include "MapPoseTree.hh"
#include "MapFeature.hh"
#include <fstream>
#include <string>

#define CURE_FEATURE_NOMINAL_RANGE 100
namespace Cure{
  class RelativePointFeature;
  class FeatureDescriptors;
 /**
   * This class provides detailed information and methods on a point
   * feature measurement.  A subclass might contain SIFT descriptors,
   * color histograms, templets.  
   * 
   *
   */
  class FeatureDescriptor
  {
  protected:
  public:
    unsigned short m_Type;
    /**
     *  Unit vector in the direction of the feature.
     */
    Cure::Vector3D m_Bearing;
    /**
     * The range to the feature.
     */  
    double m_Range;
    
    FeatureDescriptors *m_Descriptors;

    long m_FeatureKey;

    FeatureDescriptor(){
      m_Descriptors=0;
      m_Type=0;
      m_Range=CURE_FEATURE_NOMINAL_RANGE;
      m_FeatureKey=-1;
    }
    FeatureDescriptor(FeatureDescriptors *des,
		      double range=0);
    /**
     * @param b the bearing vector to the point.
     * @param des this should be connected to som MapPoseTree 
     * @param range an estimated range to the point
     */
    FeatureDescriptor(double b[3], FeatureDescriptors *des,
		      double range=0);
    
    virtual ~FeatureDescriptor();
    /**
     * @return a pointer to the sensor pose that this is relative to.
     */
    Transformation3D *sensorPose();
    /**
     * @return the CovType of the sensor pose.
     */
    unsigned short getCovType();
    /**
     * @return a pointer to the Relatic
     */
    RelativePointFeature * getRelativePointFeature();
    /**
     * Gets the bearing to the feautre.
     * @param bearing the vector pointing to the feature from m_SensorPose is
     * returned here.
     */
    void getBearing(Cure::Vector3D &bearing){
      bearing=m_Bearing;
    }
    /**
     * This weights the implied range from triangulation of two
     * descriptors.  The weight is tthe sin^2
     * angle between the bearings.
     * 
     * @param wr the range times the weight returned here. 
     * @param f the feature descriptor to compute the range in this.
     * descriptors sensor frame from
     * @return the weight or -1 if fails to find sensorpose
     */
    double getWeightedRange(double &wr,
			    FeatureDescriptor &f);

    /**
     * This weights the implied angle theta with the z-axis.  The
     * weight is the inverse range times the sin^2 angle between the
     * bearings.
     * 
     * @param wt the theta times the weight returned here. 
     * @param f the feature descriptor to compute phi in this
     * descriptors sensor frame from.
     * @return the weight or -1 if fails to find sensorpose
     */
    double getWeightedTheta(double &wt,
			    FeatureDescriptor &f);
    /**
     * This weights the implied angle phi around the x axis.  The
     * weight is the inverse range times the sin^2 angle between the
     * bearings.
     * 
     * @param wt the phi times the weight returned here. 
     * @param f the feature descriptor to compute phi in this
     * descriptors sensor frame from.
     * @return the weight or -1 if fails to find sensorpose
     */
    double getWeightedPhi(double &wt,
			  FeatureDescriptor &f);


    /**
     *  Find out if the pose is the one that this descriptor was seen from.
     *  @param p the pointer to the Transformation3D object that this->m_Sensor
     *  @return true if this descritor was not relative to p
     */
    bool notFrom(Transformation3D *p){
      return (sensorPose()!=p); 
    }
    virtual MapFeature * getMapFeature(MapBank *b){
      return Cure::getMapFeature(b,m_FeatureKey);
    } 
    /**
     * gets the integer pixel location in the image for this
     * @param pix the pixels are returned here 
     */
    virtual void getPixels(short pix[2]){}
    /**
     * gets the double pixel location in the image for this.  So 
     * fractions of a pixel are possible (ie Lucas-Kanade)
     * @param pos the pixels are returned here 
     *
     */
    virtual void getImagePosition(float pos[2]){}
    /**
     * This returns a distance from this descriptor to f.
     * @param f the descriptor to compare.
     * @return a distance from this descriptor to f.
     */
    virtual float match(FeatureDescriptor &f){return 3.0;}
    /**
     * This returns an energy of the epipolar constraint between two
     * descriptors.  It is calcualted based on the angle between the bearing
     * of the given descriptor and the bearing to the closest point on
     * the bearing ray this descriptor.  This only uses the m_Range estimate 
     * if there is not sufficient triangulation to define the closeset point.
     *
     * @param f the descriptor to compare to
     * @param bearingRadsSqErr the variance of the bearing measurments
     * @return half the angle squared / variance 0, or -1 if failed to find
     * sensorpose. 
     */
    double checkEpipolar(FeatureDescriptor &f,
			 double bearingRadsSqErr=.0004);
    /**
     * This returns an energy of the epipolar constraint between two
     * descriptors.  It is calcualted based on the angle between the
     * bearing of the given descriptor and the bearing to the closest
     * point on the bearing ray this descriptor.  It also calculates
     * the gradient and hessian wrt sensor pose covtype dimensions for
     * (f, this).  So g has rows that start with the gradient wrt
     * f.sensorPose() CovCoordinates then rows for this->sensorPose()
     * CovCoordinates.  One should accumulate these derivatives for
     * all constraints. One can then calculate a change in the pose by
     * dx = -h^-1 g.  This only uses the m_Range estimate if there is
     * not sufficient triangulation to define the closeset point.
     *
     *
     * @param f the descriptor to compare to
     * @param bearingRadsSqErr the variance of the bearing measurments
     * @param g the gradient is returned here.
     * @param h the Hessian is returned here.
     *
     * @return half the angle squared / variance or -1 if failed to find
     * sensorpose. 
     */
    double checkEpipolar(FeatureDescriptor &f, Matrix &g, Matrix &h,
			 double bearingRadsSqErr=.0004);
        
  };
 
  /**
   * This class should be associated with a PoseTree object that has
   * the m_SensorPose for the FeatureDescriptors object The
   * descriptors should be given the pose when created and deleted
   * when the pose is.
   */
  class FeatureDescriptors
  {
  public:
    /**
     * The list of all descriptors on this pose.  Some of the pointers
     * might be null.
     */
    FeatureDescriptor **m_Desc;
    /** The legth of m_Desc.*/
    unsigned short m_NumberDesc;
    
    Transformation3D *m_SensorPose;
    unsigned short m_CovType;
    MapPoseTree *m_Tree;
    


  public:
    FeatureDescriptors(){
      m_Desc=0;
      m_SensorPose=0;
      m_CovType=0;
      m_NumberDesc=0;
      m_Tree=0;
    }
    /**
     * @param n the number of Descriptors to be added
     * @param tree the PoseTree
     * @param branch the branch of th tre to the leaf of sensorpose
     */
    FeatureDescriptors(MapPoseTree *tree, short branch, int n=0){
      m_Tree=tree;
      if (m_Tree){
	m_SensorPose=m_Tree->getLeafPose(branch);
	m_CovType=m_Tree->getLeafPoseType(branch);
      }
      if (n==0){
	m_Desc=0;
	m_NumberDesc=0;
      } else {
	m_Desc=new FeatureDescriptor* [n];
	m_NumberDesc=n;
	memset(m_Desc,0,n*sizeof(FeatureDescriptor*));
      }
    }
    virtual ~FeatureDescriptors(){
      clear();
    }
    MapBank *bank(){
      if (m_Tree)return m_Tree->Bank;
      else return 0;
    }
    /**
     * This should only be used if there is no m_Tree
     */
    void setSensorPose(Transformation3D *sensorpose,unsigned short covtype){
      m_SensorPose=sensorpose;
      m_CovType=covtype;
    }
    void setMapPoseTree(MapPoseTree *tree, short branch){
      m_Tree=tree;
      if (m_Tree){
	m_SensorPose=m_Tree->getLeafPose(branch);
	m_CovType=m_Tree->getLeafPoseType(branch);
      }
    }
    /**
     * Add a descriptor to the first null position on m_Desc,
     * extending as needed.  This also sets the m_SensorPose of f.
     * @param f the descriptor to add
     *
     */
    void add (FeatureDescriptor *f){
      f->m_Descriptors=this;
      if (m_Desc){
	int again=m_NumberDesc;
	for (int i=0; i< m_NumberDesc;i++)
	  if (m_Desc[i]==f)return;
	  else if (m_Desc[i]==0)
	    if (again==m_NumberDesc)again=i;
	if (again!=m_NumberDesc){
	  m_Desc[again]=f;
	  return;
	}
	FeatureDescriptor **s=m_Desc;
	m_Desc=new FeatureDescriptor* [m_NumberDesc+1];
	memcpy(m_Desc,s,m_NumberDesc*sizeof(FeatureDescriptor*));
      }else m_Desc=new FeatureDescriptor* [1];
      m_Desc[m_NumberDesc]=f;
      m_NumberDesc++;
    }
    
    /**
    * This deletes all the FeatureDescriptors on the list.
    */
    void clear();
    /**
     * Get the ith descriptor.
     * @param i the index of the descriptor on m_Desc.
     * @return the ith Featuredescriptor stored on this object.
     */
    FeatureDescriptor * operator() (const unsigned short i){
      if (i<m_NumberDesc)
	return m_Desc[i];
      return 0;
    }
    /**
     *  Looks for f and replaces it with null.
     * @param f the descriptor to remove.
     */
    void remove(FeatureDescriptor *f){
      f->m_Descriptors=0;
      for (long i=0;i<m_NumberDesc;i++)
	if ((m_Desc[i])==f)m_Desc[i]=0;
    }
    /**
     * Finds a feature key on the list of descriptors.
     * @param key the key to look for.
     * @return 0 if not found else the first matching FeatureDescriptor.
     */
    FeatureDescriptor * findFeature(long key){
      for (int i=0;i<m_NumberDesc;i++)
	if (m_Desc[i])if (m_Desc[i]->m_FeatureKey==key)return m_Desc[i];
      return 0;
    }
    /**
     * gets the MapFeature matched to the index Descriptor 
     * @param index the index into m_Desc.
     * @return the pointer to the MapFeature;
     */
    MapFeature * getMapFeature(unsigned short index){
      MapBank *b=bank();
      if (b)
	if (index<m_NumberDesc)
	  if (m_Desc[index])
	    return m_Desc[index]->getMapFeature(b);
      return 0;
    }
    /**
     *  Looks for f aan returen its index into m_Desc.l.
     * @param f the descriptor to find.
     * @return the index into m_Descriptor of -1 if not found.
     */
    long getIndex(FeatureDescriptor *f){
      for (long i=0;i<m_NumberDesc;i++)
	if ((m_Desc[i])==f)return i;
      return -1;
    }
    /**
     * gets the bearing to the ith feature
     *
     * @param i the feature index
     * @param bearing the bearing is retured here.
     * @return 1 if fail else 0.
     */
    int getBearing(unsigned short i,Cure::Vector3D &bearing){
      if (i<m_NumberDesc)
	if (m_Desc[i]){
	  m_Desc[i]->getBearing(bearing);
	  return 0;
	}
      return 1;
    }

    
    int match(FeatureDescriptors &fd, 
	      Cure::Matrix &Errors,
	      Cure::ShortMatrix &Matches, 
	      float limit=2.0);
    bool  match(FeatureDescriptor &f, unsigned short &index,float &matchError)
    {
      matchError=10;
      index=0;
      for (int i=0; i<m_NumberDesc; i++)
   	{
	  float d=m_Desc[i]->match(f);
	  if (d<matchError)
	    {
	      matchError=d;
	      index=i;
	    }
	}
      if (matchError<10){
	return true;
      }
      return false;
    }
    /**
     * This finds a number of candidate matches to a descriptor from the 
     * set of m_Desc.
     * @param fd the FeatureDescriptro to match to
     * @param mats call with Rows set to max number of matches to find.
     * The index of the descriptors that match are returned here in order of 
     * best match.
     * @param errorlimit applied to the return of match from the
     * FeatureDescriptor
     */
   void getSomeMatches(FeatureDescriptor &fd,ShortMatrix &mats,
		       float errorlimit)
    {
      int index=0;
      Transformation3D *pose=fd.sensorPose();
      if (!pose)return;
      Matrix merr(mats.Rows,1);
      for (int i=0; i<m_NumberDesc; i++)
   	{
	  float d=m_Desc[i]->match(fd);
	  if (d<errorlimit){
	    int ind=i;
	    if (index<mats.Rows){
	      for (int j=0;j<index; j++){
		if (d<merr(j,0)){
		  double temp=merr(j,0);
		  merr(j,0)=d;
		  d=temp;
		  int tempi=mats(j,0);
		  mats(j,0)=ind;
		  ind=tempi;
		}
	      }
	      mats(index,0)=ind;
	      merr(index,0)=d;
	      index++;
	    }else{
	      if (d<merr(index-1,0)){
		for (int j=0;j<index; j++){
		  if (d<merr(j,0)){
		    double temp=merr(j,0);
		    merr(j,0)=d;
		    d=temp;
		    int tempi=mats(j,0);
		    mats(j,0)=ind;
		    ind=tempi;
		  }
		}
	      }
	    }
	  }
	}
      mats.Rows=index;
    }
    /**
     * This finds the first good matchs from a FeatureDescriptor.
     *
     * @param des the FeatureDescriptor that you are trying to match to
     * the descriptors on this.
     *
     *
     * @param errorlimit applied to the return from
     * FeatureDescriptors::match(...)
     * @param radsqerr eatimated variance in the bearing angle
     * @param coarseEnergyLimit the threshold on the epipolarconstraint.
     * @param sample the max number of matching Descriptors to try Epipolar on.
     * @return the matching descriptor
     * 
     */
    FeatureDescriptor * getEpipolarMatch(FeatureDescriptor *des,
					 float errorlimit,
					 double radsqerr,
					 double coarseEnergyLimit,
					 int sample=5)
    {
      ShortMatrix mats(sample,1);
      getSomeMatches(*des,mats,errorlimit);
      for (int j=0;j<mats.Rows; j++)
	{
	  if (des->checkEpipolar(*m_Desc[mats(j,0)],radsqerr)
				 <coarseEnergyLimit){
	    return  (*this)(mats(j,0));
	  }
	}
      return 0;
    }	

    /**
     * This finds a reverse consitent match from a FeatureDescriptor to a
     * list of this FeatureDescriptors.  
     * @param des the FeatureDescriptor that you are trying to match to
     * the descriptors on this.
     * @param index the index into des.
     *
     * @param errorlimit applied to the return from
     * FeatureDescriptors::match(...)
     * @param radsqerr eatimated variance in the bearing angle
     * @param coarseEnergyLimit the threshold on the epipolarconstraint.
     * @param sample the max number of matching Descriptors to try Epipolar on.
     * @return the matching descriptor
     * 
     */
    FeatureDescriptor * match(FeatureDescriptors &des,
			      int index,
			      double errorlimit,	
			      double radsqerr,
			      double coarseEnergyLimit,
			      int sample=5){
      FeatureDescriptor *f=des(index);
      if (f){
	FeatureDescriptor *f2=getEpipolarMatch(f,
					       radsqerr,
					       errorlimit,
					       coarseEnergyLimit,
					       sample);
	if (f2) {
	  FeatureDescriptor *f3=des.getEpipolarMatch(f2,
						     radsqerr,
						     errorlimit,
						     coarseEnergyLimit,
						     sample);
	  if (f3==f)return f2;
	  
	}
      }
      return 0;
    }    
  };
  class FeatureDescriptorList
  {
  public:
    FeatureDescriptor * Element; 
    FeatureDescriptorList *Next;
    FeatureDescriptorList();
    ~FeatureDescriptorList();
    FeatureDescriptor * operator()(unsigned short i){return get(i);}
    void clean();
    unsigned short add(FeatureDescriptor *w);
    unsigned short addUnique(FeatureDescriptor *w);
    bool remove(unsigned short k);
    unsigned short count();
    unsigned short removeDescriptor(FeatureDescriptor *sw);
    bool find(FeatureDescriptor *n); 
    FeatureDescriptor * get(unsigned short n);
    FeatureDescriptor * getRemove(unsigned short n);
  protected:
  
  };
}  
  
#endif
