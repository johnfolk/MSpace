// = AUTHOR(S)
//    John Folkesson
//    
//    Aug 1, 2007
//
//    Copyright (c) 2007 John Folkesson
//    
#ifndef CURE_MAPPOSETREE_HH
#define CURE_MAPPOSETREE_HH


#include "MapBank.hh"
#include "PoseTree.hh"
#include "Matrix.hh"

namespace Cure{
  class FeatureDescriptors;
  class FeatureDescriptor;
  /**
   * A Cure::PoseTree is a tree of tranformation links with machinary
   * to comput leaf poses and the jacobians wrt link coordinates.
   *
   * 
   *
   */
  class MapPoseTree: public MapObject,public PoseTree
  {
  public:
    /**
     * Holds the Feature Descriptor Objects for esch leaf in the tree.
     */
    Cure::FeatureDescriptors **m_Descriptors;
  protected:
  private: 
    unsigned short m_AllocatedDescriptors;
    /*
     * This is an int array of length dim initialized to 0.
     * It is meant to hold the index into the 
     * Covariance matrix or some other structures.
     */

    long  *Index;
    unsigned short m_SizeIndex;
  public:    
    
    MapPoseTree(unsigned short covtype=0):      
      MapObject(0),
      PoseTree(covtype),
      m_Descriptors(0),
      m_AllocatedDescriptors(0),
      Index(0),
      m_SizeIndex(0)
    {
      reallocate();
    }

    MapPoseTree(MapBank *b):      
      MapObject(b),
      PoseTree(0),
      m_Descriptors(0),
      m_AllocatedDescriptors(0),
      Index(0),
      m_SizeIndex(0)
    {
      reallocate();
    }



    MapPoseTree(MapPoseTree &templet):
      MapObject(templet.Bank),
      PoseTree(),
      m_Descriptors(0),
      m_AllocatedDescriptors(0),
      Index(0),
      m_SizeIndex(0)
    {
      (*this)=templet;
    }
    MapPoseTree(PoseTree &templet):
      MapObject(0),
      PoseTree(templet),
      m_Descriptors(0),
      m_AllocatedDescriptors(0),
      Index(0),
      m_SizeIndex(0)
    {
    }

    MapPoseTree(unsigned short covtype,
		double x[6]):
      MapObject(0),
      PoseTree(covtype,x), 
      m_Descriptors(0),
      m_AllocatedDescriptors(0),
      Index(0),
      m_SizeIndex(0)
    {
    }
    MapPoseTree(short covtype, double rootpose[6],
		ShortMatrix &branches, 
		ShortMatrix &posetypes,
		ShortMatrix &linktypes, 
		Matrix &initlinkposes):
      MapObject(0),
      PoseTree(covtype,rootpose,branches,posetypes,linktypes,initlinkposes),
      m_Descriptors(0),
      m_AllocatedDescriptors(0),
      Index(0),
      m_SizeIndex(0)
    {
    }
    /**
     * Deleting a MapPoseTree will delete its FeatureDesccritor objects
     * 
     */
    virtual ~MapPoseTree();  
    /**
     * The copy operator
     */
    void operator=(MapPoseTree &tree){
      clean();
      m_Time=tree.m_Time;
      m_Branches=tree.m_Branches;
      m_RootPose=tree.m_RootPose;
      m_NumberOfLinks=tree.m_NumberOfLinks;
      m_Links=new Link*[m_NumberOfLinks];
      m_CovType=tree.m_CovType;
      for (int i=0;i<m_NumberOfLinks;i++){
	m_Links[i]=new Link(*tree.m_Links[i]);
	m_Links[i]->m_FromPose=&m_RootPose;
      }
      reallocate();
      int top=m_SizeIndex;
      if (top>tree.m_SizeIndex)top=tree.m_SizeIndex;
      for (int i=0;i<top;i++)Index[i]=tree.Index[i];
    }
    /**
     * Sets the Map bank for the tree.
     * @param b the mapBank
     */
    void setBank(MapBank *b){
      Bank=b;
      if (Bank)Bank->add(this);
    }
    unsigned short getAllocatedDescriptors(){return m_AllocatedDescriptors;}
    long index(unsigned short n){
      if (n<m_SizeIndex)
	return Index[n];
      return 0;
    }
    void setIndex(unsigned short n, long value){
      if (n<m_SizeIndex)
	Index[n]=value;
    }
    virtual void reallocate();
    virtual void clean();

    virtual void init();
    /**
     * Get a pointer to a subclass
     * @param type the type of the subclass you want.
     * @return a pointer that can be cast to the subclass
     * cooresponding to the type or if that is not possible NULL.
     *
     */
    MapObject * getFrameType(long type){
      if (type==POSETREE_TYPE)return this;
      return 0;
    } 
    MapObject * getMapObjectType(short type){
      if (type==FRAME_TYPE)return this;
      return 0;
    }
    /**
     * This gets the Descriptors.
     * @param branch the branch to add it to
     * @return pointer to the branch's FeatureDescriptors object.
     */
    FeatureDescriptors * getDescriptors(unsigned short branch){
      if (branch>m_AllocatedDescriptors){
	reallocate();
	if (branch>m_AllocatedDescriptors)
	  return 0;
      }
      return m_Descriptors[branch];
    }

    /**
     * Adds a FeatureDescriptor  (This descriptor is deleted if it is 
     * still on the list when the MapPoseTree is deleted.
     * @param d the descriptor to add.
     * @param branch the branch to add it to
    */
    void addDescriptor(FeatureDescriptor *d, unsigned short branch);
    
    /**
     * Removes a FeatureDescriptor
     * @param d the descriptor to remove.
     * @param branch the branch to add it to
     * 
     */
    void removeDescriptor(FeatureDescriptor *d, unsigned short branch);
    

    /**
     * used to set all the index elements 
     * @param ind and array of length dim() that holds the values to use.
     */
    void setIndex(long *ind){
      int top=m_SizeIndex;
      for (int i=0; i<top;i++)Index[i]=ind[i]; 
    }
    /**
     *
     * Sets the index for the covDim dim for example the rows of a
     * state vector.
     * 
     * Assumes that the indecies are consecutive
     * @param ind The vaule for the first index.  
     */
    void setIndex(long ind){
        int top=m_SizeIndex;
	
	for (int i=0; i<top;i++)Index[i]=ind+i; 
    }
    
    /** 
     *get the timestamp of this posetree
     * @return the time 
     */
    Timestamp & time(){return m_Time;}
    
    /**
     * saves the data of this tree in a GenericData object.
     *
     */
    void get(GenericData &gd);
    /**
     * Restors the data of this tree from a GenericData object.
     * This will clear all Descriptors from this object.
     * @return 1 if fail else 0
     */
    int set(GenericData &gd);
    
    /**
     * Adds dx to the covdim part of the coordinates for i = 0; i<covDim().  
     * It then sets dp=0.
     * @param dx Column matrix of the change in cov dof dimensions.
     * @return number of dof dimensions if ok else -1. 
     */
    int updateX(Cure::Matrix & dx){
      Matrix x;
      getX(x);
      if (dx.Rows==x.Rows){
	x+=dx;    
	dx=0;
	setX(x);
	return covDim();
      }
      return -1;
    }

    /**
     * This is a bit smarter than updateX, it expects dxp to be a column vector
     * from the start of the array that this->Index refers to.  So this
     * will index correctly even for non consecutive indecies.
     */
    void updateIndexedX(Matrix &dx){
      Cure::Matrix d(covDim(),1);
      for (int i=0; i<dx.Rows; i++)
	{
	  d(i,0)=dx(Index[i],0);
	  dx(Index[i],0)=0;	    
	}
      updateX(d);	
    }
    /**
     * Add dx to the full coordinates x.
     *
     */
    void updateState(Cure::Matrix & dx){
      Matrix x;
      getState(x);
      if (dx.Rows==x.Rows){
	x+=dx;
	dx=0;
	setState(x);
      }
    }

    /**
     * This returns the cov coordiantes of the tree transformed to
     * some new frame along with the jacobian of the cov coodinates
     * wrt the cov coordinates of the transform to the new frame and
     * the old coc coordinates.
     *
     * @param xnew the new CovCoodinates are returned here (only the
     * base link is transformed).
     *
     * @param jac the jacobian wrt the (newnullframe cov coodinates,
     * this cov coodinates) is returned here.  
     *
     * @param newnullframe the tranform from the origin of the current
     * frame to the origin of the new frame.
     *
     * @param covtyp the CovType to use for newnullframe.
     */
    void transform(Matrix &xnew,Matrix &jac,Transformation3D &newnullframe, 
		   unsigned short covtyp){
      Transformation3D temp;
      int typ=(covtyp<<6)+m_CovType;
      temp.dominusAplusB(newnullframe,m_RootPose,jac,m_CovType,typ,1);
      getX(xnew);
      temp.getCovCoordinates(&xnew(0,0),m_CovType);
    }
    virtual MapObject * getObjectType(short type){
      if (type==POSETREE_TYPE)return this;
      return 0;
    }
    /**
     * Transforms this to a new frame.
     * @param newnullframe the tranform from the origin of the current
     * frame to the origin of the new frame.
     */
    void transform(Transformation3D &newnullframe){
      m_RootPose.leftSubtract(newnullframe);
      CalcState=(CalcState&3);
    }
    short getObjectSubType(){return POSETREE_TYPE;}
    short getObjectType(){return FRAME_TYPE;}
    void print();
    
  protected:
  
  };
  class MapPoseTreeMaker: public MapObjectMaker{
  public:
    MapPoseTreeMaker(MapBank *b):MapObjectMaker(b){
      ObjectType=FRAME_TYPE;
      ObjectSubType=POSETREE_TYPE;
    }
    MapObject * makeMapObject(){
      return new MapPoseTree(Bank);
    }
  };
  static inline MapPoseTree *getMapPoseTree(MapBank *b,long key)
  {
    if(!b)return 0;
    MapObject *mo=b->getMapObject(key);
    if (!mo) return 0;
    mo=mo->getObjectType(FRAME_TYPE);
    if (!mo) return 0;
    return (MapPoseTree *)mo->getFrameType(POSETREE_TYPE);
  }  
}
#endif
