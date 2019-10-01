// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_MAPFEATURE_H
#define CURE_MAPFEATURE_H


#include "MapBank.hh"
#include "MapPointList.hh"
#include "Matrix.hh"

namespace Cure{
  /**
   * MapFeatures are parameterized by 3 sets if coordinates {x3D's},{x2D's}
   * and {Scalars}.  Each set could contain any number of coordinate vectors.
   * The three sets differ by how the coordinates change when transformed to 
   * a sensor frame.  
   *
   * Furthermore the state vector formed by joining all the coordinates in a 
   * single column vector is refered to as X.
   * This X is then split into two parts P and the low info part.  The P part 
   * is defined by a non linear mapping about which we can 
   * know the linearization:
   * 
   * dP = B*dX, 
   *
   * And the projection of a small change in P to a small change in X is:
   *
   *  dX=Bdual*dP
   * 
   * Where B is the projection matrix, (P dimension x Full dimension), and 
   *
   * B*Bdual=I the p-dim identity matrix 
   * 
   * It is the P part that has enough information gathered to 
   * do linearization of the measurements properly, and thus EKF updates.
   *
   * The low part is the remaining dof.  These dimensions will have
   * some information but not enough.  They are set with addInfo. 
   * These updates can even use non-gaussian
   * info, such as point clouds... This kind of info is passed in a Matrix 'w'.
   *
   * The change to the low info coordinates will change the definition of the 
   * P coordinates. 
   *
   * So this all works like this:
   * First one creates the feature.
   * When one has some measurement info, one  
   * calculates the measurement info and then calls:
   * 
   *   addInfo(...)
   *  
   * That will add the measurement info and update the Low coordinates if 
   * possible.
   * 
   * After that one can try:
   *  virtual int extend(...)
   *
   * That will check if there is enough info in some dimensions to start 
   * gaussian updates. 
   * 
   * If not (return 0) continue to collect measurments and add the info.
   * 
   * When extend returns>0, it returns the matricies needed to do a pseudo
   * measurement update of the EKF.  
   * That will give you a psuedo mesurement which you can use at the 
   * current pose to initailize your gaussian filter.
   *
   *
   * After this initialization  proceed thus:
   * 1. Do matching and so on using a PosedFeature.
   * 2. calculate the info in the measurement using a PosedFeature.
   * 3. addInfo(...)
   * 4. if(extend(...))
   *      a. extend your gaussian algorithm with b_n invcov...
   * 5. calculate the gaussian update, dp, using your Gaussian Filter (and 
   *    the PosedFeature to get the jacobians).
   * 6. call  updateP(dp);
   * 7. next measurement.
   *
   *
   * 
   *
   */
  class MapFeature: public MapObject
  {
  public:
    /**
     * Here is where the 2Dand 3D points are stored, ie the state vector.
     * This will normally be set to an Array of pointers[NumberPoints]
     */
    MapPoint **Points;
    /**
     * The scalar array could contain information such as radius, hieght, 
     * width...
     */
    double *Scalars;
    /**
     * Each subclass has a unique Type.
     */
    short Type;
    /**
     * FullDim=3*Number3D+2*Number2D+NumberScalars
     * Xf has this length.
     */
    short FullDim;
    /** Length of the Points array*/
    short NumberPoints;
    /** The number of the Points that are 3D points*/
    short Number3D;
    /** The number of the Points that are 2D points*/
    short Number2D;
    /**The number of scalars*/
    short NumberScalars;
    /** Some points can be shared betweent two features */
    short NumberSharedPoints;
    /**
     * Indicates if the indices are consecutive, i.e. if the
     * P-dimensions are stored consecutive in the state vector /
     * covariance matrix. If this is the case some operations can be
     * performed easier and faster.  So if inot consecutive the
     * feature is scattered aroung in the State vector of the EKF.
     */
    short Consecutive;
    /**
     * An index into Points that correspond for shared points only.
     * Shared points have B=I for their P-Space.
     * So this pointer is null if NumberSharedPoints always was 0.
     * So PSharedIndex[0]=1 means P-Space dimension 0 is from Points[1].  
     */
    short *PSharedIndex;
    /*
     * This stores the distance of the last dense info update.
     */
    double LastDistance;  
    /*
     * This stores the 'bearing' of the last dense info update.
     * This bearing is defined to rotate as a 3D vector under 
     * coordinate transformations.
     */
    double LastBearing[3];
    /*
     * This sets the max 'distance difference' from current distance for the 
     * 'LowInfo' (cloud points).  Further points will be droped as the
     * coorelations get too weak.
     */
    double DistanceThreshold;
    /**
     * Projects from P to Xf. Xf=Bdual*P
     */
    Cure::Matrix Bdual; 
  protected: 
    /*
     * This is an int array of length FullDim initialized to 0.
     * It is meant to hold the index into the 
     * Covariance matrix for the P dimensions.
     */
    long  *Index;
    /**
     * Cumulated Information.  This ends up as inv covariance of xf for 
     * extending dim.
     */
    Cure::Matrix Info;
  public:
    MapFeature(MapBank *b=0);
    /**
     * Deleting a MapFeature will delete its points if the points end up
     * with no features.
     */
    virtual ~MapFeature();
    virtual MapFeature * copy(){return 0;}
    /**
     * Sets the Map bank for the feature.
     * @param b the mapBank
     */
    void setBank(MapBank *b){
      Bank=b;
      if (Bank)Bank->add(this);
    }
    /**
     * Get a pointer to a subclass
     * @param type the type of the subclass you want.
     * @return a pointer that can be cast to the subclass
     * cooresponding to the type or if that is not possible NULL.
     *
     */
    MapObject * getFeatureType(short type){
      if (type==Type)return this;
      return 0;
    }
    MapObject * getMapObjectType(short type){
      if (type==FEATURE_TYPE)return this;
      return 0;
    }
    short  getObjectSubType(){return Type;}
    short getObjectType(){return FEATURE_TYPE;}
    /**
     * This casts the object to a MapPoint and then removes it.
     */
    virtual int removeObject(MapObject *o);
    /**
     * Adds a point to the Points array.
     */
    virtual int addPoint(MapPoint *f,unsigned int which);
   
    /**
     * Removes the feature from f and all pointers to f are set to
     * 0.
     */
    virtual int removePoint(MapPoint *f);
    /**
     * This is better than removePoint(MapPoint *f) in that
     * It will only remove one point (in the odd situation
     * that Points[n]==f for more than one n).
     * It only removes the feature from the point if the 
     * point is only in one place in Points[n].
     *
     * @param which the index intp Points to remove.
     * @return the point removed (or 0 if none)
     */
    virtual MapPoint * removePoint(unsigned short which);
    /**
     * This is used to match the feature and to form near lists.
     * @return 1 if inside else 0
     */
    virtual int inside(double bottomleft[2], double topright[2]){return 0;}
 
    /**
     * This is used to tell if this was a one/few time observation and
     * should be deleted.  
     * 
     * @param distance if the last info was added too far from
     * (current) distance we say it hasn't been seen for a while.
     *
     * @return 1 if it hasn't been initialized and hasn't been seen
     * for a while, else 0.
     */
    virtual int partial(double distance){
      if (Bdual.Columns>0)return 0;
      double mindistance=distance-DistanceThreshold/2;
      if (LastDistance>(mindistance))return 0;
      return 1;
    }
    /**
     * This is used to tell if this is still collecting dense data.
     * 
     * @param distance if the last info was added too far from
     * (current) distance we say it hasn't been seen for a while.
     * 
     * @return false if it hasn't been seen for a while or is at full
     * P-dim, else true.
     */
    
    bool collectingData(double distance){
      if (dimCheck())return false;
      double mindistance=distance-DistanceThreshold/2;
      if (LastDistance>(mindistance))return true;
      return false;
    }
   
    /**
     * This is an aid to tracking by, for example, testing the viewing angle.. 
     * @param testValue meaning depends on Feature.
     * return 0 if can be tracked reliably, else NOT_VISABLE
     */
    virtual int testVisable(double *){return 0;}
    /**
     * This can be used to Cast the subclass.
     */
    virtual void narrow(){}  

    unsigned short pDim(){return Bdual.Columns;}
    /**
     * This is needed by getMapFeature(MapBank *b,long key)
     */
    MapFeature * narrowFeature(){return this;}

    /**
     * Get the array that holds the indices into the Covariance matrix
     * for the P dimensions of this feature
     */
    long* getIndex(){
      return Index;
    }

    /**
     * Checks if the indices are consecutive, i.e. if the P-dimensions
     * are stored consecutive in the state vector / covariance
     * matrix. If this is the case some operations can be performed
     * easier and faster.
     */
    void setConsecutive()
    {
      Consecutive=1;
      for (short i=1; i<Bdual.Columns;i++)
	if (Index[i]!=(Index[i-1]+1))Consecutive=0;
    }

    /*
     * Sets the index for the P dim for example the rows of a state vector.
     * @param ind An array of length Bdual.Columns.  
     */
    void setIndex(long *ind){
      for (short i=0; i<Bdual.Columns;i++)Index[i]=ind[i]; 
      setConsecutive();  
    }

    /**
     *
     * Sets the index for the P dim for example the rows of a state vector.
     * Assumes that the indecies are consecutive
     * @param ind The vaule for the first index.  
     */
    void setIndex(short ind){
      for (short i=0; i<Bdual.Columns;i++)Index[i]=ind+i; 
      Consecutive=1;
    }

    /**
     * Sets Index[which]=ind.
     * @param ind the new value to set Index[which] to.
     * @param which the idex into Index
     * 
     */
    void setIndex(long ind,short  which){
      Index[which]=ind; 
      setConsecutive();  
    }

    /**
     * Changed one of the indices into the covariance matrix. Will
     * look for right index and change it.
     */
    void changeIndex(long fromindex, long toindex){
      for (int i=0; i<Bdual.Columns;i++)
	if (fromindex==Index[i])Index[i]=toindex;
      setConsecutive();  
    }

    /*
     * Finds the best measurement type for this type of measurement
     * data and the current infomation collected on the map feature.
     * Measurement type 0 is no possible linearized measurement.  This
     * is for uninitailized features.
     *
     * @param type the input is a 'suggested type' 
     * @return the best MeausurementType possible given the suggestion and
     *   the P-dim.   
     */
    virtual unsigned short getMeasurementType(unsigned short type){return type;}
    /**
     *  This returns a projection onto the P space for this type of
     *  measurment.  So that this measurment should move dP in a way
     *  that project*dP=dP.
     *  
     * @param project a projection matrix on the P space into the measurement 
     * subspace of the P space.
     * @param mtype return from getMeasurmentType.
     * @return 0 if projection is the identity matrix else 1.    
     */
    virtual int  getMtypeProjection(Cure::Matrix &project, 
				     unsigned short mtype){
      project.reallocate(Bdual.Columns);
      project=1;
      return 0;
    }
    /**
     * Adds the new collected Info to the cumulation and updates the 
     * Low dimensions.  This tracks the feature low dimensions using the 
     * latest measurement.  If the distance is not > LastDistance the
     * assumption is that this info is to replace any previous info at
     * this distance.
     * This means the addInfo update is done so that the last observation 
     * agrees with the prediction as far as the LowInfo can align it.
     * @param w Meaning can vary depending on subclass, typically a point or 
     *         bearing cloud.
     * @param type Meaning can vary depending on subclass.
     * @return 0 if ok 
     */
    virtual int  addInfo(const Cure::Matrix & w,const int type=0){return 0;}

    /**
     * Adds the new collected Info to the cumulation and updates the 
     * Low dimensions.  This tracks the feature low dimensions using the 
     * latest measurement.  If the distance is not > LastDistance the
     * assumption is that this info is to replace any previous info at
     * this distance.
     * This means the addInfo update is done so that the last observation 
     * agrees with the prediction as far as the LowInfo can align it.
     * @param w Meaning can vary depending on subclass, typically a point or 
     *         bearing cloud.
     * @param map2info the transformation between the continuous info frame
     * and the map frame.
     * @param type Meaning can vary depending on subclass.
     * @return 0 if ok 
     */
    virtual int  addInfo(const Cure::Matrix & w,
			 Cure::Transformation3D & map2info,
			 const int type=0){
      return addInfo(w,type);
    }
    /**
     * @param v all the info as in addInfo except now their are many rows.
     * @param threshold meaning depends on subclass.
     *
     */
    virtual int  addFullInfo(const Cure::Matrix & v,double threshold=0)
    {return addInfo(v);}
    /**
     * Adds dp to the p part of the coordinates for i = 0; i<Bdual.Columns.  
     * It then sets dp=0. It also recenters if needed.
     * @param dp Column matrix of the change in p-space dimension.
     */
    virtual void updateP(Cure::Matrix & dp);
    /**
     * Adds dxp to the p part of the coordinates and then sets dxp[i] to 0
     * for i = 0; i<B.Rows.  Also recenters if needed.
     * @param dxp Array of the change in p-space dimension.
     * @return number of P-space dimensions 
     */
    int updateP(double *dxp)
    {
      if (Bdual.Columns)
	{
	  Cure::Matrix dp(dxp,Bdual.Columns,1); 
	  updateP(dp);
	}
      return Bdual.Columns;
    }
    /**
     * This is a bit smarter than updateP, it expects dxp to be a pointer
     * to the start of the array that this->Index refers to.  So this
     * will index correctly even for non consecutive indecies.
     * It returns 0 if the indecies might not be consecutive and
     * the number of P-dim otherwise. 
     */
    int updateIndexedP(double *dxp){
      if (Consecutive)return updateP(dxp+Index[0]);
      Cure::Matrix dp(Bdual.Columns,1);
      for (int i=0; i<Bdual.Columns; i++)
	{
	  dp(i,0)=dxp[Index[i]];
	  dxp[Index[i]]=0;	    
	}
      updateP(dp);	
      return 0;
    }
    /**
     * Add dx to the coordinates x.
     */
    void updateX(const Cure::Matrix & dx);
    /**
     * Sets the feature coordinates to x
     * @param x a column vector with 3D followed by 2D followed by scalars. 
    */
    void setX(const Cure::Matrix & x);
    /** 
     * Gets the feature coordinates in x
     * @param x a column vector with 3D followed by 2D followed by scalars. 
     */
    void getX(Cure::Matrix & x)const;
    /*
     * Projects the coordinate vector with b.
     * @param p a column vector of Pdim.
     */
    virtual void getP(Cure::Matrix & p)const;
    
    /**
     * @return the flags for M-space rows that have frames attached to the 
     * feature. 
     */
    virtual unsigned short getFramed(){return 0;}
    /*
     * Recalculate the Bdual Matrix if needed.
     *
     */
    virtual void recenter(){}
    /*
     * This is where the sub class needs to overide.  This will
     * add dimensions to P if 'enough' info has been collected. 
     * @return number of new dimensions
     */
    virtual int extend(){return 0;}
    /**
     * Normally one need not override this.
     * p becomes (Bold,b_n)^T*Xo
     * dz= 0 = b_n*dXf = Jeo*Jos*dXs + Jeo*Jof*dXf
     * ==>Jeo=b_n*Jof^-1*Jos
     * ==>Jes=b_n*Jof^-1*Jos
     * Jep=I ndim x ndim  (ndim is return value # of new dim)
     * 
     * So for this psuedo measurement:
     * define H=(Jeo(Jos,Jop))=(b_n Jof^-1 Jos, I),
     * A = H C, where C is the covariance of the state vector after extending 
     *
     * S = H C H^T + (invcov)^-1.
     *
     * The Kalman Filter psuedo measurment gives then a change in C:
     * dC = -A^T S^-1 A. 
     *
     * @param b_n matrix to project to the new p dim
     * @param invcov Returns with the (covariance)^-1 of the psuedo measurement.
     * @return number of new dimensions
     */
    virtual int extend(Cure::Matrix & b_n,Cure::Matrix & invcov);
    /*
     * Prints info to the screen.
     */
    virtual void print(int level=0);
    /**
     * gets the B matrix.
     */
    virtual void getB(Cure::Matrix &b)const
    {b.transpose_(Bdual);}

    /**
     * The B matrix has rows that are orthogonal vectors with magnitude
     * 1/(scale factor).  That scale factor will either be 1 or a length 
     * in some direction in the Xf space.  For instance the length of a 
     * wall or line feature.
     * This method returns the matrix needed to calculate the current scale
     * factor from the current feature coordinates.  If the scale factors are
     * 1 then that row is left off of scales.  To find the coorespondence
     * between rows of scales and rows of B, examine the return value.
     * It is a list of binary flags that are 1 if the corresponding B row is
     * scaled and 0 otherwise.  Thus if zero is returned no scaling is done.
     * if 1 is returned the first row of b is scaled. If 4 is returned then
     * the third row of b is scaled, (5 -> 1st and 3rd scaled) and so on.
     *
     * @param scales formed into the Matrix that multiplies Xf to get the 
     *               scale vector for the rows of the B matrix.  It will have 
     *               rows for each scaled P dim.
     * @return binary flags to scaled rows of B.
     */
    virtual unsigned short getScales(Cure::Matrix &scales){
      scales.reallocate(0,FullDim);
      return 0;
    }
    /**
     * Write the MapFeature information to a text file.
     * @param fs the file to write to. It should be opened with std::ios::out
     */
    virtual void write(std::fstream &fs );
      
    /**
     * This initializes a MapFeture from the information in a text file.
     * @param version 1 only for now.
     * @param fs the file to read from. It should be opened with std::ios::in
     * @param readKey set this to false if key has been read allready.
     * @return 1 if fail else 0
    */
    virtual int read(int version,std::fstream &fs, bool readKey=true);
    /**
     * Stores the feature in a GenericData object
     * @param gd store this feature in.
     */
    virtual void get(GenericData &gd);
    /**
     * Restores the feature from a GenericData object.  If the MapPoints
     * of gd already exists on the Bank then the data is simply copied to
     * the existing object.   By already exists we mean that a MapPoint
     * is associated with gd's BankID, Key combination.  In other cases 
     * the Points on this feature are removed and the new ones added.  If the
     * removed points have no features they are deleted.
     *
     * @param gd the Data to restore this feature from.
     *
     * @return 0 if ok,1 if ds is formated wrong, MAP_OBJECT_INVALID
     * if the point's keys are associated with the wrong type MapObject.
     * 
     */
    virtual int set(GenericData &gd);
    
  protected:
    /**
     * Check if we might be using dense data.
     */
    virtual bool dimCheck(){
      if (Bdual.Columns==FullDim)return true;
      return false;
    }
    void init();
    void setNumberSharedPoints();
  };
  static inline MapFeature *getMapFeature(MapBank *b,long key)
  {
    if(!b)return 0;
    MapObject *mo=b->getMapObject(key);
    if (!mo) return 0;
    return mo->narrowFeature();
  }  
}
#endif
