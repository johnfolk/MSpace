// = AUTHOR(S)
//    John Folkesson
//    
//    August 1, 2007
//
//    Copyright (c) 2007 John Folkesson
//    
#ifndef CURE_RELATIVEPOINTFEATURE_HH
#define CURE_RELATIVEPOINTFEATURE_HH

#include "FeatureDescriptors.hh"
#include "FeatureFun.hh"
#include "MapPoseTree.hh"
#include "GenericData.hh"
namespace Cure {
  /**
   * Class managing relative point features in the map These features
   * are scalars which are functions of (phi, theta,range) the frame
   * of the point being an additional parameter of the feature. Here
   * phi is the angle of the bearing vector with positive z, theta is
   * the angle of projection of the bearing vector to the xy plane
   * with the positive x axis.
   *
   * The functions can then be various depending on m_Function 
   *
   * m_Function:
   *
   *  0 - Just Identity, funtion=(phi, theta,range)
   *
   *  1 - Inverse Range,  The function = (phi,theta, 1/range) 
   * 
   *  2 - Non-linear theta function=(phi,m_Scale(theta^3+m_Alpha theta), range)
   *
   * @author John Folkesson
   */
  class RelativePointFeature: public MapFeature
  {
  public:
    /**
     *  0 - Just Identity, funtion=(phi, theta,range)
     *
     *  1 - Inverse Range,  The function = (phi,theta, 1/range) 
     * 
     *  2 - Non-linear theta function=(phi,m_Scale(theta^3+m_Alpha
     *  theta), range)
     *  
     *  256  - BearingOnly
     *  512  - phiRange
     *  1024 - thetaRange
     */
    unsigned short m_Function;
    RelativePointFeature **CastPtr;
    long AllocatedIndex[3];
    double  AllocatedScalars[3];
    /** A list of FeatureDescriptors for this Feature*/
    FeatureDescriptorList m_Descriptors;
    unsigned short m_Bad;
    double m_Scale;
    double m_Alpha;
    long m_ReferenceKey;
    unsigned short m_RefBranch;
    
    /** 
     * flags 1=phi is weal 2 theta is weak 4 range is weak, and 8 is
     * phi not ready for initialization, 16 is theta not ready and 32
     * is range not ready.  So if the m_WeakInfo is 0 all dimensions are
     * determined by one single measurement good enough for initialization.
     */
    unsigned short m_WeakInfo;;
  protected:
  public:
    RelativePointFeature();

    /**
     *  @param wp The templete feature, CastPtr, WeightThreshold, 
     *  TriangleThreshold and DistanceThreshold are copied from the templete.
     */
    RelativePointFeature( RelativePointFeature *wp,MapBank *b=0);
    virtual ~RelativePointFeature(){
      cleanDescriptors();
    }
    
    virtual void narrow(){*CastPtr=this;}
    virtual MapFeature *copy();
    void setReference(long key,unsigned short branch){
      m_ReferenceKey=key;
      m_RefBranch=branch;
    }
    bool bearingOnly(){
      return (m_Function&0x100);
    }
    bool phiRange(){
      return (m_Function&0x200);
    }
    bool thetaRange(){
      return (m_Function&0x400);
    }
    /**
     * This re-assignes this feature a new reference frame.
     * 
     * @param key the key of the new reference MapPoseTree object.
     * @param branch the new branch to the sensor leaf
     */
    void resetReference(long key,unsigned short branch){
      double x[3];
      int r=getMapFrameCartesian(x);
      m_ReferenceKey=key;
      m_RefBranch=branch;
      if (r)return;
      Transformation3D *tran=referencePose();
      if (tran){
	tran->transform(x,x);
	setCenter(x);
      }
    }
    /**
     * This allows one to set the MSpace to any subspace of the 3 dof.
     * @return the number of new MSpace dimensions.
    */
    int setMSpace(bool usephi=false, bool usetheta=false, bool userange=false)
    {
      if (usephi)
	m_WeakInfo=(m_WeakInfo&0x37);
      if (usetheta)
	m_WeakInfo=(m_WeakInfo&0x2F);
      if (userange)
	m_WeakInfo=(m_WeakInfo&0x1F);
      return extend();
    }
    /**
     *  @param jac Jacobian between an earlier  Mspace and the current one
     *  @param bold the earlier B matrix.
     */
    void getMSpaceJacobian(Matrix &jac,Matrix &bold){      
      jac.multiply_(bold,Bdual);
      
    }
    /**
     *  @param jac Jacobian between an earlier  Mspace and the current one
     *  @param bold the earlier B matrix.
     */
    void getMSpaceJacobian(Matrix &jac,unsigned short oldmspace){      
      Matrix bold;
      getoldB(oldmspace,bold);
      getMSpaceJacobian(jac,bold);     
    }
    unsigned short getMSpaceFlags(){
      if (Bdual.Columns==0)return 0;
      unsigned short t=0;
      if (Bdual(0,0)==1)t=1;
      else if (Bdual(1,0)==1)t=2;
      else if (Bdual(3,0)==1)t=4;
      if (Bdual.Columns==1)return t;
      if (Bdual(1,1)==1)t=(t|2);
      else if (Bdual(3,1)==1)t=(t|4);
      if (Bdual.Columns==2)return t;
      if (Bdual(1,2)==1)t=(t|4);
      return t;
    }
    void  getoldB(unsigned short mSpaceFlags, Matrix &bold){
      int dim=0;
      if (mSpaceFlags&1)dim++;
      if (mSpaceFlags&2)dim++;
      if (mSpaceFlags&4)dim++;
      bold.reallocateZero(dim,3);
      dim=0;
      if (mSpaceFlags&1){
	bold(dim,0)=1;
	dim++;
      }
      if (mSpaceFlags&2){
	bold(dim,1)=1;
	dim++;
      }
      if (mSpaceFlags&4){
	bold(dim,2)=1;
	dim++;
      }
    }

    /**
     * @param xrel the carteasian coordinates relative to the sensor frame.
     */
    void setCenter(double xrel[3]){
      double rho2=xrel[0]*xrel[0]+xrel[1]*xrel[1];
      double rho=sqrt(rho2);
      Scalars[0]=atan2(rho,xrel[2]);
      Scalars[1]=atan2(xrel[1],xrel[0]);
      double r2=(rho2+xrel[2]*xrel[2]);
      Scalars[2]=sqrt(r2);
      applyFunction();
    }
    void applyFunction(){
      if (Scalars[1]>M_PI_2){
	Scalars[0]=-Scalars[0];
	Scalars[1]-=M_PI;
      }else if (Scalars[1]<-M_PI_2){
	Scalars[0]=-Scalars[0];
	Scalars[1]+=M_PI;
      }
       if (Scalars[2]<0){
	 Scalars[0]=M_PI-Scalars[0];
	 Scalars[1]+=M_PI;
	 Scalars[2]=-Scalars[2];
	 applyFunction();
       }
       if (m_Function&1){
	 if (Scalars[2]<1E-6)
	   Scalars[2]=1E6;
	 else
	   Scalars[2]=1/Scalars[2];
       } 
       if (m_Function&2){
	 Scalars[1]=(m_Scale*Scalars[1]*(m_Alpha+Scalars[1]*Scalars[1]));
       }
     }
    /*
     * Projects the coordinate vector with b.
     * @param p a column vector of Pdim.
     */
    void getP(Cure::Matrix & p){
      int dim=0;
      if (!(m_WeakInfo&8))dim++;
      if (!(m_WeakInfo&16))dim++;
      if (!(m_WeakInfo&32))dim++;
      p.reallocate(dim,1);
      dim=0;
      if (!(m_WeakInfo&8)){
	p(dim,0)=Scalars[0];
	dim++;
      }
      if (!(m_WeakInfo&16)){
	p(dim,0)=Scalars[1];
	dim++;
      }
      if (!(m_WeakInfo&32)){
	p(dim,0)=Scalars[2];
      }
    }
    void updateP(Cure::Matrix & dp){
      int dim=0;
      if (!(m_WeakInfo&8)){
	Scalars[0]+=dp(dim,0);
	dim++;
      }
      if (!(m_WeakInfo&16)){
	Scalars[1]+=dp(dim,0);
	dim++;
      }
      if (!(m_WeakInfo&32))
	Scalars[2]+=dp(dim,0);
      dp=0;
    }
    double iterateTheta(double &d, double u){
      double d2=d*d;
      double s=((3*d2+m_Alpha)*m_Scale);
      double du=(u-(m_Scale*d*(m_Alpha+d2)))/s;
      d+=(du);
      if (du<0)du=-du;
      return du;
    }
    /**
     * gets (phi,theta,range) values from Scalars.
     * @param s the values are returned here.
     */
    void getPolar(double s[3]){
      if (!m_Function&3){
	s[0]=Scalars[0];
	s[1]=Scalars[1];
	s[2]=Scalars[2];
      }else  {
	if (m_Function&1){
	  s[0]=Scalars[0];
	  s[1]=Scalars[1];
	  s[2]=1/Scalars[2];
	} 
	if (m_Function&2){
	  s[0]=Scalars[0];
	  s[2]=Scalars[2];
	  //cbrt is cubic root
	  s[1]=cbrt(Scalars[1]/m_Scale);
	  int i=10;
	  while (iterateTheta(s[1],Scalars[1])>1E-5){
	    i--;
	    if (i<0)break;
	  }
	}
      }    
    }
    /**
     * returns the xyz relative to the reference frame
     *
     */
    void getRelativeCartesian(double xrel[3]){
      double s[3];
      s[0]=0;
      s[1]=0;
      s[2]=0;
      getPolar(s);
      xrel[2]=s[2]*cos(s[0]);
      xrel[0]=s[2]*sin(s[0])*cos(s[1]);
      xrel[1]=s[2]*sin(s[0])*sin(s[1]);
    }
    /**
     * returns the xyz relative to the reference frame
     *
     */
    void getBearing(Vector3D b){
      double t=getTheta();
      double p=getPhi();
      b(2)=cos(p);
      b(0)=sin(p);
      b(1)=(b(0)*sin(t));
      b(0)*=(b(0)*cos(t));
    }
    Transformation3D * referencePose(){
      MapPoseTree *pt=getMapPoseTree(Bank,m_ReferenceKey);
      if (pt) return pt->getLeafPose(m_RefBranch);
      return 0;
    }
    int getReferencePose(Pose3D &p){
      MapPoseTree *pt=getMapPoseTree(Bank,m_ReferenceKey);
      if (pt){
	Transformation3D trans;
	unsigned short type;
	pt->getLeafPose(trans,type,m_RefBranch);
	p.Transformation3D::operator=(trans);
	p.setCovType(type);
	p.Time=pt->m_Time;
	return 0;
      }
      return 1;
    }

    int getMapFrameCartesian(double x[3]){
      getRelativeCartesian(x);
      Transformation3D *p=referencePose();
      if (p) p->invTransform(x,x);
      else return 1;
      return 0;
    }
    /**
     * @return 1 if inside else 0
     */
    int inside(double phirange[2], double thetarange[0],
		double rangerange[2]){
      double s[3];
      getPolar(s);
      if (s[0]<phirange[0])return 0;
      if (s[0]>phirange[1])return 0;
      if (s[1]<thetarange[0])return 0;
      if (s[1]>thetarange[1])return 0;
      if (s[2]<rangerange[0])return 0;
      if (s[2]>rangerange[1])return 0;
      return 1;
    }

    /**
     * @return 1 if inside else 0
     */
    int inside(double lowleft[2], double topright[2]){
      double s[3];
      getMapFrameCartesian(s);
      if (s[0]<lowleft[0])return 0;
      if (s[1]<lowleft[1])return 0;
      if (s[1]>topright[0])return 0;
      if (s[1]>topright[1])return 0;
      return 1;
    }
 
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
      if (LastDistance<0)return 0; 
      if (Bdual.Columns>0)return 0;
      double mindistance=distance-DistanceThreshold/2;
      if (LastDistance>(mindistance))return 0;
      return 1;
    }
    unsigned short getMeasurementType(unsigned short type){
      if ((type!=32)){
	if (Bdual.Columns==0)return 0;
	if (type==1)//phi theta Range 
	  return type;
	if (type==3)//phi theta^3 Range 
	  return type;
	if (type==5)//phi theta Range  (
	  return type;
      	if (type==7)//phi theta^3 Range  (
	  return type;
	return 0;
      }
      if (Bdual.Columns==3) type=1;
      else type=0;
      return type;  
    }
    
    
    void addDescriptor(FeatureDescriptor *fd)
    {
      m_Descriptors.addUnique(fd);
    }
    void cleanDescriptors()
    {
      for (FeatureDescriptorList *flist=&m_Descriptors;flist->Next;
	   flist=flist->Next){
	if (flist->Element->m_FeatureKey==Key)
	  flist->Element->m_FeatureKey=-1;
      }
      m_Descriptors.clean();
    }
    /**
     * This checks the angle between the bearings of the last descriptor
     * to all a the others
     * @param threshold a threshold on the sin^2 angle 
     * @return true if sin^2 angle exceeds the threshold.
     */
    bool checkLastBearingInfo(double threshold)
    {
      int row=m_Descriptors.count();
      if (row>1){
	Matrix v(row,3);
	int rw=0;
	for (FeatureDescriptorList *dlist=&m_Descriptors;
	     dlist->Next; dlist=dlist->Next)      
	  {
	    FeatureDescriptor *des=dlist->Element;
	    Transformation3D *pn=des->sensorPose();
	    if (pn){
	      pn->invRotate(des->m_Bearing.X,&v(rw,0));
	      rw++;
	    }
	  }
	v.Rows=rw;
	int i=rw-1;
	for (int j=0; j<i; j++){
	  double temp=v(i,0)*v(j,0)+v(i,1)*v(j,1)+v(i,2)*v(j,2);
	  temp=1-temp*temp;
	  if (temp>threshold){
	    return true;
	  }
	}
      }
      return false;
    }
   void setRange()
    {
      Scalars[2]=calcRange();
      applyFunction();
    }
    void setTheta()
    {
      Scalars[1]=calcTheta();
      applyFunction();
    }

    void setPhi()
    {
      Scalars[0]=calcPhi();
      applyFunction();
    }
   void setRange(double r)
    {
      Scalars[2]=r;
      applyFunction();
    }
    void setTheta(double t)
    {
      Scalars[1]=t;
      applyFunction();
    }

    void setPhi(double f)
    {
      Scalars[0]=f;
      applyFunction();
    }
    double getPhi(){
      return Scalars[0];
    }
    double getTheta(){
      if (!(m_Function&2))return Scalars[1];
      double s=cbrt(Scalars[1]/m_Scale);
	int i=10;
	while (iterateTheta(s,Scalars[1])>1E-5){
	  i--;
	  if (i<0)break;
	}
	return s;
    }
    double getRange(){
      if (!(m_Function&1))return Scalars[2];
      return 1/Scalars[2];
    }
    /**
     * This calculates an estimated range from the first descriptor's
     * sensor frame.
     * @return the range.
     */
    double calcTheta(){
      double sw=0;
      double swr=0; 
      FeatureDescriptor *fd=m_Descriptors.Element;
      if (m_Descriptors.Next)
	for (FeatureDescriptorList *dlist2=m_Descriptors.Next; 
	     dlist2->Next; dlist2=dlist2->Next)      
	  {
	    double wr=0;
	    sw+=fd->getWeightedTheta(wr,*dlist2->Element);
	    swr+=wr;
	  }
      if (sw>0)swr/=sw;
      return swr;
        }
    /**
     * This calculates an estimated range from the first descriptor's
     * sensor frame.
     * @return the range.
     */
    double calcPhi(){
      return 0;
    }
    /**
     * This calculates an estimated range from the first descriptor's
     * sensor frame.
     * @return the range.
     */
    double calcRange(){
      double sw=0;
      double swr=0; 
      FeatureDescriptor *fd=m_Descriptors.Element;
      if (m_Descriptors.Next)
	for (FeatureDescriptorList *dlist2=m_Descriptors.Next; 
	     dlist2->Next; dlist2=dlist2->Next)      
	  {
	    double wr=0;
	    sw+=fd->getWeightedRange(wr,*dlist2->Element);
	    swr+=wr;
	  }
      if (sw>0)swr/=sw;
      return swr;
    }
    /**
     * Finds the jacobian of polar (phi,theta,range) wrt the Scalars.
     *
     * @param jpolar diagonal elements of jacobian of polar returned here.
     * @return flags for which elements are not unity so 1 is jpolar[0]!=1
     *         4 is jpolar[2]!=1 and so on. 
     */
    unsigned short getPolarJac(double jpolar[3]){
      jpolar[0]=1;
      if (m_Function&1){
	 jpolar[1]=1;
	 if (Scalars[2]>1E-12)
	   jpolar[2]=-1/(Scalars[2]*Scalars[2]);
	 else jpolar[2]=1E24;
	 return 4;
      } 
      if (m_Function&2){
	jpolar[2]=1;
	double t=getTheta();
	jpolar[1]=1/(m_Scale*(m_Alpha+3*t*t));
	return 2;
      }
      jpolar[1]=1;
      jpolar[2]=1;
      return 0;
    
    }
    /**
     * This gets the world xyz and its jacobian wrt reference CovCoordinates
     *  and Scalars.
     * @param jac (dx/xref, dx/Scalars) where xref is the
     *        coordinate vector of the pose tree (covDim)
     * @param x the map frame xyz is returned here
     * @return 0 if ok else 1
     */
    int getCartesian(Matrix &jac,Matrix &x);
    virtual int extend();
    /*
     * Use this to unconditionlly extend to full dimension.
     * The Index[0]=-1 is the flag that this has been done.
     */
    void forceExtend()
    {
      if (Bdual.Columns==3)return;
      m_WeakInfo=(m_WeakInfo&7);
      extend();
    }
    bool testMatch(RelativePointFeature *mw, double tolerance);
    bool getC(RelativePointFeature *mw,Matrix &c)const;
   
    virtual void getB(Matrix &b)const;

    void remove(FeatureDescriptor* f){
      m_Descriptors.removeDescriptor(f);
      f->m_FeatureKey=0;
    }
    void add(FeatureDescriptor* f){
      m_Descriptors.addUnique(f);
      f->m_FeatureKey=this->Key;
    }
    bool empty(){
      if (m_Descriptors.Next)return  false;
      return true;
    }
    bool isBad(){return (m_Bad&1);}
    bool isHard(){return (m_Bad&2);}
    void setBad(){m_Bad=(m_Bad|1);} 
    void setHard(){m_Bad=(m_Bad|2);} 
    void setBad(unsigned short k);
    void addBad(unsigned short k){m_Bad=(m_Bad|k);}
    bool checkInfo(double threshold=.0004);
    void eat(RelativePointFeature *im){
      addBad(im->m_Bad);
      for (FeatureDescriptorList  *ai  = &im->m_Descriptors;ai->Next; ai=ai->Next)
	this->add(ai->Element);
      delete im;
    }
    FeatureDescriptor * operator()(unsigned short i){
      return m_Descriptors.get(i);
    }

    /**
     * @param fd the descriptor with the info
     * @param threshold the threshold for init.
     * @return 1 if Feature can now have extend succeed else 0.
     */
    int addInfo(FeatureDescriptor *fd, double bearingthreshold,
		int numberThreshold=0){
      addDescriptor(fd);
      if (Bdual.Columns<3)
	if (m_Descriptors.count()>numberThreshold)
	  if (checkLastBearingInfo(bearingthreshold)){
	    if (m_WeakInfo&4)
	      setRange();
	    if (m_WeakInfo&2)
	      setTheta();
	    if (m_WeakInfo&1)
	      setPhi();
	    m_WeakInfo=(m_WeakInfo&7);
	    return 1;
	  }
      return 0;
    }
    
    /**
     * Adds the new collected Info to the cumulation and updates the 
     * LowInfo dimensions.
     * @param v a row matrix(1,8)  with this structure:
     * distance, w, x_i, s_i  (0..7)
     * w is a relative weight estimated to approximate the amount of
     * infomation contained in this observation (inverse variance)
     * x_i is the xyz of the camera when the observation was made.
     * s_i is a vector from x_i towards the point (magnitude ignored)
     * @param type 0.  
     * @return 0
     */
    int  addInfo(FeatureDescriptors *des,
		 const Cure::Matrix & v,double bearingthreshold,
		 int numberThreshold=0,
		 double range=5,
		 const int type=0)
    {  
      if (Bdual.Columns==3)return 0;
      if (v.Columns<8)return 1;
      if ( LastDistance<=v(0,0))
	LastDistance=v(0,0);
      else return 0;
      double b[3];
      double d=v(0,7)*v(0,7)+v(0,4)*v(0,4)+v(0,5)*v(0,5);
      d=sqrt(d);
      b[0]=v(0,5)/d;
      b[1]=v(0,6)/d;
      b[2]=v(0,7)/d;
      FeatureDescriptor *fd=new FeatureDescriptor(b,des,range);
      return addInfo(fd,bearingthreshold, numberThreshold);
    }
    int tryInitialize(double matchthreshold,
		      double radsqerr,LongList *feats=0);
 
    /**
     * This initializes based on a pixel and focal length
     * @param centerpixels the xz of the pixel 
     * @param focallength the camera focal length
     * @param distanceguess  relative to camera origin.
     */
    void initializeFromPixels(const double centerpixels[2],
			      double focallength,
			      double distanceguess);

    int merge(RelativePointFeature *mf, unsigned short type=1);
  
    /**
     * Write the MapPoint information to a text file.
     * @param fs the file to write to. It should be opened with std::ios::out
     */
    void write(std::fstream &fs ){
      if (Bdual.Columns<3)return;
      MapFeature::write(fs);
      fs<<m_Function<<" "<<m_Scale<<" "<<m_Alpha<<" "<<m_WeakInfo<<" "<<m_Bad<<"\n";
    }
      
    /**
     * This initializes a MapFeture from the information in a text file.
     * @param version 1 only for now.
     * @param fs the file to read from. It should be opened with std::ios::in
     * @param readKey set this to false if key has been read allready.
     * @return 1 if fail else 0
    */
    virtual int read(int version, std::fstream &fs, bool readKey=true){
      if (version!=1)return 1;
      int ret=MapFeature::read(1,fs,readKey);
      if (ret)return ret;
      fs>>m_Function>>m_Scale>>m_Alpha>>m_WeakInfo>>m_Bad;
      recenter();
      return 0;
    }
    /**
     * Stores the feature in a DataSet object
     * @param gd the Data to store this feature in.
     */
    virtual void get(GenericData &gd);
    /**
     * Restores the feature from a DataSet object.  It looks for the
     * ds's BankID, ReferenceKey combination and if found uses the Key
     * from this objects's Bank as ReferenceKey.  Otherwise it makes a 
     * new MapPoseTree with and takes that as the Reference.  This leaves
     * the new PoseTree uninitialized and the RefBranch might point to a
     * non-existing leaf.  This also will clean the descriptor list.
     *
     * @param gd the Data to restore this feature from.
     * @return 1 if fails else 0;
     */
    virtual int set(GenericData &gd);
    void getGenericData(GenericData &gd);
  protected:
  };

} // namespace Cure

#endif
