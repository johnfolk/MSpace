// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2007 John Folkesson
//    

#include "RelativePointFeatureHelper.hh"
#include "SensorData.hh"
#include "CureDebug.hh"

#ifndef DEPEND
#include <sstream>
#endif

using namespace Cure;

int RelativePointFeatureHelper::addDistance(Match & mat,double distance)
{
  //Must set the bounding box of the measurments a la ImagePlane
  //We must set Metric.
  Measurement &m=*mat.Measure; 
  mat.PathDistance=distance;
  if ((m.MeasurementType==32)||
      (m.MeasurementType==33)){
    double b[3];
    b[0]=-m.V(0,0);
    b[1]=m.V(2,0);
    b[2]=-m.V(1,0);
    double rho=b[0]*b[0]+b[1]*b[1];
    double d=rho+b[2]*b[2];
    d=sqrt(d);
    rho=sqrt(rho);
    rho/=d;
    m.V(0,0)=b[0]/d;
    m.V(1,0)=b[1]/d;
    m.V(2,0)=b[2]/d;
    
    m.Z.reallocate(2,1);
    m.Z(0,0)=atan2(rho,m.V(2,0));
      m.Z(1,0)=atan2(m.V(1,0),
		     m.V(0,0));
      m.CovV/=d;
      m.CovV.swapRows(2,1);
      m.CovV.swapColumns(2,1);
      m.CovV(1,0)*=-1;
      m.CovV(1,2)*=-1;
      m.CovV(0,1)*=-1;
      m.CovV(2,1)*=-1;
      m.MeasurementType+=2;
  }
  //Z=(phi,theta,range,
    //V=(phi,f(theta),range).
  
  if ((m.MeasurementType==1)||(m.MeasurementType==5)){
    m.Z=m.V;
    if (TempletFeature.m_Function&2){
      m.MeasurementType=(m.MeasurementType|2);
	{
	  m.Z(1,0)*=m.Z(1,0);
	  m.V(1,0)=(TempletFeature.m_Scale*m.V(1,0)*(TempletFeature.m_Alpha+
						     m.Z(1,0)));
	  double b=2;//  (theta cutoff)^2 = b cov_theta  b is in (1,3)
	  double cv=m.CovV(1,1);
	  if ((m.W.Columns<10)||(m.W.Rows<1))m.W.grow(1,10);
	  m.W(0,8)=TempletFeature.m_Scale;
	  m.W(0,8)=TempletFeature.m_Alpha;
	  double scale=(m.W(0,8)*(b*cv));
	  m.CovV(1,1)*=scale*scale;
	  m.CovV(1,0)*=scale;
	  m.CovV(1,2)*=scale;
	  m.CovV(0,1)*=scale;
	  m.CovV(2,1)*=scale;
	}
    }
  }
  if ((m.MeasurementType==3)||(m.MeasurementType==7))
    {//theta
      if (m.W.Columns>9){
	  TempletFeature.m_Alpha=m.W(0,9);
	  TempletFeature.m_Scale=m.W(0,8);
      }
      m.Z=m.V;
      double s=cbrt(m.Z(1,0)/TempletFeature.m_Scale);
      int i=10;
      while (TempletFeature.iterateTheta(s,m.V(1,0))>1E-5){
	i--;
	if (i<0)break;
      }
      m.Z(1,0)=s;
    }
  if (TempletFeature.bearingOnly()){
    mat.Thresholds.reallocate(1,2);
    mat.Thresholds(0,0)=m_ZBox[0];
    mat.Thresholds(0,1)=m_ZBox[1];
    mat.Metric.reallocate(2);
    mat.Metric=m_BearingOnlyMetric[0];
    mat.Metric(1,1)=m_BearingOnlyMetric[1];
    m.BoundingBox.reallocate(2);
    m.BoundingBox(0,0)=m.Z(0,0)
      -m_ZBox[0];
    m.BoundingBox(1,0)=m.Z(0,0)
      +m_ZBox[0];
    m.BoundingBox(0,1)=m.Z(1,0)
      -m_ZBox[1];
    m.BoundingBox(1,1)=m.Z(1,0)
      +m_ZBox[1];
  }else if (TempletFeature.phiRange()){
    mat.Thresholds.reallocate(1,2);
    mat.Thresholds(0,0)=m_ZBox[2];
    mat.Thresholds(0,1)=m_ZBox[3];
    mat.Metric.reallocate(3);
    mat.Metric=m_PhiRangeMetric[0];
    mat.Metric(1,1)=m_PhiRangeMetric[1];
    mat.Metric(2,2)=m_PhiRangeMetric[2];
    m.BoundingBox.reallocate(2);
    m.BoundingBox(0,0)=m.Z(0,0)
      -m_ZBox[2];
    m.BoundingBox(1,0)=m.Z(0,0)
      +m_ZBox[2];
    m.BoundingBox(0,1)=m.Z(2,0)
      -m_ZBox[3];
      m.BoundingBox(1,1)=m.Z(2,0)
      +m_ZBox[3];
  }else if (TempletFeature.thetaRange()){
    mat.Thresholds.reallocate(1,2);
    mat.Thresholds(0,0)=m_ZBox[4];
    mat.Thresholds(0,1)=m_ZBox[5];
    mat.Metric.reallocate(3);
    mat.Metric=m_ThetaRangeMetric[0];
    mat.Metric(1,1)=m_ThetaRangeMetric[1];
    mat.Metric(2,2)=m_ThetaRangeMetric[2];
    m.BoundingBox.reallocate(2);
    m.BoundingBox(0,0)=m.Z(0,0)
      -m_ZBox[4];
    m.BoundingBox(1,0)=m.Z(0,0)
      +m_ZBox[4];
    m.BoundingBox(0,1)=m.Z(1,0)
      -m_ZBox[5];
    m.BoundingBox(1,1)=m.Z(1,0)
      +m_ZBox[5];
  }
  return 0;
}

int RelativePointFeatureHelper::findMatches(FeatureDescriptors &currentdes,
					    FeatureDescriptors **des,
					    int n,
					    float errorlimit, 
					    double radsqerr, 
					    double coarseEnergyLimit)
{
  if (currentdes.m_NumberDesc==0)return 0;
  if (n<1)return 0;
  Transformation3D* poses[n];
  for (int i=0;i<n;i++){
    poses[i]=des[i]->m_SensorPose;
    if (!poses[i])return 0;
  }
  int num=0;
  Transformation3D *pose=currentdes.m_SensorPose;
  if (!pose)return 0;
  LongList keys;
  for (int k=0;k<n;k++)
    {
      for (int i=0; i<currentdes.m_NumberDesc; i++)
	{
	  if (currentdes.m_Desc[i]->m_FeatureKey==-1){
	    FeatureDescriptor *f=
	      des[k]->match(currentdes, i,
			   errorlimit,radsqerr,coarseEnergyLimit,5);
	    if (f){
	      long lg=f->m_FeatureKey;
	      if (lg==-1){
		RelativePointFeature *mf
		  =makeRelativePointFeature(f);
		lg=mf->Key;
		NearKeys.add(lg);
		Near.add(mf);
	      }
	      if (keys.addUnique(lg)){
		currentdes.m_Desc[i]->m_FeatureKey=lg;
		RelativePointFeature *mf=
		  getRelativePointFeature(lg);
		mf->m_Descriptors.add(currentdes.m_Desc[i]);
		num++;
	      }
	    }
	  }
	}
    }
  return num;
}

int RelativePointFeatureHelper::lucasKanadeTrack(FeatureDescriptors **des, 
						 int n,
						 unsigned short windowsize,
						 int offset)
{ 
  
  float threshold=(2*windowsize+1)*1*(1024);
  if (n<3)return 0;
  if(des[0]->m_NumberDesc==0)return 0;
  if(des[n-1]->m_NumberDesc==0)return 0;
  Transformation3D * poses[n];
  poses[0]=des[0]->m_SensorPose;
    if (!poses[0])return 0;
  for (int i=1;i<n;i++){
    poses[i]=des[i]->m_SensorPose;
    if (!poses[i])return 0;
    /*
    if ((poses[i-1]->Time>poses[i]->Time)){
      std::cerr<<"I do not handle tracking forward\n";
      return 0;
    }
    */
  }
  int ret=0;
  for (int i=0;i<des[0]->m_NumberDesc;i++)
    {
      FeatureDescriptor *f1 =des[0]->m_Desc[i];
      if (f1)
	{
	  long key=f1->m_FeatureKey;
	  RelativePointFeature *im=
	    getRelativePointFeature(key);
	  if (im)
	    if (!im->isBad()){
	      FeatureDescriptor *f2 =des[n-1]->findFeature(key);
	      if (f2)
		{
		  short pix[2*n+2];
		  float pos[2*n+2];
		  f1->getPixels(pix+2*n);
		  f2->getPixels(pix);
		  f1->getImagePosition(pos+2*n);
		  f2->getImagePosition(pos);
		  float factor[2];
		  factor[0]=((float)(pos[2*n]-pos[0]))/(float)(n-1);
		  factor[1]=((float)(pos[2*n+1]-pos[1]))/(float)(n-1);
		  bool trytrack=true;
		  for (int j=1;j<n;j++)
		    {
		      pos[2*j]=pos[0]+factor[0]*(float)j;
		      pos[2*j+1]=pos[1]+factor[1]*(float)j;
		      float d=pos[2*j]+.5;
		      pix[2*j]=(short)d;
		      d=pos[2*j+1]+.5;
		      pix[2*j+1]=(short)d;
		      if (!isWithin(pix+2*j,windowsize+2))trytrack=false;
		    }
		  Cure::Matrix displace[n-1];
		  for (int j=1 ;((j<n)&&trytrack);j++){
		    
		  float delta=lktrack(j-n-offset,j-n-offset+1,
				    pix+2*j-2,pix+2*j,
				    displace[j-1],windowsize);
		  if ((delta<0)||(delta>threshold))
		    trytrack=false;
		  else{
		    pos[2*j]=pix[2*j]+displace[j-1](0,0);
		    pos[2*j+1]=pix[2*j+1]+displace[j-1](1,0);
		    if (j<(n-1)){
		      pos[2*j+2]=pos[2*j]+(pos[2*n]-pos[2*j])/(float)(n-1-j);
		      pos[2*j+3]=pos[2*j+1]+(pos[2*n+1]-pos[2*j+1])/(float)(n-1-j);
		      float d=pos[2*j+2]+.5;
		      pix[2*j+2]=(short)d;
		      d=pos[2*j+3]+.5;
		      pix[2*j+3]=(short)d;
		    }
		  }
		}
		if (trytrack){
		  float tst=pos[2*n]-pos[2*n-2];
		  if ((tst>1.5)||(tst<-1.5))trytrack=false;
		}
		if (trytrack){
		  float tst=pos[2*n+1]-pos[2*n-1];
		  if ((tst>1.5)||(tst<-1.5))trytrack=false;
		}
		if (trytrack){
		  ret++;
		  for (int j=1;j<(n-1);j++)
		    {
		      FeatureDescriptor *fd=new FeatureDescriptor();
		      des[n-j-1]->add(fd);
		      im->add(fd);
		      fd->m_Bearing(1)=1;
		      fd->m_Bearing(0)=(pos[2*j]
				      -m_BearingZero[0])/m_FocalLength[0];
		      fd->m_Bearing(2)=-(pos[2*j+1]/
				       -m_BearingZero[1])/m_FocalLength[1];
		      double db=fd->m_Bearing*fd->m_Bearing;
		      db=sqrt(db);
		      fd->m_Bearing/=db;
		    }
		}
	      }
	  }
	}
    }
  return ret;
}
float RelativePointFeatureHelper::lktrack(int originalindex, int trackedindex, 
					  short originalpix[2],
					  short trackedpix[2],
					  Cure::Matrix & displace,
					  unsigned short windowsize){
  if (!isWithin(originalpix,windowsize+1))return -1;
  displace.reallocate(2,1);
  originalindex+=m_Current;
  trackedindex+=m_Current;
  originalindex=originalindex%m_Images.getNumberOfImages();
  trackedindex=trackedindex%m_Images.getNumberOfImages();
  if (originalindex<0)originalindex+=m_Images.getNumberOfImages();
  if (trackedindex<0)trackedindex+=m_Images.getNumberOfImages();
  //  IplImage *img=m_Undistorted[originalindex];
  //IplImage *imgt=m_Undistorted[trackedindex];
  //if ((!img)||(!imgt))return 1;
  Cure::ShortMatrix templet(2*windowsize+3);
  getTemplet(originalindex,templet,originalpix[0],
	     originalpix[1]);;
  float it=3;
  while (it>0){
    float delta=getDisplacement(trackedindex,displace,
				templet,
				trackedpix);
    if (delta<0)return -1;
    bool redo=false;
    if (displace(0,0)>.5){
      if (displace(0,0)>(double)windowsize)return -1;
      float dp=displace(0,0)+.5;
      short sdp=(short)dp;
      trackedpix[0]+=sdp;
      displace(0,0)-=sdp;
      redo=true;
    }
    else if (displace(0,0)<-.5){
      if (displace(0,0)<-(double)windowsize)return -1;
      float dp=displace(0,0)-.5;
      short sdp=(short)dp;
      trackedpix[0]+=sdp;
      displace(0,0)-=sdp;
      redo=true;
    } 
    if (displace(1,0)>.5){
      if (displace(1,0)>(double)windowsize)return -1;
      float dp=displace(1,0)+.5;
      short sdp=(short)dp;
      trackedpix[1]+=sdp;
      displace(1,0)-=sdp;
      redo=true;
    }
    else if (displace(1,0)<-.5){
      if (displace(1,0)<-(double)windowsize)return -1;
      float dp=displace(1,0)-.5;
      short sdp=(short)dp;
      trackedpix[1]+=sdp;
      displace(1,0)-=sdp;
      redo=true;
    } 
    if (!redo)
      return delta;
    it--;
  }
  return -1;
} 

float RelativePointFeatureHelper::getDisplacement(int imgindex,
						  Cure::Matrix &displace,
						  Cure::ShortMatrix &templet,
						  short center2[2])
{
  int off=(templet.Rows-3)/2;
  if (!isWithin(center2,off+1))return -1;
  int iw2=(center2[0]-off);
  Matrix e(2,1);
  float delta=0;
  Matrix m(2);
  e=0;
  int top=2*off+1;
  int dx=1;
  if (displace(0,0)<0)dx=0;
  int dy=1;
  if (displace(1,0)<0)dy=0;
  for (int ti=0;ti<top;ti++,iw2++)
    {
      int j2=center2[1]-off;
      for (int tj=0;tj<top;tj++,j2++)
	{
	  short t=templet(ti+1,tj+1);
	  t-=getPixelValue(iw2,j2,imgindex);
	  double gx=(double)(templet(ti+1+dx,tj+1)-templet(ti+dx,tj+1));
	  double gy=(double)(templet(ti+1,tj+1+dy)-templet(ti+1,tj+dy));
	  e(0,0)+=(double)t*gx;
	  e(1,0)+=(double)t*gy;
	  m(0,0)+=gx*gx;
	  m(0,1)+=gx*gy;
	  m(1,1)+=gy*gy;
	  delta+=((float)t*(float)t);
	}
    }
  m(1,0)=m(0,1);
  bool redo=false;
  if (m.invert()){
    Matrix lambda, ev;
    m.symmetricEigen(lambda,ev);
    m(0,0)=ev(0,0)*ev(0,0)*lambda(0,0);
    m(1,0)=ev(0,0)*ev(1,0)*lambda(0,0);
    m(0,1)=m(1,0);
    m(1,1)=ev(1,0)*ev(1,0)*lambda(0,0);
    redo=true;
  }
  displace.multiply_(m,e);
      
  if ((((dx==1)&&(displace(0,0)<0))||((dx==0)&&(displace(0,0)>0)))
      ||(((dy==1)&&(displace(1,0)<0))||((dy==0)&&(displace(1,0)>0))))redo=true;
 
  if (redo){
    Matrix dis=displace;
    iw2=(center2[0]-off);
    e=0;
    m=0;
    dx=1;
    if (displace(0,0)<0)dx=0;
    dy=1;
    if (displace(1,0)<0)dy=0;
    for (int ti=0;ti<top;ti++,iw2++)
	{
	  int j2=center2[1]-off;
	  for (int tj=0;tj<top;tj++,j2++)
	    {
	      short t=templet(ti+1,tj+1);
	      t-=getPixelValue(iw2,j2,imgindex);
	      double gx=(double)(templet(ti+1+dx,tj+1)-templet(ti+dx,tj+1));
	      double gy=(double)(templet(ti+1,tj+1+dy)-templet(ti+1,tj+dy));
	      e(0,0)+=(double)t*gx;
	      e(1,0)+=(double)t*gy;
	      m(0,0)+=gx*gx;
	      m(0,1)+=gx*gy;
	      m(1,1)+=gy*gy;
	    }
	}
    m(1,0)=m(0,1);
    if (m.invert()){
      Matrix lambda, ev;
      m.symmetricEigen(lambda,ev);
      m(0,0)=ev(0,0)*ev(0,0)*lambda(0,0);
      m(1,0)=ev(0,0)*ev(1,0)*lambda(0,0);
      m(0,1)=m(1,0);
      m(1,1)=ev(1,0)*ev(1,0)*lambda(0,0);
    }
    displace.multiply_(m,e);

    if (((dx==1)&&(displace(0,0)<0))||((dx==0)&&(displace(0,0)>0)))
      {
	displace(0,0)+=dis(0,0);
	displace(0,0)/=2;
      }  
    if (((dy==1)&&(displace(1,0)<0))||((dy==0)&&(displace(1,0)>0)))
      {
	displace(1,0)+=dis(1,0);
	displace(1,0)/=2;
      } 
  }
  delta-=(float)(e(0,0)*displace(0,0)+e(1,0)*displace(1,0));
  if (delta<0)delta=0;
  return delta;  
}
    

MapFeature * RelativePointFeatureHelper::makeMapFeature
(Cure::Pose3D &sensorpose,
 Cure::Measurement &m,
 bool addToVisable)
{
  double d=0;
  m.print();
  sensorpose.print();
  if (m.MeasurementType==33) d=m.W(0,5);

  std::cerr<<"make point ";
  RelativePointFeature *mf=makeRelativePointFeature();
  if (m.MeasurementType==1){
    mf->setPhi(m.V(0,0));
    mf->setTheta(m.V(1,0));
    mf->setRange(m.V(2,0));
  }
  if (m.MeasurementType==33)
    {
      std::cerr<<"forceextend ";
      CureCERR(40) << "Forcing extend because meastype 33\n";
      mf->forceExtend();
      m.MeasurementType=32;
      addTrackedKey(m.Key,mf->Key);
    }
  mf->print();
  std::cerr<<"Add Visable"<<addToVisable<<" ";
  getchar();
  if (addToVisable){
    PosedRelativePointFeature *pf=makePosedRelativePointFeature(mf);
    if (pf){
      pf->print();
      pf->transform(sensorpose, sensorpose.getCovType());
      PosedFeatures.add(pf);
      Visable.add(pf);
    }
  }
  return mf;
}

bool 
RelativePointFeatureHelper::supportsSubconfig(int sc)
{
  return (1 <= sc && sc <= 4);
}
void 
RelativePointFeatureHelper::printConfiguration()
{
  std::cerr << "Configured RelativePointFeatureHelper with "
	    << std::endl;
}
int
RelativePointFeatureHelper::config(const std::string &arglist)
{
  std::istringstream str(arglist);
  int version = -1;

  if ( !(str >> version)) {
    CureCERR(20) << "Failed to read version number for config params list\n";
    return 1;
  }

  int ret = 0;
  switch (version) {
  case 1:
    ret = configVer1(arglist);
    break;
  default:
    CureCERR(20) << "Cannot handle config version " << version << std::endl;
    return 1;
  }
  return ret;
}
int 
RelativePointFeatureHelper::configVer1(const std::string &arglist) 
{
  std::istringstream str(arglist);
  int version = 0;
  int subcfg = -1;
      if ( !(str >> version>> subcfg)) {
    CureCERR(20) << "Failed to read subcfg number for config params list\n";
    return 1;
  }
  switch (subcfg) {
  case 0:
    {
      double sw, sh, f, pW, pH;
      if (str >> RoughSearchRange
	  >> NearSearchRange
	  >> NewFeatureRange
	  >> MatchThreshold
	  >> TempletFeature.DistanceThreshold
	  >> sw >> sh >> f >> pW >> pH) {
	double bl[2], tr[2];
	bl[0] = -sw;
	bl[1] = -sh;
	tr[0] = sw;
	tr[1] = sh;
	setImage(bl, tr, f, pW, pH);
	return 0;
      }
    }
      break;
  case 1:
    if (str >> RoughSearchRange
        >> NearSearchRange
        >> NewFeatureRange
        >> MatchThreshold) {
      return 0;
    }
    break;
  case 2:
    if (str >> 
        TempletFeature.DistanceThreshold) {
        return 0;
    }
    break;
  case 4:
    double sw, sh, f, pW, pH;
    if (str >> sw >> sh >> f >> pW >> pH) {
      double bl[2], tr[2];
      bl[0] = -sw;
      bl[1] = -sh;
      tr[0] = sw;
      tr[1] = sh;
      setImage(bl, tr, f, pW, pH);
      return 0;
    }
    break;
  default:
    CureCERR(20) << "Cannot handle subcfg " << subcfg << std::endl;
    return 1;
  }
  // Could not read all parameters but the ones it did read will
  // change and the rest will have thier default values
  return 1;
}


//merge starts by match to a expected measurement
//then extra critera are tested 
//So for a wall that might be overlap (RoughSearchRange<0),
//for a HLine it might be similar hieght or hieght not known.
//Then the SLAM will need to deal with enforcing the 
//Constraints.  We need CP=b for the P's
/*
 * try finding a merge candidates from search list
 * If found remove them from searchlist and return them in 
 * pl1 and pl2.
 * @return 1 after finding one 0 otherwise
 */
unsigned short  RelativePointFeatureHelper::findMergeFeatures
(PosedFeature * pfp[2],
 PosedFeatureList * searchList)
{
  PosedRelativePointFeature *pl[2];
  if (searchList==0)searchList=&Visable;
  double thresholds[2];
  thresholds[0]=.005;
  thresholds[1]= 2*PixelInfo[1];
  for (PosedFeatureList *mlist=searchList; mlist->Next; 
       mlist=mlist->Next)
    {
      PosedFeature *pf=mlist->Element;
      if (pf)
	{
	  pl[0]=castPosedRelativePointFeature(pf);
	  if (pl[0])
	    {  
	      for (PosedFeatureList *mlist2=mlist->Next; mlist2->Next; 
		   mlist2=mlist2->Next)
		{
		  pf=mlist2->Element;
		}
	    }
	}
    } 
  return 0;
}


int  RelativePointFeatureHelper::getMergeContstraint
(Matrix &a, Matrix &e,Matrix & cov,
 double distance,
 unsigned short typ, PosedFeature * pl[2])
{
  a.reallocate(0);
  e.reallocate(0,1);
  cov.reallocate(0);
  return 1;
}
void 
RelativePointFeatureHelper::insertPoseTopath(Cure::Matrix &path, 
                                     Cure::Matrix &deltapath,
                                     const Transformation3D  &sensorpose,
                                     int index)
{
  int n=index+1;
  if (n>path.Rows)
    {
      path.grow(n,3);
      deltapath.grow(n,3*n);
    }
  n=path.Rows;
  int r3=3*index;
  int k=0;
  sensorpose.getXYZ(&path(index,0));
  for (int j=0; j<n; j++,k+=3)
    {
      deltapath(index,k)=path(index,0)-path(j,0);
      deltapath(index,k+1)=path(index,1)-path(j,1);
      deltapath(index,k+2)=path(index,2)-path(j,2);
      deltapath(j,r3)=-deltapath(index,k);
      deltapath(j,r3+1)=-deltapath(index,k+1);
      deltapath(j,r3+2)=-deltapath(index,k+2);
    }
}

void 
RelativePointFeatureHelper::makeBearings(Cure::Matrix &bearings,
				 bool *frames,
				 Cure::Matrix &v,
				 Cure::Transformation3D  *sensorposes,
				 const double distanceguess)
{
  int n=v.Rows;
  bearings.reallocate(n,3);
  // Unit Vector that points at the detected feature in the sensor frame 
  for (int i=0; i<n; i++)
    { 
      if (v(i,2)>0){

	bearings(i,0)=-v(i,0);
	bearings(i,1)=v(i,2)/(1-v(i,2)/distanceguess);
	bearings(i,2)=-v(i,1);
	// Rotate the direction vector to the info frame
	sensorposes[i].invRotate(&bearings(i,0),&bearings(i,0));
	double d=bearings(i,0)*bearings(i,0)+bearings(i,1)*bearings(i,1)
	  +bearings(i,2)*bearings(i,2);
	if (d<1E-40)frames[i]=false;
	else{    
	  frames[i]=true;
	  d=sqrt(d);
	  bearings(i,0)/=d;
	  bearings(i,1)/=d;
	  bearings(i,2)/=d;
	}
      }
      else
	frames[i]=false;
    }
} 
void 
RelativePointFeatureHelper::predictPixels(Cure::Matrix &vhat,
				  Cure::Transformation3D  *sensorposes,
				  Cure::Matrix &meanpoint)
{  
  int n=vhat.Rows;
  double pt[3];
  for (int i=0; i<n; i++)
    { 
      if (vhat(i,2)>0){
	sensorposes[i].transform(&meanpoint(0,0),pt);
	double d=(pt[1]-vhat(i,2));
	vhat(i,0)=-pt[0]*vhat(i,2)/d;
	vhat(i,1)=-pt[2]*vhat(i,2)/d;
      }
    }
} 

void 
RelativePointFeatureHelper::calcInnovations(Cure::Matrix &innovation,
				    Cure::Matrix &vhat,
				    Cure::Matrix &v, int currentframe)
{  

  int r=0;
  int n=vhat.Rows;
  innovation.reallocate(n,3);
  for (int i=currentframe; i<n; i++)
    { 
      if (vhat(i,2)>0){
	innovation(r,0)=v(i,0)-vhat(i,0);
	innovation(r,1)=v(i,1)-vhat(i,1);
	innovation(r,2)=i-currentframe;	
	r++;
      }
    }
  for (int i=0;i<currentframe; i++)
    { 
      if (vhat(i,2)>0){
	innovation(r,0)=v(i,0)-vhat(i,0);
	innovation(r,1)=v(i,1)-vhat(i,1);
	innovation(r,2)=i-currentframe+n;	
	r++;
      }
    }
  innovation.offset(0,0,r,3);
} 

void 
RelativePointFeatureHelper::calcTimeDerivative(Cure::Matrix &dx_dt,
				       Cure::Matrix &x)
{  
  int n=x.Rows;
  dx_dt.reallocate(n-1,3);
  for (int i=1; i<n; i++)
    { 
      double dt=x(i,2)-x(i-1,2);
      if (dt>0){
	dx_dt(i-1,0)=x(i,0)-x(i-1,0)/dt;
	dx_dt(i-1,1)=x(i,1)-x(i-1,1)/dt;
	dx_dt(i-1,2)=x(i,2)-(dt/2);
      }
    }
} 


double 
RelativePointFeatureHelper::triangulateBuffer(Cure::Matrix &meanpoint,
				      Cure::Matrix &v,
				      Cure::Transformation3D  *sensorposes,
				      const double weightThreshold, 
				      const double sqDistanceThreshold,
				      const double mseThreshold,
                                      const double secDerThreshold,
				      const double distanceguess)
{
  Cure::Matrix path,deltapath;
  for (int i=0; i<v.Rows;i++)
    insertPoseTopath(path, 
		     deltapath,
		     sensorposes[i],
		     i);
  //  calcPath(path,deltapath,sensorposes, v.Rows);
  Matrix bearings;
  bool frames[v.Rows];
  makeBearings(bearings,frames,v,sensorposes,distanceguess);
  double d= triangulateBuffer(meanpoint,path,deltapath,bearings,frames,
			   weightThreshold,sqDistanceThreshold,
			   mseThreshold);
  if (d>0){
    for (int i=0; i<v.Rows; i++)
      if (!frames[i])
        if (v(i,2)>0){
	  v(i,2)=-1;
	  d=-2;
	}

    if (secDerThreshold > 0) {
      double sd = makeDerivateCheck(meanpoint, v, sensorposes);
      if (sd > secDerThreshold) return -4;
    }
  }
  

  return d;
}

double 
RelativePointFeatureHelper::triangulateBuffer(Cure::Matrix &meanpoint,
				      Cure::Matrix &path, 
				      Cure::Matrix &deltapath, 
				      Cure::Matrix &v,
				      Cure::Transformation3D  *sensorposes,
				      const double weightThreshold, 
				      const double sqDistanceThreshold,
				      const double mseThreshold,
                                      const double secDerThreshold,
				      const double distanceguess,
                                      double *retWeight,
                                      double *retSqDist,
                                      double *retMse,
                                      double *retSecDer)
{
  Matrix bearings;
  bool frames[v.Rows];
  makeBearings(bearings,frames,v,sensorposes,distanceguess);
  double d=triangulateBuffer(meanpoint,path,deltapath,bearings,frames,
			     weightThreshold,sqDistanceThreshold,
			     mseThreshold, retWeight, retSqDist, retMse);
  if (d>0) {
    for (int i=0; i<v.Rows; i++)
      if (!frames[i])
        if (v(i,2)>0){
	  v(i,2)=-1;
	  d=-2;
	}
    if (secDerThreshold > 0) {
      double sd = makeDerivateCheck(meanpoint, v, sensorposes);
      if (retSecDer) {
        *retSecDer = sd;
        CureCERR(60) << "Storing secder\n";
      }
      if (sd > secDerThreshold) return -4;
    }
  }
  
  return d;
}

double 
RelativePointFeatureHelper::triangulateBuffer(Cure::Matrix &meanpoint,
				      Cure::Matrix &path, 
				      Cure::Matrix &deltapath, 
				      Cure::Matrix &bearings,
				      bool *frames,
				      const double weightThreshold, 
				      const double sqDistanceThreshold,
				      const double mseThreshold,
                                      double *retWeight,
                                      double *retSqDist,
                                      double *retMse)
{

  if (retSqDist) {
    *retSqDist = 0;
  }

  int n=bearings.Rows;
  Matrix r(n),w(n),rw(n),tw(n,1);
  double totalw=0;
  
  for (int i=0; i<n; i++)
    { 
      if (frames[i]){
	int j3=3*(i+1);
	for (int j=i+1; j<n; j++,j3+=3)
	  if (frames[j]){
	    double dot=(bearings(i,0)*bearings(j,0)+
			bearings(i,1)*bearings(j,1)+
			bearings(i,2)*bearings(j,2));
	    w(i,j)=1-dot*dot;
	    if (w(i,j)>1E-4)
	      {
		w(j,i)=w(i,j);
		tw(i,0)+=w(i,j);
		tw(j,0)+=w(j,i);
		double b1=bearings(i,0)*deltapath(i,j3)+
		  bearings(i,1)*deltapath(i,j3+1)+
		bearings(i,2)*deltapath(i,j3+2);
		double b2=-(bearings(j,0)*deltapath(i,j3)+
			    bearings(j,1)*deltapath(i,j3+1)+
			    bearings(j,2)*deltapath(i,j3+2));
		rw(i,j)=-(b1+dot*b2);
		r(i,j)=rw(i,j)/w(i,j);
		rw(j,i)=-(dot*b1+b2);
		r(j,i)=rw(j,i)/w(j,i);
	      }	 else w(i,j)=0;
	  }
	totalw+=tw(i,0);
      }
    }
  double mse=0;
  bool loop=true;
  while (loop){
    CureCERR(60) << "totalw=" << totalw 
                 << " (" << weightThreshold << ")\n";
    
    if (retWeight) {
      CureCERR(60) << "Storing weight\n";
      *retWeight = totalw;
    }

    if (totalw<weightThreshold)return -1;
    meanpoint.reallocate(3,1);
    meanpoint=0;
    int len=0;
    for (int i=0; i<n; i++)
      {
	if (frames[i])
	  if (tw(i,0)>0){
	    len++;
	    double range=0;
	    for (int j=0; j<n; j++)
	      if (frames[j])
		range+=rw(i,j);
	    if (tw(i,0)>0)
	      range/=tw(i,0);
	    double p[3];
	    if (tw(i,0)>0){
	      p[0]=range*bearings(i,0)+path(i,0);
	      p[1]=range*bearings(i,1)+path(i,1);
	      p[2]=range*bearings(i,2)+path(i,2);
	      meanpoint(0,0)+=(tw(i,0)*p[0]);
	      meanpoint(1,0)+=(tw(i,0)*p[1]);
	      meanpoint(2,0)+=(tw(i,0)*p[2]);
	    }
	  }
      }
    
    meanpoint(0,0)/=totalw;
    meanpoint(1,0)/=totalw;
    meanpoint(2,0)/=totalw;
    //Now check consistancy of whole set
    mse=0;
    loop=false;
    CureCERR(40) << "totalw=" << totalw 
          << " (" << weightThreshold << ")\n";
    for (int i=0; i<n; i++)
      {
	if (frames[i]){
	  double m[3];
	  m[0]=meanpoint(0,0)-path(i,0);
	  m[1]=meanpoint(1,0)-path(i,1);
	  m[2]=meanpoint(2,0)-path(i,2);
	  double test=(m[0]*bearings(i,0)+
		       m[1]*bearings(i,1)+
		       m[2]*bearings(i,2));
	  test*=test;
	  test=(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]-test);
	  //	  std::cerr<<test<<" ";
	  if (test>sqDistanceThreshold) {
	    frames[i]=false;
	    totalw-=tw(i,0);
	    tw(i,0)=0;
	    for (int j=0; j<n;j++)
	      if (frames[j]){
		totalw-=w(j,i);
		tw(j,0)-=w(j,i);
		w(i,j)=0;
		w(j,i)=0;
		rw(j,i)=0;
		rw(i,j)=0;
	      }
	    loop=true;
	  } 

          if (retSqDist && test > *retSqDist) {
            CureCERR(60) << "Storing sqDist\n";
            *retSqDist = test;
          }

	  mse+=test;
	}
      }
    mse/=len;
    CureCERR(40)<<" \n";
    CureCERR(40) <<"MSE "<<mse<<" "<<len<<" "<< "totalw=" << totalw 
		 << " (" << weightThreshold << ")\n"<<
      meanpoint(0,0)<<" "<<meanpoint(1,0)<<" "<<meanpoint(2,0)<<"\n";
  }
  
  if (retWeight) {
    CureCERR(60) << "Storing weight\n";
    *retWeight = totalw;
  }
  if (retMse) {
    CureCERR(60) << "Storing mse\n";
    *retMse = mse;
  }
  if (mse>mseThreshold)return -3;
  return totalw;
}

double
RelativePointFeatureHelper::makeDerivateCheck(Cure::Matrix &meanpoint,
				      Cure::Matrix &v,
                                      Cure::Transformation3D *sensorposes)
{
  Matrix vhat(v);
  predictPixels(vhat,sensorposes,meanpoint);
  Matrix innov;
  calcInnovations(innov,vhat,v,0);
  
  Matrix di;
  calcTimeDerivative(di,innov);
  Matrix d2i;
  calcTimeDerivative(d2i,di);
  double maxxi=0;
  double minxi=0;
  double maxyi=0;
  double minyi=0;

  /*
  for (int i=0;i<innov.Rows; i++) {
    if (innov(i,0)>maxxi)maxxi=innov(i,0);
    if (innov(i,0)<minxi)minxi=innov(i,0);
    if (innov(i,1)>maxyi)maxyi=innov(i,1);
    if (innov(i,1)<minyi)minyi=innov(i,1);
  } 
  CureCERR(30)<<"\nInovation: "
              <<minxi<<" "<<maxxi<<" "<<minyi<<" "<<maxyi<<"\n";

  maxxi=0;
  minxi=0;
  maxyi=0;
  minyi=0;
  for (int i=0;i<di.Rows; i++) {
    if (di(i,0)>maxxi)maxxi=di(i,0);
    if (di(i,0)<minxi)minxi=di(i,0);
    if (di(i,1)>maxyi)maxyi=di(i,1);
    if (di(i,1)<minyi)minyi=di(i,1);

  } 
  CureCERR(30)<<"di_dt: "<<minxi<<" "<<maxxi<<" "<<minyi<<" "<<maxyi<<"\n";
  */

  maxxi=0;
  minxi=0;
  maxyi=0;
  minyi=0;
  for (int i=0;i<d2i.Rows; i++) {
    if (d2i(i,0)>maxxi)maxxi=d2i(i,0);
    if (d2i(i,0)<minxi)minxi=d2i(i,0);
    if (d2i(i,1)>maxyi)maxyi=d2i(i,1);
    if (d2i(i,1)<minyi)minyi=d2i(i,1);
  } 
  CureCERR(40)<<"d^2i_dt^2: "<<minxi<<" "<<maxxi<<" "<<minyi<<" "<<maxyi<<"\n";

  // Return the biggest (absolute) second derivative
  double ret = fabs(minxi);
  if (maxxi > ret) ret = maxxi;
  if (fabs(minyi) > ret) ret = fabs(minyi);
  if (maxyi > ret) ret = maxyi;

  return ret;
}

/*
double 
RelativePointFeatureHelper::fullTriangulate(Cure::Matrix &meanpoint,
				      Cure::Matrix &path, 
				      Cure::Matrix &deltapath, 
				      Cure::Matrix &bearings,
				      bool *frames,
				      const double weightThreshold, 
				      const double sqDistanceThreshold,
				      const double mseThreshold)
{
  

  //Cost=r^T c r - 2 b^T r 
  // c= (n-1)I - matrix of bearings dot bearings
  // b= bearing dot deltapatah
  //r=c^-1 b
  for (int i=0
  Matrix c(n), b(n,1), r(n,1),w(n,1);
  w=0;
  c=(n-1);
  b=0;
  int k=0;
  double totalw=0;
  for (int i=1; i<n; i++)
    for (int j=0; j<i; j++)
      {
	double d=(bearings(i,0)*bearings(j,0)+bearings(i,1)*bearings(j,1)+
		 bearings(i,2)*bearings(j,2));
	c(i,j)-=d;
	c(j,i)=c(i,j);
	d=1-d;
	w(i,0)+=d;
	w(j,0)+=d;
	totalw+=d;
	b(j,0)+=(bearings(j,0)*deltapath(k,0)+bearings(j,1)*deltapath(k,1)+
		 bearings(j,2)*deltapath(k,2));
	b(i,0)-=(bearings(i,0)*deltapath(k,0)+bearings(i,1)*deltapath(k,1)+
		 bearings(i,2)*deltapath(k,2));
	k++;
      }
  if (totalw<threshold) {
    return -1;
  }
  totalw*=2;
  c.invert();
  r.multiply_(c,b);
  meanpoint.reallocate(3,1);
  meanpoint=0;
  for (int i=0; i<n; i++)
  {
  if (w(i,0)>1E-10)
  {
  double p[3];
  p[0]=r(i,0)*bearings(i,0)+path(i,0);
  p[1]=r(i,0)*bearings(i,1)+path(i,1);
  p[2]=r(i,0)*bearings(i,2)+path(i,2);
  meanpoint(0,0)+=(w(i,0)*p[0]);
  meanpoint(1,0)+=(w(i,0)*p[1]);
	    meanpoint(2,0)+=(w(i,0)*p[2]);
	    }
      else totalw-=w(i,0);
      }
      meanpoint(0,0)/=totalw;
      meanpoint(1,0)/=totalw;
      meanpoint(2,0)/=totalw;
      //Now check consistancy of whole set
      double mse=0;
  for (int i=0; i<n; i++)
  {
  if (frames[i]){
    double m[3];
    m[0]=meanpoint(0,0)-path(i,0);
    m[1]=meanpoint(1,0)-path(i,1);
    m[2]=meanpoint(2,0)-path(i,2);
    double test=(m[0]*bearings(i,0)+m[1]*bearings(i,1)+m[2]*bearings(i,2));
    test*=test;
    
    test=(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]-test);
    mse+=test;
    if (test>consistancytest) {
    frames[i]=false;
      //CureCERR(30) << "test=" << test
      //             << " > " << consistancytest << std::endl;        
      return -2;
      } else {
      //CureCERR(30) << "test=" << test
      //           << " < " << consistancytest << std::endl;        
      }
      }
      }
  
  mse/=n;
  CureCERR(30) <<"MSE "<<mse<<" "<<n<<" "<< "totalw=" << totalw 
		 << " (" << threshold << ")\n";
    
  if (mse>consistancytest/sqrt(n))return -3;
  return totalw/2;
*/
