// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "PointFeatureHelper.hh"
#include "SensorData.hh"
#include "CureDebug.hh"

#ifndef DEPEND
#include <sstream>
#endif

using namespace Cure;

MapFeature * PointFeatureHelper::makeMapFeature(Cure::Pose3D &sensorpose,
					Cure::Measurement &m,
					bool addToVisable)
{
  double d=0;
  if (m.MeasurementType==33) d=m.W(0,5);
  MapPointFeature *mf=makeMapPointFeature(sensorpose,m.V.Element,d);
  if (m.MeasurementType==33)
    {
      CureCERR(40) << "Forcing extend because meastype 33\n";
      mf->forceExtend();
      m.MeasurementType=32;
      addTrackedKey(m.Key,mf->Key);
    }
  if (addToVisable){
    PosedFeature *pf=makePosedFeature(mf);
    if (pf){
      pf->transform(sensorpose, sensorpose.getCovType());
      PosedFeatures.add(pf);
      Visable.add(pf);
    }
  }
  return mf;
}

bool 
PointFeatureHelper::supportsSubconfig(int sc)
{
  return (1 <= sc && sc <= 4);
}
void 
PointFeatureHelper::printConfiguration()
{
  std::cerr << "Configured PointFeatureHelper with "
	    << " TriangleThreshold=" << TempletFeature.TriangleThreshold
	    << " WeightThreshold=" << TempletFeature.WeightThreshold
	    << " DistanceThreshold=" << TempletFeature.DistanceThreshold
	    << " TrackThreshold=" << TempletFeature.TrackThreshold
	    << " Image={"
               << " bottomLeft[0]=" << ImagePlane[0]
	    << " bottomLeft[1]=" << ImagePlane[1]
	    << " topRight[0]=" << ImagePlane[2]
	    << " topRight[1]=" << ImagePlane[3]
               << " focalLength=" << PixelInfo[0]
	    << " pixwidth=" << PixelInfo[1]
	    << " pixheight=" << PixelInfo[2] << "}"
	    << " RoughSearchRange=" << RoughSearchRange
               << " NearSearchRange=" << NearSearchRange
	    << " MatchThreshold=" << MatchThreshold
	    << " DistanceGuess=" << DistanceGuess
	    << std::endl;
}
int
PointFeatureHelper::config(const std::string &arglist)
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
PointFeatureHelper::configVer1(const std::string &arglist) 
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
	  >> DistanceGuess
	  >> TempletFeature.TrackThreshold
	  >> TempletFeature.TriangleThreshold
	  >> TempletFeature.WeightThreshold
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
        >> MatchThreshold
        >> DistanceGuess
	>> TempletFeature.TrackThreshold) {
      return 0;
    }
    break;
  case 2:
    if (str >> TempletFeature.TriangleThreshold
	>> TempletFeature.WeightThreshold
	>> TempletFeature.DistanceThreshold) {
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

int PointFeatureHelper::addDistance(Match & mat,double distance)
{
  Measurement &m=*mat.Measure; 
  if (m.SensorType!=SensorData::SENSORTYPE_CAMERA)return 1;
  if ((m.MeasurementType!=32)&&(m.MeasurementType!=33))return 1;
  if ((m.W.Columns>4)&&(m.W.Rows)){
    setImage(m.W.Element);
  }

  mat.PathDistance=distance;
  Center.setXY(m.V.Element);
  mat.Thresholds.reallocate(1,2);
  mat.Thresholds(0,1)=RoughSearchRange;
  mat.Thresholds(0,0)=RoughSearchRange;
  mat.Metric.reallocate(2);
  mat.Metric(0,0)=1/(PixelInfo[1]*PixelInfo[1]);
  mat.Metric(0,1)=0;
  mat.Metric(1,0)=0;
  mat.Metric(1,1)=1/(PixelInfo[2]*PixelInfo[2]);

  return 0;
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
unsigned short  PointFeatureHelper::findMergeFeatures(PosedFeature * pfp[2],
				PosedFeatureList * searchList)
{
  PosedPointFeature *pl[2];
  if (searchList==0)searchList=&Visable;
  double thresholds[2];
  thresholds[0]=.005;
  thresholds[1]= 2*PixelInfo[1];
  double rect[4];
  for (PosedFeatureList *mlist=searchList; mlist->Next; 
       mlist=mlist->Next)
    {
      PosedFeature *pf=mlist->Element;
      if (pf)
	{
	  pl[0]=castPosedPointFeature(pf);
	  if (pl[0])
	    {  
	      rect[0]=pl[0]->PixelCenter.X[0];
	      rect[1]=pl[0]->PixelCenter.X[1];
	      rect[2]=pl[0]->PixelCenter.X[0];
	      rect[3]=pl[0]->PixelCenter.X[1];
	      rect[0]-=PixelInfo[1]*2;
	      rect[1]-=PixelInfo[2]*2;
	      rect[2]+=PixelInfo[1]*2;
	      rect[3]+=PixelInfo[2]*2;
	      for (PosedFeatureList *mlist2=mlist->Next; mlist2->Next; 
		   mlist2=mlist2->Next)
		{
		  pf=mlist2->Element;
		 if (pf)
		   if (pf->inRectangle(rect))
		     {
		       if (!(pf->roughMatch(pl[0]->Z.Element,thresholds)))
			 {
			   unsigned short res=0;
			   pl[1]=castPosedPointFeature(pf);
			   MapPointFeature *ml0=
			     getMapPointFeature(pl[0]->FeatureKey);
			   MapPointFeature *ml1=
			     getMapPointFeature(pl[1]->FeatureKey);
			   unsigned short typ2=ml1->Bdual.Columns;
			   unsigned short typ1=ml0->Bdual.Columns;
			   if ((typ2<3)&&(typ1<3))res=1;
			   else
			     {
			       double d=pl[0]->Center(1)-
				 pl[1]->Center(1);
			       if ((typ2<3)||(typ1<3)){
				 if ((d<.04)&(d>-.04))res=1;
			       }
			       else if ((d<.02)&(d>-.02))res=1;
			     }
			   if (res)
			     {
			       if (typ2==typ1)
				 {
				   MapPointFeature *ml=
				       getMapPointFeature(pl[0]->FeatureKey);
				   int vcount1=0;
				   if (ml)vcount1=ml->VCount;
				   ml=getMapPointFeature(pl[1]->FeatureKey);
				   int vcount2=0;
				   if (ml)vcount2=ml->VCount;
				   if (vcount1>vcount2)res=2;
				   else  res=1;
				 }
			       else if (typ2==0)res=2;
			       }
			   if (res>1)
			     {
			       pfp[0]=pl[1];
			       pfp[1]=pl[0];
			       res=1;
				   typ2=typ1;
			     }
			   else{
			     pfp[0]=pl[0];
			     pfp[1]=pl[1];
			   }
			   searchList->removeFeature(pl[0]);
			   return res;
			 }
		     }
		}
	    }
	}
    } 
  return 0;
}


int  PointFeatureHelper::getMergeContstraint(Matrix &a, Matrix &e,Matrix & cov,
					     double distance,
					     unsigned short typ, PosedFeature * pl[2])
{
  MapPointFeature *ml0=getMapPointFeature(pl[0]->FeatureKey);
  MapPointFeature *ml1=getMapPointFeature(pl[1]->FeatureKey); 
  if (!ml0)return 0;
  if (!ml1)return 0;

  int dim0=ml0->Bdual.Columns;
  int dim1=ml1->Bdual.Columns;
  a.reallocate(dim0,dim0+dim1);
  e.reallocate(dim0,1);
  cov.reallocate(dim0);
  cov=0;
  if (dim0==0)return 0;
  a.Columns=dim0;
  a=1;
  a.Element+=dim0;
  a.Columns=dim1;
  a=-1;
  a.Element-=dim0;
  a.Columns=dim0+dim1;

  Matrix x(a.Columns,1);
  x.Rows=dim0;
  ml0->getP(x);
  x.Element+=dim0;
  ml1->getP(x);
  x.Element-=dim0;
  x.Rows=a.Columns;
  e.multiply_(a,x);
  x.multTranspose(e,e,1);
  cov=(x(0,0)+1E-6);
  return 0;
}
/*
void 
PointFeatureHelper::calcPath(Cure::Matrix &path, 
                             Cure::Matrix &deltapath,
                             Transformation3D  *sensorposes, 
                             int n)
{
  path.reallocate(n,3);
  deltapath.reallocate(n,3*n);
  sensorposes[0].getXYZ(&path(0,0));
  int i3=3;
  for (int i=1; i<n; i++,i3+=3){
    int k=0;
    sensorposes[i].getXYZ(&path(i,0));
    for (int j=0; j<i; j++,k+=3)
      {
	deltapath(i,k)=path(i,0)-path(j,0);
	deltapath(i,k+1)=path(i,1)-path(j,1);
	deltapath(i,k+2)=path(i,2)-path(j,2);
	deltapath(j,i3)=-deltapath(i,k);
	deltapath(j,i3+1)=-deltapath(i,k+1);
	deltapath(j,i3+2)=-deltapath(i,k+2);
      }
    deltapath(i,k)=0;
    deltapath(i,k+1)=0;
    deltapath(i,k+2)=0;
  }
}

void 
PointFeatureHelper::growPath(Cure::Matrix &path, 
                             Cure::Matrix &deltapath,
                             Transformation3D  &sensorpose,int maxsize)
{
  int r=path.Rows;
  if (r==maxsize)return slidePath(path,deltapath,sensorpose);
  int n=r+1;
  path.grow(n,3);
  sensorpose.getXYZ(&path(r,0));
  deltapath.grow(n,3*n);
  int r3=3*r;
  int k=0;
  for (int j=0; j<r; j++,k+=3)
    {
      deltapath(r,k)=path(r,0)-path(j,0);
      deltapath(r,k+1)=path(r,1)-path(j,1);
      deltapath(r,k+2)=path(r,2)-path(j,2);
      deltapath(j,r3)=-deltapath(r,k);
      deltapath(j,r3+1)=-deltapath(r,k+1);
      deltapath(j,r3+2)=-deltapath(r,k+2);
    }
  deltapath(r,k)=0;
  deltapath(r,k+1)=0;
  deltapath(r,k+2)=0;
}
void 
PointFeatureHelper::slidePath(Cure::Matrix &path, 
                             Cure::Matrix &deltapath,
                             Transformation3D  &sensorpose)
{
  int n=path.Rows;
  path.deleteRow(0);
  path.reset(n,3);
  deltapath.deleteRow(0);
  deltapath.deleteColumn(0);
  deltapath.deleteColumn(0);
  deltapath.deleteColumn(0);
  deltapath.reset(n,3*n);
  int r=n-1;
  int r3=3*r;
  int k=0;
  sensorpose.getXYZ(&path(r,0));
  for (int j=0; j<r; j++,k+=3)
    {
      deltapath(r,k)=path(r,0)-path(j,0);
      deltapath(r,k+1)=path(r,1)-path(j,1);
      deltapath(r,k+2)=path(r,2)-path(j,2);
      deltapath(j,r3)=-deltapath(r,k);
      deltapath(j,r3+1)=-deltapath(r,k+1);
      deltapath(j,r3+2)=-deltapath(r,k+2);
    }
  deltapath(r,k)=0;
  deltapath(r,k+1)=0;
  deltapath(r,k+2)=0;
}
*/
void 
PointFeatureHelper::insertPoseTopath(Cure::Matrix &path, 
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
PointFeatureHelper::makeBearings(Cure::Matrix &bearings,
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
PointFeatureHelper::predictPixels(Cure::Matrix &vhat,
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
PointFeatureHelper::calcInnovations(Cure::Matrix &innovation,
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
PointFeatureHelper::calcTimeDerivative(Cure::Matrix &dx_dt,
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
PointFeatureHelper::triangulateBuffer(Cure::Matrix &meanpoint,
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
PointFeatureHelper::triangulateBuffer(Cure::Matrix &meanpoint,
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
PointFeatureHelper::triangulateBuffer(Cure::Matrix &meanpoint,
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
PointFeatureHelper::makeDerivateCheck(Cure::Matrix &meanpoint,
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
PointFeatureHelper::fullTriangulate(Cure::Matrix &meanpoint,
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
