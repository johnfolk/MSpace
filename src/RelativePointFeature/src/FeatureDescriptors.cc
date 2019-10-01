// = RCSID
//    $Id: FeatureDescriptors.cc ,v 1.1 2007/8/1 
//
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2007 John Folkesson
//    
#include "FeatureDescriptors.hh"
#include "RelativePointFeature.hh"
using namespace Cure;



FeatureDescriptor:: FeatureDescriptor(FeatureDescriptors *des,
				      double range){
  m_Descriptors=0;
  m_Type=0;
  if (range)m_Range=range;
  else m_Range=CURE_FEATURE_NOMINAL_RANGE;
  m_FeatureKey=-1;
  des->add(this);
}
FeatureDescriptor::FeatureDescriptor(double b[3], FeatureDescriptors *des,
				     double range){
  m_Bearing=b;
  m_Descriptors=0;
  m_Type=0;
  if (range)m_Range=range;
  else m_Range=CURE_FEATURE_NOMINAL_RANGE;
  m_FeatureKey=-1;
  des->add(this);
}

FeatureDescriptor::~FeatureDescriptor(){
  if (m_Descriptors)m_Descriptors->remove(this);
}

Transformation3D *FeatureDescriptor::sensorPose(){
  if (!m_Descriptors)return 0;
  return m_Descriptors->m_SensorPose;
}
unsigned short FeatureDescriptor::getCovType(){
  if (!m_Descriptors)return 0;
  return m_Descriptors->m_CovType;
}
RelativePointFeature * FeatureDescriptor::getRelativePointFeature(){
  if (!m_Descriptors)return 0;
  MapFeature *mf=getMapFeature(m_Descriptors->bank());
  return (RelativePointFeature *) 
    mf->getFeatureType(RELATIVEPOINTFEATURE_TYPE);
}


double FeatureDescriptor::checkEpipolar(FeatureDescriptor &f, 
					Matrix &g, Matrix &h,
					double  bearingRadsSqErr){
  Matrix jrs;   
  Cure::Pose3D relativePose, a,b;
  Transformation3D *pa=f.sensorPose();
  if (!pa)return -1;
  a=*pa;
  Transformation3D *pb=sensorPose();
  pb=sensorPose();
  if (!pb)return -1;
  b=*pb;
  a.setCovType(f.getCovType());
  b.setCovType(getCovType());
  relativePose.setCovType(63);
  relativePose.minusPlus_(a,b, &jrs); 
  Vector3D rotb1;
  relativePose.rotate(f.m_Bearing,rotb1);
  Vector3D protb1(m_Bearing);
  double dot=rotb1*m_Bearing;
  protb1*=(-dot);
  protb1+=rotb1;
  double numer=protb1*relativePose.Position;
  double denom=protb1*protb1;
  double range=m_Range;
  bool rconst=false;
  if (denom>1E-6)
      range=numer/denom;
  else  
    rconst=true;
  if (range<0.2)
    {
      rconst=true;
      range=0.2;
    }  
  
  double drange[6];
  Matrix dR[5];
  relativePose.getJac(dR);
  if (!rconst){
    Vector3D pxrel(m_Bearing);
    double d=pxrel*relativePose.Position;
    pxrel*=(-d);
    pxrel+=relativePose.Position;
    d=(-2*numer/(denom*denom));
    Vector3D tmp;
    tmp=pxrel;
    tmp/=denom;
    Vector3D tmp1=protb1;
    tmp1*=d;
    tmp+=tmp1;
    drange[0]=protb1(0)/denom;
    drange[1]=protb1(1)/denom;
    drange[2]=protb1(2)/denom;
    Vector3D bearing1(f.m_Bearing);
    drange[3]=((bearing1.leftMultBy(dR[0]))*tmp);
    drange[4]=((bearing1.leftMultBy(dR[1]))*tmp);
    drange[5]=((bearing1.leftMultBy(dR[2]))*tmp);
  }
  Matrix jxrel(3,6);
  Vector3D x(f.m_Bearing);
  x*=range;
  relativePose.transform(x,x);
  jxrel.offset(0,0,3,3);
  jxrel=dR[4];  //-Rot
  jxrel.reset(3,6);
  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      {
	jxrel(j,i+3)=dR[i](j,0)*x(0)+dR[i](j,1)*x(1)+dR[i](j,2)*x(2);
	if (!rconst)
	  {
	    jxrel(j,i+3)+=rotb1(j)*(drange[i+3]);
	    jxrel(j,i)+=rotb1(j)*(drange[i]);
	  }
      }
  
  double magx=sqrt(x*x);
  if (magx>1E-6){  
    x/=(magx);  
  } else {
    x=(m_Bearing);
    magx=1E-6;
  }  
  double bdx=(m_Bearing)*x;    
  double innovation=sqrt(1-bdx);
  Matrix jdx(1,3);
  if (innovation>1E-9){
    double d=(0.5/(magx*innovation));
    jdx(0,0)=(bdx*x(0)-m_Bearing(0))*d;
    jdx(0,1)=(bdx*x(1)-m_Bearing(1))*d;
    jdx(0,2)=(bdx*x(2)-m_Bearing(2))*d;
  }
  else {
    //Set the x of a bit so we can have a hessian
    x=rotb1.cross(m_Bearing);
    x*=.01;
    x+=m_Bearing;
    magx=1.000049999;
    x/=magx;
    bdx=.999950004;    
    double d=(.5/(1.000049999*0.007070785));
    jdx(0,0)=(bdx*x(0)-m_Bearing(0))*d;
    jdx(0,1)=(bdx*x(1)-m_Bearing(1))*d;
    jdx(0,2)=(bdx*x(2)-m_Bearing(2))*d;  
  } 
  Matrix jdr(1,6);
  jdr.multiply_(jdx,jxrel);
  if (rconst){
    Vector3D xrel=relativePose.Position;
    double mr=xrel*xrel;
    if (mr>1E-5){
      mr=sqrt(mr);
      xrel/=mr;
      mr=jdr(0,0)*xrel(0)+jdr(0,1)*xrel(1)+jdr(0,2)*xrel(2);
      jdr(0,0)-=mr*xrel(0);
      jdr(0,1)-=mr*xrel(1);
      jdr(0,2)-=mr*xrel(2);
    }
  } 
  Matrix j;
  j.multiply_(jdr,jrs);
  double sigma= bearingRadsSqErr*(1+range/magx)*(1+range/magx);
  double eng=(innovation*innovation)/(2*sigma);
  g.transpose(j);
  g*=(innovation/sigma);
  Matrix bb=j;
  bb/=(sigma);
  h.multTranspose_(j,bb,1); 
  return eng;
}
double FeatureDescriptor::checkEpipolar(FeatureDescriptor &f, 
					double bearingRadsSqErr){
  Cure::Transformation3D relativePose;
  Cure::Transformation3D *a=f.sensorPose();
  Cure::Transformation3D *b=sensorPose();
  if (!a)return -1;
  if (!b)return -1;
  relativePose=((a->inverse())+(*b));
  Vector3D rotb1;
  relativePose.rotate(f.m_Bearing,rotb1);
  Vector3D protb1(m_Bearing);
  double dot=rotb1*m_Bearing;
  protb1*=(-dot);
  protb1+=rotb1;
  double numer=protb1*relativePose.Position;
  double denom=protb1*protb1;
  double range=m_Range;
  if (denom>1E-6)
    range=numer/denom;
  if (range<0.2)
    range=0.2;
  Vector3D x(f.m_Bearing);
  x*=range;
  relativePose.transform(x,x);  
  double magx=sqrt(x*x);
  if (magx>1E-6){  
    x/=(magx);  
  } else {
    x=(m_Bearing);
    magx=1E-6;
  }  
  double bdx=(m_Bearing)*x;    
  double innovation=sqrt(1-bdx);
  Matrix jdx(1,3);
  if (innovation<=1E-9)
      magx=1.000049999;
  double sigma= bearingRadsSqErr*(1+range/magx)*(1+range/magx);
  return (innovation*innovation)/(2*sigma);
}

double FeatureDescriptor::getWeightedRange(double &wr,
					   FeatureDescriptor &f){
  Cure::Transformation3D relativePose;
  Cure::Transformation3D *a=f.sensorPose();
  Cure::Transformation3D *b=sensorPose();
  if (!a)return -1;
  if (!b)return -1;
  relativePose=((a->inverse())+(*b));
  Vector3D rotb1;
  relativePose.rotate(m_Bearing,rotb1);
  Vector3D protb1(f.m_Bearing);
  double dot=rotb1*f.m_Bearing;
  protb1*=(-dot);
  protb1+=rotb1;
  double numer=protb1*relativePose.Position;
  double denom=protb1*protb1;
  wr=m_Range;
  if (denom>1E-6)
    wr=numer/denom;
  if (wr<0.2)
    wr=0.2;  
  dot*=dot;
  dot=1-dot;
  wr*=dot;
  return dot;
}
double FeatureDescriptor::getWeightedTheta(double &wt,
					   FeatureDescriptor &f){
  Cure::Transformation3D relativePose;
  Cure::Transformation3D *a=f.sensorPose();
  Cure::Transformation3D *pb=sensorPose();
  if (!a)return -1;
  if (!pb)return -1;
  relativePose=((a->inverse())+(*pb));
  Vector3D rot;
  Vector3D b(f.m_Bearing);
  b*=f.m_Range;
  relativePose.transform(b,rot);
  double rho=rot(0)*rot(0)+rot(1)*rot(1);
  wt=atan2(rot(1),rot(0));
  double w=f.m_Range;
  if (w<1E-6)w=1E6;
  else w=1/w;
  double dot=rot(0)*m_Bearing(0)+rot(1)*m_Bearing(1)+rot(2)*m_Bearing(2);
  dot*=dot;
  rho+=rot(2)*rot(2);
  dot/=rho;
  w*=(1-dot);
  wt*=w;
  return w;
}
double FeatureDescriptor::getWeightedPhi(double &wt,
					   FeatureDescriptor &f){
  Cure::Transformation3D relativePose;
  Cure::Transformation3D *a=f.sensorPose();
  Cure::Transformation3D *pb=sensorPose();
  if (!a)return -1;
  if (!pb)return -1;
  relativePose=((a->inverse())+(*pb));
  Vector3D rot;
  Vector3D b(f.m_Bearing);
  b*=f.m_Range;
  relativePose.transform(b,rot);
  double rho=rot(0)*rot(0)+rot(1)*rot(1);
  wt=atan2(rot(2),sqrt(rho));
  rho+=rot(2)*rot(2);
  double w=f.m_Range;
  if (w<1E-6)w=1E6;
  else w=1/w;
  double dot=rot(0)*m_Bearing(0)+rot(1)*m_Bearing(1)+rot(2)*m_Bearing(2);
  dot*=dot;
  dot/=rho;
  w*=(1-dot);
  wt*=w;
  return w;
}

void FeatureDescriptors::clear(){
  if (m_Desc){
    for (int i=0; i< m_NumberDesc;i++)
      if (m_Desc[i]){
	MapFeature *mf=Cure::getMapFeature(bank(),m_Desc[i]->m_FeatureKey);
	if (mf)
	  if (mf->Type==RELATIVEPOINTFEATURE_TYPE){
	    RelativePointFeature *rpf=(RelativePointFeature *)mf;
	    rpf->m_Descriptors.removeDescriptor(m_Desc[i]);
	  }
	delete m_Desc[i];
	m_Desc[i]=0;
      }
    delete [] m_Desc;
  }
  m_Desc=0;
  m_NumberDesc=0;
}


int FeatureDescriptors::match(FeatureDescriptors &fd, 
			      Cure::Matrix &errors,
			      Cure::ShortMatrix &matches, 
			      float limit)
{      
  matches.reallocate(0,2);
  float matcherr[m_NumberDesc];     
  for (int i=0; i<m_NumberDesc; i++)
    {
      matcherr[i]=limit;
      float er;
      unsigned short index;
      if (fd.match(*m_Desc[i],index,er))
	{
	  if (er<limit)
	    {
	      unsigned short index2;
	      if (match(*fd(index),index2,er))
		if (index2==i){
		  //  std::cerr<<i<<" "<<er<<" ";
		  if (er<limit)
		    {
		      matches.append((short)i,0);
		      matches.append((short)index,1);
		      matcherr[matches.Rows]=er;
		    }
		}
	    }
	}
    }
  errors.reallocate(matches.Rows,1);
  for (int i=0;i<matches.Rows; i++)
    errors(i,0)=matcherr[i];
  return matches.Rows;
}

FeatureDescriptorList::FeatureDescriptorList()
{
  Element=0;
  Next=0;
}
 
FeatureDescriptorList::~FeatureDescriptorList()
{
  if (Next!=0) delete Next;
  Next=0;
}

void FeatureDescriptorList::clean()
{
  if (Next!=0) delete Next;
  Next=0;
  Element =0;
}
unsigned short FeatureDescriptorList::add(FeatureDescriptor *w)
{
  if (Next==0)
    {
      Element=w;
      Next=new FeatureDescriptorList();
      return 1;
    }
  return (1+Next->add(w));
}

unsigned short FeatureDescriptorList::addUnique(FeatureDescriptor *w)
{
  if (Element==w)return 0;
  if (Next==0)
    {
      Next=new FeatureDescriptorList();
      Element=w;
      return 1;
    }
  return Next->addUnique(w);
}
unsigned short  FeatureDescriptorList::count()
{
  if (Next)return (1+Next->count());
  return 0;
}
bool FeatureDescriptorList::remove(unsigned short k)
{
  if (k>0)
    { 
      if (Next)return Next->remove(k-1);
      return false;
    }
  if (Next){
    if (Next->Next){
      Element=Next->Element;
      Next->remove(0);
      return true;
    }    
    delete Next;
    Next=0;
    Element=0;
    return true;
  }
  return false;
}
FeatureDescriptor * FeatureDescriptorList::getRemove(unsigned short k)
{
  if (k>0)
    { 
      if (Next)return Next->getRemove(k-1);
      return 0;
    }
  if (Next){
    FeatureDescriptor *f=Element;
    remove(0);
    return f;
  }
  return 0;
}
unsigned short FeatureDescriptorList::removeDescriptor(FeatureDescriptor *n)
{
  if (Next==0)return 0;
  unsigned short r=0;
  while (Element==n){
    r++;
    remove(0);
  }
  if (Next){
    r+=Next->removeDescriptor(n);
  }
  return r;
}
bool FeatureDescriptorList::find(FeatureDescriptor *n)
{
  if (n==Element)return true;
  if (Next)return Next->find(n);
  return false;
}

FeatureDescriptor * FeatureDescriptorList::get(unsigned short k)
{
  if (k>0)
    { 
      if (Next)return Next->get(k-1);
      return 0;
    }
  return Element;
}
