//
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2007 John Folkesson
//    

#include "RelativePointFeature.hh"
using namespace Cure;

void RelativePointFeature::get(GenericData &gd)
{
  long c=gd.ShortData.Columns;
  if (c<8)c=8;
  long r=gd.ShortData.Rows;
  gd.forceShortDataSize(r+1,c);
  gd.ShortData(r,0)=m_Function;
  gd.ShortData(r,1)=m_Bad;
  gd.ShortData(r,2)=m_WeakInfo;
  gd.ShortData(r,3)=m_RefBranch;
  gd.ShortData.setLong(m_ReferenceKey,r,4);
  MapPoseTree *pt=getMapPoseTree(Bank,m_ReferenceKey);
  if (pt)gd.ShortData.setLong(pt->ID,r,6);
  else gd.ShortData.setLong(0,1,2);
  r=gd.Data.Rows;
  c=gd.Data.Columns;
  gd.Data.grow(r+2,c);
  gd.Data(r,0)=m_Scale;
  gd.Data(r+1,0)=m_Alpha;
}
int RelativePointFeature::set(GenericData &gd)
{
  if (MapFeature::set(gd))return 1;
  if (gd.ShortData.Rows<4)return 1;
  if (gd.ShortData.Columns<8)return 1;
  if (gd.Data.Rows<5)return 1;
  if (gd.Data.Columns<1)return 1;
  cleanDescriptors();
  long r=gd.ShortData.Rows-1;
  m_Function=gd.ShortData(r,0);
  m_Bad=gd.ShortData(r,1);
  m_WeakInfo=gd.ShortData(r,2);
  m_ReferenceKey=gd.ShortData.getLong(r,4);
  m_RefBranch=gd.ShortData(r,3);
  r=gd.Data.Rows-2;
  m_Scale=gd.Data(r,0);
  m_Alpha=gd.Data(r+1,0);
  short bankid;
  long key;
  if (getBankInfo(gd,bankid,key))return 1;
  MapObject *mo =0;
  if (Bank)mo=Bank->getMapObject(bankid,m_ReferenceKey);
  if (mo){
    if (mo->getObjectSubType()!=POSETREE_TYPE)return MAP_OBJECT_INVALID;
    if (mo->getObjectType()!=FRAME_TYPE)return MAP_OBJECT_INVALID;
  } else{
    mo=new MapPoseTree(Bank);
    if (Bank)Bank->associate(bankid,m_ReferenceKey,mo->Key);
  }
  m_ReferenceKey=mo->Key;
  return 0;
}
RelativePointFeature::RelativePointFeature():MapFeature(0)
{
  m_Function=0;
  m_Scale=130;
  m_Alpha=1E-10;
  m_Bad=0;
  m_WeakInfo=0x38;
  m_ReferenceKey=-1;
  m_RefBranch=0;
  FullDim=3;
  NumberScalars=3;
  Scalars=AllocatedScalars;
  Index=AllocatedIndex;
  Type=RELATIVEPOINTFEATURE_TYPE;
  init(); 
}
RelativePointFeature::RelativePointFeature(RelativePointFeature *wp,
				   MapBank *b)
  :MapFeature(b)
{
  m_Bad=0;
  FullDim=3;
  NumberScalars=3;
  Scalars=AllocatedScalars;
  Index=AllocatedIndex;
  Type=RELATIVEPOINTFEATURE_TYPE;
  init();
  if (wp)
     {
       m_WeakInfo=wp->m_WeakInfo;
       m_Alpha=wp->m_Alpha;
       m_Scale=wp->m_Scale;
       m_Function=wp->m_Function;
       CastPtr=wp->CastPtr;
       m_ReferenceKey=wp->m_ReferenceKey;
       m_RefBranch=wp->m_RefBranch;
     }
  else
    {
      m_Function=0;
      m_Scale=130;
      m_Alpha=1E-10;
      m_WeakInfo=0x38;
      m_ReferenceKey=-1;
      m_RefBranch=0;
      CastPtr=0;
    } 

}
MapFeature *RelativePointFeature::copy()
{
  RelativePointFeature *mw=new RelativePointFeature(this, Bank);
  mw->m_Bad=m_Bad;
   mw->Consecutive=Consecutive;
  mw->LastDistance=LastDistance;  
  mw->LastBearing[0]=LastBearing[0];
  mw->LastBearing[1]=LastBearing[1];
  mw->LastBearing[2]=LastBearing[2];
  for (int i=0;i<FullDim;i++)
    mw->Index[i]=Index[i];
  mw->DistanceThreshold=DistanceThreshold;
  mw->Bdual=Bdual; 
  mw->Info=Info;
  Matrix x(3,1);
  getX(x);
  mw->setX(x);
  return mw;
}
void RelativePointFeature::initializeFromPixels(const double pixels[2],
						double focallength,
						double distanceguess)
{
 double  v[3];
  double d=(distanceguess-focallength)/focallength;
  v[0]=-(pixels[0]*d);
  v[2]=-(pixels[1]*d);
  v[1]=distanceguess;
  setCenter(v);
}


int  RelativePointFeature::extend()
{
  if (Bdual.Columns==3)return 0;
  int dim=0;
  if (!(m_WeakInfo&0x8))dim++;
  if (!(m_WeakInfo&0x10))dim++;
  if (!(m_WeakInfo&0x20))dim++;
  if (dim<=Bdual.Columns)return 0;
  int dm=Bdual.Columns;
  Bdual.reallocateZero(3,dim);
  if (dim==3){
    Bdual=1;
    return 3-dm;;
  }
  if (dim==2){
    if (!(m_WeakInfo&0x18))
	Bdual=1;
    else if (!(m_WeakInfo&0x28))
      {
	Bdual(0,0)=1;
	Bdual(2,1)=1;
      }
    else if (!(m_WeakInfo&0x30))
      {
	Bdual(1,0)=1;
	Bdual(2,1)=1;
      }
    return 2-dm;
  }
  if (!(m_WeakInfo&0x08))
    Bdual=1;
  else if (!(m_WeakInfo&0x10))
    Bdual(1,0)=1;
  else if (!(m_WeakInfo&0x20))
    Bdual(2,0)=1;
  return 1-dm;;
}

int RelativePointFeature::merge(RelativePointFeature *mf, unsigned short typ)
{
  if (!mf)return -1;
  for (FeatureDescriptorList *d=&mf->m_Descriptors;d->Next;d=d->Next)
    if (d->Element){
      d->Element->m_FeatureKey=Key;
      m_Descriptors.addUnique(d->Element);
    }
  
  mf->m_Descriptors.clean();
  return 0;
}

void RelativePointFeature::getB(Cure::Matrix &b)const
{
  b.transpose(Bdual);
}
/**
 * phi range
 * 1 0
 * 0 0
 * 0 1
 *
 * phi theta
 * 1 0
 * 0 1
 * 0 0
 */
bool RelativePointFeature::getC(RelativePointFeature *mw,Cure::Matrix &c)const
{
  if ((Bdual.Columns==0)||(mw->Bdual.Columns==0))return false;
  if (m_Function!=mw->m_Function)return false;
  int dim=2;
  if ((Bdual.Columns==3)&&(mw->Bdual.Columns==3))dim++;
  c.reallocate(dim,Bdual.Columns);
  if (dim==3)
    c=1;
  else if (Bdual(1,1)==1)
    c=1;
  else {
    c=0;
    c(0,0)=1;
    c(1,c.Columns-1)=1;
  }
  return true; 
}
bool RelativePointFeature::testMatch(RelativePointFeature *mw, double tolerance)
{
  if ((Bdual.Columns!=3)||(mw->Bdual.Columns!=3))return false;
  return true;
}

int RelativePointFeature::getCartesian(Matrix &jac,Matrix &x){
  MapPoseTree *pt=getMapPoseTree(Bank,m_ReferenceKey);
  if (!pt)return 1;
  pt->calc();
  Transformation3D trans;
  unsigned short type;
  pt->getLeafPose(trans,type,m_RefBranch);
  Matrix jsref;
  pt->getLeafJacobian(jsref,m_RefBranch);
  unsigned short pdim=jsref.Columns;
  Matrix dR[5];
  Vector3D xr,b1,drb1, invrotxr, drrinvinvrotb1;
  double r1=getRange();    
  getBearing(b1);
  xr=b1;
  xr*=r1;
  trans.getJac(dR);
  
  x.reallocate(3,1); 
  trans.getXYZ(x.Element);
  trans.invRotate(xr,invrotxr);
  
  x(0,0)+=invrotxr(0);
  x(1,0)+=invrotxr(1);
  x(2,0)+=invrotxr(2);
 


  Matrix jxrs(3,jsref.Rows);//,jxrx(3)=dR[3];
  jac.reallocateZero(3,3+pdim);  
  unsigned short t=type;
  int k=0;
  for (int i=0;i<3;i++,t=(t>>1)){
    if (t&1){
      for (int j=0;j<3;j++)
	jxrs(j,k)=dR[4](j,i);
      k++;
    }
  }
  for (int i=3;i<6;i++,t=(t>>1)){
    if (t&1){
      Vector3D temp;
      temp=xr.leftMultBy(dR[i]);
      for (int j=0;j<3;j++)
	jxrs(j,k)=temp(j);
      k++;
    }
  }
  Matrix jxr(3,pdim+3);

  jxr.offset(0,0,3,pdim);
  jxr.multiply_(jxrs,jsref);
  jxr.offset(0,pdim,3,3);
  jxr=dR[3];
  jxr.offset(0,-pdim,3,3+pdim);
 
  double r2=r1*r1;
  double rho2=r2-xr(2)*xr(2);
  double rho=sqrt(rho2);
  Matrix juxr(3);
  double a=r2*rho;
  //Really one needs to make sure rho=0 and r=0 dont happen
  // so the z axis must be not in the sensor fov and r must not get too small
  if (a<1E-12)a=1E12;
  else a=1/a;
  double b=1E12;
  if (rho2>1E-12)b=1/rho2;
  double c=1E12;
  if (r2>1E-12)c=1/r2;
  double d=1E12;
  if (r1>1E-12)d=1/r1;
  
  juxr(0,0)=xr(2)*xr(0)*a;
  juxr(0,1)=xr(2)*xr(1)*a;
  juxr(0,2)=-rho*c;
  
  juxr(1,0)=-xr(1)*b;
  juxr(1,1)=xr(0)*b;
  juxr(1,2)=0;
  
  juxr(0,0)=xr(0)*d;
  juxr(0,1)=xr(1)*d;
  juxr(0,2)=xr(2)*d;
  double jpolar[3];
  int rows=getPolarJac(jpolar);
  for (int i=0;i<3;i++,rows=(rows>>1))
    if (rows&1)
      for (int j=0;j<3;j++)
	juxr(i,j)/=jpolar[i];
  
  jac.multiply_(juxr,jxr);
  return 0;
}
