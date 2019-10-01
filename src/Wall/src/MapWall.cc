// = RCSID
//    $Id: MapWall.cc ,v 1.1 2004/04/1 
//
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "MapWall.hh"
#include "GenericData.hh"
using namespace Cure;
using namespace MatrixStuff;


MapWall::MapWall():MapFeature(0)
{
  FullDim=4;
  NumberPoints=2;
  Number2D=2;
  Points=AllocatedPoints;
  Index=AllocatedIndex;
  Type=MAPWALL_TYPE;
  StartLimit=0;
  EndLimit=0;
  CanLimit=0;
  InfoInitialized=0;
  Length=0;
  Ls=0;
  Tangent[0]=0;
  Tangent[1]=0;
  init();
  CastPtr=0;
  CountThreshold=25;
  TightnessValue=.3;
  EndThreshold=.1;
  LengthThreshold=1;
  DistanceThreshold=0;
  VarRhoThreshold=1E-3;

}
MapWall::MapWall( MapWall * wp,MapBank *b):MapFeature(b)
{

  Type=MAPWALL_TYPE;
  FullDim=4;
  NumberPoints=2;
  Number2D=2;
  Points=AllocatedPoints;
  Index=AllocatedIndex;
  Type=MAPWALL_TYPE;
  StartLimit=0;
  EndLimit=0;
  CanLimit=0;
  InfoInitialized=0;
  Length=0;
  Ls=0;
  Tangent[0]=0;
  Tangent[1]=0;
  init();
  if (wp)
     {
       CastPtr=wp->CastPtr;
       CountThreshold=wp->CountThreshold;
       TightnessValue=wp->TightnessValue;
       EndThreshold=wp-> EndThreshold;
       LengthThreshold=wp->LengthThreshold;
       VarRhoThreshold=wp->VarRhoThreshold;
       DistanceThreshold=wp->DistanceThreshold;
     }
  else
    {
      CastPtr=0;
      CountThreshold=25;
      TightnessValue=.3;
      EndThreshold=.1;
      LengthThreshold=1;
      DistanceThreshold=0;
      VarRhoThreshold=1E-3;
   }
}
MapFeature * MapWall::copy(){
  MapWall *mw=new MapWall(this,Bank);
  mw->Consecutive=Consecutive;
  mw->LastDistance=LastDistance;  
  mw->LastBearing[0]=LastBearing[0];
  mw->LastBearing[1]=LastBearing[1];
  mw->LastBearing[2]=LastBearing[2];
  for (int i=0;i<FullDim;i++)
    mw->Index[i]=Index[i];
  mw->DistanceThreshold=DistanceThreshold;
  mw->Bdual=Bdual;
  mw->StartLimit=StartLimit;
  mw->EndLimit=EndLimit;
  mw->CanLimit=CanLimit;
  mw->InfoInitialized=InfoInitialized;
  mw->Length=Length;
  mw->Ls=Ls;
  mw->Tangent[0]=Tangent[0];
  mw->Tangent[1]=Tangent[1];
  mw->Info=Info;
  mw->NumberSharedPoints=NumberSharedPoints;
  if (PSharedIndex){
    mw->PSharedIndex=new short[4];
    for (int i=0;i<4;i++)
      mw->PSharedIndex[i]=PSharedIndex[i];
  }

  mw->initailize(Points[0]->X,  Points[1]->X);
  return mw;
}
void MapWall::setStart(MapPoint *val)
{
  addPoint(val,0);
}
void MapWall::setEnd(MapPoint *val)
{
  addPoint(val,1);
}

void MapWall::calcTangent()
{
   Tangent[0]=Points[1]->X[0]-Points[0]->X[0];
   Tangent[1]=Points[1]->X[1]-Points[0]->X[1];
   double d=Tangent[0]*Tangent[0]+Tangent[1]*Tangent[1];
   Length=sqrt(d);
   if (Length<1E-9)return;
   Tangent[0]/=Length;
   Tangent[1]/=Length;
}

void MapWall::calcNormal(double x[2])
{
  calcTangent();
   x[0]=Tangent[1];
   x[1]=-Tangent[0];
   return;
}
void MapWall::initailize(double startx[3], double endx[3], 
			 unsigned short forceextend)
{
  if (!Points[0])
    {
      MapPoint *p=new MapPoint(Bank);
      setStart(p);
    }
  if (!Points[1])
    {
      MapPoint *p=new MapPoint(Bank);
      setEnd(p);
    }
 *Points[0]=(startx);
 *Points[1]=endx;
 CanLimit=(forceextend&6);
 InfoInitialized=(forceextend&1); 
 calcTangent();
 extend();
 recenter();
}

void MapWall::initailizeFromLine(Transformation3D &t,Line2D & ln)
{
  if (!Points[0])
    {
      MapPoint *p=new MapPoint(Bank);
      setStart(p);
    }
  if (!Points[1])
    {
      MapPoint *p=new MapPoint(Bank);
      setEnd(p);
    }
  t.invTransform2D(ln.StartPoint.X, Points[0]->X);
  t.invTransform2D(ln.EndPoint.X, Points[1]->X);
  recenter();
}
int  MapWall::extend()
{
  if (Bdual.Columns>3)return 0;
  if (!InfoInitialized)return 0;
  int res=Bdual.Columns;
  if (res<2)res=2;
  if ((!StartLimit))
    {
      if (CanLimit&2)
	{
	  if (EndLimit)	StartLimit=EndLimit+1;
	  else StartLimit=2;
	  res++;
	}
    }
  if ((!EndLimit))
    {
      if (CanLimit&4)
	{
	  if (StartLimit)EndLimit=StartLimit+1;
	  else EndLimit=2;
	  res++;
	}
    }
  res-=Bdual.Columns;
  Bdual.reallocate(4,Bdual.Columns+res);
  recenter();
  return res;
}
/*
void MapWall::lengthAdjust()
{
  if (Bdual.Columns==0)return;
  if (NumberSharedPoints==0)
    {
      double d=Bdual(0,0)*Bdual(0,0)+Bdual(1,0)*Bdual(1,0);
      d*=2;
      d=Length/sqrt(d);
      Bdual(0,0)*=d;
      Bdual(1,0)*=d;
      Bdual(3,0)=-Bdual(0,0);
      Bdual(4,0)=-Bdual(1,0);
      return;
    }
  if (NumberSharedPoints==1)
    {
      double d=Bdual(0,2)*Bdual(0,2)+Bdual(1,2)*Bdual(1,2);
      d+=Bdual(2,2)*Bdual(2,2)+Bdual(3,2)*Bdual(3,2);
      d=Length/(sqrt(d)*Ls); 
      Bdual(0,2)*=d;
      Bdual(1,2)*=d;
      Bdual(3,2)*=d;
      Bdual(4,2)*=d;
    }
  
}
*/
unsigned short MapWall::getMeasurementType(unsigned short type)
{
  if (Bdual.Columns==0) return 0;
  if (!StartLimit)type=(type&5);
  if (!EndLimit)type=(type&3);
  return type;
}
void MapWall::initializeToCloud(Point2DCloud & c, Transformation3D &map2info)
{
  Cloud.Rho=c.Rho;
  Cloud.Gamma=c.Gamma;
  Cloud.Orientation=c.Orientation;
  double n[2];
  n[0]=Cloud.Normal[0]*c.Normal[0]+Cloud.Normal[1]*c.Normal[1];
  n[1]=Cloud.Normal[0]*c.Normal[1]-Cloud.Normal[1]*c.Normal[0];
  Cloud.Normal[0]=n[0];
  Cloud.Normal[1]=n[1];
  Cloud.lineFrame();
  Cloud.lineOrder();
  Cloud.Normal[0]=c.Normal[0];
  Cloud.Normal[1]=c.Normal[1];
  Matrix x(4,1);
  c.getEndpoints(x);
  map2info.invTransform2D(x.Element,x.Element);
  map2info.invTransform2D(x.Element+2,x.Element+2);
  setX(x);
  if (c.SigmaRho<1E-6)c.SigmaRho=1E-6;
  double a=(1/c.SigmaRho);
  Info=(.01/c.SigmaRho);
  InfoInitialized=1;
  // return;
  calcNormal(n);
  
  Info(0,0)+=(a*n[0]*n[0]);
  Info(0,1)+=(a*n[1]*n[0]);
  Info(1,1)+=(a*n[1]*n[1]);
  a=1/(TightnessValue*TightnessValue);
  Info(0,0)+=(a*n[1]*n[1]);
  Info(0,1)-=(a*n[1]*n[0]);
  Info(0,1)+=(a*n[1]*n[1]);  
  Info(2,2)=Info(0,0);
  Info(2,3)=Info(0,1);
  Info(3,3)=Info(1,1);
  Info(1,0)=Info(0,1);
  Info(3,2)=Info(2,3);
}
int  MapWall::addInfo(const Cure::Matrix & v,Transformation3D & map2info,
		      const int type)
{
  if (Bdual.Columns==4)return 0;
  if (v.Rows==0)return 0;
  Point2DCloud c;
  if ( LastDistance<=v(0,0))
    {
      LastDistance=v(0,0);
      double mindistance=v(0,0)-DistanceThreshold;
      if (!Cloud.LineFramed)
	{
	  Cloud=v;
	  Cloud.fitLine();
	  Cloud.lineFrame();
	}
      else
	{
	  c=v;
	  Cloud.merge(c,mindistance);
	}
     if (Bdual.Columns<2)
	{
	  Cloud.fitLine();
	  int start=Cloud.getLine(c, 0, TightnessValue,CountThreshold,
				  LengthThreshold,VarRhoThreshold);
	  if (start)
	    {
	      int rem=start-c.Cloud.Rows;
	      Point2DCloud c2;
	      int res=Cloud.getLine(c2, start, TightnessValue,CountThreshold,
				    LengthThreshold,VarRhoThreshold);
	      while (res)
		{
		  start=res;
		  if(c2.Cloud.Rows>c.Cloud.Rows)
		    {
		      c=c2;
		      rem=start-c.Cloud.Rows;
		    }
		  res=Cloud.getLine(c2, start, TightnessValue,CountThreshold,
				    LengthThreshold,VarRhoThreshold);
		}
	      Cloud.remove(rem, c.Cloud.Rows);
	      initializeToCloud(c,map2info);
	      return 0;
	    }
	  return 0;
	}
      Matrix x(4,1);
      getX(x);
      map2info.transform2D(x.Element,x.Element);
      map2info.transform2D(x.Element+2,x.Element+2);
      int res=0;
      double t[2];
      t[0]=x(1,0)*Cloud.Normal[0]-x(0,0)*Cloud.Normal[1];
      t[1]=x(3,0)*Cloud.Normal[0]-x(2,0)*Cloud.Normal[1];
      double length=t[1]-t[0];
      int startrow=Cloud.findRow(t[0]);
      int endrow=Cloud.findRow(t[1])-1;
      int extends=Cloud.extend(startrow,TightnessValue,-1);
      int extende=Cloud.extend(endrow,TightnessValue,1);
      if (extends!=-1)
	{
	  if (!StartLimit)
	    {
	      double d=Cloud.Cloud(extends,1)-t[0];
	      if (d<0)
		{
		  res=1;
		  x(0,0)-=d*Cloud.Normal[1];
		  x(1,0)+=d*Cloud.Normal[0];
		  length-=d;
		}
	    }
	}
      if (extende!=-1)
	{
	  if (!EndLimit)
	    {
	      double d=Cloud.Cloud(extende,1)-t[1];
	      if (d>0)
		{
		  res=1;
		  x(2,0)-=d*Cloud.Normal[1];
		  x(3,0)+=d*Cloud.Normal[0];
		  length+=d;
		}
	    }
	}
   
      if ((startrow>0)&&(endrow<Cloud.Cloud.Rows-1))
	{
	  extends=startrow+1;
	  if (extends<Cloud.Cloud.Rows-1)
	    {
	      extende=endrow-startrow-2;
	      if (extende>0)
		Cloud.remove(extends,extende);
	    }
	}
   
      if (type&2)
	{
	  double d=x(0,0)-v(0,1);
	  d*=d;
	  d+=(x(1,0)-v(0,2))*(x(1,0)-v(0,2));
	  if (d<EndThreshold) {
	    CanLimit=(CanLimit|2);
	    double n[2];
	    calcNormal(n);
	    if (c.SigmaRho<1E-6)c.SigmaRho=1E-6;
	    double a=(1/c.SigmaRho);
	    Info(0,0)+=(a*n[0]*n[0]);
	    Info(0,1)+=(a*n[1]*n[0]);
	    Info(1,1)+=(a*n[1]*n[1]);
	    d+=(TightnessValue*TightnessValue);
	    if (d<1E-6)d=1E-6;
	    a=1/d;
	    Info(0,0)+=(a*n[1]*n[1]);
	    Info(0,1)-=(a*n[1]*n[0]);
	    Info(0,1)+=(a*n[1]*n[1]);  
	    Info(2,2)=Info(0,0);
	    Info(2,3)=Info(0,1);
	    Info(3,3)=Info(1,1);
	    Info(1,0)=Info(0,1);
	    Info(3,2)=Info(2,3);
	  }
	}
      if (type&4)
	{
	  double d=x(2,0)-v(v.Rows-1,1);
	  d*=d;
	  d+=(x(3,0)-v(v.Rows-1,2))*(x(3,0)-v(v.Rows-1,2));
	  if (d<EndThreshold){
	    CanLimit=(CanLimit|4);
	    double n[2];
	    calcNormal(n);
	    if (c.SigmaRho<1E-6)c.SigmaRho=1E-6;
	    double a=(1/c.SigmaRho);
	    Info(0,0)+=(a*n[0]*n[0]);
	    Info(0,1)+=(a*n[1]*n[0]);
	    Info(1,1)+=(a*n[1]*n[1]);
	    d+=(TightnessValue*TightnessValue);
	    if (d<1E-6)d=1E-6;
	    a=1/d;
	    Info(0,0)+=(a*n[1]*n[1]);
	    Info(0,1)-=(a*n[1]*n[0]);
	    Info(0,1)+=(a*n[1]*n[1]);  
	    Info(2,2)=Info(0,0);
	    Info(2,3)=Info(0,1);
	    Info(3,3)=Info(1,1);
	    Info(1,0)=Info(0,1);
	    Info(3,2)=Info(2,3);
	  }
	}
      if (res)
	{
	  map2info.invTransform2D(x.Element,x.Element);
	  map2info.invTransform2D(x.Element+2,x.Element+2);
	  calcTangent();
	  double dx=(x(0,0)-Points[0]->X[0])*Tangent[0]+
	    (x(1,0)-Points[0]->X[1])*Tangent[1];
	  x(0,0)=Points[0]->X[0]+dx*Tangent[0];
	  x(1,0)=Points[0]->X[1]+dx*Tangent[1];
	  dx=(x(2,0)-Points[1]->X[0])*Tangent[0]+
	    (x(3,0)-Points[1]->X[1])*Tangent[1];
	  x(2,0)=Points[1]->X[0]+dx*Tangent[0];
	  x(3,0)=Points[1]->X[1]+dx*Tangent[1];
	  setX(x);
	}   
    }
  else
    {
      c=v;
      Cloud.replace(c);
    }
  return 0;
}
int MapWall::merge(MapWall *mf, unsigned short typ)
{
  if (!mf)return -1;
  calcTangent();
  if (typ&1)
    {
      if(Bdual.Columns==4)return -1;
      Matrix x(4,1);
      getX(x);
      Matrix nx(4,1);
      mf->getX(nx);
      if(!StartLimit)
	{
	  double t;
	  t=-(x(0,0)*Tangent[0]+x(1,0)*Tangent[1]);
	  t+=(nx(0,0)*Tangent[0]+nx(1,0)*Tangent[1]);
	  if (t<0){
	    x(0,0)+=t*Tangent[0];
	    x(1,0)+=t*Tangent[1];
	  }
	  if (mf->StartLimit){
	    CanLimit=(CanLimit|2);
	    Info=1;
	  }
	}
      if(!EndLimit)
	{
	  double t;
	  t=-(x(2,0)*Tangent[0]+x(3,0)*Tangent[1]);
	  t+=(nx(2,0)*Tangent[0]+nx(3,0)*Tangent[1]);
	  if (t>0){
	    x(2,0)+=t*Tangent[0];
	x(3,0)+=t*Tangent[1];
	  }
	  if (mf->EndLimit){
	    CanLimit=(CanLimit|4);
	    Info=1;
	  }
	  
	}
      setX(x);
      return 0;
    }
  else if (typ&0x1000)
    {
      int which2=((typ>>4)&0xF);
      int which=((typ>>8)&0xF);
      Matrix x(4,1);
      getX(x);
      Matrix nx(4,1);
      mf->getX(nx);
      if(!which)
	{
	  double t;
	  
	  t=-(x(0,0)*Tangent[0]+x(1,0)*Tangent[1]);
	  if (which2==1)
	    t+=(nx(2,0)*Tangent[0]+nx(3,0)*Tangent[1]);
	  else
	    t+=(nx(0,0)*Tangent[0]+nx(1,0)*Tangent[1]);
	  x(0,0)+=t*Tangent[0];
	  x(1,0)+=t*Tangent[1];
	}
      else
	{
	  double t;
	  t=-(x(2,0)*Tangent[0]+x(3,0)*Tangent[1]);
	  if (which2==1)
	    t+=(nx(2,0)*Tangent[0]+nx(3,0)*Tangent[1]);
	  else
	    t+=(nx(0,0)*Tangent[0]+nx(1,0)*Tangent[1]);
	  x(2,0)+=t*Tangent[0];
	  x(3,0)+=t*Tangent[1];
	}
      setX(x);
      MapPoint *f=mf->removePoint(which2);
      if (f)
	if (f->numberOfFeatures()==0)
	  {
	    delete f;
	  }
      mf->addPoint(Points[which],which2);
      mf->setNumberSharedPoints();
      setNumberSharedPoints();
      mf->calcTangent();
      Ls=Length;
      mf->Ls=mf->Length;
      if (NumberSharedPoints==0)
	{
	  std::cerr<<"*MApWall NUMBER SHARED ERROR";
	}
      if (mf->NumberSharedPoints==0)
	{
	  std::cerr<<"*MApWall NUMBER SHARED ERROR";
	}
      if (NumberSharedPoints==1)
	{
	  PSharedIndex[0]=which;
	  PSharedIndex[1]=which;
	}
      else
	{
	  PSharedIndex[2]=which;
	  PSharedIndex[3]=which;
	}
      int badindex=mf->Index[0];
      if (mf->NumberSharedPoints==1)
	{
	  mf->PSharedIndex[0]=which2;
	  mf->PSharedIndex[1]=which2;
	  int k=2;
	  if (NumberSharedPoints==1)k=0;
	  mf->setIndex(Index[k],0);
	  mf->setIndex(Index[k+1],1);
	}
      else
	{
	  mf->PSharedIndex[2]=which2;
	  mf->PSharedIndex[3]=which2;
	  int k=2;
	  if (NumberSharedPoints==1)k=0;
	  badindex=mf->Index[2];
	  mf->setIndex(Index[k],2);
	  mf->setIndex(Index[k+1],3);
	}
      recenter();
      mf->recenter();
      return badindex;
    }
  return -1;
}
void MapWall::forceExtend(int typ)
 {
  Info=1;
  if (typ==0)
    CanLimit=(CanLimit|2);
  else
    CanLimit=(CanLimit|4);
 }
bool MapWall::getC(MapWall *mw,Cure::Matrix &c)const
{
  if (Bdual.Columns<2)return false;
  int dim=2;
  if ((mw->StartLimit)&&(StartLimit))
    dim++;
  if ((mw->EndLimit)&&(EndLimit))
    dim++;
  c.reallocate(dim,Bdual.Columns);
  c=0;
  double t=Tangent[0]*mw->Tangent[0]+Tangent[1]*mw->Tangent[1];
  if (t<0)return false;
  if (NumberSharedPoints==0){
    c(0,0)=1;
    c(1,1)=1;  
    int k=2;
    int j=2;
    if ((mw->StartLimit)&&(StartLimit))
      {
	c(k,j)=1;
	k++;
	j++;
      }
    else if (StartLimit)j++;
    if ((mw->EndLimit)&&(EndLimit))
	c(k,j)=1;
  }
  else if (NumberSharedPoints==1)
    {
      int k=1;
      if (((PSharedIndex[0]==0)&&(mw->StartLimit))||
	  ((PSharedIndex[0]==1)&&(mw->EndLimit))){
	c(0,0)=1;
	c(1,1)=1;
	k++;
      }
      else {
	c(0,0)=-Tangent[1];
	c(0,1)=Tangent[0];
      }
      c(k,2)=1;
      if (c.Rows>k+1)c(k+1,3)=1;
    }
  else if (NumberSharedPoints==2)
    {
      int k=1;
      if (((PSharedIndex[0]==0)&&(mw->StartLimit))||
	  ((PSharedIndex[0]==1)&&(mw->EndLimit))){
	c(0,0)=1;
	c(1,1)=1;
	k++;
      }
      else {
	c(0,0)=-Tangent[1];
	c(0,1)=Tangent[0];
      }
      if (((PSharedIndex[1]==0)&&(mw->StartLimit))||
	  ((PSharedIndex[1]==1)&&(mw->EndLimit))){
	c(k,2)=1;
      if (c.Rows>k+1)c(k+1,3)=1;
      }
      else {
	c(k,2)=-Tangent[1];
	c(k,3)=Tangent[0];
      }
    }
  return true; 
}
/*
bool MapWall::testMatch(MapWall *mw, double tolerance)
{
  Matrix x1,x2;
  getX(x1);
  mw->getX(x2);
  Matrix dx=x2;
  dx-=x1;
  if ((Bdual.Columns<2)||(mw->Bdual.Columns<2)){
    double d=Tangent[0]*mw->Tangent[1]-Tangent[0]*mw->Tangent[1];
    if (d<0)d=-d;;
    if (Length>mw->Length)d*=mw->Length;
    else d*=mw->Length;
    if (tolerance<d)  return false;
   d=(dx(0,0)+dx(2,0))*Tangent[1]-(dx(1,0)+dx(3,0))*Tangent[0];
   if (d<0)d=-1;
   d/=2;
   if (tolerance<d)  return false;
  }
  double d=dx(0,0)*Tangent[0]+dx(1,0)*Tangent[1];
  if ((!StartLimit)||(!mw->StartLimit)) 
    {
    if (d>(Length+tolerance))return false;
    if (-d>(mw->Length+tolerance))return false;
    } 
 if ((StartLimit)&&(!mw->StartLimit))
   if (d<-tolerance)return false;
 if ((mw->StartLimit)&&(!StartLimit))
   if (d>tolerance)return false;
  d=dx(2,0)*Tangent[0]+dx(3,0)*Tangent[1];
  if ((!EndLimit)||(!mw->EndLimit)){ 
    if (-d>(Length+tolerance))return false;
    if (d>(mw->Length+tolerance))return false;
  }
  if ((EndLimit)&&(!mw->EndLimit))
    if (d>tolerance)return false;
  if ((mw->EndLimit)&&(!EndLimit))
    if (d<-tolerance)return false;
  return true;  
}

*/
bool MapWall::testMatch(MapWall *mw, double tolerance)
{
  //  return false;
  Matrix x1,x2;
  getX(x1);
  mw->getX(x2);
  Matrix dx=x2;
  dx-=x1;
  double d=Tangent[0]*mw->Tangent[0]+Tangent[1]*mw->Tangent[1];
  if (d<.5)return false; 
  d=Tangent[0]*mw->Tangent[1]-Tangent[1]*mw->Tangent[0];
  if (d<0)d=-d;;
  if (Length>mw->Length)d*=mw->Length;
  else d*=Length;
  if (tolerance<d)  return false;
  d=(dx(0,0)+dx(2,0))*Tangent[1]-(dx(1,0)+dx(3,0))*Tangent[0];
  if (d<0)d=-1;
  d/=2;
  if (tolerance<d)  return false;
  d=dx(0,0)*Tangent[0]+dx(1,0)*Tangent[1];
  if ((!StartLimit)||(!mw->StartLimit)) 
    {
    if (d>(Length+tolerance))return false;
    if (-d>(mw->Length+tolerance))return false;
    } 
 if ((StartLimit)&&(!mw->StartLimit))
   if (d<-tolerance)return false;
 if ((mw->StartLimit)&&(!StartLimit))
   if (d>tolerance)return false;
  d=dx(2,0)*Tangent[0]+dx(3,0)*Tangent[1];
  if ((!EndLimit)||(!mw->EndLimit)){ 
    if (-d>(Length+tolerance))return false;
    if (d>(mw->Length+tolerance))return false;
  }
  if ((EndLimit)&&(!mw->EndLimit))
    if (d>tolerance)return false;
  if ((mw->EndLimit)&&(!EndLimit))
    if (d<-tolerance)return false;
  return true;  
}

double MapWall::judgeMatch(MapWall *mw, double tolerance)
{
  //  if (testMatch(mw,tolerance))return 1;
  Matrix x1,x2;
  getX(x1);
  mw->getX(x2);
  Matrix dx=x2;
  dx-=x1;
  double d=Tangent[0]*mw->Tangent[0]+Tangent[1]*mw->Tangent[1];
  if (d<.5)return false; 
  d=Tangent[0]*mw->Tangent[1]-Tangent[1]*mw->Tangent[0];
  if (d<0)d=-d;
  if (Length>mw->Length)d*=mw->Length;
  else d*=Length;
  if (tolerance<d)  return 0;
  d=(dx(0,0)+dx(2,0))*Tangent[1]-(dx(1,0)+dx(3,0))*Tangent[0];
  if (d<0)d=-1;
  d/=2;
  if (tolerance<d)  return 0;
  d=dx(0,0)*Tangent[0]+dx(1,0)*Tangent[1];
  double overlap=0;
  if (d>0){
    if (!(d>Length))
      overlap=Length-d;
  }
  else{
    if (!(-d>mw->Length))
      overlap=mw->Length+d;
  }
  if (overlap==0)return 0;
  if ((StartLimit)&&(!mw->StartLimit))
    if (d<-tolerance)return 0;
  if ((mw->StartLimit)&&(!StartLimit))
    if (d>tolerance)return 0;
  d=dx(2,0)*Tangent[0]+dx(3,0)*Tangent[1];
  if ((EndLimit)&&(!mw->EndLimit))
    if (d>tolerance)return 0;
  if ((mw->EndLimit)&&(!EndLimit))
    if (d<-tolerance)return 0;
  if (Length<mw->Length)d=Length;
  else d=mw->Length;
  if (d>overlap)return overlap/d;
  return 1;  
}


void MapWall::recenter()
{
  if (Bdual.Columns==0) return;
  calcTangent();

  // Assuming that the tangent direction is alpha the normal angle
  // gamma = (alpha - pi/2)
  // cos(gamma) = cos(alpha-pi/2) =-sin(alpha) =-Tangent[1]
  // sin(gamma) = sin(alpha-pi/2) = cos(alpha) = Tangent[0]

  if (NumberSharedPoints==0)
    {
      Bdual(0,0)=-Tangent[1]*M_SQRT1_2; // cos(norm_ang) / sqrt(2)
      Bdual(1,0)=Tangent[0]*M_SQRT1_2;  // sin(norm_ang) / sqrt(2)
      if (Bdual.Columns>1)
	{
	  Bdual(0,1)=Bdual(0,0);
	  Bdual(1,1)=Bdual(1,0);
	  Bdual(2,1)=Bdual(0,0);
	  Bdual(3,1)=Bdual(1,0);
	}
      int k=2;
      if (Bdual.Columns>k)
	if (StartLimit)
	  {
	    Bdual(0,StartLimit)=Tangent[0];
	    Bdual(1,StartLimit)=Tangent[1];
	    Bdual(2,StartLimit)=0;
	    Bdual(3,StartLimit)=0;
	    k++;
	  }
      if (Bdual.Columns>k)
	if (EndLimit)
	  {
	    
	    Bdual(2,EndLimit)=Tangent[0];
	    Bdual(3,EndLimit)=Tangent[1];
	    Bdual(0,EndLimit)=0;
	    Bdual(1,EndLimit)=0;
	  }
      Bdual(0,0)*=Length;
      Bdual(1,0)*=Length;
      Bdual(2,0)=-Bdual(0,0);
      Bdual(3,0)=-Bdual(1,0);
    }
  else if (NumberSharedPoints==1)
    {
      Bdual=0;

      // The shared point has the identity matrix for B
      Bdual(PSharedIndex[0]*2,0)=1;
      Bdual(PSharedIndex[0]*2+1,1)=1;

      // Index of the non-shared point
      int k=0;
      if (PSharedIndex[0]==0) k=2;

      // Set the part of the B-matrix from the unconstrained point
      Bdual(k,2)=-Tangent[1]*Length/Ls;
      Bdual(k+1,2)=Tangent[0]*Length/Ls;

      if (Bdual.Columns>3){
	Bdual(k,3)=Tangent[0];
	Bdual(k+1,3)=Tangent[1];
      }
    }
  else if (NumberSharedPoints==2)
    {
      Bdual=0;
      Bdual(PSharedIndex[0]*2,0)=1;
      Bdual(PSharedIndex[0]*2+1,1)=1;
      Bdual(PSharedIndex[2]*2,2)=1;
      Bdual(PSharedIndex[2]*2+1,3)=1;
     }
}
unsigned short MapWall::getScales(Cure::Matrix &scales)
{
  if (Bdual.Rows==0){
    scales.reallocate(0,FullDim);
    return 0;
  }
  calcTangent();
  if (NumberSharedPoints==0)
    {
      scales.reallocate(1,FullDim);
      scales(0,0)=-Tangent[0];//*M_SQRT1_2;
      scales(0,1)=-Tangent[1];//*M_SQRT1_2;
      scales(0,2)=-scales(0,0);
      scales(0,3)=-scales(0,1);
      return 1;
    }
  else if (NumberSharedPoints==1)
    {
      scales.reallocate(1,FullDim);
      scales(0,0)=-Tangent[0];//*M_SQRT1_2;
      scales(0,1)=-Tangent[1];//*M_SQRT1_2;
      scales(0,2)=-scales(0,0);
      scales(0,3)=-scales(0,1);
      scales/=Ls;
      return 4;
    }
  else  scales.reallocate(0, FullDim);
  return 0;
}
unsigned short MapWall::getFramed(){
  if (Bdual.Columns==0)return 0;
  if (Bdual.Columns==1)return 0;
  if (Bdual.Columns==2)return 3;
  if (NumberSharedPoints==0)
    {
      if (Bdual.Columns==2)return 7;
      if (Bdual.Columns==3)return 15;
    }
  else if (NumberSharedPoints==1)
    {
      unsigned short t=4;
      if (Bdual.Columns>3)
	t=12;
      return t;
    }
  return 0;
    
}
int  MapWall::getMtypeProjection(Cure::Matrix &project, 
				  unsigned short mtype)
{
  project.reallocate(Bdual.Columns);
  if (mtype==0)
    {
      project=0;
      return 1;
    }  
  project=1;
  if (Bdual.Columns<3)return 0;
  if (mtype&6)return 0;
  if (mtype==1)
    { 
      if (NumberSharedPoints==0)
	{
	  if (project.Rows<3)return 0;
	  for (int i=2;i<project.Rows; i++)
	    project(i,i)=0;
	}
      else if (NumberSharedPoints==1)
	{
	  project(0,0)=Tangent[1]*Tangent[1];
	  project(0,1)=-Tangent[1]*Tangent[0];
	  project(1,1)=Tangent[0]*Tangent[0];
	  project(1,0)=project(0,1);
	  if (Bdual.Columns>3){
	    project(3,3)=0;
	  }
	}
      else if (NumberSharedPoints==2)
	{
	  project(0,0)=Tangent[1]*Tangent[1];
	  project(0,1)=-Tangent[1]*Tangent[0];
	  project(1,1)=Tangent[0]*Tangent[0];
	  project(1,0)=project(0,1);
	  project(2,2)=project(0,0);
	  project(2,3)=project(0,1);
	  project(3,3)=project(1,1);
	  project(3,2)=project(2,3);
	}
   }
  else if (mtype&2)
    {
      if (NumberSharedPoints==0)
	{
	  if (EndLimit)
	    project(EndLimit,EndLimit)=0;
	  else return 0;
	}
      else if (NumberSharedPoints==1)
	{
	  if ((PSharedIndex[0]==0)&&(EndLimit))
	    project(3,3)=0;
	  else if (PSharedIndex[0]==1)
	    {
	      project(0,0)=Tangent[1]*Tangent[1];
	      project(0,1)=-Tangent[1]*Tangent[0];
	      project(1,1)=Tangent[0]*Tangent[0];
	      project(1,0)=project(0,1);
	    }
	}
      else if (NumberSharedPoints==2)
	{
	  int k=0;
	  if (PSharedIndex[0]==0)k=2;
	  project(k,k)=Tangent[1]*Tangent[1];
	  project(k,k+1)=-Tangent[1]*Tangent[0];
	  project(k+1,k+1)=Tangent[0]*Tangent[0];
	  project(k+1,k)=project(k,k+1);
	}

    }
  else if (mtype&4)
    {
      if (NumberSharedPoints==0)
	{
	  if (StartLimit)
	    project(StartLimit,StartLimit)=0;
	  else return 0;
	}
      else if (NumberSharedPoints==1)
	{
	  if ((PSharedIndex[0]==1)&&(StartLimit))
	    project(3,3)=0;
	  else if (PSharedIndex[0]==0)
	    {
	      project(0,0)=Tangent[1]*Tangent[1];
	      project(0,1)=-Tangent[1]*Tangent[0];
	      project(1,1)=Tangent[0]*Tangent[0];
	      project(1,0)=project(0,1);
	    }
	}
      else if (NumberSharedPoints==2)
	{
	  int k=0;
	  if (PSharedIndex[0]==1)k=2;
	  project(k,k)=Tangent[1]*Tangent[1];
	  project(k,k+1)=-Tangent[1]*Tangent[0];
	  project(k+1,k+1)=Tangent[0]*Tangent[0];
	  project(k+1,k)=project(k,k+1);
	}
    }
  return 1;
}

void MapWall::getB(Cure::Matrix &b)const
{
  b.transpose_(Bdual);
  if (b.Rows==0)return;
  if (NumberSharedPoints==0)
    {
      double d=Bdual(0,0)*Bdual(0,0)+Bdual(1,0)*Bdual(1,0);
      d*=2;
      b(0,0)/=d;
      b(0,1)/=d;
      b(0,2)/=d;
      b(0,3)/=d;
    }
  else if (NumberSharedPoints==1)
    {
      double d=Bdual(0,2)*Bdual(0,2)+Bdual(1,2)*Bdual(1,2);
      d+=Bdual(2,2)*Bdual(2,2)+Bdual(3,2)*Bdual(3,2);
      b(2,0)/=d;
      b(2,1)/=d;
      b(2,2)/=d;
      b(2,3)/=d;
    }
}
void MapWall::print(int level)
{
  Matrix x;
  getX(x);
  std::cerr<<" MapWall: "<<Key<<" "<<ID<<" "<<Bdual.Columns<<" "<<Cloud.Cloud.Rows<<" "<<NumberSharedPoints<<" "<<StartLimit<<" "<<EndLimit<<" ";
  for (int i=0; i<FullDim;i++)std::cerr<<x(i,0)<<" ";
  std::cerr<<std::endl;
  if (level>0) Bdual.print();
  if (level>1) 
    if (Index)
      std::cerr<<Index[0]<<" "<<Index[1]<<" "<<Index[2]<<" "<<Index[3]<<"\n";
}
void MapWall::write(std::fstream &fs )
{
  MapFeature::write(fs);
  fs<<StartLimit<<" "<<EndLimit<<" "<<Ls<<" "
    <<InfoInitialized<<" "<<CanLimit<<"\n";
  if (NumberSharedPoints)
    for (int i=0;i<4;i++)
      fs<<PSharedIndex[i]<<" ";
  fs<<"\n";
}
int MapWall::read(int version, std::fstream &fs, bool readKey)
{
  if (version!=1)return 1;
  int ret=MapFeature::read(1,fs,readKey);
  if (ret)return ret;
  fs>>StartLimit>>EndLimit>>Ls>>InfoInitialized>>CanLimit;
  if (NumberSharedPoints)
    {   
      if (PSharedIndex==0)
	{
	  PSharedIndex=new short[4];
	}
      for (int i=0;i<4;i++)
	fs>>PSharedIndex[i];
    }
  recenter();
  return 0;

}
MapWall * MapWall::read(std::fstream &fs, bool readKey, int version)
{
  if (version!=1)return 0;
  long key=0;
  if (readKey){
    if (!(fs>>key))return 0;
    if (key!=Key){
      if (Bank){
	MapObject *m=Bank->getMapObject(key);
	if (m){
	  std::cerr<<"MapObject with this key existes "<<key<<"\n";
	}
	Bank->remove(this);
	Bank->add(this,key);
	if (key!=Key){
	  std::cerr<<"wanted wall key "<<key<<" got "<<Key;
	}
      }
      else
	Key=key;
    }
  }

  fs>>ID;
  for (int i=0;i<2; i++){
    fs>>key;
    MapPoint *mp=0;
    if (Bank){ 
      MapObject *m=Bank->getMapObject(key);
      if (m)
	{
	  mp=m->narrowPoint();	
	}
    }
    if (!mp){
      mp=new MapPoint();
      mp->Bank=Bank;
      if (Bank){
	Bank->add(mp,key);    
	if (key!=mp->Key){
	  std::cerr<<"wanted point key "<<key<<" got "<<mp->Key;
	}
      }      
    }
    fs>>mp->ID>>mp->Index>>mp->Local_map;
    fs>>mp->X[0]>>mp->X[1]>>mp->X[2];
    if (mp)addPoint(mp,i);
  }
  fs>>NumberSharedPoints>>Consecutive>>LastDistance;
  fs>>key;
  Bdual.reallocate(4,key);
  long *indx=getIndex();
  for (int i=0; i<4; i++)
    fs>>indx[i];
  fs>>StartLimit>>EndLimit>>Ls>>InfoInitialized>>CanLimit;
  if (NumberSharedPoints)
    {   
      if (PSharedIndex==0)
	{
	  PSharedIndex=new short[4];
	}
      for (int i=0;i<4;i++)
	fs>>PSharedIndex[i];
    }
  recenter();
  return this;
}
void MapWall::get(GenericData &gd)
{
  MapFeature::get(gd);
  long r=gd.ShortData.Rows;
  long c=gd.ShortData.Columns;
  if (c<4)c=4;
  gd.forceShortDataSize(r+1,c);
  gd.ShortData(r,0)=StartLimit;
  gd.ShortData(r,1)=EndLimit;
  gd.ShortData(r,2)=InfoInitialized;
  gd.ShortData(r,3)=CanLimit;
  r=gd.Data.Rows;
  c=gd.Data.Columns;
  gd.Data.grow(r+1,c);
  gd.Data(r,0)=Ls;
 }
int MapWall::set(GenericData &gd)
{
  if (MapFeature::set(gd))return 1;
  long r=gd.ShortData.Rows;
  long c=gd.ShortData.Columns;
  if (r<6)return 1;
  if (c<4)return 1;
  if (gd.Data.Rows<7)return 1;
  if (gd.Data.Columns<1)return 1;
  StartLimit=gd.ShortData(r-1,0);
  EndLimit=gd.ShortData(r-1,1);
  InfoInitialized=gd.ShortData(r-1,2);
  CanLimit=gd.ShortData(r-1,3);
  r=gd.Data.Rows;
  Ls=gd.Data(r-1,0);
  return 0;
}
