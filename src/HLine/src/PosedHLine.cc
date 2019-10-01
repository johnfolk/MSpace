
// = RCSID
//    $Id: PosedHLine.cc ,v 1.1 2004/04/1 
//
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "PosedHLine.hh"

using namespace Cure;


PosedHLine::PosedHLine(PosedHLine **pt,MapBank *b, long fKey, long id):
  PosedFeature(b, fKey,id)
{
  init(); 
  CastPtr=pt;
  FocalLength=0;
}
PosedHLine::PosedHLine(PosedHLine **pt,double focalLength, MapHLine *cb, long id):
  PosedFeature(cb->Bank, cb->Key,id)
{
  init(); 
  CastPtr=pt;
  FocalLength=focalLength;
}
PosedHLine::PosedHLine(PosedHLine & pw)
{
  PosedHLine::operator=(pw);
}
void PosedHLine::operator =(PosedHLine & pw)
{
  init();
  equal(pw);
  CastPtr=pw.CastPtr;
  Line=pw.Line;
  PixelLine=pw.PixelLine;
  Tangent=pw.Tangent;
  Reflected=pw.Reflected;
  FocalLength=pw.FocalLength;
}
void PosedHLine::init()
{
  Number2D=0;
  Number3D=2;
  NumberScalars=0;
  FullDim=6;
  Points3D=AllocatedPoints3D;
  Points3D[0]=&Line.StartPoint;
  Points3D[1]=&Line.EndPoint;
  Reflected=0;
  Jof.reallocate(6);
  Jof=1;
  Jop.reallocate(6,0);
  Z.reallocate(2,1);
  Xo.reallocate(6,1);
} 

int PosedHLine::transform(Transformation3D& pose,unsigned short covType)
{
  if (PosedFeature::transform(pose,covType))return MAP_OBJECT_INVALID;
  Vertical(0)=0;
  Vertical(1)=0;
  Vertical(2)=1;
  pose.rotate(Vertical,Vertical);
  if ((Line.StartPoint(1)<=0)&&(Line.EndPoint(1)<=0))return NOT_VISABLE;
  Vector3D cv;
  cv=Line.StartPoint.cross(Line.EndPoint);
  double magcv=(cv*cv);
  if (magcv<1E-20)return NOT_VISABLE;
  magcv=sqrt(magcv);
  cv/=magcv;
  double c[3];
  pose.invRotate(cv.X,c);
  int r=0;
  MapFeature *mf =getFeature();
  if (mf)
    r=mf->testVisable(c);
  if (r)return r;
  return calcZ();
}
int PosedHLine::calcZ()
{
  if (Line.getPixels(PixelLine.StartPoint.X,PixelLine.EndPoint.X,
		     FocalLength))return NOT_VISABLE;
  Tangent=Line.EndPoint-Line.StartPoint;
  Tangent/=sqrt(Tangent*Tangent);  
  Reflected=0;
  PixelLine.getPolarValues(Z.Element,0);
  if (Z(0,0)>M_PI_2)
    {
      Reflected=1;
      Z(0,0)-=M_PI;
      Z(1,0)*=-1;
    }
  if (Z(0,0)<-M_PI_2)
    {
      Reflected=1;
      Z(0,0)+=M_PI;
      Z(1,0)*=-1;
    }
  return 0;
}
int PosedHLine::roughMatch(double *z,double *thresholds)
{
  if (!(CalcState&1))return -1;
  double innov[2];
  innov[0]=z[0]-Z(0,0);
  innov[1]=z[1]-Z(1,0);
  while (innov[0]>M_PI)innov[0]-=TWO_PI;
  while (innov[0]<-M_PI)innov[0]+=TWO_PI;
  if (innov[0]>M_PI_2)
    {
      innov[0]-=M_PI;
      innov[1]+=Z(1,0);
      innov[1]*=-1;
      innov[1]-=Z(1,0);
    }
  if (innov[0]<-M_PI_2)
    {
      innov[0]+=M_PI;
      innov[1]+=Z(1,0);
      innov[1]*=-1;
      innov[1]-=Z(1,0);
    }
  if ((innov[0]>thresholds[0])||(innov[0]<-thresholds[0]))return 1;
  if ((innov[1]>thresholds[1])||(innov[1]<-thresholds[1]))return 1;
  return 0;
}
double PosedHLine::fineMatch(double *z,Cure::Matrix & sigma)
{
  if (!(CalcState&1))return -1;
  Matrix innov(2,1),a,b;
  for (int i=0; i<2;i++)
    innov(i,0)=z[i]-Z(i,0);
  while (innov(0,0)>M_PI)innov(0,0)-=TWO_PI;
  while (innov(0,0)<-M_PI)innov(0,0)+=TWO_PI;
  if (innov(0,0)>M_PI_2)
    {
      innov(0,0)-=M_PI;
      innov(1,0)+=Z(1,0);
      innov(1,0)*=-1;
      innov(1,0)-=Z(1,0);
    }
  if (innov(0,0)<-M_PI_2)
    {
      innov(0,0)+=M_PI;
      innov(1,0)+=Z(1,0);
      innov(1,0)*=-1;
      innov(1,0)-=Z(1,0);
    }
  a.multiply_(sigma,innov);
  b.multTranspose_(innov,a,1);
  return b(0,0);
}
int PosedHLine::calcEta(Cure::Matrix & v,unsigned short mtype,
			unsigned short noJac)
{
  if (!(CalcState&1))return 1;
  setMeasurementType(mtype);
  if (CalcState&4)return 0;
  if (v.Rows<5)MeasurementType=0; //We don't handle this case
  if (MeasurementType&3)
    {
      //MeasurementType-=3;
      FocalLength=v(4,0);
      if (calcZ())return 1;
      Vector3D s,e,cv,t,temp,temp2;
      t=Line.EndPoint-Line.StartPoint;
      double d=(t*t);
      if (d<1E-20)return 1;
      d=sqrt(d);
      t/=d;
      if (v.Rows==7)
	{
	  s(0)=-(v(0,0)+v(5,0));
	  s(1)=FocalLength;
	  s(2)=-(v(1,0)+v(6,0));
	  e(0)=-(v(2,0)+v(5,0));
	  e(1)=FocalLength;
	  e(2)=-(v(3,0)+v(6,0));
	}
      else if (v.Rows==5)
	{
	  s(0)=-(v(0,0));
	  s(1)=FocalLength;
	  s(2)=-(v(1,0));
	  e(0)=-(v(2,0));
	  e(1)=FocalLength;
	  e(2)=-(v(3,0));
	}
      cv=s.cross(e);
      double magcv=(cv*cv);
      if (magcv<1E-20)return 1;
      magcv=sqrt(magcv);
      cv/=magcv;
      int k=1;
      if (MeasurementType&2)k=2;
      Eta.reallocate(k,1);
      Eta(0,0)=t*cv;
      Jeo.reallocate(k,6);
      Jeo(0,3)=cv(0)/d;
      Jeo(0,0)=-Jeo(0,3);
      Jeo(0,4)=cv(1)/d;
      Jeo(0,1)=-Jeo(0,4);
      Jeo(0,5)=cv(2)/d;
      Jeo(0,2)=-Jeo(0,5);
      double b=t(0)*Eta(0,0)/d;
      Jeo(0,3)-=b;
      Jeo(0,0)+=b;
      b=t(1)*Eta(0,0)/d;
      Jeo(0,4)-=b;
      Jeo(0,1)+=b;
      b=t(2)*Eta(0,0)/d;
      Jeo(0,5)-=b;
      Jeo(0,2)+=b;
      Jev.reallocate(k,v.Rows);
      temp=t.cross(e);
      temp/=magcv;
      Jev(0,0)=temp(0);
      Jev(0,1)=temp(2);
      Jev(0,4)=-temp(1);
      temp=s.cross(t);
      temp/=magcv; 
      Jev(0,2)=temp(0);
      Jev(0,3)=temp(2);
      Jev(0,4)-=(temp(1));
      
      double a=Eta(0,0)/magcv;
      Jev(0,0)-=a*(cv(1)*e(2)-cv(2)*v(4,0));
      Jev(0,1)-=a*(cv(0)*v(4,0)-cv(1)*e(0));
      Jev(0,2)+=a*(cv(1)*s(2)-cv(2)*v(4,0));
      Jev(0,3)+=a*(cv(0)*v(4,0)-cv(1)*s(0));
      Jev(0,4)+=a*(cv(0)*(v(3,0)-v(1,0))+cv(2)*(v(0,0)-v(2,0)));
      if (v.Rows==7)
	{
	  Jev(0,5)=Jev(0,0)+Jev(0,2);
	  Jev(0,6)=Jev(0,1)+Jev(0,3);
	}   
      if (k==2)
	{
	  Vector3D r;
	  r=(Line.StartPoint+Line.EndPoint);
	  double d2=(r*t);
	  temp=t;
	  temp*=d2;
	  r-=temp;
	  Eta(1,0)=r*cv;;
	  Jeo(1,0)=(cv(0)-d2*Jeo(0,0)-Eta(0,0)*(t(0)-r(0)/d));
	  Jeo(1,1)=(cv(1)-d2*Jeo(0,1)-Eta(0,0)*(t(1)-r(1)/d));
	  Jeo(1,2)=(cv(2)-d2*Jeo(0,2)-Eta(0,0)*(t(2)-r(2)/d));
	  Jeo(1,3)=(cv(0)-d2*Jeo(0,3)-Eta(0,0)*(t(0)+r(0)/d));
	  Jeo(1,4)=(cv(1)-d2*Jeo(0,4)-Eta(0,0)*(t(1)+r(1)/d));
	  Jeo(1,5)=(cv(2)-d2*Jeo(0,5)-Eta(0,0)*(t(2)+r(2)/d));
	  temp=r.cross(e);
	  temp/=magcv;
	  Jev(1,0)=temp(0);
	  Jev(1,1)=temp(2);
	  Jev(1,4)=-temp(1);
	  temp=s.cross(r);
	  temp/=magcv;
	  Jev(1,2)=temp(0);
	  Jev(1,3)=temp(2);
	  Jeo(1,4)-=temp(1);
	  a=Eta(1,0)/magcv;
	  Jev(1,0)-=a*(cv(1)*e(2)-cv(2)*v(4,0));
	  Jev(1,1)-=a*(cv(0)*v(4,0)-cv(1)*e(0));
	  Jev(1,2)-=a*(cv(1)*(-s(2))+cv(2)*v(4,0));
	  Jev(1,3)+=a*(cv(0)*v(4,0)-cv(1)*s(0));
	  Jev(1,4)+=a*(cv(0)*(v(3,0)-v(1,0))+cv(2)*(v(0,0)-v(2,0)));
	  if (v.Rows==7)
	    {
	      Jev(1,5)=Jev(1,0)+Jev(1,2);
	      Jev(1,6)=Jev(1,1)+Jev(1,3);
	    }
	}
    }  
  else
    {
      Eta.reallocate(0);
      Jev.reallocate(0,v.Rows);
      Jeo.reallocate(0,6);
    }
  CalcState=(CalcState|0x4);
  return 0;
}
/*
int PosedHLine::addInfo(Match &mat,
			Pose3D & pose, const int type)
{
  
  MapFeature *mf =getFeature();
  if (mf)
   {
     Measurement *m=mat.Measure;
     Matrix v(1,11);
     FocalLength=m->V(4,0);//info(0,3);
     v(0,0)=mat.PathDistance;//info(0,0);
     v(0,1)=mat.Weight;//info(0,1);
     pose.getXYZ(v.Element+2);
     v(0,5)=-m->V(0,0);//info(0,2);
     v(0,7)=-m->V(1,0);//info(0,4);    
     v(0,6)=m->V(4,0)/(1-m->V(4,0)/Line.StartPoint(1));
     //FocalLength/(1-FocalLength/Line.StartPoint(1));
     v(0,9)=m->V(4,0)/(1-m->V(4,0)/Line.EndPoint(1));
     v(0,8)=-m->V(2,0);//info(0,5);
     v(0,10)=-m->V(3,0);//info(0,7);    
     pose.invRotate(v.Element+5,v.Element+5);
     pose.invRotate(v.Element+8,v.Element+8);
     int res=mf->addInfo(v,type);
     transform(pose);
     return res;
    }
  return -1;
}
*/
int PosedHLine::addInfo(Match &mat,
			Transformation3D & pose, unsigned short covType, 
			Transformation3D & map2info,			
			const int type)
{   
  MapFeature *mf =getFeature();
  if (mf)
   {
     Measurement *m=mat.Measure;
     Matrix v(1,11);
     FocalLength=m->V(4,0);
     v(0,0)=mat.PathDistance;
     v(0,1)=mat.Weight;
     pose.getXYZ(v.Element+2);
     map2info.transform(v.Element+2,v.Element+2);
     v(0,5)=-m->V(0,0);
     v(0,7)=-m->V(1,0);
     v(0,6)=m->V(4,0)/(1-m->V(4,0)/Line.StartPoint(1));
     v(0,9)=m->V(4,0)/(1-m->V(4,0)/Line.EndPoint(1));
     v(0,8)=-m->V(2,0);
     v(0,10)=-m->V(3,0);
     pose.invRotate(v.Element+5,v.Element+5);
     pose.invRotate(v.Element+8,v.Element+8);
     map2info.rotate(v.Element+5,v.Element+5);
     map2info.rotate(v.Element+8,v.Element+8);
     int res=mf->addInfo(v,map2info,type);
     transform(pose,covType);
     return res;
    }
  return -1;
}
/*
int PosedHLine::pixel(PixelFragment& pix,double  bottomLeft[2],
		      double topRight[2],double pixelWidth,
		      double pixelHeight)
{
  if ((Line.StartPoint(1)<=0)&&(Line.EndPoint(1)<=0))return NOT_VISABLE;
  if(!(lineInRectangle(PixelLine.StartPoint.X,PixelLine.EndPoint.X, 
		       bottomLeft,topRight)))return NOT_VISABLE;
  pix.clean();
  int s[2],e[2];
  double r[2];
  pixelize(PixelLine.StartPoint.X,s,pixelWidth, pixelHeight);
  pixelize(PixelLine.EndPoint.X,e,pixelWidth, pixelHeight);
  double t[3];
  if ((s[1]<e[1])||((s[1]==e[1])&&(s[0]>e[0])))
    {
     int d=s[1];
      s[1]=e[1];
      e[1]=d;
      d=s[0];
      s[0]=e[0];
      e[0]=d;
      r[0]=Line.EndPoint(1);
      r[1]=Line.StartPoint(1);
    }
  else
    {
      r[0]=Line.EndPoint(1);
      r[1]=Line.StartPoint(1);
    }
  t[0]=e[0]-s[0];
  t[1]=e[1]-s[1];
  t[3]=r[1]-r[0];
  double step=(t[0]*t[0]+t[1]*t[1]);
  pix.add(s[0],s[1],r[0]);
  int last[2];
  last[0]=s[0];
  last[1]=s[1];
  if (step==0)return 0;
  step=sqrt(.5/step);
  double f=step;
  while (f<1)
    {
      if (t[0]>0)
	e[0]=s[0]+(int)(f*(t[0])+.5);
      else
	e[0]=s[0]+(int)(f*(t[0])-.5);
      if (t[1]>0)
	e[1]=s[1]+(int)(f*(t[1])+.5);
      else
	e[1]=s[1]+(int)(f*(t[1])-.5);
      if ((e[0]!=last[0])||(e[1]!=last[1]))
	{
	  r[1]=r[0]+f*t[3];
	  pix.add(e[0],e[1],r[1]);
	  last[0]=e[0];
	  last[1]=e[1];
	}
      f+=step;
    }
  return 0;
}
*/
