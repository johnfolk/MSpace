// = RCSID
//    $Id: PosedWall.cc ,v 1.1 2004/04/1 
//
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "PosedWall.hh"

using namespace Cure;

PosedWall::PosedWall(PosedWall **wp,Trig *t,MapBank *b, long fKey, long id):
  PosedFeature(b, fKey,id)
{
  init(); 
  CastPtr=wp;
  Triger=t;
}
PosedWall::PosedWall(PosedWall **wp,Trig *t,MapWall *w, long id):
  PosedFeature(w->Bank,w->Key,id)
{
  init();  
  CastPtr=wp;
  Triger=t;
}

PosedWall::PosedWall(PosedWall & pw)
{

  PosedWall::operator=(pw);
}
void PosedWall::operator =(PosedWall & pw)
{
  init(); 
  equal(pw);
  CastPtr=pw.CastPtr;
  Line=pw.Line;
  Triger=pw.Triger;
  EndPointUpdateThreshold=pw.EndPointUpdateThreshold;
}

void PosedWall::init()
{
  Number2D=2;
  Number3D=0;
  NumberScalars=0;
  FullDim=4;
  Points2D=AllocatedPoints2D;
  Points2D[0]=&Line.StartPoint;
  Points2D[1]=&Line.EndPoint;
  Jof.reallocate(4);
  Jof=1;
  Jop.reallocate(4,0);
  Z.reallocate(4,1);
  Z=0;
  Xo.reallocate(4,1);
  EndPointUpdateThreshold=.01;
} 

int PosedWall::transform(Transformation3D& pose,unsigned short covType)
{
  if (PosedFeature::transform(pose,covType))return MAP_OBJECT_INVALID;
 return calcZ();
}
int PosedWall::calcZ()
{
  Z.reset(4,1);
  Line.getPolarValues(Z.Element,3);
  while(Z(3,0)<Z(2,0))Z(3,0)+=TWO_PI;
  double d=Z(3,0)-Z(2,0);
  if (d>M_PI) return NOT_VISABLE;
  return 0;
}

int PosedWall::roughMatch(double  *z,double *thresholds)
{
  if (!(CalcState&1))return -1;
  double innov[2];
  innov[0]=z[0]-Z(0,0);
  innov[1]=z[1]-Z(1,0);
  while (innov[0]>M_PI)innov[0]-=TWO_PI;
  while (innov[0]<-M_PI)innov[0]+=TWO_PI;
  if ((innov[0]>thresholds[0])||(innov[0]<-thresholds[0]))return 1;
  if ((innov[1]>thresholds[1])||(innov[1]<-thresholds[1]))return 1;
  return 0;
}
double PosedWall::fineMatch(Cure::Measurement &m,Cure::Matrix & sigma){
  if (!(CalcState&1))return -1;
  Line2D ln(m.V.Element,m.V.Element+2);
  double z[4];
  Z.reset(4,1);
  ln.getPolarValues(z,3);
  while(z[3]<z[2])z[3]+=TWO_PI;
  while ((z[2]-Z(2,0))>TWO_PI){
    z[2]-=TWO_PI;
    z[3]-=TWO_PI;
  }
  while ((z[2]-Z(2,0))<-TWO_PI){
    z[2]+=TWO_PI;
    z[3]+=TWO_PI;
  }
  if (Z(2,0)>z[2])
    {
      if ((z[3]-Z(2,0))<.05)return -1;
    }
  else if ((Z(3,0)-z[2])<.05)return -1;
  Matrix innov(2,1),a,b;
  for (int i=0; i<2;i++)
    innov(i,0)=z[i]-Z(i,0);
  while (innov(0,0)>M_PI)innov(0,0)-=TWO_PI;
  while (innov(0,0)<-M_PI)innov(0,0)+=TWO_PI;
  a.multiply_(sigma,innov);
  b.multTranspose_(innov,a,1);
  return b(0,0);
}

double PosedWall::fineMatch(double *z,Cure::Matrix  & sigma)
{
  if (!(CalcState&1))return -1;
  Matrix innov(2,1),a,b;
  for (int i=0; i<2;i++)
    innov(i,0)=z[i]-Z(i,0);
  while (innov(0,0)>M_PI)innov(0,0)-=TWO_PI;
  while (innov(0,0)<-M_PI)innov(0,0)+=TWO_PI;
  a.multiply_(sigma,innov);
  b.multTranspose_(innov,a,1);
  return b(0,0);
}
int PosedWall::consistentPose(Cure::Measurement &m, 
			      Cure::Transformation3D &pose)
{
  double x[3];
  x[0]=0;
  x[1]=0;
  x[2]=0;
  pose.setXYTheta(x);
  transform(pose,0);
  Line2D ln(m.V.Element,m.V.Element+2);
  double z[4];
  Z.reset(4,1);
  ln.getPolarValues(z,3);
  x[2]=Z(0,0)-z[0];
  double c=cos(Z(0,0));
  double s=sin(Z(0,0));
  x[0]=c*(Z(1,0)-z[1]);
  x[1]=s*(Z(1,0)-z[1]);

  double d=c*(Line.StartPoint(1)+Line.EndPoint(1))/2-
    s*(Line.StartPoint(0)+Line.EndPoint(0))/2;
  x[0]-=s*d;
  x[1]+=c*d;
  double cv=sqrt(m.CovV(0,0)+m.CovV(2,2)+m.CovV(1,1)+m.CovV(3,3));
  d=Line.getLength()+cv;  
  d*=(rand()-RAND_MAX/2);
  d/=RAND_MAX;
  d-=z[1]*(sin(z[3]-z[0])+sin(z[2]-z[0]));
  x[0]-=s*d;
  x[1]+=c*d;
  pose.setXYTheta(x);
  transform(pose,3);
  d=0;
  if (z[3]>Z(3,0))d-=z[1]*sin(z[3]-Z(3,0));
  if (z[2]<Z(2,0))d+=z[1]*sin(z[2]-Z(2,0));
  x[0]-=s*d;
  x[1]+=c*d;
  d=(rand()-RAND_MAX/2);
  d/=RAND_MAX;
  x[0]+=(d*cv); 
  d=(rand()-RAND_MAX/2);
  d/=RAND_MAX;
  x[1]+=(d*cv); 
  d=(rand()-RAND_MAX/2);
  d/=RAND_MAX;
  x[2]+=(d*cv/z[1]); 
  pose.setXYTheta(x);
  transform(pose,0);
  return 0;
}

int PosedWall::calcEta(Cure::Matrix & v,unsigned short mtype,
		       unsigned short noJac)
{
 if (!(CalcState&1))return 1; 
  setMeasurementType(mtype);
  if (MeasurementType&0x2){
    double ends[2];
    ends[0]=(Line.StartPoint(0)-v(0,0));
    ends[0]*=ends[0];
    ends[1]= (Line.StartPoint(1)-v(1,0));
    ends[1]*=ends[1];
    ends[0]+=ends[1];
    if(ends[0]>EndPointUpdateThreshold)
      MeasurementType=MeasurementType&0xFFFD;
  }
  if (MeasurementType&0x4){
    double ends[2];
    ends[0]=(Line.EndPoint(0)-v(2,0));
    ends[0]*=ends[0];
    ends[1]= (Line.EndPoint(1)-v(3,0));
    ends[1]*=ends[1];
    ends[0]+=ends[1];
    if(ends[0]>EndPointUpdateThreshold)
      MeasurementType=MeasurementType&0xFFFB;
  }
  if (CalcState&4)return 0;
  if (v.Rows!=4)MeasurementType=0;
  if (MeasurementType==0)
    {
      Jev.reallocate(0,v.Rows);
      Jeo.reallocate(0,4);
      Eta.reallocate(0,1);
      CalcState=(CalcState|0x4);
      return 0;
    }
  else if (MeasurementType<8)
    {
      if (calcZ())return 1;
      int k=2;
      if (MeasurementType&0x2)k++;
      if (MeasurementType&0x4)k++;
      Eta.reallocate(k,1);
      Jeo.reallocate(k,4);
      Jev.reallocate(k,v.Rows);
      double to[2], tv[2], ro[2],rv[2];
      to[0]=Line.EndPoint(0)-Line.StartPoint(0);
      to[1]=Line.EndPoint(1)-Line.StartPoint(1);
      double d=to[0]*to[0]+to[1]*to[1];
      if (d<1E-20)return 1;
      d=sqrt(d);
      to[0]/=d;
      to[1]/=d;
      tv[0]=v(2,0)-v(0,0);
      tv[1]=v(3,0)-v(1,0);
      ro[0]=Line.EndPoint(0)+Line.StartPoint(0);
      ro[1]=Line.EndPoint(1)+Line.StartPoint(1);
      rv[0]=v(2,0)+v(0,0);
      rv[1]=v(3,0)+v(1,0);
      Eta(0,0)=to[0]*tv[1]-to[1]*tv[0];
      Eta(1,0)=(ro[0]-rv[0])*tv[1]-(ro[1]-rv[1])*tv[0];
      Jev(0,0)=to[1];
      Jev(0,1)=-to[0];
      Jev(0,2)=-to[1];
      Jev(0,3)=to[0];
      Jeo(0,0)=-tv[1]/d;
      Jeo(0,1)=tv[0]/d;
      Jeo(0,2)=-Jeo(0,0);
      Jeo(0,3)=-Jeo(0,1);
      double a=Eta(0,0)/d;
      double b=a*to[0];	
      a*=to[1];
      Jeo(0,0)+=b;
      Jeo(0,1)+=a;
      Jeo(0,2)-=b;
      Jeo(0,3)-=b;
    
      Jeo(1,0)=tv[1];
      Jeo(1,1)=-tv[0];
      Jeo(1,2)=tv[1];
      Jeo(1,3)=-tv[0];
      Jev(1,0)=ro[1]-2*v(3,0);
      Jev(1,1)=-ro[0]+2*v(2,0);
      Jev(1,2)=-ro[1]+2*v(1,0);
      Jev(1,3)=ro[0]-2*v(0,0);
      k=2;
      if (MeasurementType&0x2){
	 Eta(k,0)=(Line.StartPoint(0)-v(0,0))*tv[0]+
	   (Line.StartPoint(1)-v(1,0))*tv[1];
	 Jeo(k,0)=tv[0];
	 Jeo(k,1)=tv[1];
	 Jeo(k,2)=0;
	 Jeo(k,3)=0;
	 Jev(k,0)=(-tv[0]-Line.StartPoint(0)+v(0,0));
	 Jev(k,1)=(-tv[1]-Line.StartPoint(1)+v(1,0));
	 Jev(k,2)=(Line.StartPoint(0)-v(0,0));
	 Jev(k,3)=(Line.StartPoint(1)-v(1,0));
	 k++;
      }
      if (MeasurementType&0x4){
	Eta(k,0)=(Line.EndPoint(0)-v(2,0))*tv[0]+(Line.EndPoint(1)-v(3,0))*tv[1];
	Jeo(k,2)=tv[0];
	Jeo(k,3)=tv[1];
	Jeo(k,0)=0;
	Jeo(k,1)=0;
	Jev(k,2)=(-tv[0]+Line.EndPoint(0)-v(2,0));
	Jev(k,3)=(-tv[1]+Line.EndPoint(1)-v(3,0));
	Jev(k,0)=(-Line.EndPoint(0)+v(2,0));
	Jev(k,1)=(-Line.EndPoint(1)+v(3,0));
      }
    }
  CalcState=(CalcState|0x4);
  return 0;  
}

/*
int PosedWall::addInfo(Match &mat,
		       Pose3D & pose,  const int type)
{
  MapFeature *mf =getFeature();
   if (mf)
   {
     Measurement *m=mat.Measure;
     Matrix &info=m->W;
     Matrix v(info.Rows,3);
      int top=3*info.Rows;
     double dis=mat.PathDistance;
     int j=0;
     for (int i=0; i<top; i+=2, j+=2)
       {
	 v.Element[i]=dis;
	 i++;
	 pose.invTransform2D(info.Element+j,v.Element+i);
       }
     int res=mf->addInfo(v,type);
     transform(pose);
     return res;
    }
  return -1;
}
*/
int  PosedWall::addInfo(Match &mat,
			Cure::Transformation3D & pose,unsigned short covType,
			Cure::Transformation3D & map2info,
			const int type)
{
  MapFeature *mf =getFeature();
   if (mf)
   {
     Measurement *m=mat.Measure;
     Matrix &info=m->W;
     Matrix v(info.Rows,3);
     int top=3*info.Rows;
     double dis=mat.PathDistance;
     Transformation3D sens2info;
     sens2info-=pose;
     sens2info+=map2info;
     int j=0;
     for (int i=0; i<top; i+=2,j+=2)
       {
	 v.Element[i]=dis;
	 i++;
	 sens2info.transform2D(info.Element+j,v.Element+i);
       }
     int res=mf->addInfo(v,map2info,type);
     transform(pose,covType);
     return res;
    }
  return -1;
}
int PosedWall::scan(ScanFragment& s)
{

  Z.reset(4,1);
  double d=Z(3,0)-Z(2,0);
  if (d<-M_PI)d+=(2*M_PI);
  if (d<0) return NOT_VISABLE;
  s.clean();
  int j=Triger->toDegrees(Z(2,0));
  int top=Triger->toDegrees(Z(3,0))+1;
  double normal[2];
  normal[0]=cos(Z(0,0));
  normal[1]=sin(Z(0,0));
  for (;j<top; j++)
    {
      d=Z(1,0)*Triger->cosOfThetaMinusGamma(j,normal[0],normal[1]);
      s.add(j,d);
    }
  return 0;
}
