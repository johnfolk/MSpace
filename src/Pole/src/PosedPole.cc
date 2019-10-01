// = RCSID
//    $Id: posedPole.cc ,v 1.1 2004/04/1 
//
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "PosedPole.hh"

using namespace Cure;
PosedPole::PosedPole(PosedPole **tp,Trig *t,MapBank *b, long fKey, long id):
  PosedFeature(b, fKey, id)
{
  init(); 
  CastPtr=tp;
  StartAngle=0;
  EndAngle=0;
  Triger=t;
  Normal[0]=1;
  Normal[1]=0;
  A[0]=0;
  A[1]=0;
  A[2]=0;
  A[3]=0;
  U[0]=0;
  U[1]=0;
  U[2]=0;
}

PosedPole::PosedPole(PosedPole **tp,Trig *t,MapPole *tr,long id):
  PosedFeature(tr->Bank, tr->Key, id)
{
  init(); 
  CastPtr=tp;
  StartAngle=0;
  EndAngle=0;
  Triger=t;
  Normal[0]=1;
  Normal[1]=0;
  A[0]=0;
  A[1]=0;
  A[2]=0;
  A[3]=0;
  U[0]=0;
  U[1]=0;
  U[2]=0;
}
PosedPole::PosedPole(PosedPole & pt)
{
  PosedPole::operator=(pt);
}
void PosedPole::init()
{
  Number2D=1;
  Number3D=0;
  NumberScalars=1;
  FullDim=3;
  Points2D=AllocatedPoints2D;
  Points2D[0]=&Center;
  Jof.reallocate(3);
  Jof(0,2)=0;
  Jof(1,2)=0;
  Jof(2,2)=0;
  Jof(2,0)=0;
  Jof(2,1)=0;
  Z.reallocate(3,1);
  Xo.reallocate(3,1);
  Jop.reallocate(3,0);
  Z=0;
  Xo=0;
} 
void PosedPole::operator =(PosedPole & pt)
{
  init();  
  equal(pt);
  CastPtr=pt.CastPtr;
  Center=pt.Center;
  StartAngle=pt.StartAngle;
  EndAngle=pt.EndAngle;
  Triger=pt.Triger;
  Normal[0]=pt.Normal[0];
  Normal[1]=pt.Normal[1];
  A[0]=pt.A[0];
  A[1]=pt.A[1];
  A[2]=pt.A[2];
  A[3]=pt.A[3];
  U[0]=pt.U[0];
  U[1]=pt.U[1];
  U[2]=pt.U[2];
  CalcState=pt.CalcState;
  MeasurementType=pt.MeasurementType;
  Z=pt.Z;
  Jeo=pt.Jeo;
  Jep=pt.Jep;
  Jos=pt.Jos;
  Jof=pt.Jof;
  Jop=pt.Jop;
  Xo=pt.Xo;
  Eta=pt.Eta;
}

int PosedPole::transform(Transformation3D& pose,unsigned short covType)
{
  if (PosedFeature::transform(pose,covType))return MAP_OBJECT_INVALID;
  Center.getPolarValues(Z.Element);
  Z(2,0)=Xo(2,0);
  pose.getRot2D(A);
  A[0]*=Xo(2,0);
  A[1]*=Xo(2,0);
  A[2]*=Xo(2,0);
  A[3]*=Xo(2,0);
  Normal[0]=cos(Z(0,0));
  Normal[1]=sin(Z(0,0));
  double tempz[2];
  pose.invRotate2D(Normal,tempz);
  tempz[0]=-tempz[0];
  tempz[1]=-tempz[1];
  U[1]=Triger->intAtan(tempz);
  double temp=Center(0)*A[3]-Center(1)*A[1];
  double a=Center(0)*A[2]-Center(1)*A[0];
  double b=A[1]*A[2]-A[3]*A[0];

  if ((temp<.0001)&&(temp>-.0001))
    {
      if ((a<.0001)&&(a>-.0001))
	{
	  U[0]=U[1];
	  U[2]=U[1];
	}
      else 
	{
	  U[0]=Triger->intAsin((-b/a));
	  if (U[0]>90)
	    {
	      U[0]=U[1];
	      U[2]=U[1];
	    }
	  else
	    {
	      U[2]=180-U[0];
	      if (U[2]>180)U[2]-=360;
	    }
	}
    }
  else
    {
      a/=temp;
      b/=temp;
      temp=a*a;
      temp+=1;
      double root=temp-(b*b);
      if (root<0)
	{
	  U[0]=U[1];
	  U[2]=U[1];
	}
      else
	{
	  root=sqrt(root);
	  double tmp=a*b;
	  double su1=(tmp+root)/temp;
	  double su2=(tmp+root)/temp;
	  int ans[2];
	  ans[0]=Triger->intAsin(su1);
	  ans[1]=Triger->intAsin(su2);
	  if(ans[0]>90)
	    {
	      if(ans[1]>90)
		{
		  U[0]=U[1];
		  U[2]=U[1];
		}
	      else
		{
		  U[0]=ans[1];
		  U[2]=180-ans[1];
		}
	    }
	  else if(ans[1]>90)
	    {
	      U[0]=ans[0];
	      U[2]=180-ans[0];
	    }
	  else
	    {
	      double test[4];
	      test[0]=a*su1+b-Triger->cosOf(ans[0]);
	      test[1]=a*su1+b-Triger->cosOf(180-ans[0]);
	      test[2]=a*su2+b-Triger->cosOf(ans[1]);
	      test[3]=a*su2+b-Triger->cosOf(180-ans[1]);
	      for (int i=0; i<4; i++)
		if (test[i]<0)test[i]=-test[i];
	      if (test[0]<test[1])
		{
		  U[0]=ans[0];
		  U[2]=180-ans[0];
		}
	      else
		{
		  U[2]=ans[0];
		  U[0]=180-ans[0];
		  double d=test[0];
		  test[0]=test[1];
		  test[1]=d;
		}
	      if (test[2]<test[1])
		{
		  U[2]=ans[1];
		  test[1]=test[2];
		  if (test[0]>test[1])
		    {
		      int i=U[0];
		      U[0]=U[2];
		      U[2]=i;
		      double d=test[0];
		      test[0]=test[1];
		      test[1]=d;
		    }
		}
	      if (test[3]<test[1])
		U[2]=180-ans[1];
	    }
	}
    }
  int ang[3];
  for (int i=0; i<3;i++)
    ang[i]=thetaOfU(U[i]);
  int tst=ang[0]-ang[1];
  if (tst>180)ang[0]-=360;
  if (tst<-180)ang[0]+=360;
  tst=ang[2]-ang[1];
  if (tst>180)ang[2]-=360;
  if (tst<-180)ang[2]+=360;
  if (ang[0]>ang[2])
    {
      tst=ang[0];
      ang[0]=ang[2];
      ang[2]=tst;
      tst=U[0];
      U[0]=U[2];
      U[2]=tst;
    }
  EndAngle=Triger->toRads(ang[2]);
  StartAngle=Triger->toRads(ang[0]);
  CalcState=1;
  return 0;
}
int PosedPole::calcEta(Cure::Matrix & v,unsigned short mtype)
{
  if (!(CalcState&1))return -1; 
  setMeasurementType(mtype);
  if (CalcState&4)return 0;
  int dim=2;
  if (MeasurementType&1)dim=3;
  if (v.Rows!=dim)return -1;
  Eta.reallocate(dim,1);
  Jev.reallocate(dim);
  Jev=1;
  Jeo.reallocate(dim,3);
  Jeo.Columns=2;
  Jeo.Rows=2;
  Center.getPolarJacobian(Jeo);
  Jeo.Columns=3;
  Jeo(0,2)=0;
  Jeo(1,2)=1;
  Eta(0,0)=v(0,0)-Z(0,0);
  Eta(1,0)=v(1,0)-Z(1,0);
  while(Eta(0,0)>M_PI)Eta(0,0)-=TWO_PI;
  while(Eta(0,0)<-M_PI)Eta(0,0)+=TWO_PI;
  if(MeasurementType&1)
    {
      Jeo.Rows=3;
      Jeo(2,0)=0;
      Jeo(2,1)=0;
      Jeo(2,2)=1;
      Eta(2,0)=v(2,0)-Z(2,0);
    }
  CalcState=(CalcState|0x4);
  return 0;
}
int PosedPole::roughMatch(double *z,double *thresholds)
{
  if (!(CalcState&1))return -1;
  double innov=z[0]-Z(0,0);
   while (innov>M_PI)innov-=TWO_PI;
  while (innov<-M_PI)innov+=TWO_PI;
  if ((innov>thresholds[0])||(innov<-thresholds[0]))return 1;
  innov=z[1]-Z(1,0); 
  if ((innov>thresholds[1])||(innov<-thresholds[1]))return 1;
  if(thresholds[2]>0)
    {
      innov=z[2]-Z(2,0);
      if ((innov>thresholds[2])||(innov<-thresholds[2]))return 1;
    }
  return 0;
}
double PosedPole::fineMatch(double *z,Cure::Matrix  & sigma)
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

int PosedPole::thetaOfU(int u)
{
  double temp[2];
  double n[2];
  n[0]=Triger->cosOf(u);
  n[1]=Triger->sinOf(u);
  temp[0]=Center(0)+A[0]*n[0]+A[1]*n[1];
  temp[1]=Center(1)+A[2]*n[0]+A[3]*n[1];
  return Triger->intAtan(temp);
}
int PosedPole::thetaOfU(double *r,int u)
{
  double temp[2];
  double n[2];
  n[0]=Triger->cosOf(u);
  n[1]=Triger->sinOf(u);
  temp[0]=Center(0)+A[0]*n[0]+A[1]*n[1];
  temp[1]=Center(1)+A[2]*n[0]+A[3]*n[1];
  r[0]=sqrt(temp[1]*temp[1]+temp[0]*temp[0]);
  return Triger->intAtan(temp);
}
int PosedPole::scan(ScanFragment& s)
{
  s.clean();
  int n=Triger->toDegrees(EndAngle-StartAngle);
  if (n<1)n=1;
  double r[n];
  int a[n];
  a[0]=thetaOfU(r,U[0]);
  s.add(a[0],r[0]);
  n--;
  a[n]=thetaOfU(r+n,U[2]);
  double temp=0;
  int i=thetaOfU(&temp,U[1]);
  int index=i-a[0];
  if (index<-180) index+=360;
  if (index>180)index-=360;
  if (index<0)return 1;
  if (index>n)return 1;
  a[index]=i;
  r[index]=temp;
  if (index>1)
    recursiv(U[0],U[1],&a[0],&r[0],index);
  index=n-index;
  if (index>1)
    recursiv(U[1],U[2],&a[1],&r[1],index);
  n++;
  for (int i=1; i<n;i++)
    s.add(a[i],r[i]);
  return 0;
}
void PosedPole::recursiv(int u0,int u1,int *a,double *r,int length)
{
  int i=u1-u0;
  if (i<0) i+=360;
  if (i<180) i-=360;
  i=(i>>1);
  double temp=0;
  int ang=thetaOfU(&temp,i);
  int index=ang-a[0];
  if (index<-180) index+=360;
  if (index>180)index-=360;
  if ((index>0)&&(index<length))
    {
      a[index]=ang;
      r[index]=temp;
      if (index>1)recursiv(u0,i,a, r,index);
      int inx=length-index;
      if (inx>1)recursiv(i,u1,&a[index], &r[index],inx);
      return;
    }
  if (index==0)
    {
      recursiv(i,u1,a, r,length);
      return;
    }
  recursiv(u0,i,a, r,length);
}
