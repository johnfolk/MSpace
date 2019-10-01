// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "PosedRelativePointFeature.hh"


using namespace Cure;

PosedRelativePointFeature::PosedRelativePointFeature
(PosedRelativePointFeature **pcc,
 MapBank *b, long fKey, long id):
  PosedFeature(b, fKey,id)
{
  init();
  CastPtr=pcc;
  
 }
PosedRelativePointFeature::PosedRelativePointFeature
(PosedRelativePointFeature **pcc,
 RelativePointFeature *mcc, long id):
  PosedFeature(mcc->Bank, mcc->Key,id)
{
  init();
  CastPtr=pcc;
}
PosedRelativePointFeature::PosedRelativePointFeature(PosedRelativePointFeature & pw)
{
  PosedRelativePointFeature::operator=(pw);
}
void PosedRelativePointFeature::operator =(PosedRelativePointFeature & pw)
{
  init();
  equal(pw);
  CastPtr=pw.CastPtr;
  SensorPose=pw.SensorPose;

}
void PosedRelativePointFeature::init()
{
  Number2D=0;
  Number3D=0;
  NumberScalars=3;
  FullDim=3;
  Jof.reallocate(3);
  Jof=1;
  Jeo.reallocate(2,3);
  Jeo=1;
  Z.reallocate(7,1);
  Xo.reallocate(3,1);
  Jop.reallocate(3,1);
} 

int PosedRelativePointFeature::transform(Cure::Transformation3D& sensorPose,
					 unsigned short covType)
{ 

  SensorPose=sensorPose;
  SensorPose.print();
  PoseType=covType;
  RelativePointFeature *mf=getRelativePointFeature();
  if (!mf)return MAP_OBJECT_INVALID;
  Transformation3D *refPose=mf->referencePose();  
  Vector3D x;
  std::cerr<<refPose<<"\n";
  x.print();
  mf->getRelativeCartesian(x.X);
  x.print();
  refPose->invTransform(x.X,x.X);
  sensorPose.transform(x.X,Xo.Element);
  Matrix temp;
  temp.multTranspose_(Xo,Xo,1);
  Z(6,0)=sqrt(temp(0,0)-Xo(2,0)*Xo(2,0));
  double d=sqrt(temp(0,0));
  Z(5,0)=d;
  if (d<1E-15)d=1E-15;
  Z(4,0)=atan2(Xo(1,0),Xo(0,0));  
  Z(3,0)=atan2(Z(6,0),Xo(2,0));
  Z(0,0)=Xo(0,0)/d;
  Z(1,0)=Xo(1,0)/d;
  Z(2,0)=Xo(2,0)/d;
  CalcState=1;
  return 0;
}

double PosedRelativePointFeature::fineMatch(double * z,Cure::Matrix & sigma)
{
  if (!(CalcState&1))return -1;
  if ((MeasurementType==34)||(MeasurementType==35)){
    double d=Z(0,0)*z[0]+Z(1,0)*z[1]+Z(2,0)*z[2];
    d*=d;
    d=1-d;
    return d*sigma(0,0);
  }
  if ((MeasurementType==1)||(MeasurementType==2)||
      (MeasurementType==5)||(MeasurementType==7)){
    Matrix inn(3,1), a,d;
    inn(0,0)=Z(3,0)-z[0];
    if (inn(0,0)>M_PI)inn(0,0)-=2*M_PI;
    if (inn(0,0)<M_PI)inn(0,0)+=2*M_PI;
    inn(1,0)=Z(4,0)-z[1];
    if (inn(1,0)>M_PI)inn(1,0)-=2*M_PI;
    if (inn(1,0)<M_PI)inn(1,0)+=2*M_PI;
    inn(2,0)=Z(5,0)-z[2];
    a.multiply_(sigma,inn);
    d.multTranspose_(inn,a,1);
    return d(0,0);
  }
  return -1;
}
int PosedRelativePointFeature::roughMatch(double  *z,double *thresholds)
{ 
  if (!(CalcState&1))return -1;
  if ((MeasurementType==32)||(MeasurementType==33)){
    double d=Z(0,0)*z[0]+Z(1,0)*z[1]+Z(2,0)*z[2];
    d*=d;
    d=1-d;
    if (d>thresholds[0]) return 1;
    return 0;
  }
  if ((MeasurementType==1)||(MeasurementType==2)||
      (MeasurementType==5)||(MeasurementType==7)){
    
    double d=Z(3,0)-z[0];
    if (d>M_PI)d-=2*M_PI;
    if (d<M_PI)d+=2*M_PI;
    if (d>thresholds[0])return 1;
    if (-d>thresholds[0])return 1;
    d=Z(4,0)-z[1];
    if (d>M_PI)d-=2*M_PI;
    if (d<M_PI)d+=2*M_PI;
    if (d>thresholds[1])return 1;
    if (-d>thresholds[1])return 1;
    d=Z(5,0)-z[2];
    if (d>thresholds[2])return 1;
    if (-d>thresholds[2])return 1;
    return 0;
  }
  return 1;
}

//Must set SensorPose before calling this 
int PosedRelativePointFeature::calcEta(Cure::Matrix & v,
				      unsigned short mtype, 
				      unsigned short noJac)
{
  if (CalcState&1)return 0;
  setMeasurementType(mtype);
  RelativePointFeature *mf=getRelativePointFeature();
  if (!mf)return MAP_OBJECT_INVALID;
  Transformation3D *refPose,relpose;
  refPose=mf->referencePose();
  if (!refPose)return MAP_OBJECT_INVALID;

  short fulldim=3; 
  if ((MeasurementType==32)||(MeasurementType==33)) //bearing only measurements
    fulldim=2;
  double  inn[3]; 

  relpose=(refPose->inverse())+(SensorPose);
  Vector3D dx,b1,rotdx, rotb1;;
  double r1=mf->getRange();
  mf->getBearing(b1);
  relpose.getXYZ(dx.X);
  relpose.rotate(b1,rotb1);
  relpose.rotate(dx,rotdx);
  if ((MeasurementType==32)||(MeasurementType==33)){
    //bearing only measurements
    Vector3D b2;
    b2(0)=v(0,0);
    b2(1)=v(1,0);
    b2(2)=v(2,0);
    double a=b2*rotdx;
    double b=dx*b1;
    double c=b2*rotb1;
    double den=c*c;
    double num=b-a*c;
    den=1-den;
    bool rconst=false;
    double range=r1;
    if (den<1E-6)rconst=true;
    double dx2=dx*dx;
    if ((dx2-a)<1E-3)rconst=true;
    if ((dx2-b)<1E-3)rconst=true;
    if (!rconst){
      range=num/den;
      if (range<.1){
	rconst=true;
	range=0.1;
      }
    }
    if (rconst){
      if (fulldim>0)fulldim--;
    }
    Vector3D btilde(b1);
    btilde*=range;
    btilde-=dx;
    double btilde2=btilde*btilde;
    double magbtilde=sqrt(btilde2);
    Vector3D invrotb2;
    relpose.invRotate(b2,invrotb2);
    inn[0]=sqrt(2*(1-(btilde*invrotb2)/magbtilde));
    if (!rconst){
      Vector3D x1(b1);
      x1*=r1;
      x1-=dx;
      double x12=x1*x1;
      double magx1=sqrt(x12);
      inn[1]=sqrt(2*(1-(btilde*x1)/(magbtilde*magx1)));
    }
  }else {  //range and bearing measurments
    inn[0]=-v(0,0);
    inn[1]=-v(1,0);
    inn[2]=-v(2,0);
    Vector3D u(rotb1);
    u*=(r1);
    u-=rotdx;
    double u2=u*u;
    double r=sqrt(u2);
    double rho2=(u2-u(2)*u(2));
    double rho=sqrt(rho2);
    inn[0]+=atan2(rho,u(2));
    double theta=atan2(u(1),u(0));
    double adjust=1;
    if (MeasurementType&2){
      adjust=(theta*theta);
      //We should use the Measurment scale slpha 
      //scale=Mesurement->W(0,8) and alpha =W(0,9)
      inn[1]+=(mf->m_Scale*theta*(mf->m_Alpha+adjust));
      //inn[1]+=(m_W(0,8)*theta*(m_W(0,9)+adjust));
    }
    else inn[1]+=theta;
    inn[2]+=r;

  } // done range/bearing part
  Eta.reallocate(fulldim,1);
  for (short i=0; i<fulldim;i++)
    Eta(i,0)=inn[i];
  unsigned short fdim=mf->Bdual.Columns;
  if ((fdim<2)&&(fulldim<3)){
    Eta.offset(0,0,1,1);
  } 
  CalcState=1;
  return 0;
}
