// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2006 John Folkesson
//    

#include "PosedRangedPoint.hh"


using namespace Cure;

PosedRangedPoint::PosedRangedPoint(PosedRangedPoint **pcc,MapBank *b, 
				   long fKey, long id):
  PosedFeature(b, fKey,id)
{
  init();
  CastPtr=pcc;
}
PosedRangedPoint::PosedRangedPoint(PosedRangedPoint **pcc,
				       MapRangedPoint *mcc, long id):
  PosedFeature(mcc->Bank, mcc->Key,id)
{
  init();
  CastPtr=pcc;
}
PosedRangedPoint::PosedRangedPoint(PosedRangedPoint & pw)
{
  PosedRangedPoint::operator=(pw);
}
void PosedRangedPoint::operator =(PosedRangedPoint & pw)
{
  init();
  equal(pw);
  CastPtr=pw.CastPtr;
  Center=pw.Center;
}
void PosedRangedPoint::init()
{
  Number2D=0;
  Number3D=1;
  NumberScalars=0;
  FullDim=3;
  Points3D=AllocatedPoints3D;
  Points3D[0]=&Center;
  Jof.reallocate(3);
  Jeo.reallocate(3,3);
  Z.reallocate(3,1);
  Xo.reallocate(3,1);
  Jof=1;
  Jop.reallocate(3,0);
  Jev.reallocate(3);
  Jev=-1;
} 

int PosedRangedPoint::transform(Transformation3D&  pose, 
				 unsigned short covType)
{
  if (PosedFeature::transform(pose,covType))return MAP_OBJECT_INVALID;
  return calcZ();
  MapFeature *mf =getFeature();

  // Get predicted direction in the map frame from the sensor to the feature
  double s[3];
  pose.invRotate(Center.X,s);

  double d=s[0]*s[0]+s[1]*s[1]+s[2]*s[2];
  if (d<1E-25)return NOT_VISABLE;

  // Normalize s
  d=sqrt(d);
  s[0]/=d;
  s[1]/=d;
  s[2]/=d;
  int r=0;
  if (mf)
    r=mf->testVisable(s);
  if (r)return r;
  return calcZ();
}
int PosedRangedPoint::calcZ()
{
  double rho2=Center(0)*Center(0)+
    Center(2)*Center(2);
  double y2=Center(1)*Center(1);
  double r2=rho2+y2;
  double r=sqrt(r2);
  if (r<1E-2) //too close to see.
    return NOT_VISABLE;
  double rho=sqrt(rho2);
  Z(0,0)=atan2(rho,Center(1));
  Z(1,0)=atan2(Center(2),Center(0));
  Z(2,0)=r;
  if (Jeo.Rows==1)
    Jeo.offset(-2,0,3,3);
  else if (Jeo.Rows==0)
    Jeo.offset(0,0,3,3);
  Jeo(0,1)=-rho/(r2);
  Jeo(1,1)=0;
  Jeo(2,0)=Center(0)/r2;
  Jeo(2,1)=Center(1)/r2;
  Jeo(2,2)=Center(2)/r2;
  if (rho2<1e-10){
    rho2=1e-10;
    rho=1e-5;
  }
  Jeo(1,0)=-Center(2)/rho2;
  Jeo(1,2)=Center(0)/rho2;
  Jeo(0,0)=Center(0)*Center(1)/(r2*rho);
  Jeo(0,2)=Center(2)*Center(1)/(r2*rho);
  CalcState=1;
  return 0;
}



int PosedRangedPoint::roughMatch(double *z ,double *thresholds)
{
  if (!(CalcState&1))return -1;
  for(int i=0;i<3;i++)
    {
      double d=Z(i,0)-z[i];
      if ((d>thresholds[i])||(d<-thresholds[i]))return 1;
    }
  return 0;
}


int PosedRangedPoint::calcEta(Cure::Matrix & v,unsigned short mtype,
			       unsigned short noJac)
{
  // Make sure that Xo and Z have been calulated, i.e. that transform
  // was called
  if (!(CalcState&1))return 1;

  setMeasurementType(mtype);
 
  // No need to call this function twice
  if (CalcState&4)return 0;

  // Make sure that the V vector is large enough
  if (v.Rows!=3)MeasurementType=0;

  // If MeasurementType is !=0 at this point we have a map feature
  // which is initialized and we have enough data in the V vector
  if (MeasurementType==6)
    {
      Eta.reallocate(3,1);
      Eta(0,0)=Z(0,0)-v(0,0);
      Eta(1,0)=Z(1,0)-v(1,0);
      Eta(2,0)=Z(2,0)-v(2,0);
      if (Jeo.Rows==1)
	Jeo.offset(-2,0,3,3);
      else if (Jeo.Rows==0)
	Jeo.offset(0,0,3,3);
    }
  else if (MeasurementType==4) //range only
    {
      Eta.reallocate(1);
      Eta(0,0)=Z(2,0)-v(0,0);
      if (Jeo.Rows==3)
	Jeo.offset(2,0,1,3);
      else if (Jeo.Rows==0)
	Jeo.offset(2,0,1,3);
    }
  else
    {
      Eta.reallocate(0);
      if (Jeo.Rows==3)
	Jeo.offset(0,0,0,0);
      else if (Jeo.Rows==1)
	Jeo.offset(-2,0,0,0);
     
    }

  // Update the CalcState so that we do not have to perform these cals
  // again
  CalcState |= 0x4;
  return 0;
}

int PosedRangedPoint::addInfo(Match &mat,
			       Transformation3D & pose, 
			       unsigned short covType, 
			       Transformation3D & map2info,
			       const int type)
{
  MapFeature *mf =getFeature();
  if (mf)
   {
     Transformation3D sens2info;
     sens2info-=pose;
     sens2info+=map2info;
     Measurement *m=mat.Measure;
     Matrix v(1,8);
     v(0,0)=mat.PathDistance;
     v(0,1)=10;

     // Sensor pose in sensor frame 
     v(0,2)=0;
     v(0,3)=0;
     v(0,4)=0;
     // Transform sensor pose to info frame 
     sens2info.transform(v.Element+2,v.Element+2);

     // Vector that points at the detected feature in the sensor frame 
     v(0,5)=sin(m->V(0,0));
     v(0,6)=cos(m->V(0,0));
     v(0,7)=sin(m->V(1,0))*v(0,5);
     v(0,5)*=cos(m->V(1,0));
     v(0,5)*=m->V(2,0);
     v(0,6)*=m->V(2,0);
     v(0,7)*=m->V(2,0);
     // Rotate the direction vector to the info frame
     sens2info.rotate(v.Element+5,v.Element+5);
     
     // Add the information to the actual map feature
     int res=mf->addInfo(v,map2info,type);
     // Recalculate pose in the sensor frame and predict measurements
     transform(pose,covType);
     return res;
    }
  return -1;
}

int PosedRangedPoint::testVisable(double  bottomLeft[2],double topRight[2])
{
  if(!(FeatureFcn::pointInRectangle(Center.X, 
                                    bottomLeft,topRight)))return NOT_VISABLE;
  return 0;

}
