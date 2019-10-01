// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "PosedPointFeature.hh"


using namespace Cure;

PosedPointFeature::PosedPointFeature(PosedPointFeature **pcc,MapBank *b, long fKey, long id):
  PosedFeature(b, fKey,id)
{
  init();
  CastPtr=pcc;
  FocalLength=0;
}
PosedPointFeature::PosedPointFeature(PosedPointFeature **pcc,
				       double focal,
				       MapPointFeature *mcc, long id):
  PosedFeature(mcc->Bank, mcc->Key,id)
{
  init();
  CastPtr=pcc;
  FocalLength=focal;
}
PosedPointFeature::PosedPointFeature(PosedPointFeature & pw)
{

  PosedPointFeature::operator=(pw);
}
void PosedPointFeature::operator =(PosedPointFeature & pw)
{
  init();
  equal(pw);
  CastPtr=pw.CastPtr;
  Center=pw.Center;
  PixelCenter=pw.PixelCenter;
  FocalLength=pw.FocalLength;
}
void PosedPointFeature::init()
{
  Number2D=0;
  Number3D=1;
  NumberScalars=0;
  FullDim=3;
  Points3D=AllocatedPoints3D;
  Points3D[0]=&Center;
  Jof.reallocate(3);
  Jeo.reallocate(2,3);
  Z.reallocate(2,1);
  Xo.reallocate(3,1);
  Jof=1;
  Jop.reallocate(3,0);
} 

int PosedPointFeature::transform(Transformation3D&  pose, 
				 unsigned short covType)
{
  if (PosedFeature::transform(pose,covType))return MAP_OBJECT_INVALID;
  MapFeature *mf =getFeature();

  // Get predicted direction in the map frame from the camera to the feature
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
int PosedPointFeature::calcZ()
{
  if (Center.getPixels(PixelCenter.X,FocalLength))  
    return NOT_VISABLE;
  Z(0,0)=PixelCenter(0);
  Z(1,0)=PixelCenter(1);
  CalcState=1;
  return 0;
}
int PosedPointFeature::roughMatch(double *z ,double *thresholds)
{
  if (!(CalcState&1))return -1;
  for(int i=0;i<2;i++)
    {
      double d=Z(i,0)-z[i];
      if ((d>thresholds[i])||(d<-thresholds[i]))return 1;
    }
  return 0;
}


int PosedPointFeature::calcEta(Cure::Matrix & v,unsigned short mtype,
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
  if (MeasurementType)
    {
      // Get the Jacobian based on the focal length
      Center.getPixelJacobian(Jeo, v(2,0));

      // @todo Do we really want to do this already here, we are not
      // done yet and it is anyway done in the end??
      CalcState |= 0x4;

      // Calculate the innovation 
      Eta.reallocate(2,1);
      Eta(0,0)=Z(0,0)-v(0,0);
      Eta(1,0)=Z(1,0)-v(1,0);

      // Calculate Jacobian of the innovation wrt to the v=(u,v,f)
      // where (u,v) are the pixel measurement point and f is the
      // focal length.
      Jev.reallocate(2,v.Rows);
      Jev(0,0)=-1;
      Jev(0,1)=0;
      Jev(1,1)=-1;
      Jev(1,0)=0;      
      Jev(0,2)=Z(0,0)/FocalLength;  // Approx!
      Jev(1,2)=Z(1,0)/FocalLength;  // Approx!

      if (v.Rows==5)
	{
          if (Jev.Columns < 5) {
            std::cerr << "WARNING: Jev only has size has size "
                      << Jev.Rows << " x" << Jev.Columns << std::endl;
          }
	  Eta(0,0)-=v(3,0);
	  Eta(1,0)-=v(4,0);
	  Jev(0,3)=Jev(0,0);
	  Jev(0,4)=Jev(0,1);
	  Jev(1,3)=Jev(1,0);
	  Jev(1,4)=Jev(1,1);
	}
    }
  else  // if (MeasurementType)
    {
      Eta.reallocate(0);
      Jev.reallocate(0,v.Rows);
      Jeo.reallocate(0,3);
    }

  // Update the CalcState so that we do not have to perform these cals
  // again
  CalcState |= 0x4;

  return 0;
}

  //v= distance, w, x_i, s_i,  (0..)
/*
int PosedPointFeature::addInfo(Match &mat,
			Pose3D & pose, const int type)
{
  MapFeature *mf =getFeature();
  if (mf)
   {
     Measurement *m=mat.Measure;
     Matrix v(1,8);
     FocalLength=m->V(2,0);
     v(0,0)=mat.PathDistance;
     v(0,1)=10;
     pose.getXYZ(v.Element+2);
     v(0,5)=-m->V(0,0);//info(0,2);
     v(0,6)=m->V(2,0)/(1-m->V(2,0)/Center(1));
     v(0,7)=-m->V(1,0);//info(0,4);      
     pose.invRotate(v.Element+5,v.Element+5);
     int res=mf->addInfo(v,type);
     transform(pose);
     return res;
    }
  return -1;
}
*/
int PosedPointFeature::addInfo(Match &mat,
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
     FocalLength=m->V(2,0);
     v(0,0)=mat.PathDistance;
     v(0,1)=10;

     // Sensor pose in sensor frame 
     v(0,2)=0;
     v(0,3)=0;
     v(0,4)=0;
     // Transform sensor pose to info frame 
     sens2info.transform(v.Element+2,v.Element+2);

     // Vector that points at the detected feature in the sensor frame 
     v(0,5)=-m->V(0,0);
     v(0,6)=m->V(2,0)/(1-m->V(2,0)/Center(1));
     v(0,7)=-m->V(1,0);
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

int PosedPointFeature::testVisable(double  bottomLeft[2],double topRight[2])
{
  if (Center(1)<(FocalLength+1E-15))return NOT_VISABLE;
  if(!(FeatureFcn::pointInRectangle(PixelCenter.X, 
                                    bottomLeft,topRight)))return NOT_VISABLE;
  return 0;

}
/*
int PosedPointFeature::pixel(PixelFragment& pix,double  bottomLeft[2],
			     double topRight[2],double pixelWidth,
			     double pixelHeight)
{
  if (Center(1)<(FocalLength+1E-15))return NOT_VISABLE;
  if(!(pointInRectangle(PixelCenter.X, 
			bottomLeft,topRight)))return NOT_VISABLE;
  pix.clean();
  int c[2];
  double r[2];
  pixelize(PixelCenter.X,c,pixelWidth,pixelHeight);
  r[0]=Center(1);
  pix.add(c[0],c[1]+2,r[0]);
  pix.add(c[0]-1,c[1]+1,r[0]);
  pix.add(c[0],c[1]+1,r[0]);
  pix.add(c[0]+1,c[1]-1,r[0]);
  pix.add(c[0]-2,c[1],r[0]);
  pix.add(c[0]-1,c[1],r[0]);
  pix.add(c[0],c[1],r[0]);
  pix.add(c[0]+1,c[1],r[0]);
  pix.add(c[0]+2,c[1],r[0]);
  pix.add(c[0]-1,c[1]-1,r[0]);
  pix.add(c[0],c[1]-1,r[0]);
  pix.add(c[0]+1,c[1]-1,r[0]);
  pix.add(c[0],c[1]-2,r[0]);
  return 0;
}
*/


