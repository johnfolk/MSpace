// = RCSID
//    $Id: MapFeature.cc ,v 1.1 2004/04/1 
//
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "PosedFeature.hh"
using namespace Cure;

void Cure::PosedFeature::equal(PosedFeature &pw)
{
  Bank=pw.Bank;
  FeatureKey=pw.FeatureKey;
  ID=pw.ID;
  PoseType=pw.PoseType;
  CalcState=pw.CalcState;
  MeasurementType=pw.MeasurementType;
  Z=pw.Z;
  Xo=pw.Xo;
  if (CalcState&4)
    {
      Eta=pw.Eta;
      Jeo=pw.Jeo;
      Jev=pw.Jev;
    }
  if (CalcState&2)
    {
      Jos=pw.Jos;
      Jof=pw.Jof;
      Jop=pw.Jop;
    }
  if (CalcState&8)
    {
      Jep=pw.Jep;
      Jes=pw.Jes;
    }
}
int Cure::PosedFeature::calcXo(Cure::Transformation3D & trans)
{
  MapFeature *mf =getMapFeature(Bank,FeatureKey);
  if (!mf){
    CalcState=0;
    return MAP_OBJECT_INVALID;
  }  
  int irow=0;
  for(int i=0;i<mf->Number3D; i++,irow+=3)
    {
      Points3D[i]->transform(trans,*mf->Points[i]);
      Points3D[i]->getXYZ(Xo.Element+irow);
    }
  for(int i=0;i<mf->Number2D; i++,irow+=2)
    {
      Points2D[i]->transform2D(trans,*mf->Points[i+mf->Number3D]);
      Points2D[i]->getXY(Xo.Element+irow);
    }
  for(int i=0;i<mf->NumberScalars; i++,irow++)
    {
      Xo(irow,0)=mf->Scalars[i];;
    }
  CalcState=1;
  return 0;
     
} 
int Cure::PosedFeature::setXo()
{
  MapFeature *mf =getMapFeature(Bank,FeatureKey);
  if (!mf){
    CalcState=0;
    return MAP_OBJECT_INVALID;
  }  
  int irow=0;
  for(int i=0;i<mf->Number3D; i++,irow+=3)
    {
      (*Points3D[i])=(*mf->Points[i]);
      Points3D[i]->getXYZ(Xo.Element+irow);
    }
  for(int i=0;i<mf->Number2D; i++,irow+=2)
    {
      *Points2D[i]=(*mf->Points[i+mf->Number3D]);
      Points2D[i]->getXY(Xo.Element+irow);
    }
  for(int i=0;i<mf->NumberScalars; i++,irow++)
    {
      Xo(irow,0)=mf->Scalars[i];;
    }
  CalcState=1;
  return 0;
     
} 
int Cure::PosedFeature::roughMatch(double *z,double *thresholds)
{
  if (!(CalcState&1))return -1;
  for (int i=0; i<Z.Rows;i++)
    {
      double d=z[i]-Z(i,0);
    if ((d>thresholds[i])||(d<-thresholds[i]))return 1;
    }
  return 0;
}

double Cure::PosedFeature::fineMatch(double *z,Cure::Matrix  & sigma)
{
  if (!(CalcState&1))return -1;
  Matrix innov(sigma.Rows,1), a, d;

  // Calculate innovation 
  for (int i=0; i<innov.Rows;i++)
    innov(i,0)=z[i]-Z(i,0);

  // Calculate mahalanobis distance (assuming that sigma is the
  // innovation covariance matrix) as d = inov^T * inv(S) * inov
  a.multiply_(sigma,innov);
  d.multTranspose_(innov,a,1);

  return d(0,0);
}


int Cure::PosedFeature::calcJo(Matrix dR[10])
{
  if (CalcState&2)return 0;
  if (!(CalcState&1))return -1;
  MapFeature *mf =getMapFeature(Bank,FeatureKey);
  if (!mf)return MAP_OBJECT_INVALID;

  // The dR array of matrices contains parts of the Jacobian. The
  // first 5 (index 0-4) are 2D and the last 5 are 3D (index 5-9).
  // The first 3 for 2D and 3D are wrt theta, phi and psi respectively
  // and the last 2 are the rotation matrix itself and the negation
  // thereof.

  // This function calculates the Jacobian of the relation
  // x_o = R*(pos(x_f)-pos(x_s))
  // where pos(x_f) and pos(x_s) are the feature and sensor positions
  // respectively in the map frame. and R is the rotation matrix from
  // the map frame to the sensor frame. Note that R depends on the
  // angles of the sensor, i.e. R=R(ang(x_s)).

  // Jos is the Jacobian wrt to the sensor pose
  // and Jof is the Jacobian wrt to the feature pose

  // Each parameterization point (3D, 2D, 0D) results 3, 2 or 1 rows
  // respectively in the Jos matrix. The number of columns depends on
  // the number of uncertain dimensions in the pose of the sensor. If
  // it is constrained to be in the plane, i.e. x,y,theta there will
  // be three columns. 

  // Derivating x_o wrt to the position part of x_s is trivial and
  // results in -R. The derivation wrt to the angles of x_s is a bit
  // more involved.

  // Initially we assume that the sensor pose has x,y,theta, i.e. dim 3
  int tcols=3;

  if (PoseType&4)tcols++;  // has z
  if (PoseType&16)tcols++; // has phi
  if (PoseType&32)tcols++; // has psi

  // Make sure that Jos is of the right size
  Jos.reallocate(Xo.Rows,tcols);

  // LOOP THROUGH 3D POINTS

  // Insert the Jacobian wrt to pos(x_s) in the first 3 columns. We
  // might not use the 3rd column (dep on z) later
  int irow=0;
  for(int i=0;i<Number3D; i++,irow+=3)
    {
      // Copy the already calculated 3x3 Jacobian into the right place
      // in the large Jos
      Jos.offset(irow,0,3,3);
      Jos=dR[9];
      Jos.reset(Xo.Rows,tcols);
     }

  // Add the columns of the Jacobian corresponding to ang(x_s)
  int row=0;
  for (int c=0; c<Number3D; c++, row+=3)
    {
      double *x=Points3D[c]->X;
      int cols=3;  
      if (!(PoseType&4)) cols--;  // known z (no uncertainty)

      // Add derivatives wrt to ang(x_s)
      int k=8;
      for (int i=0; i<3; i++,k*=2)
	{
	  if (PoseType&k)
	    {
	      for(int j=0;j<3;j++)
		Jos(j+row,cols)=dR[5+i](j,0)*x[0]+
		  dR[5+i](j,1)*x[1]+dR[5+i](j,2)*x[2];
	      cols++;
	    }
	}
    } 

  // LOOP THROUGH 2D POINTS

  // Insert the Jacobian wrt to pos(x_s) in the first 3 columns. We
  // might not use the 3rd column (dep on z) later
  if (PoseType&4) {
    // Sensor pose has uncertain z
    for(int i=0;i<Number2D; i++,irow+=2)
      {
	Jos.offset(irow,0,2,2);
	Jos=dR[4];
	Jos.reset(Xo.Rows,tcols);
	Jos(irow,2)=0;
	Jos(irow+1,2)=0;
      }
  } else {
    // Sensor pose has known z
    for(int i=0;i<Number2D; i++,irow+=2)
      {
	Jos.offset(irow,0,2,2);
	Jos=dR[4];
	Jos.reset(Xo.Rows,tcols);
      } 
  }
  // Add the columns of the Jacobian corresponding to ang(x_s)
  for (int c=0; c<Number2D; c++, row+=2)
    {
      double *x=Points2D[c]->X;
      int cols=3;  
      if (!(PoseType&4)) cols--;  // known z (no uncertainty)

      // Add derivatives wrt to ang(x_s)
      int k=8;
      for (int i=0; i<3; i++,k*=2)
	{
	  if (PoseType&k)
	    {
	      for(int j=0;j<2;j++)
		Jos(j+row,cols)=dR[i](j,0)*x[0]+
		  dR[i](j,1)*x[1];
	      cols++;
	    }
	}
    }
  // Go through the scalars which are of course all independent of the
  // sensor pose and thus give derivative 0
  for (int i=0; i<NumberScalars;irow++,i++)
    for (int j=0; j<tcols; j++)
      Jos(irow,j)=0;

  // The Jacobian Jof, i.e. d(x_o)/d(x_f) is simply R for each
  // parameter point

  // LOOP THROUGH 3D POINTS
  irow=0;
  for(int i=0;i<Number3D; i++,irow+=3)
    {
	Jof.offset(irow,irow,3,3);
	Jof=dR[8];
	Jof.reset(Xo.Rows,Xo.Rows);
    }

  // LOOP THROUGH 2D POINTS
  for(int i=0;i<Number2D; i++,irow+=2)
    {
	Jof.offset(irow,irow,2,2);
	Jof=dR[3];
	Jof.reset(Xo.Rows,Xo.Rows);
    }
  // Create Jop which is just Jof*Bdual
  Jop.multiply_(Jof,mf->Bdual);
  // Update the CalcState to tell that we have now calculates Jos, Jof
  // and Jop
  CalcState |= 0x2;

  return 0;
} 

int  Cure::PosedFeature::extend(Transformation3D & t,
			  Matrix & jaces, 
			  Matrix & invcov)
{
  MapFeature *mf=getMapFeature(Bank,FeatureKey);
  if (!mf) return 0; 
  Matrix b_n;
  int res=mf->extend(b_n, invcov); 
  if (res)
    {
      int irow=0;
      Matrix a(b_n.Columns),dR[10];
      a=0;
      if (Number3D)
	{
	  t.getJac(dR+5);
	  Matrix b;
	  b=dR[8];
	  b.invert();
	  for(int i=0;i<Number3D; i++,irow+=3)
	    {
	      a.offset(irow,3);
	      a=b;
	      a.reset(b_n.Columns,b_n.Columns);
	    }
	}
      if (Number2D)
	{
	  t.get2DJac(dR);
	  Matrix b;
	  b=dR[3];
	  b.invert();
	  for(int i=0;i<Number2D; i++,irow+=2)
	    {
	      a.offset(irow,2);
	      a=b;
	      a.reset(b_n.Columns,b_n.Columns);
	    }
	}
      calcJo(dR);
      Matrix b;
      b.multiply_(b_n,a);
      jaces.multiply_(b,Jos);
    }     
  return res;
}

int Cure::PosedFeature::calcJe()
{
  if (CalcState&0x8)return 0;
  if (!(CalcState&0x7))return -1;
  Jep.multiply_(Jeo,Jop);
  Jes.multiply_(Jeo,Jos);
  CalcState=(CalcState|0x8);
  return 0;
}
int Cure::PosedFeature::linearizeMeasurement(Cure::Matrix  & v,
				       Cure::Pose3D & pose,
				       unsigned short mtype)
{
  // Extracts the Xo vector from the point parameterization (the 3D
  // and 2D points and the scalars) and calculates the Z, i.e. the
  // predicated measurement.
  transform(pose);

  // Calculate the innovation and the Jacobians Jeo and Jev
  if (calcEta(v,mtype))return -1;

  Matrix dR[10];
  if (Number3D)
    pose.getJac(dR+5);

  if (Number2D)
    pose.get2DJac(dR);
  if (calcJo(dR))return -1;
  calcJe();
  return 0;
}
void Cure::PosedFeature::print(int level)
{
  MapFeature *mf =getMapFeature(Bank,FeatureKey);
  if (mf)mf->print(level);
}
