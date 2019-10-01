
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "EKF.hh"
#include "GenericData.hh"


using namespace Cure;
EKF::EKF(MapBank *b, int dim){
  Bank=b;
  Keys=0;
  Check=0;
  ParameterRows=0;
  FeatureRows=0;
  Sensors=0;
  for (int i=0; i<5; i++)
    {
      SensorRows[i]=0;
      SensorCovTypes[i]=0;
      OffsetCovTypes[i]=0;
    }
  RobotCovType=0;
  RobotRows=0;
  Gain=1;
  reallocate(dim);
}
EKF::~EKF()
{
  if (Keys)
    {
      delete [] Keys;
      Keys=0;
    }
}
void EKF::reallocate(int dim)
{
  if (dim<1)dim=1;
  if (Keys)
    {
      if (dim<C.RowInc)dim=C.RowInc;

      int r=C.Rows;
      C.grow(dim);
      Dx.grow(dim,1);
      Dx.reset(r,1);
      C.reset(r,r);
      LongList *d=Keys;
      Keys=new LongList[dim];
      for (int i=0; i<C.Rows;i++)Keys[i]=d[i];
      delete [] d;
      return;
    }
  C.reallocate(dim);
  C.reset(0,0);
  Dx.reallocate(dim,1);
  Dx.reset(0,1);
  Keys=new LongList[dim];
}
void EKF::initRobotPose(Pose3D & pose)
{ 
  if (0!=C.Rows)
    {
      std::cerr
	<<"EKF: You have to initRobotPose once,before "<<
	"setting up other dimensions\n";
      return ;
    }
  if (C.RowInc<pose.Covariance.Rows)reallocate(pose.Covariance.Rows);
  Robot=pose;
  RobotCovType=Robot.getCovType();
  RobotRows=calcRows(RobotCovType);
  C.Rows=RobotRows;
  ParameterRows=C.Columns=RobotRows;
  FeatureRows=RobotRows;
  C=1E-9;
  Dx.Rows=RobotRows;
  Dx=0;
  pose.Covariance=C;
  Robot.Covariance=C;
}
int EKF::initVelocity(unsigned short veltype)
{
  Robot.setVelType(veltype);
  int r=calcRows(veltype);
  if (FeatureRows!=ParameterRows){
    std::cerr<<"Can't ste velocity type after adding features\n"; 
    return 1;
  }
  int cr=C.Rows+r;
  reallocate(cr);
  FeatureRows=cr;
  return 0;
}

int EKF::initSensorOffset(Cure::Pose3D & offsetPose, int sensorCovType)
{
   if (Sensors>4)
    {
      std::cerr<<"EKF can only handle up to 5 Sensors\n";
      return -1;
    }

  if (0==C.Rows)
    {
      std::cerr<<"EKF: You have to initRobotPose before initSensorOffset\n";
      return -1;
    }
  if (ParameterRows!=C.Rows)
    {
      std::cerr<<"EKF: You have to initSensorOffset before adding Parameters\n";
      return -1;
    }
  if (FeatureRows!=C.Rows)
    {
      std::cerr<<"EKF: You have to initSensorOffset before adding Features\n";
      return -1;
    }
  int r=C.Rows+offsetPose.Covariance.Rows;
  if (C.RowInc<r)reallocate(r);
  SensorRows[Sensors]=C.Rows;
  C.offset(SensorRows[Sensors],offsetPose.Covariance.Rows);
  C=offsetPose.Covariance;
  C.reset(r,r);
  Dx.Rows=r;
  Dx=0;
  ParameterRows=r;
  FeatureRows=r;
  Offsets[Sensors]=offsetPose;
  if (!((sensorCovType&11)==11)){
    std::cerr<<"We don't allow sensor Covariance types with less than xy and theta\n"<<sensorCovType<<" ";
    sensorCovType=(sensorCovType|11);
    std::cerr<<sensorCovType;
  }
  SensorCovTypes[Sensors]=sensorCovType;
  OffsetCovTypes[Sensors]=offsetPose.getCovType();
  Sensors++;
  return (Sensors-1);
}
void EKF::offsetPredict(Cure::Pose3D & offsetpose, int sensorIndex)
{
  Offsets[sensorIndex]=offsetpose;
  int r=C.Rows;
  int rows=calcRows(OffsetCovTypes[sensorIndex]);
  C.offset(SensorRows[sensorIndex],0,rows,r);
  C=0;
  C.reset(r,r);
  C.offset(0,SensorRows[sensorIndex],r,rows);
  C=0;
  C.reset(r,r);
  C.offset(SensorRows[sensorIndex],rows);
  C=offsetpose.Covariance;
  C.reset(r,r);
  Dx.offset(SensorRows[sensorIndex],0,rows,1);
  Dx=0;
  Dx.reset(r,1);
}
//I assume here that the acceleration is xyz and inc the Euler angles 
//The velocity type has to coorespond to the accel. update
void EKF::incrementalPredict(Pose3D & pose, Pose3D & incPart, 
			     GenericData &accPart)
{
  
  if ((accPart.ShortData.Rows>0)&&(accPart.ShortData.Columns>0))
    {
      bool u=false;
      for (int i=0;i<RobotRows;i++)
	if ((Dx(i,0)>1E-6)||(Dx(i,0)<-1E-6))u=true;
      for (int i=ParameterRows;i<FeatureRows;i++)
	if ((Dx(i,0)>1E-6)||(Dx(i,0)<-1E-6))u=true;
      if (u) updatePoses();
      unsigned short atype=Robot.getVelType();
      int arows=calcRows(atype);      
      if (arows>(FeatureRows-ParameterRows)){
	std::cerr
	  <<"You can't do an acceration predict without allocated enough ParameterRows\n";
	return;
      }
      if ((arows!=accPart.Data.Columns)||(arows!=(accPart.Data.Rows-1))){
	std::cerr
	  <<"You can't do an acceration predict without proper accPart\n";
	return;
      }
      Cure::Timestamp t=accPart.Time;
      t-=Robot.Time;
      double dt=t.getDouble();
      if (dt<1E-6)return;
      double dt2=dt*dt;
      unsigned short inctype=incPart.getCovType();
      if (inctype!=48){
	std::cerr<<"Accel predicte must use a angular velocity prediction.";
	return;
      }
      Matrix v(3,1),x(3,1),euler(3,1);
      Matrix jacR[3];
      Matrix rotinv(3);
      Matrix hessian[9];

      Robot.getXYZ(x.Element);
      Robot.getAngles(euler.Element);
      v(0,0)=Robot.Velocity[0];
      v(1,0)=Robot.Velocity[1];
      v(2,0)=Robot.Velocity[2];
      Robot.getJacobianR(jacR);
      Robot.getRinv(rotinv.Element);
      Robot.getHessianR(hessian);
      Matrix incomega(3,1);    
      incPart.getCovCoordinates(incomega.Element, inctype);
      incomega/=dt;
      Matrix covmm(6);
      covmm.offset(3,3,3,3);
      covmm=(incPart.Covariance);
      covmm/=dt2;
      covmm.offset(-3,-3,3,3);
      Matrix ar(3,1);
      Matrix absomega(3,1),atran(3),jaca[3];
      ar(0,0)=accPart.Data(0,0);
      ar(1,0)=accPart.Data(0,1);
      ar(2,0)=accPart.Data(0,2);
      accPart.Data.offset(1,0,3,3);
      covmm=(accPart.Data);
      accPart.Data.offset(-1,0,4,3);
      covmm.offset(0,0,6,6);
      if (Robot.transformedAngluarVel(incomega,absomega,atran,jaca))return;
      
      Matrix  jw[3];
      jw[0]=(jacR[0]);
      jw[0]*=(absomega(0,0));
      jw[1]=(jacR[1]);
      jw[1]*=(absomega(1,0));
      jw[2]=(jacR[2]);
      jw[2]*=(absomega(2,0));
      Matrix jdotw(jw[0]);
      jdotw+=jw[1];
      jdotw+=jw[2];
       Matrix jdotwv;
      jdotwv.multiply_(jdotw,v);
      
      Matrix rota(ar);
      rota-=jdotwv;
      Matrix a(3,1);
      a.multiply_(rotinv,rota);
      
      Matrix hwv[3];
      for (int i=0;i<3;i++){
	Matrix hw=hessian[3*i];
	hw*=absomega(0,0);
	Matrix temp(hessian[3*i+1]);
	temp*=absomega(1,0);
	hw+=temp;
	temp=hessian[3*i+2];
	temp*=absomega(2,0);
	hw+=temp;
	hwv[i].multiply_(hw,v);
      }
      Matrix jwtheta(3);
      //Matrix jev(3);
      jwtheta.offset(0,0,3,1);
      jwtheta=0;
      
      //note jaca[0] ==0 its just to not confuse people
      jwtheta.offset(0,1,3,1);
      jwtheta.multiply_(jaca[1],incomega);
      jwtheta.offset(0,1,3,1);
      jwtheta.multiply_(jaca[2],incomega);
      jwtheta.offset(0,-2,3,3);
      for (int i=1;i<3;i++)
	{
	  Matrix temp=jacR[0];
	  temp*=jwtheta(0,i);
	  Matrix temp2=jacR[1];
	  temp2*=jwtheta(1,i);
	  temp+=temp2;
	  temp2=jacR[2];
	  temp2*=jwtheta(2,i);
	  temp+=temp2;
	  hwv[i].addProduct_(temp,v);
	}
      Matrix jaeuler(3);
      jaeuler.offset(0,0,3,1);
      jaeuler.multiply_(jacR[0],rota);
      jaeuler.subtractProduct_(rotinv,hwv[0]);
      jaeuler.offset(0,1,3,1);
      jaeuler.multiply_(jacR[1],rota);
      jaeuler.subtractProduct_(rotinv,hwv[1]);
      jaeuler.offset(0,1,3,1);
      jaeuler.multiply_(jacR[2],rota);
      jaeuler.subtractProduct_(rotinv,hwv[2]);
      jaeuler.offset(0,-2,3,3);
      
      Matrix jav(3);
      jav=0;
      jav.subtractProduct_(rotinv,jdotw);
      
      Matrix jaw(3);
      
      Matrix tj(3);
      Matrix temp;
      temp.multiply_(jacR[0],v);
      tj.offset(0,0,3,1);
      tj.multiply_(rotinv,temp);
      temp.multiply_(jacR[1],v);
      tj.offset(0,1,3,1);
      tj.multiply_(rotinv,temp);
      temp.multiply_(jacR[2],v);
      tj.offset(0,1,3,1);
      tj.multiply_(rotinv,temp);
      tj.offset(0,-2,3,3);
      jaw(0,0)=-(atran(0,0)*tj(0,0)+atran(1,0)*tj(0,1)+atran(2,0)*tj(0,2));
      jaw(1,0)=-(atran(0,0)*tj(1,0)+atran(1,0)*tj(1,1)+atran(2,0)*tj(1,2));
      jaw(2,0)=-(atran(0,0)*tj(2,0)+atran(1,0)*tj(2,1)+atran(2,0)*tj(2,2));
      jaw(0,1)=-(atran(0,1)*tj(0,0)+atran(1,1)*tj(0,1)+atran(2,1)*tj(0,2));
      jaw(1,1)=-(atran(0,1)*tj(1,0)+atran(1,1)*tj(1,1)+atran(2,1)*tj(1,2));
      jaw(2,1)=-(atran(0,1)*tj(2,0)+atran(1,1)*tj(2,1)+atran(2,1)*tj(2,2));
      jaw(0,2)=-(atran(0,2)*tj(0,0)+atran(1,2)*tj(0,1)+atran(2,2)*tj(0,2));
      jaw(1,2)=-(atran(0,2)*tj(1,0)+atran(1,2)*tj(1,1)+atran(2,2)*tj(1,2));
      jaw(2,2)=-(atran(0,2)*tj(2,0)+atran(1,2)*tj(2,1)+atran(2,2)*tj(2,2));
      double newx[9];
      newx[0]=x(0,0)+v(0,0)*dt+a(0,0)*dt2/2;
      newx[1]=x(1,0)+v(1,0)*dt+a(1,0)*dt2/2;
      newx[2]=x(2,0)+v(2,0)*dt+a(2,0)*dt2/2;
      newx[3]=euler(0,0)+absomega(0,0)*dt;
      newx[4]=euler(1,0)+absomega(1,0)*dt;
      newx[5]=euler(2,0)+absomega(2,0)*dt;
      newx[6]=v(0,0)+a(0,0)*dt;
      newx[7]=v(1,0)+a(1,0)*dt;
      newx[8]=v(2,0)+a(2,0)*dt;
      Robot.setXYZ(newx);
      Robot.setAngles(newx+3);
      Robot.Velocity[0]=newx[6];
      Robot.Velocity[1]=newx[7];
      Robot.Velocity[2]=newx[8];
      
      int rows=ParameterRows+3;
      Matrix j(rows,rows);
      j=0;
      double d=(dt2/2);
      j.offset(0,3,3,3);
      j=jaeuler;
      j*=d;
      j.offset(0,rows-6,3,3);
      
      j=jav;
      j*=d;
      j(0,0)+=dt;
      j(1,1)+=dt;
      j(2,2)+=dt;
      j.offset(3,-rows+6,3,3);
      j=jwtheta;
      j*=dt;
      j.offset(rows-6,0,3,3);
      j=jaeuler;	
      j*=dt;
      j.offset(0,rows-6,3,3);
      j=jav;
      j*=dt;
      j.offset(6-rows,6-rows,rows,rows);

      Matrix k(rows,6);
      k.offset(0,0,3,3);
      
      k=rotinv;
      k*=d;
      k.offset(0,3,3,3);
      k=jaw;
      k*=d;
      k.offset(3,0,3,3);
      k=atran;
      k*=dt;
      k.offset(rows-6,-3,3,3);
      k=rotinv;
      k*=dt;
      k.offset(0,3,3,3);
      k=jaw;
      k*=dt;
      k.offset(6-rows,-3,rows,6);
      Matrix cm; 
      temp.multiply_(k,covmm);
      cm.multTranspose(temp,k,2);
      long crows=C.Rows;
      C.offset(0,0,rows,crows);
      Matrix jc;
      jc.multiply_(j,C);
      jc.offset(0,0,rows,rows);
      temp.multTranspose(jc,j,2);
      temp+=cm;
      jc.offset(0,0,rows,crows);
      C+=jc;
      C.offset(0,0,rows,rows);
      C+=temp;
      temp.transpose(jc);
      C.offset(0,0,crows,rows);
      C+=temp;
      C.offset(0,0,crows,crows);
      setRobotCov();
      Robot.Time=accPart.Time;
      pose=Robot;
      return;
    }
  //Did not include the w cov int cova and the jac of the Jrot in the jvx
  return incrementalPredict(pose,incPart);
}
void EKF::incrementalPredict(Pose3D & pose, Pose3D & inc)
{
  if (RobotCovType!=pose.getCovType())
    {
      std::cerr<<"EKF: IncrementalPredict must use pose covariance type: "
	       <<RobotCovType<<"\n";
      return ;
    }
  checkKeys(0);
  Matrix jac;
  inc.Covariance*=Gain;
    Timestamp t=inc.Time;
  t-=Robot.Time;
  double dt=t.getDouble();

  if (inc.Time>Robot.Time)
    Robot.Time=inc.Time;
  double x[6];
  Robot.getCoordinates(x);
  Robot.add(Robot,inc,&jac);
  // Remember how many rows the covariance matrix has (and thus columns)
  int r=C.Rows;

  if (Robot.getVelType()){
    if (dt>1e-6){
      dt*=dt;
      double dx[6];
      Robot.getCoordinates(dx);
      for (int i=0;i<6;i++)
	Robot.Velocity[i]=(dx[i]-x[i])/dt;
      Matrix jacV;
      Robot.separateVelocity(jac,jacV,jac,Robot, inc);
      int vr=calcRows(Robot.getVelType());
      C.offset(ParameterRows,0,vr,r);
      C=0;
      C.offset(-ParameterRows,ParameterRows,r,vr);
      C=0;
      C.offset(ParameterRows,0,
	       vr,vr);
      Matrix temp;
      temp.multTranspose_(inc.Covariance,jacV,2);
      C.multiply_(jacV,temp);
      C/=dt;
      C.offset(-ParameterRows,r);
    }
  }
  //probably should set part of C.

  // Get only the jac wrt robot state variables
  jac.Columns=RobotRows;


#if 1
  Matrix temp(RobotRows,r);
  C.offset(0, 0,RobotRows,r);
  temp.multiply_(jac,C);
  temp.offset(0,RobotRows);
  temp=Robot.Covariance;
  temp.reset(RobotRows,r);
  C=temp;
  C.reset(r,r);
  C.offset(0,r,RobotRows);
  C.transpose_(temp);
  C.reset(r,r);
#else

  // An alternative implementation to the above with a bit more
  // comments. It seems however that I have missunderstodd something
  // because it does not qive exactly the same results. Could it be
  // just numerical effects? When I run through HolosTest the final
  // pose for with the old code is

  // Robot pose: ClassType:0 SubType:704 ID:0 Packed:0 vtype:0 Ptype:11
  // Pos:[0.111109 -1.46696e-05 2.58682e-14] Rot:[0.00926754 -1.97799e-14 3.28765e-15]  v:[0 0 0 0 0 0 ]
  // P=[2.10401e-05 -7.94206e-07 6.44575e-06
  //    -7.94206e-07 1.69544e-05 3.70567e-07
  //    6.44575e-06 3.70567e-07 2.07527e-05 ]; %Rows/cols/RowInc:3 3 3

  // and with the new

  // Robot pose: ClassType:0 SubType:704 ID:0 Packed:0 vtype:0 Ptype:11
  //  Pos:[0.111212 0.00208734 2.58682e-14] Rot:[0.00887191 -1.97799e-14 3.28765e-15]  v:[0 0 0 0 0 0 ]
  // P=[2.14096e-05 -8.29394e-07 7.51571e-06
  //    -8.29394e-07 1.79003e-05 5.8945e-07
  //    7.51571e-06 5.8945e-07 2.42295e-05 ]; %Rows/cols/RowInc:3 3 3

  // Get the submatrix to hold the part of the covariance
  // matrix that contains the correlation ot the robot pose
  C.offset(0, RobotRows, RobotRows, r - RobotRows);

  // Create temporary storage
  Matrix temp(RobotRows, r - RobotRows);
  
  // Perform temp = Fr * Crf, where Fr is the jacobian of the
  // prediction function wrt to the robot state variables and Crf is
  // the correlation to the robot state.
  temp.multiply_(jac,C);

  // Copy the results into the covariance matrix
  
  // Crf: correlation part submatrix (upper rows)
  C = temp;  

  // Cfr: correlation part submatrix (left columns)
  C.reset(r,r);
  C.offset(RobotRows, 0, r-RobotRows, RobotRows);
  C.transpose_(temp);

  // Crr: robot state covariance
  C.reset(r,r);
  C.offset(0,RobotRows);
  C = Robot.Covariance;

  // Make C be the full covariance matrix again
  C.reset(r,r);
#endif

  pose=Robot;
}

int EKF::checkEigen()
{
  Matrix lambda(C.Rows);
  C.symmetricEigen(lambda);
  int test=0;
  std::cerr<<"EIGEN VALUES FOR C MATRIX:"<<std::endl;
  for (int i=0; i<C.Rows;i++)
    {
      if (lambda(i,i)<1E-16)test=1;
      std::cerr<<lambda(i,i)<<" ";
    }
  std::cerr<<std::endl;
  return test;
}
int EKF::update(Match &mat, double mahalanobisTest,
		Cure::Pose3D & pose)
{
  PosedFeature &pf=*mat.MatchedFeature;
  long *index=pf.getIndex();
  if (!index)return -1;

  // Since this function might be called when looping through a set of
  // measurements and thus this might not be the first call, we need
  // to propagate the change in P-space coordinates to the "real
  // coordinates" of this particular feature and possibly change the
  // projection matrix before we can perform the update with the new
  // measurement.
  updateP(pf);

  //Robot=pose;

  Measurement *m=mat.Measure;

  // Count the number of P-dimensions for this feature
  int dim=0;
  while ((index[dim])&&(dim<pf.FullDim))dim++;

  // Perform the update if there are P-dimensions
  if (dim)
    {
      Matrix cov;
      
      // Calculate the innovation and all Jacobians
      pf.linearizeMeasurement( m->V, Robot,m->MeasurementType);
      if (pf.Eta.Rows==0) return 0;

      // Calculate measurement covariance R=Jev*Cvv*Jev^T
      Matrix invcov;
      invcov.multTranspose_(m->CovV,pf.Jev,2);	        
      cov.multiply_(pf.Jev,invcov);

      Matrix **(jacobians)[1];
      Matrix * jac[2];
      jacobians[0]=jac;
      jac[0]=&(pf.Jes);
      jac[1]=&(pf.Jep);

      // We need to keep track of the indices of all variables that
      // our measurement depends on since these will give non-zero
      // columns if a full Jacobian were to be calculated. The
      // linearizeMeasurement call above have given us the non-zero
      // parts but we need to define where they fit in when we do the
      // matrix multiplications.
      long **(indecies)[1];
      long *(ind)[2];
      indecies[0]=ind;

      // Create an array of indices for the robot pose variables
      long poseindex[RobotRows];
      indecies[0][0]=poseindex;
      for (long i=0; i<RobotRows; i++)
	poseindex[i]=i;

      // Set the indices for the map feature that we have measured
      indecies[0][1]=index;

      Matrix Sinv,A,w; 
      // Calculate A=J*C and Sinv=(innovation covariance)^(-1)
      getMetric(A, Sinv,cov,jacobians,indecies, 1,2);

      // Check if we are performing a mahalanobis test
      if (mahalanobisTest>0)
	{
          Matrix temp0,temp1;
	  temp0.multiply_(Sinv,pf.Eta);
	  temp1.multTranspose_(pf.Eta,temp0,1);
	  if (temp1(0,0)>mahalanobisTest)return 1;
	}

      // Calculate Kalman gain
      w.multTranspose_(A,Sinv,1);

      // Update state covariance matrix C-=w*A=A^T*Sinv*A, where A=J*C
      C.subtractProduct_(w,A);

      // Calculate the change in M-space coordinates
      Dx.subtractProduct_(w,pf.Eta);

      // Propagate the change in robot pose to the absolute robot
      // pose, also update the robot pose covariance matrix
      updatePoses();

      // Update the feature coordinates based on the change in P-space
      // coordinates. Note that we 
      updateP(pf); 
    }
  setRobotCov();
  pose=Robot;
  return 0;

}
int EKF::updateP(PosedFeature & pf)
{
  long *index=pf.getIndex();
  if (!index)return -1;
  return pf.updateIndexedP(Dx.Element);
}
MapFeature * EKF::getFeature(int i)
{
  
  MapFeature *mf=getMapFeature(Bank,Keys[i].Element);
  if (mf) return mf;
  long key=0;
  while ((Keys[i].Next)&&(key==0))
    {
      mf=getMapFeature(Bank,Keys[i].Element);
      if (mf)
	key=Keys[i].Element;
      else Keys[i].remove(0);
    }
  return mf;
}
void EKF::updateP()
{

  int i=FeatureRows;
  while(i<C.Rows)
    {
      MapFeature *mf=getFeature(i);
      if (mf){
	long j=mf->updateIndexedP(Dx.Element);
	if (j)i+=j;
	else
	  {
	    j=Keys[i].Element;
	    int t=C.Rows;
	    while ((i<t)&&(Keys[i].Element==j))i++;
	  }
      }
      else removeRows(i,1);
    }
}

int EKF::removeRows(int start,int dim)
{
  for (int i=0; i<dim;i++)
    {
      C.deleteRow(start);
      C.deleteColumn(start);
      Dx.deleteRow(start);
    }  
  for (int i=start; i<C.Rows; i++)  
    {
      Keys[i]=Keys[i+dim];
      for (LongList *llist=Keys+i; llist->Next; llist=llist->Next)
	{
	  MapFeature *mf=getMapFeature(Bank,llist->Element);
	  if (mf)mf->changeIndex(i+dim,i);
	}
    }
  int top1=C.Rows+dim;
  for(int i=C.Rows; i<top1;i++)
    Keys[i].clear();
  return C.Rows;
}
void EKF::checkKeys(int start)
{
  if (start<FeatureRows)start=FeatureRows;
  int i=start; 
  while (i<C.Rows)
    {

      MapFeature *mf=getFeature(i);
      if (mf) {
	long ky=Keys[i].Element;
	while ((i<C.Rows)&&(Keys[i].Element==ky))i++;
      }
      else removeRows(i,1);
    }
}
    
/**
 *	do S= A*(jac1,jac2)^T + cov
 *	(xs,xp)->(xs,xp)+W*de
 *	C->C-W*A =(I-WJ)C=(I-A^TS^-1J)C
 *	A=(0,...0,jac1,0,...0,jac2,0...0)*C
 *	W=A^T*S^-1
 *
 *
 * @param A =jacobians*C
 * @param Sinv =Inverse(A*jacobians^T+cov)
 * @param cov The covariance of the measurments
 * @param jacobians[i][m] i runs 1 to numblocks m and runs from 1 to numJac
 *        for each i all the jacobians have the same number of rows while the
 *        columns are indexed int C by index. 
 * @param index  The index[i][m][j] gives the C row of the jth 
 *        column of jacobian[i][m].
 * @param numblocks This is the number of measurement blocks,
 * @param numJac  This is the number of jacobians in each block.
 */
void EKF::getMetric(Matrix & A, Matrix & Sinv, 
		    Matrix & cov, Matrix *** jacobians,
		    long *** index,  int numblocks,int numJac)
{
  int etadim=0;
  for (int i=0; i<numblocks;i++)
    etadim+=jacobians[i][0]->Rows;
  A.reallocate(etadim,C.Rows);
  A=0;
  for (int i=0; i<numblocks;i++)
    {
      Matrix *jac=jacobians[i][0];
      A.Rows=jac->Rows;
      for (int m=0; m<numJac; m++)
	{
	  long *indx=index[i][m];
	  jac=jacobians[i][m];
	  int oldCol=jac->Columns;
	  int k=0;
	  while (k<oldCol)
	    {
	      int rowoff=indx[k]*C.RowInc;
	      C.Element+=rowoff;
	      int col=1;
	      int ii=indx[k]+1;
	      for (int j=k+1;((j<oldCol)&&(indx[j]==ii)); ii++,j++)col++;
	      jac->Columns=col;
	      C.Rows=col; 
	      A.addProduct_(*jac,C);
	      C.Element-=rowoff;
	      jac->Element+=col;
	      k+=col;
	    }
	  jac->Columns=oldCol;
	  jac->Element-=oldCol;
	}
      A.Element+=(A.RowInc*jac->Rows);
    }
  A.Rows=etadim;
  A.Element-=etadim*A.RowInc;
  C.Rows=C.Columns;

  //C<->Atranspose
  //A<->Sinv
  Sinv.reallocate(etadim);
  Sinv=0;
  Matrix b;
  for (int i=0; i<numblocks;i++)
    {
      Matrix *jac=jacobians[i][0];
      Sinv.Rows=jac->Rows;
      for (int m=0; m<numJac; m++)
	{
	  long *indx=index[i][m];
	  jac=jacobians[i][m];
	  int oldCol=jac->Columns;
	  int k=0;
	  while (k<oldCol)
	    {
	      A.Element+=indx[k];
	      int col=1;
	      int ii=indx[k]+1;
	      for (int j=k+1;((j<oldCol)&&(indx[j]==ii)); ii++,j++)col++;
	      jac->Columns=col;
	      A.Columns=col; 
	      b.multTranspose_(*jac,A,2);
	      Sinv+=b;
	      A.Element-=indx[k];
	      jac->Element+=col;
	      k+=col;
	    }
	  jac->Columns=oldCol;
	  jac->Element-=oldCol;
	}
      Sinv.Element+=(Sinv.RowInc*jac->Rows);
    }
  Sinv.Rows=etadim;
  Sinv.Element-=etadim*Sinv.RowInc;
  A.Columns=C.Columns;
  Sinv+=cov;
  Sinv.symetrize();
  Sinv.invert();
}
int EKF::extend(Match &mat)
{ 
  PosedFeature &pf=*mat.MatchedFeature;
  int dim=extend(pf);
  setRobotCov();
  return dim;
}
int EKF::extend(Pose3D &pose, PosedFeature &pf)
{
  int dim=extend(pf);
   pose=Robot;
  return dim;
}
int EKF::extend(PosedFeature &pf)
{ 
  long *index=pf.getIndex();
  if (!index)return -1; 
  updateP(pf);
  int dim=0;
  while ((index[dim])&&(dim<pf.FullDim))dim++;
  Matrix jaces,invcov;
  int  res=pf.extend(Robot,jaces,invcov);
  if (res)
    {
      if (res<0)return res;
      insertRows(dim,res,index, pf.FeatureKey);
      dim+=res;
      int top=C.Rows-dim;      
      MapFeature *mf=getMapFeature(Bank,pf.FeatureKey);
      if (mf)mf->setIndex(top);
      invcov.invert();
      Matrix r_n(res,res);
      r_n=1;
      Matrix **(jacobians)[1];
      Matrix * jac[2];
      jacobians[0]=jac;
      long **(indecies)[1];
      long *(ind)[2];
      indecies[0]=ind;
      jac[0]=&(jaces);
      jac[1]=&r_n;
      long poseindex[RobotRows];
      indecies[0][0]=poseindex;
      for (long i=0; i<RobotRows; i++)
	poseindex[i]=i;
      indecies[0][1]=index+dim-res;
      Cure::Matrix  w,A, Sinv;
      getMetric(A, Sinv,invcov,jacobians,indecies, 1,2);
      w.multTranspose_(A,Sinv,1);
      C.subtractProduct_(w,A);
      setRobotCov();
    }
  return dim;
}
int EKF::extend(PosedFeature &pf,int sensorIndex )
{
  long *index=pf.getIndex();
  if (!index)return -1; 
  updateP(pf);
  int dim=0;
  while ((index[dim])&&(dim<pf.FullDim))dim++;
  Matrix jaces,invcov;
  int cols=(OffsetCovTypes[sensorIndex]<<6);
  cols+=RobotCovType;
  Transformation3D sensorpose;
  Matrix jsro(6,12);
  sensorpose.doAplusB(Robot,Offsets[sensorIndex], jsro,
		      SensorCovTypes[sensorIndex], cols, 1);
  int  res=pf.extend(sensorpose,jaces,invcov);
  if (res)
    {
      if (res<0)return res;
      insertRows(dim,res,index, pf.FeatureKey);
      dim+=res;
      int top=C.Rows-dim;      
      MapFeature *mf=getMapFeature(Bank,pf.FeatureKey);
      if (mf)mf->setIndex(top);
      invcov.invert();
      Matrix r_n(res,res);
      r_n=1;
      Matrix **(jacobians)[1];
      Matrix * jac[3];
      jacobians[0]=jac;
      long **(indecies)[1];
      long *(ind)[3];
      indecies[0]=ind;   
      Matrix jer,jeoff;
      jsro.offset(0,0,jsro.Rows,RobotRows);
      jer.multiply_(jaces,jsro);
      jsro.offset(0,RobotRows,jsro.Rows,calcRows(OffsetCovTypes[sensorIndex]));
      jeoff.multiply_(jaces,jsro);
      jac[0]=&jer;
      jac[1]=&jeoff;
      jac[2]=&r_n;
      long poseindex[RobotRows];
      indecies[0][0]=poseindex;
      for (long i=0; i<RobotRows; i++)
	poseindex[i]=i;
      long offindex[jeoff.Columns];
      indecies[0][1]=offindex;
      for (int i=0; i<jeoff.Columns; i++)
	offindex[i]=SensorRows[sensorIndex]+i;
      indecies[0][2]=index+dim-res;
      Cure::Matrix  w,A, Sinv; 
      getMetric(A, Sinv,invcov,jacobians,indecies, 1,3);
      w.multTranspose_(A,Sinv,1);
      C.subtractProduct_(w,A);
      setRobotCov();
    }
  return dim;
}
  /**
   * constraint is a*(pl[0].P,pl[1].P)=b 
   * dX=Kb;
   * C'=(I-Ka)C(I-Ka)
   * K=Ca^T(aCa^T)^-1
   * 
   */
int EKF::merge(PosedFeature *pf[2],Matrix &a,Matrix &e,Matrix & cov,
		double mahalanobisTest, Pose3D &pose)
{
  long *index0=pf[0]->getIndex();
  long *index1=pf[1]->getIndex();
  if (!index0)return 1;
  if (!index1)return 1;
  int dim0=0;
  while ((index0[dim0])&&(dim0<pf[0]->FullDim))dim0++;
  int dim1=0;
  while ((index1[dim1])&&(dim1<pf[1]->FullDim))dim1++;
  // Robot=pose;
 if (dim0)
    {
      if (a.Columns!=(dim0+dim1))
	{
	  std::cerr<<"EKF::merge A has wrong columns";
	  a.print();
	  return 0;;
	}
      Matrix  jac1, jac0;
      a.Columns=dim0;
      jac0=a;
      a.Columns=dim1;
      a.Element+=dim0;
      jac1=a;
      a.Columns+=dim0;
      a.Element-=dim0;
      Matrix Sinv,A,w;
      Matrix **(jacobians)[1];
      Matrix * jac[2];
      jacobians[0]=jac;
      jacobians[0][0]=&jac0;
      jacobians[0][1]=&jac1;
      long **(indecies)[1];
      long *(ind)[2];
      indecies[0]=ind;
      indecies[0][0]=index0;
      indecies[0][1]=index1;
      getMetric(A, Sinv,cov,jacobians,indecies, 1,2);
      if (mahalanobisTest>0)
	{
	  jac0.multiply_(Sinv,e);
	  jac1.multTranspose(e,jac0,1);
	  if ((jac1(0,0)>(mahalanobisTest*e.Rows)))return 0;
	}
      w.multTranspose_(A,Sinv,1);
      C.subtractProduct_(w,A);
      Dx.subtractProduct_(w,e);
      updatePoses();
      pf[0]->transform(Robot);
      pf[1]->transform(Robot);
    }
 setRobotCov();
 updateP(*pf[0]);
 updateP(*pf[1]);
 pose=Robot;
 return 1; 
}

int EKF::updateSensor(Match &mat, double mahalanobisTest, int sensorIndex)
{
  PosedFeature &pf=*mat.MatchedFeature;
  long *index=pf.getIndex();
  if (!index)return -1;
  Measurement *m=mat.Measure;
  //Check if this has been extended without adding rows to C
  if (index[0]==-1){
    MapFeature *mf=getMapFeature(Bank,pf.FeatureKey);
    if (!mf)return -1;
    int res=mf->Bdual.Columns;
    insertRows(0,res,index, pf.FeatureKey);
    mf->setIndex(C.Rows-res); 
    for (int i=C.Rows-res; i<C.Rows; i++)
      C(i,i)=100;
  }
  updateP(pf);
  int dim=0;
  while ((index[dim])&&(dim<pf.FullDim))dim++;
  if (dim)
    {
      Matrix Sinv,A,w, temp0,temp1;
      Matrix **(jacobians)[1];
      Matrix * jac[3];
      jacobians[0]=jac;
      long **(indecies)[1];
      long *(ind)[3];
      indecies[0]=ind;   
      int cols=(OffsetCovTypes[sensorIndex]<<6);
      cols+=RobotCovType;
      Pose3D sensorpose;
      Matrix jsro(6,12);
      sensorpose.setCovType(SensorCovTypes[sensorIndex]);
      sensorpose.doAplusB(Robot,Offsets[sensorIndex], jsro,
			  SensorCovTypes[sensorIndex], cols, 1);
      Matrix cov;
      pf.linearizeMeasurement( m->V, sensorpose,m->MeasurementType);
      if (pf.Eta.Rows==0) return 0;
      Matrix invcov;
      invcov.multTranspose_(m->CovV,pf.Jev,2);	        
      cov.multiply_(pf.Jev,invcov);
      Matrix jer,jeoff;
      jsro.offset(0,0,jsro.Rows,RobotRows);
      jer.multiply_(pf.Jes,jsro);
      jsro.offset(0,RobotRows,jsro.Rows,calcRows(OffsetCovTypes[sensorIndex]));
      jeoff.multiply_(pf.Jes,jsro);
      jac[0]=&jer;
      jac[1]=&jeoff;
      jac[2]=&pf.Jep;
      long poseindex[RobotRows];
      indecies[0][0]=poseindex;
      for (long i=0; i<RobotRows; i++)
	poseindex[i]=i;
      long offindex[jeoff.Columns];
      indecies[0][1]=offindex;
      for (int i=0; i<jeoff.Columns; i++)
	offindex[i]=SensorRows[sensorIndex]+i;
      indecies[0][2]=index;
      getMetric(A, Sinv,cov,jacobians,indecies,1,3);
      if (mahalanobisTest>0)
	{
	  temp0.multiply_(Sinv,pf.Eta);
	  temp1.multTranspose_(pf.Eta,temp0,1);
	  if (temp1(0,0)>mahalanobisTest)return 1;
	}
      w.multTranspose_(A,Sinv,1);
      C.subtractProduct_(w,A);
      Dx.subtractProduct_(w,pf.Eta);
      updatePoses();
      updateP(pf); 
    }
  return 0;
}

void EKF::updatePoses()
{
  double x[6];
  Robot.getCovCoordinates(x,RobotCovType);
  int i=0;
  while(i<RobotRows)
    {
      x[i]+=Dx.Element[i];
      Dx.Element[i]=0;
      i++;
    }
  Robot.setCovCoordinates(x,RobotCovType);
  int r=C.Rows;
  C.offset(0,RobotRows);
  Robot.Covariance=C;
  C.offset(RobotRows,0);
  for (int k=0; k<Sensors; k++)
    {
      int top=calcRows(OffsetCovTypes[k]);
      Offsets[k].getCovCoordinates(x,OffsetCovTypes[k]);
      C.offset(0,top);
      for (int j=0; j<top; j++,i++)
	{
	  x[j]+=Dx.Element[i];
	  Dx.Element[i]=0;
	
	}
      Offsets[k].setCovCoordinates(x,OffsetCovTypes[k]);
      Offsets[k].Covariance=C;
      C.offset(top,0);
    }
  C.reset(r,r);
  int j=0;
  unsigned short t=Robot.getVelType();
  while ((t)&&(i<FeatureRows)){
    if (t&1){
      Robot.Velocity[j]+=Dx.Element[i];
      Dx.Element[i]=0;
      i++;
      j++;
    }
    t=(t>>1);
  }
}

void  EKF::updateSensor(Match *matches, int n, Pose3D &probot,
			Transformation3D & map2info, 
			double manhanalobisTest,int sensorIndex)
{
  if (sensorIndex<0)return update(matches,n,probot,
				 map2info, 
				 manhanalobisTest);
  double best=100000;
  for (int i=0; i<n; i++)
    {
      if (matches[i].MatchDistance>0)
	if (matches[i].MatchDistance<best)
	  best=matches[i].MatchDistance;
    }
  Transformation3D startpose=Robot+Offsets[sensorIndex];
  Transformation3D psensor=startpose;
  best*=1.1;
  for (int i=0; i<n; i++)
    if (matches[i].MatchDistance>-1)
      if (matches[i].MatchDistance<best)
	{
	  if(updateSensor(matches[i],manhanalobisTest,sensorIndex)==0)
	    {
	      psensor=Robot;
	      psensor+=Offsets[sensorIndex];
	      map2info=((psensor-startpose)+map2info);
	      startpose=psensor;
	      addinfo(matches[i],psensor,SensorCovTypes[sensorIndex],map2info);
	      extend(matches[i],sensorIndex);
	    }
	  matches[i].clear();  //this will delete the PosedFeature
	}  
  for (int i=0; i<n; i++)
    {
      if (matches[i].MatchDistance>-1)
	{
	  if(updateSensor(matches[i],manhanalobisTest,sensorIndex)==0)
	    {	
	      psensor=Robot;
	      psensor+=Offsets[sensorIndex];
	      map2info=((psensor-startpose)+map2info);
	      startpose=psensor;
	      addinfo(matches[i],psensor,SensorCovTypes[sensorIndex],map2info);
	      extend(matches[i],sensorIndex);
	    }
	}  
      matches[i].clear();//this will delete the PosedFeature
    }
  probot=Robot;
}

void  EKF::update(Match *matches, int n, Pose3D &psensor,
		    Transformation3D & map2info, 
		    double manhanalobisTest)
{
  //Robot=psensor;

  // Look for the best match
  double best=100000;
  for (int i=0; i<n; i++)
    {
      if (matches[i].MatchDistance>0)
	if (matches[i].MatchDistance<best)
	  best=matches[i].MatchDistance;
    }

  // Variable used to be able to correctly update the map2info
  // tranformation as we update the robot pose in world coordinates.
  Transformation3D startpose=Robot;

  // Perform update using the matches that are close to best
  // (incl. the best one).
  best*=1.1;;
  for (int i=0; i<n; i++)
    if (matches[i].MatchDistance>-1)
      if (matches[i].MatchDistance<best)
	{
	  if(update(matches[i],manhanalobisTest,Robot)==0)
	    {
	      map2info=((Robot-startpose)+map2info);
	      startpose=Robot;
	      addinfo(matches[i],map2info);
	      extend(matches[i]);
	    }
	  matches[i].clear();  //this will delete the PosedFeature
	}  

  // Update with the rest of the features
  for (int i=0; i<n; i++)
    {
      if (matches[i].MatchDistance>-1)
	{
	  if(update(matches[i],manhanalobisTest,Robot)==0)
	    {	
	      map2info=((Robot-startpose)+map2info);
	      startpose=Robot;
	      addinfo(matches[i],map2info);
	      extend(matches[i]);
	    }
	}  
      matches[i].clear();//this will delete the PosedFeature
    }

  // NOTE that this assumes that this function assumes that the robot
  // pose is equal to the sensor pose before this function is
  // called. This is typically achieved with incrementalPredict.

  // Return the updated position
  psensor=Robot;
}
void  EKF::update(Match *matches, int n, Pose3D &psensor, 
		    double manhanalobisTest)
{
  //Robot=psensor;
  double best=100000;
  for (int i=0; i<n; i++)
    {
      if (matches[i].MatchDistance>0)
	if (matches[i].MatchDistance<best)
	  best=matches[i].MatchDistance;
    }
  Transformation3D map2info;
  best*=1.1;;

  for (int i=0; i<n; i++)
    if (matches[i].MatchDistance>-1)
      if (matches[i].MatchDistance<best)
	{
	  if(update(matches[i],manhanalobisTest,Robot)==0)
	    {
	      addinfo(matches[i],map2info);
	      extend(matches[i]);
	    }
	  matches[i].clear();  //this will delete the PosedFeature
	}  
  for (int i=0; i<n; i++)
    {
      if (matches[i].MatchDistance>-1)
	{
	  if(update(matches[i],manhanalobisTest,Robot)==0)
	    {	
	      addinfo(matches[i],map2info);
	      extend(matches[i]);
	    }
	}  
      matches[i].clear();//this will delete the PosedFeature
    }
  psensor=Robot;
}
void EKF::permuteRows(Cure::Matrix & permute,PosedFeature **pf,long number)
{
  updateP();
  long index[permute.Rows];
  int i=0;  
  for (long j=0; j<number; j++)
    {
      long *index1=pf[j]->getIndex();
      if (!index1)return;
      long k=0;
      while ((index1[k])&&(k<pf[j]->FullDim)){
	index[i]=index1[k];  
	i++;
	k++;
      }
    }
  if (i!=permute.Rows)
    {
       for (int j=0; j<number; j++)
	 pf[j]->print();
       permute.print();
	std::cerr<<"EKF ERROR PERMUTE ROWS "<<i<<std::endl;
      for (int j=0; j<i; j++)
	std::cerr<<index[j]<<" ";
    }
  permuteRows(permute,index);
  //  permute.invert();
  //permuteRows(permute,index);
}
void EKF::permuteRows(Cure::Matrix & permute,long *index)
{
  Matrix temp(permute.Rows,C.Rows);
  temp=0; 
  for (int i=0; i<permute.Rows; i++)
    for (int k=0; k<C.Rows; k++)
      for (int j=0;j<permute.Rows;j++)
	temp(i,k)+=permute(i,j)*C(index[j],k);
  Matrix temp2(permute.Rows);
  temp2=0;
  for (int i=0; i<permute.Rows; i++)
    for (int k=0; k<permute.Rows; k++)
      for (int j=0;j<permute.Rows;j++)
	temp2(i,k)+=temp(i,index[j])*permute(k,j);
  for (int i=0; i<permute.Rows; i++)
    for (int k=0; k<C.Rows; k++)
      {
	C(k,index[i])=temp(i,k);
	C(index[i],k)=temp(i,k);
      }
  for (int i=0; i<permute.Rows; i++)
    for (int k=0; k<permute.Rows; k++)
      C(index[i],index[k])=temp2(i,k);
}

void EKF::checkIndex(long ind)
{
  if (ind<FeatureRows) return;
  for (int i=ind; i<ind+3; i++)
      for (LongList *llist=&Keys[i]; llist->Next; )
	{
	  MapFeature *mf=getMapFeature(Bank,Keys[i].Element);
	  long *index=mf->getIndex();
	  if (mf)
	    {
	      int tst=1;
	      for (int j=0; j<mf->Bdual.Columns; j++)
		{
		  int k=index[j];
		  if (k==i) tst=0;
		  else
		    {
		      int tst2=1;
		      for (LongList *llist2=&Keys[k]; llist2->Next;
			   llist2=llist2->Next)
			if (llist2->Element==mf->Key)tst2=0;
		      if (tst2)Keys[k].add(mf->Key);
		    }
		}
	      if (tst)
		{
		  llist->remove(0);
		}
	      else llist=llist->Next;
	    }
	  else llist->remove(0);
      }
  for (int i=ind+2; i>ind-1; i--)
    if(Keys[i].Next==0)
      removeRows(i,1);
}
void EKF::insertRows(int currentDim,int addedDim, long *index, long featureKey)
{
  if (addedDim)
    {
      while (C.Rows+addedDim>C.RowInc)reallocate(C.RowInc+100);
       if (currentDim!=0)
	{
	  LongList tempKeys[currentDim];
	  Matrix temp(currentDim,C.Rows);
	  double tempdx[currentDim];
	  temp.Rows=1;
	  C.Rows=1;
	  long ind[currentDim],ind2[currentDim];
	  for (long i=0; i<currentDim; i++)
	    {
	      ind[i]=index[i];
	      ind2[i]=ind[i];
	      int rowoff=C.RowInc*ind[i];
	      C.Element+=rowoff;
	      temp=C;
	      C.Element-=rowoff;
	      temp.Element+=temp.RowInc;
	      tempdx[i]=Dx(ind[i],0);
	      tempKeys[i]=Keys[ind[i]];
	      for (LongList *llist=tempKeys+i; llist->Next; llist=llist->Next)
		{
		  MapFeature *mf=getMapFeature(Bank,llist->Element);
		  if (mf)mf->changeIndex(ind[i],-ind[i]);
		}
	    } 
	  C.Rows=C.Columns;
	  temp.Element-=(temp.RowInc*currentDim);
	  temp.Rows=currentDim;
	  Matrix tempC(currentDim);
	  tempC.Columns=1;
	  temp.Columns=1;
	  for (int i=0;i<currentDim;i++){
	    temp.Element+=ind[i];
	    tempC=temp;
	    temp.Element-=ind[i];
	    tempC.Element++;
	  }
	  tempC.Element-=currentDim;
	  tempC.Columns=currentDim;
	  temp.Columns=temp.RowInc;
	  for (int i=0;i<currentDim;i++){
	    temp.deleteColumn(ind2[i]);
	    removeRows(ind2[i],1);
	    for (int j=i+1; j<currentDim; j++)
	      if (ind2[j]>ind2[i])ind2[j]--;
	  }	  
	  int top=C.Rows; 
	  C.Rows+=currentDim;
	  C.Columns+=currentDim;      
	  Dx.Rows+=currentDim;
	  int r=C.Rows;
	  int c=C.Columns;

	  C.offset(top,0,currentDim,top);
	  temp.offset(0,currentDim,top);
	  C=temp;
	  C.offset(0,top,currentDim,currentDim);
	  tempC.offset(0,currentDim,currentDim);
	  C=tempC;
	  C.reset(r,c);
    
	  for (int i=0; i<top;i++)
	    for (int j=top; j<C.Rows; j++)
	      C(i,j)=C(j,i);
	  for (int i=0;i<currentDim;i++){
	    Dx(top+i,0)=tempdx[i];;
	    Keys[top+i]=tempKeys[i];
	    for (LongList *llist=Keys+i+top; llist->Next; llist=llist->Next)
	      {
		MapFeature *mf=getMapFeature(Bank,llist->Element);
		if (mf)mf->changeIndex(-ind[i],top+i);
	      }
	  }
	}
      C.Rows+=addedDim;
      C.Columns=C.Rows;
      Dx.Rows=C.Rows;
      for (int i=C.Rows-addedDim; i<C.Rows; i++)
	{
	  Dx(i,0)=0;
	  for (int j=0; j<C.Rows; j++)
	    C(i,j)=C(j,i)=0;
	  C(i,i)=10000;
	  Keys[i].add(featureKey);
	}
    }
}
void EKF::setRobotCov()
{ 
  int r=C.Rows;
  C.offset(0,RobotRows);
  Robot.Covariance=C;
  C.reset(r,r);
}

int EKF::getFeatureCov(Matrix &fCov,MapFeature *mf)
{
  long *ind=mf->getIndex();
  long dim=0;
  if (ind)
    while ((ind[dim]>0)&&(dim<mf->FullDim))dim++;      
  fCov.reallocate(dim);
  for(int i=0; i<dim; i++)
    for(int j=0; j<dim; j++)
      fCov(i,j)=C(ind[i],ind[j]);
  return dim;
}

int EKF::getFeatureRobotCorr(Matrix &corr,MapFeature *mf)
{
  long *ind=mf->getIndex();
  long dimF=0;

  if (ind) {
    // Incremenet/count number of feature dimensions as long as there
    // is something nonzero in the index vector that stores the
    // indices into the covariance matrix.
    while ( (ind[dimF] > 0) && (dimF < mf->FullDim)) {
      dimF++;      
    }
  }

  int dimR=Robot.Covariance.Rows;

  corr.reallocate(dimR, dimF);
  for(int i=0; i<dimR; i++) {
    for(int j=0; j<dimF; j++) {
      corr(i,j)=C(i,ind[j]);
    }
  }

  return dimF;
}
