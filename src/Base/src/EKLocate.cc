
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "EKLocate.hh"


using namespace Cure;
EKLocate::EKLocate(MapBank *b, int dim){
  Bank=b;
  Check=0;
  ParameterRows=0;
  Sensors=0;
  for (int i=0; i<5; i++)
    {
      SensorRows[i]=0;
      SensorCovTypes[5]=0;
      OffsetCovTypes[0]=0;
    }
  RobotCovType=0;
  RobotRows=0;
  Gain=1;
  reallocate(dim);
}
EKLocate::~EKLocate()
{
}
void EKLocate::reallocate(int dim)
{
  if (dim<1)dim=1;
  if (C.RowInc==0)
    {
      if (dim<C.RowInc)dim=C.RowInc;
      int r=C.Rows;
      C.grow(dim);
      Dx.grow(dim,1);
      Dx.reset(r,1);
      C.reset(r,r);
      return;
    }
  C.reallocate(dim);
  C.reset(0,0);
  Dx.reallocate(dim,1);
  Dx.reset(0,1);
}
void EKLocate::initRobotPose(Pose3D & pose)
{ 
  if (0!=C.Rows)
    {
      std::cerr<<"EKLocate: You have to initRobotPose before anything\n";
      return ;
    }
  if (C.RowInc<pose.Covariance.Rows)reallocate(pose.Covariance.Rows);
  Robot=pose;
  RobotCovType=Robot.getCovType();
  RobotRows=calcRows(RobotCovType);
  C.Rows=RobotRows;
  ParameterRows=C.Columns=RobotRows;
  // FeatureRows=RobotRows;
  C=1E-9;
  Dx.Rows=RobotRows;
  Dx=0;
  pose.Covariance=C;
  Robot.Covariance=C;
}

int EKLocate::resetRobotPose(const Pose3D & pose)
{ 
  if (0==C.Rows) {
    std::cerr<<"EKLocate: You have to call initRobotPose before anything\n";
    return 1;
  }
  /*
//This test is no longer valid
  if (C.Rows != pose.Covariance.Rows) {
    std::cerr << "EKLocate: You need to have same size for pose cov:"
	      << C.Rows << "!=" << pose.Covariance.Rows << std::endl; 
    return 1;
  }
  */
  if (RobotCovType != pose.getCovType()) {
    std::cerr << "EKLocate: You need to have the same cov type, "
              << "got " << pose.getCovType() 
              << " wanted " << RobotCovType << "\n";
    return 1;
  }
  Robot=pose;
  //We now Need to zero out the Sensor cros covariance too.
  int r=C.Rows;
  C.offset(0,0,pose.Covariance.Rows,r);
  C=0;
  C.offset(0,0,r,pose.Covariance.Rows);
  C=0;
  C.offset(0,0,pose.Covariance.Rows,pose.Covariance.Rows);
  C=pose.Covariance;
  C.reset(r,r);
  return 0;
}

int EKLocate::initSensorOffset(Cure::Pose3D & offsetPose, int sensorCovType)
{
   if (Sensors>4)
    {
      std::cerr<<"EKLocate can only handle up to 5 Sensors\n";
      return -1;
    }

  if (0==C.Rows)
    {
      std::cerr<<"EKLocate: You have to initRobotPose before initSensorOffset\n";
      return -1;
    }
  if (ParameterRows!=C.Rows)
    {
      std::cerr<<"EKLocate: You have to initSensorOffset before adding Parameters\n";
      return -1;
    }
  if (!((sensorCovType&11)==11)){
    std::cerr<<"We don't allow sensor Covariance types with less than xy and theta\n";
    sensorCovType=(sensorCovType|11);
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
  //  FeatureRows=r;
  Offsets[Sensors]=offsetPose;
  SensorCovTypes[Sensors]=sensorCovType;
  OffsetCovTypes[Sensors]=offsetPose.getCovType();
  Sensors++;
  return (Sensors-1);
}
void EKLocate::offsetPredict(Cure::Pose3D & offsetpose, int sensorIndex)
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
void EKLocate::incrementalPredict(Pose3D & pose, Pose3D & inc)
{
  if (RobotCovType!=pose.getCovType())
    {
      std::cerr<<"EKLocate: IncrementalPredict must use pose covariance type: "
	       <<RobotCovType<<"\n";
      return ;
    }
  Matrix jac;
  inc.Covariance*=Gain;
  if (inc.Time>Robot.Time)
    Robot.Time=inc.Time;
  Robot.add(Robot,inc,&jac);

  // Get only the jac wrt robot state variables
  jac.Columns=RobotRows;

  // Remember how many rows the covariance matrix has (and thus columns)
  int r=C.Rows;

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
  pose=Robot;
}

int EKLocate::checkEigen()
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
int EKLocate::update(Match &mat, double mahalanobisTest,
		Cure::Pose3D & pose)
{
  PosedFeature &pf=*mat.MatchedFeature;

  //int *index=pf.getIndex();
  //if (!index)return -1;

  // Since this function might be called when looping through a set of
  // measurements and thus this might not be the first call, we need
  // to propagate the change in P-space coordinates to the "real
  // coordinates" of this particular feature and possibly change the
  // projection matrix before we can perform the update with the new
  // measurement.
  //  updateP(pf);

  //Robot=pose;

  Measurement *m=mat.Measure;

  // Count the number of P-dimensions for this feature
  //int dim=0;
  //while ((index[dim])&&(dim<pf.FullDim))dim++;

  // Perform the update if there are P-dimensions
  if (1)//(dim)
    {
      Matrix cov;
      // Calculate the innovation and all Jacobians
      if (pf.linearizeMeasurement( m->V, Robot,m->MeasurementType))return -1;
      if (pf.Eta.Rows==0) return 0;

      // Calculate measurement covariance R=Jev*Cvv*Jev^T
      Matrix invcov;
      invcov.multTranspose_(m->CovV,pf.Jev,2);	        
      cov.multiply_(pf.Jev,invcov);

      Matrix **(jacobians)[1];
      Matrix * jac[1];
      jacobians[0]=jac;
      jac[0]=&(pf.Jes);
      //      jac[1]=&(pf.Jep);

      // We need to keep track of the indices of all variables that
      // our measurement depends on since these will give non-zero
      // columns if a full Jacobian were to be calculated. The
      // linearizeMeasurement call above have given us the non-zero
      // parts but we need to define where they fit in when we do the
      // matrix multiplications.
      int **(indecies)[1];
      int *(ind)[1];
      indecies[0]=ind;

      // Create an array of indices for the robot pose variables
      int poseindex[RobotRows];
      indecies[0][0]=poseindex;
      for (int i=0; i<RobotRows; i++)
	poseindex[i]=i;

      // Set the indices for the map feature that we have measured
      //  indecies[0][1]=index;

      Matrix Sinv,A,w; 
      // Calculate A=J*C and Sinv=(innovation covariance)^(-1)
      getMetric(A, Sinv,cov,jacobians,indecies, 1,1);

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
      updatePose();

      // Update the feature coordinates based on the change in P-space
      // coordinates. Note that we 
      //updateP(pf); 
    }
  setRobotCov();
  pose=Robot;
  return 0;

}
void EKLocate::updatePose()
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
  C.reset(r,r);
}
/*
int EKLocate::removeRows(int start,int dim)
{
  for (int i=0; i<dim;i++)
    {
      C.deleteRow(start);
      C.deleteColumn(start);
      Dx.deleteRow(start);
    }
  return C.Rows;
}
*/
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
void EKLocate::getMetric(Matrix & A, Matrix & Sinv, 
		    Matrix & cov, Matrix *** jacobians,
		    int *** index,  int numblocks,int numJac)
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
	  int *indx=index[i][m];
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
	  int *indx=index[i][m];
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
double EKLocate::getloglikelihood(Measurement *m, PosedFeature *pf,
				  Transformation3D &sensorpose)
{
  if (pf->transform(sensorpose,0))return -1E32;//feature should not be visable
  // Calculate the innovation and the Jacobians Jeo and Jev
  if (pf->calcEta(m->V,m->MeasurementType))return -1E32;
  int n=pf->Eta.Rows;
  if (n==0) return -1E32;
  Matrix a,b,invcov;
  a.multTranspose_(m->CovV,pf->Jev,2);	        
  invcov.multiply_(pf->Jev,a);
  if (invcov.invert())return -1E32;    
  a.multiply_(invcov,pf->Eta);
  b.multTranspose_(pf->Eta,a,1);
  return -b(0,0)/2;
}
double  EKLocate::getloglikelihood(Match *matches, int n,
				   Transformation3D &probot,
				   int sensorIndex,
				   double nomatchlikelihood,
				   double matchtest)
{
  if (sensorIndex<0)return 0;
  Transformation3D sensorpose(probot);
  sensorpose+=Offsets[sensorIndex];
  double r=0;
  for (int i=0; i<n; i++)
    {
      double d=nomatchlikelihood; 
      if (matches[i].MatchDistance>-1)
	{
	  d=getloglikelihood(matches[i],sensorpose);
	}  
      if (d<matchtest)d=nomatchlikelihood;
      r+=d;
    }
  return r;
}
double EKLocate::getloglikelihood(Match &mat,
				  Cure::Transformation3D &sensorpose)
{
  PosedFeature *pf=mat.MatchedFeature;
  Measurement *m=mat.Measure;
  return getloglikelihood(m, pf,sensorpose);
}

int EKLocate::updateSensor(Match &mat, double mahalanobisTest, int sensorIndex)
{
  PosedFeature &pf=*mat.MatchedFeature;
  //  int *index=pf.getIndex();
  //if (!index)return -1;
  // updateP(pf);
  //  int dim=0;
  Measurement *m=mat.Measure;
  //while ((index[dim])&&(dim<pf.FullDim))dim++;
  //if (dim)
    {
      Matrix Sinv,A,w, temp0,temp1;
      Matrix **(jacobians)[1];
      Matrix * jac[2];
      jacobians[0]=jac;
      int **(indecies)[1];
      int *(ind)[2];
      indecies[0]=ind;   
      int cols=(OffsetCovTypes[sensorIndex]<<6);
      cols+=RobotCovType;
      Pose3D sensorpose;
      Matrix jsro(6,12);
      sensorpose.setCovType(SensorCovTypes[sensorIndex]);
      sensorpose.doAplusB(Robot,Offsets[sensorIndex], jsro,
			  SensorCovTypes[sensorIndex], cols, 1);
      Matrix cov;
      if ( pf.linearizeMeasurement( m->V, sensorpose,m->MeasurementType))
	return -1;
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
      //      jac[2]=&pf.Jep;
      int poseindex[RobotRows];
      indecies[0][0]=poseindex;
      for (int i=0; i<RobotRows; i++)
	poseindex[i]=i;
      int offindex[jeoff.Columns];
      indecies[0][1]=offindex;
      for (int i=0; i<jeoff.Columns; i++)
	offindex[i]=SensorRows[sensorIndex]+i;
      //      indecies[0][2]=index;
      getMetric(A, Sinv,cov,jacobians,indecies,1,2);
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
      //updateP(pf); 
    }
  return 0;
}
void EKLocate::updatePoses()
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
}

void  EKLocate::updateSensor(Match *matches, int n, Pose3D &probot,
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
	      //addinfo(matches[i],psensor,SensorCovTypes[sensorIndex],map2info);
	      //extend(matches[i],sensorIndex);
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
	      // addinfo(matches[i],psensor,SensorCovTypes[sensorIndex],map2info);
	      //extend(matches[i],sensorIndex);
	    }
	}  
      matches[i].clear();//this will delete the PosedFeature
    }
  probot=Robot;
}

void  EKLocate::update(Match *matches, int n, Pose3D &psensor,
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
	      // addinfo(matches[i],map2info);
	      // extend(matches[i]);
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
	      //addinfo(matches[i],map2info);
	      // extend(matches[i]);
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
void  EKLocate::update(Match *matches, int n, Pose3D &psensor, 
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
	      //addinfo(matches[i],map2info);
	      //extend(matches[i]);
	    }
	  matches[i].clear();  //this will delete the PosedFeature
	}  
  for (int i=0; i<n; i++)
    {
      if (matches[i].MatchDistance>-1)
	{
	  if(update(matches[i],manhanalobisTest,Robot)==0)
	    {	
	      //addinfo(matches[i],map2info);
	      // extend(matches[i]);
	    }
	}  
      matches[i].clear();//this will delete the PosedFeature
    }
  psensor=Robot;
}
void EKLocate::setRobotCov()
{ 
  int r=C.Rows;
  C.offset(0,RobotRows);
  Robot.Covariance=C;
  C.reset(r,r);
}
