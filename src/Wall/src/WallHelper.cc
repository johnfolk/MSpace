// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "WallHelper.hh"
#include "CureDebug.hh"

#ifndef DEPEND
#include <sstream>
#endif

using namespace Cure;

MapFeature * WallHelper::makeMapFeature(Cure::Pose3D &sensorpose,
					Cure::Measurement &m,
					bool addToVisable)
{
  MapFeature *mf=makeMapWall(sensorpose,m.V.Element);
  if (addToVisable){
    PosedFeature *pf=makePosedFeature(mf);
    if (pf){
      pf->transform(sensorpose, sensorpose.getCovType());
      PosedFeatures.add(pf);
      Visable.add(pf);
    }
  }
  return mf;
}
bool
WallHelper::supportsSubconfig(int sc)
{
  return (1 <= sc && sc <= 3);
}
void 
WallHelper::printConfiguration()
{
  std::cerr << "\nConfigured WallHelper with MatchSettings:\n" 
	    << "  RoughSearchRange=" << RoughSearchRange
	    << " NearSearchRange=" << NearSearchRange
	    << " NewFeatureRange=" << NewFeatureRange
	    << " MatchThreshold=" << MatchThreshold
	    << " EndPointUpdateThreshold=" << EndPointUpdateThreshold
	    << "\nInialization:\n"
	    << "  CountThreshold=" << TempletFeature.CountThreshold
	    << " TightnessValue=" << TempletFeature.TightnessValue
	    << " EndThreshold=" << TempletFeature.EndThreshold
	    << " LengthThreshold=" << TempletFeature.LengthThreshold
	    << " VarRhoThreshold=" << TempletFeature.VarRhoThreshold
	    << " DistanceThreshold=" << TempletFeature.DistanceThreshold
	    << " UseEndpoints=" << UseEndpoints
	    << "\nMerge Settings:\n"
	    << "  MergeOverLap=" << MergeOverLap
	    << " MergeMaxGammaError=" << MergeMaxGammaError
	    << " MergeMaxRhoError=" << MergeMaxRhoError
	    << " CornerSeparation=" << CornerSeparation
	    << std::endl;
}
int
WallHelper::config(const std::string &arglist)
{
  std::istringstream str(arglist);
  int version = -1;
  if ( !(str >> version)) {
    CureCERR(20) << "Failed to read version number for config params list\n";
    return 1;
  }

  int ret = 0;
  switch (version) {
  case 1:
    ret = configVer1(arglist);
    break;
  default:
    CureCERR(20) << "Cannot handle config version " << version << std::endl;
    return 1;
  }
  
  return ret;
}
  
int 
WallHelper::configVer1(const std::string &arglist) 
  {
  std::istringstream str(arglist);
  int version = 0;
  int subcfg = -1;
  if ( !(str >> version>> subcfg)) {
    CureCERR(20) << "Failed to read subcfg number for config params list\n";
    return 1;
  }
  switch (subcfg) {
  case 0:
    if (str >> RoughSearchRange
	>> NearSearchRange
	>> NewFeatureRange
	>> MatchThreshold
	>> EndPointUpdateThreshold 
	>> TempletFeature.CountThreshold
	>> TempletFeature.TightnessValue
	>> TempletFeature.EndThreshold
	>> TempletFeature.LengthThreshold
	>> TempletFeature.VarRhoThreshold
	>> TempletFeature.DistanceThreshold
	>> UseEndpoints
	>> MergeOverLap
	>> MergeMaxGammaError
	>> MergeMaxRhoError
	>> CornerSeparation) {
      return 0;
    }
    break;
  case 1:
    if (str >> RoughSearchRange
	>> NearSearchRange
	>> NewFeatureRange
	>> MatchThreshold
	>> EndPointUpdateThreshold) {
      return 0;
    }
    break;
  case 2:
    if (str >> TempletFeature.CountThreshold
	>> TempletFeature.TightnessValue
	>> TempletFeature.EndThreshold
	>> TempletFeature.LengthThreshold
	>> TempletFeature.VarRhoThreshold
	>> TempletFeature.DistanceThreshold
	>> UseEndpoints) {
      return 0;
    }
    break;
  case 3:
    if (str >> MergeOverLap
	>> MergeMaxGammaError
	>> MergeMaxRhoError
	>> CornerSeparation) {
      return 0;
    }
    break;
  default:
    CureCERR(20) << "Cannot handle subcfg " << subcfg << std::endl;
    return 1;
  }

  // Could not read all parameters but the ones it did read will
  // change and the rest will have thier default values
  return 1;
}

void WallHelper::loosenMerging(double factor)
{
  if ( MergeOverLap>0) MergeOverLap/=factor;
  else MergeOverLap*=factor;
  CornerSeparation*=factor;
  MergeMaxGammaError*=factor;
  MergeMaxRhoError*=factor; 
}
int WallHelper::addDistance(Match & mat,double distance)
{
  Measurement &m=*mat.Measure;

  mat.PathDistance=distance;
  if (m.SensorType!=1)return 1;
  if (!(m.MeasurementType<8))return 1;
  if (!UseEndpoints)
    if(m.MeasurementType>1)
      m.MeasurementType=1;
  
  if (m.Z(1,0)<0){
    m.Z(1,0)=-m.Z(1,0);
    if (m.Z(0,0)<0)m.Z(0,0)+=M_PI;
    else m.Z(0,0)-=M_PI;
  }
  Line.StartPoint.setXY(m.V.Element);
  Line.EndPoint.setXY(m.V.Element+2);
  mat.Thresholds.reallocate(1,2);
  mat.Thresholds(0,1)=RoughSearchRange;
  double d=Line.StartPoint.getDistance(Line.EndPoint);
  mat.Thresholds(0,0)=10*RoughSearchRange/d;
  d*=d;
  if (mat.Thresholds(0,0)>.5)mat.Thresholds(0,0)=.5;
  mat.Metric.reallocate(2);
  mat.Metric=2/(mat.Measure->CovV(0,0)+mat.Measure->CovV(1,1));
  mat.Metric(0,0)*=d;
  return 0;
}
/*
void WallHelper::makeMatch(Match & mat,Matrix & v,double distance)
{
  makeMeasurement(*mat.Measure,v, distance);
  addDistance(mat, distance);
}
void WallHelper::makeMeasurement(Measurement & m,Matrix & v,double distance)
{
  m.W.reallocate(v.Rows,3);
  m.W.Element++;
  m.W.Columns=2;
  m.W=v;
  m.W.Element--;
  m.W.Columns=3;
  for (int i=0; i<v.Rows;i++)m.W(i,0)=distance;
  if (v.Rows<3)
    {
      return;
    }
  Point2DCloud cloud;
  cloud=m.W;
  cloud.fitLine();
  m.V.reallocate(4,1);
  m.Z.reallocate(2,1);
  m.Z(0,0)=cloud.Gamma;
  m.Z(1,0)=cloud.Rho;
  m.CovV.reallocate(4);
  m.V(0,0)=v(0,0);
  m.V(1,0)=v(0,1);
  m.V(2,0)=v(v.Rows-1,0);
  m.V(3,0)=v(v.Rows-1,1);
  m.CovV=cloud.SigmaRho;
  Line.StartPoint.setXY(m.V.Element);
  Line.EndPoint.setXY(m.V.Element+2);
  m.BoundingBox.reallocate(2);
  double *rect=m.BoundingBox.Element;
  Line.StartPoint.getXY(rect);
  Line.EndPoint.getXY(rect+2);
  if (rect[0]>rect[2])
    {
      double d=rect[0];
      rect[0]=rect[2];
      rect[2]=d;
    }
  if (rect[1]>rect[3])
    {
      double d=rect[1];
      rect[1]=rect[3];
      rect[3]=d;
    }
}


*/
//merge starts by match to a expected measurement
//then extra critera are tested 
//So for a wall that might be overlap (RoughSearchRange<0),
//for a HLine it might be similar hieght or hieght not known.
//Then the SLAM will need to deal with enforcing the 
//Constraints.  We need CP=b for the P's
/*
 * try finding a merge candidates from search list
 * If found remove them from searchlist and return them in 
 * pl1 and pl2.
 * @return 1 after finding one 0 otherwise
 */
unsigned short  WallHelper::findMergeFeatures(PosedFeature * pfp[2],
				PosedFeatureList * searchList)
{
  PosedWall *pl[2];
  if (searchList==0)searchList=&Visable;
  double thresholds[2];
  thresholds[0]=MergeMaxGammaError*10; 
  thresholds[1]= MergeMaxRhoError*10;
  double rect[4],r2[4];
  for (PosedFeatureList *mlist=searchList; mlist->Next; 
       mlist=mlist->Next)
    {
      PosedFeature *pf=mlist->Element;
      if (pf)
	{
	  pl[0]=castPosedWall(pf);
	  MapWall *ml0=getMapWall(pl[0]->FeatureKey);
	  if (pl[0])
	    {      
	      r2[0]=rect[0]=pl[0]->Line.StartPoint.X[0];
	      r2[1]=rect[1]=pl[0]->Line.StartPoint.X[1];
	      r2[2]=rect[2]=pl[0]->Line.EndPoint.X[0];
	      r2[3]=rect[3]=pl[0]->Line.EndPoint.X[1];
	      
	      if (rect[0]>rect[2])
		{
		  rect[0]=r2[2];
		  rect[2]=r2[0];
		  r2[0]=rect[0];
		  r2[2]=rect[2];
		}
	     if (rect[1]>rect[3])
		{
		  rect[1]=r2[3];
		  rect[3]=r2[1];
		  r2[1]=rect[1];
		  r2[3]=rect[3];
		}
	     double ov=MergeOverLap;
	     r2[0]+=ov;
	     r2[1]+=ov;
	     r2[2]-=ov;
	     r2[3]-=ov;
	     if (r2[0]>r2[2])
	       {
		 double d=r2[0];
		 r2[0]=r2[2];
		 r2[2]=d;
	       }
	     if (r2[1]>r2[3])
	       {
		 double d=r2[1];
		 r2[1]=r2[3];
		 r2[3]=d;
	       }
	     if (ov<0)ov=-ov;
	     ov+=(thresholds[1]+.01);
	     ov+=.5;
	     rect[0]-=ov;
	     rect[1]-=ov;
	     rect[2]+=ov;
	     rect[3]+=ov;
	    
	       // for (PosedFeatureList *mlist2=mlist->Next; mlist2->Next; 
	     for (PosedFeatureList *mlist2=&PosedFeatures; mlist2->Next; 
		  mlist2=mlist2->Next)
	       {
		 pf=mlist2->Element;
		 if (pf)
		   if (pf->inRectangle(rect))
		     if(pf->inRectangle(r2))
		       {
			 pl[1]=castPosedWall(pf);
			 MapWall *ml1=getMapWall(pl[1]->FeatureKey); 
			 if ((ml1)&&(ml1!=ml0))
			   {
			     ml0->calcTangent();
			     ml1->calcTangent();
			     double d=ml0->Tangent[0]*ml1->Tangent[0]+
			       ml0->Tangent[1]*ml1->Tangent[1];
			     if (d>0)
			       if (!(pf->
				   roughMatch(pl[0]->Z.Element,thresholds)))
			       {
				 unsigned short res=1;
				 unsigned short typ2=ml1->Bdual.Columns;
				 unsigned short typ1=ml0->Bdual.Columns;
				 MapWall *ml=getMapWall(pl[0]->FeatureKey);
				 int vcount0=0;
				 if (ml)vcount0=ml->Cloud.Cloud.Rows;
				 ml=getMapWall(pl[1]->FeatureKey);
				 int vcount1=0;
				 if (ml)vcount1=ml->Cloud.Cloud.Rows;
			     
				 if (typ1==0){
				   if (typ2==0)
				     {
				       if (vcount0>vcount1)res=2;
				       else  res=1;
				     }
				   else res=1;
				 }
				 else if (typ2==0)res=2;
				 else if (typ1==2){
				   if (typ2==2)
				     {
				       if (vcount0>vcount1)res=2;
				       else  res=1;
				     }
				   else res=1;
				 }
				 else if (typ2==2)res=2;
				 else if (typ1==3){
				   if (typ2==3)
				     {
				       if (vcount0>vcount1)res=2;
				       else  res=1;
				     }
				   else res=1;
				 }
				 else if (typ2==3)res=2;
				 else 
				   {
				     if (vcount0>vcount1)res=2;
				     else  res=1;
				   }
				 if (res>1)
				   {
				     pfp[0]=pl[1];
				     pfp[1]=pl[0];
				     res-=1;
				   }
				 else{
				   pfp[0]=pl[0];
				   pfp[1]=pl[1];
				 }
				 searchList->removeFeature(pl[0]);
				 return res;
			       }
			 }
		 }
	       }
	     rect[0]+=ov;
	     rect[1]+=ov;
	     rect[2]-=ov;
	     rect[3]-=ov;
	     ov=CornerSeparation;
	     if (ov<0)ov=-ov;
	     rect[0]-=ov;
	     rect[1]-=ov;
	     rect[2]+=ov;
	     rect[3]+=ov;
	      if (rect[0]>rect[2])
		{
		  double d=rect[0];
		  rect[0]=rect[2];
		  rect[2]=d;
		}
	     if (rect[1]>rect[3])
		{
		  double d=rect[1];
		  rect[1]=rect[3];
		  rect[3]=d;
		}

	     if (UseEndpoints)
	       //    for (PosedFeatureList *mlist2=mlist->Next; mlist2->Next; 
	       for (PosedFeatureList *mlist2=&PosedFeatures; mlist2->Next; 
		    mlist2=mlist2->Next){
		 pf=mlist2->Element;
		 if (pf)
		   if (pf->inRectangle(rect))
		     {
		       pl[1]=castPosedWall(pf);
		       if (pl[1])
			 {

			   MapWall *ml0=getMapWall(pl[0]->FeatureKey);
			   MapWall *ml1=getMapWall(pl[1]->FeatureKey);
			   if (((ml0)&&(ml1))&&(ml0!=ml1))
			     {
			       ml0->calcTangent();
			       ml1->calcTangent();
			       double d=ml0->Tangent[0]*ml1->Tangent[0]+
				 ml0->Tangent[1]*ml1->Tangent[1];
			       if ((d<.75)&&(d>-.75))
				 {
				   Vector3D v01=(*ml1->Points[0])-(*ml0->Points[1]);
				   Vector3D v10=(*ml1->Points[1])-(*ml0->Points[0]);
			       
				   unsigned short res2=0;
				   double d01=v01*v01;
				   double d10=v10*v10;
				   if (d01<CornerSeparation)//(d01<.04)
				     {
				       ml0->calcTangent();
				       ml1->calcTangent();
				       d=v01(0)*ml0->Tangent[0]+
					 v01(1)*ml0->Tangent[1];
				       v01(0)-=d*ml0->Tangent[0];
				       v01(1)-=d*ml0->Tangent[1];
				       d=v01(0)*ml1->Tangent[0]+
					 v01(1)*ml1->Tangent[1];
				       v01(0)-=d*ml1->Tangent[0];
				       v01(1)-=d*ml1->Tangent[1];
				       d01=v01*v01;
				       if (d01<.0001)//(d01<.0001)
					 {
					   if ((ml1->Bdual.Columns>1)&&
					       (ml0->Bdual.Columns>1))
					     if (ml1->Points[0]!=ml0->Points[1])
					       res2=0x1010;
					   if(res2)				  
					     if ((ml0->Points[1]->numberOfFeatures()>1)||(ml1->Points[0]->numberOfFeatures()>1))
					       {
						 res2=0;//0x1100;
						 //pfp[0]=pl[1];
						 //pfp[1]=pl[0];
						 return res2;
					       }
					 }
				     }
				   else if (d10<CornerSeparation)//(d10<.01)
				     {
				       ml0->calcTangent();
				       ml1->calcTangent();
				       double d=v10(0)*ml0->Tangent[0]+
					 v10(1)*ml0->Tangent[1];
				       v10(0)-=d*ml0->Tangent[0];
				       v10(1)-=d*ml0->Tangent[1];
				       d=v10(0)*ml1->Tangent[0]+
					 v10(1)*ml1->Tangent[1];
				       v10(0)-=d*ml1->Tangent[0];
				       v10(1)-=d*ml1->Tangent[1];
				       d10=v10*v10;
				       if (d10<.0001)//(d10<.0001)	 
					 if ((ml1->Bdual.Columns>1)&&
					     (ml0->Bdual.Columns>1))
					   if (ml1->Points[1]!=ml0->Points[0])
					     res2=0x1100;
				       if(res2)				  
					 if ((ml0->Points[0]->
					      numberOfFeatures()>1)||
					     (ml1->Points[1]->numberOfFeatures()>1))
					   {
					     return 0;
					     res2=0x1010;
					     pfp[0]=pl[1];
					     pfp[1]=pl[0];
					     return res2;
					   }
				     }
				   if (res2) {
				     pfp[0]=pl[0];
				     pfp[1]=pl[1];
				     return res2;
				   }
				 }
			     }
			 }
		     }
	       }
	    }
	}
    }
  return 0;
}

int WallHelper::getMergeContstraint(Matrix &a, Matrix &e,Matrix & cov,
				     double distance,
				     unsigned short typ, PosedFeature * pl[2])
{
  int ret=0;
  MapWall *ml0=getMapWall(pl[0]->FeatureKey);
  MapWall *ml1=getMapWall(pl[1]->FeatureKey); 
  if (!ml0)return ret;
  if (!ml1)return ret;


  int dim0=ml0->Bdual.Columns;
  int dim1=ml1->Bdual.Columns;
  if (typ&1)
    {
      int erows=2;
      int limcase=0;
      if (dim0>2)
	{
	  if ((ml0->StartLimit)&&(ml1->StartLimit))limcase=1;
	  if ((ml0->EndLimit)&&(ml1->EndLimit))limcase+=2;
	  if (limcase&1)erows++; 
	  if (limcase&2)erows++; 
	}
      a.reallocate(erows,dim0+dim1);
      e.reallocate(erows,1);
      cov.reallocate(erows);
      a=0;
      cov=0;
      if (dim0==0)return ret;
      if (dim1==0)return ret;
      Matrix b, c1,c2,xi,xj,pi,pj;
      ml0->getB(b);
      if (!getCommon(ml0,ml1,c1,c2))return 0;
      ml0->getX(xi);
      ml1->getX(xj);
      pi.multiply_(b,xi);
      pj.multiply_(b,xj);
      xi.multiply_(c1,pi);
      xj.multiply_(c2,pj);
      xi-=xj;
      e(0,0)=xi(0,0);
      e(1,0)=xi(1,0);
      a(0,0)=1;
      a(1,1)=1;
      a(0,dim0)=-1;
      a(1,dim0+1)=-1;
      /*
      a(0,0)=c1(0,0);
      a(1,0)=c1(1,0);
      a(0,1)=c1(0,1);
      a(1,1)=c1(1,1);
      a(0,dim0)=-c2(0,0);
      a(1,dim0)=-c2(1,0);
      a(0,dim0+1)=-c2(0,1);
      a(1,dim0+1)=-c2(1,1);
      */
      cov(0,0)=e(0,0)*e(0,0);
      cov(1,1)=e(1,0)*e(1,0);
      if (cov(0,0)<1E-8)cov(0,0)=1E-8;
      if (cov(1,1)<1E-6)cov(1,1)=1E-6;      

      if ((e(0,0)<(MergeMaxGammaError))&&(e(0,0)>-(MergeMaxGammaError)))
	if (((e(1,0)<(MergeMaxRhoError)))&&(e(1,0)>-(MergeMaxRhoError)))
	  ret=1;
      if (e.Rows==2)return ret;     

      Cure::Matrix p0;
      Cure::Matrix p1;
      ml0->getP(p0);
      ml1->getP(p1);
      int k=2;
      if (limcase&1)
	{
	  int c=ml0->StartLimit;
	  a(k,c)=1;
	  e(k,0)=p0(c,0);
	  c=ml1->StartLimit;
	  a(k,c+dim0)=-1;
	  e(k,0)-=p1(c,0);
	  cov(k,k)=e(k,0)*e(k,0);
	  if (cov(k,k)<1E-6)cov(k,k)=1E-6;
	  k++;
	}
      if (limcase&2)
	{
	  int c=ml0->EndLimit;
	  a(k,c)=1;
	  e(k,0)=p0(c,0);
	  c=ml1->EndLimit;
	  a(k,c+dim0)=-1;
	  e(k,0)-=p1(c,0);
	  cov(k,k)=e(k,0)*e(k,0);
	  if (cov(k,k)<1E-6)cov(k,k)=1E-6;
	}
      return ret;     
    }
  else   if (typ&0x1000)
    {
      if (dim0<2)return 0;
      a.reallocate(2,dim0+dim1);
      e.reallocate(2,1);
      cov.reallocate(2);
      ml0->calcTangent();
      ml1->calcTangent();
      ml0->recenter();
      ml1->recenter();
      a.Columns=dim0;
      ml0->Bdual.Rows=2;

      if (typ&0x10)
	{
	  ml0->Bdual.Element+=2*ml0->Bdual.RowInc;
	  a=ml0->Bdual;
	  ml0->Bdual.Element-=2*ml0->Bdual.RowInc;
	}
      else
	{
	  a=ml0->Bdual;
	}
      ml0->Bdual.Rows=4;
      a.Columns=dim1;
      a.Element+=dim0;;
      ml1->Bdual.Rows=2;
      Matrix b;
      if (typ&0x100)
	{
	  ml1->Bdual.Element+=2*ml1->Bdual.RowInc;
	  b=ml1->Bdual;
	  ml1->Bdual.Element-=2*ml1->Bdual.RowInc;
	}
      else
	{
	  b=ml1->Bdual;
	}
      b*=-1;
      a=b;
      ml1->Bdual.Rows=4;
      a.Columns+=dim0;
      a.Element-=dim0;;
      Matrix x(a.Columns,1);
      x.Rows=dim0;
      ml0->getP(x);
      x.Element+=dim0;
      x.Rows=dim1;
      ml1->getP(x);
      x.Element-=dim0;
      x.Rows=a.Columns;
       e.multiply_(a,x);
      cov(0,0)=e(0,0)*e(0,0);
      cov(1,1)=e(1,0)*e(1,0);
      cov(0,1)=0;
      cov(1,0)=0;
      if (cov(0,0)<1E-6)cov(0,0)=1E-6;
      if (cov(1,1)<1E-6)cov(1,1)=1E-6;
    }
  return 0;
}

void WallHelper::getPermutation(Cure::Matrix & permute,
				PosedFeature *pf[2],
				unsigned short typ)
{
  MapWall *ml0=getMapWall(pf[0]->FeatureKey);
  MapWall *ml1=getMapWall(pf[1]->FeatureKey); 
  if (!ml0)return;
  if (!ml1)return;
  int dim0=ml0->Bdual.Columns;
  int dim1=ml1->Bdual.Columns;

  permute.reallocate(dim0+dim1);
  if (typ&1){  
    permute=1;
    return;
  }
  permute=0; 
  if ((dim0<3)||(dim1<3))
    {
      std::cerr<<"WALLHELPER error";
      return;
    }
  if (typ&0x1000)
    {
      if (ml0->NumberSharedPoints==0){
	if (typ&0x10)
	  {
	    permute(0,0)=ml0->Bdual(2,0);
	    permute(0,1)=ml0->Bdual(2,1);
	    permute(1,0)=ml0->Bdual(3,0);
	    permute(1,1)=ml0->Bdual(3,1);
	    int k=ml0->EndLimit;
	    if (k<2)
	      {
		std::cerr<<"WALLHELPER 2error"<<k;
		return;
	      }
	    permute(0,k)=ml0->Bdual(2,k);
	    permute(1,k)=ml0->Bdual(3,k);
	    permute(2,0)=ml0->Bdual(0,0);
	    permute(2,1)=ml0->Bdual(0,1);
	    if (dim0==4)
	      {
		k=ml0->StartLimit;
		if (k<2)
		  {
		    std::cerr<<"WALLHELPER 3error"<<k;
		    return;
		  }
		permute(3,k)=1;
	      }
	  } 
	else
	  {
	    permute(0,0)=ml0->Bdual(0,0);
	    permute(0,1)=ml0->Bdual(0,1);
	    permute(1,0)=ml0->Bdual(1,0);
	    permute(1,1)=ml0->Bdual(1,1);
	    int k=ml0->StartLimit;
	    if (k<2)
	      {
		std::cerr<<"WALLHELPER 4error"<<k;
		return;
	      }
	
	    permute(0,k)=ml0->Bdual(0,k);
	    permute(1,k)=ml0->Bdual(1,k);
	    permute(2,0)=ml0->Bdual(2,0);
	    permute(2,1)=ml0->Bdual(2,1);
	    if (dim0==4)
	      {
		k=ml0->EndLimit;
		if (k<2)
		  {
		    std::cerr<<"WALLHELPER 5error"<<k;
		    return;
		  }
		permute(3,k)=1;
	      }

	  }
      }
      else {
	permute.Rows=dim0;
	permute.Columns=dim0;
	if (permute.Rows==ml0->Bdual.Rows)
	  permute=ml0->Bdual;
	else
	  permute=1;
	permute.Rows+=dim1;
	permute.Columns=permute.Rows;
      }
    
      permute.Element+=(permute.RowInc+1)*dim0;
      permute.Rows=dim1;
      permute.Columns=dim1;
	
      if (ml1->NumberSharedPoints==0){
	if (typ&0x100)
	  {
	    permute(0,0)=ml1->Bdual(2,0);
	    permute(0,1)=ml1->Bdual(2,1);
	    permute(1,0)=ml1->Bdual(3,0);
	    permute(1,1)=ml1->Bdual(3,1);
	    int k=ml1->EndLimit;
	    if (k<2)
	      {
		std::cerr<<"WALLHELPER 6error"<<k;
		return;
	      }
	    permute(0,k)=ml1->Bdual(2,k);
	    permute(1,k)=ml1->Bdual(3,k);
	    permute(2,0)=ml1->Bdual(0,0);
	    permute(2,1)=ml1->Bdual(0,1);
	    if (dim1==4)
	      {
		k=ml1->StartLimit;
		if (k<2)
		  {
		    std::cerr<<"WALLHELPER 7error"<<k;
		    return;
		  }
		permute(3,k)=1;
	      }
	  } 
	else
	  {
	    permute(0,0)=ml1->Bdual(0,0);
	    permute(0,1)=ml1->Bdual(0,1);
	    permute(1,0)=ml1->Bdual(1,0);
	    permute(1,1)=ml1->Bdual(1,1);
	    int k=ml1->StartLimit;
	    if (k<2)
	      {
		std::cerr<<"WALLHELPER 8error"<<k;
		return;
	      }
	    permute(0,k)=ml1->Bdual(0,k);
	    permute(1,k)=ml1->Bdual(1,k);
	    permute(2,0)=ml1->Bdual(2,0);
	    permute(2,1)=ml1->Bdual(2,1);
	    if (dim1==4)
	      {
		k=ml1->EndLimit;
		if (k<2)
		  {
		    std::cerr<<"WALLHELPER 9error"<<k;
		    return;
		  }
		permute(3,k)=1;
	      }
	  }
      }
      else {
	if (permute.Rows==ml1->Bdual.Rows)
	  permute=ml1->Bdual;
	else
	  permute=1;
      }
      permute.Rows+=dim0;
      permute.Columns=permute.Rows;
      permute.Element-=(permute.RowInc+1)*dim0;
    }
  if (permute.Rows!=dim0+dim1)
    {
      std::cerr<<"WALLHELPER 10error";
      return;
    }
}
