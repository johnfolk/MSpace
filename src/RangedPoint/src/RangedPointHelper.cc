// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2006 John Folkesson
//    

#include "RangedPointHelper.hh"
#include "SensorData.hh"
#include "CureDebug.hh"

#ifndef DEPEND
#include <sstream>
#endif

using namespace Cure;

MapFeature * RangedPointHelper::makeMapFeature(Cure::Pose3D &sensorpose,
					       Cure::Measurement &m,
					       bool addToVisable)
{
  double v[3],rho;
  rho=m.V(2,0)*sin(m.V(0,0));
  v[1]=m.V(2,0)*cos(m.V(0,0));
  v[0]=rho*cos(m.V(1,0));
  v[2]=rho*sin(m.V(1,0));
  MapRangedPoint *mf=makeMapRangedPoint(sensorpose,v);
  if (m.Key>0)
    {
      if (m.MeasurementType&1)
	mf->forceExtend();
      addTrackedKey(m.Key,mf->Key);
    }
  if (addToVisable){
    PosedFeature *pf=makePosedRangedPoint(mf);
    if (pf){
      pf->transform(sensorpose, sensorpose.getCovType());
      PosedFeatures.add(pf);
      Visable.add(pf);
    }
  }
  return mf;
}

bool 
RangedPointHelper::supportsSubconfig(int sc)
{
  return (1 <= sc && sc <= 4);
}
void 
RangedPointHelper::printConfiguration()
{
  std::cerr << "Configured RangedPointHelper with "
	    << " DistanceThreshold=" << TempletFeature.DistanceThreshold
	    << " RoughSearchRange=" << RoughSearchRange
	    << " NearSearchRange=" << NearSearchRange
	    << " MatchThreshold=" << MatchThreshold
	    << std::endl;
}
int
RangedPointHelper::config(const std::string &arglist)
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
RangedPointHelper::configVer1(const std::string &arglist) 
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
    {
      if (str >> RoughSearchRange
	  >> NearSearchRange
	  >> NewFeatureRange
	  >> MatchThreshold
	  >> TempletFeature.DistanceThreshold) {
	return 0;
      }
    }
    break;
  case 1:
    if (str >> RoughSearchRange
        >> NearSearchRange
        >> NewFeatureRange
        >> MatchThreshold) {
      return 0;
    }
    break;
  case 2:
    if (str >> TempletFeature.DistanceThreshold) {
      return 0;
    }
    break;
  case 4:
    return 0;
    break;
  default:
    CureCERR(20) << "Cannot handle subcfg " << subcfg << std::endl;
    return 1;
  }   // Could not read all parameters but the ones it did read will
  // change and the rest will have their default values
  return 1;
}

int RangedPointHelper::addDistance(Match & mat,double distance)
{
 mat.PathDistance=distance;
 mat.Thresholds.reallocate(1,3);
 mat.Thresholds(0,0)=.5;
 mat.Thresholds(0,1)=.5;
 mat.Thresholds(0,2)=RoughSearchRange;
 mat.Metric.reallocate(3);
 mat.Metric=25;
 mat.Metric(2,2)=1;
 return 0;
}

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
unsigned short  RangedPointHelper::findMergeFeatures(PosedFeature * pfp[2],
				PosedFeatureList * searchList)
{
  PosedRangedPoint *pl[2];
  if (searchList==0)searchList=&Visable;
  double thresholds[3];
  thresholds[0]=.1;
  thresholds[1]=.1;
  thresholds[2]=.3;
  double rect[4];
  for (PosedFeatureList *mlist=searchList; mlist->Next; 
       mlist=mlist->Next)
    {
      PosedFeature *pf=mlist->Element;
      if (pf)
	{
	  pl[0]=castPosedRangedPoint(pf);
	  if (pl[0])
	    {  
	      rect[0]=pl[0]->Center.X[0];
	      rect[1]=pl[0]->Center.X[1];
	      rect[2]=pl[0]->Center.X[0];
	      rect[3]=pl[0]->Center.X[1];
	      rect[0]-=thresholds[2];
	      rect[1]-=thresholds[2];
	      rect[2]+=thresholds[2];
	      rect[3]+=thresholds[2];
	      for (PosedFeatureList *mlist2=mlist->Next; mlist2->Next; 
		   mlist2=mlist2->Next)
		{
		  pf=mlist2->Element;
		 if (pf)
		   if (pf->inRectangle(rect))
		     {
		       if (!(pf->roughMatch(pl[0]->Z.Element,thresholds)))
			 {
			   unsigned short res=0;
			   pl[1]=castPosedRangedPoint(pf);
			   MapRangedPoint *ml0=
			     getMapRangedPoint(pl[0]->FeatureKey);
			   MapRangedPoint *ml1=
			     getMapRangedPoint(pl[1]->FeatureKey);
			   unsigned short typ2=ml1->Bdual.Columns;
			   unsigned short typ1=ml0->Bdual.Columns;
			   if ((typ2<3)&&(typ1<3))res=1;
			   else
			     {
			       double d=pl[0]->Center(1)-
				 pl[1]->Center(1);
			       if ((typ2<3)||(typ1<3)){
				 if ((d<.04)&(d>-.04))res=1;
			       }
			       else if ((d<.02)&(d>-.02))res=1;
			     }
			   if (res)
			     {
			       if (typ2==typ1)
				 {
				   if (pl[0]->FeatureKey<pl[1]->FeatureKey)
				     res=2;
				   else  res=1;
				 }
			       else if (typ2==0)res=2;
			       }
			   if (res>1)
			     {
			       pfp[0]=pl[1];
			       pfp[1]=pl[0];
			       res=1;
			       typ2=typ1;
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
	}
    } 
  return 0;
}


int  RangedPointHelper::getMergeContstraint(Matrix &a, Matrix &e,Matrix & cov,
					     double distance,
					     unsigned short typ, 
					    PosedFeature * pl[2])
{
  MapRangedPoint *ml0=getMapRangedPoint(pl[0]->FeatureKey);
  MapRangedPoint *ml1=getMapRangedPoint(pl[1]->FeatureKey); 
  if (!ml0)return 0;
  if (!ml1)return 0;
  int dim0=ml0->Bdual.Columns;
  int dim1=ml1->Bdual.Columns;
  a.reallocate(dim0,dim0+dim1);
  e.reallocate(dim0,1);
  cov.reallocate(dim0);
  cov=0;
  if (dim0==0)return 0;
  a.Columns=dim0;
  a=1;
  a.Element+=dim0;
  a.Columns=dim1;
  a=-1;
  a.Element-=dim0;
  a.Columns=dim0+dim1;
  Matrix x(a.Columns,1);
  x.Rows=dim0;
  ml0->getP(x);
  x.Element+=dim0;
  ml1->getP(x);
  x.Element-=dim0;
  x.Rows=a.Columns;
  e.multiply_(a,x);
  x.multTranspose(e,e,1);
  cov=(x(0,0)+1E-6);
  return 0;
}
