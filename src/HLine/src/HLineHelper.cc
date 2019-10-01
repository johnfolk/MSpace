// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "MSpace/HLine/HLineHelper.hh"
#include "Utils/CureDebug.hh"

#ifndef DEPEND
#include <sstream>
#endif

using namespace Cure;
MapFeature * 
HLineHelper::makeMapFeature(Cure::Pose3D &sensorpose,
					Cure::Measurement &m,
					bool addToVisable)
{
  MapFeature *mf=makeMapHLine(sensorpose,m.V.Element,DistanceGuess);
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
HLineHelper::supportsSubconfig(int sc)
{
  return (1 <= sc && sc <= 4);
}
void 
HLineHelper::printConfiguration()
{
  std::cerr << "Configured HLineHelper with "
	    << " TanThreshold=" << TempletFeature.TanThreshold
	    << " RhoThreshold=" << TempletFeature.RhoThreshold
               << " DistanceThreshold=" << TempletFeature.DistanceThreshold
	    << " TrackThreshold=" << TempletFeature.TrackThreshold
	    << " Image={"
	    << " bottomLeft[0]=" << ImagePlane[0]
	    << " bottomLeft[1]=" << ImagePlane[1]
	    << " topRight[0]=" << ImagePlane[2]
               << " topRight[1]=" << ImagePlane[3]
	    << " focalLength=" << PixelInfo[0]
	    << " pixwidth=" << PixelInfo[1]
	    << " pixheight=" << PixelInfo[2] << "}"
	    << " RoughSearchRange=" << RoughSearchRange
	    << " NearSearchRange=" << NearSearchRange
	    << " MatchThreshold=" << MatchThreshold
	    << " DistanceGuess=" << DistanceGuess
	    << std::endl;
  
}
int
HLineHelper::config(const std::string &arglist)
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
HLineHelper::configVer1(const std::string &arglist) 
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
      double sw, sh, f, pW, pH;
      if (str >> RoughSearchRange
	  >> NearSearchRange
        >> NewFeatureRange
	  >> MatchThreshold
	  >> DistanceGuess
	  >> TempletFeature.TrackThreshold
	  >> TempletFeature.TanThreshold
	  >> TempletFeature.RhoThreshold
	  >> TempletFeature.DistanceThreshold
	  >> sw >> sh >> f >> pW >> pH) {
	double bl[2], tr[2];
	bl[0] = -sw;
	bl[1] = -sh;
	tr[0] = sw;
	tr[1] = sh;
	setImage(bl, tr, f, pW, pH);
	return 0;
      }
    }
      break;
  case 1:
    if (str >> RoughSearchRange
        >> NearSearchRange
        >> NewFeatureRange
        >> MatchThreshold
        >> DistanceGuess
	>> TempletFeature.TrackThreshold) {
      return 0;
    }
    break;
  case 2:
    if (str >> TempletFeature.TanThreshold
	>> TempletFeature.RhoThreshold
	>> TempletFeature.DistanceThreshold) {
        return 0;
    }
    break;
  case 4:
    double sw, sh, f, pW, pH;
    if (str >> sw >> sh >> f >> pW >> pH) {
      double bl[2], tr[2];
      bl[0] = -sw;
      bl[1] = -sh;
      tr[0] = sw;
      tr[1] = sh;
      setImage(bl, tr, f, pW, pH);
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


void HLineHelper::setWidth(Match & m1,Match & m2)
{

  MapHLine *ml=getMapHLine(m1.MatchedFeature->FeatureKey);
  if (ml)
    {
      double d=m1.Measure->Z(1,0)-
	m2.Measure->Z(1,0);
      d/=PixelInfo[0];
      d*=(castPosedHLine(m1.MatchedFeature)->Line.StartPoint(1));
      d*=d;
      if (d>ml->SqWidth)ml->SqWidth=d;
    }
}
void HLineHelper::setCovV(Match & m1)
{
  MapHLine *ml=getMapHLine(m1.MatchedFeature->FeatureKey);
  if (ml)
    {
       if (ml->SqWidth>0)
	{
	  double d=PixelInfo[0]/
	    (castPosedHLine(m1.MatchedFeature)->Line.StartPoint(1));
	  d*=d;
	  d*=ml->SqWidth;
	  double s=m1.Measure->CovV(0,0);
	  if (s>0)
	    {
	      d+=s;
	      d/=s;
	      s=m1.Measure->CovV(4,4);
	      m1.Measure->CovV*=d;
	      if (s>0)m1.Measure->CovV(4,4)=s;
	      else m1.Measure->CovV(4,4)=PixelInfo[0]*PixelInfo[0]/10;
	    }
	  else 
	    {
	      s=m1.Measure->CovV(4,4);
	      m1.Measure->CovV=d;
	      if (s>0)m1.Measure->CovV(4,4)=s;
	      else m1.Measure->CovV(4,4)=PixelInfo[0]*PixelInfo[0]/10;
	    }
	}
    }
}
int HLineHelper::addDistance(Match & mat,double distance)
{
  Measurement &m=*mat.Measure;
  mat.PathDistance=distance;
  //  m.W(0,0)=distance;
  // Line.StartPoint.setXY(m.V.Element);
  // Line.EndPoint.setXY(m.V.Element+2);
  //    Line.StartPoint.getDistance(Line.EndPoint);
  if (m.SensorType!=2)return 1;
  if (!(m.MeasurementType<8))return 1;
  if ((m.W.Columns>4)&&(m.W.Rows<2)){
    setImage(m.W.Element);
  }
  double lenx=m.V(0,0)-m.V(2,0);
  lenx*=lenx;
  double leny=m.V(1,0)-m.V(3,0);
  leny*=leny;
  lenx+=leny;
  mat.Weight=lenx*5E4;
  mat.Thresholds.reallocate(1,2);
  mat.Thresholds(0,1)=RoughSearchRange;
  mat.Thresholds(0,0)=10*RoughSearchRange/sqrt(lenx);
  if (mat.Thresholds(0,0)>.5)mat.Thresholds(0,0)=.5;
  //  double lenx=mat.Measure->W(0,1)/5E4;

  mat.Metric.reallocate(2);
  mat.Metric(0,1)=0;
  mat.Metric(1,0)=0;
  mat.Metric(1,1)=1/(PixelInfo[1]*PixelInfo[2]);
  mat.Metric(0,0)=mat.Metric(1,1)*lenx/10;
  return 0;
  m.CovV*=1;

  //  return;
  //m.CovV.reallocate(5);
   // m.CovV=0;
  //this is to fix pixel center error
  // m.CovV(0,0)=(PixelInfo[1]*PixelInfo[1])*4000;
  //m.CovV(1,1)=(PixelInfo[2]*PixelInfo[2])*4000;
  // m.CovV(2,2)=m.CovV(0,0);
  // m.CovV(3,3)=m.CovV(1,1);
  m.CovV(0,2)=m.CovV(0,0)*(1/5);
  m.CovV(1,3)=m.CovV(1,1)*(1/5);
  m.CovV(2,0)=m.CovV(0,2);
  m.CovV(3,1)=m.CovV(1,3);
  m.CovV(4,4)=(PixelInfo[0]*PixelInfo[0]*.05);
  return 0;
}
/*
void HLineHelper::makeMeasurement(Measurement & m,Matrix & v,double distance)
{

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
unsigned short  HLineHelper::findMergeFeatures(PosedFeature * pfp[2],
				PosedFeatureList * searchList)
{
  PosedHLine *pl[2];
  if (searchList==0)searchList=&Visable;
  double thresholds[2];
  thresholds[0]=.005;
  thresholds[1]= 2*PixelInfo[1];
  double rect[4];
  for (PosedFeatureList *mlist=searchList; mlist->Next; 
       mlist=mlist->Next)
    {
      PosedFeature *pf=mlist->Element;
      if (pf)
	{
	  pl[0]=castPosedHLine(pf);
	  if (pl[0])
	    {      
	      rect[0]=pl[0]->PixelLine.StartPoint.X[0];
	      rect[1]=pl[0]->PixelLine.StartPoint.X[1];
	      rect[2]=pl[0]->PixelLine.EndPoint.X[0];
	      rect[3]=pl[0]->PixelLine.EndPoint.X[1];
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
	     rect[0]+=PixelInfo[1];
	     rect[1]+=PixelInfo[2];
	     rect[2]-=PixelInfo[1];
	     rect[3]-=PixelInfo[2];
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
			   pl[1]=castPosedHLine(pf);
			   MapHLine *ml0=getMapHLine(pl[0]->FeatureKey);
			   MapHLine *ml1=getMapHLine(pl[1]->FeatureKey);
			   unsigned short typ2=ml1->Bdual.Columns;
			   unsigned short typ1=ml0->Bdual.Columns;
			   if ((typ2<3)&&(typ1<3))res=1;
			   else
			     {
			       double d=pl[0]->Line.StartPoint(1)-
				 pl[1]->Line.StartPoint(1);
			       if ((typ2<3)||(typ1<3)){
				 if ((d<.04)&(d>-.04))res=1;
			       }
			       else if ((d<.02)&(d>-.02))res=1;
			     }
			   if (res)
			     {
			        double d=pl[0]->Tangent*pl[1]->Tangent;
			       if (typ1==0){
				 if (typ2==0)
				   {
				     MapHLine *ml=getMapHLine(pl[0]->FeatureKey);
				     int vcount1=0;
				     if (ml)vcount1=ml->VCount;
				     ml=getMapHLine(pl[1]->FeatureKey);
				     int vcount2=0;
				     if (ml)vcount2=ml->VCount;
				     if (vcount1>vcount2)res=5;
				     else  res=1;
				   }
				 else res=1;
			       }
			       else if (typ2==0)res=5;
			       else if (typ1==1){
				 if (typ2==1)
				   {
				     MapHLine *ml=getMapHLine(pl[0]->FeatureKey);
				     int vcount1=0;
				     if (ml)vcount1=ml->VCount;
				     ml=getMapHLine(pl[1]->FeatureKey);
				     int vcount2=0;
				     if (ml)vcount2=ml->VCount;
				     if (vcount1>vcount2)res=6;
				     else  res=2;
				   }
				 else res=2;
			       }
			       else if (typ2==1)res=6;
			       else {
				 MapHLine *ml=getMapHLine(pl[0]->FeatureKey);
				 int vcount1=0;
				 if (ml)vcount1=ml->VCount;
				 ml=getMapHLine(pl[1]->FeatureKey);
				 int vcount2=0;
				 if (ml)vcount2=ml->VCount;
				 if (vcount1>vcount2)res=8;
				 else  res=4;
			       }
			       if (res>4)
				 {
				   pfp[0]=pl[1];
				   pfp[1]=pl[0];
				   res-=4;
				   typ2=typ1;
				 }
			       else{
				 pfp[0]=pl[0];
				 pfp[1]=pl[1];
			       }
			       if (typ2==0)res+=8;
			       else if (typ2==1)res+=16;
			       else res+=32;
			       if (d<0)res+=64;
			       searchList->removeFeature(pl[0]);
			       return res;
			     }
			 }
		     }
	       }
	    }
	}
    }
  return 0;
}


int HLineHelper::getMergeContstraint(Matrix &a, Matrix &e,Matrix & cov,
				      double distance,
				      unsigned short typ, PosedFeature * pl[2])
{
  int dim0=0;
  if (typ&2)dim0=1;
  else if (typ&4)dim0=3;
  int dim1=0;
  if (typ&16)dim1=1;
  else if (typ&32)dim1=3;  
  a.reallocate(dim0,dim0+dim1);
  e.reallocate(dim0,1);
  cov.reallocate(dim0);
  a=0;
  cov=0;
  if (dim0==0)return 0;
  MapHLine *ml0=getMapHLine(pl[0]->FeatureKey);
  MapHLine *ml1=getMapHLine(pl[1]->FeatureKey);
  ml0->calcTangent();
  ml1->calcTangent();
  a(0,0)=1;///ml0->Length;;
  e(0,0)=atan2(-ml0->Tangent[0],ml0->Tangent[1]);
  if (typ&64)
  {
    a(0,dim0)=1;///ml1->Length;
    e(0,0)+=atan2(-ml1->Tangent[0],ml1->Tangent[1]);
  } 
 else {
   a(0,dim0)=-1;///ml1->Length;
   e(0,0)-=atan2(-ml1->Tangent[0],ml1->Tangent[1]);
 }
  e(0,0)*=-M_SQRT1_2;
  cov(0,0)=e(0,0)*e(0,0);
  if (dim0==1)return 0;
  Cure::Matrix p0;
  Cure::Matrix p1;
  ml0->getP(p0);
  ml1->getP(p1);
  a(1,1)=1;
  a(2,2)=1;
  e(1,0)=p0(1,0);
  e(2,0)=p0(2,0);
  if (typ&64)
  {
    a(1,dim0+1)=1;
    a(2,dim0+2)=1;
    e(1,0)+=p1(1,0);
    e(2,0)+=p1(2,0);
  } 
 else {
   a(1,dim0+1)=-1;
   a(2,dim0+2)=-1;
   e(1,0)-=p1(1,0);
   e(2,0)-=p1(2,0);
 }
  cov(1,1)=e(1,0)*e(1,0);
  cov(2,2)=e(2,0)*e(2,0);
  //  e(1,0)*=2;
  //e(2,0)*=2;
  return 0;
}

void HLineHelper::match(Match *matches,int n,
			double distance,Cure::Pose3D &pcam, LongList &mapFeats)
{
  //  visFeatures.clear();
  MapHelper::match(matches,n,distance,pcam,mapFeats);
  //THis is special for HLINES
  //Here we increase CovV if line has SqWidth>0 and we set Width if two
  //measurments match the line. 
  for (int i=0; i<n; i++)
    {
      if (!(matches[i].MatchDistance<0))
	if (matches[i].MatchedFeature)
	  {
	    long key=matches[i].MatchedFeature->FeatureKey;
		for (int j=i+1;j<n;j++)
		  if (!(matches[j].MatchDistance<0))
		    if (matches[j].MatchedFeature)
		      if (key==matches[j].MatchedFeature->FeatureKey)
			setWidth(matches[i], matches[j]);
		setCovV(matches[i]);
	  }
    }
 
}
