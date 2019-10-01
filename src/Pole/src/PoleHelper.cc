// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "PoleHelper.hh"


using namespace Cure;
void PoleHelper::matchPoles(PosedFeatureList & matches, Cure::Point2D & pt,
			    double searchlimit,	PosedFeatureList * poles)
{

  double rect[4];
  rect[0]=pt(0)-searchlimit;
  rect[1]=pt(1)-searchlimit;
  rect[2]=pt(0)+searchlimit;
  rect[3]=pt(1)+searchlimit;
  double z[2],thresholds[3];
  if (poles==0)poles=&Visable;
  pt.getPolarValues(z);
  thresholds[0]=searchlimit/z[1];
  thresholds[1]=searchlimit;
  thresholds[2]=0;
  for (PosedFeatureList *mlist=poles; mlist->Next; 
       mlist=mlist->Next)
    {
      PosedPole *pf=castPosedPole(mlist->Element);
      if (pf)
	{
	  if(FeatureFcn::pointInRectangle(pf->Center.X, rect,rect+2))
	    {
	      if (!(pf->roughMatch(z,thresholds)))
		matches.add(pf);
	    }
	}
    }
}

double PoleHelper::matchPole(PosedPole & pp, Cure::Point2D & pt,
			       Matrix & sigma,
			       PosedFeatureList * poles)
{
  if (poles==0)poles=&Visable;
  double z[2];
  pt.getPolarValues(z);
  int rows=sigma.Rows;
  if (rows<2)
    {
      rows=2;
      sigma.reallocate(2);
      sigma=1;
      sigma(0,0)=.05;
    }
  else
    {
      sigma.Rows=2;
      sigma.Columns=2;
    }
  double dist=1000;
  for (PosedFeatureList *mlist=poles; mlist->Next; 
       mlist=mlist->Next)
    {
      PosedPole *pf=castPosedPole(mlist->Element);
      if (pf)
	{
	  {
	    
	    double d=pf->fineMatch(z,sigma);
	    if ((d<dist)&&(d!=-1))
	      {
		pp=*pf;
		dist=d;
	      }
	  }
	}
    }
  sigma.Rows=rows;
  sigma.Columns=rows;
  if (dist<1000)return dist;
  return -1;
}
double PoleHelper::makeMeasurement(double pixels[3],double distance,
				    Measurement & m)
{
  return 0;
}
