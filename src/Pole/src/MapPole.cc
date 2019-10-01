// = RCSID
//    $Id: MapPole.cc ,v 1.1 2004/04/1 
//
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "MapPole.hh"
using namespace Cure;
using namespace MatrixStuff;
MapPole::MapPole():MapFeature(0)
{
  FullDim=3;
  NumberPoints=1;
  Number2D=1;
  NumberScalars=1;
  Scalars=&Radius;
  Points=AllocatedPoints;
  Index=AllocatedIndex;
  Type=MAPPOLE_TYPE;
  Radius=0;
  InfoInitialized=0;
  init();
  CastPtr=0;
  CountThreshold=25;
  RadiusThreshold=.5;
  DistanceThreshold=0;
}
MapPole::MapPole(MapPole * wp,MapBank *b):MapFeature(b)
{
  FullDim=3;
  NumberPoints=1;
  Number2D=1;
  NumberScalars=1;
  Scalars=&Radius;
  Points=AllocatedPoints;
  Index=AllocatedIndex;
  Type=MAPPOLE_TYPE;
  Radius=0;
  InfoInitialized=0;
  init();
  if (wp)
     {
       CastPtr=wp->CastPtr;
       CountThreshold=wp->CountThreshold;
       RadiusThreshold=wp->RadiusThreshold;
       DistanceThreshold=wp->DistanceThreshold;
     }
  else
    {
      CastPtr=0;
      CountThreshold=10;
      RadiusThreshold=.5;
      DistanceThreshold=0;
    }
}
void MapPole::setCenter(MapPoint *val)
{
  addPoint(val,0);
}
void MapPole::initailizeFromPoint(Transformation3D &t,Point2D & pt, double r)
{
  if (!Points[0])
    {
      MapPoint *p=new MapPoint(Bank);
      setCenter(p);
    }
  t.invTransform(pt, *Points[0]);
  setRadius(r);
}

void MapPole::setRadius(double r)
{
  Radius=r;
}
int  MapPole::extend()
{
  if (Bdual.Columns==3)return 0;
  if (Bdual.Columns==0)
    {
      if (InfoInitialized)
	{
	  Bdual.reallocate(3,3);
	  Bdual=1;
	  return 1;
	}
    }
  return 0;  
}

int  MapPole::updateInfo(const Cure::Matrix & v, 
			 Cure::Matrix & covAdjust,const int type)
{
  if (Bdual.Columns==3)return 0;
  if (v.Rows==0)return 0;
  double mindistance=v(0,2)-DistanceThreshold;
  Cloud.merge(v,mindistance);
  if (Bdual.Columns==0)
    {
      Matrix x(3,1);
     
      if (Cloud.fitCircle())
      {
	x(0,0)=Cloud.Center[0];
	x(1,0)=Cloud.Center[1];
	x(2,0)=Cloud.Radius;
	if ((Cloud.CloudCount>CountThreshold)&&
	    (Cloud.VarRadius<RadiusThreshold))
	  {
	    Info(0,0)=Cloud.Cov[0];
	    Info(0,1)=Cloud.Cov[1];
	    Info(0,2)=Cloud.Cov[2];
	    Info(1,0)=Cloud.Cov[3];
	    Info(1,1)=Cloud.Cov[4];
	    Info(1,2)=Cloud.Cov[5];
	    Info(2,0)=Cloud.Cov[6];
	    Info(2,1)=Cloud.Cov[7];
	    Info(2,2)=Cloud.Cov[8];
	    if (!Info.invert())
	      InfoInitialized=1;
	  }
      }
      else
	{
	  x=0;
	  for (int i=0; i<v.Rows; i++)
	    {
	      x(0,0)+=v(i,0);
	      x(1,0)+=v(i,1);
	    }
	  x(0,0)/=v.Rows;
	  x(1,0)+=v.Rows;
	}
      setX(x);
    }
 return 0;
}
  
void MapPole::getB(Cure::Matrix &b)const
{
  b.reallocate(Bdual.Columns);
  b=1; 
}
