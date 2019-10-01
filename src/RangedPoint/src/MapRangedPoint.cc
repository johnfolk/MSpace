//
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2006 John Folkesson
//    

#include "MapRangedPoint.hh"
using namespace Cure;

MapRangedPoint::MapRangedPoint():MapFeature(0)
{
  FullDim=3;
  NumberPoints=1;
  Number3D=1;
  Points=AllocatedPoints;
  Index=AllocatedIndex;
  Type=MAPRANGEDPOINT_TYPE;
  init();
  CastPtr=0;
}
MapRangedPoint::MapRangedPoint(MapRangedPoint *wp,
				   MapBank *b):MapFeature(b)
{
  FullDim=3;
  NumberPoints=1;
  Number3D=1;
  Points=AllocatedPoints;
  Index=AllocatedIndex;
  Type=MAPRANGEDPOINT_TYPE;
  init();
  if (wp)
     {
       CastPtr=wp->CastPtr;
       DistanceThreshold=wp->DistanceThreshold;
     }
  else
    {
      CastPtr=0;
      DistanceThreshold=20;
    } 
}

int  MapRangedPoint::extend()
{
  if (Bdual.Columns==3)return 0;
  if (1)//we need some test
    {
      Bdual.reallocate(3,3);
      Bdual=1;
      return 3;
    }
  return 0;  
}

void  MapRangedPoint::forceExtend()
{
  if (Bdual.Columns==3)return;
  Bdual.reallocate(3,3);
  Index[0]=-1;
  Bdual=1;
}

unsigned short MapRangedPoint::getMeasurementType(unsigned short type)
{
  if (!type)return type;
  if (!(type&0x10))type=0;
  type=type&0x6;
  if (Bdual.Columns==3) {
    return type;
  }
  else type=0;
  return type;  
}

int  MapRangedPoint::addInfo(const Cure::Matrix & v,
			      Transformation3D & map2info,
			      const int type)
{
  if (Bdual.Columns==3)return 0;
  return 0;
}

int  MapRangedPoint::addInfo(const Cure::Matrix & v,
			      const int type)
{
  return 0;
}
int MapRangedPoint::merge(MapRangedPoint *mf, unsigned short typ)
{
  if (!mf)return -1;
  return 0;
}

void MapRangedPoint::getB(Cure::Matrix &b)const
{
  b.reallocate(Bdual.Columns,Bdual.Rows);
  b=1; 
}
bool MapRangedPoint::getC(MapRangedPoint *mw,Cure::Matrix &c)const
{
  if ((Bdual.Columns!=3)||(mw->Bdual.Columns!=3))return false;
  c.reallocate(3);
  c=1;
  return true; 
}
bool MapRangedPoint::testMatch(MapRangedPoint *mw, double tolerance)
{
  if ((Bdual.Columns!=3)||(mw->Bdual.Columns!=3))return false;
  return true;
}

