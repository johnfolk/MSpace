// = RCSID
//    $Id: MapPoint.cc ,v 1.1 2004/04/1 
//
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "MapPoint.hh"

using namespace Cure;

MapPoint::MapPoint(MapBank *b): MapObject(b),Point3D()
{
  Index=-1;
  Local_map=-1;
}
MapPoint::MapPoint(const MapPoint& p):MapObject(p.Bank)
{
  Index=-1;
  Local_map=-1;
  Vector3D::operator=(p);
} 
MapPoint::MapPoint(const Point3D& p, MapBank *b): MapObject(b)
{

  Index=-1;
  Local_map=-1;
  Vector3D::operator=(p);
} 
MapPoint::~MapPoint()
{
  while(Features.Next)
    {
      MapObject *f=Features.Element;
      f->removeObject(this);
      Features.remove(f);
    }
}
void MapPoint::write(std::fstream &fs )
{
  fs<<Key<<" "<<ID<<" ";
  fs<<Index<<" "<<Local_map<<" ";
  fs<<X[0]<<" "<<X[1]<<" "<<X[2]<<"\n";
}
void MapPoint::read(std::fstream &fs )
{
  fs>>ID>>Index>>Local_map;
  fs>>X[0]>>X[1]>>X[2];
}
int MapPoint::addFeature(MapObject *f)
{
  return Features.add(f);
}
int MapPoint::removeFeature(MapObject *f)
{
  return Features.remove(f);
}
void MapPoint::initializeFromPixels(Cure::Transformation3D &t,
				    const double pixels[2],
				    double focallength,
				    double distanceguess)
{
  Cure::Vector3D v;
  double d=(distanceguess-focallength)/focallength;
  v(0)=-(pixels[0]*d);
  v(2)=-(pixels[1]*d);
  v(1)=distanceguess;
  t.invTransform(v, *this);
}
