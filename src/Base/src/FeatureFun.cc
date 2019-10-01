// = RCSID
//    $Id: FeatureFun.cc ,v 1.1 2004/04/1 
//
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "FeatureFun.hh"

using namespace Cure;

namespace Cure {
namespace FeatureFcn {

int lineInRectangle(double start[2],double end[2], double  bottomLeft[2],
	   double topRight[2])
{
  if (((start[0]) < topRight[0])&& 
      ((start[0]) > bottomLeft[0])&&
      ((start[1]) < topRight[1])&&
      ((start[1]) > bottomLeft[1]))return 1;
  if (((end[0]) < topRight[0])&& 
      ((end[0]) > bottomLeft[0])&&
      ((end[1]) < topRight[1])&&
      ((end[1]) > bottomLeft[1]))return 1;
  double m, b, temp;
  m=(((start[1])- (end[1])))/
    ((((start[0])- (end[0]))));
  b=((start[1]))-m*((start[0]));
  if (((topRight[0]<(start[0]))&&
       (topRight[0]>(end[0])))||
      ((topRight[0]>(start[0]))&&
       (topRight[0]<(end[0]))))
    { 
      temp=b+m*topRight[0];
      if ((temp <topRight[1])&&(temp>bottomLeft[1])) return 1;
    }
  if (((bottomLeft[0]<(start[0]))&&
       (bottomLeft[0]>(end[0])))||
      ((bottomLeft[0]>(start[0]))&&
       (bottomLeft[0]<(end[0]))))
    { 
      temp=b+m*bottomLeft[0];
      if ((temp <topRight[1])&&(temp>bottomLeft[1])) return 1;
    }
  if (m==0)  m=.0001;
  if (((topRight[1]<(start[1]))&&
       (topRight[1]>(end[1])))||
      ((topRight[1]>(start[1]))&&
       (topRight[1]<(end[1]))))
    { 
      temp=(topRight[1]-b)/m;
      if ((temp <topRight[0])&&(temp>bottomLeft[0])) return 1;
    }
  if (((bottomLeft[1]<(start[1]))&&
       (bottomLeft[1]>(end[1])))||
      ((bottomLeft[1]>(start[1]))&&
       (bottomLeft[1]<(end[1]))))
    { 
      temp=(bottomLeft[1]-b)/m;
      if ((temp <topRight[0])&&(temp>bottomLeft[0])) return 1;
    }
  return 0;
}

int pointInRectangle(double center[2], double  bottomLeft[2],
	   double topRight[2])
{
  if (center[0] > topRight[0]) return 0; 
  if (center[0] < bottomLeft[0]) return 0;
  if (center[1] < bottomLeft[1]) return 0;
  if (center[1] > topRight[1]) return 0;
  return 1;
}
int near(MapFeatureList& near,MapFeatureList& features,
	       const double distance,const double center[2])
{
  near.clean();
  double bottomLeft[2],topRight[2];
  bottomLeft[0]= (center[0]-distance);
  topRight[0]= (center[0]+distance);
  bottomLeft[1]=(center[1]-distance);
  topRight[1]= (center[1]+distance);
  for (MapFeatureList *flist=&features;flist->Next;flist=flist->Next)
    {
      MapFeature *f=flist->Element;  
      if (f->inside(bottomLeft,topRight)) {
	near.add(f);
      }
    }
  return near.Count;
}
void listFeatures(MapFeatureList& features,MapBank& b, LongList &llist)
{
  features.clean();
  for (LongList *l=&llist; l->Next; l=l->Next)
    {
      MapFeature *mf =getMapFeature(&b,l->Element);
      if (mf)features.add(mf);
    }
}
void listLongs(MapFeatureList & features, LongList& llist)
{
  llist.clear();
  for (MapFeatureList *l=&features; l->Next; l=l->Next)
    {
      llist.add(l->Element->Key);
    }
}

} // namespace FeatureFcn
} // namespace Cure
