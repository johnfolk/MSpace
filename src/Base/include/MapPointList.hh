// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_MAPPOINTLIST_H
#define CURE_MAPPOINTLIST_H

#include "MapPoint.hh"

namespace Cure{

  /**
   * List for map points, i.e. 2D or 3D points that are being
   * estimated in the map.
   *
   * @author John Folkesson
   * @see
   */
  class MapPointList
  {
  public:
    int count;
    MapPoint * element; 
    MapPointList *next;
    MapPointList *prev;  
  public:
    MapPointList();
    ~MapPointList();
    void clean();
    int add(MapPoint *mf);
    int addUnique(MapPoint *mf);
    int remove(long k);
    int removePoint(MapPoint *mf);
    MapPoint * get(long n);
  };
}
#endif
