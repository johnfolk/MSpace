// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_MAPFEATURELIST_H
#define CURE_MAPFEATURELIST_H

#include "MapObject.hh"
#include "MapBank.hh"
#include "LongList.hh"
namespace Cure{

  // Forward declaration
  class MapFeature;

  /**
   * List for map features
   *
   * @author John Folkesson
   * @see
   */
  class MapFeatureList
  {
  public:
    int Count;
    MapFeature * Element; 
    MapFeatureList *Next;
    MapFeatureList *Prev;  
  public:
    MapFeatureList();
    ~MapFeatureList();
    /**
     * Forget ervything on the list without deleting the features.
     */
    void clean();
    int add(MapFeature *mf);
    int addUnique(MapFeature *mf);
    int remove(long k);
    int removeFeature(MapFeature *mf);
    MapFeature * get(long n);
  };  // class MapFeatureList

} // namespace Cure

#endif
