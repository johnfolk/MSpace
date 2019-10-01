// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_MAPOBJECTLIST_H
#define CURE_MAPOBJECTLIST_H

#include "MapObject.hh"

namespace Cure{

  /**
   * List of map objects 
   * 
   * @author John Folkesson
   */
  class MapObjectList
  {
  public:
    MapObject *Element;
    MapObjectList *Next;
  public:
    MapObjectList(){Next=0;}
    ~MapObjectList(){
      if (Next)delete Next;
    }
    int add(MapObject *mo);
    int remove(MapObject *mo);
    int removeKey(long key);
    MapObject * get(long key){
      if (Next)
	{
	  if (Element->Key==key)
	    return Element;
	  else return Next->get(key);
	}
      return 0; 
    }

    int count()
    {
      if (Next)return (1+Next->count());
      return 0; 
    }
  };

} // namespace Cure

#endif
