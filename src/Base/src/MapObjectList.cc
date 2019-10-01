// = RCSID
//    $Id: MapObjectList.cc ,v 1.1 2004/04/1 
//
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "MapObjectList.hh"
using namespace Cure;
int MapObjectList::add(MapObject *pa)
{
  if (Next==0)
    {
      Element=pa;
      Next=new MapObjectList();
      return 1; 
    }
  return Next->add(pa);
}
int MapObjectList::remove(MapObject *pa)
{
  if (Next==0) return 1;
  if (Element==pa)
    {
      Element=Next->Element;
      MapObjectList *l=Next;
      Next=Next->Next;
      l->Next=0;
      delete l;
      return 0;
    }
  return Next->remove(pa);
}
int MapObjectList::removeKey(long key)
{
  if (Next)
    {
      if (Element->Key==key)
	return remove(Element);
      else return Next->removeKey(key);
    }
  return 1; 
}


