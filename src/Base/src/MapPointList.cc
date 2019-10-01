// = RCSID
//    $Id: MapPointList.cc ,v 1.1 2004/04/1 
//
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "MapPointList.hh"

using namespace Cure; 

MapPointList::MapPointList()
{
  element=0;
  next=0;
  count=0;
  prev=0;
}

MapPointList::~MapPointList()
{
  if (next!=0) delete next;
}

void MapPointList::clean()
{
  if (next!=0) delete next;
  next=0;
  prev=0;
  count=0; 
  element =0;
}

int MapPointList::add(MapPoint *w)
{
  if (next==0)
    {
      next=new MapPointList();
      next->prev=this;
  }
  count++;
  if (count>1)
    {
      return next->add(w);
    }
  element=w;
  return 1;
}

int MapPointList::addUnique(MapPoint *w)
{
  if (element==w)return 0;
  if (next==0)
    {
      next=new MapPointList();
      next->prev=this;
    }
  
  if (count>0)
    {
      if (next->addUnique(w))
	{
	  count++;
	  return 1;
    	}
      return 0;
    }
  count++;
  element=w;
  return 1;
}
int MapPointList::remove(long k)
{
  if (k<count)
    { 
      if (k>0) 
	{
	  int i=next->remove(k-1);
	  count-=i;
	  return i;
	}
     else 
	{
	  if(count>1)
	    {
	      element=(next->get(0));
	      next->remove(0);
	    }
	  else 
	    {
	      delete next;
	      element=0;
	      next=0;
	    }
	}
      count--;
      return 1;
    }
  return 0;
}
int MapPointList::removePoint(MapPoint *n)
{
  if (count==0)return 0;
  int m=1;
  while (m>0)
    {
      m=-1;
      int i=0;
      for (MapPointList *nList=this; nList->next; i++,nList=nList->next) 
	{
	  if (nList->element==n)
	    {
	      m=i;
	    }
	}
      if (m>-1)  remove (m);
    }
  return 1;
}
MapPoint * MapPointList::get(long n)
{
  if(n<count) 
    {
      if (n==0) return element;
      if(n>0)
	return next->get(n-1);
    }
  return 0; 
}
