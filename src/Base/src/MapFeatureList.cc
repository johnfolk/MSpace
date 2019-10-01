// = RCSID
//    $Id: MapFeatureList.cc ,v 1.1 2004/04/1 
//
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "MapFeatureList.hh"

using namespace Cure;

MapFeatureList::MapFeatureList()
{
  Element=0;
  Next=0;
  Count=0;
  Prev=0;
}

MapFeatureList::~MapFeatureList()
{
  if (Next!=0) delete Next;
}

 void MapFeatureList::clean()
{
  if (Next!=0) delete Next;
  Next=0;
  Prev=0;
  Count=0; 
  Element =0;
}


int MapFeatureList::add(MapFeature *w)
{
  if (Next==0)
    {
      Next=new MapFeatureList();
      Next->Prev=this;
  }
  Count++;
  if (Count>1)
    {
      return Next->add(w);
    }
  Element=w;
  return 1;
}

int MapFeatureList::addUnique(MapFeature *w)
{
  if (Element==w)return 0;
  if (Next==0)
    {
      Next=new MapFeatureList();
      Next->Prev=this;
    }
  
  if (Count>0)
    {
      if (Next->addUnique(w))
	{
	  Count++;
	  return 1;
    	}
      return 0;
    }
  Count++;
  Element=w;
  return 1;
}
int MapFeatureList::remove(long k)
{
  if (k<Count)
    { 
      if (k>0) 
	{
	  int i=Next->remove(k-1);
	  Count-=i;
	  return i;
	}
     else 
	{
	  if(Count>1)
	    {
	      Element=(Next->get(0));
	      Next->remove(0);
	    }
	  else 
	    {
	      delete Next;
	      Element=0;
	      Next=0;
	    }
	}
      Count--;
      return 1;
    }
  return 0;
}
int MapFeatureList::removeFeature(MapFeature *n)
{
  if (Count==0)return 0;
  int m=1;
  while (m>0)
    {
      m=-1;
      int i=0;
      for (MapFeatureList *nList=this; nList->Next; i++,nList=nList->Next) 
	{
	  if (nList->Element==n)
	    {
	      m=i;
	    }
	}
      if (m>-1)  remove (m);
    }
  return 1;
}
Cure::MapFeature * MapFeatureList::get(long n)
{
  if(n<Count) 
    {
      if (n==0) return Element;
      if(n>0)
	return Next->get(n-1);
    }
  return 0; 
}
