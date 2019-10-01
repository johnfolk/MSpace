// = RCSID
//    $Id: PosedFeatureList.cc ,v 1.1 2004/04/1 
//
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "PosedFeatureList.hh"
 
using namespace Cure;

PosedFeatureList::PosedFeatureList()
{
  Element=0;
  Next=0;
  Count=0;
  Prev=0;
}

PosedFeatureList::~PosedFeatureList()
{
  if (Next!=0) delete Next;
}

void PosedFeatureList::clean()
{
  if (Next!=0) delete Next;
  Next=0;
  Prev=0;
  Count=0; 
  Element =0;
}
void PosedFeatureList::clear()
{
  if (Next!=0) 
    {
      delete Element;
      Next->clear();
      delete Next;
    }
  Next=0;
  Prev=0;
  Count=0; 
  Element =0;
}

int PosedFeatureList::add(PosedFeature *w)
{
  if (Next==0)
    {
      Next=new PosedFeatureList();
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

int PosedFeatureList::addUnique(PosedFeature *w)
{
  if (Element==w)return 0;
  if (Next==0)
    {
      Next=new PosedFeatureList();
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
int PosedFeatureList::remove(long k)
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
int PosedFeatureList::removeFeature(PosedFeature *n)
{
  if (Count==0)return 0;
  int m=1;
  while (m>0)
    {
      m=-1;
      int i=0;
      for (PosedFeatureList *nList=this; nList->Next; i++,nList=nList->Next) 
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
PosedFeature *  PosedFeatureList::keyGet(long key)
{
  if (Next)
    {
      if (Element->FeatureKey==key)return Element;
      return Next->keyGet(key);
    }
  return 0;
}
PosedFeature * PosedFeatureList::get(long n)
{
  if(n<Count) 
    {
      if (n==0) return Element;
      if(n>0)
	return Next->get(n-1);
    }
  return 0; 
}
void PosedFeatureList::listType(int typ, PosedFeatureList & result)
{
  if (Next)
    {
      MapFeature * mf=Element->getFeature();
      if (mf)
	if (mf->Type==typ)result.add(Element);
      Next->listType(typ,result);
    }
}
