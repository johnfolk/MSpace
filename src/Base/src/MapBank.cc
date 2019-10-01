// = RCSID
//    $Id: MapObjectList.cc ,v 1.1 2004/04/1 
//
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "MapBank.hh"

//using namespace Cure;

namespace Cure {

MapBank::MapBank(int depth, unsigned long id)
{
  BankID=id;
  NumberOfBanks=0;
  Banks=0;
  InvBanks=0;
  NextIndex=1;
  Size=1;
  Size=Size<<depth;
  Mask=Size-1;
  Objects=new MapObjectList[Size];
  for (int i=0; i<Size;i++)
    Objects[i].Element=0;
}
MapBank::~MapBank()
{
  bool warn=true;
  for (int i=0; i<Size;i++)
    if(Objects[i].Next!=0)
      {
	if (warn)
	  {
	    long c=0;
	    for (int j=0; j<Size;j++)
	      c+=Objects[j].count();
	    warn=false;
	    std::cerr<<"MEMORY LEAK: MapBank being deleted before its Objects. ";
	    std::cerr<<"The Objects should be deleted first."<<std::endl;  
   	    std::cerr<<"Here is one of the "<<c<<" that were left:"<<std::endl;  
	    Objects[i].Element->print();
	  } 
      }
  delete [] Objects;
  for (unsigned short i=0;i<NumberOfBanks;i++)
    if (Banks[i])delete Banks[i];
  if (Banks)delete []Banks;
  for (unsigned short i=0;i<NumberOfBanks;i++)
    if (InvBanks[i])delete InvBanks[i];
  if (InvBanks)delete []InvBanks;
  NumberOfBanks=0;
  Banks=0;
  InvBanks=0;
  std::cerr<<"MapBank deleted."<<std::endl;
}
void MapBank::add(MapObject *pa, long index)
{
  
  if (index>-1){
    if (index>=NextIndex)
      NextIndex=index;
    else if (!(getMapObject(index)))
      {
	long i=(index&Mask);
	pa->Key=index;
	Objects[i].add(pa);
	return;
      }
  }
  index=(NextIndex&Mask);
  pa->Key=NextIndex;
  Objects[index].add(pa);
  NextIndex++;
}
void MapBank::setDataSet(DataSet &ds, MapObjectMaker **makers, short n, LongList *keys)
{
  long k=0;
  TimestampedData  *td=ds.getTPointer(k);
  while (td){
    k++;
    GenericData *gd=td->narrowGenericData();
    td=ds.getTPointer(k);
    if (gd){
      short objectType;
      short objectSubType;
      if (!getObjectInfo(*gd,objectType,objectSubType))
	for (int i=0;i<n;i++)
	  if (makers[i]->getObjectType()==objectType)
	    if (makers[i]->getObjectSubType()==objectSubType)
	      {
		MapObject *mo=makers[i]->makeMapObject(*gd);
		if ((keys)&&(mo))keys->addUnique(mo->Key);
		i=n;
	      }
    }
    
  }
}
}
