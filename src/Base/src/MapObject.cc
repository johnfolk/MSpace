// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#include "MapObject.hh"
#include "MapBank.hh"
using namespace Cure;

MapObject::MapObject(MapBank *b){
  Key=-1;
  ID=-1;
  Bank=b;
  if (Bank)Bank->add(this);
}
MapObject::~MapObject()
{
  if (Bank) Bank->remove(this);
  Bank=0;
}
void MapObject::get(GenericData &gd){
  gd.GenericType=getObjectType();
  gd.setShortDataSize(1, 6);
  gd.ShortData(0,0)=getObjectSubType();
  gd.ShortData.setLong(Key,0,2);
  gd.ShortData.setLong(ID,0,4);
  gd.ShortData(0,1)=Bank->BankID;
}
int MapObject::set(GenericData &gd)
{
   if (gd.GenericType!=getObjectType())return 1;
   if (gd.ShortData.Rows<1)return 1;
   if (gd.ShortData.Columns<6)return 1;
   if (gd.ShortData(0,0)!=getObjectSubType())return 1;
   long k=gd.ShortData.getLong(0,2);
   ID=gd.ShortData.getLong(0,4);
   if (Bank)Bank->associate(gd.ShortData(0,1),k,Key);
   return 0;
}
MapObject * MapObjectMaker::makeMapObject(GenericData &gd)
{
  if (!Bank)return 0;
  short bankid;
  long key;
  if (getBankInfo(gd, bankid,key))return 0;
  MapObject *mo=Bank->getMapObject(bankid,key);
  short objectType, objectSubType;
  if (getObjectInfo(gd,objectType,objectSubType))return 0;
  if (objectType!=getObjectType())return 0;
  if (objectSubType!=getObjectSubType())return 0;
  if (mo){
    if (mo->getObjectType()!=objectType)return 0;
    else if (mo->getObjectSubType()!=objectSubType)return 0;
  }else mo=makeMapObject();
  if (mo) mo->set(gd);
  Bank->associate(bankid,key,mo->Key);
  return mo;
}
