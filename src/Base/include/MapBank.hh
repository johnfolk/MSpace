// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_MAPBANK_HH
#define CURE_MAPBANK_HH

#include "iostream"
#include "MapObjectList.hh"
#include "LongList.hh"
#include "mapDefines.hh"
#include "LongArray.hh"
#include "DataSet.hh"
namespace Cure{

  
  /**
   * Object that manages map objects. The idea is that you should not
   * store pointers to the objects in the map but instead indices that
   * you use to look the feature pointer up. This way you can delete
   * objects from the map without having to tell everyone about it. If
   * your index does not work anymore it means that the feature does
   * not exist anymore. Can help avoid segmentation faults.
   *
   * @author John Folkesson
   * @see AddressBank
   */
  class MapBank
  {
  public:
    /**
     * This identifies the bank when transfering objects between banks.
     * these should be assigned unique starting from 0 in ascending order.
     */
    unsigned short BankID;
    MapObjectList *Objects;
    int Size;
    long Mask;
    long NextIndex;
  protected:
    /** 
     * These hash tables associate a key from some other external bank
     * with a key from this bank.  No attempt need be made to keep the
     * banks in sync but the mapobjects need to be added and
     * associated before trying to reference them.  One can associate
     * more than one external key to the same internal key but if more
     * than one key with the same bankid is associated to the same key
     * here then the inverse hash will only find the last one added.
     * Thus some keys on Banks could end up associated with deleted
     * features.  This might not matter as one will just get null when
     * looking them up. However the Banks array might end up with lots
     * of garbage.
     */
    LongArray **Banks,**InvBanks;
    unsigned short NumberOfBanks;
  public:
    MapBank(int depth=7, unsigned long id=0);
    virtual ~MapBank();
    void add(MapObject *pa, long key=-1);
    void remove(MapObject *pa){
      long key=pa->Key;
      long inx=((key)&Mask);
      Objects[inx].remove(pa);
      pa->Key=-1;
      long ind=-1;
      for (int i=0;i<NumberOfBanks;i++){
	if (InvBanks[i]->getHash(key,ind)){
	  InvBanks[i]->removeHash(key);
	  Banks[i]->removeHash(ind);
	}
      }
    }
    MapObject * getMapObject(long key){
      if (key<0)return 0;
      long inx=(key&Mask);
      return Objects[inx].get(key);
    }   
    void associate(unsigned short bankid, long key, long thiskey){
      if (bankid==BankID)return;
      if (key==-1)return;
      if (bankid>=NumberOfBanks){
	if (bankid>0x7FFF)return;
	else addBankHash(bankid);
      }
      Banks[bankid]->putHash(key,thiskey);
      InvBanks[bankid]->putHash(thiskey,key);
    }
    /**
     * This allows associating MapBank Keys from other banks with
     * those from this one.
     *
     * @bankid the id of the bank to be associated with use numbers 0,1,2...
     * in order.
     * @param depth the depth of the new hash tables, (ie the size is
     * 2^depth).  default depth=-1 sets the depth to that of the MapBank.
     */
    void addBankHash(unsigned short bankid,short depth=-1){
      if (bankid>=NumberOfBanks){
	LongArray **b=Banks;
	Banks=new LongArray*[bankid+1];
	memcpy(Banks,b,NumberOfBanks*sizeof(LongArray*));
	for (unsigned short i=NumberOfBanks;i<(bankid+1);i++)
	  Banks[i]=0;
	delete[]b;
	b=InvBanks;
	InvBanks=new LongArray*[bankid+1];
	memcpy(InvBanks,b,NumberOfBanks*sizeof(LongArray*));
	for (unsigned short i=NumberOfBanks;i<(bankid+1);i++)
	  InvBanks[i]=0;
	NumberOfBanks=bankid+1;
	delete[]b;
      }
  
      if (!Banks[bankid]){
	Banks[bankid]=new LongArray();
	InvBanks[bankid]=new LongArray();
	if (depth==-1){
	  depth=-1;
	  unsigned long s=Size;
	  while (s){
	    s=(s>>1);
	    depth++;
	  }
	}
	if (depth<0)depth=0;
	Banks[bankid]->setHashDepth((long)depth);
	InvBanks[bankid]->setHashDepth((long)depth);
      }
    }
    
    MapObject * getMapObject(unsigned short bankid, long key){
      if (bankid==BankID)return getMapObject(key);
      if (bankid>=NumberOfBanks)return 0;
      if (Banks[bankid])
	if (Banks[bankid]->getHash(key,key))return getMapObject(key);  
      return 0;
    }
    
    MapObject * getObjectType(long key, unsigned short type){
      MapObject *mo=getMapObject(key);
      if (mo) return mo->getObjectType(type);
    }
    MapObject * getObjectType(unsigned short bankid,long key, unsigned short type){
      MapObject *mo=getMapObject(bankid,key);
      if (mo) return mo->getObjectType(type);
    }
    /**
     * This will store the data from a list of map objects in a
     * DataSet object that can be sent to another module using a
     * different MapBank.
     *
     * @param keys the keys of the mapobjects to store
     * @param ds the MapObjects data is put in this
     */
    void getDataSet(LongList &keys, DataSet &ds){
      ds.setSetSize(keys.count());
      long k=0;
      for (LongList *llist=&keys;llist->Next; llist=llist->Next){
	MapObject *mo=getMapObject(llist->Element);
	if (mo){
	  ds.setup(k,GENERIC_TYPE);
	  GenericData *gd=ds.getTPointer(k)->narrowGenericData();
	  mo->get(*gd);
	  k++;
	}
      }
      ds.setSetSize(k);
    }
    /**
     * This creates MapObjects for all objects stored on a DataSet and 
     * if they have not already been set.
     * 
     * 
     * @param ds the objects are stored as DataSets within ds.
     * @param makers the MapObjects are created by these.
     * @param n the length of the makers array.
     * @param keys optionally returns the keys of the objects here
     */
    void setDataSet(DataSet &ds, MapObjectMaker **makers, short n, 
		    LongList *keys=0);
  };
}
#endif
