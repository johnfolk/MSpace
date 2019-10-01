// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_MAPOBJECT_H
#define CURE_MAPOBJECT_H
#include "GenericData.hh"
#include "mapDefines.hh"
namespace Cure{
// Forward declarations
  class MapBank;
  class MapFeature;
  class MapPoint;

  //namespace RelativeBase{

  //}
  
  /**
   * Base class for objects in the map
   *
   * @author John Folkesson
   * @see
   */
  class MapObject
  {
  public:
    long ID;
    long Key;
    MapBank * Bank;
  public:
    MapObject(MapBank *b=0);
    virtual ~MapObject();
    virtual int removeObject(MapObject*){return 0;}
    virtual void narrow(){}  
    virtual MapFeature * narrowFeature(){return 0;}
    virtual short getObjectSubType(){return 0;}
    virtual short getObjectType(){return 0;}
    /**
     * Subclass can cast pointers with this by overriddding this 
     * like
     *
     * virtual MapObject *getFeatureType(long type){
     *     
     * if (type==MAPRELATIVEPOINT_TYPE)return this;
     *
     * return MapRelativeFeature::getObjectType(type);
     *}
     *
     * Then below the class definiation in the header you can do:
     *
     *static inline MapRelativePoint *getMapRelativePoint(MapBank *b,long key)
     *{
     *
     * if(b){
     *
     * MapObject *mo=b->getMapObject(key);
     * 
     *  if (mo) 
     *    return (MapRelativePoint*)mo->getFeatureType(MAPRELATIVEPOINT_TYPE);
     * 
     * }
     *
     * return 0;
     * }  
     *
     *
     */
    virtual MapObject * getObjectType(short type){
      if (type==getObjectType())return this;
      return 0;
    }
    virtual MapObject * getStateNodeType(short ){return 0;}
    virtual MapObject * getEnergyNodeType(short){return 0;}
    virtual MapObject * getFeatureType(short type){
      MapObject *mo=getObjectType(FEATURE_TYPE);
      if (!mo) return 0;
      if (mo->getObjectSubType()==type)
	return mo;
      return 0;
    }
    virtual MapObject * getLocalMapType(short ){return 0;}
    virtual MapObject * getFrameType(long type){
      MapObject *mo=getObjectType(FRAME_TYPE);
      if (!mo) return 0;
      if (mo->getObjectSubType()==type)
	return mo;
      return 0;
    }
    virtual MapPoint * narrowPoint(){return 0;}
    virtual void print(){
      std::cerr<<"Map Object "<<Key<<" "<<ID<<" \n";
    }

    /**
     * This packs all the informatiuon on this object into a GenericData
     * GenericType is the ObjectType 
     * The ObjectSubType, BankID, Key, and ID are stored in ShortData row 0
     * The subclass fills Data and extends ShortData
     *
     * @param gd Will be filled with this object's data
     */
    virtual void get(GenericData &gd);

    /**
     * This will set this object to the values stored in a GenericData.
     * It also associates the stored BankID/Key with this object's key
     * @param gd the data to unpack
     *
     * @return 1 if failed else 0, MAP_OBJECT_INVALID if gd is not a
     * right format
     */
    virtual int set(GenericData &gd);
  };
  /**
   * The base class of a factory object for MapObjects.
   *
   */
  class MapObjectMaker
  {
  public:
    /** 
     * Pointer to bank that manages all objects in the map and makes
     * it possible to delete features from application programs
     * without letting the SLAM program know
     */
    MapBank *Bank;
    /** The type of objects this makess */
    short ObjectType;
    /** The Subtype of objects this makess */
    short ObjectSubType;
    
 
    /** The subclass should set ObjectType and ObjectSubtype*/
    MapObjectMaker(MapBank *b){
      Bank=b;
      ObjectType=0;
      ObjectSubType=0;
    }
    virtual ~MapObjectMaker(){}
    /**
     * This looks for an object already on the Bank with ds's bankid
     * and key combination if theres is one but its the wrong type
     * then this returns 0. If the type stored in ds is not the type
     * of this maker then return 0.  Otherwise the object is setto ds.
     * If the object is not on the bank then a new one is made and
     * associated on the bank.
     *
     * @param ds the data is restored from this.
     * @return a MapObject set to ds or 0 if ds is not consistent.
     */
    MapObject * makeMapObject(GenericData &ds);
    /** 
     * The ObjectType is FRAME_TYPE, FEATURE_TYPE and so on @see mapDefines.hh 
     * @return the ObjectType produced by this class.
     */
    short getObjectType(){return ObjectType;}
    /** 
     * The ObjectSubType is MAPWALL_TYPE, POSETREE_TYPE, ... 
     * @see mapDefines.hh 
     * @return the ObjectSubType produced by this class.
     */
    short getObjectSubType(){return ObjectSubType;}
    
    /** This is the only virtual function*/
    virtual MapObject * makeMapObject()=0;
  };
   /**
     * Get the key and bank ID from an object stored in a GenericData.
     * @param gd holds the MapObject data.
     * @param bankid the BankId of the object is returned here.
     * @param key the key of the object is returned here.
     * @return 0 if ok else 1.
     */
  static inline int getBankInfo(GenericData &gd, short &bankid,long &key){
     if (gd.ShortData.Rows<1)return 1;
    if (gd.ShortData.Columns<6)return 1;
    key=gd.ShortData.getLong(0,2);
    bankid=gd.ShortData(0,1);
    return 0;
  }
  
   /**
     * Get the key and bank ID from an object stored in a DataSet.
     * @param gd holds the MapObject data.
     * @param objectType the type of object is returned here.
     * @param objectSubType the FeatureType etc. is returned here.
     * @return 0 if ok else 1.
     */
  static inline int getObjectInfo(GenericData &gd, short &objectType,
				  short &objectSubType){
    objectType=gd.GenericType;
    objectSubType=gd.ShortData(0,0);
    return 0;
  }
  
  
} // namespace Cure

#endif
