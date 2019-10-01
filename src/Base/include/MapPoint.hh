// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_MAPPOINT_H
#define CURE_MAPPOINT_H

#include "Point3D.hh"
#include "MapBank.hh"
#include "MapObject.hh"

namespace Cure{

  /**
   * A MapPoint is either a 3D or 2D feature coordinate coordinate.
   * These can be part of more than one feature.   
   */
  class  MapPoint: public MapObject, public Cure::Point3D
  {
  public:
    //The index is for indexing into a state vector but 
    //it can be used for something else.
    long  Index; 

    //The local_map is for indexing into a set of local_maps but 
    //it can be used for something else.
    long Local_map;

    MapObjectList Features;
  public:
    MapPoint(MapBank *b=0);
    MapPoint(const MapPoint& p);
    MapPoint(const Point3D& p, MapBank *b=0);
    ~MapPoint();
    void operator = (const Point3D& p){Vector3D::operator=(p);} 
    void operator = (const Point2D& p){Vector3D::operator=(p);} 
    void operator = (const Vector3D& p){Vector3D::operator=(p);} 
    void operator = (const double x[3]){Vector3D::operator=(x);}
    /**Cast from MapObject to MapPoint*/ 
    MapPoint* narrowPoint(){return this;}
    /** 
     * Add a feature to the feature list.  This list contains the
     * features that this point is attached to.  This does nothing to f.
     * @param f the MapFeature to attach the point to.
     */
    int addFeature(MapObject *f);
    /** 
     * Removes a feature from the feature list.  This list contains the
     * features that this point is attached to.  This does nothing to f.
     * @param f the MapFeature to remove.
     */
    int removeFeature(MapObject *f);
    /**
     * @return the number of features this point is attached to.
     */
    int numberOfFeatures(){return Features.count();}
    /**
     * Initializes the point location relative to some sensor.
     * @param sensorpose The transformation from the map to the sensor frame.
     * @param relx the xyz in the sensor frame.
     */
    void initializeFromRelative(Cure::Transformation3D &sensorpose,
				const double relx[3]){
      sensorpose.invTransform(relx,X);
    }

    /**
     * Initializes the point location based on a camera image.
     * @param t The transformation from the map to the camera frame.
     * @param pixels these are centered and right handed,
     * so normal image (u,v) becomes pixels[0] = u-centerpixel_u
     * and pixels[1]=centerpixel_v-v  
     * @param focallength the camera focallegth.
     * @param distanceguess an inital estimate of depth.
     */
    void initializeFromPixels(Cure::Transformation3D &t,
			      const double pixels[2],
			      double focallength,
			      double distanceguess);
    
    void get(GenericData &gd){
      MapObject::get(gd);
      gd.setShortDataSize(2,6);
      gd.ShortData.setLong(Index,1,0);
      gd.ShortData.setLong(Local_map,1,2);
      gd.Data.grow(1,3);
      gd.Data(0,0)=X[0];
      gd.Data(0,1)=X[1];
      gd.Data(0,2)=X[2];
    }
    int set(GenericData &gd){
      if (MapObject::set(gd))return 1;
      if (gd.ShortData.Rows<2)return 1;
      if (gd.Data.Columns<3)return 1;
      if (gd.Data.Rows<1)return 1;    
      Index=gd.ShortData.getLong(1,0);
      Local_map=gd.ShortData.getLong(1,2);
      X[0]=gd.Data(0,0);
      X[1]=gd.Data(0,1);
      X[2]=gd.Data(0,2);
      return 0;
    }
    void write(std::fstream &fs );
    void read(std::fstream &fs );
  };

} // namespace Cure

#endif
