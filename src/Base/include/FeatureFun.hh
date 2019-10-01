// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_FEATUREFUN_H
#define CURE_FEATUREFUN_H


#include "MapFeature.hh"
#include "MapFeatureList.hh"
namespace Cure {
namespace FeatureFcn {

  /**
   * This checks to see if a line passes thru a rectangle that is aligned 
   * with the coordinate axis. 
   *
   * @param start the xy of the start endpoint for the line
   * @param end the xy of the end endpoint for the line
   * @param bottomLeft the xy of the bottom left corner of the rectangle
   * @param topRight the xy of the topRight corner of the rectangle
   * @return 1 if line passes thru the rectangle. 
   */
  int lineInRectangle(double start[2],double end[2], double  bottomLeft[2],
		      double topRight[2]);
  /**
   * This checks to see if a point is in a rectangle that is aligned 
   * with the coordinate axis. 
   *
   * @param center the xy of the point.
   * @param bottomLeft the xy of the bottom left corner of the rectangle
   * @param topRight the xy of the topRight corner of the rectangle
   * @return 1 if center is in the rectangle. 
   */
  int pointInRectangle(double center[2], double  bottomLeft[2],
		       double topRight[2]);

  /**
   * This looks thru a list of features to find those inside a square
   * with a given center point and edges a given distance away, aligned
   * with the xy axis.
   * So the square has sides 2*distance and center is in the middle of it.
   */
  int near(MapFeatureList& near,MapFeatureList & features,
	   const double distance,const double center[2]);
  /**
   * Turns a list of longs into a list of valid features.
   */
  void listFeatures( MapFeatureList & features,MapBank &b, LongList &llist);
  /**
   * Turns a  list of valid features into a list of longs.
   */
  void listLongs(MapFeatureList& features, LongList &llist);

}; // namespace FeatureFcn

}; // namespace Cure

#endif
