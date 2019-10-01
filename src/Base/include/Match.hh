// = AUTHOR(S)
//    John Folkesson
//    
//    August, 2004
//
//    Copyright (c) 2004 John Folkesson
//    

#ifndef CURE_MATCHES_HH
#define CURE_MATCHES_HH


#include <iostream>
#include <stdlib.h> 
#include "PosedFeature.hh"
#include "Measurement.hh"
#include "Pose3D.hh"

namespace Cure{

  /**
   * Class that manages matches between measurements and the map in
   * the form of a posed feature.
   *
   * @author John Folkesson
   */
  class Match
  {
  public:
    /** 
     * Pointer to the measurement with high level measurement data Z,
     * the dense information W, etc
     */
    Cure::Measurement *Measure;

    /** @todo comment */
    Cure::Matrix Metric;

    /** @todo comment */
    Cure::Matrix Thresholds;

    /** The feature that this measurement ha sbeen matched to */
    PosedFeature *MatchedFeature;

    /** @todo comment */
    double PathDistance;

    /** @todo comment */
    double Weight;

    /** @todo comment */
    double MatchDistance;
  protected:
  public:
    Match(){
      Measure=0;
      MatchedFeature=0;
      MatchDistance=-1;
      Weight=0;
      PathDistance=0;
    }
    ~Match(){
      clear();
    }
    Match(Match &m)
    {
      (*this)=m;
    }

    void operator = (Match &m)
    {
      MatchedFeature=m.MatchedFeature->copy();
      MatchDistance=m.MatchDistance;
      Measure=m.Measure;
      Metric=m.Metric;
      Thresholds=m.Thresholds;
      Weight=m.Weight;
      PathDistance=m.PathDistance;
    } 
    void clear(){
      if (MatchedFeature)delete MatchedFeature;
      MatchedFeature=0;
      MatchDistance=-1;
      Weight=0;
      PathDistance=0;
    }
    void clean(){
      MatchedFeature=0;
      MatchDistance=-1;
      Weight=0;
      PathDistance=0;
    }
    void print(int level=0);
  };
}
#endif
