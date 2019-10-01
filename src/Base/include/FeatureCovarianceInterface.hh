//
// = FILENAME
//    FeatureCovarianceInterface.hh
//
// = FUNCTION
//
// = AUTHOR(S)
//    Patric Jensfelt
//
// = COPYRIGHT
//    Copyright (c) 2006 Patric Jensfelt
//
/*----------------------------------------------------------------------*/

#ifndef Cure_FeatureCovarianceInterface_hh
#define Cure_FeatureCovarianceInterface_hh

#include "Matrix.hh"
#include "MapFeature.hh"

namespace Cure {

class FeatureCovarianceInterface {
public:

  /**
   * Destructor that is needed because there are virtual functions that
   * want a virtual destructor.
   */
  virtual ~FeatureCovarianceInterface() {}

  /**
   * Get the covariance of a feature's M-Space.
   * @param returns witht the Covariance of just this feature's
   *        M-Space.
   * @param mf the MapFeature to find the covaraince of
   * @return M-Space dimension of mf 
   */
 virtual int getFeatureCov(Matrix &fCov, MapFeature *mf) = 0;

};

};

#endif
