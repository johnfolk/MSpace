// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_MAPWALL_H
#define CURE_MAPWALL_H

#include <string.h>  // memcpy
#include <fstream>  
#include "MapFeature.hh"
#include "FeatureFun.hh"
#include "MatrixStuff.hh"
#include "Point2DCloud.hh"
#include "Line2D.hh"

namespace Cure {

  /**
   * Class that manages features of type wall in the map
   *
   * @author Patric Jensfelt
   * @see
   */
  class MapWall: public MapFeature
  {
  public:
    /** narrow() will put this in here*/
    MapWall **CastPtr;
    MapPoint * AllocatedPoints[2];
    long AllocatedIndex[4];
    double Tangent[2];
    /* 
     * This says that the Start endpoint is measured.
     */
    int StartLimit;
    /* 
     * This says that the End endpoint is measured.
     */
    int EndLimit;


    /**
     * This is the min length for starting the P-Space.
     */
    double LengthThreshold;
    /**
     * This is the min scan points for starting the P-Space.
     */
    int CountThreshold;
    /**
     * This is the max variance in rho for starting the P-Space.
     */
    double VarRhoThreshold;
    /**
     * This is how close the points in the cloud must be to be included in wall. 
     */
    double TightnessValue;
    /**
     * If the squared distance from the Measured end point and the 
     * Current endpoint is less then this the end is limited and 
     * P-dim increased.
     */
    double EndThreshold;


    /**
     * The distance between the end points, updated every time that
     * calcTangent is called
     */
    double Length;

    Cure::Point2DCloud Cloud;
 
  protected:
    /**
     * This is the length of the line at the time that the endpoint
     * was merged.  It is used to scale the M dim for the case of a  single
     * shared point.  Thus the scale of the M-space continues to change with 
     * length if the unshared endpoint move tangent to the line.
     * This is important for angle measurements made after the merge.  
     * Thus there is a leverage of these measurements based on the ratio
     * of the current length to Ls.  
     *
     * This is a compromise as we throw away some information on the normal 
     * component of the unshared point in order to be able to make the line
     * longer.  We see that if we do not let the covariance of the point 
     * increase as the line gets longer we will end up with a overly precise
     * estimate of the angle of the wall.  By scaling with length we slightly
     * underestimate the precision on the position (rho) of the line, while 
     * correctly estimating the angle.
     *
     * The estimate is exact if the line never changes length after merging 
     * the endpoint.
     *
     * An alternative would be to not allow the line to extend the unshared 
     * point tangent to the line as new dense info is added.  We would then 
     * need to start a new line for that information.  We could later merge 
     * the two lines for an exact treatment.  This is easier.
     *
     */
    double Ls;

    int InfoInitialized;
    int CanLimit;
  public:
    MapWall();
    MapWall(MapWall * wp,MapBank *b=0);
    /** This is used to cast this is cound in _CastPtr after call,*/
    virtual void  narrow(){*CastPtr=this;}
    virtual MapFeature * copy();

    void setStart(MapPoint *val);
    void setEnd(MapPoint *val);  
 
    int inside(double bottomleft[2], double topright[2]){
      return FeatureFcn::lineInRectangle(Points[0]->X, 
                                         Points[1]->X,bottomleft, topright);
    }
    virtual int  extend();
    virtual unsigned short getMeasurementType(unsigned short type);
    virtual int  getMtypeProjection(Cure::Matrix &project, 
				     unsigned short mtype);

    /**
     * Adds the new collected Info to the cumulation and updates the 
     * LowInfo dimensions.
     * Have not implemented type!=0 yet
     * @param v A cloud of (distance,x,y) (N x 3) used to set endpoints.
     * @param covAdjust for EKF do cov -> covAdjust*cov*covAdjust.
     * @param type 1 for start measured, 2 for end measured, 3 both, 0 niether.  
     * @return 0 if ok 
     */
    virtual int  addInfo(const Cure::Matrix & v,const int type=0){
      Cure::Transformation3D I;
      return addInfo(v,I,type);
    }
    virtual int  addInfo(const Cure::Matrix & v,
			 Cure::Transformation3D & map2info,
			 const int type=0);
    void recenter();

    bool getC(MapWall *mw,Cure::Matrix &c)const;
    bool testMatch(MapWall *mw, double tolerance);
    /**
     * @return number between 0 and 1 that gives a percentage 
     *         overlap of the walls.
     */
    double judgeMatch(MapWall *mw, double tolerance);
    virtual void getB(Cure::Matrix &b)const;
    /**
     * The B matrix has rows that are orthogonal vectors with magnitude
     * 1/(scale factor).  That scale factor will either be 1 or a length 
     * in some direction in the Xf space.  For instance the length of the 
     * wall.
     * This method returns the matrix needed to calculate the current scale
     * factor from the current feature coordinates.  If the scale factors are
     * 1 then that row is left off of scales.  To find the coorespondence
     * between rows of scales and rows of B, examine the return value.
     * It is a list of binary flags that are 1 if the corresponding B row is
     * scaled and 0 otherwise.  Thus if zero is returned no scaling is done.
     * if 1 is returned the first row of b is scaled. If 4 is returned then
     * the third row of b is scaled.
     *
     * @param scales formed into the Matrix that multiplies Xf to get the 
     *               scale vector for the rows of the B matrix.  It will have 
     *               rows for each scaled P dim.
     * @return binary flags to scaled rows of B.
     */
    virtual unsigned short getScales(Cure::Matrix &scales);
    
    /**
     * @return the flags for M-space rows that have frames attached to the 
     * feature. 
     */
    virtual unsigned short getFramed();

    /**
     * This is to set the endpoints only and init the M-space for 
     * known walls.
     * @param startx the xyz of the start point.
     * @param endx the xyz of the end point.
     * @param forceextend bit flags to indicate the M-space, 1 for rho
     *        and gamma, 2 for start, and 4 for end.
     */
    void initailize(double startx[3], double endx[3], 
		    unsigned short forceextend=0);
    void initailizeFromLine(Cure::Transformation3D &t,Cure::Line2D & ln); 
    /**
     * Merges all info from mf to this.  
     * If merging a point it replaces mf's point and deletes 
     * the point that is not not left attached to 
     * any feature.  It then changes mf->Index to the index of the common
     * point.  It also sets of SharedIndex for this and mf so that 
     * Points[SharedIndex[k]] is the point that P space dim k is: 
     * P(k,0)=Points[SharedIndex[k]].X(k-n,0) where n=2*SharedIndex[k].  
     */
    int merge(MapWall *mf, unsigned short type=1);
    void forceExtend(int typ);
    virtual void print(int level);
    virtual void print(){print(0);}
    /**
     * Calculates the unit tangent vector (Tangent) based on the current 
     * endpoints.
     */
    void calcTangent();
    /**
     * Calculates the unit normal vector  based on the current 
     * endpoints.  Pointing from the origin perpendicular to the line 
     */
    void calcNormal(double x[2]);
    void initializeToCloud(Cure::Point2DCloud & c,Cure::Transformation3D & map2info);
    
    /**
     * Depreciated use int read(version...)
     * This initializes a MapWall from the information in a text file.
     * @param fs the file to read from. It should be opened with std::ios::in
     * @param readKey set this to false if key has been read allready.
     * @return 0 if fail else this
    */
    MapWall * read(std::fstream &fs, bool readKey=true, int version=1);
    /**
     * This initializes a MapFeture from the information in a text file.
      * @param version 1 only for now.
     * @param fs the file to read from. It should be opened with std::ios::in
     * @param readKey set this to false if key has been read allready.
     * @return 1 if fail else 0
     */
    virtual int read(int version, std::fstream &fs, bool readKey=true);

    /**
     * Write the MapWall information to a text file.
     * @param fs the file to write to. It should be opened with std::ios::out
     */
    void write(std::fstream &fs );
    /**
     * Stores the feature in a Data object
     * @param gd the Data to store this feature in.
     */
    virtual void get(GenericData &gd);
    /**
     * Restores the feature from a GenericData object
     * @param gd the Data to restore this feature from.
     * @return 1 if fails else 0;
     */
    virtual int set(GenericData &gd);
  protected:
    bool dimCheck(){
      if (Bdual.Columns==2)return true;
      return false;
    }
    //    void lengthAdjust();
  };

} // namespace Cure

#endif
