// = AUTHOR(S)
//    John Folkesson
//    
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_MAPHLINE_HH
#define CURE_MAPHLINE_HH


#include "MapFeature.hh"
#include "FeatureFun.hh"
#include "MatrixStuff.hh"
#include "LinkedArray.hh"

namespace Cure{

/**
 * Class used for horizontal lines in the map
 *
 * @author John Folkesson
 */
class MapHLine: public MapFeature
  {
  public:
    MapHLine **CastPtr;
    MapPoint * AllocatedPoints[2];
    long AllocatedIndex[6];
    double Rho[3];
    double Tangent[3];
    double InfoTangent[3];
    double Length;

    int VCount;
    Cure::LinkedArray Vectors;
    double TotalWgt;
    double SqWidth;
    /**
     * The 'weights' are summed for each pair of planes that 
     * contribute to the weighted triangulation.  These weights are
     * the product of the weight supplied with each pair of bearing vectors
     * and the sin^2(angle between the planes).
     * This sum is then compared to TanThreshold and if it is > the 
     * P dim are expanded.
     */
    double TanThreshold;
    double RhoThreshold;
    double TrackThreshold;
  protected:
 
  public:
    /**
     *  @param wp The templete feature, CastPtr, WeightThreshold and
     *         DistanceThreshold  are copied from the templete.
     */
    MapHLine();
    MapHLine(MapHLine * wp,MapBank *b=0);
    virtual void narrow(){*CastPtr=this;}
    void setStart(MapPoint *val);
    void setEnd(MapPoint *val);  
    /**
     * Write the MapPoint information to a text file.
     * @param fs the file to write to. It should be opened with std::ios::out
     */
    void write(std::fstream &fs ){
      if (Bdual.Columns<1)return;
      MapFeature::write(fs);
      fs<<LastBearing[0]<<" "<<LastBearing[1]<<" "<<LastBearing[2]<<"\n";
    }
      
    /**
     * This initializes a MapFeture from the information in a text file.
     * @param version 1 only for now.
     * @param fs the file to read from. It should be opened with std::ios::in
     * @param readKey set this to false if key has been read allready.
     * @return 1 if fail else 0
    */
    virtual int read(int version, std::fstream &fs, bool readKey=true){
      if (version!=1)return 1;
      int ret=MapFeature::read(1,fs,readKey);
      if (ret)return ret;
      fs>>LastBearing[0]>>LastBearing[1]>>LastBearing[2];
      recenter();
      return 0;
    }
    /**
     * @param testValue=z compenent of c=(s x e)/|s x e| in the map frame
     * @return 0 if this can be reliably matched to based on testValue
     */
    int testVisable(double *c){
      if (Bdual.Columns>1)return 0;
      if (Bdual.Columns>0)return 0;
      double testValue=1-(LastBearing[0]*c[0]+LastBearing[1]*c[1]+LastBearing[2]*c[2]);
      if (testValue>TrackThreshold)
	return NOT_VISABLE;
      return 0;
    }
    int merge(MapHLine *mf,unsigned short typ=1);
    int inside(double bottomleft[2], double topright[2]){
      return FeatureFcn::lineInRectangle(Points[0]->X, 
                                         Points[1]->X,bottomleft, topright);
    }
    virtual unsigned short getFramed(){
      if (Bdual.Columns==1)return 1;
      if (Bdual.Columns==3)return 7;
      return 0;
    }
    virtual int extend();
    virtual unsigned short getMeasurementType(unsigned short type);
    /**
     * Adds the new collected Info to the cumulation and updates the 
     * LowInfo dimensions.  This tracks the line using the latest measurment
     * if the line is not fully initialized.  This means the LowInfo update is 
     * done so that the last observation agrees with the prediction as far as 
     * the LowInfo can align it.
     * If distance is<LastDistance the new info will replace
     * the previous info at this distance.
     *
     * @param v Row Matrix(1,11)(distance, w, x, y, z, sx, sy, sz, ex, ey, ez) 
     * weight to give this measurment, xyz are the camera position, sx sy sz
     * is a vector from the camera towards the startpoint, (length ignored)... 
     * @param covAdjust for EKF do cov -> covAdjust*cov*covAdjust.
     * @param type 0.  
     * @return 0 
     */
    virtual int  addInfo(const Cure::Matrix & v,const int type=0){
      Cure::Transformation3D I;
      return addInfo(v,I,type);
    }
    virtual int  addInfo(const Cure::Matrix & v, 
			 Cure::Transformation3D & map2info,
			 const int type=0);
    void recenter();
    virtual void getB(Cure::Matrix &b)const;
    /**
     * The B matrix has rows that are orthogonal vectors with magnitude
     * 1/(scale factor).  That scale factor will either be 1 or a length 
     * in some direction in the Xf space.  For instance the length of the 
     * line.
     * This method returns the matrix needed to calculate the current scale
     * factor from the current feature coordinates.  If the scale factors are
     * 1 then that row is left off of scales.  To find the coorespondence
     * between rows of scales and rows of B, examine the return value.
     * It is a list of binary flags that are 1 if the corresponding B row is
     * scaled and 0 otherwise.  Thus if zero is returned no scaling is done.
     * if 1 is returned the first row of b is scaled.
     * 
     *
     * @param scales formed into the Matrix that multiplies Xf to get the 
     *               scale vector for the rows of the B matrix.  It will have 
     *               rows for each scaled P dim.
     * @return binary flags to scaled rows of B.
     */
    virtual unsigned short getScales(Cure::Matrix &scales);
    bool testMatch(MapHLine *mw, double tolerance);
    bool getC(MapHLine *mw,Cure::Matrix &c)const;
    /**
     * 
     * @param distanceGuess distance relative to camera.
     */
    void initializeFromPixels(Cure::Transformation3D &t,
			      const double startpixels[2],
			      const double endpixels[2],
			      double focallength,
			      double distanceGuess);

    virtual void print(int level=0);
    /**
     * Calculates the unit tangent vector (Tangent) based on the current 
     * endpoints.
     */
    void calcTangent();
    /**
     * Calculates the vector, Rho, based at the origin, perpendicular 
     * to Tangent and intersecting the infinite extention of the line. 
     * This is based on the current endpoints and Tangent.
     * This then sets Rho to map frame.
     */
    void calcRho();
    /**
     * Restores the feature from a GenericData object.  If the MapPoints
     * of gd already exists on the Bank then the data is simply copied to
     * the existing object.   By already exists we mean that a MapPoint
     * is associated with gd's BankID, Key combination.  In other cases 
     * the Points on this feature are removed and the new ones added.  If the
     * removed points have no features they are deleted.
     * No Initialization info is restored.
     * @param gd the Data to restore this feature from.
     *
     * @return 0 if ok,1 if ds is formated wrong, MAP_OBJECT_INVALID
     * if the point's keys are associated with the wrong type MapObject.
     * 
     */
    virtual int set(GenericData &gd);

  protected:
    virtual bool dimCheck(){
      if (Bdual.Columns==3)return true;
      return false;
    }
    /**
     * Converts from raw info to cooked.
     * Tangent and Rho must be in the info frame before calling this.
     * @param v the new dense info matrix
     * @param the array to store the made dense info 
     */
    void makeInfo(const Matrix &v, double m[28]);
    /**
     * Removes the old information from Vectors. 
     */
    void  prune(double minDistance);

    //  virtual void recenter();
    /**
     * This adds the Rho/Tangent dependent parts of the derived 
     * dense information,
     * such as dot products with Tangent, and the best LSQ fit distances,
     * (r and q).
     *
     * Be sure that Rho and Tangent are in info frame when calling this.
     *
     * @param v this fills in v[5,16..19,24..27] and uses the rest. 
     * v=(0:distance, 1:weight, 2..4:x=camera position, 5:Tangent_dot_x, 
     * 6..8:ĉ=(ŝ x ê)/|(ŝ x ê)|, 9..11:(t_i weighted implied tangent),
     * 12..14:ŝ, 15:ŝ_dot_x, 16:Tangent_dot_ŝ,
     * 17:Rho_dot_ŝ, 18:r_s, 19:q_s, 20..22:ê, 23:ê_dot_x, 24:Tangent_dot_ê,
     * 25:Rho_dot_ê, 26:r_e, 27:q_e)
     */ 
    void fit_it(double v[28]);

    int getEllipsiod(Cure::Matrix &cov, double wsums[4]);
    // void cloud_it(Cure::Matrix & cloud);
    /**
     * This fits all observations to a line then takes the resulting cloud
     * and recalcualtes the line parameters: Tangent and Rho. it then returns
     * the covariance matrix as three eigen vectors and eigen values, (the
     * first, largest, one is the tangent direction.  
     * InfoTangent must be current when calling this.
     * Will set InfoTangent and Rho (in info frame).
     * @return 0 if ok else 1
     */
    int eigen_it(double ev[9], double lambda[3],  double wsums[4]);
    /**
     * Endpoints are slide tangent to the line as far as is consistant 
     * with all the meassurements in Vectors.  This then changes Length.
     */
    void setEnds(Cure::Transformation3D & map2info);
    /**
     * Sets the endpoints to be consistant to the measurment m, normally 
     * the latest measurment.  If Tangent is in P space then the line
     * is only moved such that the P space is unchanged. ie. normal to
     * the line. 
     */
    void trackLine(Cure::Matrix &m,Cure::Transformation3D & map2info);
    /**
     * Uses the current map frame endpoints to estimate the Tangent Rho vector
     * where the equation of the line is x=Rho+t*Tangent for scalar t in 
     * some interval and X, Rho and Tangent 3D vectors. 
     * Ends up with both the Tangent and Rho trandfromed to the info frame.
     */
    void calcInfoTanRho(Cure::Transformation3D & map2info);
    /*
     * This estimates the Rho in the info frame based on dense info and
     * Tangent.
     * Tangent must be in infor frame when calling this
     */
    int estimateRho();
    void lengthAdjust();
    
  };

} // namespace Cure

#endif
