// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "MapHelper.hh"
#include "CureDebug.hh"
#include <iostream>

using namespace Cure;
bool MapHelper::getCommon(MapFeature  *mf1,MapFeature  *mf2,
			  Cure::Matrix &c1, Cure::Matrix &c2)
{
  if (!getC(mf1,mf2,c1))return false;
  Matrix scales;
  c2=c1;
  unsigned short temp=mf1->getScales(scales);
  if (temp)
    {
      Matrix xi,xj;
      mf1->getX(xi);
      mf2->getX(xj);
      Matrix ti, tj;
      ti.multiply_(scales,xi);
      tj.multiply_(scales,xj);
      int row=0;
      for (int j=0;temp>0;j++)
	{
	  if (temp&0x1)
	    {
	      double d=(ti(row,0)/tj(row,0));
	      for (int k=0;k<c2.Rows; k++)c2(k,j)*=d;
	      row++;
		  }  
	  temp=(temp>>1);
	}
    }
  return true;     
}
void MapHelper::loosenMatching(double factor)
{
  RoughSearchRange*=factor;
  NewFeatureRange*=factor;
  MatchThreshold*=factor;
}
void MapHelper::match(Match *matches,int n,
		      double distance,Pose3D &psensor, LongList &mapFeats)
		      
{
  Visable.clean();
  // Delete features that have not been observed for long and have not
  // accumulated enough info yet to have any P-dimensions.
  deletePartialFeatures(distance,&mapFeats);
  // Find map features that are near the robot (stored in Near list)
  makeNear(psensor.Position.X,&mapFeats);
  CureDO(50) {
    int nn = 0;
    for (MapFeatureList *mlist=&Near; mlist->Next; mlist=mlist->Next) {
      nn++;
    }
    CureCERR(0) << "Found " << nn << " features near\n";
  }
  
  // Clear prev list of PosedFeatures and add new ones based on Near list
  makePosedFeatures();


  CureDO(50) {
    int np = 0;
    for (PosedFeatureList *mlist=&PosedFeatures; mlist->Next; 
         mlist=mlist->Next) {
      np++;
    }
    CureCERR(0) << "Made " << np << " posed features\n";
  }
  // Clear list of  visible features from last iteration  and look for
  // the currently visible features (from PosedFeatures list).
  makeVisable(psensor,SensorType);
  CureDO(50) {
    int nv = 0;
    for (PosedFeatureList *mlist=&Visable; mlist->Next; mlist=mlist->Next) {
      nv++;
    }
    CureCERR(0) << "Found " << nv << " visible features\n";
  }

  for (int i=0; i<n;i++) {
    if (!(addDistance(matches[i],distance))) {
      if(matches[i].Measure->Key==0){
	// List of features in sensor coordinates that are potential
	// matches to this measurement.
	PosedFeatureList potmatches;
      
	// Look for rough matches between a measurement,
	// matches[i].Measure, and the Visable features. Basically
	// checking which feature in the map are predicted to be close
	// enough to the measurement.
	matchFeatures(potmatches,matches[i],&Visable);
	
      
	// Make finer match between the potential matches and the
	// measurement. Map features that do not match the measurement
	// but are quite close will be marked with MatchDistance==-2
	// so that they are not added as new features.
	matchFeature(matches[i],&potmatches);
      
      }  else {
	double dist=0;
	long fkey=findTrackedKey(matches[i].Measure->Key);
	if (fkey==-1)matches[i].MatchDistance=-1;  
	// Not in the tracked table -> new feat
	else{
	  matches[i].MatchDistance=-2;
	  for (PosedFeatureList *mlist=&Visable; mlist->Next; 
	       mlist=mlist->Next)
	    {
	      PosedFeature *pf=mlist->Element;
	      if (pf)
		{
		  if (pf->FeatureKey==fkey){
		    dist=pf->fineMatch(*matches[i].Measure,matches[i].Metric);
		    if (dist<MatchThreshold)
		      {
			matches[i].MatchedFeature=pf->copy();
			matches[i].MatchDistance=dist;
		      }
		  }
		}
	      
	    }
	}
      }
      // A measurement that does not match anything and creates a
      // new map feature assuming that it is not too close to an
      // existing map feature (MatchDistance==-2)
      if (matches[i].MatchDistance==-1) {
        MapFeature *tempfeat=makeMapFeature(psensor,
                                            *matches[i].Measure,
                                            true);
        if (tempfeat) {
          mapFeats.add(tempfeat->Key);
	  matches[i].MatchedFeature=makePosedFeature(tempfeat);
          matches[i].MatchedFeature->transform(psensor);
          matches[i].MatchDistance=0;
        }
      }
    }
  }
}

void MapHelper::simpleMatch(Match *matches,int n,
			    Pose3D &psensor, double threshold)
{
  if (threshold==0)threshold=MatchThreshold;
  threshold=MatchThreshold;
  MatchThreshold=threshold;
  for (int i=0; i<n;i++) {
    if(matches[i].Measure->Key!=0)
      {
	double dist=0;
	long fkey=findTrackedKey(matches[i].Measure->Key);
	if (fkey==-1)matches[i].Measure->Key=0;  
	else{
	  matches[i].MatchDistance=-2;
	  for (PosedFeatureList *mlist=&PosedFeatures; mlist->Next; 
	       mlist=mlist->Next)
	    {
	      PosedFeature *pf=mlist->Element;
	      if (pf)
		{
		  if (pf->FeatureKey==fkey){
		    int r=pf->transform(psensor);
		    if (!r){
		      dist=pf->fineMatch(*matches[i].Measure,matches[i].Metric);
		      if (dist<MatchThreshold)
			{
			  matches[i].MatchedFeature=pf->copy();
			  matches[i].MatchDistance=dist;
			}
		    else matches[i].MatchedFeature=0; 
		    }
		  }
	      }
	      
	    }
	}
      }
    if(matches[i].Measure->Key==0){
      PosedFeatureList potmatches;      
      Measurement *m=matches[i].Measure;
      double rect[4],rect4[8];
      rect[0]=m->BoundingBox(0,0)-RoughSearchRange;
      rect[1]=m->BoundingBox(0,1)-RoughSearchRange;
      rect[2]=m->BoundingBox(1,0)+RoughSearchRange;
      rect[3]=m->BoundingBox(1,1)+RoughSearchRange;
      if (rect[0]>rect[2]){
	double d=rect[0];
	rect[0]=rect[2];
	rect[2]=d;
      }
      if (rect[1]>rect[3]){
	double d=rect[1];
	rect[1]=rect[3];
	rect[3]=d;
      }
      psensor.invTransform2D(rect,rect4);
      psensor.invTransform2D(rect+2,rect4+2);
      double d=rect[2];
      rect[2]=rect[0];
      rect[0]=d;
      psensor.invTransform2D(rect,rect4+4);
      psensor.invTransform2D(rect+2,rect4+6);
      rect[0]=rect4[0];
      if (rect[0]>rect4[2])rect[0]=rect4[2];
      if (rect[0]>rect4[4])rect[0]=rect4[4];
      if (rect[0]>rect4[6])rect[0]=rect4[6];
      rect[2]=rect4[0];
      if (rect[2]<rect4[2])rect[2]=rect4[2];
      if (rect[2]<rect4[4])rect[2]=rect4[4];
      if (rect[2]<rect4[6])rect[2]=rect4[6];
      rect[1]=rect4[1];
      if (rect[1]>rect4[3])rect[1]=rect4[3];
      if (rect[1]>rect4[5])rect[1]=rect4[5];
      if (rect[1]>rect4[7])rect[1]=rect4[7];
      rect[3]=rect4[1];
      if (rect[3]<rect4[3])rect[3]=rect4[3];
      if (rect[3]<rect4[5])rect[3]=rect4[5];
      if (rect[3]<rect4[7])rect[3]=rect4[7];
      Visable.clean();
      for (PosedFeatureList *mlist=&PosedFeatures; mlist->Next;
	   mlist=mlist->Next )
	{
	  PosedFeature *pf=mlist->Element;
	  pf->setXo();
	  if (pf)
	    if (pf->inRectangle(rect))
	      {
		int r=pf->transform(psensor);
		if (!r)
		  {
		    if (SensorType==CAMERA_TYPE){
		      if (pf->inRectangle(ImagePlane)) 
			Visable.add(pf);
		    } else {
		      Visable.add(pf);
		    }
		  }
	      }
	}
      // matchFeatures(potmatches,matches[i],&Visable);
      matchFeature(matches[i],&Visable);//&potmatches);
    }
  }
}
void MapHelper::makeVisable(Pose3D& pose, int sensortype)
{
  Visable.clean();
  for (PosedFeatureList *mlist=&PosedFeatures; mlist->Next;) {
    PosedFeature *pf=mlist->Element;
    int r=pf->transform(pose);
    if (r) {
      if (r&MAP_OBJECT_INVALID) {
        PosedFeatures.removeFeature(pf);
        delete pf;
      } else {
        mlist=mlist->Next;
      }
    } else {
      if ((m_UseImagePlane)||(sensortype==CAMERA_TYPE)){
	if (pf->inRectangle(ImagePlane)) {
          Visable.add(pf);
        }
      } else {
	Visable.add(pf);
      }
      mlist=mlist->Next;
    }
  }
  if (sensortype==SCAN_TYPE) {
    double scan[360];
    makeFeatureScan(Visable,scan);
  }
}

void MapHelper::makeFeatureScan(PosedFeatureList &visFeatures,
				 double scan[360])
{
  ScanFragment sf;    
  PosedFeature *featurescan[360];
  for (int i=0;i<360;i++)
    {
      featurescan[i]=0;
      scan[i]=MAXSCAN;
    }
  for (PosedFeatureList *mlist=&Visable; mlist->Next; 
       mlist=mlist->Next)
    {
      PosedFeature *pf=mlist->Element;
      int r=pf->scan(sf);
      if (!(r&NOT_VISABLE))
	{
	  for (ScanFragment *s=&sf;s->Next; s=s->Next)
	    {
	      while(s->Index<0)s->Index+=360;
	      while(s->Index>359)s->Index-=360;
	      if (scan[s->Index]>s->Range)
		{
		  scan[s->Index]=s->Range;
		  featurescan[s->Index]=pf;
		}
	    }
	}
    }
  long lastkey=-1;
  visFeatures.clean();
  for (int i=0;i<360;i++)
    {
      if (featurescan[i])
	{
	  long key=featurescan[i]->FeatureKey;
	  if (key!=lastkey)
	    {
	      lastkey=key;
	      if (!(visFeatures.keyGet(key)))
		{
		  visFeatures.add(featurescan[i]);
		}
	    }
	}
    }
}
void MapHelper::matchFeatures(PosedFeatureList & results, Match &mat,
				PosedFeatureList * searchList)
{
  Measurement *m=mat.Measure;
  results.clean();
  if (searchList==0)searchList=&Visable;

  // Create a rectangle for initial maching
  double rect[4];
  rect[0]=m->BoundingBox(0,0)-RoughSearchRange;
  rect[1]=m->BoundingBox(0,1)-RoughSearchRange;
  rect[2]=m->BoundingBox(1,0)+RoughSearchRange;
  rect[3]=m->BoundingBox(1,1)+RoughSearchRange;

  // Make sure that the corners in the initial match rectangle are
  // defined correctly, i.e. that the bounding box is correct
  if (rect[0]>rect[2])
    {
      double d=rect[0];
      rect[0]=rect[2];
      rect[2]=d;
    }
  if (rect[1]>rect[3])
    {
      double d=rect[1];
      rect[1]=rect[3];
      rect[3]=d;
    }
  for (PosedFeatureList *mlist=searchList; mlist->Next; 
       mlist=mlist->Next)
    {
      PosedFeature *pf=mlist->Element;
      if (pf)
	if (pf->inRectangle(rect))
	  {
	    if (!(pf->roughMatch(m->Z.Element,mat.Thresholds.Element)))
	      {
		results.add(pf);
	      }
	  }
    }
}
void MapHelper::matchFeature(Match & mat,PosedFeatureList * searchlist)
{
  Measurement &m=*mat.Measure;
  if (searchlist==0)searchlist=&Visable;
  double dist=10*NewFeatureRange;
  PosedFeature *mpf=0;
  if (mat.MatchedFeature)
    {
      delete mat.MatchedFeature;
      mat.MatchedFeature=0;
    }

  // Lopp through the features and look for the "nearest" one
  if(mat.Measure->Key==0){
    mat.MatchDistance=-1;
    for (PosedFeatureList *mlist=searchlist; mlist->Next; 
	 mlist=mlist->Next)
      {
	PosedFeature *pf=mlist->Element;
	if (pf)
	  {
	    double d=pf->fineMatch(m,mat.Metric);
	    if ((d<dist)&&(d!=-1))
	      {
		mpf=pf;
		dist=d;
	      }
	  }
	
      }
    // A match that is too far away to be a match but too close to an
    // existing feature will result in MatchDistance=-2 and can thus be
    // discarded. We do not want features that are to close to each other.
    //
    // Here we assume that the NewFeatureRange is larger than
    // MatchThreshold, i..e if the dist it is closer than
    // NewFeatureRange it is also closer than MatchThreshold.
    
    if (dist<NewFeatureRange)
      {
	if ((mpf)&&(dist<MatchThreshold))
	  {
	    mat.MatchedFeature=mpf->copy();
	    mat.MatchDistance=dist;
	  }
	else 
	  {
	    mat.MatchDistance=-2;
	  }
      }
  }
  else {
    dist=0;
    long fkey=findTrackedKey(mat.Measure->Key);
    if (fkey==-1)mat.MatchDistance=-1;  // Not in the tracked table -> new feat
    else{
      mat.MatchDistance=-2;
      for (PosedFeatureList *mlist=searchlist; mlist->Next; 
	   mlist=mlist->Next)
	{
	  PosedFeature *pf=mlist->Element;
	  if (pf)
	    {
	      if (pf->FeatureKey==fkey){

		dist=pf->fineMatch(m,mat.Metric);
		if (dist<MatchThreshold)
		  {
		    mat.MatchedFeature=pf->copy();
		    mat.MatchDistance=dist;
		  }
	      }
	    }
	}
    }
  }
}

int MapHelper::writeTrackKeys(std::fstream &fs )
{
  if (TrackHashSize==0)return 0;
  if (!fs.good())return 1;
  fs<<FeatureType<<" "<<TrackHashSize<<" "<<"\n";
  for (int i=0; i<TrackHashSize;i++)
    {
      for (LongPairList *lplist=&TrackKeys[i]; lplist->Next; lplist=lplist->Next)
	fs<<lplist->Element<<" "<<lplist->UnsignedElement<<"\n";
    }
  fs<<0<<" "<<0<<"\n";
  return 0;
}
void MapHelper::readTrackKeys(std::fstream &fs )
{
  unsigned short ths=0;
  fs>>ths;
  setTrackHashSize(ths);
  unsigned long tkey=1;
  long fkey=1;
  while(tkey!=0)
    {
      fs>>fkey>>tkey;
      addTrackedKey(tkey,fkey);
    }
}
int MapHelper::writeFeature(long key,std::fstream &fs )
{
  if (!fs.good())return 1;
  MapFeature *mf=getMapFeature(key);
  if (!mf)return -1;
  mf->write(fs);
  return 0;
}

MapFeature * MapHelper::readFeature(std::fstream &fs, 
				    int version,int *type,bool readtype)
			       
{
  if (version!=1){
    std::cerr<<"MapHelper::readFeature Vesion not supported\n";
    return 0;
  }
  if (!fs.good())return 0;
  if (readtype){ 
    if (!(fs>>*type))return 0;
    if (*type!=FeatureType)return 0;
  }
  long key=0;
  if (!(fs>>key))return 0;
   MapObject *m=Bank->getMapObject(key);
   if (m){
    std::cerr<<"MapObject with this key existes "<<key<<"\n";
  }
  MapFeature *mw=makeMapFeature();
   if (!mw){
    std::cerr<<"read Could not make mapFeature";
  }
  mw->Bank=Bank;
  Bank->add(mw,key);
  if (key!=mw->Key){
    std::cerr<<"wanted Feature key "<<key<<" got "<<mw->Key;
  }
   int ret= mw->read(version,fs, false);
   if (!ret)return mw;
  delete mw;
  return 0;
}
