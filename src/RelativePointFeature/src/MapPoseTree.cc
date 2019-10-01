// = RCSID
//    $Id: MapPoseTree.cc ,v 1.1 2007/08/1 
//
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2007 John Folkesson
//    

#include "MapPoseTree.hh"
#include "FeatureDescriptors.hh"
#include "RelativePointFeature.hh"
using namespace Cure;

MapPoseTree::~MapPoseTree()
{
  if (Index)delete[]Index;
  Index=0;
  m_SizeIndex=0;
  MapPoseTree::clean();
  if (Index)delete[]Index;
  Index=0;
  m_SizeIndex=0;
}
void MapPoseTree::reallocate(){
  CalcState=0;
  int d=m_SizeIndex;
  int dm=dim();
  long *ind=Index;
  
  if (d!=dm){
    m_SizeIndex=dm;
    if (m_SizeIndex>0){
      Index=new long[m_SizeIndex];
      if ((ind)&&(d>0))memcpy(Index,ind,d*sizeof(long));
      for (int i=d;i<m_SizeIndex;i++)
	  Index[i]=0;
    }else Index=0;
    if (ind)delete []ind;
  }
  FeatureDescriptors **fd=m_Descriptors;
  int k=m_AllocatedDescriptors;
  m_AllocatedDescriptors=m_Branches.Rows;
  if (k<m_AllocatedDescriptors){
    if (m_AllocatedDescriptors>0)
      m_Descriptors=new FeatureDescriptors*[m_AllocatedDescriptors];
    else m_Descriptors=0;  
    for (int i=0;i<k;i++)
      m_Descriptors[i]=fd[i];
    for (short i=k;i<m_AllocatedDescriptors;i++){
      m_Descriptors[i]=new FeatureDescriptors(this, i);
    }
    if (fd)delete []fd;
  }else if (k>m_AllocatedDescriptors){
    if (m_AllocatedDescriptors>0)
      m_Descriptors=new FeatureDescriptors*[m_AllocatedDescriptors];
    else m_Descriptors=0;  
    for(int i=m_AllocatedDescriptors;i<k;i++){
      FeatureDescriptors *fds=fd[i];
      for (int j=0;j<fds->m_NumberDesc;j++){
	FeatureDescriptor *des=fds->m_Desc[j];
	if (des){
	  RelativePointFeature *mf=des->getRelativePointFeature();	   
	  if(mf)mf->remove(des);
	}
      }
      delete fds;
    } 
    for (int i=0;i<m_AllocatedDescriptors;i++)
      m_Descriptors[i]=fd[i];
    if (fd)delete []fd;
  }
  covDim();
}
void MapPoseTree::addDescriptor(FeatureDescriptor *d, unsigned short branch)
{ 
  if (branch>m_AllocatedDescriptors){
    reallocate();
    if (branch>m_AllocatedDescriptors)
      return;
  }
  m_Descriptors[branch]->add(d);
}
void MapPoseTree::removeDescriptor(FeatureDescriptor *d, 
				   unsigned short branch)
{ 
  if (branch>m_AllocatedDescriptors){
    reallocate();
    if (branch>m_AllocatedDescriptors)
      return;
  }
  m_Descriptors[branch]->remove(d);
}
void MapPoseTree::init(){
  Index=0;
  m_SizeIndex=0;
  m_Descriptors=0;
  m_AllocatedDescriptors=0;
  PoseTree::init(); 
}
  
 void MapPoseTree::clean(){
  for (int i=0;i<m_AllocatedDescriptors;i++)
     { 
       FeatureDescriptors *fd=m_Descriptors[i];
       for (int j=0;j<fd->m_NumberDesc;j++){
	 FeatureDescriptor *des=fd->m_Desc[j];
	 if (des){
	   RelativePointFeature *mf=des->getRelativePointFeature();	   
	   if(mf)mf->remove(des);
	 }
       }
       delete m_Descriptors[i];
     }
   if (m_Descriptors)
     delete []m_Descriptors;
   m_Descriptors=0;
   m_AllocatedDescriptors=0;
   PoseTree::clean();
 }
void MapPoseTree::get(GenericData &gd)
{
  MapObject::get(gd);
  calcPoses();
  int c=dim();
  int dm=c;
  c*=2;
  if (c<7)c=7;
  if (c<gd.ShortData.Columns)c=gd.ShortData.Columns;
  gd.forceShortDataSize(2, c);
  for (int i=0; i<dm;i++){
    //Index[i]=0;
    if (Index)  gd.ShortData.setLong(Index[i],1,2*i);
    else {
      gd.ShortData(1,2*i)=0;
      gd.ShortData(1,2*i+1)=0;
    }
  }  
  ShortMatrix temp;
  PoseTree::save(temp,gd.Data);
  gd.ShortData(0,6)=temp.Columns;
  if (temp.Columns>c)c=temp.Columns;
  gd.forceShortDataSize(2+temp.Rows,c);
  gd.ShortData.offset(2,0,temp.Rows,temp.Columns);
  gd.ShortData=temp;
  gd.ShortData.offset(-2,0,temp.Rows+2,c);
  gd.setTime(time());
}
int MapPoseTree::set(GenericData &gd)
{
  if (MapObject::set(gd))return 1;
  for (int i=0;i<m_AllocatedDescriptors;i++)
    { 
       FeatureDescriptors *fd=m_Descriptors[i];
       for (int j=0;j<fd->m_NumberDesc;j++){
	 FeatureDescriptor *des=fd->m_Desc[j];
	 if (des){
	   RelativePointFeature *mf=des->getRelativePointFeature();	   
	   if(mf)mf->remove(des);
	 }
       }
       delete m_Descriptors[i];
    }
  delete []m_Descriptors;
  m_Descriptors=0;
  m_AllocatedDescriptors=0;
  PoseTree::clean();
  long r=gd.ShortData.Rows;
  long c=gd.ShortData.Columns;
  if (c<7)return 1;
  if (r<2)return 1;
  gd.ShortData.offset(2,0,r-2,gd.ShortData(0,6));
  if (PoseTree::restore(gd.ShortData,gd.Data,gd.Time)){
    gd.ShortData.offset(-2,0,r,c);
    return 1;
  }
  gd.ShortData.offset(-2,0,r,c);
  reallocate();
  int dm=dim();
  if (dm)
    for (int i=0; i<dm;i++){
      if (Index)  Index[i]=gd.ShortData.getLong(1,2*i);
    }  
  return 0;
}
void MapPoseTree::print(){
  std::cerr<<"MapPoseTree "<<Key<<" "<<m_CovDim<<"\n";
  PoseTree::print();

}
