//
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "MapPointFeature.hh"
using namespace Cure;

MapPointFeature::MapPointFeature():MapFeature(0)
{
  FullDim=3;
  NumberPoints=1;
  Number3D=1;
  Points=AllocatedPoints;
  Index=AllocatedIndex;
  Type=MAPPOINTFEATURE_TYPE;
  TotalWgt=0;
  VCount=0;
  init();
  CastPtr=0;
  TriangleThreshold=.05;
  WeightThreshold=1;
  DistanceThreshold=20;
  TrackThreshold=.1;
}
MapPointFeature::MapPointFeature(MapPointFeature *wp,
				   MapBank *b):MapFeature(b)
{
  FullDim=3;
  NumberPoints=1;
  Number3D=1;
  Points=AllocatedPoints;
  Index=AllocatedIndex;
  Type=MAPPOINTFEATURE_TYPE;
  TotalWgt=0;
  VCount=0;
  init();
  if (wp)
     {
       TrackThreshold=wp->TrackThreshold;
       CastPtr=wp->CastPtr;
       TriangleThreshold=wp->TriangleThreshold;
       WeightThreshold=wp->WeightThreshold;
       DistanceThreshold=wp->DistanceThreshold;
     }
  else
    {
      CastPtr=0;
      TriangleThreshold=.05;
      WeightThreshold=1;
      DistanceThreshold=20;
      TrackThreshold=.1;
    } 
}

void MapPointFeature::initializeFromPixels(Transformation3D &t,
					    const double centerpixels[2],
					    double focallength,
					    double distanceguess)
{
  if (!Points[0])
    {
      MapPoint *p=new MapPoint(Bank);
      setCenter(p);
    }
  Points[0]->initializeFromPixels(t,centerpixels,focallength,distanceguess);
}


int  MapPointFeature::extend()
{
  if (Bdual.Columns==3)return 0;
  if (TotalWgt>WeightThreshold)
    {
      Bdual.reallocate(3,3);
      Bdual=1;
      return 3;
    }
  return 0;  
}

void  MapPointFeature::forceExtend()
{
  if (Bdual.Columns==3)return;
  Bdual.reallocate(3,3);
  Index[0]=-1;
  Bdual=1;
}
void MapPointFeature::prune(double minDistance)
{
  LinkedArray *la=&Vectors;
  if (!la->Next)return;
  if((*la)(0)<minDistance)
    {
      while(la->Next)
	{
	  if((*la)(0)<minDistance)la=la->Next;
	  else break;
	}
      if (la->Next)
	{
	  Vectors.reallocate(la->Length);
	  for (int i=0; i<Vectors.Length; i++)Vectors(i)=(*la)(i);
	  la->cut();
	  delete la;
	  VCount--;
	  la=&Vectors;
	}
      else 
	{
	  Vectors.clear();
	  VCount=0;
	  return;
	}
    }
  while(la->Next)
    {
      if((*la)(0)<minDistance)
	{
	  LinkedArray *ta=la;
	  la=la->Next;
	  ta->cut();
	  delete ta;
	  VCount--;
	}
      else la=la->Next;
    }
}

unsigned short MapPointFeature::getMeasurementType(unsigned short type)
{
  if ((type!=32))type=0;
  if (!type)return type;
  if (Bdual.Columns==3) type=1;
  else type=0;
  return type;  
}

int  MapPointFeature::addInfo(const Cure::Matrix & v,
			      Transformation3D & map2info,
			      const int type)
{

  if (Bdual.Columns==3)return 0;
  Matrix m(1,11); 
  //distance, w, x_i, s_i, s*x_i, r_i, sum_w  (0..10)
  m(0,0)=v(0,0);
  m(0,1)=v(0,1);
  m(0,2)=v(0,2);
  m(0,3)=v(0,3);
  m(0,4)=v(0,4);
  //  map2info.transform(m.Element+2,m.Element+2);
  double d=sqrt(v(0,7)*v(0,7)+v(0,5)*v(0,5)+v(0,6)*v(0,6));
  m(0,5)=v(0,5)/d;
  m(0,6)=v(0,6)/d;;
  m(0,7)=v(0,7)/d;
  m(0,8)=m(0,2)*m(0,5)+m(0,3)*m(0,6)+m(0,4)*m(0,7);
  m(0,9)=0;
  m(0,10)=0;
  if ( LastDistance<=v(0,0))
    {
      LastDistance=v(0,0);
      double mindistance=v(0,0)-DistanceThreshold;
      prune(mindistance);
      map2info.invRotate(m.Element+5, LastBearing);
    }
  else
    {
      LinkedArray *la=&Vectors;
      while(la->Next)
	{
	  if((*la)(0)!=v(0,0))la=la->Next;
	  else break;
	}
      if((*la)(0)!=v(0,0))return 1;
      while(la->Next)
	{
	  if((*la)(0)==v(0,0))
	    {
	      LinkedArray *ta=la;
	      la=la->Next;
	      ta->cut();
	      delete ta;
	      VCount--;
	    }
	  else la=la->Next;
	}    
    }
  TotalWgt=0;
  for (LinkedArray *vlist=&Vectors;vlist->Next;vlist=vlist->Next)
    {
      double w=m(0,1)*(*vlist)(1);
      double temp[3];
      temp[0]=m(0,2)-(*vlist)(2);
      temp[1]=m(0,3)-(*vlist)(3);
      temp[2]=m(0,4)-(*vlist)(4);
      double d=temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2];
      if (d>TriangleThreshold)
	{
	  double trig[2];
	  trig[0]=(*vlist)(5)*m(0,5)+(*vlist)(6)*m(0,6)+(*vlist)(7)*m(0,7);
	  trig[1]=temp[0]*(*vlist)(5)+temp[1]*(*vlist)(6)+
	    temp[2]*(*vlist)(7);
	  trig[0]=(1-trig[0]*trig[0]);
	  (*vlist)(9)+=(w*sqrt(trig[0]*(d-trig[1]*trig[1])));
	  (*vlist)(10)+=(w*trig[0]);
	}
      TotalWgt+=(*vlist)(10);
    }
  VCount++;
  Vectors.add(m);
  if ( LastDistance!=v(0,0))return 0;
  if (TotalWgt>WeightThreshold)
    {
      double oldx[4],newx[4];
      double oldev[9],oldlambda[3],ev[9],lambda[3];
      double oldquality=estimate_it(oldx,oldev,oldlambda);
      //      std::cerr<<"oldx "<<oldx[0]<<" "<<oldx[1]<<" "<<oldx[2]<<" ";
      // Points[0]->setXYZ(oldx);
      
      double newquality=2*oldquality;//reestimate_it(oldx,newx,ev,lambda);
      //std::cerr<<"oldx "<<newx[0]<<" "<<newx[1]<<" "<<newx[2]<<" ";
      if (newquality>oldquality)
	for (int i=0;i<3;i++)
	  {
	    lambda[i]=oldlambda[i];
	    newx[i]=oldx[i];
	    ev[3*i]=oldev[3*i];
	    ev[3*i+1]=oldev[3*i+1];
	    ev[3*i+2]=oldev[3*i+2];
	  }
      map2info.invTransform(newx,newx);
      Points[0]->setXYZ(newx);
      map2info.invRotate(ev,ev);
      map2info.invRotate(ev+3,ev+3);
      map2info.invRotate(ev+6,ev+6);
      Info=1;
      lambda[0]*=2;
      if (lambda[0]<1E-2)lambda[0]=1E-2;
      double d=lambda[0]*1E-6;
      if (lambda[1]<d)lambda[1]=d;
      if (lambda[2]<d)lambda[2]=d;
      
      for (int i=0; i<3; i++)
	for (int j=0;j<3;j++)
	  {
	    Info(i,j)=ev[3+i]*ev[3+j]/lambda[1]+
	      ev[6+i]*ev[6+j]/lambda[2]+ev[i]*ev[j]/lambda[0];
	  }
      return 0;
    }
  map2info.invTransform(m.Element+2,m.Element+2);
  map2info.invRotate(m.Element+5,m.Element+5);
  m(0,8)=m(0,2)*m(0,5)+m(0,3)*m(0,6)+m(0,4)*m(0,7);
  trackPoint(m);
  return 0;
}

int  MapPointFeature::addInfo(const Cure::Matrix & v,
			      const int type)
{

  if (Bdual.Columns==3)return 0;
  Matrix m(1,11); 
  //distance, w, x_i, s_i, s*x_i, r_i, sum_w  (0..10)
  m(0,0)=v(0,0);
  m(0,1)=v(0,1);
  m(0,2)=v(0,2);
  m(0,3)=v(0,3);
  m(0,4)=v(0,4);
  double d=sqrt(v(0,7)*v(0,7)+v(0,5)*v(0,5)+v(0,6)*v(0,6));
  m(0,5)=v(0,5)/d;
  m(0,6)=v(0,6)/d;;
  m(0,7)=v(0,7)/d;
  m(0,8)=m(0,2)*m(0,5)+m(0,3)*m(0,6)+m(0,4)*m(0,7);
  m(0,9)=0;
  m(0,10)=0;
  if ( LastDistance<=v(0,0))
    {
      LastDistance=v(0,0);
      double mindistance=v(0,0)-DistanceThreshold;
      prune(mindistance);
      //  map2info.invRotate(m.Element+5, LastBearing);
      LastBearing[0]=m.Element[5];
      LastBearing[1]=m.Element[6];
      LastBearing[2]=m.Element[7];
    }
  else
    {
      LinkedArray *la=&Vectors;
      while(la->Next)
	{
	  if((*la)(0)!=v(0,0))la=la->Next;
	  else break;
	}
      if((*la)(0)!=v(0,0))return 1;
      while(la->Next)
	{
	  if((*la)(0)==v(0,0))
	    {
	      LinkedArray *ta=la;
	      la=la->Next;
	      ta->cut();
	      delete ta;
	      VCount--;
	    }
	  else la=la->Next;
	}    
    }
  TotalWgt=0;
  for (LinkedArray *vlist=&Vectors;vlist->Next;vlist=vlist->Next)
    {
      double w=m(0,1)*(*vlist)(1);
      double temp[3];
      temp[0]=m(0,2)-(*vlist)(2);
      temp[1]=m(0,3)-(*vlist)(3);
      temp[2]=m(0,4)-(*vlist)(4);
      double d=temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2];
      if (d>TriangleThreshold)
	{
	  double trig[2];
	  trig[0]=(*vlist)(5)*m(0,5)+(*vlist)(6)*m(0,6)+(*vlist)(7)*m(0,7);
	  trig[1]=temp[0]*(*vlist)(5)+temp[1]*(*vlist)(6)+
	    temp[2]*(*vlist)(7);
	  trig[0]=(1-trig[0]*trig[0]);
	  (*vlist)(9)+=(w*sqrt(trig[0]*(d-trig[1]*trig[1])));
	  (*vlist)(10)+=(w*trig[0]);
	}
      TotalWgt+=(*vlist)(10);
    }
  VCount++;
  Vectors.add(m);
  if ( LastDistance!=v(0,0))return 0;
  if (TotalWgt>WeightThreshold)
    {
      double oldx[4],newx[4];
      double oldev[9],oldlambda[3],ev[9],lambda[3];
      double oldquality=estimate_it(oldx,oldev,oldlambda);
      //      std::cerr<<"oldx "<<oldx[0]<<" "<<oldx[1]<<" "<<oldx[2]<<" ";
      // Points[0]->setXYZ(oldx);
      
      double newquality=2*oldquality;//reestimate_it(oldx,newx,ev,lambda);
      //std::cerr<<"oldx "<<newx[0]<<" "<<newx[1]<<" "<<newx[2]<<" ";
      if (newquality>oldquality)
	for (int i=0;i<3;i++)
	  {
	    lambda[i]=oldlambda[i];
	    newx[i]=oldx[i];
	    ev[3*i]=oldev[3*i];
	    ev[3*i+1]=oldev[3*i+1];
	    ev[3*i+2]=oldev[3*i+2];
	  }
      //map2info.invTransform(newx,newx);
      Points[0]->setXYZ(newx);
      // map2info.invRotate(ev,ev);
      // map2info.invRotate(ev+3,ev+3);
      // map2info.invRotate(ev+6,ev+6);
      Info=1;
      lambda[0]*=2;
      if (lambda[0]<1E-2)lambda[0]=1E-2;
      double d=lambda[0]*1E-6;
      if (lambda[1]<d)lambda[1]=d;
      if (lambda[2]<d)lambda[2]=d;
      
      for (int i=0; i<3; i++)
	for (int j=0;j<3;j++)
	  {
	    Info(i,j)=ev[3+i]*ev[3+j]/lambda[1]+
	      ev[6+i]*ev[6+j]/lambda[2]+ev[i]*ev[j]/lambda[0];
	  }
      return 0;
    }
  //  map2info.invTransform(m.Element+2,m.Element+2);
  //map2info.invRotate(m.Element+5,m.Element+5);
  m(0,8)=m(0,2)*m(0,5)+m(0,3)*m(0,6)+m(0,4)*m(0,7);
  trackPoint(m);
  return 0;
}
/**
 * v (5..7) UNIT bearing vector in map frame (magnitude must be one.
 */
int  MapPointFeature::addFullInfo(const Cure::Matrix & v, double minAngle)
{
  Vectors.clear();
  VCount=0;
  if (Bdual.Columns==3)return 0;
  int rows=v.Rows;
  LastDistance=0;
  for (int rw=0;rw<rows;rw++)
    if ( LastDistance<=v(rw,0))
      {
	LastDistance=v(rw,0);
	LastBearing[0]=v(rw,5);
	LastBearing[1]=v(rw,6);
	LastBearing[2]=v(rw,7);
      }
  double mindistance=LastDistance-DistanceThreshold;  
  bool newway=true;
  if (newway){
    double testtri=0;
  for (int rw=0;rw<rows;rw++)
    {
      if (v(rw,0)>=mindistance){
	Matrix m(1,11); 
	//distance, w, x_i, s_i, s*x_i, r_i, sum_w  (0..10)
	m(0,0)=v(rw,0);
	m(0,1)=v(rw,1);
	m(0,2)=v(rw,2);
	m(0,3)=v(rw,3);
	m(0,4)=v(rw,4);
	m(0,5)=v(rw,5);
	m(0,6)=v(rw,6);;
	m(0,7)=v(rw,7);
	m(0,8)=m(0,2)*m(0,5)+m(0,3)*m(0,6)+m(0,4)*m(0,7);
	m(0,9)=0;
	m(0,10)=0;
	
	for (LinkedArray *vlist=&Vectors;vlist->Next;vlist=vlist->Next)
	  {
	    double w=m(0,1)*(*vlist)(1);
	    double temp[3];
	    temp[0]=m(0,2)-(*vlist)(2);
	    temp[1]=m(0,3)-(*vlist)(3);
	    temp[2]=m(0,4)-(*vlist)(4);
	    double d=temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2];
	    if (d>TriangleThreshold)
	      {
		double trig[2];
		trig[0]=(*vlist)(5)*m(0,5)+(*vlist)(6)*m(0,6)+(*vlist)(7)*m(0,7);
		trig[1]=temp[0]*(*vlist)(5)+temp[1]*(*vlist)(6)+
		  temp[2]*(*vlist)(7);
		trig[0]=(1-trig[0]*trig[0]);
		(*vlist)(9)+=(w*sqrt(trig[0]*(d-trig[1]*trig[1])));
		(*vlist)(10)+=(w*trig[0]);
		if (trig[0]>testtri)testtri=trig[0];
	      }
	  }
	VCount++;
	Vectors.add(m);
      }
    }
  TotalWgt=0;
  for (LinkedArray *vlist=&Vectors;vlist->Next;vlist=vlist->Next)
    TotalWgt+=(*vlist)(10);
  if (TotalWgt>WeightThreshold)
    {
      double oldx[4],newx[4];
      double oldev[9],oldlambda[3],ev[9],lambda[3];
      double oldquality=estimate_it(oldx,oldev,oldlambda);
      //      std::cerr<<"oldx "<<oldx[0]<<" "<<oldx[1]<<" "<<oldx[2]<<" ";
      // Points[0]->setXYZ(oldx);
      
      double newquality=2*oldquality;//reestimate_it(oldx,newx,ev,lambda);
      //std::cerr<<"oldx "<<newx[0]<<" "<<newx[1]<<" "<<newx[2]<<" ";
      if (newquality>oldquality)
	for (int i=0;i<3;i++)
	  {
	    lambda[i]=oldlambda[i];
	    newx[i]=oldx[i];
	    ev[3*i]=oldev[3*i];
	    ev[3*i+1]=oldev[3*i+1];
	    ev[3*i+2]=oldev[3*i+2];
	  }
      Points[0]->setXYZ(newx);
      Info=1;
      lambda[0]*=2;
      if (lambda[0]<1E-2)lambda[0]=1E-2;
      double d=lambda[0]*1E-6;
      if (lambda[1]<d)lambda[1]=d;
      if (lambda[2]<d)lambda[2]=d;
      
      for (int i=0; i<3; i++)
	for (int j=0;j<3;j++)
	  {
	    Info(i,j)=ev[3+i]*ev[3+j]/lambda[1]+
	      ev[6+i]*ev[6+j]/lambda[2]+ev[i]*ev[j]/lambda[0];
	  }
      if (testtri<minAngle)TotalWgt=0;
      return 0;
    }
  return 0;
}else{
  for (int rw=0;rw<rows;rw++)
    {
      if (v(rw,0)>=mindistance){
	Matrix m(1,11); 
	//distance, w, x_i, s_i, s*x_i, r_i, sum_w  (0..10)
	m(0,0)=v(rw,0);
	m(0,1)=v(rw,1);
	m(0,2)=v(rw,2);
	m(0,3)=v(rw,3);
	m(0,4)=v(rw,4);
	m(0,5)=v(rw,5);
	m(0,6)=v(rw,6);;
	m(0,7)=v(rw,7);
	m(0,8)=m(0,2)*m(0,5)+m(0,3)*m(0,6)+m(0,4)*m(0,7);
	m(0,9)=0;
	m(0,10)=0;
	
	for (LinkedArray *vlist=&Vectors;vlist->Next;vlist=vlist->Next)
	  {
	    double w=m(0,1)*(*vlist)(1);
	    double temp[3];
	    temp[0]=m(0,2)-(*vlist)(2);
	    temp[1]=m(0,3)-(*vlist)(3);
	    temp[2]=m(0,4)-(*vlist)(4);
	    double d=temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2];
	    if (d>TriangleThreshold)
	      {
		double trig[2];
		trig[0]=(*vlist)(5)*m(0,5)+(*vlist)(6)*m(0,6)+(*vlist)(7)*m(0,7);
		trig[1]=temp[0]*(*vlist)(5)+temp[1]*(*vlist)(6)+
		  temp[2]*(*vlist)(7);
		trig[0]=(1-trig[0]*trig[0]);
		(*vlist)(9)+=(w*sqrt(trig[0]*(d-trig[1]*trig[1])));
		(*vlist)(10)+=(w*trig[0]);
	      }
	  }
	VCount++;
	Vectors.add(m);
      }
    }
  TotalWgt=0;
  for (LinkedArray *vlist=&Vectors;vlist->Next;vlist=vlist->Next)
    TotalWgt+=(*vlist)(10);
  if (TotalWgt>WeightThreshold)
    {
      double oldx[4],newx[4];
      double oldev[9],oldlambda[3],ev[9],lambda[3];
      double oldquality=estimate_it(oldx,oldev,oldlambda);
      //      std::cerr<<"oldx "<<oldx[0]<<" "<<oldx[1]<<" "<<oldx[2]<<" ";
      // Points[0]->setXYZ(oldx);
      
      double newquality=2*oldquality;//reestimate_it(oldx,newx,ev,lambda);
      //std::cerr<<"oldx "<<newx[0]<<" "<<newx[1]<<" "<<newx[2]<<" ";
      if (newquality>oldquality)
	for (int i=0;i<3;i++)
	  {
	    lambda[i]=oldlambda[i];
	    newx[i]=oldx[i];
	    ev[3*i]=oldev[3*i];
	    ev[3*i+1]=oldev[3*i+1];
	    ev[3*i+2]=oldev[3*i+2];
	  }
      Points[0]->setXYZ(newx);
      Info=1;
      lambda[0]*=2;
      if (lambda[0]<1E-2)lambda[0]=1E-2;
      double d=lambda[0]*1E-6;
      if (lambda[1]<d)lambda[1]=d;
      if (lambda[2]<d)lambda[2]=d;
      
      for (int i=0; i<3; i++)
	for (int j=0;j<3;j++)
	  {
	    Info(i,j)=ev[3+i]*ev[3+j]/lambda[1]+
	      ev[6+i]*ev[6+j]/lambda[2]+ev[i]*ev[j]/lambda[0];
	  }
      return 0;
    }
  return 0;
}
}
  //distance, w, x_i, s_i, s*x_i, r_i, sum_w  (0..10)
void MapPointFeature::trackPoint(Matrix &m)
{

  if (Bdual.Columns==3)return;
  Matrix x(3,1);
  getX(x);
  double sxf=x(0,0)*m(0,5)+x(1,0)*m(0,6)+x(2,0)*m(0,7);
  double ss=m(0,5)*m(0,5)+m(0,6)*m(0,6)+m(0,7)*m(0,7);
  if (ss<1E-20)return;
  double t=sxf-m(0,8);
  t/=ss;
  x(0,0)=m(0,2)+m(0,5)*t;
  x(1,0)=m(0,3)+m(0,6)*t;
  x(2,0)=m(0,4)+m(0,7)*t;
  setX(x);
}

int MapPointFeature::eigen_it(Matrix & cloud, double wsums[4],
			  double ev[9], double lambda[3])
{
  wsums[0]=0;
  wsums[1]=0;
  wsums[2]=0;
  wsums[3]=0;
  Matrix cov(3);
  cov=0;
  for (int i=0; i<cloud.Rows; i++)
    {
      wsums[0]+=cloud(i,0)*cloud(i,3);
      wsums[1]+=cloud(i,1)*cloud(i,3);
      wsums[2]+=cloud(i,2)*cloud(i,3);
      wsums[3]+=cloud(i,3);
      cov(0,0)+=cloud(i,0)*cloud(i,0)*cloud(i,3);
      cov(0,1)+=cloud(i,0)*cloud(i,1)*cloud(i,3);
      cov(0,2)+=cloud(i,0)*cloud(i,2)*cloud(i,3);
      cov(1,1)+=cloud(i,1)*cloud(i,1)*cloud(i,3);
      cov(1,2)+=cloud(i,1)*cloud(i,2)*cloud(i,3);
      cov(2,2)+=cloud(i,3)*cloud(i,3)*cloud(i,3);
    }
  
  wsums[0]/=wsums[3];
  wsums[1]/=wsums[3];
  wsums[2]/=wsums[3];
  for (int i=0; i<3; i++)
    for (int j=i;j<3;j++)
      {
	cov(i,j)/=wsums[3];	   
	cov(i,j)-=(wsums[i]*wsums[j]);
	if (i!=j)cov(j,i)=cov(i,j);
      }
  Matrix lam,evec;
  cov.symmetricEigen(lam,evec);
  for (int i=0; i<3; i++)
    {
      lambda[i]=lam(i,i);
      for (int j=0; j<3; j++)
	  ev[3*i+j]=evec(j,i);
    }  
  return 0;
}
double MapPointFeature::reestimate_it(double oldx[4],double newx[4],double ev[9],
			      double lambda[3])
{
  int i=0;
  Matrix cloud(VCount,4);
  for (LinkedArray *vlist=&Vectors;vlist->Next;vlist=vlist->Next, i++)
    {
      cloud(i,3)=(*vlist)(1);
      if (cloud(i,3)>1E-9)
	{
	  double t=(*vlist)(5)*oldx[0]+(*vlist)(6)*oldx[1]+(*vlist)(7)*oldx[2];
	  t-=(*vlist)(8);
	  cloud(i,0)=(*vlist)(2)+(*vlist)(5)*t;
	  cloud(i,1)=(*vlist)(3)+(*vlist)(6)*t;
	  cloud(i,2)=(*vlist)(4)+(*vlist)(7)*t;
	}
      else
	{
	  cloud(i,3)=0;
	  cloud(i,0)=0;
	  cloud(i,1)=0;
	  cloud(i,2)=0;
	}
    }
  if (eigen_it(cloud,newx,ev,lambda))
    {
      Info=1E-10;
      return 1E30;
    }	
  if (lambda[0]<1E-4)lambda[0]=1E-4;
  if (lambda[1]<1E-4)lambda[1]=1E-4;
  if (lambda[2]<1E-4)lambda[2]=1E-4;
  lambda[0]+=1E-4; //sensor error
  lambda[1]+=1E-4; //sensor error
  lambda[2]+=1E-4;
  return (lambda[0]*lambda[1]*lambda[2]);
}
double MapPointFeature::estimate_it(double wsums[4],double ev[9],double lambda[3])
{	
  int i=0;
  Matrix cloud(VCount,4);
  for (LinkedArray *vlist=&Vectors;vlist->Next;vlist=vlist->Next, i++)
    {
      cloud(i,3)=(*vlist)(10);
      if (cloud(i,3)>1E-9)
	{
	  cloud(i,0)=(*vlist)(9)*(*vlist)(5)/(*vlist)(10)+(*vlist)(2);
	  cloud(i,1)=(*vlist)(9)*(*vlist)(6)/(*vlist)(10)+(*vlist)(3);
	  cloud(i,2)=(*vlist)(9)*(*vlist)(7)/(*vlist)(10)+(*vlist)(4);
	}
      else
	{
	  cloud(i,3)=0;
	  cloud(i,0)=0;
	  cloud(i,1)=0;
	  cloud(i,2)=0;
	}
    }

  if (eigen_it(cloud,wsums,ev,lambda))
    {
      for (int i=0; i<9;i++)ev[i]=0;
      ev[0]=1;
      lambda[0]=4E30;
      ev[4]=1;
      lambda[1]=2E30;
      ev[8]=1;
      lambda[2]=1E30;
      return 8E90;
    }	
  if (lambda[0]<1E-4)lambda[0]=1E-4;
  if (lambda[1]<1E-4)lambda[1]=1E-4;
  if (lambda[2]<1E-4)lambda[2]=1E-4;
  lambda[0]+=1E-4; //sensor error
  lambda[1]+=1E-4; //sensor error
  lambda[2]+=1E-4;
  return (lambda[0]*lambda[1]*lambda[2]);
}
int MapPointFeature::merge(MapPointFeature *mf, unsigned short typ)
{
  if (!mf)return -1;
  return 0;
}

void MapPointFeature::getB(Cure::Matrix &b)const
{
  b.reallocate(Bdual.Columns,Bdual.Rows);
  b=1; 
}
bool MapPointFeature::getC(MapPointFeature *mw,Cure::Matrix &c)const
{
  if ((Bdual.Columns!=3)||(mw->Bdual.Columns!=3))return false;
  c.reallocate(3);
  c=1;
  return true; 
}
bool MapPointFeature::testMatch(MapPointFeature *mw, double tolerance)
{
  if ((Bdual.Columns!=3)||(mw->Bdual.Columns!=3))return false;
  return true;
}
int MapPointFeature::set(GenericData &gd)
{
  if (MapFeature::set(gd))return 1;
  TotalWgt=0;
  VCount=0;
  TriangleThreshold=.05;
  WeightThreshold=1;
  DistanceThreshold=20;
  TrackThreshold=.1;
  return 0;
}

