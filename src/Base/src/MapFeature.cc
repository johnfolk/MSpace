// = RCSID
//    $Id: MapFeature.cc ,v 1.1 2004/04/1 
//
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "MapFeature.hh"
using namespace Cure;
Cure::MapFeature::MapFeature(MapBank *b):MapObject(b)
{
  LastBearing[0]=0;
  LastBearing[1]=0;
  LastBearing[2]=1;
  LastDistance=-1;
  DistanceThreshold=0;
  FullDim=0;
  NumberPoints=0;
  Number3D=0;
  Number2D=0;
  NumberScalars=0;
  NumberSharedPoints=0;
  Consecutive=1;
  Scalars=0;
  Points=0;
  PSharedIndex=0;
 }
void Cure::MapFeature::init()
{
  for (int i=0; i<NumberPoints; i++)
    Points[i]=0;
  for (int i=0; i<NumberScalars; i++)
    Scalars[i]=0;
  for (int i=0; i<FullDim;i++)Index[i]=0;
  Info.reallocate(FullDim);
  Bdual.reallocate(FullDim,0);
  Info=1;
}
Cure::MapFeature::~MapFeature()
{
  if (PSharedIndex)delete [] PSharedIndex;
  PSharedIndex=0;
  for (int i=0;i<(Number3D+Number2D);i++)
    if (Points[i]!=0)
      {
	Points[i]->removeFeature(this);
	if (!(Points[i]->Features.Next))delete Points[i];
	Points[i]=0;
      }
}
void MapFeature::write(std::fstream &fs )
{
  fs<<Type<<" "<<Key<<" "<<ID<<"\n";
  for (int i=0; i<NumberPoints; i++){
    Points[i]->write(fs);
  }
  for (int i=0; i<NumberScalars;i++)
    fs<<Scalars[i]<<"\n";
  fs<<NumberSharedPoints<<" "<<Consecutive<<" "<<LastDistance<<" "
    <<Bdual.Columns<<"\n";
  long *indx=getIndex();
  for (int i=0; i<FullDim; i++)
    fs<<indx[i]<<" ";
  fs<<"\n";
}
int  MapFeature::read(int version,std::fstream &fs, bool readKey)
{
  if (version!=1)return 1;
  long key=0;
  if (readKey){
    if (!(fs>>key))return 1;
    if (key!=Key){
      if (Bank){
	MapObject *m=Bank->getMapObject(key);
	if (m){
	  std::cerr<<"MapObject with this key existes "<<key<<"\n";
	}
	Bank->remove(this);
	Bank->add(this,key);
	if (key!=Key){
	  std::cerr<<"wanted Feature key "<<key<<" got "<<Key;
	}
      }
      else
	Key=key;
    }
  }

  fs>>ID;
  for (int i=0;i<NumberPoints; i++){
    fs>>key;
    MapPoint *mp=0;
    if (Bank){ 
      MapObject *m=Bank->getMapObject(key);
      if (m)
	{
	  mp=m->narrowPoint();	
	}
    }
    if (!mp){
      mp=new MapPoint();
      mp->Bank=Bank;
      if (Bank){
	Bank->add(mp,key);    
	if (key!=mp->Key){
	  std::cerr<<"wanted point key "<<key<<" got "<<mp->Key;
	}
      }      
    }
    mp->read(fs);
    if (mp)addPoint(mp,i);
  }
  for (int i=0; i<NumberScalars;i++)
    fs>>Scalars[i];
  fs>>NumberSharedPoints>>Consecutive>>LastDistance;
  fs>>key;
  Bdual.reallocate(FullDim,key);
  long *indx=getIndex();
  for (int i=0; i<FullDim; i++)
    fs>>indx[i];
  return 0;
}

int Cure::MapFeature::addPoint(MapPoint *f,unsigned int which)
{
  if (which>(unsigned int)(NumberPoints-1))return 1;
  if (Points[which])
    {
      Points[which]->removeFeature(this);
      Points[which]=0;
    }
  f->addFeature(this);
  Points[which]=f;
  setNumberSharedPoints();
  return 0;
}

void Cure::MapFeature::setNumberSharedPoints()
  {
    NumberSharedPoints=0;
    for (int i=0; i<NumberPoints; i++)
      if (Points[i])
	if (Points[i]->numberOfFeatures()>1)NumberSharedPoints++;
   if (NumberSharedPoints) 
     if (PSharedIndex==0)
       {
	 PSharedIndex=new short[FullDim];
	 for (int i=0; i<FullDim; i++)
	   PSharedIndex[i]=-1;
       }
  }
int Cure::MapFeature::removePoint(MapPoint *f)
{
  int r=0;
  if (!f)return r;
  f->removeFeature(this);
  for (int i=0; i<NumberPoints; i++)
    if (Points[i]==f)
      {
	Points[i]=0;
	r++;
      }
  setNumberSharedPoints();
  return r;
}
MapPoint * Cure::MapFeature::removePoint(unsigned short which)
{
  if (which>NumberPoints)return 0;
  MapPoint *f=Points[which];
  if (!f)return f;
  Points[which]=0;
  int r=0;
  for (int i=0; i<NumberPoints; i++)
    if (Points[i]==f)r++;
  if (r==0)
    f->removeFeature(this);
  setNumberSharedPoints();
  return f;
}
int Cure::MapFeature::removeObject(MapObject *f)
{
  MapPoint *p=f->narrowPoint();
  return removePoint(p);
}
void Cure::MapFeature::updateP(Matrix & dp)
{
  // Calculate the change in the "real" feature coordinates based on
  // the change in the P-space
  Matrix dx(FullDim,1);
  dx.multiply_(Bdual,dp);

  // Add this change to the points that are used to parameterize this
  // feature
  updateX(dx);

  dp=0;
}

void Cure::MapFeature::updateX(const Matrix & dx)
{
  Vector3D newdx;
  int ind=0;
  int i=0;

  // Loop through the 3D point coordinates
  for (; i<Number3D; i++,ind+=3)
    {
      double *d=&dx.Element[ind*dx.RowInc];
      newdx.setXYZ(d);
      Points[i]->Vector3D::operator+=(newdx);
    }

  // Loop through the 2D points and make sure that the Z-coord is 0
  // before doing so
  newdx.X[2]=0;
  for (; i<NumberPoints; i++,ind+=2)
    {
      double *d=&dx.Element[ind*dx.RowInc];
      newdx.setXY(d);
      Points[i]->Vector3D::operator+=(newdx);
     }

  // Loop through the scalar coordinates
  for (int i=0; i<NumberScalars;i++,ind++)
    Scalars[i]+=dx(ind,0);

  // Recenter, i.e. recalculate Bdual if needed
  recenter();
}

void Cure::MapFeature::setX(const  Matrix & x)
{
  int ind=0;
  int i=0;
  for (; i<Number3D; i++,ind+=3)
    Points[i]->setXYZ(&x.Element[ind*x.RowInc]);
  for (; i<NumberPoints; i++,ind+=2)
    Points[i]->setXY(&x.Element[ind*x.RowInc]);
  for (i=0; i<NumberScalars;i++,ind++)
    Scalars[i]=x(ind,0);
  recenter();
}
void Cure::MapFeature::getX(Matrix & x)const
{
  x.reallocate(FullDim,1);
  int ind=0;
  int i=0;
  for (; i<Number3D; i++,ind+=3)
    Points[i]->getXYZ(&x.Element[ind*x.RowInc]);
  for (; i<NumberPoints; i++,ind+=2)
    Points[i]->getXY(&x.Element[ind*x.RowInc]);
  for (i=0; i<NumberScalars;i++,ind++)
    x(ind,0)=Scalars[i];
}
void Cure::MapFeature::getP(Matrix  & p)const
{
  Matrix x;
  getX(x);
  p.reallocate(Bdual.Columns,1);
  Matrix b(Bdual.Columns,Bdual.Rows);
  getB(b);
  p.multiply_(b,x);
}
int  Cure::MapFeature::extend(Cure::Matrix & b_n,   
			Cure::Matrix & invcov)
{
  int res=extend();
  Matrix b(Bdual.Columns,Bdual.Rows);
  getB(b);
  if (res==b.Rows)
    {
      b_n .multTranspose_(Info,b,2);
      invcov.multiply_(b,b_n);
      b_n=b;
    }
  else if (res>0)
    {
      int oldr=b.Rows-res;
      int temp=b.RowInc*oldr;
      b.Element+=temp;
      b.Rows=res;
      b_n.reallocate(res,b.Columns);
      b_n=b;
      b.Element-=temp;
      b.Rows+=oldr; 
      Matrix a;
      Matrix lam;
      Info.symmetricEigen(lam);
      int test=0;
      for (int i=0; i<Info.Rows;i++)if(lam(i,i)<1E-15) test=1;
      invcov.reallocate(res);
      if (test)
	invcov=1;
      else{
	a.multTranspose_(Info,b_n,2);
	invcov.multiply_(b_n,a);     
      }
    }
  return res;  
}
void Cure::MapFeature::print(int level)
{
  Matrix x;
  getX(x);
  std::cerr<<"Cure::MapFeature: "<<Key<<" ID: "<<ID<<" P-dim: "<<Bdual.Columns
	   <<"Shared Points: "<<NumberSharedPoints<<" Xf: ";
  for (int i=0; i<FullDim;i++)std::cerr<<x(i,0)<<" ";
  std::cerr<<std::endl;
  if (level>0)Bdual.print();
}

void MapFeature::get(GenericData &gd){
  MapObject::get(gd);
  int c=2*FullDim;
  if (c<8)c=8;
  if (c<gd.ShortData.Columns)c=gd.ShortData.Columns;
  gd.forceShortDataSize(3+NumberPoints, c);
  gd.ShortData(1,0)=NumberPoints;
  gd.ShortData(1,1)=Number3D;
  gd.ShortData(1,2)=Number2D;
  gd.ShortData(1,3)=NumberScalars;
  gd.ShortData(1,4)=NumberSharedPoints;
  for (int i=0; i<FullDim;i++)
    if (Index)  gd.ShortData.setLong(Index[i],2,2*i);
    else gd.ShortData(2,i)=0;
  Matrix b;
  getB(b);
  c=3*NumberPoints+NumberScalars;
  gd.Data.reallocateZero(c,Bdual.Columns+1);
  gd.Data.offset(0,1,Bdual.Rows,Bdual.Columns);
  gd.Data=(Bdual);
  gd.Data.offset(0,-1,c,Bdual.Columns+1);
  for (int i=0; i<(NumberPoints);i++){  
    if (Points[i]){
      gd.ShortData.setLong(Points[i]->Key,i+3,0);
      gd.ShortData.setLong(Points[i]->ID,i+3,2);
      gd.ShortData.setLong(Points[i]->Index,i+3,4);
      gd.ShortData.setLong(Points[i]->Local_map,i+3,6);
      for (int j=0;j<3;j++)gd.Data(j+3*i,0)=Points[i]->X[j];
    }else {
      for (int j=0;j<8;j++)
	gd.ShortData(i+3,j)=0;     
      for (int j=0;j<3;j++)gd.Data(j+3*i,0)=0;
    }
  }
  for (int i=0; i<(NumberScalars);i++)
    gd.Data(i+3*NumberPoints,0)=Scalars[i];
}
int MapFeature::set(GenericData &gd){
  if (MapObject::set(gd))return 1;
  short bankid;
  long key;
  if (getBankInfo(gd,bankid,key))return 1;
  int c=2*FullDim;
  if (c<5)c=8;
  if (gd.ShortData.Rows<3+NumberPoints)return 1;
  if (gd.ShortData.Columns<c)return 1;
  c=3*NumberPoints+NumberScalars;
  if (gd.Data.Rows<c)return 1;
  if (gd.Data.Columns<1)return 1;
  for (int i=0;i<(NumberPoints);i++){
    long k=gd.ShortData.getLong(i+3,0);
    MapObject *mo =0;
    if (Bank)mo=Bank->getMapObject(bankid,k);
    MapPoint *pt=0;
    if (mo){
      pt=mo->narrowPoint();
      if (!pt) return MAP_OBJECT_INVALID;
    }
    if (!pt)pt=new MapPoint(Bank);
    MapPoint *pt2=Points[i];
    addPoint(pt,i);
    Points[i]->ID=gd.ShortData.getLong(i+3,2);
    Points[i]->Index=gd.ShortData.getLong(i+3,4);
    Points[i]->Local_map=gd.ShortData.getLong(i+3,6);
    Points[i]->X[0]=gd.Data(3*i,0);
    Points[i]->X[1]=gd.Data(3*i+1,0);
    Points[i]->X[2]=gd.Data(3*i+2,0);
    if (Bank)Bank->associate(bankid,k,Points[i]->Key);
    if (pt2)
      if (!(pt2->Features.Next))delete pt2;
  }
  key=gd.Data.Rows;
  gd.Data.offset(0,1,FullDim,gd.Data.Columns-1);
  Bdual=(gd.Data);
  gd.Data.offset(0,-1,key,gd.Data.Columns+1);
  for (int i=0; i<FullDim;i++)
    if (Index)  Index[i]=gd.ShortData.getLong(2,2*i);
  for (int i=0; i<(NumberScalars);i++)
    Scalars[i]=gd.Data(3*NumberPoints+i,0); 
  return 0;
}
