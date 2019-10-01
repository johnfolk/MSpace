// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "MapHLine.hh"

using namespace Cure;
using namespace MatrixStuff;

MapHLine::MapHLine():MapFeature(0)
{
  FullDim=6;
  NumberPoints=2;
  Number3D=2;
  Points=AllocatedPoints;
  Index=AllocatedIndex;
  Type=MAPHLINE_TYPE;
  VCount=0;
  TotalWgt=0;
  Tangent[0]=0;
  Tangent[1]=0;
  Tangent[2]=0;
  InfoTangent[0]=0;
  InfoTangent[1]=0;
  InfoTangent[2]=0;
  Length=0;
  init();
  CastPtr=0;
  TanThreshold=1;
  RhoThreshold=1;
  DistanceThreshold=20;
  TrackThreshold=.1;
  SqWidth=0;
}
MapHLine::MapHLine(MapHLine *wp,MapBank *b):MapFeature(b)
{
  FullDim=6;
  NumberPoints=2;
  Number3D=2;
  Points=AllocatedPoints;
  Index=AllocatedIndex;
  Type=MAPHLINE_TYPE;
  VCount=0;
  TotalWgt=0;
  Tangent[0]=0;
  Tangent[1]=0;
  Tangent[2]=0;
  InfoTangent[0]=0;
  InfoTangent[1]=0;
  InfoTangent[2]=0;
  Length=0;
  SqWidth=0;
  init();
  if (wp)
    {
      CastPtr=wp->CastPtr;
      TanThreshold=wp->TanThreshold;
      RhoThreshold=wp->RhoThreshold;
      DistanceThreshold=wp->DistanceThreshold;
      TrackThreshold=wp->TrackThreshold;
    }
  else
    {
      CastPtr=0;
      TanThreshold=1;
      RhoThreshold=1;
      DistanceThreshold=20;
      TrackThreshold=.1;
    }
}
void MapHLine::setStart(MapPoint *val)
{
  addPoint(val,0);
}
void MapHLine::setEnd(MapPoint *val)
{
  addPoint(val,1);
}
void   MapHLine::calcInfoTanRho(Transformation3D & map2info)
{
  double x[6];
  map2info.transform(Points[0]->X,x);
  map2info.transform(Points[1]->X,x+3);
  InfoTangent[0]=x[3]-x[0];
  InfoTangent[1]=x[4]-x[1];
  double d=InfoTangent[0]*InfoTangent[0]+InfoTangent[1]*InfoTangent[1];
  Length=sqrt(d);
  if (Length<1E-9)return;
  InfoTangent[0]/=Length;
  InfoTangent[1]/=Length;
  Rho[0]=(x[3]+x[0])/2;
  Rho[1]=(x[4]+x[1])/2;
  Rho[2]=(x[5]+x[2])/2;
  d=Rho[0]*InfoTangent[0]+Rho[1]*InfoTangent[1];
  Rho[0]-=InfoTangent[0]*d;
  Rho[1]-=InfoTangent[1]*d;
}

void MapHLine::calcTangent()
{
  // First calc dx and dy between line end points
  // dx
  Tangent[0]=Points[1]->X[0]-Points[0]->X[0];
  // dy
  Tangent[1]=Points[1]->X[1]-Points[0]->X[1];

  // Get distance between points
  double d=Tangent[0]*Tangent[0]+Tangent[1]*Tangent[1];
  Length=sqrt(d);
  if (Length<1E-9)return;

  // Normalize so that Tangent becomes cos(alpha) and sin(alpha)
  Tangent[0]/=Length;
  Tangent[1]/=Length;
}
void MapHLine::calcRho()
{
  Rho[0]=(Points[1]->X[0]+Points[0]->X[0])/2;
  Rho[1]=(Points[1]->X[1]+Points[0]->X[1])/2;
  Rho[2]=(Points[1]->X[2]+Points[0]->X[2])/2;
  double d=Rho[0]*Tangent[0]+Rho[1]*Tangent[1];
  Rho[0]-=Tangent[0]*d;
  Rho[1]-=Tangent[1]*d;
}

void MapHLine::initializeFromPixels(Transformation3D &t,
				    const double startpixels[2],
				    const double endpixels[2],
				    double focallength,
				    double heightguess)
{
  if (!Points[0])
    {
      MapPoint *p=new MapPoint(Bank);
      setStart(p);
    }
  if (!Points[1])
    {
      MapPoint *p=new MapPoint(Bank);
      setEnd(p);
    }
  Cure::Vector3D v;
  double y=heightguess;
  double d=(y-focallength)/focallength;
  if ((d<1E-9)&&(d>-1E-9))d=1E-9;
  v(0)=-(startpixels[0]*d);
  v(2)=-(startpixels[1]*d);
  v(1)=heightguess;
  t.invTransform(v, *Points[0]);
  // d=((*Points[0])(2)-t.Position(2));
  v(0)=-(endpixels[0]*d);
  v(2)=-(endpixels[1]*d);
  v(1)=heightguess;
  t.invTransform(v, *Points[1]);
  /*
  Cure::Vector3D v;

  double d=(1-focallength);
  v(0)=-(startpixels[0]*d);
  v(2)=-(startpixels[1]*d);
  v(1)=focallength;
  t.invTransform(v, *Points[0]);
  d=((*Points[0])(2)-t.Position(2));
  if ((d<1E-9)&&(d>-1E-9))d=1E-9;
  d=heightguess/d;
  (*Points[0])(0)+=d*((*Points[0])(0)-t.Position(0));
  (*Points[0])(1)+=d*((*Points[0])(1)-t.Position(1));
  (*Points[0])(2)=heightguess+t.Position(2);
  d=(2-focallength);
  v(0)=-(endpixels[0]*d);
  v(2)=-(endpixels[1]*d);
  v(1)=focallength;
  t.invTransform(v, *Points[1]);
  d=((*Points[1])(2)-t.Position(2));
  if ((d<1E-9)&&(d>-1E-9))d=1E-9;    
  d=heightguess/d;
  (*Points[1])(0)+=d*((*Points[1])(0)-t.Position(0));
  (*Points[1])(1)+=d*((*Points[1])(1)-t.Position(1));
  (*Points[1])(2)=heightguess+t.Position(2);
  */
  recenter();
  calcRho();
}

int  MapHLine::extend()
{
  if (Bdual.Columns>2)return 0;
  if (Bdual.Columns==0)
    {
      if (TotalWgt>TanThreshold)
	{
	  Bdual.reallocate(6,1);
	  recenter();
	  return 1;
	}
    }
  else if (Bdual.Columns==1)
    {
      if (TotalWgt>RhoThreshold)
	{
	  Bdual.reallocate(6,3);
	  recenter();
	  return 2;
	}
    }
  return 0;  
}
void MapHLine::lengthAdjust()
{
  if (Bdual.Columns==0)return;
  double d=Bdual(0,0)*Bdual(0,0)+Bdual(1,0)*Bdual(1,0);
  d*=2;
  d=Length/sqrt(d);
  Bdual(0,0)*=d;
  Bdual(1,0)*=d;
  Bdual(3,0)=-Bdual(0,0);
  Bdual(4,0)=-Bdual(1,0);
}
unsigned short MapHLine::getMeasurementType(unsigned short type)
{
  if (!(type<8))type=0;
  if (!type)return type;
  if (Bdual.Columns==3) type=2;
  else if (Bdual.Columns==1) type=1;
  else type=0;
  return type;  
}

int MapHLine::estimateRho()
{
  double d=0;
  int i=0;
  int k=0;
  double wgt[VCount];
  int top=3*VCount;
  double n[top];
  double a[top];
  for (LinkedArray *vlist=&Vectors;vlist->Next;vlist=vlist->Next,i+=3,k++)
    {	
      wgt[k]=(*vlist)(1);
      a[i]=(*vlist)(2);
      a[i+1]=(*vlist)(3);
      a[i+2]=(*vlist)(4);
      n[i]=(*vlist)(12)+(*vlist)(20);
      n[i+1]=(*vlist)(13)+(*vlist)(21);
      n[i+2]=(*vlist)(14)+(*vlist)(22);
      d=n[i]*InfoTangent[0]+n[i+1]*InfoTangent[1];
      n[i]-=d*InfoTangent[0];
      n[i+1]-=d*InfoTangent[1];
      d=sqrt(n[i]*n[i]+n[i+1]*n[i+1]+n[i+2]*n[i+2]);
      n[i]/=d;
      n[i+1]/=d;
      n[i+2]/=d;
      d=a[i]*InfoTangent[0]+a[i+1]*InfoTangent[1];
      a[i]-=d*InfoTangent[0];
      a[i+1]-=d*InfoTangent[1];
    }
  double sw=0;
  double maxw=0;
  double temp[3];
  temp[0]=0;
  temp[1]=0;
  temp[2]=0;
  k=0;
  for (i=0;i<top; i+=3,k++)
    {
      int h=k+1; 
      for (int j=i+3; j<top;j+=3,h++)
	{
	  double dij[3];
	  dij[0]=a[j]-a[i];
	  dij[1]=a[j+1]-a[i+1];
	  dij[2]=a[j+2]-a[i+2];
	  d=sqrt(dij[0]*dij[0]+dij[1]*dij[1]+dij[2]*dij[2]);
	  if (d>.05)
	    {
	      dij[0]/=d;
	      dij[1]/=d;
	      dij[2]/=d;
	      double trig[4];
	      trig[0]=n[i]*dij[0]+n[i+1]*dij[1]+n[i+2]*dij[2];
	      trig[2]=-n[j]*dij[0]-n[j+1]*dij[1]-n[j+2]*dij[2];
	      trig[1]=sqrt(1-trig[0]*trig[0]);
	      trig[3]=sqrt(1-trig[2]*trig[2]);
	      double t=trig[0]*trig[3]+trig[2]*trig[1];
	      double w=t;
	      w*=(wgt[k]*wgt[h]);
	      d*=w;
	      trig[3]*=d;
	      trig[1]*=d;
	      temp[0]+=(n[i]*trig[3]+n[j]*trig[1]);
	      temp[1]+=(n[i+1]*trig[3]+n[j+1]*trig[1]);
	      temp[2]+=(n[i+2]*trig[3]+n[j+2]*trig[1]);
	      w*=t;
	      temp[0]+=w*(a[i]+a[j]);
	      temp[1]+=w*(a[i+1]+a[j+1]);
	      temp[2]+=w*(a[i+2]+a[j+2]);
	      sw+=w;
	      if (maxw<w)maxw=w;
	    }
	}
    }
  maxw*=64;
  sw*=2;
  if (sw<1E-35)return 1;
  Rho[0]=temp[0]/sw;
  Rho[1]=temp[1]/sw;
  Rho[2]=temp[2]/sw;
  d=(Rho[0]*InfoTangent[0]+Rho[1]*InfoTangent[1]);
  Rho[0]-=d*InfoTangent[0];
  Rho[1]-=d*InfoTangent[1];
  if ((sw<RhoThreshold)||(maxw<RhoThreshold))
    {
      // calcRho();
      return 1;
    }
  for (LinkedArray *vlist=&Vectors;vlist->Next;vlist=vlist->Next)
    fit_it(vlist->Element);
  return 0;
}

int MapHLine::getEllipsiod(Matrix &cov, double wsums[4])
{
  Matrix cloud(2*VCount,4);
  int j=0;
  for (LinkedArray *vlist=&Vectors;vlist->Next;vlist=vlist->Next,j+=2)
    {	
      fit_it(vlist->Element);
      cloud(j,0)=(*vlist)(2)+(*vlist)(19)*(*vlist)(12);
      cloud(j,1)=(*vlist)(3)+(*vlist)(19)*(*vlist)(13);
      cloud(j,2)=(*vlist)(4)+(*vlist)(19)*(*vlist)(14);
      cloud(j,3)=(*vlist)(1);
      cloud(j+1,0)=(*vlist)(2)+(*vlist)(27)*(*vlist)(20);
      cloud(j+1,1)=(*vlist)(3)+(*vlist)(27)*(*vlist)(21);
      cloud(j+1,2)=(*vlist)(4)+(*vlist)(27)*(*vlist)(22);
      cloud(j+1,3)=(*vlist)(1);
    }
  wsums[0]=0;
  wsums[1]=0;
  wsums[2]=0;
  wsums[3]=0;
  cov=0;
  for (int i=0; i<cloud.Rows; i++)
    {
      if (cloud(i,3))
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
	  cov(2,2)+=cloud(i,2)*cloud(i,2)*cloud(i,3);
	}
    }
  if (wsums[3]<1E-127) return 1;
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
  return 0;
}
int   MapHLine::eigen_it(double ev[9], double lambda[3],  double wsums[4])
{
  Matrix cov(3);
  getEllipsiod(cov,wsums);
  Matrix lam,evec;
  int res=cov.symmetricEigen(lam,evec);
  for (int i=0; i<3; i++)
    {
      lambda[i]=lam(i,i);
      for (int j=0; j<3; j++)
	  ev[3*i+j]=evec(j,i);
    }  
  if (res) return res;
  lambda[0]+=1E-4;
  lambda[1]+=1E-4;
  lambda[2]+=1E-4;
  
  if (VCount>2)
    {
      double d=(((double)VCount*1.5)/(((double)VCount)-2.5));
      lambda[0]*=d;
      lambda[1]*=d;
      lambda[2]*=d;
    }
  else
    {
      lambda[0]*=100;
      lambda[1]*=100;
      lambda[2]*=100;
    }
  double d=ev[0]*InfoTangent[0]+ev[1]*InfoTangent[1];
  if (lambda[1]<1E-4)lambda[1]=1E-4;
  if (lambda[2]<1E-4)lambda[2]=1E-4;
  if (d>.5)
    {
      d=sqrt(ev[0]*ev[0]+ev[1]*ev[1]);
      InfoTangent[0]=ev[0]/d;
      InfoTangent[1]=ev[1]/d;
    }
  else if (d<-.5)
    {
      d=sqrt(ev[0]*ev[0]+ev[1]*ev[1]);
      InfoTangent[0]=-ev[0]/d;
      InfoTangent[1]=-ev[1]/d;
    }
  else  return 1;
  lambda[0]*=2;  //give half of info to each end
  lambda[1]*=2; 
  lambda[2]*=2;
  d=wsums[0]*InfoTangent[0]+wsums[1]*InfoTangent[1];
  if(estimateRho())
    {
      Rho[0]=wsums[0]-d*InfoTangent[0];
      Rho[1]=wsums[1]-d*InfoTangent[1];
      Rho[2]=wsums[2];
    }
  return 0;
}
    
/**
 * distance, w, x_i, T*x_i, c_i, t_i, (0..11)
 * s_i, s_i*x_i, T*s_i, Rho*s_i,  r_s, q_s (12..19)
 * e_i, e_i*x_i, T*e_i, Rho*e_i, r_e, q_e (20..27)
 * t_i=sum c_i x c_j for j>i
 */
void MapHLine::fit_it(double v[28])
{
  v[5]=InfoTangent[0]*v[2]+InfoTangent[1]*v[3];
  v[16]=InfoTangent[0]*v[12]+InfoTangent[1]*v[13];
  v[17]=Rho[0]*v[12]+Rho[1]*v[13]+Rho[2]*v[14];
  v[24]=InfoTangent[0]*v[20]+InfoTangent[1]*v[21];
  v[25]=Rho[0]*v[20]+Rho[1]*v[21]+Rho[2]*v[22];
  double temp=(1-v[16]*v[16]);
  if (temp<1E-35)temp=1E35;
  else temp=1/temp;
  v[18]=temp*(v[5]-v[16]*(v[15]-v[17]));
  v[19]=temp*(v[17]-v[15]+(v[16]*v[5]));
  temp=(1-v[24]*v[24]);
  if (temp<1E-35)temp=1E35;
  else temp=1/temp;
  v[26]=temp*(v[5]-v[24]*(v[23]-v[25]));
  v[27]=temp*(v[25]-v[23]+(v[24]*v[5]));
}
void MapHLine::prune(double minDistance)
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

int  MapHLine::addInfo(const Cure::Matrix & v, Transformation3D &map2info,
		       const int type)
{
  if (Bdual.Columns==6)return 0;
  if ( LastDistance<=v(0,0))
    {
      LastDistance=v(0,0);
      double minDistance=(v(0,0)-DistanceThreshold); 
      prune(minDistance);
      calcInfoTanRho(map2info);  
      Matrix m(1,28);
      makeInfo(v, m.Element);
      map2info.invRotate(m.Element+6,LastBearing);
      VCount++;
      Vectors.add(m);
      if (Bdual.Columns==0)
	{
	  double wtan[2];
	  wtan[0]=0;
	  wtan[1]=0;
	  for (LinkedArray *vlist=&Vectors;vlist->Next;vlist=vlist->Next)
	    {
	      wtan[0]+=(*vlist)(9);
	      wtan[1]+=(*vlist)(10);
	    }
	  TotalWgt=sqrt(wtan[0]*wtan[0]+wtan[1]*wtan[1]);
	  if (TotalWgt>TanThreshold)
	    {
	      InfoTangent[0]=wtan[0]/TotalWgt;
	      InfoTangent[1]=wtan[1]/TotalWgt;
	      double ev[9],lambda[3];
	      int tst=10;
	      double quality=1E32;
	      double oldtan[2],oldRho[3];
	      Info=1;
	      while (tst>0)
		{
		  oldRho[0]=Rho[0];
		  oldRho[1]=Rho[1];
		  oldRho[2]=Rho[2];
		  oldtan[0]=InfoTangent[0];
		  oldtan[1]=InfoTangent[1];
		  double wsums[4];
		  if (!(eigen_it(ev, lambda,wsums)))
		    {
		      double d=lambda[1]*lambda[2];
		      if (d<quality)
			{
			  quality=d-1E-11; 	  
			  map2info.invRotate(ev,ev);
			  map2info.invRotate(ev+3,ev+3);
			  map2info.invRotate(ev+6,ev+6);
			  for (int i=0; i<2; i++)
			    for (int j=0;j<2;j++)
			      {
				Info(i,j)=ev[3+i]*ev[3+j]/lambda[1]+
				  ev[6+i]*ev[6+j]/lambda[2];
				Info(i+3,j+3)=Info(i,j);
			      }
			  Info(2,5)=Info(2,2);
			  Info(5,2)=Info(2,5);
			  oldRho[0]=Rho[0];
			  oldRho[1]=Rho[1];
			  oldRho[2]=Rho[2];
			} 
		      else
			tst=-1;
		    }
		  else
		    {
		      TotalWgt=0;
		      InfoTangent[0]=oldtan[0];
		      InfoTangent[1]=oldtan[1];
		      Rho[0]=oldRho[0];
		      Rho[1]=oldRho[1];
		      Rho[2]=oldRho[2];
		      trackLine(m,map2info);	
		      return 0;
		    }
		  tst--;
		}
	      if (tst<0)
		{
		  InfoTangent[0]=oldtan[0];
		  InfoTangent[1]=oldtan[1];
		  Rho[0]=oldRho[0];
		  Rho[1]=oldRho[1];
		  Rho[2]=oldRho[2];
		}
	      fit_it(m.Element);
	      trackLine(m,map2info);
	      return 0;
	    }
	  trackLine(m,map2info);
	  return 0;  
	}
      if (Bdual.Columns==1)
	{

	  double wtan[2];
	  wtan[0]=0;
	  wtan[1]=0;
	  for (LinkedArray *vlist=&Vectors;vlist->Next;vlist=vlist->Next)
	    {
	      wtan[0]+=(*vlist)(9);
	      wtan[1]+=(*vlist)(10);
	    }
	  TotalWgt=sqrt(wtan[0]*wtan[0]+wtan[1]*wtan[1]);
	  if (TotalWgt>RhoThreshold)
	    {
	      
	      if(estimateRho())
		{
		  TotalWgt=0;
		  trackLine(m,map2info);
		  return 0;
		}
	      int tst=10;
	      double quality=1E32;
	      double oldtan[2],oldRho[3];

	      Info=1;
	  
	      while (tst>0)
		{
		  oldRho[0]=Rho[0];
		  oldRho[1]=Rho[1];
		  oldRho[2]=Rho[2];
		  oldtan[0]=InfoTangent[0];
		  oldtan[1]=InfoTangent[1];
		  double ev[9],lambda[3];
		  double wsums[4];
		  if (!(eigen_it(ev, lambda,wsums)))
		    {
		      double d=lambda[1]*lambda[2];
		      if (d<quality)
			{
			  double dot=wsums[0]*oldtan[0]+wsums[1]*oldtan[1];
			  oldRho[0]=wsums[0]-oldtan[0]*dot;
			  oldRho[1]=wsums[1]-oldtan[1]*dot;
			  oldRho[2]=wsums[2];
			  quality=d-1E-14; 	     
			  map2info.invRotate(ev,ev);
			  map2info.invRotate(ev+3,ev+3);
			  map2info.invRotate(ev+6,ev+6);
			  for (int i=0; i<3; i++)
			    for (int j=0;j<3;j++)
			      {
				Info(i,j)=ev[3+i]*ev[3+j]/lambda[1]+
				  ev[6+i]*ev[6+j]/lambda[2];
				Info(i+3,j+3)=Info(i,j);
			      }
			} 
		      else
			tst=-1;
		    }
		  else
		    {
		      TotalWgt=0;
		      InfoTangent[0]=oldtan[0];
		      InfoTangent[1]=oldtan[1];
		      Rho[0]=oldRho[0];
		      Rho[1]=oldRho[1];
		      Rho[2]=oldRho[2];
		      map2info.invRotate(InfoTangent,Tangent);
		      calcRho();
		      lengthAdjust();
		      return 0;
		    }
		  tst--;
		}
	      InfoTangent[0]=oldtan[0];
	      InfoTangent[1]=oldtan[1];
	      //	  if (tst<0)
	      {
		Rho[0]=oldRho[0];
		Rho[1]=oldRho[1];
		Rho[2]=oldRho[2];
	      }
	      Matrix x;
	      getX(x);
	      map2info.transform(x.Element,x.Element);
	      map2info.transform(x.Element+3,x.Element+3);
	      double t=(x(0,0)-Rho[0])*InfoTangent[0]+
		(x(1,0)-Rho[1])*InfoTangent[1];
	      x(0,0)=Rho[0]+InfoTangent[0]*t;
	      x(1,0)=Rho[1]+InfoTangent[1]*t;
	      x(2,0)=Rho[2];

	      t=(x(3,0)-Rho[0])*InfoTangent[0]+(x(4,0)-Rho[1])*InfoTangent[1];
	      x(3,0)=Rho[0]+InfoTangent[0]*t;
	      x(4,0)=Rho[1]+InfoTangent[1]*t;
	      x(5,0)=Rho[2];
	      map2info.invTransform(x.Element,x.Element);
	      map2info.invTransform(x.Element+3,x.Element+3);
	      setX(x);
	      setEnds(map2info);
	      return 0;
	    }
	  trackLine(m,map2info);
	  return 0;  
	}
      double min=m(0,18);
      double max=m(0,26);
      Matrix x;
      getX(x);
      map2info.transform(x.Element,x.Element);
      map2info.transform(x.Element+3,x.Element+3);
      double dl[2];
      dl[0]=0;
      dl[1]=0;
      double d=x(0,0)*InfoTangent[0]+x(1,0)*InfoTangent[1];
      if (min<d)dl[0]=min-d;
      d=x(3,0)*InfoTangent[0]+x(4,0)*InfoTangent[1];
      if (max>d)dl[1]=max-d; 
      if ((dl[0]<0)||(dl[1]>0))
	{
	  d=((x(0,3)-x(0,0))*(x(3,0)-x(0,0))+(x(4,0)-x(1,0))*(x(4,0)-x(1,0)));
	  d=sqrt(d);
	  d=(d-dl[0]+dl[1])/d;
	  x(0,0)+=dl[0]*InfoTangent[0];
	  x(1,0)+=dl[0]*InfoTangent[1];
	  x(3,0)+=dl[1]*InfoTangent[0];
	  x(4,0)+=dl[1]*InfoTangent[1];
	  map2info.invTransform(x.Element,x.Element);
	  map2info.invTransform(x.Element+3,x.Element+3);      
	  setX(x);
	  lengthAdjust();
	  return 0;
	}
    }
  else
    {
      LinkedArray *la=&Vectors;
      while(la->Next)
	{
	  if((*la)(0)!=v(0,0))la=la->Next;
	  else break;
	}
      if (la->Next)
	{
	  if ((*la)(0)!=v(0,0))return 1;
	}
      else return 1;
      calcInfoTanRho(map2info);  
      makeInfo(v, la->Element);
    }
  return 0;
}
void MapHLine::makeInfo(const Matrix &v, double m[28])
{
  m[0]=v(0,0);
  m[1]=v(0,1);
  m[2]=v(0,2);
  m[3]=v(0,3);
  m[4]=v(0,4);
  // distance, w, x_i, T*x_i, c_i, t_i, (0..11)
  //s_i, s_i*x_i, T*s_i, Rho*s_i,  r_s, q_s (12..19)
  //e_i, e_i*x_i, T*e_i, Rho*e_i, r_e, q_e (20..27)
  //t_i weighted tangent 
  double d=sqrt(v(0,7)*v(0,7)+v(0,5)*v(0,5)+v(0,6)*v(0,6));
  m[12]=v(0,5)/d;
  m[13]=v(0,6)/d;;
  m[14]=v(0,7)/d;
  d=sqrt(v(0,10)*v(0,10)+v(0,8)*v(0,8)+v(0,9)*v(0,9));
  m[20]=v(0,8)/d;
  m[21]=v(0,9)/d;
  m[22]=v(0,10)/d;
  m[15]=m[12]*m[2]+m[13]*m[3]+m[14]*m[4];
  m[23]=m[20]*m[2]+m[21]*m[3]+m[22]*m[4];
  m[6]=m[13]*m[22]-m[14]*m[21];
  m[7]=m[14]*m[20]-m[12]*m[22];
  m[8]=m[12]*m[21]-m[13]*m[20];
  m[9]=m[7]*m[1];
  m[10]=-m[6]*m[1];
  m[11]=0;
  d=m[9]*InfoTangent[0]+m[10]*InfoTangent[1];
  if (d<0)
    {
      m[9]=-m[9];
      m[10]=-m[10];
    }
  d=sqrt(m[6]*m[6]+m[7]*m[7]+m[8]*m[8]);
  m[6]/=d;
  m[7]/=d;
  m[8]/=d;
  fit_it(m);
  if (m[18]>m[26])
    {
      for (int i=0; i<8;i++){
	double d=m[12+i];
	m[12+i]=m[20+i];
	m[20+i]=d;
      }
      m[6]=-m[6];
      m[7]=-m[7];
      m[8]=-m[8];
    }
}

void MapHLine::trackLine(Matrix &m, Transformation3D & map2info)
{
  Matrix x(6,1);
  if (Bdual.Columns==0)
    {
      x(0,0)=m(0,2)+m(0,12)*m(0,19);
      x(1,0)=m(0,3)+m(0,13)*m(0,19);
      x(2,0)=m(0,4)+m(0,14)*m(0,19);
      x(3,0)=m(0,2)+m(0,20)*m(0,27);
      x(4,0)=m(0,3)+m(0,21)*m(0,27);
      x(5,0)=m(0,4)+m(0,22)*m(0,27);
      x(2,0)+=x(5,0);
      x(2,0)/=2;
      x(5,0)=x(2,0);
      map2info.invTransform(x.Element,x.Element);
      map2info.invTransform(x.Element+3,x.Element+3);
      setX(x);
      setEnds(map2info);
      return;
    }
  getX(x);
  map2info.transform(x.Element,x.Element);
  map2info.transform(x.Element+3,x.Element+3);
  x(2,0)=m(0,4)+m(0,14)*m(0,19);
  x(5,0)=m(0,4)+m(0,22)*m(0,27);
  x(2,0)+=x(5,0);
  x(2,0)/=2;
  x(5,0)=x(2,0);
  double dx[2];
  dx[0]=m(0,2)+m(0,12)*m(0,19)-x(0,0);
  dx[1]=m(0,3)+m(0,13)*m(0,19)-x(1,0);
  dx[0]+=(m(0,2)+m(0,20)*m(0,27)-x(3,0));
  dx[1]+=(m(0,3)+m(0,21)*m(0,27)-x(4,0));
  double d=-dx[0]*InfoTangent[1]+dx[1]*InfoTangent[0];
  d/=2;
  x(0,0)-=d*InfoTangent[1];
  x(1,0)+=d*InfoTangent[0];
  x(3,0)-=d*InfoTangent[1];
  x(4,0)+=d*InfoTangent[0];
  map2info.invTransform(x.Element,x.Element);
  map2info.invTransform(x.Element+3,x.Element+3);
  setX(x);
  setEnds(map2info);
  return;

}
void MapHLine::setEnds(Transformation3D & map2info)
{
  calcTangent();
  Matrix x;
  getX(x);
  map2info.transform(x.Element,x.Element);
  map2info.transform(x.Element+3,x.Element+3);
   map2info.rotate(Tangent,InfoTangent);
  double min=x(0,0)*InfoTangent[0]+x(1,0)*InfoTangent[1];
  double max=x(3,0)*InfoTangent[0]+x(4,0)*InfoTangent[1];
  double oldRho[3];
  oldRho[0]=Rho[0];
  oldRho[1]=Rho[1];
  oldRho[2]=Rho[2];
  Rho[0]=(x.Element[3]+x.Element[0])/2;
  Rho[1]=(x.Element[4]+x.Element[1])/2;
  Rho[2]=(x.Element[5]+x.Element[2])/2;
  double tmp=Rho[0]*InfoTangent[0]+Rho[1]*InfoTangent[1];
  Rho[0]-=InfoTangent[0]*tmp;
  Rho[1]-=InfoTangent[1]*tmp;
  for (LinkedArray *vlist=&Vectors;vlist->Next;vlist=vlist->Next)
    {
      fit_it(vlist->Element);
      if ( (*vlist)(18)<min)min=(*vlist)(18);
      if ( (*vlist)(26)>max)max=(*vlist)(26);
     }
  double dl[2];
  dl[0]=0;
  dl[1]=0;
  double d=x(0,0)*InfoTangent[0]+x(1,0)*InfoTangent[1];
  if (min<d)dl[0]=min-d;
  d=x(3,0)*InfoTangent[0]+x(4,0)*InfoTangent[1];
  if (max>d)dl[1]=max-d; 
  if ((dl[0]<0)||(dl[1]>0))
    {
      d=((x(0,3)-x(0,0))*(x(3,0)-x(0,0))+(x(4,0)-x(1,0))*(x(4,0)-x(1,0)));
      d=sqrt(d);
      d=(d-dl[0]+dl[1])/d;
      x(0,0)+=dl[0]*InfoTangent[0];
      x(1,0)+=dl[0]*InfoTangent[1];
      x(3,0)+=dl[1]*InfoTangent[0];
      x(4,0)+=dl[1]*InfoTangent[1];
      map2info.invTransform(x.Element,x.Element);
      map2info.invTransform(x.Element+3,x.Element+3);
      setX(x);
      lengthAdjust();
    }
  if(Bdual.Columns<3)
    {
      Rho[0]=oldRho[0];
      Rho[1]=oldRho[1];
      Rho[2]=oldRho[2];
    }
}
int MapHLine::merge(MapHLine *mf, unsigned short typ)
{
  if (!mf)return -1;
  calcTangent();
  Matrix x(6,1);
  getX(x);
  Matrix nx(6,1);
  mf->getX(nx);
  double t=-(x(0,0)*Tangent[0]+x(1,0)*Tangent[1]);
  t+=(nx(0,0)*Tangent[0]+nx(1,0)*Tangent[1]);
  if (t<0){
    x(0,0)+=t*Tangent[0];
    x(1,0)+=t*Tangent[1];
  }
  t=-(x(3,0)*Tangent[0]+x(4,0)*Tangent[1]);
  t+=(nx(3,0)*Tangent[0]+nx(4,0)*Tangent[1]);
  if (t>0){
    x(3,0)+=t*Tangent[0];
    x(4,0)+=t*Tangent[1];
  }
  setX(x);
  return 0;
}
void MapHLine::recenter()
{
  if (Bdual.Columns==0)return;

  // Update the Tangent[0] and Tangent[1] variables
  calcTangent();

  // Assuming that the tangent direction is alpha the normal angle
  // gamma = (alpha - pi/2)
  // cos(gamma) = cos(alpha-pi/2) =-sin(alpha) =-Tangent[1]
  // sin(gamma) = sin(alpha-pi/2) = cos(alpha) = Tangent[0]
  Bdual(0,0)=-Tangent[1]*M_SQRT1_2; // cos(norm_ang) / sqrt(2)
  Bdual(1,0)=Tangent[0]*M_SQRT1_2;  // sin(norm_ang) / sqrt(2)
  if (Bdual.Columns==3)
    {
      Bdual(0,1)=Bdual(0,0);
      Bdual(1,1)=Bdual(1,0);
      Bdual(2,1)=0;
      Bdual(3,1)=Bdual(0,0);
      Bdual(4,1)=Bdual(1,0);
      Bdual(5,1)=0;
      Bdual(0,2)=0;
      Bdual(1,2)=0;
      Bdual(2,2)=M_SQRT1_2;
      Bdual(3,2)=0;
      Bdual(4,2)=0;
      Bdual(5,2)=M_SQRT1_2;
    }
  Bdual(0,0)*=Length;
  Bdual(1,0)*=Length;
  Bdual(2,0)=0;
  Bdual(3,0)=-Bdual(0,0);
  Bdual(4,0)=-Bdual(1,0);
  Bdual(5,0)=0;
}
unsigned short MapHLine::getScales(Cure::Matrix &scales)
{
  if (Bdual.Rows==0){
    scales.reallocate(0,FullDim);
    return 0;
  }
  calcTangent();
  scales.reallocate(1,FullDim);
  scales(0,0)=-Tangent[0];//*M_SQRT1_2;
  scales(0,1)=-Tangent[1];//*M_SQRT1_2;
  scales(0,2)=0;
  scales(0,3)=-scales(0,0);
  scales(0,4)=-scales(0,1);
  scales(0,5)=0;
  return 1;
}

void MapHLine::getB(Cure::Matrix &b)const
{
  b.transpose_(Bdual);
  if (b.Rows==0)return;
  double d=Bdual(0,0)*Bdual(0,0)+Bdual(1,0)*Bdual(1,0);
  d*=2;
  b(0,0)/=d;
  b(0,1)/=d;
  b(0,3)/=d;
  b(0,4)/=d;
}
bool MapHLine::getC(MapHLine *mw,Cure::Matrix &c)const
{
  if ((Bdual.Columns<1)||(mw->Bdual.Columns<1))return false;
  int dim=1;
  if ((Bdual.Columns==3)&&(mw->Bdual.Columns==3))dim=3;
  c.reallocate(dim,Bdual.Columns);
  c=0;
  for (int i=0; i<dim; i++)c(i,i)=1;
  return true; 
}
bool MapHLine::testMatch(MapHLine *mw, double tolerance)
{
  Matrix x1,x2;
  getX(x1);
  mw->getX(x2);
  Matrix dx=x1;
  calcTangent();
  mw->calcTangent();
  double d=Tangent[0]*mw->Tangent[1]-Tangent[0]*mw->Tangent[1];
  if (d<0)d=-d;
  if (d>M_PI_2){
    x2.offset(3,0,3,1);
    dx.offset(0,0,3,1);
    dx-=x2;
    dx.offset(3,0,3,1);
    x2.offset(-3,0,3,1);
    dx-=x2;
    dx.reset(6,1);
    d-=M_PI;
  }
  if (d<0)d=-d;
  else dx-=x2;
  if ((Bdual.Columns<1)||(mw->Bdual.Columns<1)){
    if (Length>mw->Length)d*=mw->Length;
    else d*=mw->Length;
    if (tolerance<d)  return false;
  }
  if ((Bdual.Columns==1)||(mw->Bdual.Columns==1)){
    d=(dx(0,0)+dx(3,0))*Tangent[1]-(dx(1,0)+dx(4,0))*Tangent[0];
    if (d<0)d=-d;
    d/=2;
    if (tolerance<d)  return false;
    d=dx(2,0)/10;
    if (tolerance<d)  return false;
  }
  tolerance*=2;
  d=dx(0,0)*Tangent[0]+dx(1,0)*Tangent[1];
  if (d<-(mw->Length+tolerance))return false;
  if (-d>(Length+tolerance))return false;
  d=dx(3,0)*Tangent[0]+dx(4,0)*Tangent[1];
  if (-d>(mw->Length+tolerance))return false;
  if (d>(Length+tolerance))return false;
  return true;  

}

void MapHLine::print(int level)
{
  Matrix x;
  getX(x);
  std::cerr<<"MapHLine: "<<Key<<" "<<ID<<" "<<Bdual.Columns<<" "<<VCount<<" ";
  for (int i=0; i<FullDim;i++)std::cerr<<x(i,0)<<" ";
  std::cerr<<std::endl;
  if (level>0)Bdual.print();
}

int MapHLine::set(GenericData &gd)
{
  if (MapFeature::set(gd))return 1;
  VCount=0;
  TotalWgt=0;
  Tangent[0]=0;
  Tangent[1]=0;
  Tangent[2]=0;
  InfoTangent[0]=0;
  InfoTangent[1]=0;
  InfoTangent[2]=0;
  Length=0;
  TanThreshold=1;
  RhoThreshold=1;
  DistanceThreshold=20;
  TrackThreshold=.1;
  SqWidth=0;
  return 0;
}
