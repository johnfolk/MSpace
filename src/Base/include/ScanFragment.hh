// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef SCANFRAGMENT_H
#define SCANFRAGMENT_H
#define MAXSCAN 10000


/**
 * A class that hold the data from a single scan point and can be used
 * to form lists of connected scan points
 *
 * @author John Folkesson
 * @see
 */
class ScanFragment
{
public:
  int Index;
  double Range;
  ScanFragment *Next;
public:
  ScanFragment()
    {
      Index=0;
      Range=0;
      Next=0;
    }
  ~ScanFragment()
    {
      clean();
    }
  void clean()
    {
      if (Next)
	{
	  delete Next;
	  Next=0;
	}
    }
  ScanFragment *add(int indx,double r)
    {
      if (Next)return Next->add(indx,r);
      Next=new ScanFragment();
      Range=r;
      Index=indx;
      return this;
    }
};
#endif
