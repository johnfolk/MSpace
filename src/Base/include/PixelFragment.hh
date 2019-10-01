// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_PIXELFRAGMENT_H
#define CURE_PIXELFRAGMENT_H

namespace Cure{
  
  /**
   * Helper class for visual slam
   *
   * @author John Folkesson
   * @see
   */
  class PixelFragment
  {
  public:
    int X;
    int Y;
    double Range;
    PixelFragment *Next;
  public:
    PixelFragment()
    {
      X=0;
      Y=0;
      Range=0;
      Next=0;
    }
    ~PixelFragment()
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

    PixelFragment *add(int x,int y,double r)
    {
      if (Next)return Next->add(x,y,r);
      Next=new PixelFragment();
      Range=r;
      X=x;
      Y=y;
      return this;
    }
  };

} // namespace Cure

#endif
