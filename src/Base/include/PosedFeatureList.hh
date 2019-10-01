// = AUTHOR(S)
//    John Folkesson
//    
//    March 11, 2004
//
//    Copyright (c) 2004 John Folkesson
//    
#ifndef CURE_POSEDFEATURELIST_H
#define CURE_POSEDFEATURELIST_H

#include "PosedFeature.hh"
namespace Cure{

  /**
   * List for posed features
   *
   * @author John Folkesson
   * @see
   */
  class PosedFeatureList
  {
  public:
    int Count;
    PosedFeature * Element; 
    PosedFeatureList *Next;
    PosedFeatureList *Prev;  
  public:
    PosedFeatureList();
    ~PosedFeatureList();
    /*
     * Deletes all the links but not the PosedFeatures (Element)
     */
    void clean();
    /**
     * This is clean but it also deletes Element
     */
    void clear();
    int add(PosedFeature *mf);
    int addUnique(PosedFeature *mf);
    int remove(long k);
    int removeFeature(PosedFeature *mf);
    PosedFeature * get(long n);
    PosedFeature *  keyGet(long key);
    void listType(int typ, PosedFeatureList & result);
  };

} // namespace Cure

#endif
