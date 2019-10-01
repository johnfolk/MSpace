
// = AUTHOR(S)
//    John Folkesson
//    Copyright (c) 2004 John Folkesson
//    

#include "Match.hh"

void Cure::Match::print(int level)
{
  std::cerr<<" Match distance: "<<MatchDistance
	   <<" Weight "<<Weight
	   <<" PathDistance "<<PathDistance<<"\n";
  if (level&0x20)
    if (MatchedFeature)MatchedFeature->print();
  if (level&40)
   if (Measure) Measure->print(level&0x1F);
}
