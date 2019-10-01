//    
//
// = AUTHOR(S)
//    John Folkesson
//    
//    Copyright (c) 2007 John Folkesson
//    

#ifndef CURE_MAPDEFINES_HH
#define CURE_MAPDEFINES_HH

// Error code 
#define MAP_OBJECT_INVALID -1
#define NOT_VISABLE  -2 

//Different types of Objects:
#define LOCALMAP_TYPE  1
#define FRAME_TYPE 2
#define ENERGYNODE_TYPE 3
#define FEATURE_TYPE  4  
#define STATENODE_TYPE  5  

// Different types of features
#define MAPWALL_TYPE          1
#define MAPRANGEDPOINT_TYPE  2
#define MAPPOLE_TYPE         3
#define MAPPOINTFEATURE_TYPE 64
#define MAPHLINE_TYPE        66
#define RELATIVEPOINTFEATURE_TYPE 128


// Different types of StateNodes
#define MAPRELATIVEFEATURE_TYPE 1
#define MAPRELATIVEPOINT_TYPE   2
#define MAPMATCHEDPOINT_TYPE 3
#define POSENODE_TYPE           128


// Different types of LocalMaps
#define EKFFILTER_TYPE     1
#define TRIANGLEFILTER_TYPE     2


//Different types of EnergyNodes
#define FEATUREMEASURMENT_TYPE     1

//Different types of Frames
#define POSEDFRAME_TYPE     1
#define POSETREE_TYPE       2

//Measurements
#define CAMERA_TYPE 1
#define SCAN_TYPE 2
#define RANGEBEARING_TYPE 4


namespace Cure{


} // namespace Cure

#endif 
