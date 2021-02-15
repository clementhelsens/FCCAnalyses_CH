#ifndef  VERTEXINGACTS_ANALYZERS_H
#define  VERTEXINGACTS_ANALYZERS_H

#include <cmath>
#include <vector>
#include "edm4hep/TrackState.h"
#include "ROOT/RVec.hxx"



namespace VertexingACTS{

bool initialize(ROOT::VecOps::RVec<edm4hep::TrackState> tracks);

}

#endif