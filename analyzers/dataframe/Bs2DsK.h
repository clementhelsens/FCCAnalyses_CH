
#ifndef BS2DSK_ANALYZERS_H
#define BS2DSK_ANALYZERS_H

#include <cmath>
#include <vector>

#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"
#include "podio/ObjectID.h"
#include "TLorentzVector.h"


#include "MCParticle.h"
#include "ReconstructedParticle2MC.h"
#include "Vertex.h"

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  selChargedRP_KfromBs ( ROOT::VecOps::RVec<int> RP2MC_indices, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec< std::vector<int> > ndecays  );

ROOT::VecOps::RVec< ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> > selChargedRP_DsfromBs( ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int>  ind, ROOT::VecOps::RVec< std::vector<int> > ndecays  ) ;

ROOT::VecOps::RVec< int > n_DsTracks( ROOT::VecOps::RVec< ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> > all_Ds_tracks );

ROOT::VecOps::RVec< FCCAnalysesVertex >  DsVertex( ROOT::VecOps::RVec< ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> > all_Ds_tracks,
                                                     ROOT::VecOps::RVec<edm4hep::TrackState> tracks  ) ;

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  ReconstructedDs( ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> ind, ROOT::VecOps::RVec< std::vector<int> > ndecays  ) ;

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  ReconstructedDs_atVertex( ROOT::VecOps::RVec<FCCAnalysesVertex> DsVertices,
                ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind,
                ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) ;

ROOT::VecOps::RVec<edm4hep::TrackState>  ReconstructedDs_atVertex_TrackState( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> theDss)  ;

ROOT::VecOps::RVec<  ROOT::VecOps::RVec<edm4hep::TrackState> > tracks_for_fitting_the_Bs_vertex
                        ( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> RecoDs_atVertex,
                          ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  BachelorK,
                          ROOT::VecOps::RVec<edm4hep::TrackState> all_tracks  ) ;

ROOT::VecOps::RVec< FCCAnalysesVertex >  BsVertex( ROOT::VecOps::RVec< ROOT::VecOps::RVec<edm4hep::TrackState> > all_Bs_tracks,
                                                   ROOT::VecOps::RVec<edm4hep::TrackState> tracks  ) ;

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> RecoPiplus( ROOT::VecOps::RVec< ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> > thelegs) ;
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> RecoKplus( ROOT::VecOps::RVec< ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> > thelegs) ;
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> RecoKminus( ROOT::VecOps::RVec< ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> > thelegs) ;


#endif
