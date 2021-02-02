
#ifndef  RECONSTRUCTEDPARTICLE2MC_ANALYZERS_H
#define  RECONSTRUCTEDPARTICLE2MC_ANALYZERS_H

#include <cmath>
#include <vector>

#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"
#include "podio/ObjectID.h"
#include "TLorentzVector.h"




/// Return the D0 of a track to a reconstructed particle
ROOT::VecOps::RVec<float> getRP2MC_p (ROOT::VecOps::RVec<int> recin,
				      ROOT::VecOps::RVec<int> mcin,
				      ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
				      ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);

/// Return the D0 of a track to a reconstructed particle
ROOT::VecOps::RVec<float> getRP2MC_px (ROOT::VecOps::RVec<int> recin,
				       ROOT::VecOps::RVec<int> mcin,
				       ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
				       ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);

/// Return the D0 of a track to a reconstructed particle
ROOT::VecOps::RVec<float> getRP2MC_py (ROOT::VecOps::RVec<int> recin,
				       ROOT::VecOps::RVec<int> mcin,
				       ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
				       ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);

/// Return the D0 of a track to a reconstructed particle
ROOT::VecOps::RVec<float> getRP2MC_pz (ROOT::VecOps::RVec<int> recin,
				       ROOT::VecOps::RVec<int> mcin,
				       ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
				       ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);

/// Return the D0 of a track to a reconstructed particle
ROOT::VecOps::RVec<float> getRP2MC_mass (ROOT::VecOps::RVec<int> recin,
					 ROOT::VecOps::RVec<int> mcin,
					 ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
					 ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);

/// Return the D0 of a track to a reconstructed particle
ROOT::VecOps::RVec<float> getRP2MC_charge (ROOT::VecOps::RVec<int> recin,
					   ROOT::VecOps::RVec<int> mcin,
					   ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
					   ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);

/// Return the D0 of a track to a reconstructed particle
ROOT::VecOps::RVec<float> getRP2MC_pdg (ROOT::VecOps::RVec<int> recin,
					ROOT::VecOps::RVec<int> mcin,
					ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
					ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);

/// Return the D0 of a track to a reconstructed particle
ROOT::VecOps::RVec<TLorentzVector> getRP2MC_tlv (ROOT::VecOps::RVec<int> recin,
						 ROOT::VecOps::RVec<int> mcin,
						 ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
						 ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);

/// Return the D0 of a track to a reconstructed particle
ROOT::VecOps::RVec<int> getRP2MC_index (ROOT::VecOps::RVec<int> recin,
					ROOT::VecOps::RVec<int> mcin,
					ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco);

/// Return the D0 of a track to a reconstructed particle
ROOT::VecOps::RVec<int> getRP2MC_index_test (ROOT::VecOps::RVec<int> recin,
					     ROOT::VecOps::RVec<int> mcin,
					     ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
					     ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
					     ROOT::VecOps::RVec<int> parents);


/// Return the D0 of a track to a reconstructed particle
ROOT::VecOps::RVec<int> getRP2MC_parentid (ROOT::VecOps::RVec<int> recin,
					   ROOT::VecOps::RVec<int> mcin,
					   ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
					   ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
					   ROOT::VecOps::RVec<int> parents);


/// select ReconstructedParticles with transverse momentum greater than a minimum value [GeV]
struct getRP2MC_p_func {
  ROOT::VecOps::RVec<float>  operator() (ROOT::VecOps::RVec<int> recin,
					 ROOT::VecOps::RVec<int> mcin,
					 ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
					 ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);
};

struct selMuons_JPsimatch {
  selMuons_JPsimatch( int pdg_mother, int pdg_daughter1, int pdg_daughter2) ;
  int m_pdg_mother = 443;
  int m_pdg_daughter1 = 13;
  int m_pdg_daughter2 = -13;
  std::vector<edm4hep::ReconstructedParticleData> operator() (ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> mcdaughters) ;
};

struct selRP_PDG {
  selRP_PDG(int arg_PDG, bool arg_chargedOnly);
  int m_PDG = 13 ;
  bool m_chargedOnly = true;
  std::vector<edm4hep::ReconstructedParticleData>  operator() (ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) ;
};

std::vector<edm4hep::ReconstructedParticleData> selRP_ChargedHadrons ( ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) ;

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> selRP_Bs2JPsiPhi( ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> mcdaughters) ;

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> selChargedRP_MCmatch_daughtersOf( int index, ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> ind ) ;


int getTrack2MC_index (  int track_index,
                         ROOT::VecOps::RVec<int> recind,
                         ROOT::VecOps::RVec<int> mcind,
                         ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco);

#endif

