#include "ReconstructedParticle2MC.h"

#include "MCParticle.h"




ROOT::VecOps::RVec<float> getRP2MC_p(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<float> result;
  result.resize(reco.size(),-1.);

  for (unsigned int i=0; i<recind.size();i++) {
    TLorentzVector tlv;
    tlv.SetXYZM(mc.at(mcind.at(i)).momentum.x,mc.at(mcind.at(i)).momentum.y,mc.at(mcind.at(i)).momentum.z,mc.at(mcind.at(i)).mass);
    result[recind.at(i)]=tlv.P();
  }
  return result;
}

ROOT::VecOps::RVec<TLorentzVector> getRP2MC_tlv(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<TLorentzVector> result;
  result.resize(reco.size(),TLorentzVector());

  for (unsigned int i=0; i<recind.size();i++) {
    TLorentzVector tlv;
    tlv.SetXYZM(mc.at(mcind.at(i)).momentum.x,mc.at(mcind.at(i)).momentum.y,mc.at(mcind.at(i)).momentum.z,mc.at(mcind.at(i)).mass);
    result[recind.at(i)]=tlv;
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP2MC_px(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<float> result;
  result.resize(reco.size(),-1.);
  for (unsigned int i=0; i<recind.size();i++) {
    result[recind.at(i)]=mc.at(mcind.at(i)).momentum.x;
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP2MC_py(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<float> result;
  result.resize(reco.size(),-1.);
  for (unsigned int i=0; i<recind.size();i++) {
    result[recind.at(i)]=mc.at(mcind.at(i)).momentum.y;
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP2MC_pz(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<float> result;
  result.resize(reco.size(),-1.);
  for (unsigned int i=0; i<recind.size();i++) {
    result[recind.at(i)]=mc.at(mcind.at(i)).momentum.z;
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP2MC_pdg(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<float> result;
  result.resize(reco.size(),-1.);
  for (unsigned int i=0; i<recind.size();i++) {
    result[recind.at(i)]=mc.at(mcind.at(i)).PDG;
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP2MC_charge(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<float> result;
  result.resize(reco.size(),-1.);
  for (unsigned int i=0; i<recind.size();i++) {
    result[recind.at(i)]=mc.at(mcind.at(i)).charge;
  }
  return result;
}

ROOT::VecOps::RVec<float> getRP2MC_mass(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<float> result;
  result.resize(reco.size(),-1.);
  for (unsigned int i=0; i<recind.size();i++) {
    result[recind.at(i)]=mc.at(mcind.at(i)).mass;
  }
  return result;
}

ROOT::VecOps::RVec<int> getRP2MC_index(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco) {
  ROOT::VecOps::RVec<int> result;
  result.resize(reco.size(),-1.);
  for (size_t i=0; i<recind.size();i++) {
    result[recind.at(i)]=mcind.at(i);
  }
  return result;
}


ROOT::VecOps::RVec<int> getRP2MC_index_test(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> parents) {
  std::cout <<"=========NEW EVENT========="<<std::endl;
  ROOT::VecOps::RVec<int> result;
  result.resize(reco.size(),-1.);
  for (size_t i=0; i<recind.size();i++) {
    if (result[recind.at(i)]>-1){
      auto & p_prev = mc.at(result[recind.at(i)]);
      auto & p_now  = mc.at(mcind.at(i));
      auto & p_reco = reco.at(recind.at(i));
      TLorentzVector tlv_prev;
      TLorentzVector tlv_now;
      TLorentzVector tlv_reco;
      tlv_prev.SetXYZM(p_prev.momentum.x, p_prev.momentum.y, p_prev.momentum.z, p_prev.mass);
      tlv_now.SetXYZM(p_now.momentum.x, p_now.momentum.y, p_now.momentum.z, p_now.mass);
      tlv_reco.SetXYZM(p_reco.momentum.x, p_reco.momentum.y, p_reco.momentum.z, p_reco.mass);

      std::cout << "reco energy " << tlv_reco.E() << " eta " << tlv_reco.Eta() << " phi " << tlv_reco.Phi() << " previous PDG " << p_prev.PDG << " energy " << tlv_prev.E() << " eta " << tlv_prev.Eta() << " phi " << tlv_prev.Phi() << " dR reco " << tlv_reco.DeltaR(tlv_prev) << "  new PDG  " << p_now.PDG << " energy " << tlv_now.E() << " eta " << tlv_now.Eta() << " phi " << tlv_now.Phi() << " dR reco " << tlv_reco.DeltaR(tlv_now) <<std::endl;
      for (unsigned j = mc.at(result[recind.at(i)]).parents_begin; j != mc.at(result[recind.at(i)]).parents_end; ++j){
	std::cout << "   prev==index " << j <<" parents " << parents.at(j) << "  PDGID "<< mc.at(parents.at(j)).PDG << " px " << mc.at(parents.at(j)).momentum.x << "  status  " << mc.at(parents.at(j)).generatorStatus << std::endl;
	for (unsigned k = mc.at(parents.at(j)).parents_begin; k != mc.at(parents.at(j)).parents_end; ++k)
	  std::cout << "   prev==index " << k <<" grandparents " << parents.at(k) << "  PDGID "<< mc.at(parents.at(k)).PDG << " px " << mc.at(parents.at(k)).momentum.x << "  status  " << mc.at(parents.at(k)).generatorStatus << std::endl;
      }

      for (unsigned j = mc.at(mcind.at(i)).parents_begin; j != mc.at(mcind.at(i)).parents_end; ++j){
	std::cout << "   now==index  " << j <<" parents " << parents.at(j) << "  PDGID "<< mc.at(parents.at(j)).PDG << " px " << mc.at(parents.at(j)).momentum.x << "  status  " << mc.at(parents.at(j)).generatorStatus << std::endl;
	for (unsigned k = mc.at(parents.at(j)).parents_begin; k != mc.at(parents.at(j)).parents_end; ++k)
	  std::cout << "   now==index " << k <<" grandparents " << parents.at(k) << "  PDGID "<< mc.at(parents.at(k)).PDG << " px " << mc.at(parents.at(k)).momentum.x << "  status  " << mc.at(parents.at(k)).generatorStatus << std::endl;
      }
    }
    result[recind.at(i)]=mcind.at(i);
  }
  return result;
}



ROOT::VecOps::RVec<int> getRP2MC_parentid (ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> parents){
  ROOT::VecOps::RVec<int> result;
  result.resize(reco.size(),-1.);
  for (unsigned int i=0; i<recind.size();i++) {
    if (mc.at(mcind.at(i)).parents_begin!=mc.at(mcind.at(i)).parents_end){
      result[recind.at(i)]=parents.at(mc.at(mcind.at(i)).parents_begin);
    }
  }

  /*  if (recind.size()>reco.size()){ 
    std::cout << recind.size() <<"========="<<reco.size()<<std::endl;
    for (unsigned int i=0; i<recind.size();i++) {
      if (i<recind.size()-1 && recind[i]==recind[i+1]){
	
	TLorentzVector tlv;
	tlv.SetXYZM(mc.at(mcind.at(i)).momentum.x,mc.at(mcind.at(i)).momentum.y,mc.at(mcind.at(i)).momentum.z,mc.at(mcind.at(i)).mass);
	TLorentzVector tlv2;
	tlv2.SetXYZM(reco.at(recind.at(i)).momentum.x,reco.at(recind.at(i)).momentum.y,reco.at(recind.at(i)).momentum.z,reco.at(recind.at(i)).mass);
	std::cout << "n mc " << mc.size() << " rec ind " << recind.at(i) << " reco P "<< tlv2.P()<< "  mc ind " << mcind.at(i) << " truth P " << tlv.P() << " pdg_id " << mc.at(mcind.at(i)).PDG  << "  parent id " << parents.at(mc.at(mcind.at(i)).parents_begin) << " parent pdg id " << mc.at(parents.at(mc.at(mcind.at(i)).parents_begin)).PDG << std::endl;

	tlv.SetXYZM(mc.at(mcind.at(i+1)).momentum.x,mc.at(mcind.at(i+1)).momentum.y,mc.at(mcind.at(i+1)).momentum.z,mc.at(mcind.at(i+1)).mass);
	tlv2.SetXYZM(reco.at(recind.at(i+1)).momentum.x,reco.at(recind.at(i+1)).momentum.y,reco.at(recind.at(i+1)).momentum.z,reco.at(recind.at(i+1)).mass);
	std::cout << "n mc " << mc.size() << " rec ind " << recind.at(i+1) << " reco P "<< tlv2.P()<< "  mc ind " << mcind.at(i+1) << " truth P " << tlv.P() << " pdg_id " << mc.at(mcind.at(i+1)).PDG  << "  parent id " << parents.at(mc.at(mcind.at(i+1)).parents_begin) << " parent pdg id " << mc.at(parents.at(mc.at(mcind.at(i+1)).parents_begin)).PDG << std::endl;
	}
    }
    }*/
  return result;
}


ROOT::VecOps::RVec<float>  getRP2MC_p_func::operator() (ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<float> result;
  result.resize(reco.size(),-1.);

  for (unsigned int i=0; i<recind.size();i++) {
    TLorentzVector tlv;
    tlv.SetXYZM(mc.at(mcind.at(i)).momentum.x,mc.at(mcind.at(i)).momentum.y,mc.at(mcind.at(i)).momentum.z,mc.at(mcind.at(i)).mass);
    result[recind.at(i)]=tlv.P();
  }

  if (recind.size()>reco.size()){ 
    std::cout << recind.size() <<"========="<<reco.size()<<std::endl;
     for (unsigned int i=0; i<recind.size();i++) {
       TLorentzVector tlv;
       tlv.SetXYZM(mc.at(mcind.at(i)).momentum.x,mc.at(mcind.at(i)).momentum.y,mc.at(mcind.at(i)).momentum.z,mc.at(mcind.at(i)).mass);
       TLorentzVector tlv2;
       tlv2.SetXYZM(reco.at(recind.at(i)).momentum.x,reco.at(recind.at(i)).momentum.y,reco.at(recind.at(i)).momentum.z,reco.at(recind.at(i)).mass);
       std::cout << "n mc " << mc.size() << " rec ind " << recind.at(i) << " reco P "<< tlv2.P()<< "  mc ind " << mcind.at(i) << " truth P " << tlv.P() << " pdg_id " << mc.at(mcind.at(i)).PDG << " parent_begin " <<  mc.at(mcind.at(i)).parents_begin << " parent_end " <<  mc.at(mcind.at(i)).parents_end << " daut_begin " <<  mc.at(mcind.at(i)).daughters_begin << " daut_end " <<  mc.at(mcind.at(i)).daughters_end <<std::endl;
     }
  }
  return result;
}


// -------------------------------------------------------------------------------------------------

// -- select RecoParticles associated with MC muons
// -- ( for muons from JPsi, can not use the Muon collection because it oontains
// -- only the isolated muons)

selRP_PDG::selRP_PDG( int arg_pdg, bool arg_chargedOnly ): m_PDG(arg_pdg), m_chargedOnly(arg_chargedOnly)  {} ;
std::vector<edm4hep::ReconstructedParticleData>  selRP_PDG::operator() (ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {

  std::vector<edm4hep::ReconstructedParticleData> result;

  for (int i=0; i<recind.size();i++) {
      int reco_idx = recind.at(i);
      int mc_idx = mcind.at(i);
      int pdg = mc.at(mc_idx).PDG ;
      if ( m_chargedOnly ) {
        if ( reco.at( reco_idx ).charge ==0 ) continue;
      }
      if ( std::abs( pdg ) == std::abs( m_PDG)  ) {
         result.push_back( reco.at( reco_idx ) ) ;
      }
  }
  return result;
}

// -------------------------------------------------------------------------------------------------

// -- select RecoParticles associated with a charged hadron :
// -- take all charged RecoParticles that are not associated with  a MC lepton

std::vector<edm4hep::ReconstructedParticleData> selRP_ChargedHadrons (ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {

  std::vector<edm4hep::ReconstructedParticleData> result;

  for (int i=0; i<recind.size();i++) {
      int reco_idx = recind.at(i);
      int mc_idx = mcind.at(i);
      int pdg = mc.at(mc_idx).PDG ;
      if ( reco.at( reco_idx ).charge == 0 ) continue;
      if ( std::abs( pdg ) == 11 || std::abs( pdg ) == 13 || std::abs( pdg ) == 15 ) continue ;
      result.push_back( reco.at( reco_idx ) ) ;
  }

  return result;
}


// -------------------------------------------------------------------------------------------------

// -- select the reco'ed particles associated with MC muons that come from the
// -- decay of a J/Psi

/// selMuons_JPsimatch::selMuons_JPsimatch( int arg_dum ) : m_dummy(arg_dum) {};

selMuons_JPsimatch::selMuons_JPsimatch( int pdg_mother, int pdg_daughter1, int pdg_daughter2) {
  m_pdg_mother = pdg_mother;
  m_pdg_daughter1 = pdg_daughter1;
  m_pdg_daughter2 = pdg_daughter2;
} ;

    
std::vector<edm4hep::ReconstructedParticleData> selMuons_JPsimatch::operator() (ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> mcdaughters) {

    
  std::vector<edm4hep::ReconstructedParticleData> result;
  std::vector< std::array<int, 2> >  MCmuons_from_JPsis = get_MC_legs_from_mothers( m_pdg_mother, m_pdg_daughter1, m_pdg_daughter2, true)( mc, mcdaughters) ;
  int nJPsis = MCmuons_from_JPsis.size() ;
  //std::cout << " number of mothers in this event  = " << nJPsis << std::endl;

  if ( nJPsis <1) return result ;


  for( int ijpsi=0; ijpsi < nJPsis; ijpsi ++) {

      std::array<int, 2>  MCmuons_JPsi = MCmuons_from_JPsis[ ijpsi ] ;
      if (MCmuons_JPsi[0]<0 || MCmuons_JPsi[1]<0) continue;

      // std::cout << " indices of the MC legs = " << MCmuons_JPsi[0] << " " << MCmuons_JPsi[1] << std::endl;

      int nlegs =0;
      for (int i=0; i<recind.size();i++) {
          int reco_idx = recind.at(i);
          // keep only charged particles
          if ( reco.at( reco_idx ).charge == 0 ) continue;
          int mc_idx = mcind.at(i);
          if ( mc_idx == MCmuons_JPsi[0] || mc_idx == MCmuons_JPsi[1]) {
             result.push_back( reco.at( reco_idx ) ) ;	 
             nlegs ++;
          }
      }

      if ( nlegs == 1 ) {   // one of the muons from this JPsi didnot make a RecoParticle
			    // e.g. outside of the tracker acceptance
       			    // in which case, remove the other leg if it was found
       	result.pop_back() ;
      }

  }  // loop over the JPsis

  return result;
}

// -------------------------------------------------------------------------------------------------

// -- select the reco'ed particles associated with Bs -> JPsi(mumu) phi(KK)

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> selRP_Bs2JPsiPhi( ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> mcdaughters) {

  std::vector<edm4hep::ReconstructedParticleData> result;

  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> muons_from_JPsi = selMuons_JPsimatch( 531, 13, -13)( recind, mcind, reco, mc, mcdaughters);
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> kaons_from_Phi = selMuons_JPsimatch( 531, 321, -321)( recind, mcind, reco, mc, mcdaughters);

  //if ( muons_from_JPsi.size() != 2 ) std::cout << " -- problem with this Bs, #muons from JPsi = " << muons_from_JPsi.size() << std::endl; 
  //if ( kaons_from_Phi.size() != 2 )  std::cout << " -- problem with this Bs, #kaons from Phi = " << kaons_from_Phi.size() << std::endl;

  if ( muons_from_JPsi.size() == 2 && kaons_from_Phi.size() == 2 ) {
  // keep only the first Bs, anyway cases with 2 Bs are rare
    result.push_back( muons_from_JPsi[0] );
    result.push_back( muons_from_JPsi[1] );
    result.push_back(kaons_from_Phi[0]);
    result.push_back(kaons_from_Phi[1]);
  }

 return ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>( result);
}


// -------------------------------------------------------------------------------------------------

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> selChargedRP_MCmatch_daughtersOf( int index, ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> ind ) {

  std::vector<edm4hep::ReconstructedParticleData> result;

  std::vector<int>  stable_daughters = list_of_stable_particles_from_decay( index, mc, ind );

  //for ( int j=0; j < stable_daughters.size(); j++) {
      //std::cout << "      Ds decayed into PDG = " << mc[stable_daughters[j]].PDG << std::endl;
  //}

      int nlegs =0;
      for (int i=0; i<recind.size();i++) {
          int reco_idx = recind.at(i);
          // keep only charged particles
          if ( reco.at( reco_idx ).charge == 0 ) continue;
          int mc_idx = mcind.at(i);
      
          if ( std::find( stable_daughters.begin(), stable_daughters.end(), mc_idx ) != stable_daughters.end() ) {
		result.push_back( reco.at( reco_idx ) ) ;
		nlegs ++;
	  }

      }

   return ROOT::VecOps::RVec( result );
  
}

// -------------------------------------------------------------------------------------------------

int getTrack2MC_index (  int track_index,
                                        ROOT::VecOps::RVec<int> recind,
                                        ROOT::VecOps::RVec<int> mcind,
                                        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco) {

  int mc_index = -1;

      for (int i=0; i<recind.size();i++) {
          int reco_idx = recind.at(i);
          // keep only charged particles
          if ( reco.at( reco_idx ).charge == 0 ) continue;
          mc_index = mcind.at(i);
          if ( reco.at( reco_idx ).tracks_begin == track_index ) return mc_index;
      }
 return mc_index;
}


	  









