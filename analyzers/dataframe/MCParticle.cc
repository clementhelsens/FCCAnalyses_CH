#include "MCParticle.h"


ROOT::VecOps::RVec<float> getMC_pt(ROOT::VecOps::RVec<edm4hep::MCParticleData> in){
 ROOT::VecOps::RVec<float> result;
 for (size_t i = 0; i < in.size(); ++i) {
   result.push_back(sqrt(in[i].momentum.x * in[i].momentum.x + in[i].momentum.y * in[i].momentum.y));
 }
 return result;
}

ROOT::VecOps::RVec<edm4hep::MCParticleData> mergeParticles(ROOT::VecOps::RVec<edm4hep::MCParticleData> x, ROOT::VecOps::RVec<edm4hep::MCParticleData> y) {
  //to be keept as std::vector
  std::vector<edm4hep::MCParticleData> result;
  result.reserve(x.size() + y.size());
  result.insert( result.end(), x.begin(), x.end() );
  result.insert( result.end(), y.begin(), y.end() );
  return ROOT::VecOps::RVec(result);
}


ROOT::VecOps::RVec<float> getMC_time(ROOT::VecOps::RVec<edm4hep::MCParticleData> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.time);
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_pdg(ROOT::VecOps::RVec<edm4hep::MCParticleData> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.PDG);
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_genStatus(ROOT::VecOps::RVec<edm4hep::MCParticleData> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.generatorStatus);
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_simStatus(ROOT::VecOps::RVec<edm4hep::MCParticleData> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.simulatorStatus);
  }
  return result;
}

ROOT::VecOps::RVec<edm4hep::Vector3d> getMC_vertex(ROOT::VecOps::RVec<edm4hep::MCParticleData> in){
  ROOT::VecOps::RVec<edm4hep::Vector3d> result;
  for (auto & p: in) {
    result.push_back(p.vertex);
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_vertex_x(ROOT::VecOps::RVec<edm4hep::MCParticleData> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.vertex.x);
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_vertex_y(ROOT::VecOps::RVec<edm4hep::MCParticleData> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.vertex.y);
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_vertex_z(ROOT::VecOps::RVec<edm4hep::MCParticleData> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.vertex.z);
  }
  return result;
}

ROOT::VecOps::RVec<edm4hep::Vector3d> getMC_endPoint(ROOT::VecOps::RVec<edm4hep::MCParticleData> in){
  ROOT::VecOps::RVec<edm4hep::Vector3d> result;
  for (auto & p: in) {
    result.push_back(p.endpoint);
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_endPoint_x(ROOT::VecOps::RVec<edm4hep::MCParticleData> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.endpoint.x);
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_endPoint_y(ROOT::VecOps::RVec<edm4hep::MCParticleData> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.endpoint.y);
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_endPoint_z(ROOT::VecOps::RVec<edm4hep::MCParticleData> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.endpoint.z);
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_mass(ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.mass);
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_eta(ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.Eta());
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_phi(ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.Phi());
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_e(ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.E());
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_p(ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.P());
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_px(ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.momentum.x);
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_py(ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.momentum.y);
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_pz(ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.momentum.z);
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_charge(ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.charge);
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_y(ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.Rapidity());
  }
  return result;
}

ROOT::VecOps::RVec<float> getMC_theta(ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.Theta());
  }
  return result;
}

ROOT::VecOps::RVec<TLorentzVector> getMC_tlv(ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  ROOT::VecOps::RVec<TLorentzVector> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv);
  }
  return result;
}

int getMC_n(ROOT::VecOps::RVec<edm4hep::MCParticleData> x) {
  int result =  x.size();
  return result;
}



selMC_pT::selMC_pT(float arg_min_pt) : m_min_pt(arg_min_pt) {};

ROOT::VecOps::RVec<edm4hep::MCParticleData>  selMC_pT::operator() (ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  ROOT::VecOps::RVec<edm4hep::MCParticleData> result;
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    if (std::sqrt(std::pow(p.momentum.x,2) + std::pow(p.momentum.y,2)) > m_min_pt) {
      result.emplace_back(p);
    }
  }
  return result;
}

selMC_genStatus::selMC_genStatus(int arg_status) : m_status(arg_status) {};
ROOT::VecOps::RVec<edm4hep::MCParticleData>  selMC_genStatus::operator() (ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  ROOT::VecOps::RVec<edm4hep::MCParticleData> result;
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    if (p.generatorStatus == m_status) {
      result.emplace_back(p);
    }
  }
  return result;
}


selMC_PDG::selMC_PDG(int arg_pdg, bool arg_chargeconjugate) : m_pdg(arg_pdg), m_chargeconjugate( arg_chargeconjugate )  {};

std::vector<edm4hep::MCParticleData>  selMC_PDG::operator() (ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  std::vector<edm4hep::MCParticleData> result;
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    if ( m_chargeconjugate ) {
      	if ( std::abs( p.PDG ) == std::abs( m_pdg)  ) result.emplace_back(p);
    }
    else {
	if ( p.PDG == m_pdg ) result.emplace_back(p);
    }
  }
  return result;
}


getMC_decay::getMC_decay(int arg_mother, int arg_daughters, bool arg_inf){m_mother=arg_mother; m_daughters=arg_daughters; m_inf=arg_inf;};
bool getMC_decay::operator() (ROOT::VecOps::RVec<edm4hep::MCParticleData> in,  ROOT::VecOps::RVec<int> ind){

  bool result=false;
  for (size_t i = 0; i < in.size(); ++i) {
    if (in[i].PDG!=m_mother)continue;
    // std::cout << "Here a mother with PDG = " << m_mother << std::endl;
    int ndaughters=0;
    for (unsigned j = in.at(i).daughters_begin; j != in.at(i).daughters_end; ++j) {
      //std::cout << "       daughter : " << in[ind.at(j)].PDG << std::endl;
      if (std::abs(in[ind.at(j)].PDG)==m_daughters && m_inf==false)ndaughters+=1;
      else if (std::abs(in[ind.at(j)].PDG)<=m_daughters && m_inf==true)ndaughters+=1;
    }
    //if (ndaughters>1){
    if (ndaughters>=1){
      result=true;
      return result;
    }
  }
  return result;
}


ROOT::VecOps::RVec<int> getMC_parentid(ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> parents){
  ROOT::VecOps::RVec<int> result;
  /*std::cout <<"================== Full Truth=================" <<std::endl;
  for (size_t i = 0; i < mc.size(); ++i) {
    std::cout << "i= " << i << "  PDGID "<< mc.at(i).PDG  <<  "  status  " << mc.at(i).generatorStatus << std::endl;
    for (unsigned j = mc.at(i).parents_begin; j != mc.at(i).parents_end; ++j) 
      std::cout << "   ==index " << j <<" parents " << parents.at(j) << "  PDGID "<< mc.at(parents.at(j)).PDG << "  status  " << mc.at(parents.at(j)).generatorStatus << std::endl;

  }*/
    
  //std::cout <<"================== NEW EVENT=================" <<std::endl;
  for (size_t i = 0; i < mcind.size(); ++i) {

    if (mcind.at(i)<0){
      result.push_back(-999);
      continue;
    }
    //std::cout << "mc ind " << mcind.at(i) << "  PDGID "<< mc.at(mcind.at(i)).PDG  << "  status  " << mc.at(mcind.at(i)).generatorStatus << std::endl;
    for (unsigned j = mc.at(mcind.at(i)).parents_begin; j != mc.at(mcind.at(i)).parents_end; ++j) {
      //std::cout << "   ==index " << j <<" parents " << parents.at(j) << "  PDGID "<< mc.at(parents.at(j)).PDG << "  status  " << mc.at(parents.at(j)).generatorStatus << std::endl;
      // result.push_back(parents.at(j));
    }
    //std::cout << mc.at(mcind.at(i)).parents_begin <<"---"<< mc.at(mcind.at(i)).parents_end<< std::endl;	
    if (mc.at(mcind.at(i)).parents_end - mc.at(mcind.at(i)).parents_begin>1) {
      //std::cout << "-999" << std::endl;
      result.push_back(-999);
    }
    else {
      //std::cout << "not -999 "<< parents.at(mc.at(mcind.at(i)).parents_begin) << std::endl;		    
      result.push_back(parents.at(mc.at(mcind.at(i)).parents_begin));
    }
  }
  return result;
}


getMC_tree::getMC_tree(int arg_index) : m_index(arg_index) {};
ROOT::VecOps::RVec<int> getMC_tree::operator() (ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind){
  ROOT::VecOps::RVec<int> result;
  auto & particle = in[m_index];
  
  //for (unsigned j = in.at(i).parents_begin; j != in.at(i).parents_end; ++j){
  //  if 
  //  result.push_back(ind.at(j));
  
  
  std::cout << "Thomas logic"<<std::endl;
  
  for (size_t i = 0; i < in.size(); ++i) {
    // all the other cout
    std::cout << i  << " status " << in[i].generatorStatus << " pdg " << in[i].PDG << " p beg "<< in.at(i).parents_begin << " p end " <<in.at(i).parents_end << "  mc size " << in.size() << "  ind size "<<ind.size() << std::endl;
    for (unsigned j = in.at(i).parents_begin; j != in.at(i).parents_end; ++j) {
      std::cout << "   ==index " << j <<" parents " << ind.at(j) << std::endl;
    }
  }
  //std::cout << "END Thomas logic"<<std::endl;

  /*  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    std::cout <<  "here" << std::endl;
    
    if (p.generatorStatus != m_index) continue;
    ROOT::VecOps::RVec<int> tree;
    tree.push_back(in.at(ind.at(i)).parents_begin);
    while(true){
      std::cout <<  "tree back " << tree.back() << std::endl;
      //      std::cout << 
      tree.push_back(in.at(ind.at(tree.back())).parents_begin);
    }
    result.push_back(tree);
  }
  return result;*/
  return result;
}




filterMC_pdgID::filterMC_pdgID(int arg_pdgid, bool arg_abs){m_pdgid = arg_pdgid; m_abs = arg_abs;};
bool  filterMC_pdgID::operator() (ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    if ((m_abs && abs(p.PDG) == m_pdgid) || (p.PDG == m_pdgid)) return true;
  }
  return false;
}


getMC_EventPrimaryVertex::getMC_EventPrimaryVertex( int arg_genstatus) { m_genstatus = arg_genstatus; };
TVector3 getMC_EventPrimaryVertex::operator() ( ROOT::VecOps::RVec<edm4hep::MCParticleData> in )  {
  TVector3 result(-1e12,-1e12,-1e12);
  for (auto & p: in) {
     if ( p.generatorStatus == m_genstatus ) {   // generator status code for the incoming particles of the hardest subprocess
       TVector3 res( p.vertex.x, p.vertex.y, p.vertex.z );
       result = res;
       break;
     }
   }

  return result;
}

// ----------------------------------------------------------------------------------------------------------------------------------

std::vector<int> list_of_stable_particles_from_decay( int i, ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind) {
  std::vector<int> res;

	// could maybe use getMC_tree above, but for the while, the latter
	   	// looks like work in progress

  // i = index of a MC particle in the Particle block
  // in = the Particle collection
  // ind = the block with the indices for the daughters, Particle#1.index

  // returns a vector with the indices (in the Particle block) of the stable daughters of the particle i

  int db = in.at(i).daughters_begin ;
  int de = in.at(i).daughters_end;
  ///std::cout << " in list_of_stable_particles_from_decay  i = " <<  i << " db " << db << " de " << de << std::endl;

  if ( db != de ) {// particle is unstable
    int d1 = ind[db] ;
    int d2 = ind[de-1];
    //std::cout << " d1 d2 " << d1 << " " << d2 << std::endl;
    for (int idaughter = d1; idaughter <= d2; idaughter++) {
      std::vector<int> rr = list_of_stable_particles_from_decay( idaughter, in, ind) ;
      res.insert( res.end(), rr.begin(), rr.end() );
    }
  }
  else {    // particle is stable
     res.push_back( i ) ;
     return res ;
  }
  return res;
}

// ----------------------------------------------------------------------------------------------------------------------------------

std::vector<int> list_of_particles_from_decay( int i, ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind) {
  std::vector<int> res;

  // i = index of a MC particle in the Particle block
  // in = the Particle collection
  // ind = the block with the indices for the daughters, Particle#1.index

  // returns a vector with the indices (in the Particle block) of the daughters of the particle i

  int db = in.at(i).daughters_begin ;
  int de = in.at(i).daughters_end;
  if  ( db == de ) return res;   // particle is stable
  int d1 = ind[db] ;
  int d2 = ind[de-1];
  for (int idaughter = d1; idaughter <= d2; idaughter++) {
     res.push_back( idaughter);
  }
  return res;
}

// ----------------------------------------------------------------------------------------------------------------------------------

get_MC_legs_from_mothers::get_MC_legs_from_mothers( int pdg_mother, int pdg_daughter1, int pdg_daughter2, bool stableDaughters) {
  m_pdg_mother = pdg_mother;
  m_pdg_daughter1 = pdg_daughter1;
  m_pdg_daughter2 = pdg_daughter2;
  m_stableDaughters = stableDaughters;
} ;


std::vector< std::array<int, 2> > get_MC_legs_from_mothers::operator() ( ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind) {

// Example : pdg_mother = 443 ( mother = JPsi )
//           pdg_daughter1 = 13, pdg_daughter2 = -13 (muons)

// Returns a vector of array<int, 2>  :
//    size of the vector = number of JPsis in this event
//    the array contains the indices of the mu- and the mu+,in this order

   std::vector< std::array<int, 2> >  result;

   std::vector<int> theMothers ;	  // contains the indices of the mothers (e.g. the JPsis) in the MCParticleblock

   for ( int i=0; i < in.size(); i++){
     int pdg = in[i].PDG ;
     if ( pdg == m_pdg_mother ){	
        theMothers.push_back( i );
     }
   }
   int nMothers = theMothers.size();	   // number of JPsis in this event
   bool Mothers  = ( nMothers >= 1);

//	std::cout << "    nMothers = " << nMothers << std::endl;

   if ( ! Mothers ) return result;

   ///if ( nMothers > 1 ) std::cout << " -- more than 1 JPsi : " << nMothers << std::endl;

   // debug :
                  /*
     for( int i=0; i < in.size(); i++){
           std::string here="";
           if ( in[i].PDG == 443 ) here ="   ---   here is a JPsi " ;
           int db = in.at(i).daughters_begin ;
           int de = in.at(i).daughters_end;
           if ( db != de ) {
                 std::cout << i << " pdg: " << in[i].PDG << " decay products from " << ind[db] << " to " << ind[de-1] << here << std::endl;
           }
           else {
                 std::cout << i << " pdg: " << in[i].PDG << " is stable " << here << std::endl ;
           }
     }
    
 */

   std::vector<int> allMuonsFromJPsis;  // that is used to handle potential duplicates
					// (e.g. a JPsi that comes itself from a  JPsi with a different status code)

   for( int i=0; i < theMothers.size(); i++) {   // loop over the JPsis

     std::array<int, 2> resu;
     resu[0] = -1;
     resu[1] = -1;

     int ijp = theMothers[i] ;

     std::vector<int> products ;
     if ( m_stableDaughters ) {
	products = list_of_stable_particles_from_decay( ijp, in, ind ) ;
     }
     else {
	products = list_of_particles_from_decay( ijp, in, ind ) ;
     }

	//debug:
	  std::cout << " result of list_of_stable_particles_from_decay " << std::endl;
          for (int idec=0; idec< products.size(); idec++) {
	    std::cout << "  pdg = " << in.at(products[idec]).PDG << std::endl;
	  }
     //if ( products.size() != 2 ) continue;
     if ( products.size() < 2 ) continue;

     
     //std::cout << " a D0 that decayed into " << in [ products[0]].PDG << " and " << in [ products[1]].PDG << std::endl;
     int n1 = 0;
     int n2 = 0;

     for ( int j=0; j < products.size(); ++j) {	   // did the JPsi decay into muons ?
        int idx = products[j] ;
        //std::cout << "     decay product : " << in[idx].PDG << std::endl ;
                if (in[idx].PDG == m_pdg_daughter1) {
		resu[0] = idx ;
		n1 ++;
	}
        if (in[idx].PDG == m_pdg_daughter2) {
		resu[1] = idx;
		n2 ++;
	}
        if ( in[idx].PDG == m_pdg_daughter1 || in[idx].PDG == m_pdg_daughter2) {
	   TLorentzVector leg1 ;
           leg1.SetXYZM( in[idx].momentum.x, in[idx].momentum.y, in[idx].momentum.z, in[idx].mass );
           //std::cout << " a MC leg mass = " << in[idx].mass << " pT = " << leg1.Pt() << " theta = " << leg1.Theta() << " eta = " << leg1.Eta() << " phi = " << leg1.Phi() << std::endl;
        }
     }

     if ( n1 == 0 || n2 == 0) continue;

     // to remove potential duplicated J/Psis :
     if(std::find(allMuonsFromJPsis.begin(), allMuonsFromJPsis.end(), resu[0]) != allMuonsFromJPsis.end()) {
	continue;  	// skip this JPsi
     }
     else {

	allMuonsFromJPsis.push_back( resu[0] );
        allMuonsFromJPsis.push_back( resu[1] );
        result.push_back( resu );
     }

   }  /// end loop over the JPsis

    //std::cout  << "    exit get_MC_muons_from_JPsis " << std::endl;
         return result;
}



// --------------------------------------------------------------------------------------------------
//   decay vertex of a MCparticle (because the "endpoints"  are not filled in the  files)


edm4hep::Vector3d  getMC_EndPoint( int index,  ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind )  {

 edm4hep::Vector3d result(-1e12, -1e12, -1e12);

 if (index >= in.size()) {
    std::cout << " invalid index in getMC_EndPoint" << std::endl;
    return result;   // invalid index
 }

     int db = in.at(index).daughters_begin ;
     int d1 = ind[db] ;   // first daughter
     result = in.at(d1).vertex ;

 return result;
  
}

// --------------------------------------------------------------------------------------------------
//   decay vertices of all particles with a given PDG ID

getMC_decayVertex::getMC_decayVertex( int pdg_mother, bool chargeconjugate) {
 m_pdg_mother = pdg_mother;
 m_chargeconjugate = chargeconjugate;
} ;

ROOT::VecOps::RVec<edm4hep::Vector3d> getMC_decayVertex::operator() ( ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind ) {

  std::vector< edm4hep::Vector3d> result;

  for ( int i=0; i< in.size(); i++){

     bool keep  =  in.at(i).PDG == m_pdg_mother ;
     if ( m_chargeconjugate) {
	keep = abs( in.at(i).PDG) == m_pdg_mother ;
     }
     if ( ! keep ) continue;
 

/*
     std::vector<int>  stable_daughters = list_of_stable_particles_from_decay( i, in, ind );
     if ( stable_daughters.size() < 2 ) continue;  // something weird..
     for ( int j=0 ; j< stable_daughters.size();  j++) {
	int idx = stable_daughters[j];
        result.push_back(in.at(idx).vertex);
	break;	//could actually add  a check that all daughters have the same vertex
     }
*/

     // int db = in.at(i).daughters_begin ;
     // int d1 = ind[db] ;   // first daughter
     // result.push_back(in.at(d1).vertex) ;

     result.push_back(  getMC_EndPoint(i, in, ind) );

  }
  return ROOT::VecOps::RVec(result);

}


// --------------------------------------------------------------------------------------------------

// specific to Bs -> Ds K :

ROOT::VecOps::RVec< std::vector<int> > getMC_indices_Bs2DsK( ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind)  {
  std::vector< std::vector<int> > results;

  for (int i=0; i < in.size(); i++) {
    if ( in.at(i).PDG == 531) {	// a Bs
	//std::cout << " --- here a Bs " << std::endl;
	// by construction, all Bs's should decay to Ds + K
	//	it could oscillate though, in which case there is only one daughter with PDG = -531
        std::vector<int> the_daughters = list_of_particles_from_decay( i, in, ind ) ;
	std::vector<int> daughters_ordered;    // reorder if needed: put the Ds first, then the bachelor K
	daughters_ordered.push_back( -1);
        daughters_ordered.push_back( -1);
        bool hasDs = false;
	bool hasK = false;
        for (int id = 0; id < the_daughters.size(); id++) {
	   int idx = the_daughters[id];
	   if ( in.at(idx).PDG == 431 ) {
		hasDs = true;
		daughters_ordered[0] = idx;
	   }
	   if ( in.at(idx).PDG == -321) {   // the bachelor K
		hasK = true;
		daughters_ordered[1] = idx;
	   }
	}
        if ( hasDs) {
		//results.push_back( the_daughters  );
		if ( the_daughters.size() != 2 ) std::cout << "   weird in getMC_indices_Bs2DsK: #daughters from Bs decay is " << the_daughters.size() << std::endl;
		if (! hasK) std::cout << "   weird in getMC_indices_Bs2DsK: no K found " << std::endl;
		results.push_back( daughters_ordered );  
	}
    } // end isBs
   
  }
 return ROOT::VecOps::RVec( results) ;
}

// --------------------------------------------------------------------------------------------------

/// number of decays Bs -> Ds K in this event  (could be 0 because not every event has produced a Bs,
///			and could be 0 or 2 because of Bs oscillations)
int MC_nBs2DsK( ROOT::VecOps::RVec< std::vector<int> >   ndecays ) {
   // ROOT::VecOps::RVec< std::vector<int> > ndecays = getMC_indices_Bs2DsK( in, ind );
   return ndecays.size();
}

// --------------------------------------------------------------------------------------------------

ROOT::VecOps::RVec<edm4hep::MCParticleData> getMC_Ds( ROOT::VecOps::RVec<edm4hep::MCParticleData> in,  ROOT::VecOps::RVec< std::vector<int> >   ndecays)  {
   std::vector< edm4hep::MCParticleData> results;
   for (int idecay = 0; idecay < ndecays.size(); idecay++) {
     std::vector<int> indices = ndecays[idecay];
     ///for (int id = 0; id < indices.size(); id++) {
      //int idx = indices[id];
      ///if ( in.at(idx).PDG  == 431 ) {   // the Ds
         int idx = indices[0] ;   // the Ds
         results.push_back( in.at(idx) );
      ///}
     ///}
   }
   return ROOT::VecOps::RVec<edm4hep::MCParticleData>( results );
}

ROOT::VecOps::RVec<edm4hep::MCParticleData> getMC_BachelorK( ROOT::VecOps::RVec<edm4hep::MCParticleData> in,  ROOT::VecOps::RVec< std::vector<int> >   ndecays)  {
   std::vector< edm4hep::MCParticleData> results;
   for (int idecay = 0; idecay < ndecays.size(); idecay++) {
     std::vector<int> indices = ndecays[idecay];
     //for (int id = 0; id < indices.size(); id++) {
      //int idx = indices[id];
      //if ( in.at(idx).PDG  == -321 ) {   // the bachelor K
         int idx = indices[1];     // the bachelor K
         results.push_back( in.at(idx) );
      //}
     //}
   }
   return ROOT::VecOps::RVec<edm4hep::MCParticleData>( results );
}

ROOT::VecOps::RVec<edm4hep::MCParticleData> getMC_Kplus( ROOT::VecOps::RVec<edm4hep::MCParticleData> in,  ROOT::VecOps::RVec<int> ind, ROOT::VecOps::RVec< std::vector<int> >   ndecays)  {
   std::vector< edm4hep::MCParticleData> results;
   edm4hep::MCParticleData dummy;
   for (int idecay = 0; idecay < ndecays.size(); idecay++) {
     std::vector<int> indices = ndecays[idecay];
     bool found = false;
     int idx = indices[0] ;   // the Ds 
     std::vector<int> Ds_daughters = list_of_stable_particles_from_decay( idx, in, ind ) ; 
     for (int ida=0; ida < Ds_daughters.size(); ida++) {
        int k = Ds_daughters[ida];
	if ( in.at(k).PDG == 321) {
		results.push_back( in.at(k) );
		found = true;
		break;
	}
     }
     if (! found) results.push_back( dummy );
   }
  return results;
}


ROOT::VecOps::RVec<edm4hep::MCParticleData> getMC_Kminus( ROOT::VecOps::RVec<edm4hep::MCParticleData> in,  ROOT::VecOps::RVec<int> ind, ROOT::VecOps::RVec< std::vector<int> >   ndecays)  {
   std::vector< edm4hep::MCParticleData> results;
   edm4hep::MCParticleData dummy;
   for (int idecay = 0; idecay < ndecays.size(); idecay++) {
     std::vector<int> indices = ndecays[idecay];
     bool found = false;
     int idx = indices[0] ;   // the Ds 
     std::vector<int> Ds_daughters = list_of_stable_particles_from_decay( idx, in, ind ) ; 
     for (int ida=0; ida < Ds_daughters.size(); ida++) {
        int k = Ds_daughters[ida];
        if ( in.at(k).PDG == -321) {
                results.push_back( in.at(k) );
                found = true;
                break;
        }
     }
     if (! found) results.push_back( dummy );
   }
  return results;
}

ROOT::VecOps::RVec<edm4hep::MCParticleData> getMC_Piplus( ROOT::VecOps::RVec<edm4hep::MCParticleData> in,  ROOT::VecOps::RVec<int> ind, ROOT::VecOps::RVec< std::vector<int> >   ndecays)  {
   std::vector< edm4hep::MCParticleData> results;
   edm4hep::MCParticleData dummy;
   for (int idecay = 0; idecay < ndecays.size(); idecay++) {
     std::vector<int> indices = ndecays[idecay];
     bool found = false;
     int idx = indices[0] ;   // the Ds 
     std::vector<int> Ds_daughters = list_of_stable_particles_from_decay( idx, in, ind ) ; 
     for (int ida=0; ida < Ds_daughters.size(); ida++) {
        int k = Ds_daughters[ida];
        if ( in.at(k).PDG == 211) {
                results.push_back( in.at(k) );
                found = true;
                break;
        }
     }
     if (! found) results.push_back( dummy );
   }
  return results;
}




ROOT::VecOps::RVec< edm4hep::Vector3d > MC_Ds_DecayVertex( ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind, ROOT::VecOps::RVec< std::vector<int> >   ndecays)  {

  std::vector< edm4hep::Vector3d > results;

  //ROOT::VecOps::RVec< std::vector<int> > ndecays = getMC_indices_Bs2DsK( in, ind );

  for (int idecay = 0; idecay < ndecays.size(); idecay++) {
    std::vector<int> indices = ndecays[idecay];
    bool hasDs = false;
    for (int id = 0; id < indices.size(); id++) {
      int idx = indices[id];
      if ( in.at(idx).PDG  == 431 ) {   // the Ds
        hasDs = true;
	edm4hep::Vector3d vertex  = getMC_EndPoint( idx, in, ind) ;
        results.push_back( vertex );
      }
    }
    if ( ! hasDs) std::cout << " weird.. no Ds  found in MC_Ds_DecayVertex " << std::endl;
  }
  
  return ROOT::VecOps::RVec( results );
}


ROOT::VecOps::RVec< edm4hep::Vector3d > MC_Bs_DecayVertex( ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind, ROOT::VecOps::RVec< std::vector<int> >   ndecays)  {

  std::vector< edm4hep::Vector3d > results;

  //ROOT::VecOps::RVec< std::vector<int> > ndecays = getMC_indices_Bs2DsK( in, ind );

  for (int idecay = 0; idecay < ndecays.size(); idecay++) {
    std::vector<int> indices = ndecays[idecay];
    bool hasDs = false;
    for (int id = 0; id < indices.size(); id++) {
      int idx = indices[id];
      if ( in.at(idx).PDG  == 431 ) {   // the Ds
        hasDs = true;
        edm4hep::Vector3d vertex  = in.at(idx).vertex;   // take the vertex of the Ds
        results.push_back( vertex );
      }
    }
    if ( ! hasDs) std::cout << " weird.. no Ds  found in MC_Bs_DecayVertex " << std::endl;
  }

  return ROOT::VecOps::RVec( results );
}




