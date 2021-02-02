#include "Bs2DsK.h"

// -- select the reco'ed particles associated with the bachelor K in Bs -> Ds K 


//ROOT::VecOps::RVec< edm4hep::ReconstructedParticleData>  selChargedRP_KfromBs ( ROOT::VecOps::RVec<int> RP2MC_indices, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco, ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec< std::vector<int> > ndecays  ){

ROOT::VecOps::RVec< edm4hep::ReconstructedParticleData>  selChargedRP_KfromBs ( ROOT::VecOps::RVec<int> RP2MC_indices, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco, ROOT::VecOps::RVec< std::vector<int> > ndecays  ){

 ROOT::VecOps::RVec< edm4hep::ReconstructedParticleData>   results;

 edm4hep::ReconstructedParticleData dummy;
 dummy.energy = -99;

 //ROOT::VecOps::RVec< std::vector<int> > ndecays = getMC_indices_Bs2DsK( mc, mcind );

  for (int idecay = 0; idecay < ndecays.size(); idecay++) {
    std::vector<int> indices = ndecays[idecay];

    bool Found_K = false;
    //for (int id = 0; id < indices.size(); id++) {
      //int idx = indices[id];
      int idx = indices[1];  // the Kaon
      //if ( mc.at(idx).PDG  == -321 ) {   // the K-

       	for ( int j=0; j < reco.size(); j++) {
		if ( RP2MC_indices[j] == idx ) {
		    Found_K = true;
		    results.push_back( reco.at( j) ) ;
		    break;
		}
 	}

        if ( ! Found_K)  // did not find a reco Particle matched tp the Kaon (e.g. K was outside the acceptance)
	   results.push_back( dummy );

      //}
    //}
  }

 return results;
}



// -- select the reco'ed particles associated with Ds (from Bs) -> Phi(KK) Pi+

ROOT::VecOps::RVec< ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> > selChargedRP_DsfromBs( ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,  ROOT::VecOps::RVec<int>  ind, ROOT::VecOps::RVec< std::vector<int> > ndecays  ){

 std::vector<  ROOT::VecOps::RVec< edm4hep::ReconstructedParticleData> >  results;

 //ROOT::VecOps::RVec< std::vector<int> > ndecays = getMC_indices_Bs2DsK( mc, mcind );

  for (int idecay = 0; idecay < ndecays.size(); idecay++) {
    std::vector<int> indices = ndecays[idecay];

    //for (int id = 0; id < indices.size(); id++) {
      //int idx = indices[id];
      //if ( mc.at(idx).PDG  == 431 ) {   // the Ds

      int idx = indices[0];   // the Ds

       //std::cout << " ... Ds found, enter in selChargedRP_MCmatch_daughtersOf  " << std::endl;
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> the_Ds_tracks = selChargedRP_MCmatch_daughtersOf( idx, recind, mcind, reco, mc, ind);
       //std::cout << "    ... number of matched particles = " << the_Ds_tracks.size() << std::endl;
        results.push_back(  the_Ds_tracks );

      //}
    //}
  }

 return ROOT::VecOps::RVec( results );
}

ROOT::VecOps::RVec< int > n_DsTracks( ROOT::VecOps::RVec< ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> > all_Ds_tracks ) {
  std::vector<int> res;
  for (int i=0; i < all_Ds_tracks.size(); i++) {
    int ntracks = all_Ds_tracks[i].size();
    res.push_back( ntracks );
  }
 return ROOT::VecOps::RVec( res );
}

// ----------------------------------------------------------------------------------------------

// the RecoParticles associated with the Ds daughters, Ds+ -> K+ K- Pi+

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> RecoKplus( ROOT::VecOps::RVec< ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> > thelegs) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> results;
  	edm4hep::ReconstructedParticleData dummy;
  	bool found = false;
  	for (int idecay=0; idecay < thelegs.size(); idecay++) {
    	ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs = thelegs[idecay];
    	for (int ileg = 0; ileg < legs.size(); ileg ++) {
		if ( legs[ileg].charge > 0 && legs[ileg].mass > 0.45 && legs[ileg].mass < 0.55) {  // the K+
			found = true;
			results.push_back( legs[ileg] );
			break;
		}
    	}
    	if (! found) results.push_back( dummy );
  	}
   return results;
}

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> RecoKminus( ROOT::VecOps::RVec< ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> > thelegs) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> results;
  	edm4hep::ReconstructedParticleData dummy;
  	bool found = false;
  	for (int idecay=0; idecay < thelegs.size(); idecay++) {
    	ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs = thelegs[idecay];
    	for (int ileg = 0; ileg < legs.size(); ileg ++) {
        	if ( legs[ileg].charge < 0 && legs[ileg].mass > 0.45 && legs[ileg].mass < 0.55) {  // the K-
                	found = true;
                	results.push_back( legs[ileg] );
                	break;
        	}
    	}
    	if (! found) results.push_back( dummy );
  	}
  return results;
}

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> RecoPiplus( ROOT::VecOps::RVec< ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> > thelegs) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> results;
        edm4hep::ReconstructedParticleData dummy;
        bool found = false;
        for (int idecay=0; idecay < thelegs.size(); idecay++) {
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs = thelegs[idecay];
        for (int ileg = 0; ileg < legs.size(); ileg ++) {
                if ( legs[ileg].charge > 0 && legs[ileg].mass > 0.1 && legs[ileg].mass < 0.2) {  // the Piplus
                        found = true;
                        results.push_back( legs[ileg] );
                        break;
                }
        }
        if (! found) results.push_back( dummy );
        }
  return results;
}


ROOT::VecOps::RVec< FCCAnalysesVertex >  DsVertex( ROOT::VecOps::RVec< ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> > all_Ds_tracks,
                                                   ROOT::VecOps::RVec<edm4hep::TrackState> tracks  ) {
   std::vector< FCCAnalysesVertex > vertices;
   for (int i=0; i < all_Ds_tracks.size(); i++) {
      ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> theTracks = all_Ds_tracks[i];
      FCCAnalysesVertex vertex = VertexFitter(3, theTracks, tracks) ;
      vertices.push_back( vertex) ;
   }
  return ROOT::VecOps::RVec< FCCAnalysesVertex >( vertices );
}


// ----------------------------------------------------------------------------------------------------------------------------

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  ReconstructedDs( ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> ind, ROOT::VecOps::RVec< std::vector<int> > ndecays  ) {

	// reconstruct the Ds from the pre-fit tracks
	// no updated parameters from the vertex fit, and (px, py) correspond to the dca
	// and not to the vertex

 std::vector< edm4hep::ReconstructedParticleData>   result;

   for (int idecay = 0; idecay < ndecays.size(); idecay++) {
    std::vector<int> indices = ndecays[idecay];
    edm4hep::ReconstructedParticleData theDs;

    for (int id = 0; id < indices.size(); id++) {
      int idx = indices[id];
      if ( mc.at(idx).PDG  == 431 ) {   // the Ds

      TLorentzVector theDs4v;
      float Ds_charge = 0;

      std::vector<int>  stable_daughters = list_of_stable_particles_from_decay( idx, mc, ind );

      for (int i=0; i<recind.size();i++) {
          int reco_idx = recind.at(i);
          // keep only charged particles
          if ( reco.at( reco_idx ).charge == 0 ) continue;
          int mc_idx = mcind.at(i);

          if ( std::find( stable_daughters.begin(), stable_daughters.end(), mc_idx ) != stable_daughters.end() ) {
	        TLorentzVector leg;
		leg.SetXYZM( reco.at( reco_idx ).momentum.x, reco.at( reco_idx ).momentum.y, reco.at( reco_idx ).momentum.z, mc.at(mc_idx).mass );
		theDs4v +=  leg;
		Ds_charge += mc.at(mc_idx).charge;
          }
      }

      theDs.charge =  Ds_charge;
      theDs.momentum.x = theDs4v.Px();
      theDs.momentum.y = theDs4v.Py();
      theDs.momentum.z = theDs4v.Pz();
      theDs.mass = theDs4v.M();

      }
    }
    result.push_back( theDs );
    }
 return ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>( result );
}

// ----------------------------------------------------------------------------------------------------------------------------


ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  ReconstructedDs_atVertex( ROOT::VecOps::RVec<FCCAnalysesVertex> DsVertices,
		ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, 	
		ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {


        // reconstruct the Ds from the post-fit tracks

 std::vector< edm4hep::ReconstructedParticleData>   result;

   for (int idecay = 0; idecay < DsVertices.size(); idecay++) {		// loop over the Ds decays (in case > 1 Ds->KKPi in the event)
    FCCAnalysesVertex DsVertex = DsVertices[idecay];

    edm4hep::ReconstructedParticleData theDs;

    ROOT::VecOps::RVec<int> the_tracks_indices = DsVertex.reco_ind ;	// indices of the TrackStates that made this vertex
    ROOT::VecOps::RVec< TVector3 >  updated_track_momentum_at_vertex = DsVertex.updated_track_momentum_at_vertex ;

    TLorentzVector theDs4v;
    float Ds_charge = 0;
    //std::cout << " here a Ds decay " << std::endl;
    //std::cout << " Number of tracks at the vertex " << the_tracks_indices.size() << std::endl;
    //std::cout << " updated_track_momentum_at_vertex.size() = " << updated_track_momentum_at_vertex.size() << std::endl;

    for ( int itra = 0; itra < the_tracks_indices.size(); itra ++) {
	
	int index = the_tracks_indices[itra];
	//std::cout << " track index " << index << std::endl;
	TVector3 track_momentum = updated_track_momentum_at_vertex[itra];
        //std::cout << " track momentum : " ;  track_momentum.Print();

	int mc_index = getTrack2MC_index( index , recind, mcind, reco ) ;   // ssociation of the track to a MC Particle
        //std::cout << " ... a track matched to PDG = " << mc.at( mc_index ).PDG << std::endl;
	TLorentzVector leg;
	leg.SetXYZM( track_momentum.Px(), track_momentum.Py(), track_momentum.Pz(), mc.at( mc_index ).mass );
	theDs4v +=  leg;
	Ds_charge += mc.at(mc_index).charge;
        
    } // end loop over the tracks

    theDs.charge =  Ds_charge;
    theDs.momentum.x = theDs4v.Px();
    theDs.momentum.y = theDs4v.Py();
    theDs.momentum.z = theDs4v.Pz();
    theDs.mass = theDs4v.M();
    theDs.referencePoint = DsVertex.vertex.position ;   // the Ds decay vertex

    result.push_back( theDs );

    }   // end loop over the Ds's

 return ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>( result );
}


// ----------------------------------------------------------------------------------------------------------------------------


ROOT::VecOps::RVec<edm4hep::TrackState>  ReconstructedDs_atVertex_TrackState( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> theDss) {

// return a TrackState corresponding to the reco'ed Ds

  ROOT::VecOps::RVec<edm4hep::TrackState > result;

  float norm = 1e-3;   // to convert from mm to meters 
  for (int iDs = 0; iDs < theDss.size(); iDs++) {

	edm4hep::ReconstructedParticleData theDs = theDss[iDs];
      	TVector3 vertex( theDs.referencePoint.x * norm, theDs.referencePoint.y * norm, theDs.referencePoint.z * norm ); 
	TVector3 momentum ( theDs.momentum.x, theDs.momentum.y, theDs.momentum.z)  ;
     	// get the corresponding track parameters using Franco's code 
     	TVectorD track_param = XPtoPar( vertex, momentum, theDs.charge );

	edm4hep::TrackState track;
	track.D0 	= track_param[0] * 1e3 ; // from meters to mm
	track.phi	= track_param[1];
	track.omega	= track_param[2] / ( 0.5*1e3 ) ; // C from Franco = rho/2, and convert from m-1 to mm-1
	track.Z0	= track_param[3] * 1e3  ;   // from meters to mm
	track.tanLambda = track_param[4];

	track.referencePoint.x = vertex[0];
        track.referencePoint.y = vertex[1];
        track.referencePoint.z = vertex[2];

	// assume here that the parameters are perfectly measured
	std::array<float, 15> covMatrix;
	for (int icov=0; icov < 15; icov++) {
	   covMatrix[icov] = 1e-15;
	}
	track.covMatrix = covMatrix;

	result.push_back( track );

  }

  return ROOT::VecOps::RVec<edm4hep::TrackState>( result );
}

// ----------------------------------------------------------------------------------------------------------------------------


ROOT::VecOps::RVec<  ROOT::VecOps::RVec<edm4hep::TrackState> > tracks_for_fitting_the_Bs_vertex
			( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> RecoDs_atVertex, 
			  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  BachelorK,
			  ROOT::VecOps::RVec<edm4hep::TrackState> all_tracks  ) {

 
 ROOT::VecOps::RVec<  ROOT::VecOps::RVec<edm4hep::TrackState> > results;

 ROOT::VecOps::RVec<edm4hep::TrackState> DsPseudoTrackState = ReconstructedDs_atVertex_TrackState( RecoDs_atVertex );

 int ndecays = DsPseudoTrackState.size();
 if ( BachelorK.size() != ndecays) {
   std::cout << " --- problem in tracks_for_fitting_the_Bs_vertex :  DsPseudoTrackState.size() = " << ndecays << " but BachelorK.size() = " << BachelorK.size() << std::endl;
   return results;
 }

 for (int idecay = 0; idecay < ndecays; idecay ++) {
      	std::vector< edm4hep::TrackState> thetracks;
	thetracks.push_back( DsPseudoTrackState[idecay] );    // the pseudo-Ds track
	int track_index_BachelorK = BachelorK[idecay].tracks_begin ;

	if ( track_index_BachelorK < all_tracks.size() ) {
	    thetracks.push_back(  all_tracks.at( track_index_BachelorK ) ) ;
	}
	else {
	    std::cout << " --- problem in tracks_for_fitting_the_Bs_vertex: track index of the bachelor K = " << track_index_BachelorK << " but size of the track collection is " << all_tracks.size() << std::endl;
	}

	results.push_back( ROOT::VecOps::RVec<edm4hep::TrackState> (thetracks) );
 }

 return results;
}


ROOT::VecOps::RVec< FCCAnalysesVertex >  BsVertex( ROOT::VecOps::RVec< ROOT::VecOps::RVec<edm4hep::TrackState> > all_Bs_tracks,
                                                   ROOT::VecOps::RVec<edm4hep::TrackState> tracks  ) {
   std::vector< FCCAnalysesVertex > vertices;
   for (int i=0; i < all_Bs_tracks.size(); i++) {
      ROOT::VecOps::RVec<edm4hep::TrackState> theTracks = all_Bs_tracks[i];
      FCCAnalysesVertex vertex = VertexFitter_Tk(2, theTracks, tracks) ;
      vertices.push_back( vertex) ;
   }
  return ROOT::VecOps::RVec< FCCAnalysesVertex >( vertices );
}







// ---------------------------------------------------------------------------------------------------------------------------

/*
edm4hep::TrackState Track_from_RP( edm4hep::ReconstructedParticleData> particle) {

  edm4hep::TrackState track;
  float pt = sqrt( pow( particle.momentum.x,2) + pow( particle.momentum.y,2) );
  track.omega = ;
  
}


// ---------------------------------------------------------------------------------------------------------------------------

ROOT::VecOps::RVec<edm4hep::TrackState> ReconstructedDsTrackState( ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> ind, ROOT::VecOps::RVec< std::vector<int> > ndecays  ) {

 std::vector< edm4hep::TrackState>   result;

   for (int idecay = 0; idecay < ndecays.size(); idecay++) {
    std::vector<int> indices = ndecays[idecay];
    edm4hep::ReconstructedParticleData theDs;

    for (int id = 0; id < indices.size(); id++) {
      int idx = indices[id];
      if ( mc.at(idx).PDG  == 431 ) {   // the Ds

      TLorentzVector theDs4v;
      float Ds_charge = 0;

      std::vector<int>  stable_daughters = list_of_stable_particles_from_decay( idx, mc, ind );

      for (int i=0; i<recind.size();i++) {
          int reco_idx = recind.at(i);
          // keep only charged particles
          if ( reco.at( reco_idx ).charge == 0 ) continue;
          int mc_idx = mcind.at(i);

          if ( std::find( stable_daughters.begin(), stable_daughters.end(), mc_idx ) != stable_daughters.end() ) {
                TLorentzVector leg;
                leg.SetXYZM( reco.at( reco_idx ).momentum.x, reco.at( reco_idx ).momentum.y, reco.at( reco_idx ).momentum.z, mc.at(mc_idx).mass );
                theDs4v +=  leg;
                Ds_charge += mc.at(mc_idx).charge;
          }
      }

      theDs.charge =  Ds_charge;
      theDs.momentum.x = theDs4v.Px();
      theDs.momentum.y = theDs4v.Py();
      theDs.momentum.z = theDs4v.Pz();
      theDs.mass = theDs4v.M();

      }
    }
    result.push_back( theDs );
    }
 return ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>( result );
}

*/


