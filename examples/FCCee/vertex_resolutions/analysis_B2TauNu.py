import sys
import ROOT

PDGID=541
Filter=""
if PDGID==541:
    Filter="filterMC_pdgID(541, true)(Particle)==true"
elif PDGID==521:
    Filter="(filterMC_pdgID(521, false)(Particle)==true && filterMC_pdgID(-521, false)(Particle)==false) || (filterMC_pdgID(521, false)(Particle)==false && filterMC_pdgID(-521, false)(Particle)==true)"

print ("Load cxx analyzers ... ",)
ROOT.gSystem.Load("libedm4hep")
ROOT.gSystem.Load("libpodio")
ROOT.gSystem.Load("libFCCAnalyses")
ROOT.gErrorIgnoreLevel = ROOT.kFatal
_edm  = ROOT.edm4hep.ReconstructedParticleData()
_pod  = ROOT.podio.ObjectID()
_fcc  = ROOT.getMC_px

print ('edm4hep  ',_edm)
print ('podio    ',_pod)
print ('fccana   ',_fcc)

ROOT.gInterpreter.Declare("""
edm4hep::Vector3d MyMCDecayVertex(ROOT::VecOps::RVec<edm4hep::Vector3d> in1, ROOT::VecOps::RVec<edm4hep::Vector3d> in2) {

   edm4hep::Vector3d vertex(1e12, 1e12, 1e12); 
   if ( in1.size() == 0 && in2.size()==0) {
      std::cout <<"no vtx " <<std::endl;
      return vertex;
   }

   if ( in1.size() == 1 && in2.size()==0) vertex=in1[0];
   else if ( in1.size() == 0 && in2.size()==1) vertex=in2[0];
   else{
std::cout << "in1.size() " << in1.size() << "  in2.size() " <<in1.size()<< std::endl;
   }
   return vertex;
}

""")

ROOT.gInterpreter.Declare("""
float MyMinEnergy(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {

   float min=999999.;
   for (auto & p: in) {
    if (p.energy<min && p.energy>0) min=p.energy;
  }
  return min;
}

""")

class analysis():

    #__________________________________________________________
    def __init__(self, inputlist, outname, ncpu):
        self.outname = outname
        if ".root" not in outname:
            self.outname+=".root"

        #ROOT.ROOT.EnableImplicitMT(ncpu)

        self.df = ROOT.RDataFrame("events", inputlist)
        print (" done")
    #__________________________________________________________
    def run(self):
        #df2 = (self.df.Range(1000)
        df2 = (self.df
               .Filter(Filter)

               .Alias("Particle1", "Particle#1.index")
               .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
               .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")

               # MC event primary vertex
               .Define("MC_PrimaryVertex",  "getMC_EventPrimaryVertex(21)( Particle )" )

               # number of tracks
               .Define("ntracks","get_nTracks(EFlowTrack_1)")

               # Retrieve the decay vertex of all MC particles
               .Define("MC_DecayVertices",  "getMC_endPoint( Particle, Particle1)" )

               # MC indices of the decay Bc/Bu -> nu nu pi+ pi- pi+
               #     - if the event contains > 1 such decays, the first one is kept
               #     - the cases Bs -> Bsbar -> mu mu K K are included here
               #   first boolean: if true, look at the stable daughters, otherwise at the intermediate daughters
               #   second boolean: if true, include the charge conjugate decays 
               .Define("B2NuNuPiPiPi_indices",   "getMC_indices_ExclusiveDecay( %s, { 16, -16, 211, -211, 211 }, true, false)( Particle, Particle1)"%(PDGID))
               .Define("Bbar2NuNuPiPiPi_indices","getMC_indices_ExclusiveDecay( -%s, { -16, 16, -211, 211, -211 }, true, false)( Particle, Particle1)"%(PDGID))

               .Define("Piminus", "selMC_leg(4) ( B2NuNuPiPiPi_indices , Particle)" )
               .Define("Piplus",  "selMC_leg(4) ( Bbar2NuNuPiPiPi_indices , Particle)" )
               .Define("TauMCDecayVertex",  "MyMCDecayVertex(getMC_vertex(Piminus),getMC_vertex(Piplus))")

               # the MC Bc or Bcbar:
               .Define("B",  "if (B2NuNuPiPiPi_indices.size()>0) return selMC_leg(0) (B2NuNuPiPiPi_indices , Particle ); else return selMC_leg(0) (Bbar2NuNuPiPiPi_indices , Particle );")

                # Kinematics of the B :
               .Define("B_theta", "getMC_theta( B )")
               .Define("B_phi",   "getMC_phi( B )")
               .Define("B_e",     "getMC_e( B )")
               .Define("B_charge","getMC_charge( B )")

                # and the MC legs of the B :
               .Define("Nu1_vec",  "if (B2NuNuPiPiPi_indices.size()>0) return selMC_leg(1) (B2NuNuPiPiPi_indices , Particle ); else return selMC_leg(1) (Bbar2NuNuPiPiPi_indices , Particle );")
               .Define("Nu2_vec",  "if (B2NuNuPiPiPi_indices.size()>0) return selMC_leg(2) (B2NuNuPiPiPi_indices , Particle ); else return selMC_leg(2) (Bbar2NuNuPiPiPi_indices , Particle );")
               .Define("Pion1_vec",  "if (B2NuNuPiPiPi_indices.size()>0) return selMC_leg(3) (B2NuNuPiPiPi_indices , Particle ); else return selMC_leg(3) (Bbar2NuNuPiPiPi_indices , Particle );")
               .Define("Pion2_vec",  "if (B2NuNuPiPiPi_indices.size()>0) return selMC_leg(4) (B2NuNuPiPiPi_indices , Particle ); else return selMC_leg(4) (Bbar2NuNuPiPiPi_indices , Particle );")
               .Define("Pion3_vec",  "if (B2NuNuPiPiPi_indices.size()>0) return selMC_leg(5) (B2NuNuPiPiPi_indices , Particle ); else return selMC_leg(5) (Bbar2NuNuPiPiPi_indices , Particle );")

               .Define("Nu1e_vec",     "getMC_e( Nu1_vec )")
               .Define("Nu2e_vec",     "getMC_e( Nu2_vec )")
               .Define("Pion1e_vec",   "getMC_e( Pion1_vec )")
               .Define("Pion2e_vec",   "getMC_e( Pion2_vec )")
               .Define("Pion3e_vec",   "getMC_e( Pion3_vec )")

               .Define("Nu1",  "if (Nu1e_vec.size()>0 && Nu2e_vec.size()>0 && Nu1e_vec[0] > Nu2e_vec[0]) return Nu1_vec; else return Nu2_vec;")
               .Define("Nu2",  "if (Nu1e_vec.size()>0 && Nu2e_vec.size()>0 && Nu1e_vec[0] < Nu2e_vec[0]) return Nu1_vec; else return Nu2_vec;")
               .Define("Pion1", "if (Pion1e_vec.size()>0 && Pion2e_vec.size()>0 && Pion3e_vec.size()>0 && Pion1e_vec[0] > Pion2e_vec[0] && Pion1e_vec[0] > Pion3e_vec[0]) return Pion1_vec; else if (Pion1e_vec.size()>0 && Pion2e_vec.size()>0 && Pion3e_vec.size()>0 && Pion2e_vec[0]>Pion1e_vec[0] && Pion2e_vec[0]>Pion3e_vec[0]) return Pion2_vec; else return Pion3_vec;")
               .Define("Pion2", "if (Pion1e_vec.size()>0 && Pion2e_vec.size()>0 && Pion3e_vec.size()>0 && Pion1e_vec[0] < Pion2e_vec[0] && Pion1e_vec[0] > Pion3e_vec[0]) return Pion1_vec; else if (Pion1e_vec.size()>0 && Pion2e_vec.size()>0 && Pion3e_vec.size()>0 && Pion2e_vec[0]<Pion1e_vec[0] && Pion2e_vec[0]>Pion3e_vec[0]) return Pion2_vec; else return Pion3_vec;")
               .Define("Pion3", "if (Pion1e_vec.size()>0 && Pion2e_vec.size()>0 && Pion3e_vec.size()>0 && Pion1e_vec[0] < Pion2e_vec[0] && Pion1e_vec[0] < Pion3e_vec[0]) return Pion1_vec; else if (Pion1e_vec.size()>0 && Pion2e_vec.size()>0 && Pion3e_vec.size()>0 && Pion2e_vec[0]<Pion1e_vec[0] && Pion2e_vec[0]<Pion3e_vec[0]) return Pion2_vec; else return Pion3_vec;")

               
               # Kinematics of the B decay:
               .Define("Nu1_theta", "getMC_theta( Nu1 )")
               .Define("Nu1_phi",   "getMC_phi( Nu1 )")
               .Define("Nu1_e",     "getMC_e( Nu1 )")
               .Define("Nu1_charge","getMC_charge( Nu1 )")

               .Define("Nu2_theta", "getMC_theta( Nu2 )")
               .Define("Nu2_phi",   "getMC_phi( Nu2 )")
               .Define("Nu2_e",     "getMC_e( Nu2 )")
               .Define("Nu2_charge","getMC_charge( Nu2 )")

               .Define("Pion1_theta", "getMC_theta( Pion1 )")
               .Define("Pion1_phi",   "getMC_phi( Pion1 )")
               .Define("Pion1_e",     "getMC_e( Pion1 )")
               .Define("Pion1_charge","getMC_charge( Pion1 )")

               .Define("Pion2_theta", "getMC_theta( Pion2 )")
               .Define("Pion2_phi",   "getMC_phi( Pion2 )")
               .Define("Pion2_e",     "getMC_e( Pion2 )")
               .Define("Pion2_charge","getMC_charge( Pion2 )")

               .Define("Pion3_theta", "getMC_theta( Pion3 )")
               .Define("Pion3_phi",   "getMC_phi( Pion3 )")
               .Define("Pion3_e",     "getMC_e( Pion3 )")
               .Define("Pion3_charge","getMC_charge( Pion3 )")

               # Returns the RecoParticles associated with the 5 Bc decay products.
               # The size of this collection is always 5 provided that Bc2TauNuNuPiPiPi_indices is not empty,
               # possibly including "dummy" particles in case one of the leg did not make a RecoParticle.
               # This is on purpose, to maintain the mapping with the indices - i.e. the 1st particle in 
               # the list is the mu+, then the mu-, etc.
               .Define("BRecoParticles",  "if (B2NuNuPiPiPi_indices.size()>0) return selRP_matched_to_list( B2NuNuPiPiPi_indices, MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle); else return selRP_matched_to_list( Bbar2NuNuPiPiPi_indices, MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle);")

               #Returns the pion with minimum energy
               .Define("minPionE", "MyMinEnergy(BRecoParticles)")

               # the corresponding tracks - here, dummy particles, if any, are removed
               .Define("BTracks",   "getRP2TRK( BRecoParticles, EFlowTrack_1)" )

               # number of tracks used to reconstruct the Bc vertex 
               .Define("n_BTracks", "getTK_n( BTracks )")

               # the reco'ed vertex :
               .Define("BVertexObject",   "VertexFitter_Tk( 2, BTracks)" )
               .Define("BVertex",  "get_VertexData( BVertexObject )")

               # Angular separation between the tracks from the Bc decays
               .Define("deltaAlpha_max","angular_separation(0)( BRecoParticles )")
               .Define("deltaAlpha_min","angular_separation(1)( BRecoParticles )")
               .Define("deltaAlpha_ave","angular_separation(2)( BRecoParticles )")










               
               
                # and the MC legs of the Bs :
               #.Define("Muplus",  " selMC_leg(1)( Bs2MuMuKK_indices, Particle )")
               #.Define("Muminus",  " selMC_leg(2)( Bs2MuMuKK_indices, Particle )")
               #.Define("Kplus",  " selMC_leg(3)( Bs2MuMuKK_indices, Particle )")
               #.Define("Kminus",  " selMC_leg(4)( Bs2MuMuKK_indices, Particle )")

                # Kinematics of the Bs legs (MC) :
               #.Define("Muplus_theta",  "getMC_theta( Muplus )")
               #.Define("Muplus_phi",  "getMC_phi( Muplus )")
               #.Define("Muplus_e",  "getMC_e( Muplus )")
               #.Define("Muminus_theta",  "getMC_theta( Muminus )")
               #.Define("Muminus_phi",  "getMC_phi( Muminus )")
               #.Define("Muminus_e",  "getMC_e( Muminus )")
               #.Define("Kplus_theta",  "getMC_theta( Kplus )")
               #.Define("Kplus_phi",  "getMC_phi( Kplus )")
               #.Define("Kplus_e",  "getMC_e( Kplus )")
               #.Define("Kminus_theta",  "getMC_theta( Kminus )")
               #.Define("Kminus_phi",  "getMC_phi( Kminus )")
               #.Define("Kminus_e",  "getMC_e( Kminus )")

               #.Define("Bc_theta",   "getMC_theta( Bc )")
               #.Define("Bc_phi",   "getMC_phi( Bc )")
               #.Define("Bc_e",   "getMC_e( Bc )")

               
               # Decay vertex of theBs
		# does not work...  the endpoint is not filled in the files
               #.Define("BsDecayVertex",  "getMC_endPoint( Bs )")

               # Careful with the following: if Bs -> Bsbar, this returns the prod vertex of the Bsbar !
               #.Define("BsDecayVertex",   "getMC_decayVertex(531, false)( Particle, Particle1)")

               # Better use a custom method in Bs2JPsiPhi :
               #.Define("BsMCDecayVertex",   "BsMCDecayVertex( Bs2MuMuKK_indices, Particle )")

               # Returns the RecoParticles associated with the 4 Bs decay products.
               # The size of this collection is always 4 provided that Bs2MuMuKK_indices is not empty,
               # possibly including "dummy" particles in case one of the leg did not make a RecoParticle.
               # This is on purpose, to maintain the mapping with the indices - i.e. the 1st particle in 
               # the list is the mu+, then the mu-, etc.
               #.Define("BsRecoParticles",  "selRP_matched_to_list( Bs2MuMuKK_indices, MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")


               # the reco'ed Bs legs
               #.Define("RecoMuplus",   "selRP_leg(0)( BsRecoParticles )")
               #.Define("RecoMuminus",  "selRP_leg(1)( BsRecoParticles )")
               #.Define("RecoKplus",    "selRP_leg(2)( BsRecoParticles )")
               #.Define("RecoKminus",   "selRP_leg(3)( BsRecoParticles )")
               # and the kinematics :
               #.Define("RecoMuplus_theta",  "getRP_theta( RecoMuplus )")
               #.Define("RecoMuplus_phi",  "getRP_phi( RecoMuplus )")
               #.Define("RecoMuplus_e",  "getRP_e( RecoMuplus )")
               #.Define("RecoMuminus_theta",  "getRP_theta( RecoMuminus )")
               #.Define("RecoMuminus_phi",  "getRP_phi( RecoMuminus )")
               #.Define("RecoMuminus_e",  "getRP_e( RecoMuminus )")
               #.Define("RecoKplus_theta",  "getRP_theta( RecoKplus )")
               #.Define("RecoKplus_phi",  "getRP_phi( RecoKplus )")
               #.Define("RecoKplus_e",  "getRP_e( RecoKplus )")
               #.Define("RecoKminus_theta",  "getRP_theta( RecoKminus )")
               #.Define("RecoKminus_phi",  "getRP_phi( RecoKminus )")
               #.Define("RecoKminus_e",  "getRP_e( RecoKminus )")

               # the reco'ed legs, with the momenta at the Bs decay vertex
	       # This is not needed for the vertexing here. This was mostly to check
     	       # my propagation of the tracks from their point of dca to the Bs vertex.
               #.Define("RecoMuplus_atVertex",  "selRP_leg_atVertex(0) ( BsRecoParticles, BsVertexObject, EFlowTrack_1 )")
               #.Define("RecoMuplus_atVertex_theta",   "getRP_theta( RecoMuplus_atVertex )")
               #.Define("RecoMuplus_atVertex_phi",   "getRP_phi( RecoMuplus_atVertex )")
               #.Define("RecoMuminus_atVertex",  "selRP_leg_atVertex(1) ( BsRecoParticles, BsVertexObject, EFlowTrack_1 )")
               #.Define("RecoMuminus_atVertex_theta",   "getRP_theta( RecoMuminus_atVertex )")
               #.Define("RecoMuminus_atVertex_phi",   "getRP_phi( RecoMuminus_atVertex )")
               #.Define("RecoKplus_atVertex",  "selRP_leg_atVertex(2) ( BsRecoParticles, BsVertexObject, EFlowTrack_1 )")
               #.Define("RecoKplus_atVertex_theta",   "getRP_theta( RecoKplus_atVertex )")
               #.Define("RecoKplus_atVertex_phi",   "getRP_phi( RecoKplus_atVertex )")
               #.Define("RecoKminus_atVertex",  "selRP_leg_atVertex(3) ( BsRecoParticles, BsVertexObject, EFlowTrack_1 )")
               #.Define("RecoKminus_atVertex_theta",   "getRP_theta( RecoKminus_atVertex )")
               #.Define("RecoKminus_atVertex_phi",   "getRP_phi( RecoKminus_atVertex )")


        )


        # select branches for output file
        branchList = ROOT.vector('string')()
        for branchName in [
                #"RP_MC_parentindex",
                "MC_PrimaryVertex",
                "ntracks",
                "BVertex",
                "TauMCDecayVertex",
                "n_BTracks",
                "deltaAlpha_max",
                "deltaAlpha_min",
                "deltaAlpha_ave",
                "minPionE",
                "B_theta",
                "B_phi",
                "B_e",
                "B_charge",
                
                "Nu1_theta",
                "Nu1_phi",
                "Nu1_e",
                "Nu1_charge",
                
                "Nu2_theta",
                "Nu2_phi",
                "Nu2_e",
                "Nu2_charge",

                "Pion1_theta",
                "Pion1_phi",
                "Pion1_e",
                "Pion1_charge",
 
                "Pion2_theta",
                "Pion2_phi",
                "Pion2_e",
                "Pion2_charge",

                "Pion3_theta",
                "Pion3_phi",
                "Pion3_e",
                "Pion3_charge",

                ]:
            branchList.push_back(branchName)
        df2.Snapshot("events", self.outname, branchList)

# python examples/FCCee/vertex_resolutions/analysis_Bc2TauNu.py "/eos/experiment/fcc/ee/generation/DelphesEvents/fcc_tmp_v02/p8_ee_Zbb_ecm91_EvtGen_Bc2TauNuTAUHADNU/events_13*" events_Bc2TauNuTAUHADNU.root
# python examples/FCCee/vertex_resolutions/analysis_Bc2TauNu.py "/eos/experiment/fcc/ee/generation/DelphesEvents/fcc_tmp_v02/p8_ee_Zbb_ecm91_EvtGen_Bu2TauNuTAUHADNU/events_130*" events_Bu2TauNuTAUHADNU_test.root

if __name__ == "__main__":

    if len(sys.argv)==1:
        print ("usage:")
        print ("python ",sys.argv[0]," file.root")
        print ("python ",sys.argv[0]," dir/*.root")
        sys.exit(3)


    import glob
    filelist = glob.glob(sys.argv[1])
    print ("Create dataframe object from ", )
    fileListRoot = ROOT.vector('string')()
    for fileName in filelist:
        fileListRoot.push_back(fileName)
        print (fileName, " ",)
        print (" ...")
        
    outDir = sys.argv[0].replace(sys.argv[0].split('/')[0],'outputs/').replace('analysis_Bc2TauNu.py','')+'/'
    import os
    os.system("mkdir -p {}".format(outDir))
    
    outfile = outDir+'events.root'
    if len(sys.argv)==3:
        outfile = outDir+sys.argv[2]
    ncpus = 8
    analysis = analysis(fileListRoot, outfile, ncpus)
    analysis.run()

    

