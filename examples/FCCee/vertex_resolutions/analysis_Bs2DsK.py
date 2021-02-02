import sys
import ROOT

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

               .Alias("Particle1", "Particle#1.index")

               # MC event primary vertex
               .Define("MC_PrimaryVertex",  "getMC_EventPrimaryVertex(21)( Particle )" )

               # number of tracks
               .Define("ntracks","get_nTracks(EFlowTrack_1)")

               # Select primary tracks based on the matching to MC
		  # This can be used  to select primary tracks when the
		  # gen-level primary vertex  is not  (0,0,0)
               .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
               .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
               .Define("PrimaryTracks",  "SelPrimaryTracks(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle, MC_PrimaryVertex)" )
               .Define("nPrimaryTracks", "getRP_n(PrimaryTracks)")
               # Reconstruct the vertex from these primary tracks :
               .Define("Vertex_primaryTracks",  "VertexFitter ( 1, PrimaryTracks, EFlowTrack_1) ")   # primary vertex, in mm

               # Reco to MC indices :
               .Define( "RP2MC_indices",  "getRP2MC_index( MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles )")

               # Bs to Ds K decay : indices
               .Define("Bs2DsK_indices",  "getMC_indices_Bs2DsK(  Particle,  Particle1  )")

	       # Number of decays Bs -> Ds K
               .Define("nBs",  "MC_nBs2DsK( Bs2DsK_indices )")

               # the MC Ds from the Bs  decay
               .Define("Ds",  "getMC_Ds( Particle, Bs2DsK_indices )")
               # the Ds   energy   
               .Define("Ds_E",  "getMC_e( Ds )")
               .Define("Ds_pt", "getMC_pt( Ds ) ")
               .Define("Ds_theta", "getMC_theta( Ds )")
               .Define("Ds_phi", "getMC_phi( Ds )")

               # The MC legs from the Ds decay
               .Define("Kplus",   "getMC_Kplus( Particle, Particle1, Bs2DsK_indices)")
               .Define("Kminus",   "getMC_Kminus( Particle, Particle1, Bs2DsK_indices)")
               .Define("Piplus",   "getMC_Piplus( Particle, Particle1, Bs2DsK_indices)")

              
               # the MC bachelor K from the Ds decay
               .Define("BachelorK",  "getMC_BachelorK( Particle, Bs2DsK_indices )")
               .Define("BachelorK_E",  "getMC_e( BachelorK )")
               .Define("BachelorK_theta",  "getMC_theta( BachelorK )")

               # Decay vertex of the Ds
               .Define("DsDecayVertex", " MC_Ds_DecayVertex( Particle,  Particle1, Bs2DsK_indices  )")
               # Decay vertex of the Bs
               .Define("BsDecayVertex", " MC_Bs_DecayVertex( Particle,  Particle1, Bs2DsK_indices  )")

	       # Tracks associated with Ds (from Bs) decays :
               .Define("DsTracks",  "selChargedRP_DsfromBs( MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle, Particle1, Bs2DsK_indices )" )
	       # Number of such tracks :
               .Define("nDsTracks",  "n_DsTracks( DsTracks )")

               # the RecoParticles associated with the K+, K- and Pi+ of teh Ds decay
               .Define("RecoKplus",  "RecoKplus( DsTracks) ")
               .Define("RecoKminus", "RecoKminus( DsTracks)")
               .Define("RecoPiplus", "RecoPiplus( DsTracks)")

               # Reco'ed vertex of the Ds
               .Define("DsVertexObject",  "DsVertex( DsTracks, EFlowTrack_1) ")
               .Define("DsVertex",  " get_VertexData( DsVertexObject )")
             
	
               # Reconstructed Ds
                .Define("RecoDs", " ReconstructedDs( MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle, Particle1, Bs2DsK_indices )" )
                .Define("RecoDs_atVertex",  "ReconstructedDs_atVertex( DsVertexObject, MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles, Particle)" )

                .Define("RecoDs_pt",  "getRP_pt( RecoDs )")
                .Define("RecoDs_theta",  "getRP_theta( RecoDs ) ")
                .Define("RecoDs_phi",   "getRP_phi( RecoDs )")

                .Define("RecoDs_atVertex_pt",  "getRP_pt( RecoDs_atVertex )")
                .Define("RecoDs_atVertex_theta",  "getRP_theta( RecoDs_atVertex )")
                .Define("RecoDs_atVertex_phi",   "getRP_phi( RecoDs_atVertex )")


               # the RecoParticles associated with the bachelor K's
               .Define("RecoBachelorK",  "selChargedRP_KfromBs( RP2MC_indices, ReconstructedParticles, Bs2DsK_indices )" )

            
               # the list of tracks to reconstruct the Bs vertex
               .Define("BsTracks",  "tracks_for_fitting_the_Bs_vertex( RecoDs_atVertex, RecoBachelorK, EFlowTrack_1 )")

               # Reco'ed Bs vertex
               .Define("BsVertexObject",  "BsVertex( BsTracks, EFlowTrack_1) ")
               .Define("BsVertex",  "get_VertexData( BsVertexObject ) ")


               # Bs to Ds(KKpi) K  decay ?
               #.Define("Bsdecay",  "getMC_decay(531, 443, false)(Particle, Particle1)")
               #.Define("Bsbardecay",  "getMC_decay(-531, 443, false)(Particle, Particle1)")

               # MC Bs
               #.Define("Bs",  "selMC_PDG(531,  false)(Particle)")
               #.Define("Bsbar",  "selMC_PDG(-531,  false)(Particle)")
               #.Define("nBs",  "getMC_n( Bs)")
               #.Define("nBsbar",  "getMC_n( Bsbar)")

               # Decay vertex of theBs
		# does not work...  the endpoint is not filled in the files
               #.Define("BsDecayVertex",  "getMC_endPoint( Bs )")

               #.Define("BsDecayVertex",   "getMC_decayVertex(531, false)( Particle, Particle1)")

               #.Define("BsTracks",   "selRP_Bs2JPsiPhi(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle,Particle1)")
               #.Define("BsVertex",   "VertexFitter( 1, BsTracks, EFlowTrack_1)" )
        )


        # select branches for output file
        branchList = ROOT.vector('string')()
        for branchName in [
                #"MC_PrimaryVertex",
                #"ntracks",
                #"nPrimaryTracks",
                #"Vertex_primaryTracks",     # on Zuds: both track selections lead to very similar results for the vertex
                #"Bs2DsK_indices",
                "nBs",
                "Ds",
                "Ds_E",
                "Ds_pt",
                "Ds_theta",
                "Ds_phi",
                "BachelorK",
                "BachelorK_E",
                "BachelorK_theta",
                "Kplus",
                "Kminus",
                "Piplus",
                "DsDecayVertex",
                "BsDecayVertex",
                #"DsTracks",
                "nDsTracks",
                "RecoKplus",
                "RecoKminus",
                "RecoPiplus",
                "DsVertex",
                "RecoDs",
                "RecoDs_atVertex",
                "RecoDs_pt",
                "RecoDs_theta",
                "RecoDs_phi",
                "RecoDs_atVertex_pt",
                "RecoDs_atVertex_theta",
                "RecoDs_atVertex_phi",
                "RecoBachelorK",
                "BsVertex"

                ]:
            branchList.push_back(branchName)
        df2.Snapshot("events", self.outname, branchList)



if __name__ == "__main__":

    if len(sys.argv)==1:
        print ("usage:")
        print ("python ",sys.argv[0]," file.root")
        sys.exit(3)
    infile = sys.argv[1]
    outDir = 'FCCee/'+sys.argv[0].split('/')[1]+'/'
    import os
    os.system("mkdir -p {}".format(outDir))
    outfile = outDir+infile.split('/')[-1]
    ncpus = 0
    analysis = analysis(infile, outfile, ncpus)
    analysis.run()

    tf = ROOT.TFile(infile)
    entries = tf.events.GetEntries()
    p = ROOT.TParameter(int)( "eventsProcessed", entries)
    outf=ROOT.TFile(outfile,"UPDATE")
    p.Write()
