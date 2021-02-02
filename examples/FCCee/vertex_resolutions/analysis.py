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
               #.Define("Vertex_primaryTracks",  "VertexFitter ( 1, PrimaryTracks, EFlowTrack_1) ")   # primary vertex, in mm

               # Bs to JPsi decay ?
               .Define("Bsdecay",  "getMC_decay(531, 443, false)(Particle, Particle1)")
               .Define("Bsbardecay",  "getMC_decay(-531, 443, false)(Particle, Particle1)")

               # MC Bs
               .Define("Bs",  "selMC_PDG(531,  false)(Particle)")
               .Define("Bsbar",  "selMC_PDG(-531,  false)(Particle)")
               .Define("nBs",  "getMC_n( Bs)")
               .Define("nBsbar",  "getMC_n( Bsbar)")

               # Decay vertex of theBs
		# does not work...  the endpoint is not filled in the files
               #.Define("BsDecayVertex",  "getMC_endPoint( Bs )")

               .Define("BsDecayVertex",   "getMC_decayVertex(531, false)( Particle, Particle1)")

               .Define("BsTracks",   "selRP_Bs2JPsiPhi(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle,Particle1)")
               .Define("BsVertexObject",   "VertexFitter( 2, BsTracks, EFlowTrack_1)" )
               .Define("BsVertex",  "get_oneVertexData( BsVertexObject )")

               # number of tracks used to reconstruct the Bs vertex 
               .Define("n_BsTracks", "getRP_n( BsTracks )")

               # Angular separation between the tracks from the Bs decays
               .Define("deltaAlpha_max","angular_separation(0)( BsTracks )")
               .Define("deltaAlpha_min","angular_separation(1)( BsTracks )")
               .Define("deltaAlpha_ave","angular_separation(2)( BsTracks )")

        )


        # select branches for output file
        branchList = ROOT.vector('string')()
        for branchName in [
                "MC_PrimaryVertex",
                "ntracks",
                "nPrimaryTracks",
                #"Vertex_primaryTracks",     # on Zuds: both track selections lead to very similar results for the vertex
                "Bsdecay",
                "Bsbardecay",
                "BsDecayVertex",
                "nBs",
                "nBsbar",
                "BsVertex",
                "BsTracks",
                "n_BsTracks",
                "deltaAlpha_max",
                "deltaAlpha_min",
                "deltaAlpha_ave",

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
