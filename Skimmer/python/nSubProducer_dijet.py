import ROOT
import math, os, sys
import numpy as np
import fastjet
ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *

class nSubProdDijet(Module):
    def __init__(self, minJetPt=200., maxJetEta=2.4, minLeadJetPt=450, sysSource=[]):
        self.writeHistFile=True
        self.verbose = False
        self.beta = 0.0
        self.zcut = 0.1
        self.R = 0.8
        self.maxTau = 10

        ### Kinematics Cuts Jets ###
        self.minJetPt = minJetPt
        self.maxJetEta = maxJetEta
        self.minLeadJetPt = minLeadJetPt

        ### Kinenatic Cuts Muons ###
        self.minTightMuonPt = 53.
        self.maxTightMuonEta = 2.1
        self.minMuonPt = 20.
        self.maxMuonEta = 2.4

        ### Kinenatic Cuts Electrons ###
        self.minTightElectronPt = 120.
        self.minElectronPt = 35.
        self.range1ElectronEta = [0,1.442]
        self.range2ElectronEta = [1.56,2.5]

        self.sysSource = ['_nom'] + [ isys+i for i in [ 'Up', 'Down' ] for isys in sysSource if not isys.endswith('nom') ]

	self.kinematic_labels = ["_pt", "_eta", "_phi", "_mass"]
	self.nSub_labels = ["_tau_0p5_", "_tau_1_", "_tau_2_"]
        ## JES from https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources#Main_uncertainties_2016_80X
        self.JESLabels = [ "AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias", "Fragmentation", "SinglePionECAL", "SinglePionHCAL", "FlavorQCD", "TimePtEta", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeBal", "RelativeSample", "RelativeFSR", "RelativeStatFSR", "RelativeStatEC", "RelativeStatHF", "PileUpDataMC", "PileUpPtRef", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF", "PileUpMuZero", "PileUpEnvelope", "SubTotalPileUp", "SubTotalRelative", "SubTotalPt", "SubTotalScale", "SubTotalAbsolute", "SubTotalMC", "Total", "TotalNoFlavor", "TotalNoTime", "TotalNoFlavorNoTime" ]

        ### defining variables (not sure why here)
        self.nSub0p5 = ROOT.NsubjettinessWrapper( 0.5, 0.8, 0, 0 ) #beta, cone size, measureDef 0=Normalize, axesDef 0=KT_axes
        self.nSub1 = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 0 )
        self.nSub2 = ROOT.NsubjettinessWrapper( 2, 0.8, 0, 0 )

        self.nSub0p5_WTA_kT = ROOT.NsubjettinessWrapper( 0.5, 0.8, 0, 3 ) #beta, cone size, measureDef 0=Normalize, axesDef 0=WTA_kT_axes
        self.nSub1_WTA_kT = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 3 )
        self.nSub2_WTA_kT = ROOT.NsubjettinessWrapper( 2, 0.8, 0, 3 )

        self.nSub0p5_OP_kT = ROOT.NsubjettinessWrapper( 0.5, 0.8, 0, 6 ) #beta, cone size, measureDef 0=Normalize,axesDef 6=onepass_kT_axes
        self.nSub1_OP_kT = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 6 )
        self.nSub2_OP_kT = ROOT.NsubjettinessWrapper( 2, 0.8, 0, 6 )
        self.sd = ROOT.SoftDropWrapper(self.beta, self.zcut, self.R, self.minJetPt)

        print ("Load C++ Recluster worker module")
        ROOT.gSystem.Load("libPhysicsToolsNanoAODJMARTools.so")



    def beginJob(self, histFile, histDirName):
        Module.beginJob(self, histFile, histDirName)

	self.nSub_labels = ["_tau_0p5_", "_tau_1_", "_tau_2_"]
        ## JES from https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources#Main_uncertainties_2016_80X
        #self.JESLabels = [ "AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias", "Fragmentation", "SinglePionECAL", "SinglePionHCAL", "FlavorQCD", "TimePtEta", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeBal", "RelativeSample", "RelativeFSR", "RelativeStatFSR", "RelativeStatEC", "RelativeStatHF", "PileUpDataMC", "PileUpPtRef", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF", "PileUpMuZero", "PileUpEnvelope", "SubTotalPileUp", "SubTotalRelative", "SubTotalPt", "SubTotalScale", "SubTotalAbsolute", "SubTotalMC", "Total", "TotalNoFlavor", "TotalNoTime", "TotalNoFlavorNoTime" ]


        ### Booking histograms
        self.addObject( ROOT.TH1F('PUweight',   ';PUWeight',   20, 0, 2) )
        #### general selection
        for isel in [ '_noSelnoWeight', '_noSel', '_dijetSel' ]:
            self.addObject( ROOT.TH1F('nPVs'+isel,   ';number of PVs',   100, 0, 100) )
            self.addObject( ROOT.TH1F('nleps'+isel,   ';number of leptons',   20, 0, 20) )
            self.addP4Hists( 'muons', isel )
            self.addP4Hists( 'eles', isel )
            self.addObject( ROOT.TH1F('nAK8jets'+isel,   ';number of AK8 jets',   20, 0, 20) )
            self.addP4Hists( 'AK8jets', isel )
            self.addObject( ROOT.TH1F('METPt'+isel,   ';MET (GeV)',   200, 0, 2000) )
        self.addP4Hists( 'leadAK8jet', '_dijetSel' )

        for sysUnc in self.sysSource:
            self.addObject( ROOT.TH1F('recoJetPt'+sysUnc, ';AK8 reco jet pt [GeV]', 500, 0, 5000) )
            self.addObject( ROOT.TH1F('recoJetSDmass'+sysUnc, ';AK8 reco jet SD mass [GeV]', 500, 0, 500) )
            self.addObject( ROOT.TH1F('recoJetEta'+sysUnc, ';AK8 reco jet eta', 40, -5, 5) )
            self.addObject( ROOT.TH1F('recoJetTau21'+sysUnc, ';AK8 reco jet #tau_{21}', 40, 0, 1) )
            self.addObject( ROOT.TH1F('genJetPt'+sysUnc, ';AK8 gen jet pt [GeV]', 500, 0, 5000) )
            self.addObject( ROOT.TH1F('genJetEta'+sysUnc, ';AK8 gen jet eta', 40, -5, 5) )
            self.addObject( ROOT.TH1F('genJetTau21'+sysUnc, ';AK8 gen jet #tau_{21}', 40, 0, 1) )

            self.addObject( ROOT.TH2F('respJetTau21'+sysUnc, ';AK8 gen jet #tau_{21};AK8 reco jet #tau_{21}', 40, 0, 1, 40, 0, 1) )
            for tauN in range(self.maxTau):
                for x in self.nSub_labels:

                    self.addObject( ROOT.TH1F("recoJet"+x+str(tauN)+sysUnc, ";AK8 jet #tau", 100, 0, 1 ) )
                    self.addObject( ROOT.TH1F("genJet"+x+str(tauN)+sysUnc, ";AK8 jet #tau", 100, 0, 1 ) )
                    self.addObject( ROOT.TH2F('respJet'+x+str(tauN)+sysUnc, ';AK8 gen jet #tau_{21};AK8 reco jet #tau_{21}', 100, 0, 1, 100, 0, 1) )

                    self.addObject( ROOT.TH1F("recoJet"+x+str(tauN)+"_WTA_kT"+sysUnc, ";AK8 jet #tau", 100, 0, 1 ) )
                    self.addObject( ROOT.TH1F("genJet"+x+str(tauN)+"_WTA_kT"+sysUnc, ";AK8 jet #tau", 100, 0, 1 ) )
                    self.addObject( ROOT.TH2F('respJet'+x+str(tauN)+"_WTA_kT"+sysUnc, ';AK8 gen jet #tau_{21};AK8 reco jet #tau_{21}', 100, 0, 1, 100, 0, 1) )

                    self.addObject( ROOT.TH1F("recoJet"+x+str(tauN)+"_OP_kT"+sysUnc, ";AK8 jet #tau", 100, 0, 1 ) )
                    self.addObject( ROOT.TH1F("genJet"+x+str(tauN)+"_OP_kT"+sysUnc, ";AK8 jet #tau", 100, 0, 1 ) )
                    self.addObject( ROOT.TH2F('respJet'+x+str(tauN)+"_OP_kT"+sysUnc, ';AK8 gen jet #tau_{21};AK8 reco jet #tau_{21}', 100, 0, 1, 100, 0, 1) )

    def addP4Hists(self, s, t ):
        self.addObject( ROOT.TH1F(s+'_pt'+t,  s+';p_{T} (GeV)',   200, 0, 2000) )
        self.addObject( ROOT.TH1F(s+'_eta'+t, s+';#eta', 100, -4.0, 4.0 ) )
        self.addObject( ROOT.TH1F(s+'_phi'+t, s+';#phi', 100, -3.14259, 3.14159) )
        self.addObject( ROOT.TH1F(s+'_mass'+t,s+';mass (GeV)', 100, 0, 1000) )



    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        for ijet in ["0"]: #,1] :
            for x in self.kinematic_labels:
                self.out.branch("goodrecojet" + ijet + "%s"%x,  "F")
                if not x.startswith(('_eta', '_phi')):
                    for u in [ 'Up', 'Down' ]:
                        self.out.branch("goodrecojet" + ijet + "%s"%x+'_jer'+u,  "F")
                        for jes in self.JESLabels:
                            self.out.branch("goodrecojet" + ijet + "%s"%x+'_jesUnc'+jes+u,  "F")

            self.out.branch("ngoodrecojet" + ijet,  "I")  ### dummy for nanoAOD Tools
            self.out.branch("goodrecojet" + ijet + '_pt_raw',  "F")
            self.out.branch("goodrecojet" + ijet + '_mass_raw',  "F")
            self.out.branch("goodrecojet" + ijet + '_corr_JEC',  "F")
            self.out.branch("goodrecojet" + ijet + '_corr_JER',  "F")
            self.out.branch("goodrecojet" + ijet + '_corr_JMS',  "F")
            self.out.branch("goodrecojet" + ijet + '_corr_JMR',  "F")
            self.out.branch("goodrecojet" + ijet + '_mass_jmrUp',  "F")
            self.out.branch("goodrecojet" + ijet + '_mass_jmsUp',  "F")
            self.out.branch("goodrecojet" + ijet + '_mass_jmrDown',  "F")
            self.out.branch("goodrecojet" + ijet + '_mass_jmsDown',  "F")
            self.out.branch("goodrecojet" + ijet + "_softdrop_mass",  "F")
            self.out.branch("goodrecojet" + ijet + "_tau21",  "F")
            self.out.branch("goodrecojet" + ijet + "_N21",  "F")
            self.out.branch("goodrecojet" + ijet + "_nConst",  "I")

            for x in self.kinematic_labels:
                self.out.branch("goodgenjet" + ijet + "%s"%x,  "F")
            self.out.branch("goodgenjet" + ijet + "_tau21", "F")
            self.out.branch("goodgenjet" + ijet + "_nConst",  "I")

            self.out.branch("dR_gen_reco_AK8", "F")
            self.out.branch("dR_genW_genAK8", "F")
            self.out.branch("genEventNo_taus_are_0",  "I")
            self.out.branch("recoEventNo_taus_are_0",  "I")

            for x in self.kinematic_labels:
                self.out.branch("sd_goodrecojet" + ijet + "%s"%x,  "F")
                if not x.startswith(('_eta', '_phi')):
                    for u in [ 'Up', 'Down' ]:
                        self.out.branch("sd_goodrecojet" + ijet + "%s"%x+'_jer'+u,  "F")
                        for jes in self.JESLabels:
                            self.out.branch("sd_goodrecojet" + ijet + "%s"%x+'_jesUnc'+jes+u,  "F")
            self.out.branch("nsd_goodrecojet" + ijet,  "I")  ### dummy for nanoAOD Tools
            self.out.branch("sd_goodrecojet" + ijet + "_nConst",  "I")
            self.out.branch("sd_goodrecojet" + ijet + '_pt_raw',  "F")
            self.out.branch("sd_goodrecojet" + ijet + '_mass_raw',  "F")
            self.out.branch("sd_goodrecojet" + ijet + '_corr_JEC',  "F")
            self.out.branch("sd_goodrecojet" + ijet + '_corr_JER',  "F")
            self.out.branch("sd_goodrecojet" + ijet + '_corr_JMS',  "F")
            self.out.branch("sd_goodrecojet" + ijet + '_corr_JMR',  "F")
            self.out.branch("sd_goodrecojet" + ijet + '_mass_jmrUp',  "F")
            self.out.branch("sd_goodrecojet" + ijet + '_mass_jmsUp',  "F")
            self.out.branch("sd_goodrecojet" + ijet + '_mass_jmrDown',  "F")
            self.out.branch("sd_goodrecojet" + ijet + '_mass_jmsDown',  "F")

            for x in self.kinematic_labels:
                self.out.branch("sd_goodgenjet" + ijet + "%s"%x,  "F")
            self.out.branch("sd_goodgenjet" + ijet + "_nConst",  "I")


            for tauN in range(self.maxTau):
                for x in self.nSub_labels:

                    self.out.branch("goodrecojet" + ijet + "%s"%x +str(tauN),  "F")
                    self.out.branch("goodgenjet" + ijet + "%s"%x +str(tauN),  "F")

                    self.out.branch("goodrecojet" + ijet + "%s"%x +str(tauN) + "_WTA_kT",  "F")
                    self.out.branch("goodgenjet" + ijet + "%s"%x +str(tauN) + "_WTA_kT",  "F")

                    self.out.branch("goodrecojet" + ijet + "%s"%x +str(tauN) + "_OP_kT",  "F")
                    self.out.branch("goodgenjet" + ijet + "%s"%x +str(tauN) + "_OP_kT",  "F")

                    self.out.branch("sd_goodrecojet" + ijet + "%s"%x +str(tauN),  "F")
                    self.out.branch("sd_goodgenjet" + ijet + "%s"%x +str(tauN),  "F")

                    self.out.branch("sd_goodrecojet" + ijet + "%s"%x +str(tauN) + "_WTA_kT",  "F")
                    self.out.branch("sd_goodgenjet" + ijet + "%s"%x +str(tauN) + "_WTA_kT",  "F")

                    self.out.branch("sd_goodrecojet" + ijet + "%s"%x +str(tauN) + "_OP_kT",  "F")
                    self.out.branch("sd_goodgenjet" + ijet + "%s"%x +str(tauN) + "_OP_kT",  "F")

        pass


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        '''process event, return True (go to next module) or False (fail, go to next event)'''

        isMC = event.run == 1

        jets = list(Collection(event, 'FatJet' ))
        pfCands = list(Collection(event, 'PFCandsAK8' ))
        electrons = list(Collection(event, 'Electron'))
        muons = list(Collection(event, 'Muon'))
        bjets = list(Collection(event, 'Jet'))
        met = Object(event, 'MET')
        PV = Object(event, 'PV')


        #### Lepton cleaning
        recoElectrons  = [x for x in electrons if x.pt>self.minElectronPt and x.cutBased_HEEP and ((self.range1ElectronEta[0]<abs(x.p4().Eta())<self.range1ElectronEta[1]) or (self.range2ElectronEta[0]<abs(x.p4().Eta())<self.range2ElectronEta[1]))]
        #recoElectrons.sort(key=lambda x:x.pt, reverse=True)

        recoMuons = [ x for x in muons if x.pt > self.minMuonPt and x.highPtId > 1 and abs(x.p4().Eta()) < self.maxMuonEta and x.pfRelIso03_all < 0.1]
        #recoMuons.sort(key=lambda x:x.pt, reverse=True)
        nleptons = len(recoMuons)+len(recoElectrons)
        ##################################################

        #### MET (not sure if needed)
        MET = ROOT.TLorentzVector()
        MET.SetPtEtaPhiE(met.pt, 0., met.phi, met.sumEt)
        ##################################################

        #### Basic jet selection
        recojets = [ x for x in jets if x.pt_nom > self.minJetPt and abs(x.p4().Eta()) < self.maxJetEta]
        recojets.sort(key=lambda x:x.pt_nom,reverse=True)
        ##################################################

        #### Weight
        if isMC: weight = event.puWeight * event.genWeight
        else: weight = 1
        ##################################################

        #### Checking no selection without weights
        getattr( self, 'nPVs_noSelnoWeight' ).Fill( PV.npvsGood )
        getattr( self, 'nleps_noSelnoWeight' ).Fill( nleptons )
        for imuon in recoMuons:
            getattr( self, 'muons_pt_noSelnoWeight' ).Fill( imuon.pt )
            getattr( self, 'muons_eta_noSelnoWeight' ).Fill( imuon.eta )
            getattr( self, 'muons_phi_noSelnoWeight' ).Fill( imuon.phi )
        for iele in recoElectrons:
            getattr( self, 'eles_pt_noSelnoWeight' ).Fill( iele.pt )
            getattr( self, 'eles_eta_noSelnoWeight' ).Fill( iele.eta )
            getattr( self, 'eles_phi_noSelnoWeight' ).Fill( iele.phi )
        getattr( self, 'nAK8jets_noSelnoWeight' ).Fill( len(recojets) )
        for ijet in recojets:
            getattr( self, 'AK8jets_pt_noSelnoWeight' ).Fill( ijet.pt )
            getattr( self, 'AK8jets_eta_noSelnoWeight' ).Fill( ijet.eta )
            getattr( self, 'AK8jets_phi_noSelnoWeight' ).Fill( ijet.phi )
            getattr( self, 'AK8jets_mass_noSelnoWeight' ).Fill( ijet.msoftdrop )
        getattr( self, 'METPt_noSelnoWeight' ).Fill( MET.Pt() )

        #### Checking no selection with weights
        getattr( self, 'nPVs_noSel' ).Fill( PV.npvsGood, weight )
        getattr( self, 'nleps_noSel' ).Fill( nleptons, weight )
        for imuon in recoMuons:
            getattr( self, 'muons_pt_noSel' ).Fill( imuon.pt, weight )
            getattr( self, 'muons_eta_noSel' ).Fill( imuon.eta, weight )
            getattr( self, 'muons_phi_noSel' ).Fill( imuon.phi, weight )
        for iele in recoElectrons:
            getattr( self, 'eles_pt_noSel' ).Fill( iele.pt, weight )
            getattr( self, 'eles_eta_noSel' ).Fill( iele.eta, weight )
            getattr( self, 'eles_phi_noSel' ).Fill( iele.phi, weight )
        getattr( self, 'nAK8jets_noSel' ).Fill( len(recojets), weight )
        for ijet in recojets:
            getattr( self, 'AK8jets_pt_noSel' ).Fill( ijet.pt, weight )
            getattr( self, 'AK8jets_eta_noSel' ).Fill( ijet.eta, weight )
            getattr( self, 'AK8jets_phi_noSel' ).Fill( ijet.phi, weight )
            getattr( self, 'AK8jets_mass_noSel' ).Fill( ijet.msoftdrop, weight )
        getattr( self, 'METPt_noSel' ).Fill( MET.Pt(), weight )

        ##### Applying selection
        if (nleptons==0) and (len(recojets)>1):

            passAK8jet = {}          ### Storing goodreco jet as list for later use
            passAK8jet['jet'] = recojets[0]

            ##### Computing quantities
            if (recojets[0].pt > self.minLeadJetPt*0.8):    #### store jets in the range of the pt selection

                #### Run calculations of NSub bases and store for ungroomed recojets

                #### Applying PUPPI weights to the PF candidates
                constituents = ROOT.vector("TLorentzVector")()
                recoCandsPUPPIweightedVec = ROOT.vector("TLorentzVector")()
                for p in pfCands :
                    t = ROOT.TLorentzVector( p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
                    pw =  p.puppiWeight
                    tp = ROOT.TLorentzVector(p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
                    tp = tp * pw
                    #if pw > 0. :
                    #print "Applying PUPPI weights"
                    recoCandsPUPPIweightedVec.push_back(tp)

                #### Storing only the PF candidates that are close to the leadAK8jet (constituents)
                for x in recoCandsPUPPIweightedVec:
                    if abs(recojets[0].p4().DeltaR( x )) < 0.8: constituents.push_back(x)

                #### Computing n-subjetiness basis from PF PUPPI constituents
                nsub0p5 = self.nSub0p5.getTau( self.maxTau, constituents )
                nsub1 = self.nSub1.getTau( self.maxTau, constituents )
                nsub2 = self.nSub2.getTau( self.maxTau, constituents )

                nsub0p5_WTA_kT = self.nSub0p5_WTA_kT.getTau( self.maxTau, constituents )
                nsub1_WTA_kT = self.nSub1_WTA_kT.getTau( self.maxTau, constituents )
                nsub2_WTA_kT = self.nSub2_WTA_kT.getTau( self.maxTau, constituents )

                nsub0p5_OP_kT = self.nSub0p5_OP_kT.getTau( self.maxTau, constituents )
                nsub1_OP_kT = self.nSub1_OP_kT.getTau( self.maxTau, constituents )
                nsub2_OP_kT = self.nSub2_OP_kT.getTau( self.maxTau, constituents )

                try: passAK8jet['tau21'] = recojets[0].tau2/recojets[0].tau1
                except ValueError: passAK8jet['tau21'] = -1
                try: passAK8jet['tau32'] = recojets[0].tau3/recojets[0].tau2
                except ValueError: passAK8jet['tau32'] = -1

                #### Filling branch with passAK8jet info after selection
                self.out.fillBranch("ngoodrecojet0", 1 )  ### dummy for nanoAOD Tools
                self.out.fillBranch("goodrecojet0_pt",  recojets[0].pt_nom )
                self.out.fillBranch("goodrecojet0_eta",  recojets[0].p4().Eta() )
                self.out.fillBranch("goodrecojet0_phi",  recojets[0].p4().Phi() )
                self.out.fillBranch("goodrecojet0_mass",  recojets[0].p4().M() )
                self.out.fillBranch("goodrecojet0_softdrop_mass", recojets[0].msoftdrop_nom)
                self.out.fillBranch("goodrecojet0_pt_raw",  recojets[0].pt_raw )
                self.out.fillBranch("goodrecojet0_mass_raw",  recojets[0].mass_raw )
                self.out.fillBranch("goodrecojet0_tau21", passAK8jet['tau21'] )
                self.out.fillBranch("goodrecojet0_N21",  recojets[0].n2b1)
                self.out.fillBranch("goodrecojet0_corr_JEC",  recojets[0].corr_JEC )
                self.out.fillBranch("goodrecojet0_corr_JER",  ( recojets[0].corr_JER if isMC else -99999 ) )
                self.out.fillBranch("goodrecojet0_corr_JMS",  ( recojets[0].corr_JMS if isMC else -99999 ) )
                self.out.fillBranch("goodrecojet0_corr_JMR",  ( recojets[0].corr_JMR if isMC else -99999 ) )
                self.out.fillBranch("goodrecojet0_mass_jmrUp",  ( recojets[0].mass_jmrUp if isMC else -99999 ) )
                self.out.fillBranch("goodrecojet0_mass_jmsUp",  ( recojets[0].mass_jmsUp if isMC else -99999 ) )
                self.out.fillBranch("goodrecojet0_mass_jmrDown",  ( recojets[0].mass_jmrDown if isMC else -99999 ) )
                self.out.fillBranch("goodrecojet0_mass_jmsDown",  ( recojets[0].mass_jmsDown if isMC else -99999 ) )
                for q in ['pt', 'mass']:
                    for u in [ 'Up', 'Down' ]:
                        self.out.fillBranch("goodrecojet0_"+q+'_jer'+u, ( getattr( recojets[0], q+'_jer'+u ) if isMC else -99999 ) )
                        for jes in self.JESLabels:
                            self.out.fillBranch("goodrecojet0_"+q+'_jesUnc'+jes+u, ( getattr( recojets[0], q+'_jes'+jes+u ) if isMC else -99999 ) )


                #### filling histos and branches with nsub basis
                for tauN in range(self.maxTau):
                    if nsub0p5[tauN]==0. or nsub1[tauN]==0. or nsub2[tauN]==0.:
                        self.out.fillBranch("recoEventNo_taus_are_0", event.event)
                    self.out.fillBranch("goodrecojet0_tau_0p5_"+str(tauN),  nsub0p5[tauN]  )
                    self.out.fillBranch("goodrecojet0_tau_1_"+str(tauN),  nsub1[tauN]  )
                    self.out.fillBranch("goodrecojet0_tau_2_"+str(tauN),  nsub2[tauN]  )
                    passAK8jet['0p5'+str(tauN)] = nsub0p5[tauN]
                    passAK8jet['1'+str(tauN)] = nsub1[tauN]
                    passAK8jet['2'+str(tauN)] = nsub2[tauN]

                    self.out.fillBranch("goodrecojet0_tau_0p5_"+str(tauN) + "_WTA_kT",  nsub0p5_WTA_kT[tauN]  )
                    self.out.fillBranch("goodrecojet0_tau_1_"+str(tauN) + "_WTA_kT",  nsub1_WTA_kT[tauN]  )
                    self.out.fillBranch("goodrecojet0_tau_2_"+str(tauN) + "_WTA_kT",  nsub2_WTA_kT[tauN]  )
                    passAK8jet['0p5WTAkT'+str(tauN)] = nsub0p5_WTA_kT[tauN]
                    passAK8jet['1WTAkT'+str(tauN)] = nsub1_WTA_kT[tauN]
                    passAK8jet['2WTAkT'+str(tauN)] = nsub2_WTA_kT[tauN]

                    self.out.fillBranch("goodrecojet0_tau_0p5_"+str(tauN) + "_OP_kT",  nsub0p5_OP_kT[tauN]  )
                    self.out.fillBranch("goodrecojet0_tau_1_"+str(tauN) + "_OP_kT",  nsub1_OP_kT[tauN]  )
                    self.out.fillBranch("goodrecojet0_tau_2_"+str(tauN) + "_OP_kT",  nsub2_OP_kT[tauN]  )
                    passAK8jet['0p5OPkT'+str(tauN)] = nsub0p5_OP_kT[tauN]
                    passAK8jet['1OPkT'+str(tauN)] = nsub1_OP_kT[tauN]
                    passAK8jet['2OPkT'+str(tauN)] = nsub2_OP_kT[tauN]


                #### For MC only, loop over the genJet collection
                if isMC:
                    genParticles = Collection(event, 'GenPartAK8')
                    genPartAll = Collection(event, 'GenPart')
                    genjets  = list(Collection(event, 'GenJetAK8' ))

                    #### Matching recojet and genJet
                    tmpGenInd = -1
                    tmpGenDeltaR = 99999
                    for i in xrange(0,len(genjets)):
                        genDeltaR = abs(genjets[i].p4().DeltaR(recojets[0].p4()))
                        if genDeltaR < tmpGenDeltaR:
                            tmpGenDeltaR = genDeltaR
                            tmpGenInd = i

                    goodgenjet = genjets[tmpGenInd]
                    passAK8jet['genjet'] = genjets[tmpGenInd]

                    #### Checking with GenParticles are constituents of the matched genjet
                    genConstituents = ROOT.vector("TLorentzVector")()
                    for x in genParticles:
                        if abs(passAK8jet['genjet'].p4().DeltaR( x.p4() )) < 0.8: genConstituents.push_back(x.p4())

                    #### Computing nsub basis for genJets
                    genNsub0p5 = self.nSub0p5.getTau( self.maxTau, genConstituents )
                    genNsub1 = self.nSub1.getTau( self.maxTau, genConstituents )
                    genNsub2 = self.nSub2.getTau( self.maxTau, genConstituents )

                    genNsub0p5_WTA_kT = self.nSub0p5_WTA_kT.getTau( self.maxTau, genConstituents )
                    genNsub1_WTA_kT = self.nSub1_WTA_kT.getTau( self.maxTau, genConstituents )
                    genNsub2_WTA_kT = self.nSub2_WTA_kT.getTau( self.maxTau, genConstituents )

                    genNsub0p5_OP_kT = self.nSub0p5_OP_kT.getTau( self.maxTau, genConstituents )
                    genNsub1_OP_kT = self.nSub1_OP_kT.getTau( self.maxTau, genConstituents )
                    genNsub2_OP_kT = self.nSub2_OP_kT.getTau( self.maxTau, genConstituents )

                    self.out.fillBranch("goodgenjet0_pt",  passAK8jet['genjet'].pt)
                    self.out.fillBranch("goodgenjet0_eta",  passAK8jet['genjet'].p4().Eta() )
                    self.out.fillBranch("goodgenjet0_phi",  passAK8jet['genjet'].p4().Phi() )
                    self.out.fillBranch("goodgenjet0_mass",  passAK8jet['genjet'].p4().M() )
                    self.out.fillBranch("dR_gen_reco_AK8",  tmpGenDeltaR)

                    try: passAK8jet['gentau21'] = genNsub1[1]/genNsub1[0]
                    except ValueError: passAK8jet['gentau21'] = -1
                    try: passAK8jet['gentau32'] = genNsub1[2]/genNsub1[1]
                    except ValueError: passAK8jet['gentau32'] = -1
                    self.out.fillBranch("goodgenjet0_tau21", passAK8jet['gentau21'] )

                    ###### filling histos

                    for tauN in range(self.maxTau):
                        self.out.fillBranch("goodgenjet0_tau_0p5_"+str(tauN),  genNsub0p5[tauN]  )
                        self.out.fillBranch("goodgenjet0_tau_1_"+str(tauN),  genNsub1[tauN]  )
                        self.out.fillBranch("goodgenjet0_tau_2_"+str(tauN),  genNsub2[tauN]  )
                        passAK8jet['gen0p5'+str(tauN)] = genNsub0p5[tauN]
                        passAK8jet['gen1'+str(tauN)] = genNsub1[tauN]
                        passAK8jet['gen2'+str(tauN)] = genNsub2[tauN]

                        self.out.fillBranch("goodgenjet0_tau_0p5_"+str(tauN) + "_WTA_kT",  genNsub0p5_WTA_kT[tauN]  )
                        self.out.fillBranch("goodgenjet0_tau_1_"+str(tauN) + "_WTA_kT",  genNsub1_WTA_kT[tauN]  )
                        self.out.fillBranch("goodgenjet0_tau_2_"+str(tauN) + "_WTA_kT",  genNsub2_WTA_kT[tauN]  )
                        passAK8jet['gen0p5WTAkT'+str(tauN)] = genNsub0p5_WTA_kT[tauN]
                        passAK8jet['gen1WTAkT'+str(tauN)] = genNsub1_WTA_kT[tauN]
                        passAK8jet['gen2WTAkT'+str(tauN)] = genNsub2_WTA_kT[tauN]

                        self.out.fillBranch("goodgenjet0_tau_0p5_"+str(tauN) + "_OP_kT",  genNsub0p5_OP_kT[tauN]  )
                        self.out.fillBranch("goodgenjet0_tau_1_"+str(tauN) + "_OP_kT",  genNsub1_OP_kT[tauN]  )
                        self.out.fillBranch("goodgenjet0_tau_2_"+str(tauN) + "_OP_kT",  genNsub2_OP_kT[tauN]  )
                        passAK8jet['gen0p5OPkT'+str(tauN)] = genNsub0p5_OP_kT[tauN]
                        passAK8jet['gen1OPkT'+str(tauN)] = genNsub1_OP_kT[tauN]
                        passAK8jet['gen2OPkT'+str(tauN)] = genNsub2_OP_kT[tauN]
            ##################### end of computing quantities

            ##### FIlling histograms
            for sysUnc in self.sysSource:
                if( getattr(passAK8jet['jet'], 'pt'+sysUnc) > self.minLeadJetPt ):

                    #### Checking selection with weights
                    if sysUnc.endswith('nom'):
                        getattr( self, 'nPVs_dijetSel' ).Fill( PV.npvsGood, weight )
                        getattr( self, 'nleps_dijetSel' ).Fill( nleptons, weight )
                        for imuon in recoMuons:
                            getattr( self, 'muons_pt_dijetSel' ).Fill( imuon.pt, weight )
                            getattr( self, 'muons_eta_dijetSel' ).Fill( imuon.eta, weight )
                            getattr( self, 'muons_phi_dijetSel' ).Fill( imuon.phi, weight )
                        for iele in recoElectrons:
                            getattr( self, 'eles_pt_dijetSel' ).Fill( iele.pt, weight )
                            getattr( self, 'eles_eta_dijetSel' ).Fill( iele.eta, weight )
                            getattr( self, 'eles_phi_dijetSel' ).Fill( iele.phi, weight )
                        getattr( self, 'nAK8jets_dijetSel' ).Fill( len(recojets), weight )
                        for ijet in recojets:
                            getattr( self, 'AK8jets_pt_dijetSel' ).Fill( ijet.pt, weight )
                            getattr( self, 'AK8jets_eta_dijetSel' ).Fill( ijet.eta, weight )
                            getattr( self, 'AK8jets_phi_dijetSel' ).Fill( ijet.phi, weight )
                            getattr( self, 'AK8jets_mass_dijetSel' ).Fill( ijet.msoftdrop, weight )
                        getattr( self, 'METPt_dijetSel' ).Fill( MET.Pt(), weight )
                        getattr( self, 'leadAK8jet_pt_dijetSel' ).Fill( getattr(passAK8jet['jet'], 'pt'), weight )
                        getattr( self, 'leadAK8jet_eta_dijetSel' ).Fill( getattr(passAK8jet['jet'], 'eta'), weight )
                        getattr( self, 'leadAK8jet_phi_dijetSel' ).Fill( getattr(passAK8jet['jet'], 'phi'), weight )
                        getattr( self, 'leadAK8jet_mass_dijetSel' ).Fill( getattr(passAK8jet['jet'], 'msoftdrop'), weight )
                        if isMC:
                            getattr( self, 'genJetPt'+sysUnc ).Fill( passAK8jet['genjet'].pt, weight )
                            getattr( self, 'genJetEta'+sysUnc ).Fill( passAK8jet['genjet'].eta, weight )
                            getattr( self, 'genJetTau21'+sysUnc ).Fill( passAK8jet['gentau21'], weight )
                            for tauN in range(self.maxTau):
                                getattr( self, "genJet_tau_0p5_"+str(tauN)+sysUnc ).Fill( passAK8jet['gen0p5'+str(tauN)], weight )
                                getattr( self, "genJet_tau_1_"+str(tauN)+sysUnc ).Fill( passAK8jet['gen1'+str(tauN)], weight )
                                getattr( self, "genJet_tau_2_"+str(tauN)+sysUnc ).Fill( passAK8jet['gen2'+str(tauN)], weight )
                                getattr( self, "genJet_tau_0p5_"+str(tauN)+'_WTA_kT'+sysUnc ).Fill( passAK8jet['gen0p5WTAkT'+str(tauN)], weight )
                                getattr( self, "genJet_tau_1_"+str(tauN)+'_WTA_kT'+sysUnc ).Fill( passAK8jet['gen1WTAkT'+str(tauN)], weight )
                                getattr( self, "genJet_tau_2_"+str(tauN)+'_WTA_kT'+sysUnc ).Fill( passAK8jet['gen2WTAkT'+str(tauN)], weight )
                                getattr( self, "genJet_tau_0p5_"+str(tauN)+'_OP_kT'+sysUnc ).Fill( passAK8jet['gen0p5OPkT'+str(tauN)], weight )
                                getattr( self, "genJet_tau_1_"+str(tauN)+'_OP_kT'+sysUnc ).Fill( passAK8jet['gen1OPkT'+str(tauN)], weight )
                                getattr( self, "genJet_tau_2_"+str(tauN)+'_OP_kT'+sysUnc ).Fill( passAK8jet['gen2OPkT'+str(tauN)], weight )

                    #### filling histos
                    getattr( self, 'recoJetPt'+sysUnc ).Fill( getattr(passAK8jet['jet'], 'pt'+sysUnc ), weight )
                    getattr( self, 'recoJetEta'+sysUnc ).Fill( getattr(passAK8jet['jet'], 'eta'), weight )
                    getattr( self, 'recoJetSDmass'+sysUnc ).Fill( getattr(passAK8jet['jet'], 'msoftdrop'+sysUnc ), weight )
                    getattr( self, 'recoJetTau21'+sysUnc ).Fill( passAK8jet['tau21'], weight )
                    getattr( self, 'respJetTau21'+sysUnc ).Fill( passAK8jet['gentau21'], passAK8jet['tau21'], weight )
                    for tauN in range(self.maxTau):
                        getattr( self, "recoJet_tau_0p5_"+str(tauN)+sysUnc ).Fill( passAK8jet['0p5'+str(tauN)], weight )
                        getattr( self, 'respJet_tau_0p5_'+str(tauN)+sysUnc ).Fill( passAK8jet['gen0p5'+str(tauN)], passAK8jet['0p5'+str(tauN)], weight )
                        getattr( self, "recoJet_tau_1_"+str(tauN)+sysUnc ).Fill( passAK8jet['1'+str(tauN)], weight )
                        getattr( self, 'respJet_tau_1_'+str(tauN)+sysUnc ).Fill( passAK8jet['gen1'+str(tauN)], passAK8jet['1'+str(tauN)], weight )
                        getattr( self, "recoJet_tau_2_"+str(tauN)+sysUnc ).Fill( passAK8jet['2'+str(tauN)], weight )
                        getattr( self, 'respJet_tau_2_'+str(tauN)+sysUnc ).Fill( passAK8jet['gen2'+str(tauN)], passAK8jet['2'+str(tauN)], weight )
                        getattr( self, "recoJet_tau_0p5_"+str(tauN)+'_WTA_kT'+sysUnc ).Fill( passAK8jet['0p5WTAkT'+str(tauN)], weight )
                        getattr( self, 'respJet_tau_0p5_'+str(tauN)+'_WTA_kT'+sysUnc ).Fill( passAK8jet['gen0p5WTAkT'+str(tauN)], passAK8jet['0p5WTAkT'+str(tauN)], weight )
                        getattr( self, "recoJet_tau_1_"+str(tauN)+'_WTA_kT'+sysUnc ).Fill( passAK8jet['1WTAkT'+str(tauN)], weight )
                        getattr( self, 'respJet_tau_1_'+str(tauN)+'_WTA_kT'+sysUnc ).Fill( passAK8jet['gen1WTAkT'+str(tauN)], passAK8jet['1WTAkT'+str(tauN)], weight )
                        getattr( self, "recoJet_tau_2_"+str(tauN)+'_WTA_kT'+sysUnc ).Fill( passAK8jet['2WTAkT'+str(tauN)], weight )
                        getattr( self, 'respJet_tau_2_'+str(tauN)+'_WTA_kT'+sysUnc ).Fill( passAK8jet['gen2WTAkT'+str(tauN)], passAK8jet['2WTAkT'+str(tauN)], weight )
                        getattr( self, "recoJet_tau_0p5_"+str(tauN)+'_OP_kT'+sysUnc ).Fill( passAK8jet['0p5OPkT'+str(tauN)], weight )
                        getattr( self, 'respJet_tau_0p5_'+str(tauN)+'_OP_kT'+sysUnc ).Fill( passAK8jet['gen0p5OPkT'+str(tauN)], passAK8jet['0p5OPkT'+str(tauN)], weight )
                        getattr( self, "recoJet_tau_1_"+str(tauN)+'_OP_kT'+sysUnc ).Fill( passAK8jet['1OPkT'+str(tauN)], weight )
                        getattr( self, 'respJet_tau_1_'+str(tauN)+'_OP_kT'+sysUnc ).Fill( passAK8jet['gen1OPkT'+str(tauN)], passAK8jet['1OPkT'+str(tauN)], weight )
                        getattr( self, "recoJet_tau_2_"+str(tauN)+'_OP_kT'+sysUnc ).Fill( passAK8jet['2OPkT'+str(tauN)], weight )
                        getattr( self, 'respJet_tau_2_'+str(tauN)+'_OP_kT'+sysUnc ).Fill( passAK8jet['gen2OPkT'+str(tauN)], passAK8jet['2OPkT'+str(tauN)], weight )



#        #### Computing Softdrop jets
#        sdrecojets = self.sd.result( recoCandsPUPPIweightedVec )
#        ##################################################
#
#
#        ### Storing goodreco jet as list for later use
#        goodsdrecojet = []
#        #Run calculations of NSub bases and store for groomed recojets
#        if len(sdrecojets)>0:
#            sdrecojet = sdrecojets[0]
#            #### Checking which recojet is matched to the sdrecojets
#            sdrecojetFull = False   ## dummy
#            tmpsdrecojet = ROOT.TLorentzVector( )
#            tmpsdrecojet.SetPtEtaPhiM( sdrecojet.perp(), sdrecojet.eta(), sdrecojet.phi(), sdrecojet.m() )
#            deltaR = 99999
#            for ireco in jets:  ### to be safe, loop over all the jets
#                tmpDeltaR = tmpsdrecojet.DeltaR( ireco.p4() )
#                if (tmpDeltaR < deltaR):
#                    deltaR=tmpDeltaR
#                    if (tmpDeltaR<0.3): sdrecojetFull = ireco
#
#            # Cluster only the particles near the appropriate jet to save time
#            sd_constituents = ROOT.vector("TLorentzVector")()
#
#            for x in sdrecojet.constituents():
#                sd_constits = ROOT.TLorentzVector( x.px(), x.py(), x.pz(), x.E())
#                if abs(sdrecojet.delta_R( x )) < 0.8:
#                    sd_constituents.push_back(sd_constits)
#            sd_nsub0p5 = self.nSub0p5.getTau( self.maxTau, sd_constituents )
#            sd_nsub1 = self.nSub1.getTau( self.maxTau, sd_constituents )
#            sd_nsub2 = self.nSub2.getTau( self.maxTau, sd_constituents )
#
#            sd_nsub0p5_WTA_kT = self.nSub0p5_WTA_kT.getTau( self.maxTau, sd_constituents )
#            sd_nsub1_WTA_kT = self.nSub1_WTA_kT.getTau( self.maxTau, sd_constituents )
#            sd_nsub2_WTA_kT = self.nSub2_WTA_kT.getTau( self.maxTau, sd_constituents )
#
#            sd_nsub0p5_OP_kT = self.nSub0p5_OP_kT.getTau( self.maxTau, sd_constituents )
#            sd_nsub1_OP_kT = self.nSub1_OP_kT.getTau( self.maxTau, sd_constituents )
#            sd_nsub2_OP_kT = self.nSub2_OP_kT.getTau( self.maxTau, sd_constituents )
#
#            self.out.fillBranch("nsd_goodrecojet0", 1 )  ## dummy
#            self.out.fillBranch("sd_goodrecojet0_pt",  sdrecojet.perp() )
#            self.out.fillBranch("sd_goodrecojet0_eta",  sdrecojet.eta() )
#            self.out.fillBranch("sd_goodrecojet0_phi",  sdrecojet.phi() )
#            self.out.fillBranch("sd_goodrecojet0_mass",  sdrecojet.m() )
#            self.out.fillBranch("sd_goodrecojet0_pt_raw",  sdrecojetFull.pt_raw )
#            self.out.fillBranch("sd_goodrecojet0_mass_raw",  sdrecojetFull.mass_raw )
#            self.out.fillBranch("sd_goodrecojet0_corr_JEC",  sdrecojetFull.corr_JEC )
#            if isMC:
#                self.out.fillBranch("sd_goodrecojet0_corr_JER",  sdrecojetFull.corr_JER )
#                self.out.fillBranch("sd_goodrecojet0_corr_JMS",  sdrecojetFull.corr_JMS )
#                self.out.fillBranch("sd_goodrecojet0_corr_JMR",  sdrecojetFull.corr_JMR )
#                self.out.fillBranch("sd_goodrecojet0_mass_jmrUp",  sdrecojetFull.mass_jmrUp )
#                self.out.fillBranch("sd_goodrecojet0_mass_jmsUp",  sdrecojetFull.mass_jmsUp )
#                self.out.fillBranch("sd_goodrecojet0_mass_jmrDown",  sdrecojetFull.mass_jmrDown )
#                self.out.fillBranch("sd_goodrecojet0_mass_jmsDown",  sdrecojetFull.mass_jmsDown )
#                for q in ['pt', 'mass']:
#                    for u in [ 'Up', 'Down' ]:
#                        self.out.fillBranch("sd_goodrecojet0_"+q+'_jer'+u,  getattr( sdrecojetFull, q+'_jer'+u ) )
#                        for jes in self.JESLabels:
#                            self.out.fillBranch("sd_goodrecojet0_"+q+'_jesUnc'+jes+u,  getattr( sdrecojetFull, q+'_jes'+jes+u ) )
#
#            for tauN in range(self.maxTau):
#
#                self.out.fillBranch("sd_goodrecojet0_tau_0p5_"+str(tauN),  sd_nsub0p5[tauN]  )
#                self.out.fillBranch("sd_goodrecojet0_tau_1_"+str(tauN),  sd_nsub1[tauN]  )
#                self.out.fillBranch("sd_goodrecojet0_tau_2_"+str(tauN),  sd_nsub2[tauN]  )
#
#                self.out.fillBranch("sd_goodrecojet0_tau_0p5_"+str(tauN) + "_WTA_kT",  sd_nsub0p5_WTA_kT[tauN]  )
#                self.out.fillBranch("sd_goodrecojet0_tau_1_"+str(tauN) + "_WTA_kT",  sd_nsub1_WTA_kT[tauN]  )
#                self.out.fillBranch("sd_goodrecojet0_tau_2_"+str(tauN) + "_WTA_kT",  sd_nsub2_WTA_kT[tauN]  )
#
#                self.out.fillBranch("sd_goodrecojet0_tau_0p5_"+str(tauN) + "_OP_kT",  sd_nsub0p5_OP_kT[tauN]  )
#                self.out.fillBranch("sd_goodrecojet0_tau_1_"+str(tauN) + "_OP_kT",  sd_nsub1_OP_kT[tauN]  )
#                self.out.fillBranch("sd_goodrecojet0_tau_2_"+str(tauN) + "_OP_kT",  sd_nsub2_OP_kT[tauN]  )
#
#
#        ### ALE: move this part here, not need to compute this if does not pass selection
#        if isMC:
#            self.realW = 0
#            self.isgenW=0
#            genParticles = Collection(event, 'GenPartAK8')
#            genPartAll = Collection(event, 'GenPart')
#            gen_ak8 = list(Collection(event, 'GenJetAK8' ))
#
#
#            #list real W's, match gen-level hadronically daughting W's to selected W candidate jets at reco-level
#            genWmoms = []
#            genWdaughts = []
#
#
#            Tmom = [x for x in genPartAll if x.pt>10 and abs(x.pdgId)==6]
#            Tdaughts = [x for x in genPartAll if x.pt>1 and (abs(x.pdgId)==5 or abs(x.pdgId)==24)]
#
#            Wdaughts = [x for x in genPartAll if x.pt>1. and 0<abs(x.pdgId)<9]
#
#            Wmom = [x for x in genPartAll if x.pt>10. and abs(x.pdgId)==24]
#
#
#            if len(Wdaughts)>0 and len(Wmom)>0 and len(Tdaughts)>0 and len(Tmom)>0 :
#                for x in Wdaughts:
#                    for y in Wmom:
#                        try:
#                            if y==Wmom[x.genPartIdxMother] and (y in Tdaughts):
#                                self.out.fillBranch("PID_W",y.pdgId)
#                                #print y, "W from top decaying to qq"
#                                #self.out.fillBranch("GPIdx_W",x.genPartIdxMother)
#                                genWmoms.append(y)
#                                genWdaughts.append(x)
#
#                        except:
#                            continue
#
#
#            genjets = [x for x in gen_ak8 if x.pt > 0.8*self.minJetPt and  abs(x.p4().Eta()) < self.maxJetEta]   ### ALE: not sure that we want to apply a selection in genJets. What happends if you remove this?
#            genjets.sort(key=lambda x:x.pt,reverse=True)
#
#
#            if len(genjets)<1:
#                print "#exit if no gen-jets in event", event.event
#                return False #exit if no gen-jets in event
#
#
#            dRmin=[0.2,0]
#
#            for i in xrange(0,len(genjets)):
#                genjet_4v = ROOT.TLorentzVector()
#                genjet_4v.SetPtEtaPhiM(genjets[i].pt, genjets[i].eta, genjets[i].phi, genjets[i].mass)
#                if abs(genjet_4v.DeltaR(recoAK8))<dRmin[0]:
#                    dRmin[0] = abs(genjet_4v.DeltaR(recoAK8))
#                    dRmin[1] = i
#
#            goodgenjet = genjets[dRmin[1]]
#            genlevelAK8 = ROOT.TLorentzVector()
#            genlevelAK8.SetPtEtaPhiM(goodgenjet.pt, goodgenjet.eta, goodgenjet.phi, goodgenjet.mass)
#
#            ###Matching recojets to gen-level W's and W daughters
#            for W in genWmoms:
#                genW_4v = ROOT.TLorentzVector()
#                genW_4v.SetPtEtaPhiM(W.pt,W.eta,W.phi,W.mass)
#
#                if abs(recoAK8.DeltaR(genW_4v))<0.8:
#
#                    ndec = 0
#                    for Wdec in genWdaughts:
#
#                        gendec_4v = ROOT.TLorentzVector()
#                        gendec_4v.SetPtEtaPhiM(Wdec.pt,Wdec.eta,Wdec.phi,Wdec.mass)
#
#                        if abs(recoAK8.DeltaR(gendec_4v))<0.6:
#                            ndec +=1
#
#                    if ndec>1:
#                        self.out.fillBranch("dR_AK8motherW_mW", abs(recoAK8.DeltaR(genW_4v)))
#                        self.out.fillBranch("dR_AK8qi_mW", abs(recoAK8.DeltaR(gendec_4v)))
#                        self.out.fillBranch("PID_mW", W.pdgId)
#                        self.realW = 1
#                    else:
#                        self.out.fillBranch("PID_umW", W.pdgId)
#                        self.realW = 0
#
#            for W in genWmoms:
#                genW_4v = ROOT.TLorentzVector()
#                genW_4v.SetPtEtaPhiM(W.pt,W.eta,W.phi,W.mass)
#
#                if abs(genlevelAK8.DeltaR(genW_4v))<0.8:
#
#                    ndec = 0
#                    for Wdec in genWdaughts:
#                        gendec_4v = ROOT.TLorentzVector()
#                        gendec_4v.SetPtEtaPhiM(Wdec.pt,Wdec.eta,Wdec.phi,Wdec.mass)
#
#                        if abs(genlevelAK8.DeltaR(gendec_4v))<0.6:
#                            ndec +=1
#
#                    if ndec>1:
#                        #self.out.fillBranch("genmatchedgenAK8", 1)
#                        self.isgenW=1
#                        self.out.fillBranch("dR_genW_genAK8", abs(genlevelAK8.DeltaR(genW_4v)))
#                    else:
#                        #self.out.fillBranch("genmatchedgenAK8", 0)
#                        self.isgenW=0
#
#            genpart_CandsVec = ROOT.vector("TLorentzVector")()
#
#            #Store the mMDT groomed (leading) genjet
#            for p in genParticles :
#                genpart_CandsVec.push_back( ROOT.TLorentzVector( p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E()) )
#
#            sdgenjets = self.sd.result( genpart_CandsVec )
#
#            if len(genpart_CandsVec) == 0 :
#                return False
#
#            #Run calculations of NSub bases and store for ungroomed genjets
#            # Cluster only the particles near the appropriate jet to save time
#            constituents = ROOT.vector("TLorentzVector")()
#
#
#            #Run calculations of NSub bases and store for groomed genjets
#            for i_sd_genjet,sd_genjet in enumerate(sdgenjets):
#
#                # Cluster only the particles near the appropriate jet to save time
#                constituents = ROOT.vector("TLorentzVector")()
#
#                for x in sd_genjet.constituents():
#
#                    sd_constits = ROOT.TLorentzVector( x.px(), x.py(), x.pz(), x.E())
#                    if abs( sd_genjet.delta_R( x )) < 0.8:
#                        constituents.push_back(sd_constits)
#
#                nsub0p5 = self.nSub0p5.getTau( self.maxTau, constituents )
#                nsub1 = self.nSub1.getTau( self.maxTau, constituents )
#                nsub2 = self.nSub2.getTau( self.maxTau, constituents )
#
#                nsub0p5_WTA_kT = self.nSub0p5_WTA_kT.getTau( self.maxTau, constituents )
#                nsub1_WTA_kT = self.nSub1_WTA_kT.getTau( self.maxTau, constituents )
#                nsub2_WTA_kT = self.nSub2_WTA_kT.getTau( self.maxTau, constituents )
#
#                nsub0p5_OP_kT = self.nSub0p5_OP_kT.getTau( self.maxTau, constituents )
#                nsub1_OP_kT = self.nSub1_OP_kT.getTau( self.maxTau, constituents )
#                nsub2_OP_kT = self.nSub2_OP_kT.getTau( self.maxTau, constituents )
#
#                if (i_sd_genjet < 1 ): #to extract only the leading jet
#
#                    self.out.fillBranch("sd_goodgenjet" + str(i_sd_genjet) + "_pt",   sd_genjet.perp() )
#                    self.out.fillBranch("sd_goodgenjet" + str(i_sd_genjet) + "_eta",  sd_genjet.eta() )
#                    self.out.fillBranch("sd_goodgenjet" + str(i_sd_genjet) + "_phi",  sd_genjet.phi() )
#                    self.out.fillBranch("sd_goodgenjet" + str(i_sd_genjet) + "_mass", sd_genjet.m() )
#
#                    for tauN in range(self.maxTau):
#
#                        self.out.fillBranch("sd_goodgenjet" + str(i_sd_genjet) + "_tau_0p5_"+str(tauN),  nsub0p5[tauN]  )
#                        self.out.fillBranch("sd_goodgenjet" + str(i_sd_genjet) + "_tau_1_"+str(tauN),  nsub1[tauN]  )
#                        self.out.fillBranch("sd_goodgenjet" + str(i_sd_genjet) + "_tau_2_"+str(tauN),  nsub2[tauN]  )
#
#                        self.out.fillBranch("sd_goodgenjet" + str(i_sd_genjet) + "_tau_0p5_"+str(tauN) + "_WTA_kT",  nsub0p5_WTA_kT[tauN]  )
#                        self.out.fillBranch("sd_goodgenjet" + str(i_sd_genjet) + "_tau_1_"+str(tauN) + "_WTA_kT",  nsub1_WTA_kT[tauN]  )
#                        self.out.fillBranch("sd_goodgenjet" + str(i_sd_genjet) + "_tau_2_"+str(tauN) + "_WTA_kT",  nsub2_WTA_kT[tauN]  )
#
#                        self.out.fillBranch("sd_goodgenjet" + str(i_sd_genjet) + "_tau_0p5_"+str(tauN) + "_OP_kT",  nsub0p5_OP_kT[tauN]  )
#                        self.out.fillBranch("sd_goodgenjet" + str(i_sd_genjet) + "_tau_1_"+str(tauN) + "_OP_kT",  nsub1_OP_kT[tauN]  )
#                        self.out.fillBranch("sd_goodgenjet" + str(i_sd_genjet) + "_tau_2_"+str(tauN) + "_OP_kT",  nsub2_OP_kT[tauN]  )
#
#

        return True

    ##################### Helpful functions
    def printP4( self, c ):
        if hasattr( c, "p4"):
            s = ' %6.2f %5.2f %5.2f %6.2f ' % ( c.p4().Perp(), c.p4().Eta(), c.p4().Phi(), c.p4().M() )
        elif hasattr( c, "Perp"):
            s = ' %6.2f %5.2f %5.2f %6.2f ' % ( c.Perp(), c.Eta(), c.Phi(), c.M() )
        else:
            s = ' %6.2f %5.2f %5.2f %6.2f ' % ( c.perp(), c.eta(), c.phi(), c.m() )
        return s
    def printCollection(self,coll):
        for ic,c in enumerate(coll):
            s = self.printP4( c )
            print ' %3d : %s' % ( ic, s )

    def getSubjets(self, p4, subjets, dRmax=0.8):
        ret = []
        for subjet in subjets :
            if p4.DeltaR(subjet.p4()) < dRmax and len(ret) < 2 :
                ret.append(subjet.p4())
        return ret

