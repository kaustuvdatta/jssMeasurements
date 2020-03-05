import ROOT
import math, os, sys
import numpy as np
import fastjet
ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *

class nSubProd(Module):
    def __init__(self, selection='dijet', sysSource=[]):
        self.writeHistFile=True

        self.selection = selection

        ### Kinematics Cuts AK8Jets ###
        self.minAK8JetPt = 170  ### this is the basic minimum, not the final
        self.maxJetAK8Eta = 2.4

        ### Kinematics Cuts Jets ###
        self.minJetPt = 30
        self.maxJetEta = 2.4
        self.minBDisc = 0.3093  ### L: 0.0614, M: 0.3093, T: 07221

        ### Kinenatic Cuts Muons ###
        self.minMuonPt = 20.
        self.maxMuonEta = 2.4

        ### Kinenatic Cuts Electrons ###
        self.minElectronPt = 35.
        self.range1ElectronEta = [0,1.442]
        self.range2ElectronEta = [1.56,2.5]

        ### Defining nsubjetiness basis
        self.maxTau = 10
	self.nSub_labels = ["_tau_0p5_", "_tau_1_", "_tau_2_"]
        self.nSub0p5 = ROOT.NsubjettinessWrapper( 0.5, 0.8, 0, 0 ) #beta, cone size, measureDef 0=Normalize, axesDef 0=KT_axes
        self.nSub1 = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 0 )
        self.nSub2 = ROOT.NsubjettinessWrapper( 2, 0.8, 0, 0 )

        self.nSub0p5_WTA_kT = ROOT.NsubjettinessWrapper( 0.5, 0.8, 0, 3 ) #beta, cone size, measureDef 0=Normalize, axesDef 0=WTA_kT_axes
        self.nSub1_WTA_kT = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 3 )
        self.nSub2_WTA_kT = ROOT.NsubjettinessWrapper( 2, 0.8, 0, 3 )

        self.nSub0p5_OP_kT = ROOT.NsubjettinessWrapper( 0.5, 0.8, 0, 6 ) #beta, cone size, measureDef 0=Normalize,axesDef 6=onepass_kT_axes
        self.nSub1_OP_kT = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 6 )
        self.nSub2_OP_kT = ROOT.NsubjettinessWrapper( 2, 0.8, 0, 6 )

        ### Softdrop quantities
        self.beta = 0.0
        self.zcut = 0.1
        self.R = 0.8
        self.sd = ROOT.SoftDropWrapper(self.beta, self.zcut, self.R, self.minAK8JetPt)

        print ("Load C++ Recluster worker module")
        ROOT.gSystem.Load("libPhysicsToolsNanoAODJMARTools.so")

        ### Helpers
	self.kinematic_labels = ["_pt", "_eta", "_phi", "_mass"]

        ### Uncerstinties
        self.sysSource = ['_nom'] + [ isys+i for i in [ 'Up', 'Down' ] for isys in sysSource if not isys.endswith('nom') ]
        ## JES from https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources#Main_uncertainties_2016_80X
        self.JESLabels = [ "AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias", "Fragmentation", "SinglePionECAL", "SinglePionHCAL", "FlavorQCD", "TimePtEta", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeBal", "RelativeSample", "RelativeFSR", "RelativeStatFSR", "RelativeStatEC", "RelativeStatHF", "PileUpDataMC", "PileUpPtRef", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF", "PileUpMuZero", "PileUpEnvelope", "SubTotalPileUp", "SubTotalRelative", "SubTotalPt", "SubTotalScale", "SubTotalAbsolute", "SubTotalMC", "Total", "TotalNoFlavor", "TotalNoTime", "TotalNoFlavorNoTime" ]


    def beginJob(self, histFile, histDirName):
        Module.beginJob(self, histFile, histDirName)

        tauBins = 100

        ### Booking histograms
        self.addObject( ROOT.TH1F('PUweight',   ';PUWeight',   20, 0, 2) )
        #### general selection
        selList = (['_dijetSel' ] if self.selection.startswith('dijet') else [ '_WSel', '_topSel' ] )
        for isel in [ '_noSelnoWeight', '_noSel'] + selList:
            self.addObject( ROOT.TH1F('nPVs'+isel,   ';number of PVs',   100, 0, 100) )
            self.addObject( ROOT.TH1F('nleps'+isel,   ';number of leptons',   20, 0, 20) )
            self.addP4Hists( 'muons', isel )
            self.addP4Hists( 'eles', isel )
            self.addObject( ROOT.TH1F('nAK8jets'+isel,   ';number of AK8 jets',   20, 0, 20) )
            self.addP4Hists( 'AK8jets', isel )
            self.addObject( ROOT.TH1F('nAK4jets'+isel,   ';number of AK4 jets',   20, 0, 20) )
            self.addP4Hists( 'AK4jets', isel )
            self.addObject( ROOT.TH1F('METPt'+isel,   ';MET (GeV)',   200, 0, 2000) )

        for iSel in selList:
            self.addP4Hists( 'leadAK8jet', iSel )

            for sysUnc in self.sysSource:
                self.addObject( ROOT.TH1F('recoJetPt'+sysUnc+iSel, ';AK8 reco jet pt [GeV]', 500, 0, 5000) )
                self.addObject( ROOT.TH1F('recoJetSDmass'+sysUnc+iSel, ';AK8 reco jet SD mass [GeV]', 500, 0, 500) )
                self.addObject( ROOT.TH1F('recoJetEta'+sysUnc+iSel, ';AK8 reco jet eta', 40, -5, 5) )
                self.addObject( ROOT.TH1F('recoJetTau21'+sysUnc+iSel, ';AK8 reco jet #tau_{21}', tauBins, 0, 1) )
                self.addObject( ROOT.TH1F('genJetPt'+sysUnc+iSel, ';AK8 gen jet pt [GeV]', 500, 0, 5000) )
                self.addObject( ROOT.TH1F('genJetEta'+sysUnc+iSel, ';AK8 gen jet eta', 40, -5, 5) )
                self.addObject( ROOT.TH1F('genJetTau21'+sysUnc+iSel, ';AK8 gen jet #tau_{21}', tauBins, 0, 1) )
                self.addObject( ROOT.TH2F('respJetTau21'+sysUnc+iSel, ';AK8 gen jet #tau_{21};AK8 reco jet #tau_{21}', tauBins, 0, 1, tauBins, 0, 1) )

                self.addObject( ROOT.TH1F('recoSDJetPt'+sysUnc+iSel, ';AK8 reco SD jet pt [GeV]', 500, 0, 5000) )
                self.addObject( ROOT.TH1F('recoSDJetSDmass'+sysUnc+iSel, ';AK8 reco SD jet SD mass [GeV]', 500, 0, 500) )
                self.addObject( ROOT.TH1F('recoSDJetEta'+sysUnc+iSel, ';AK8 reco SD jet eta', tauBins, -5, 5) )
                self.addObject( ROOT.TH1F('recoSDJetTau21'+sysUnc+iSel, ';AK8 reco SD jet #tau_{21}', tauBins, 0, 1) )
                self.addObject( ROOT.TH2F('respSDJetTau21'+sysUnc+iSel, ';AK8 gen jet #tau_{21};AK8 reco SD jet #tau_{21}', tauBins, 0, 1, tauBins, 0, 1) )

                for tauN in range(self.maxTau):
                    for x in self.nSub_labels:

                        if sysUnc.endswith('nom'):
                            self.addObject( ROOT.TH1F("genJet"+x+str(tauN)+sysUnc+iSel, ";AK8 jet #tau", tauBins, 0, 1 ) )
                            self.addObject( ROOT.TH1F("genJet"+x+str(tauN)+"_WTA_kT"+sysUnc+iSel, ";AK8 jet #tau", tauBins, 0, 1 ) )
                            self.addObject( ROOT.TH1F("genJet"+x+str(tauN)+"_OP_kT"+sysUnc+iSel, ";AK8 jet #tau", tauBins, 0, 1 ) )
                        self.addObject( ROOT.TH1F("recoJet"+x+str(tauN)+sysUnc+iSel, ";AK8 jet #tau", tauBins, 0, 1 ) )
                        self.addObject( ROOT.TH2F('respJet'+x+str(tauN)+sysUnc+iSel, ';AK8 gen jet #tau_{21};AK8 reco jet #tau_{21}', tauBins, 0, 1, tauBins, 0, 1) )

                        self.addObject( ROOT.TH1F("recoJet"+x+str(tauN)+"_WTA_kT"+sysUnc+iSel, ";AK8 jet #tau", tauBins, 0, 1 ) )
                        self.addObject( ROOT.TH2F('respJet'+x+str(tauN)+"_WTA_kT"+sysUnc+iSel, ';AK8 gen jet #tau_{21};AK8 reco jet #tau_{21}', tauBins, 0, 1, tauBins, 0, 1) )

                        self.addObject( ROOT.TH1F("recoJet"+x+str(tauN)+"_OP_kT"+sysUnc+iSel, ";AK8 jet #tau", tauBins, 0, 1 ) )
                        self.addObject( ROOT.TH2F('respJet'+x+str(tauN)+"_OP_kT"+sysUnc+iSel, ';AK8 gen jet #tau_{21};AK8 reco jet #tau_{21}', tauBins, 0, 1, tauBins, 0, 1) )

                        self.addObject( ROOT.TH1F("recoSDJet"+x+str(tauN)+sysUnc+iSel, ";AK8 jet #tau", tauBins, 0, 1 ) )
                        self.addObject( ROOT.TH2F('respSDJet'+x+str(tauN)+sysUnc+iSel, ';AK8 gen jet #tau_{21};AK8 reco SD jet #tau_{21}', tauBins, 0, 1, tauBins, 0, 1) )

                        self.addObject( ROOT.TH1F("recoSDJet"+x+str(tauN)+"_WTA_kT"+sysUnc+iSel, ";AK8 jet #tau", tauBins, 0, 1 ) )
                        self.addObject( ROOT.TH2F('respSDJet'+x+str(tauN)+"_WTA_kT"+sysUnc+iSel, ';AK8 gen jet #tau_{21};AK8 reco SD jet #tau_{21}', tauBins, 0, 1, tauBins, 0, 1) )

                        self.addObject( ROOT.TH1F("recoSDJet"+x+str(tauN)+"_OP_kT"+sysUnc+iSel, ";AK8 jet #tau", tauBins, 0, 1 ) )
                        self.addObject( ROOT.TH2F('respSDJet'+x+str(tauN)+"_OP_kT"+sysUnc+iSel, ';AK8 gen jet #tau_{21};AK8 reco SD jet #tau_{21}', tauBins, 0, 1, tauBins, 0, 1) )

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

            for x in self.kinematic_labels:
                self.out.branch("goodgenjet" + ijet + "%s"%x,  "F")
            self.out.branch("goodgenjet" + ijet + "_tau21", "F")

            self.out.branch("dR_gen_reco_AK8", "F")
            self.out.branch("dR_genW_genAK8", "F")
            self.out.branch("genEventNo_taus_are_0",  "I")
            self.out.branch("recoEventNo_taus_are_0",  "I")

            for x in self.kinematic_labels:
                self.out.branch("sd_goodrecojet" + ijet + "%s"%x,  "F")
            self.out.branch("nsd_goodrecojet" + ijet,  "I")  ### dummy for nanoAOD Tools
            self.out.branch("sd_goodrecojet" + ijet + "_tau21",  "F")

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

        AK8jets = list(Collection(event, 'FatJet' ))
        electrons = list(Collection(event, 'Electron'))
        muons = list(Collection(event, 'Muon'))
        jets = list(Collection(event, 'Jet'))
        met = Object(event, 'MET')

        #### Lepton cleaning
        recoElectrons  = [x for x in electrons if x.pt>self.minElectronPt and x.cutBased_HEEP and ((self.range1ElectronEta[0]<abs(x.eta)<self.range1ElectronEta[1]) or (self.range2ElectronEta[0]<abs(x.eta)<self.range2ElectronEta[1]))]
        recoElectrons.sort(key=lambda x:x.pt, reverse=True)

        recoMuons = [ x for x in muons if x.pt > self.minMuonPt and x.highPtId > 1 and abs(x.p4().Eta()) < self.maxMuonEta and x.pfRelIso03_all < 0.1]
        recoMuons.sort(key=lambda x:x.pt, reverse=True)
        nleptons = len(recoMuons)+len(recoElectrons)
        ##################################################

        #### MET (not sure if needed)
        MET = ROOT.TLorentzVector()
        MET.SetPtEtaPhiE(met.pt, 0., met.phi, met.sumEt)
        ##################################################

        #### Basic AK8 jet selection
        recoAK8jets = [ x for x in AK8jets if x.pt_nom > self.minAK8JetPt and abs(x.eta) < self.maxJetAK8Eta ] # and x.jetId>4]
        recoAK8jets.sort(key=lambda x:x.pt_nom,reverse=True)
        ##################################################

        #### Basic AK4 bjet selection
        recoAK4bjets = [ x for x in jets if x.pt > self.minJetPt and abs(x.eta) < self.maxJetEta and x.btagDeepB > self.minBDisc ] # and x.jetId>4 ]
        recoAK4bjets.sort(key=lambda x:x.pt,reverse=True)
        ##################################################

        #### Weight
        if isMC: weight = event.puWeight * event.genWeight
        else: weight = 1
        ##################################################

        #### Checking no selection without weights
        getattr( self, 'nPVs_noSelnoWeight' ).Fill( getattr( event, 'PV_npvsGood') )
        getattr( self, 'nleps_noSelnoWeight' ).Fill( nleptons )
        for imuon in recoMuons:
            getattr( self, 'muons_pt_noSelnoWeight' ).Fill( imuon.pt )
            getattr( self, 'muons_eta_noSelnoWeight' ).Fill( imuon.eta )
            getattr( self, 'muons_phi_noSelnoWeight' ).Fill( imuon.phi )
        for iele in recoElectrons:
            getattr( self, 'eles_pt_noSelnoWeight' ).Fill( iele.pt )
            getattr( self, 'eles_eta_noSelnoWeight' ).Fill( iele.eta )
            getattr( self, 'eles_phi_noSelnoWeight' ).Fill( iele.phi )
        getattr( self, 'nAK8jets_noSelnoWeight' ).Fill( len(recoAK8jets) )
        for ijet in recoAK8jets:
            getattr( self, 'AK8jets_pt_noSelnoWeight' ).Fill( ijet.pt )
            getattr( self, 'AK8jets_eta_noSelnoWeight' ).Fill( ijet.eta )
            getattr( self, 'AK8jets_phi_noSelnoWeight' ).Fill( ijet.phi )
            getattr( self, 'AK8jets_mass_noSelnoWeight' ).Fill( ijet.msoftdrop )
        getattr( self, 'nAK4jets_noSelnoWeight' ).Fill( len(recoAK4bjets) )
        for ijet in recoAK4bjets:
            getattr( self, 'AK4jets_pt_noSelnoWeight' ).Fill( ijet.pt )
            getattr( self, 'AK4jets_eta_noSelnoWeight' ).Fill( ijet.eta )
            getattr( self, 'AK4jets_phi_noSelnoWeight' ).Fill( ijet.phi )
        getattr( self, 'METPt_noSelnoWeight' ).Fill( MET.Pt() )

        #### Checking no selection with weights
        getattr( self, 'nPVs_noSel' ).Fill( getattr( event, 'PV_npvsGood'), weight )
        getattr( self, 'nleps_noSel' ).Fill( nleptons, weight )
        for imuon in recoMuons:
            getattr( self, 'muons_pt_noSel' ).Fill( imuon.pt, weight )
            getattr( self, 'muons_eta_noSel' ).Fill( imuon.eta, weight )
            getattr( self, 'muons_phi_noSel' ).Fill( imuon.phi, weight )
        for iele in recoElectrons:
            getattr( self, 'eles_pt_noSel' ).Fill( iele.pt, weight )
            getattr( self, 'eles_eta_noSel' ).Fill( iele.eta, weight )
            getattr( self, 'eles_phi_noSel' ).Fill( iele.phi, weight )
        getattr( self, 'nAK8jets_noSel' ).Fill( len(recoAK8jets), weight )
        for ijet in recoAK8jets:
            getattr( self, 'AK8jets_pt_noSel' ).Fill( ijet.pt, weight )
            getattr( self, 'AK8jets_eta_noSel' ).Fill( ijet.eta, weight )
            getattr( self, 'AK8jets_phi_noSel' ).Fill( ijet.phi, weight )
            getattr( self, 'AK8jets_mass_noSel' ).Fill( ijet.msoftdrop, weight )
        getattr( self, 'nAK4jets_noSel' ).Fill( len(recoAK4bjets), weight )
        for ijet in recoAK4bjets:
            getattr( self, 'AK4jets_pt_noSel' ).Fill( ijet.pt, weight )
            getattr( self, 'AK4jets_eta_noSel' ).Fill( ijet.eta, weight )
            getattr( self, 'AK4jets_phi_noSel' ).Fill( ijet.phi, weight )
        getattr( self, 'METPt_noSel' ).Fill( MET.Pt(), weight )

        ##### Applying selection
        if self.selection.startswith('dijet'): passSelection = (nleptons==0) and (len(recoAK8jets)>1)
        elif self.selection.startswith('Wtop'):
            if (len(recoMuons)==1) and (len(recoAK8jets)>0) and (len(recoAK4bjets)>1) and (met.pt>40):

                ### removing ak4 jets inside leadAK8 jet
                for bjet in recoAK4bjets:
                    if recoAK8jets[0].p4().DeltaR( bjet.p4() )<0.8 : recoAK4bjets.remove(bjet)

                ### defining muon isolation and leptonic top
                muonIso = []
                leptonicTop = []
                for bjet in recoAK4bjets:
                    if recoMuons[0].p4().DeltaR( bjet.p4() )<0.4: muonIso.append( False )
                    if (recoMuons[0].p4().DeltaR( bjet.p4() )>0.4) and (recoMuons[0].p4().DeltaR( bjet.p4() )<1.5): leptonicTop.append( True )

                passSelection = all(muonIso) and ((MET+recoMuons[0].p4()).Pt()>200) and any(leptonicTop) and (abs(recoAK8jets[0].eta)<1.5) and (recoAK8jets[0].p4().DeltaR(recoMuons[0].p4())>0.8)
            else: passSelection = False
        else: passSelection = False

        if passSelection:
            if self.selection.startswith('dijet'):
                self.createHistosTrees( '_dijetSel', recoAK8jets[0], 450., isMC, event, weight, recoElectrons, recoMuons, MET, recoAK8jets, recoAK4bjets )
            else:
                if (recoAK8jets[0].msoftdrop<100) and (recoAK8jets[0].msoftdrop>60):
                    self.createHistosTrees( '_WSel', recoAK8jets[0], 200., isMC, event, weight, recoElectrons, recoMuons, MET, recoAK8jets, recoAK4bjets )
                elif (recoAK8jets[0].msoftdrop>140):
                    self.createHistosTrees( '_topSel', recoAK8jets[0], 350., isMC, event, weight, recoElectrons, recoMuons, MET, recoAK8jets, recoAK4bjets )

        return True

    def createHistosTrees(self, iSel, recoAK8jet, minLeadAK8JetPt, isMC, event, weight, recoElectrons, recoMuons, MET, recoAK8jets, recoAK4bjets):

        pfCands = list(Collection(event, 'PFCandsAK8' ))
        passAK8jet = {}          ### Storing goodreco jet as list for later use
        passAK8jet['jet'] = recoAK8jet

        ##### Computing quantities
        if (recoAK8jet.pt > minLeadAK8JetPt*0.8):    #### store jets in the range of the pt selection

            #### Run calculations of NSub bases and store for ungroomed recoAK8jets (default in CMS)

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
                if abs(recoAK8jet.p4().DeltaR( x )) < 0.8: constituents.push_back(x)

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

            try: passAK8jet['tau21'] = recoAK8jet.tau2/recoAK8jet.tau1
            except ZeroDivisionError: passAK8jet['tau21'] = -1
            try: passAK8jet['tau32'] = recoAK8jet.tau3/recoAK8jet.tau2
            except ZeroDivisionError: passAK8jet['tau32'] = -1

            #### Filling branch with passAK8jet info after selection
            self.out.fillBranch("ngoodrecojet0", 1 )  ### dummy for nanoAOD Tools
            self.out.fillBranch("goodrecojet0_pt",  recoAK8jet.pt_nom )
            self.out.fillBranch("goodrecojet0_eta",  recoAK8jet.p4().Eta() )
            self.out.fillBranch("goodrecojet0_phi",  recoAK8jet.p4().Phi() )
            self.out.fillBranch("goodrecojet0_mass",  recoAK8jet.p4().M() )
            self.out.fillBranch("goodrecojet0_softdrop_mass", recoAK8jet.msoftdrop_nom)
            self.out.fillBranch("goodrecojet0_pt_raw",  recoAK8jet.pt_raw )
            self.out.fillBranch("goodrecojet0_mass_raw",  recoAK8jet.mass_raw )
            self.out.fillBranch("goodrecojet0_tau21", passAK8jet['tau21'] )
            self.out.fillBranch("goodrecojet0_N21",  recoAK8jet.n2b1)
            self.out.fillBranch("goodrecojet0_corr_JEC",  recoAK8jet.corr_JEC )
            self.out.fillBranch("goodrecojet0_corr_JER",  ( recoAK8jet.corr_JER if isMC else -99999 ) )
            self.out.fillBranch("goodrecojet0_corr_JMS",  ( recoAK8jet.corr_JMS if isMC else -99999 ) )
            self.out.fillBranch("goodrecojet0_corr_JMR",  ( recoAK8jet.corr_JMR if isMC else -99999 ) )
            self.out.fillBranch("goodrecojet0_mass_jmrUp",  ( recoAK8jet.mass_jmrUp if isMC else -99999 ) )
            self.out.fillBranch("goodrecojet0_mass_jmsUp",  ( recoAK8jet.mass_jmsUp if isMC else -99999 ) )
            self.out.fillBranch("goodrecojet0_mass_jmrDown",  ( recoAK8jet.mass_jmrDown if isMC else -99999 ) )
            self.out.fillBranch("goodrecojet0_mass_jmsDown",  ( recoAK8jet.mass_jmsDown if isMC else -99999 ) )
            for q in ['pt', 'mass']:
                for u in [ 'Up', 'Down' ]:
                    self.out.fillBranch("goodrecojet0_"+q+'_jer'+u, ( getattr( recoAK8jet, q+'_jer'+u ) if isMC else -99999 ) )
                    for jes in self.JESLabels:
                        self.out.fillBranch("goodrecojet0_"+q+'_jesUnc'+jes+u, ( getattr( recoAK8jet, q+'_jes'+jes+u ) if isMC else -99999 ) )


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
            ################################################## end of ungroomed jets

            ##################################################
            #### Run calculations of NSub bases and store for groomed recoAK8jets

            #### Computing Softdrop jets
            sdrecoAK8jets = self.sd.result( recoCandsPUPPIweightedVec )

            ### Storing goodreco jet as list for later use
            if len(sdrecoAK8jets)>0:

                deltaR = 99999
                sdIndex = -1
                for isdInd, isdjet in enumerate(sdrecoAK8jets):

                    #### Checking which recojet is matched to the sdrecoAK8jets
                    tmpisdjet = ROOT.TLorentzVector( )
                    tmpisdjet.SetPtEtaPhiM( isdjet.perp(), isdjet.eta(), isdjet.phi(), isdjet.m() )
                    tmpDeltaR = tmpisdjet.DeltaR( passAK8jet['jet'].p4() )
                    if (tmpDeltaR < deltaR):
                        deltaR=tmpDeltaR
                        sdIndex= isdInd

                passAK8jet['sdjet'] = sdrecoAK8jets[sdIndex]

                # Cluster only the particles near the appropriate jet to save time
                sd_constituents = ROOT.vector("TLorentzVector")()

                for x in passAK8jet['sdjet'].constituents():
                    sd_constits = ROOT.TLorentzVector( x.px(), x.py(), x.pz(), x.E())
                    if abs(passAK8jet['sdjet'].delta_R( x )) < 0.8:
                        sd_constituents.push_back(sd_constits)
                sd_nsub0p5 = self.nSub0p5.getTau( self.maxTau, sd_constituents )
                sd_nsub1 = self.nSub1.getTau( self.maxTau, sd_constituents )
                sd_nsub2 = self.nSub2.getTau( self.maxTau, sd_constituents )

                sd_nsub0p5_WTA_kT = self.nSub0p5_WTA_kT.getTau( self.maxTau, sd_constituents )
                sd_nsub1_WTA_kT = self.nSub1_WTA_kT.getTau( self.maxTau, sd_constituents )
                sd_nsub2_WTA_kT = self.nSub2_WTA_kT.getTau( self.maxTau, sd_constituents )

                sd_nsub0p5_OP_kT = self.nSub0p5_OP_kT.getTau( self.maxTau, sd_constituents )
                sd_nsub1_OP_kT = self.nSub1_OP_kT.getTau( self.maxTau, sd_constituents )
                sd_nsub2_OP_kT = self.nSub2_OP_kT.getTau( self.maxTau, sd_constituents )

                self.out.fillBranch("nsd_goodrecojet0", 1 )  ## dummy
                self.out.fillBranch("sd_goodrecojet0_pt",  passAK8jet['sdjet'].perp() )
                self.out.fillBranch("sd_goodrecojet0_eta",  passAK8jet['sdjet'].eta() )
                self.out.fillBranch("sd_goodrecojet0_phi",  passAK8jet['sdjet'].phi() )
                self.out.fillBranch("sd_goodrecojet0_mass",  passAK8jet['sdjet'].m() )

                try: passAK8jet['sdtau21'] = sd_nsub1[1]/sd_nsub1[0]
                except ZeroDivisionError: passAK8jet['sdtau21'] = -1
                try: passAK8jet['sdtau32'] = sd_nsub1[2]/sd_nsub1[1]
                except ZeroDivisionError: passAK8jet['sdtau32'] = -1
                self.out.fillBranch("sd_goodrecojet0_tau21", passAK8jet['sdtau21'] )

                for tauN in range(self.maxTau):
                    self.out.fillBranch("sd_goodrecojet0_tau_0p5_"+str(tauN),  sd_nsub0p5[tauN]  )
                    self.out.fillBranch("sd_goodrecojet0_tau_1_"+str(tauN),  sd_nsub1[tauN]  )
                    self.out.fillBranch("sd_goodrecojet0_tau_2_"+str(tauN),  sd_nsub2[tauN]  )
                    passAK8jet['sd0p5'+str(tauN)] = sd_nsub0p5[tauN]
                    passAK8jet['sd1'+str(tauN)] = sd_nsub1[tauN]
                    passAK8jet['sd2'+str(tauN)] = sd_nsub2[tauN]

                    self.out.fillBranch("sd_goodrecojet0_tau_0p5_"+str(tauN) + "_WTA_kT",  sd_nsub0p5_WTA_kT[tauN]  )
                    self.out.fillBranch("sd_goodrecojet0_tau_1_"+str(tauN) + "_WTA_kT",  sd_nsub1_WTA_kT[tauN]  )
                    self.out.fillBranch("sd_goodrecojet0_tau_2_"+str(tauN) + "_WTA_kT",  sd_nsub2_WTA_kT[tauN]  )
                    passAK8jet['sd0p5WTAkT'+str(tauN)] = sd_nsub0p5_WTA_kT[tauN]
                    passAK8jet['sd1WTAkT'+str(tauN)] = sd_nsub1_WTA_kT[tauN]
                    passAK8jet['sd2WTAkT'+str(tauN)] = sd_nsub2_WTA_kT[tauN]

                    self.out.fillBranch("sd_goodrecojet0_tau_0p5_"+str(tauN) + "_OP_kT",  sd_nsub0p5_OP_kT[tauN]  )
                    self.out.fillBranch("sd_goodrecojet0_tau_1_"+str(tauN) + "_OP_kT",  sd_nsub1_OP_kT[tauN]  )
                    self.out.fillBranch("sd_goodrecojet0_tau_2_"+str(tauN) + "_OP_kT",  sd_nsub2_OP_kT[tauN]  )
                    passAK8jet['sd0p5OPkT'+str(tauN)] = sd_nsub0p5_OP_kT[tauN]
                    passAK8jet['sd1OPkT'+str(tauN)] = sd_nsub1_OP_kT[tauN]
                    passAK8jet['sd2OPkT'+str(tauN)] = sd_nsub2_OP_kT[tauN]
            ##########################################################

            #### For MC only, loop over the genJet collection
            if isMC:
                genParticles = Collection(event, 'GenPartAK8')
                genjets  = list(Collection(event, 'GenJetAK8' ))

                #### Matching recojet and genJet
                tmpGenInd = -1
                tmpGenDeltaR = 99999
                for i in xrange(0,len(genjets)):
                    genDeltaR = abs(genjets[i].p4().DeltaR(recoAK8jet.p4()))
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
                except ZeroDivisionError: passAK8jet['gentau21'] = -1
                try: passAK8jet['gentau32'] = genNsub1[2]/genNsub1[1]
                except ZeroDivisionError: passAK8jet['gentau32'] = -1
                self.out.fillBranch("goodgenjet0_tau21", passAK8jet['gentau21'] )

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
            ################################################## end of gen jets

        ##################### end of computing quantities

        ##### Filling histograms
        for sysUnc in self.sysSource:
            if( getattr(passAK8jet['jet'], 'pt'+(sysUnc if not sysUnc.startswith('_pu') else '_nom')) > minLeadAK8JetPt ):

                #### Checking nominal selection with weights
                if sysUnc.endswith('nom'):
                    getattr( self, 'nPVs'+iSel ).Fill( getattr( event, 'PV_npvsGood'), weight )
                    getattr( self, 'nleps'+iSel ).Fill( len(recoMuons)+len(recoElectrons), weight )
                    for imuon in recoMuons:
                        getattr( self, 'muons_pt'+iSel ).Fill( imuon.pt, weight )
                        getattr( self, 'muons_eta'+iSel ).Fill( imuon.eta, weight )
                        getattr( self, 'muons_phi'+iSel ).Fill( imuon.phi, weight )
                    for iele in recoElectrons:
                        getattr( self, 'eles_pt'+iSel ).Fill( iele.pt, weight )
                        getattr( self, 'eles_eta'+iSel ).Fill( iele.eta, weight )
                        getattr( self, 'eles_phi'+iSel ).Fill( iele.phi, weight )
                    getattr( self, 'nAK8jets'+iSel ).Fill( len(recoAK8jets), weight )
                    for ijet in recoAK8jets:
                        getattr( self, 'AK8jets_pt'+iSel ).Fill( ijet.pt, weight )
                        getattr( self, 'AK8jets_eta'+iSel ).Fill( ijet.eta, weight )
                        getattr( self, 'AK8jets_phi'+iSel ).Fill( ijet.phi, weight )
                        getattr( self, 'AK8jets_mass'+iSel ).Fill( ijet.msoftdrop, weight )
                    getattr( self, 'nAK4jets'+iSel ).Fill( len(recoAK4bjets), weight )
                    for ijet in recoAK4bjets:
                        getattr( self, 'AK4jets_pt'+iSel ).Fill( ijet.pt, weight )
                        getattr( self, 'AK4jets_eta'+iSel ).Fill( ijet.eta, weight )
                        getattr( self, 'AK4jets_phi'+iSel ).Fill( ijet.phi, weight )
                    getattr( self, 'METPt'+iSel ).Fill( getattr(event,'MET_pt'), weight )
                    getattr( self, 'leadAK8jet_pt'+iSel ).Fill( getattr(passAK8jet['jet'], 'pt'), weight )
                    getattr( self, 'leadAK8jet_eta'+iSel ).Fill( getattr(passAK8jet['jet'], 'eta'), weight )
                    getattr( self, 'leadAK8jet_phi'+iSel ).Fill( getattr(passAK8jet['jet'], 'phi'), weight )
                    getattr( self, 'leadAK8jet_mass'+iSel ).Fill( getattr(passAK8jet['jet'], 'msoftdrop'), weight )
                    if isMC:
                        getattr( self, 'genJetPt'+sysUnc+iSel ).Fill( passAK8jet['genjet'].pt, weight )
                        getattr( self, 'genJetEta'+sysUnc+iSel ).Fill( passAK8jet['genjet'].eta, weight )
                        getattr( self, 'genJetTau21'+sysUnc+iSel ).Fill( passAK8jet['gentau21'], weight )
                        for tauN in range(self.maxTau):
                            getattr( self, "genJet_tau_0p5_"+str(tauN)+sysUnc+iSel ).Fill( passAK8jet['gen0p5'+str(tauN)], weight )
                            getattr( self, "genJet_tau_1_"+str(tauN)+sysUnc+iSel ).Fill( passAK8jet['gen1'+str(tauN)], weight )
                            getattr( self, "genJet_tau_2_"+str(tauN)+sysUnc+iSel ).Fill( passAK8jet['gen2'+str(tauN)], weight )
                            getattr( self, "genJet_tau_0p5_"+str(tauN)+'_WTA_kT'+sysUnc+iSel ).Fill( passAK8jet['gen0p5WTAkT'+str(tauN)], weight )
                            getattr( self, "genJet_tau_1_"+str(tauN)+'_WTA_kT'+sysUnc+iSel ).Fill( passAK8jet['gen1WTAkT'+str(tauN)], weight )
                            getattr( self, "genJet_tau_2_"+str(tauN)+'_WTA_kT'+sysUnc+iSel ).Fill( passAK8jet['gen2WTAkT'+str(tauN)], weight )
                            getattr( self, "genJet_tau_0p5_"+str(tauN)+'_OP_kT'+sysUnc+iSel ).Fill( passAK8jet['gen0p5OPkT'+str(tauN)], weight )
                            getattr( self, "genJet_tau_1_"+str(tauN)+'_OP_kT'+sysUnc+iSel ).Fill( passAK8jet['gen1OPkT'+str(tauN)], weight )
                            getattr( self, "genJet_tau_2_"+str(tauN)+'_OP_kT'+sysUnc+iSel ).Fill( passAK8jet['gen2OPkT'+str(tauN)], weight )

                #### filling reco histos
                if sysUnc.startswith('_pu'):
                    WEIGHT = event.genWeight + getattr( event, 'puWeight'+sysUnc.split('pu')[1] )
                    getattr( self, 'recoJetPt'+sysUnc+iSel ).Fill( getattr(passAK8jet['jet'], 'pt_nom' ), WEIGHT )
                    getattr( self, 'recoJetSDmass'+sysUnc+iSel ).Fill( getattr(passAK8jet['jet'], 'msoftdrop_nom' ), WEIGHT )
                else:
                    WEIGHT =  weight
                    getattr( self, 'recoJetPt'+sysUnc+iSel ).Fill( getattr(passAK8jet['jet'], 'pt'+sysUnc ), WEIGHT )
                    getattr( self, 'recoJetSDmass'+sysUnc+iSel ).Fill( getattr(passAK8jet['jet'], 'msoftdrop'+sysUnc ), WEIGHT )
                getattr( self, 'recoJetEta'+sysUnc+iSel ).Fill( getattr(passAK8jet['jet'], 'eta'), WEIGHT )
                getattr( self, 'recoJetTau21'+sysUnc+iSel ).Fill( passAK8jet['tau21'], WEIGHT )
                getattr( self, 'respJetTau21'+sysUnc+iSel ).Fill( passAK8jet['gentau21'], passAK8jet['tau21'], WEIGHT )
                for tauN in range(self.maxTau):
                    getattr( self, "recoJet_tau_0p5_"+str(tauN)+sysUnc+iSel ).Fill( passAK8jet['0p5'+str(tauN)], WEIGHT )
                    getattr( self, 'respJet_tau_0p5_'+str(tauN)+sysUnc+iSel ).Fill( passAK8jet['gen0p5'+str(tauN)], passAK8jet['0p5'+str(tauN)], WEIGHT )
                    getattr( self, "recoJet_tau_1_"+str(tauN)+sysUnc+iSel ).Fill( passAK8jet['1'+str(tauN)], WEIGHT )
                    getattr( self, 'respJet_tau_1_'+str(tauN)+sysUnc+iSel ).Fill( passAK8jet['gen1'+str(tauN)], passAK8jet['1'+str(tauN)], WEIGHT )
                    getattr( self, "recoJet_tau_2_"+str(tauN)+sysUnc+iSel ).Fill( passAK8jet['2'+str(tauN)], WEIGHT )
                    getattr( self, 'respJet_tau_2_'+str(tauN)+sysUnc+iSel ).Fill( passAK8jet['gen2'+str(tauN)], passAK8jet['2'+str(tauN)], WEIGHT )
                    getattr( self, "recoJet_tau_0p5_"+str(tauN)+'_WTA_kT'+sysUnc+iSel ).Fill( passAK8jet['0p5WTAkT'+str(tauN)], WEIGHT )
                    getattr( self, 'respJet_tau_0p5_'+str(tauN)+'_WTA_kT'+sysUnc+iSel ).Fill( passAK8jet['gen0p5WTAkT'+str(tauN)], passAK8jet['0p5WTAkT'+str(tauN)], WEIGHT )
                    getattr( self, "recoJet_tau_1_"+str(tauN)+'_WTA_kT'+sysUnc+iSel ).Fill( passAK8jet['1WTAkT'+str(tauN)], WEIGHT )
                    getattr( self, 'respJet_tau_1_'+str(tauN)+'_WTA_kT'+sysUnc+iSel ).Fill( passAK8jet['gen1WTAkT'+str(tauN)], passAK8jet['1WTAkT'+str(tauN)], WEIGHT )
                    getattr( self, "recoJet_tau_2_"+str(tauN)+'_WTA_kT'+sysUnc+iSel ).Fill( passAK8jet['2WTAkT'+str(tauN)], WEIGHT )
                    getattr( self, 'respJet_tau_2_'+str(tauN)+'_WTA_kT'+sysUnc+iSel ).Fill( passAK8jet['gen2WTAkT'+str(tauN)], passAK8jet['2WTAkT'+str(tauN)], WEIGHT )
                    getattr( self, "recoJet_tau_0p5_"+str(tauN)+'_OP_kT'+sysUnc+iSel ).Fill( passAK8jet['0p5OPkT'+str(tauN)], WEIGHT )
                    getattr( self, 'respJet_tau_0p5_'+str(tauN)+'_OP_kT'+sysUnc+iSel ).Fill( passAK8jet['gen0p5OPkT'+str(tauN)], passAK8jet['0p5OPkT'+str(tauN)], WEIGHT )
                    getattr( self, "recoJet_tau_1_"+str(tauN)+'_OP_kT'+sysUnc+iSel ).Fill( passAK8jet['1OPkT'+str(tauN)], WEIGHT )
                    getattr( self, 'respJet_tau_1_'+str(tauN)+'_OP_kT'+sysUnc+iSel ).Fill( passAK8jet['gen1OPkT'+str(tauN)], passAK8jet['1OPkT'+str(tauN)], WEIGHT )
                    getattr( self, "recoJet_tau_2_"+str(tauN)+'_OP_kT'+sysUnc+iSel ).Fill( passAK8jet['2OPkT'+str(tauN)], WEIGHT )
                    getattr( self, 'respJet_tau_2_'+str(tauN)+'_OP_kT'+sysUnc+iSel ).Fill( passAK8jet['gen2OPkT'+str(tauN)], passAK8jet['2OPkT'+str(tauN)], WEIGHT )


                if passAK8jet['sdjet']:

                    if sysUnc.startswith('_nom'):
                        getattr( self, 'recoSDJetPt'+sysUnc+iSel ).Fill( passAK8jet['sdjet'].perp(), WEIGHT )
                        getattr( self, 'recoSDJetSDmass'+sysUnc+iSel ).Fill( passAK8jet['sdjet'].m(), WEIGHT )
                        getattr( self, 'recoSDJetEta'+sysUnc+iSel ).Fill( passAK8jet['sdjet'].eta(), WEIGHT )
                        getattr( self, 'recoSDJetTau21'+sysUnc+iSel ).Fill( passAK8jet['sdtau21'], WEIGHT )
                        getattr( self, 'respSDJetTau21'+sysUnc+iSel ).Fill( passAK8jet['gentau21'], passAK8jet['sdtau21'], WEIGHT )

                    for tauN in range(self.maxTau):
                        getattr( self, "recoSDJet_tau_0p5_"+str(tauN)+sysUnc+iSel ).Fill( passAK8jet['sd0p5'+str(tauN)], WEIGHT )
                        getattr( self, 'respSDJet_tau_0p5_'+str(tauN)+sysUnc+iSel ).Fill( passAK8jet['gen0p5'+str(tauN)], passAK8jet['sd0p5'+str(tauN)], WEIGHT )
                        getattr( self, "recoSDJet_tau_1_"+str(tauN)+sysUnc+iSel ).Fill( passAK8jet['sd1'+str(tauN)], WEIGHT )
                        getattr( self, 'respSDJet_tau_1_'+str(tauN)+sysUnc+iSel ).Fill( passAK8jet['gen1'+str(tauN)], passAK8jet['sd1'+str(tauN)], WEIGHT )
                        getattr( self, "recoSDJet_tau_2_"+str(tauN)+sysUnc+iSel ).Fill( passAK8jet['sd2'+str(tauN)], WEIGHT )
                        getattr( self, 'respSDJet_tau_2_'+str(tauN)+sysUnc+iSel ).Fill( passAK8jet['gen2'+str(tauN)], passAK8jet['sd2'+str(tauN)], WEIGHT )
                        getattr( self, "recoSDJet_tau_0p5_"+str(tauN)+'_WTA_kT'+sysUnc+iSel ).Fill( passAK8jet['sd0p5WTAkT'+str(tauN)], WEIGHT )
                        getattr( self, 'respSDJet_tau_0p5_'+str(tauN)+'_WTA_kT'+sysUnc+iSel ).Fill( passAK8jet['gen0p5WTAkT'+str(tauN)], passAK8jet['sd0p5WTAkT'+str(tauN)], WEIGHT )
                        getattr( self, "recoSDJet_tau_1_"+str(tauN)+'_WTA_kT'+sysUnc+iSel ).Fill( passAK8jet['sd1WTAkT'+str(tauN)], WEIGHT )
                        getattr( self, 'respSDJet_tau_1_'+str(tauN)+'_WTA_kT'+sysUnc+iSel ).Fill( passAK8jet['gen1WTAkT'+str(tauN)], passAK8jet['sd1WTAkT'+str(tauN)], WEIGHT )
                        getattr( self, "recoSDJet_tau_2_"+str(tauN)+'_WTA_kT'+sysUnc+iSel ).Fill( passAK8jet['sd2WTAkT'+str(tauN)], WEIGHT )
                        getattr( self, 'respSDJet_tau_2_'+str(tauN)+'_WTA_kT'+sysUnc+iSel ).Fill( passAK8jet['gen2WTAkT'+str(tauN)], passAK8jet['sd2WTAkT'+str(tauN)], WEIGHT )
                        getattr( self, "recoSDJet_tau_0p5_"+str(tauN)+'_OP_kT'+sysUnc+iSel ).Fill( passAK8jet['sd0p5OPkT'+str(tauN)], WEIGHT )
                        getattr( self, 'respSDJet_tau_0p5_'+str(tauN)+'_OP_kT'+sysUnc+iSel ).Fill( passAK8jet['gen0p5OPkT'+str(tauN)], passAK8jet['sd0p5OPkT'+str(tauN)], WEIGHT )
                        getattr( self, "recoSDJet_tau_1_"+str(tauN)+'_OP_kT'+sysUnc+iSel ).Fill( passAK8jet['sd1OPkT'+str(tauN)], WEIGHT )
                        getattr( self, 'respSDJet_tau_1_'+str(tauN)+'_OP_kT'+sysUnc+iSel ).Fill( passAK8jet['gen1OPkT'+str(tauN)], passAK8jet['sd1OPkT'+str(tauN)], WEIGHT )
                        getattr( self, "recoSDJet_tau_2_"+str(tauN)+'_OP_kT'+sysUnc+iSel ).Fill( passAK8jet['sd2OPkT'+str(tauN)], WEIGHT )
                        getattr( self, 'respSDJet_tau_2_'+str(tauN)+'_OP_kT'+sysUnc+iSel ).Fill( passAK8jet['gen2OPkT'+str(tauN)], passAK8jet['sd2OPkT'+str(tauN)], WEIGHT )


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

