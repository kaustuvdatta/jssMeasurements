import FWCore.ParameterSet.Config as cms
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJets
from RecoJets.Configuration.RecoPFJets_cff import ak4PFJets

process = cms.Process("Demo")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1)  )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring([
'/store/mc/RunIISummer16MiniAODv2/WprimeToWZToWhadZinv_narrow_M-1000_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/00000/FE6534F4-7D02-E711-BC7C-02163E014328.root', '/store/mc/RunIISummer16MiniAODv2/WprimeToWZToWhadZinv_narrow_M-1000_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/00000/FE5ACB8A-E801-E711-AEAB-02163E01A443.root', '/store/mc/RunIISummer16MiniAODv2/WprimeToWZToWhadZinv_narrow_M-1200_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/00000/FEB8CE02-EB02-E711-BFF7-E41D2D08DFB0.root', '/store/mc/RunIISummer16MiniAODv2/WprimeToWZToWhadZinv_narrow_M-1200_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/00000/D0E4D11B-EA02-E711-A679-848F69FD288F.root', '/store/mc/RunIISummer16MiniAODv2/WprimeToWZToWhadZinv_narrow_M-1200_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/00000/CC323E51-EA02-E711-9852-20CF307C9944.root', '/store/mc/RunIISummer16MiniAODv2/WprimeToWZToWhadZinv_narrow_M-1200_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/00000/A88E7438-EA02-E711-88D2-00266CFFCB7C.root', '/store/mc/RunIISummer16MiniAODv2/WprimeToWZToWhadZinv_narrow_M-1200_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/00000/A485B54A-EA02-E711-A8F4-00259021A526.root', '/store/mc/RunIISummer16MiniAODv2/WprimeToWZToWhadZinv_narrow_M-1200_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/00000/96129FD7-7202-E711-90A8-002590D0AFB4.root', '/store/mc/RunIISummer16MiniAODv2/WprimeToWZToWhadZinv_narrow_M-1200_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/00000/88F62F85-EA02-E711-8169-FA163E91EB47.root', '/store/mc/RunIISummer16MiniAODv2/WprimeToWZToWhadZinv_narrow_M-1200_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/00000/60EE974D-EA02-E711-A185-002590E7DFA2.root', '/store/mc/RunIISummer16MiniAODv2/WprimeToWZToWhadZinv_narrow_M-1200_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/00000/3C4C901B-EA02-E711-BE44-A4BF0107E164.root', '/store/mc/RunIISummer16MiniAODv2/WprimeToWZToWhadZinv_narrow_M-1200_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/00000/00A6ABE4-8802-E711-B6E9-02163E011B62.root', '/store/mc/RunIISummer16MiniAODv2/WprimeToWZToWhadZinv_narrow_M-1600_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/F85EB5C6-C802-E711-9CB0-047D7B881D74.root', '/store/mc/RunIISummer16MiniAODv2/WprimeToWZToWhadZinv_narrow_M-1600_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/8828EBFA-C802-E711-AF5E-001C23C0A63C.root', '/store/mc/RunIISummer16MiniAODv2/WprimeToWZToWhadZinv_narrow_M-3500_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/D8DDDA4C-FA02-E711-AD8C-A4BF0108B4E2.root', '/store/mc/RunIISummer16MiniAODv2/WprimeToWZToWhadZinv_narrow_M-3500_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/AE315850-FA02-E711-9385-FA163EE355DD.root', '/store/mc/RunIISummer16MiniAODv2/WprimeToWZToWhadZinv_narrow_M-3500_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/90000/9E0A736B-4603-E711-9700-FA163EE55D2F.root',
            ])
        )

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

####### CMS way to save output root files
process.TFileService = cms.Service("TFileService",
        fileName = cms.string("skimOutputFile_Wprime_M-1000_1200_1600_3500.root"),
        closeFileFast = cms.untracked.bool(True)
        )


####### Selecting genparticles without muons
process.GenParticlesNoNu = cms.EDFilter("CandPtrSelector",
        src = cms.InputTag("packedGenParticles"),
        cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16")
        )

###### CLusteting genJets, and applying pt cut
process.ak8GenJets = ak4GenJets.clone(
#process.ak8GenJets = ak4PFJets.clone(
        #src = "packedPFCandidates",
        src = "GenParticlesNoNu",
        rParam = 0.8,
        jetAlgorithm = 'AntiKt',
        jetPtMin = 200,
        )


###### Calculating Nsubjetiness for GenJet collection
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
process.newNsub0p5 = Njettiness.clone(
        src = cms.InputTag( 'ak8GenJets' ),
        Njets=cms.vuint32(range(1,8)),         # compute 1-, 2-, 3-, 4- subjettiness
        # variables for measure definition :
        #measureDefinition = cms.uint32( 0 ), # CMS default is normalized measure
        beta = cms.double(0.5),              # CMS default is 1
        R0 = cms.double( 0.8 ),              # CMS default is jet cone size
        axesDefinition = cms.uint32( 0 ),    # CMS default is 1-pass KT axes (6), KT_Axes (0)
        )
process.newNsub1 = process.newNsub0p5.clone( beta = 1.0 )
process.newNsub2 = process.newNsub0p5.clone( beta = 2.0 )

###### Running analyzer, which creates a tree with variables
process.addingNsub = cms.EDAnalyzer("Skimmer",
        GenjetNsub0p5 = cms.string( 'newNsub0p5' ) ,
        GenjetNsub1 = cms.string( 'newNsub1' ) ,
        GenjetNsub2 = cms.string( 'newNsub2' ) ,
        srcGenJet = cms.InputTag( 'ak8GenJets'),
        )


process.path = cms.Path(
        process.GenParticlesNoNu
        * process.ak8GenJets
        * process.newNsub0p5
        * process.newNsub1
        * process.newNsub2
        * process.addingNsub
         )
