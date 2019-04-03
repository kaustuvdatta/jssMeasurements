import FWCore.ParameterSet.Config as cms
process = cms.Process('NANO')
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(),
)
process.source.fileNames = [

    #  2017
    # data
    #'/store/user/rappocc/JetHT/JetHT_Run2017B-17Nov2017-v1/180525_135243/0000/test_data_94X_NANO_90.root'
    #'/store/user/rappocc/JetHT/JetHT_Run2017D-17Nov2017-v1/180330_195538/0000/test_data_94X_NANO_85.root'
    # mc:
    #'/store/user/rappocc/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/180330_193826/0000/test94X_NANO_69.root'
    #'root://cms-xrd-global.cern.ch//store/user/asparker/NanoAODJMARTools-skims/nanoskim-JetsandLepton-94XMC-TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8-1of3-trees.root'
    #'root://cmseos.fnal.gov//store/user/srappocc/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8_RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2/180330_194407/0000/test94X_NANO_397.root'
    #'file:nano102x_on_mini94x_2016_mc_NANO_98.root',
    'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/BulkGravToWW_narrow_M-3000_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/A6F736EF-1515-E711-8C8A-48FD8E2824DB.root'
    # 2016
    # data
    #'/store/user/srappocc/JetHT/JetHT_Run2016H-07Aug17-v1/180329_210546/0000/test_data_80X_NANO_101.root'
]
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

#process.options = cms.untracked.PSet()

process.output = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('jetObservables_nanoskim.root'),
        fakeNameForCrab =cms.untracked.bool(True))
process.out = cms.EndPath(process.output)
