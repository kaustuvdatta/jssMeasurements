#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor


from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

#this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles,runsAndLumis

from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeight_2016, puAutoWeight_2016 #puWeight_2016
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import btagSF2016
#from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import jecUncertAll #jetRecalib2016BCD, jetRecalib2016BCDAK8Puppi

from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import jetmetUncertainties2016AK8PuppiAll

# our module
from Analysis.jetObservables.nSubProducer import nsubjettinessProducer
#from Analysis.jetObservables.nSubproducer_WP import nsubjettinessProducer
#from Analysis.jetObservables.nSubProducer_matched_TTSemilept import nsubjettinessProducer
#from Analysis.jetObservables.nSubProducer_matched_bbveto_TTSemilept import nsubjettinessProducer
# Output file name, and decide if this is MC or data


haddname = "jetObservables_nanoskim.root"


#if any( ['JetHT' in i for i in inputFiles() ] ):
#    isData = True
#else :
#    isData = False

print '---------------------------------------------------'
print 'Input files:'
print inputFiles()


modulesToRun = [ puWeight_2016(),
        jetmetUncertainties2016AK8PuppiAll(),
        btagSF2016(),
        nsubjettinessProducer()
        ]
p1=PostProcessor( ".", inputFiles(), "", "keep_and_drop.txt",
                    modules=modulesToRun,
                    provenance=True,
                    fwkJobReport=True,
                    #jsonInput=runsAndLumis(),
                    #prefetch=True,
                    #maxEntries=20000,
                    haddFileName=haddname )

p1.run()
print "DONE"
os.system("ls -lR")

