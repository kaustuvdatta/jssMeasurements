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

from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeight_2016, puAutoWeight_2016
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import btagSF2016
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *

# our module
from Analysis.jetObservables.nSubProducer_dijet import nsubjettinessProducerDijet
from Analysis.jetObservables.nSubProducer_boostedW_axdefScan import nsubjettinessProducer

import argparse

parser = argparse.ArgumentParser(description='Runs MEAnalysis')
parser.add_argument(
    '--sample',
    action="store",
    help="Sample to process",
    default='ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8'
)
parser.add_argument(
    '--numEvents',
    action="store",
    type=int,
    help="Number of events to process",
    default=1000000000000,
)
parser.add_argument(
    '--iFile',
    action="store",
    help="Input file (for condor)",
    default=""
)
parser.add_argument(
    '--oFile',
    action="store",
    help="Output file (for condor)",
    default=""
)
parser.add_argument(
    '--local',
    action="store_true",
    help="Run local or condor/crab"
)
parser.add_argument(
    '--year',
    action="store",
    help="year of data",
    choices=["2016", "2017", "2018"],
    default="2016",
    required=False
)
parser.add_argument(
    '--selection',
    action="store",
    help="Event selection",
    choices=["dijet", "W", "top"],
    default="W",
    required=False
)
args = parser.parse_args(sys.argv[1:])
if args.sample.startswith(('/EGamma', '/Single', '/JetHT' )) or ('EGamma' in args.iFile or 'Single' in args.iFile or ('JetHT' in args.iFile)):
    isMC = False
    print "sample is data"
else: isMC = True

### General selections:
PV = "(PV_npvsGood>0)"
METFilters = "( (Flag_goodVertices==1) && (Flag_globalSuperTightHalo2016Filter==1) && (Flag_HBHENoiseFilter==1) && (Flag_HBHENoiseIsoFilter==1) && (Flag_EcalDeadCellTriggerPrimitiveFilter==1) && (Flag_BadPFMuonFilter==1) )"
if not isMC: METFilters = METFilters + ' && (Flag_eeBadScFilter==1)'

if args.selection.startswith('dijet'):
    Triggers = '(HLT_PFHT900)'  #### need to include other triggers
else:
    Triggers = '(HLT_Mu50==1)'
#if args.year.startswith('2016'): Triggers = ...
#elif args.year.startswith('2017'): Triggers =  ...
#elif args.year.startswith('2018'): Triggers = ...

cuts = PV + " && " + METFilters + " && " + Triggers

#### Modules to run
jetmetCorrector = createJMECorrector(isMC=isMC, dataYear=args.year, jesUncert="All", redojec=True)
fatJetCorrector = createJMECorrector(isMC=isMC, dataYear=args.year, jesUncert="All", redojec=True, jetType = "AK8PFPuppi")

modulesToRun = []
if isMC: modulesToRun.append( puWeight_2016() )
modulesToRun.append( jetmetCorrector() )
modulesToRun.append( fatJetCorrector() )
#if isMC: modulesToRun.append( btagSF2016() )
if args.selection.startswith('dijet'): modulesToRun.append( nsubjettinessProducerDijet() )
else: modulesToRun.append( nsubjettinessProducer() )


#### Make it run
p1=PostProcessor(
        '.', (inputFiles() if not args.iFile else [args.iFile]),
        cut=cuts,
        #branchsel="keep_and_drop.txt",
        modules=modulesToRun,
        provenance=True,
        fwkJobReport=True,
        #jsonInput=runsAndLumis(),
        maxEntries=args.numEvents,
        prefetch=args.local,
        longTermCache=args.local,
        haddFileName= "jetObservables_"+args.selection+"_nanoskim.root",
        histFileName = "jetObservables_"+args.selection+"_histograms.root",
        histDirName = 'jetObservables',
        )
p1.run()
print "DONE"
if not args.local: os.system("ls -lR")


