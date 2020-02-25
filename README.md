# jetObservables


This analysis uses extended NanoAOD format. More info about this format [here](https://twiki.cern.ch/twiki/bin/view/CMS/JetMET/JMARNanoAODv1)

## Instructions

### Basic setup 
To set up the code:
```bash
cmsrel  CMSSW_10_6_5
cd  CMSSW_10_6_5/src
cmsenv
git cms-addpkg GeneratorInterface/Core
git clone https://github.com/UBParker/NanoAODJMAR.git PhysicsTools/NanoAODJMAR
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
git clone https://github.com/cms-jet/NanoAODJMARTools.git PhysicsTools/NanoAODJMARTools
git clone ssh://git@gitlab.cern.ch:7999/asparker/QJetMass.git Analysis/QJetMass    ### This is just as example
git clone https://github.com/alefisico/jetObservables.git -b 106X Analysis/jetObservables
ln -s $CMSSW_BASE/src/PhysicsTools/NanoAODTools/scripts/haddnano.py Analysis/jetObservables/test/
scram b -j 6
```

### How to run skimmer

To run local:
```bash
cd Analysis/jetObservables/test/
python jetObservables_nSubProducer.py --sample TTJets --local
```

To run in crab:
```bash
cd Analysis/jetObservables/test/
python multicrab_nSubProducer.py --datasets TTJets -v 106X_v01
```

### How to compute cross sections

We are using a modified version of [these instructions](https://twiki.cern.ch/twiki/bin/viewauth/CMS/HowToGenXSecAnalyzer#Automated_scripts_to_compute_the). To set up the code, if you did not include it before:

```bash
cd $CMSSW_BASE/src/
git cms-addpkg GeneratorInterface/Core
scram b -j 8
```

To run:
```bash
cd Analysis/jetObservables/test/
git clone https://github.com/cms-sw/genproductions.git
mv genproductions/test/calculateXSectionAndFilterEfficiency/ .
rm -rf genproductions
```

To run:
```bash
cd $CMSSW_BASE/src/Analysis/jetObservables/test/calculateXSectionAndFilterEfficiency/
./calculateXSectionAndFilterEfficiency.sh -f datasets.txt -c Moriond17 -d MINIAODSIM -n 1000000  ## run using list of dataset names mode
```
where the input parameters are:
 * -f wants the input file containing the list of dataset names (default) or McM prepID (requires -m)
 * -c specifies the campaign, i.e. the string to be used to search for the secondary dataset name /.../*Moriond17*/*
 * -d specifies the datatier to be used, i.e.  /.../*/MINIAODSIM
 * -n number of events to be used for each dataset to compute the cross section
 * -m use the McM prepID instead of the dataset names

