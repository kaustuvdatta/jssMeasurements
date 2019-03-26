# jetObservables

This analysis uses extended NanoAOD format. More info about this format [here](https://twiki.cern.ch/twiki/bin/view/CMS/JetMET/JMARNanoAODv1)

To set up the code:
```bash
cmsrel  CMSSW_10_2_10
cd  CMSSW_10_2_10/src
cmsenv
git cms-merge-topic cms-nanoAOD:master-102X
git cms-merge-topic  cmantill:lsfinNanoAOD
git cms-merge-topic jmhogan:master-102X
git clone https://github.com/UBParker/NanoAODJMAR.git PhysicsTools/NanoAODJMAR
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
git clone https://github.com/cms-jet/NanoAODJMARTools.git PhysicsTools/NanoAODJMARTools
git clone ssh://git@gitlab.cern.ch:7999/asparker/QJetMass.git Analysis/QJetMass
git clone https://github.com/alefisico/jetObservables.git Analysis/jetObservables
ln -s $CMSSW_BASE/src/PhysicsTools/NanoAODTools/scripts/haddnano.py Analysis/jetObservables/test/
scram b -j 6
```

To run local:
```bash
cd Analysis/jetObservables/test/
python jetObservables_crab_extNanoAOD.py 1
```

To run in crab:
```bash
cd Analysis/jetObservables/test/
python multicrab.py -d SingleMuon -v v01
```
