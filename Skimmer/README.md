## Skimmer step

This step creates skim files for different selections: dijet, boosted W and top. The main files are [nSubProducer_withAllSel.py](python/nSubProducer_withAllSel.py) and [multicrab_nSubProducer.py](test/multicrab_nSubProducer.py)

### How to run skimmer

To run local:
```bash
cd jetObservables/Skimmer/test/
python jetObservables_nSubProducer.py --sample TTJets --local --selection dijet
```
The local option is to run in your local machine. Selection can be `dijet` or `Wtop`.

To run in crab:
```bash
cd jetObservables/Skimmer/test/
python multicrab_nSubProducer.py --datasets TTJets -v 106X_v01 --selection dijet
```
