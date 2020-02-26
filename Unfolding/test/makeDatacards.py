#####################################
#####################################

#!/usr/bin/env python
import argparse, os, shutil, sys
from ROOT import *
import ROOT
from datasets import *
from DrawHistogram import plotSimpleComparison
#import CMS_lumi as CMS_lumi
#18 import tdrstyle as tdrstyle
####gReset()
gROOT.SetBatch()
#gROOT.ForceStyle()
#tdrstyle.setTDRStyle()
#gStyle.SetOptStat(0)


def createDatacards( dataFile, sigFiles, bkgFiles, variables ):
    """docstring for createDataCards"""

    ### Getting input histos
    dataHistos = loadHistograms( dataFile, variables, isMC=False )
    dataHisto = dataHistos['datarecoJet']

    MCHisto = loadHistograms( sigFiles, variables )
    MCHisto['ttbarMadgraphrecoJet'].Scale(  dataHisto.Integral()/MCHisto['ttbarMadgraphrecoJet'].Integral())
    MCHisto['ttbarMadgraphgenJet'].Scale(  dataHisto.Integral()/MCHisto['ttbarMadgraphgenJet'].Integral())
    MCHisto['ttbarMadgraphrespJet'].Scale(  dataHisto.Integral()/MCHisto['ttbarMadgraphrespJet'].Integral())

    #print '----> Using RooUnfold'
    #R = ROOT.RooUnfoldResponse( MCHisto['recoJet'].Clone(), MCHisto['genJet'].Clone(), MCHisto['respJet'].Clone()  )
    #unfold = ROOT.RooUnfoldInvert( R, dataHisto )
    #hUnf = unfold.Hreco()

    ######## Cross check
    plotSimpleComparison( dataHisto, 'ttbarPower', MCHisto['ttbarMadgraphrecoJet'], 'ttbarMadgraph', 'tau21', rebinX=variables['Tau21'][0], version=args.version  )

    #can2D = TCanvas('can2D', 'can2D', 750, 800 )
    #MCHisto['respJet'].Draw("colz")
    #can2D.SaveAs('test2D.png')

    for ivar in variables:
        datacardName = 'datacard_'+ivar
        datacard = open( datacardName+'.txt', 'w')
        datacard.write("* imax\n")
        datacard.write("* jmax\n")
        datacard.write("* kmax\n")
        datacard.write("----------------\n")
        datacard.write("bin ")
        for ireco in range(1, dataHisto.GetNbinsX()+1 ):
            if dataHisto.GetBinContent(ireco)<1: continue
            #if ireco<5 or ireco>15: continue
            datacard.write("Reco_"+str(ireco)+" ")
        datacard.write("\n")

        datacard.write("observation ")
        for ireco in range(1, dataHisto.GetNbinsX()+1 ):
            if dataHisto.GetBinContent(ireco)<1: continue
            #if ireco<5 or ireco>15: continue
            datacard.write( str( int(dataHisto.GetBinContent(ireco)))+" " )
        datacard.write("\n")
        datacard.write("----------------\n")

        cleanup=True
        datacard.write("bin ")
        for ireco in range(1, MCHisto['ttbarMadgraphrecoJet'].GetNbinsX()+1):
            if dataHisto.GetBinContent(ireco)<1: continue
            #if ireco<5 or ireco>15: continue
            for igen in range(1, MCHisto['ttbarMadgraphgenJet'].GetNbinsX()+1):
                # remove un-necessary processes
                if cleanup and MCHisto['ttbarMadgraphrespJet'].GetBinContent(igen,ireco)<1: continue
                datacard.write("Reco_%d "%ireco) ##sig igen, in reco ireco
            #datacard.write("Reco_%d "%ireco)## bkg
        datacard.write("\n")

        datacard.write("process ")
        for ireco in range(1, MCHisto['ttbarMadgraphrecoJet'].GetNbinsX()+1):
            if dataHisto.GetBinContent(ireco)<1: continue
            #if ireco<5 or ireco>15: continue
            for igen in range(1, MCHisto['ttbarMadgraphgenJet'].GetNbinsX()+1):
                if cleanup and MCHisto['ttbarMadgraphrespJet'].GetBinContent(igen,ireco)<1: continue
                datacard.write("Gen_%d "%igen) ##sig igen, in reco ireco
            #datacard.write("Bkg ")## bkg
        datacard.write("\n")

        datacard.write("process ")
        for ireco in range(1, MCHisto['ttbarMadgraphrecoJet'].GetNbinsX()+1):
            if dataHisto.GetBinContent(ireco)<1: continue
            #if ireco<5 or ireco>15: continue
            for igen in range(1, MCHisto['ttbarMadgraphgenJet'].GetNbinsX()+1):
                if cleanup and MCHisto['ttbarMadgraphrespJet'].GetBinContent(igen,ireco)<1: continue
                datacard.write("%d "%(-igen)) ## 0 -1, -2 --> for signal
            #datacard.write("1 ")## bkg >0 for bkg
        datacard.write("\n")

        datacard.write("rate ")
        for ireco in range(1, MCHisto['ttbarMadgraphrecoJet'].GetNbinsX()+1):
            if dataHisto.GetBinContent(ireco)<1: continue
            #if ireco<5 or ireco>15: continue
            for igen in range(1, MCHisto['ttbarMadgraphgenJet'].GetNbinsX()+1):
                if cleanup and MCHisto['ttbarMadgraphrespJet'].GetBinContent(igen,ireco)<1: continue
                datacard.write("%.2f "% ( MCHisto['ttbarMadgraphrespJet'].GetBinContent(igen,ireco)))
            #datacard.write("1 ")
            #datacard.write("%.2f "%(MCHisto['recoJetScaled'].GetBinContent(ireco) if MCHisto['recoJetScaled'].GetBinContent(ireco)>0 else 1))##
        datacard.write("\n")
        datacard.write("----------------\n")

        po=' '.join(["--PO map='.*Gen_%d:r_bin%d[1,-1,20]'"%(igen,igen) for igen in range(1, dataHisto.GetNbinsX()+1) if dataHisto.GetBinContent(igen)>0 ])
        cmdText2Workspace = "text2workspace.py --X-allow-no-background "+datacardName+".txt -o "+datacardName+".root -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel "+po
        datacard.write("### RUN WITH COMMANDS: ####\n")
        datacard.write("# "+cmdText2Workspace+"\n")
        print cmdText2Workspace
        cmdCombine = 'combine -M MultiDimFit --algo singles -d "+datacardName+".root -t 0'
        datacard.write("# Run with command: "+cmdCombine+" \n")
        datacard.write("############################\n")

##########################################################################
def loadHistograms( samples, variables, isMC=True ):
    """docstring for loadHistograms"""

    allHistos = {}
    for var in variables:
        for isam in samples:
            for ih in [ 'recoJet', 'genJet', 'respJet' ]:
                allHistos[isam+ih] = samples[isam][1].Get( 'jetObservables/'+ih+var )
                #if isMC:
                MCScale = checkDict( samples[isam][0], dictSamples )['XS'] * args.lumi / checkDict( samples[isam][0], dictSamples )['2016']['nevents']
                allHistos[isam+ih].Scale( MCScale )

                if len(variables[var])==1:
                    if not ih.startswith('resp'): allHistos[isam+ih].Rebin( variables[var][0] )
                    else: allHistos[isam+ih].Rebin2D( variables[var][0], variables[var][0] )

    return allHistos



###########################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-p", "--process", action='store', dest="process", default="datacard", help="Process: datacard or plot." )
    parser.add_argument("-f", "--inputFolder", action='store', dest="inputFolder", default="", help="Path of the folder where the root files are located (wihtout the root://blah//)." )
    parser.add_argument("-o", "--outputFolder", action='store', dest="outputFolder", default="test", help="Name of output folder" )
    parser.add_argument('-l', '--lumi', action='store', type=float, default=41530., help='Luminosity, example: 1.' )
    parser.add_argument("-v", "--version", action='store', dest="version", default="v00", help="Version" )

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    dataFile = {}
    #dataFile['data'] = [ 'data', TFile( 'Rootfiles/jetObservables_histograms_SingleMuon2016ALL.root' ), 'Data', 'kBlack' ]
    dataFile['data'] = [ 'TT_TuneCUETP8M2T4_13TeV-powheg-pythia8', TFile( 'Rootfiles/jetObservables_histograms_TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root' ), 'Data', 'kBlack' ]

    sigFiles = {}
    sigFiles['ttbarMadgraph'] = [ 'TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8', TFile( 'Rootfiles/jetObservables_histograms_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root'), 'ttbar (madgraph)', 'kBlue' ]
    #sigFiles['ttbarPowheg'] = [ 'TT_TuneCUETP8M2T4_13TeV-powheg-pythia8', TFile( 'Rootfiles/jetObservables_histograms_TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root'), 'ttbar (madgraph)', 'kBlue' ]

    bkgFiles = {}
    bkgFiles['ST_s'] = [ 'ST_s-channel_4f_InclusiveDecays_13TeV-amcatnlo-pythia8', TFile( 'Rootfiles/jetObservables_histograms_ST_s-channel_4f_InclusiveDecays_13TeV-amcatnlo-pythia8.root' ), 'Single top', 'kMagenta' ]
    bkgFiles['ST_t_antitop'] = [ 'ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1', TFile( 'Rootfiles/jetObservables_histograms_ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root' ), 'Single top', 'kMagenta' ]
    bkgFiles['ST_t_top'] = [ 'ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1', TFile( 'Rootfiles/jetObservables_histograms_ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root' ), 'Single top', 'kMagenta' ]
    bkgFiles['ST_tW_antitop'] = [ 'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4', TFile( 'Rootfiles/jetObservables_histograms_ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4.root' ), 'Single top', 'kMagenta' ]
    bkgFiles['ST_tW_top'] = [ 'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4', TFile( 'Rootfiles/jetObservables_histograms_ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4.root' ), 'Single top', 'kMagenta' ]
    bkgFiles['WJets'] = [ 'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8', TFile( 'Rootfiles/jetObservables_histograms_WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root' ), 'WJets', 'kCyan' ]
    ##bkgFiles[] = [ '', TFile( 'Rootfiles/jetObservables_histograms_'+ibkg+'.root' ), '', 'kMagenta' ]

    variables = {}
    variables[ 'Tau21' ] = [ 4 ]

    createDatacards( dataFile, sigFiles, bkgFiles, variables )

