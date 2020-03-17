import ROOT
import math, os, sys
import array
import glob
import numpy as np
import root_numpy as rtnpy
import h5py


#This class aims to make loading of datasets from samples, application of further cuts after the nSub skimmer, etc. 
#more convenient 
class nSubExtractor:
	#M-body variable encodes max N for tau_N, 
	#selections available are only for boosted W, but more to be added soon


	def __init__(self, sample="TTbar/1", isMC=1, selection='W'):

		self.selection = selection
		self.nSub_labels = ["_tau_0p5_", "_tau_1_", "_tau_2_"]
		self.maxTau = 7
		self.axisdef = "Echeme_excl_kT"
		
		#self.cuts = []
		self.sample_names = ["TTbar/1", "ST/1", "ST/2", "ST/3", "ST/4", "ST/5", "Wjets/2",
							 "Data/B", "Data/C", "Data/D", "Data/E", "Data/F", "Data/G", "Data/H"]

		self.genweights = 1.
		self.lumi = 5.75+2.57+4.24+4.03+3.11+7.57+8.65 #B+C+D+E+F+G+H
		
		# 'k-factor'*pb-1->fb-1*lumi(2016)*list of cross-sections in pb^-1
		# for "TTbar/1", "ST/1", "ST/2", "ST/3", "ST/4", "ST/5", "Wjets/2"
		self.xs_weights = 0.8*1000*lumi*np.array([831.76, (80.95*0.322), (136.02*0.322), 35.6, 35.6, 10.12, 60781.5]) 
		self.nevents = np.array([ 76915549, 38811017, 66960888, 998276, 992024, 2989199, 158307515])
				

		if self.selection=="W":	
			self.ptmin = 200.	
			self.mSDmin = 65.
			self.mSDmax = 105.
			self.leptWpTmin = 200.
		

		#self.sample = sample# ["TTbar/1", "ST/1", "ST/2", "ST/3", "ST/4", "ST/5", "Wjets/2", Data/B", "Data/C", "Data/D", "Data/E", "Data/F", "Data/G", "Data/H"]
		
		if isMC: self.sample_list = sample

		else: self.sample_list = [i for i in self.sample_names if i.startswith(sample)]

		self.filepaths = ['/work/kadatta/Unfolding_samples/Data/%s/jetObservables_nanoskim_*.root'%sample for sample in self.sample_list]

	def create_var_sel_list():
		
		tau_reco = []
		if isMC: tau_gen = []

		for tauN in range(self.maxTau):
			for x in self.nSub_labels:

		        if tau_N==0 and x=="_tau_0p5_": 
		        	tau_gen[0] ="goodgenjet0"+x+str(tauN)
		        	tau_reco[0] ="goodgrecojet0"+x+str(tauN)
		        	
		        else:
		        	tau_gen.append("goodgenjet0"+x+str(tauN))
		        	tau_reco.append("goodgrecojet0"+x+str(tauN))


		if isMC:
		    var_list = ['goodrecojet0_softdrop_mass', 'goodrecojet0_mass', 
		    			'goodrecojet0_pt', 'goodrecojet0_eta', 'goodrecojet0_phi',
	                    'leptonicW_pT', 'lepton_pT', 'puWeight', 'PV_npvsGood', 'btagWeight_CSVV2',
	                    'dr_LepJet', 'dphi_MetJet', 'dphi_WJet','genmatchedrecoAK8', 'passedMETfilters',
	                    'goodgenjet0_mass',
	                    tau_reco.flatten(),
	                    tau_gen.flatten(),
	                    'goodgenjet0_pt', 'goodgenjet0_eta', 'goodgenjet0_phi',]

	        if 'Wjets' in sample or 'Wj' in sample:
	        	var_list.append('genWeight')

	        selection_indices = [var_list.index('goodrecojet0_softdrop_mass'),
                    			 var_list.index('goodrecojet0_pt'),
                  				 var_list.index('leptonicW_pT'),
                   				 var_list.index('dr_LepJet'),
                    			 var_list.index('dphi_MetJet'),
                     			 var_list.index('dphi_WJet'),
			        			 var_list.index('passedMETfilters')]
			sel_i = selection_indices
			tau_reco_ind = var_list.index('goodrecojet0_tau_0p5_0')
			tau_gen_ind = var_list.index('goodgenjet0_tau_0p5_0')
	
			return var_list, sel_i, tau_reco_ind, tau_gen_ind

	    else: 
	        var_list = ['goodrecojet0_softdrop_mass', 'goodrecojet0_mass', 
	        			'goodrecojet0_pt', 'goodrecojet0_eta', 'goodrecojet0_phi', 
	                    'leptonicW_pT', 'lepton_pT', 'PV_npvsGood',
	                    'dr_LepJet', 'dphi_MetJet', 'dphi_WJet', 'passedMETfilters',
	                    tau_reco.flatten(),]
	                    
	        selection_indices = [var_list.index('goodrecojet0_softdrop_mass'),
                    			 var_list.index('goodrecojet0_pt'),
                  				 var_list.index('leptonicW_pT'),
                   				 var_list.index('dr_LepJet'),
                    			 var_list.index('dphi_MetJet'),
                     			 var_list.index('dphi_WJet'),
			        			 var_list.index('passedMETfilters')]
			sel_i = selection_indices
			tau_reco_ind = var_list.index('goodrecojet0_tau_0p5_0')

	    	return var_list, sel_i, tau_reco_ind
	    return -1

	def sample_loader_MC():

		filelist=[]

		for path in self.filepaths:
			files=glob.glob(path)
			for f in files:
				filelist.append(f)

		dataset = [] 
		c=0
		x = 0 

		var_list, sel_i, tau_reco_ind, tau_gen_ind = create_var_sel_list()

		for f in filelist:
		    dataset = []

		    F = ROOT.TFile.Open(f, 'read')
		    
		    T = F.Get("Events");
		    
		    dataset = rtnpy.tree2array(T, branches=var_list)
		    
		    evt_list = []
		    for i in range(0,dataset.shape[0]):

		        if mSDmin<dataset[i][sel_i[0]]<=mSDmax and ptmin<dataset[i][sel_i[1]]<ptmax and leptWpTmin<dataset[i][sel_i[2]] and dataset[i][sel_i[3]]>np.pi/2 and dataset[i][sel_i[4]]>2. and dataset[i][sel_i[5]]>2. and dataset[i][sel_i[6]]==1:
		            evt_list.append(i)
		    
		    if c==0:
		        reco_nSub_basis = np.ones((len(evt_list),21))
		        gen_nSub_basis = np.ones((len(evt_list),21))
		        
		    else:
		        reco_nSub_basis = np.concatenate((reco_nSub_basis, np.ones((len(evt_list),21))))
		        gen_nSub_basis = np.concatenate((gen_nSub_basis, np.ones((len(evt_list),21))))

		    #print ("on file %d"%c), x, len(evt_list), reco_nSub_basis.shape[0]
		    for i in evt_list:
		        y=0
		        for z in range(tau_reco_ind,tau_reco_ind+21):
		            reco_nSub_basis[x][y] = dataset[i][z]
		            y = y+1
		        y=0
		        for z in range(tau_gen_ind,tau_gen_ind+21):
		            gen_nSub_basis[x][y] = dataset[i][z]   
		            y = y+1
		        x = x+1
		    if c==0:
		        data = dataset[evt_list]
		    else:
		        data = np.concatenate((data, dataset[evt_list]))
		        
		    c=c+1

		#if isMC: 
		return data, reco_nSub_basis, gen_nSub_basis

		#else: return data, reco_nSub_basis