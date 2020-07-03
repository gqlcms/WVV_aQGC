from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName   = '16_WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8_3'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'

tempfilelist = config.JobType.inputFiles = ['Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFchs.txt','Summer16_07Aug2017_V11_MC_L2Relative_AK4PFchs.txt','Summer16_07Aug2017_V11_MC_L3Absolute_AK4PFchs.txt','Summer16_07Aug2017_V11_MC_L1FastJet_AK8PFchs.txt','Summer16_07Aug2017_V11_MC_L2Relative_AK8PFchs.txt','Summer16_07Aug2017_V11_MC_L3Absolute_AK8PFchs.txt','Summer16_07Aug2017_V11_MC_L1FastJet_AK8PFPuppi.txt','Summer16_07Aug2017_V11_MC_L2Relative_AK8PFPuppi.txt','Summer16_07Aug2017_V11_MC_L3Absolute_AK8PFPuppi.txt','Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFPuppi.txt','Summer16_07Aug2017_V11_MC_L2Relative_AK4PFPuppi.txt','Summer16_07Aug2017_V11_MC_L3Absolute_AK4PFPuppi.txt']
filepath ='/afs/cern.ch/work/q/qiguo/B2G/work/04/04_17/16code/CMSSW_10_2_18/src/ExoDiBosonResonances/EDBRTreeMaker/test/'
filelist = []
for i in tempfilelist :
   filelist.append(filepath+i)

config.JobType.inputFiles = filelist
#config.JobType.inputFiles = ['Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs.txt','Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs.txt','Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFchs.txt','Fall17_17Nov2017_V32_MC_L1FastJet_AK8PFchs.txt','Fall17_17Nov2017_V32_MC_L2Relative_AK8PFchs.txt','Fall17_17Nov2017_V32_MC_L3Absolute_AK8PFchs.txt','Fall17_17Nov2017_V32_MC_L1FastJet_AK8PFPuppi.txt','Fall17_17Nov2017_V32_MC_L2Relative_AK8PFPuppi.txt','Fall17_17Nov2017_V32_MC_L3Absolute_AK8PFPuppi.txt','Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFPuppi.txt','Fall17_17Nov2017_V32_MC_L2Relative_AK4PFPuppi.txt','Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFPuppi.txt','L1PrefiringMaps_new.root']

config.JobType.psetName    = 'analysis.py'

config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob =1
config.Data.totalUnits = -1
config.Data.publication = False

config.Data.outLFNDirBase='/store/user/qiguo/crab/VVV/16sample'

config.Data.outputDatasetTag = '16_WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8_3'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERNBOX'
