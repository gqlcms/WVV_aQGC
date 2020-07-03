import FWCore.ParameterSet.Config as cms

import os
tempfilelist = ['Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFchs.txt','Summer16_07Aug2017_V11_MC_L2Relative_AK4PFchs.txt','Summer16_07Aug2017_V11_MC_L3Absolute_AK4PFchs.txt','Summer16_07Aug2017_V11_MC_L1FastJet_AK8PFchs.txt','Summer16_07Aug2017_V11_MC_L2Relative_AK8PFchs.txt','Summer16_07Aug2017_V11_MC_L3Absolute_AK8PFchs.txt','Summer16_07Aug2017_V11_MC_L1FastJet_AK8PFPuppi.txt','Summer16_07Aug2017_V11_MC_L2Relative_AK8PFPuppi.txt','Summer16_07Aug2017_V11_MC_L3Absolute_AK8PFPuppi.txt','Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFPuppi.txt','Summer16_07Aug2017_V11_MC_L2Relative_AK4PFPuppi.txt','Summer16_07Aug2017_V11_MC_L3Absolute_AK4PFPuppi.txt']
for i in tempfilelist:
  os.system("cp JECfile/"+i+" .")



process = cms.Process( "TEST" )
#process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True),allowUnscheduled=cms.untracked.bool(True))
#,
#				     SkipEvent = cms.untracked.vstring('ProductNotFound'))             
 
######## Sequence settings ##########
corrJetsOnTheFly = True
runOnMC  = True
runOnSig = False
#****************************************************************************************************#

#process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
if runOnMC:
   process.GlobalTag.globaltag = '94X_mcRun2_asymptotic_v3'#'MCRUN2_74_V9::All'
elif not(runOnMC):
   process.GlobalTag.globaltag = '94X_dataRun2_v10'

hltFiltersProcessName = 'RECO'
if runOnMC:
   hltFiltersProcessName = 'PAT' #'RECO'


process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False)
process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
                                                    inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
                                                    reverseDecision = cms.bool(False)
                                                    )
process.ApplyBaselineHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
                                                       inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
                                                       reverseDecision = cms.bool(False)
                                                       )
process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

baddetEcallist = cms.vuint32(
    [872439604,872422825,872420274,872423218,
     872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,
     872437052,872420649,872422436,872421950,
     872437185,872422564,872421566,872421695,
     872421955,872421567,872437184,872421951,
     872421694,872437056,872437057,872437313])


process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
    "EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallist, 
    taggingMode = cms.bool(True),
    debug = cms.bool(False)
    )


# DeepAK8: set up TransientTrackBuilder
process.load('Configuration.StandardSequences.MagneticField_cff')
process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName=cms.string('TransientTrackBuilder')
)
# ---------------------------------------------------------


option = 'RECO'

process.load("ExoDiBosonResonances.EDBRCommon.goodMuons_cff")
process.load("ExoDiBosonResonances.EDBRCommon.goodElectrons_cff")
process.load("ExoDiBosonResonances.EDBRCommon.goodJets_cff")
process.load("ExoDiBosonResonances.EDBRCommon.leptonicW_cff")
process.load("ExoDiBosonResonances.EDBRCommon.hadronicW_cff")
process.load("ExoDiBosonResonances.EDBRCommon.goodPuppi_cff")

from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
if not runOnMC:JETCorrLevels = ['L2Relative', 'L3Absolute', 'L2L3Residual']
if runOnMC:JETCorrLevels = ['L2Relative', 'L3Absolute']
jetToolbox(process, 'ak8', 'dummySeqAK8', 'noOutput',
           PUMethod='Puppi', JETCorrPayload='AK8PFPuppi', JETCorrLevels=JETCorrLevels,
           Cut='pt > 170.0 && abs(rapidity()) < 2.4',
           dataTier='miniAOD', runOnMC=runOnMC,
           addSoftDrop=True, addNsub=True,addSoftDropSubjets=True, subJETCorrPayload='AK4PFPuppi', subJETCorrLevels=JETCorrLevels   # must add soft-drop
)

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
from RecoBTag.MXNet.pfDeepBoostedJet_cff import _pfDeepBoostedJetTagsProbs,_pfDeepBoostedJetTagsMetaDiscrs,_pfMassDecorrelatedDeepBoostedJetTagsProbs,_pfMassDecorrelatedDeepBoostedJetTagsMetaDiscrs
#use the v2 training :  https://twiki.cern.ch/twiki/bin/viewauth/CMS/DeepAKXTagging#Additional_instructions_to_run_D
from RecoBTag.MXNet.pfDeepBoostedJet_cff import pfDeepBoostedJetTags, pfMassDecorrelatedDeepBoostedJetTags
from RecoBTag.MXNet.Parameters.V02.pfDeepBoostedJetPreprocessParams_cfi import pfDeepBoostedJetPreprocessParams as pfDeepBoostedJetPreprocessParamsV02
from RecoBTag.MXNet.Parameters.V02.pfMassDecorrelatedDeepBoostedJetPreprocessParams_cfi import pfMassDecorrelatedDeepBoostedJetPreprocessParams as pfMassDecorrelatedDeepBoostedJetPreprocessParamsV02
pfDeepBoostedJetTags.preprocessParams = pfDeepBoostedJetPreprocessParamsV02
pfDeepBoostedJetTags.model_path = 'RecoBTag/Combined/data/DeepBoostedJet/V02/full/resnet-symbol.json'
pfDeepBoostedJetTags.param_path = 'RecoBTag/Combined/data/DeepBoostedJet/V02/full/resnet-0000.params'
pfMassDecorrelatedDeepBoostedJetTags.preprocessParams = pfMassDecorrelatedDeepBoostedJetPreprocessParamsV02
pfMassDecorrelatedDeepBoostedJetTags.model_path = 'RecoBTag/Combined/data/DeepBoostedJet/V02/decorrelated/resnet-symbol.json'
pfMassDecorrelatedDeepBoostedJetTags.param_path = 'RecoBTag/Combined/data/DeepBoostedJet/V02/decorrelated/resnet-0000.params'

updateJetCollection(
    process,
    jetSource=cms.InputTag('packedPatJetsAK8PFPuppiSoftDrop'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    rParam=0.8,
    jetCorrections = ('AK8PFPuppi', cms.vstring(JETCorrLevels), 'None'),
    btagDiscriminators = _pfDeepBoostedJetTagsProbs + _pfDeepBoostedJetTagsMetaDiscrs+_pfMassDecorrelatedDeepBoostedJetTagsProbs + _pfMassDecorrelatedDeepBoostedJetTagsMetaDiscrs,
    postfix='AK8WithPuppiDaughters',   # !!! postfix must contain "WithPuppiDaughter" !!!
    printWarning = False
)

if option == 'RECO':
    process.goodMuons.src = "slimmedMuons"
    process.goodElectrons.src = "slimmedElectrons"
    process.Wtoenu.MET  = "slimmedMETs"
    process.Wtomunu.MET = "slimmedMETs"
    process.goodPuppiAK4.src = "slimmedJets"
    if runOnMC:
	process.goodPuppi.src = "JetUserData"
	met_recalculation = "JetUserDataak4"
    if  not runOnMC:
	process.goodPuppi.src = "selectedUpdatedPatJetsAK8WithPuppiDaughters"
	met_recalculation = "cleanPuppiAK4"









from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.leptonSequence = cms.Sequence(process.muSequence +
                                      process.egmGsfElectronIDSequence*process.eleSequence +
                                      process.leptonicVSequence )

process.jetSequence = cms.Sequence(
                                   process.fatPuppiSequence)






process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.load("RecoMET.METFilters.BadChargedCandidateFilter_cfi")
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.metfilterSequence = cms.Sequence(process.BadPFMuonFilter+process.BadChargedCandidateFilter)

######### JEC ########

jetsAK8puppi = "cleanPuppi"
 
if runOnMC:
   jecLevelsAK8puppi = [
                                   'Summer16_07Aug2017_V11_MC_L1FastJet_AK8PFPuppi.txt',
                                   'Summer16_07Aug2017_V11_MC_L2Relative_AK8PFPuppi.txt',
                                   'Summer16_07Aug2017_V11_MC_L3Absolute_AK8PFPuppi.txt'
     ]
   jecLevelsAK4chs = [
                                   'Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFchs.txt',
                                   'Summer16_07Aug2017_V11_MC_L2Relative_AK4PFchs.txt',
                                   'Summer16_07Aug2017_V11_MC_L3Absolute_AK4PFchs.txt'
    ]
else:
   jecLevelsAK8puppi = [
                                   'Summer16_07Aug2017BCD_V11_DATA_L1FastJet_AK8PFPuppi.txt',
                                   'Summer16_07Aug2017BCD_V11_DATA_L2Relative_AK8PFPuppi.txt',
                                   'Summer16_07Aug2017BCD_V11_DATA_L3Absolute_AK8PFPuppi.txt',
                                   'Summer16_07Aug2017BCD_V11_DATA_L2L3Residual_AK8PFPuppi.txt'
     ]
   jecLevelsAK4chs = [
                                   'Summer16_07Aug2017BCD_V11_DATA_L1FastJet_AK4PFPuppi.txt',
                                   'Summer16_07Aug2017BCD_V11_DATA_L2Relative_AK4PFPuppi.txt',
                                   'Summer16_07Aug2017BCD_V11_DATA_L3Absolute_AK4PFPuppi.txt',
                                   'Summer16_07Aug2017BCD_V11_DATA_L2L3Residual_AK4PFPuppi.txt'
     ]
if runOnMC:
	jLabel = "slimmedJetsAK8"
	jLabel = "selectedUpdatedPatJetsAK8WithPuppiDaughters"
	jetAlgo    = 'AK8PFPuppi'
	jer_era = "Summer16_07Aug2017_V11_MC"
	triggerResultsLabel      = "TriggerResults"
	triggerSummaryLabel      = "hltTriggerSummaryAOD"
	hltProcess = "HLT"
	process.JetUserData = cms.EDProducer(
	                                     'JetUserData',
	                                     jetLabel          = cms.InputTag(jLabel),
	                                     rho               = cms.InputTag("fixedGridRhoFastjetAll"),
	                                     coneSize          = cms.double(0.8),
	                                     getJERFromTxt     = cms.bool(False),
	                                     jetCorrLabel      = cms.string(jetAlgo),
	                                     jerLabel          = cms.string(jetAlgo),
	                                     resolutionsFile   = cms.string(jer_era+'_PtResolution_'+jetAlgo+'.txt'),
	                                     scaleFactorsFile  = cms.string(jer_era+'_SF_'+jetAlgo+'.txt'),
	                                     ### TTRIGGER ###
	                                     triggerResults = cms.InputTag(triggerResultsLabel,"",hltProcess),
	                                     triggerSummary = cms.InputTag(triggerSummaryLabel,"",hltProcess),
	                                     hltJetFilter       = cms.InputTag("hltPFHT"),
	                                     hltPath            = cms.string("HLT_PFHT800"),
	                                     hlt2reco_deltaRmax = cms.double(0.2),
	                                     candSVTagInfos         = cms.string("pfInclusiveSecondaryVertexFinder"),
	                                     jecAk8chsPayloadNames_jetUserdata = cms.vstring( jecLevelsAK8puppi ),
	                                     vertex_jetUserdata = cms.InputTag("offlineSlimmedPrimaryVertices"),
	                                     )
	jLabelak4 = "slimmedJets"
	jetAlgoak4    = 'AK4PFchs'
	jer_era = "Summer16_07Aug2017_V11_MC"
	triggerResultsLabel      = "TriggerResults"
	triggerSummaryLabel      = "hltTriggerSummaryAOD"
	hltProcess = "HLT"
	process.JetUserDataak4 = cms.EDProducer(
	                                        'JetUserDataak4',
	                                        jetLabel          = cms.InputTag(jLabelak4),
	                                        rho               = cms.InputTag("fixedGridRhoFastjetAll"),
	                                        coneSize          = cms.double(0.4),
	                                        getJERFromTxt     = cms.bool(False),
	                                        jetCorrLabel      = cms.string(jetAlgoak4),
	                                        jerLabel          = cms.string(jetAlgoak4),
	                                        resolutionsFile   = cms.string(jer_era+'_PtResolution_'+jetAlgoak4+'.txt'),
	                                        scaleFactorsFile  = cms.string(jer_era+'_SF_'+jetAlgoak4+'.txt'),
	                                        ### TTRIGGER ###
	                                        triggerResults = cms.InputTag(triggerResultsLabel,"",hltProcess),
	                                        triggerSummary = cms.InputTag(triggerSummaryLabel,"",hltProcess),
	                                        hltJetFilter       = cms.InputTag("hltPFHT"),
	                                        hltPath            = cms.string("HLT_PFHT800"),
	                                        hlt2reco_deltaRmax = cms.double(0.2),
	                                        candSVTagInfos         = cms.string("pfInclusiveSecondaryVertexFinder"),
	                                        jecAK4chsPayloadNames_JetUserData = cms.vstring( jecLevelsAK4chs ),
	                                        vertex_JetUserData = cms.InputTag("offlineSlimmedPrimaryVertices"),
	                                        )
#L1 Prefiring
from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
    DataEra = cms.string("2016BtoH"),   #("2016BtoH"), #Use 2016BtoH for 2016
    UseJetEMPt = cms.bool(False),
    PrefiringRateSystematicUncty = cms.double(0.2),
    SkipWarnings = False)
process.treeDumper = cms.EDAnalyzer("VVVTreeMaker",
                                    EDBRChannel = cms.string("VW_CHANNEL"),
                                    lhe =  cms.InputTag("externalLHEProducer"),
                                    isGen = cms.bool(False),
                                    isJEC = cms.bool(corrJetsOnTheFly),
                                    RunOnMC  = cms.bool(runOnMC),
                                    RunOnSig = cms.bool(runOnSig),
                                    generator =  cms.InputTag("generator"),
                                    genSrc =  cms.InputTag("prunedGenParticles"),
                                    pileup  =   cms.InputTag("slimmedAddPileupInfo"),
                                    leptonicVSrc = cms.InputTag("leptonicV"),
                                    t1jetSrc_userak4 = cms.InputTag(met_recalculation),
                                    looseMuonSrc = cms.InputTag("looseMuons"),
                                    looseElectronSrc = cms.InputTag("looseElectrons"),
                                    goodMuSrc = cms.InputTag("goodMuons"),
                                    MuSrc = cms.InputTag("slimmedMuons"),
                                    EleSrc = cms.InputTag("slimmedElectrons"),
                                    metSrc = cms.InputTag("slimmedMETs"),
                                    ak4jetsSrc = cms.InputTag("cleanPuppiAK4"), 
                                    jets = cms.InputTag("slimmedJets"),
                                    puppijets = cms.InputTag(jetsAK8puppi),
                                    jecAK4chsPayloadNames = cms.vstring( jecLevelsAK4chs ),					
                                    jecAK8puppiPayloadNames = cms.vstring( jecLevelsAK8puppi ),
                                    jecpath = cms.string(''),
                                    rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                    electronIDs = cms.InputTag("heepElectronID-HEEPV50-CSA14-25ns"),
                                    muons = cms.InputTag("slimmedMuons"),
                                    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                    hltToken    = cms.InputTag("TriggerResults","","HLT"),
                                    elPaths1     = cms.vstring("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v*"),#EXO-15-002#1
                                    elPaths2     = cms.vstring("HLT_Ele35_WPTight_Gsf_v*"), #B2G-15-005#2
                                    elPaths3     = cms.vstring("HLT_Ele45_WPLoose_Gsf_v*"),
                                    elPaths4     = cms.vstring("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"),#("HLT_Ele35_WPLoose_Gsf_v*"),#3
                                    elPaths5     = cms.vstring("HLT_Photon175_v*"),#("HLT_Ele35_WPLoose_Gsf_v*"),
                                    elPaths6     = cms.vstring("HLT_Ele38_WPTight_Gsf_v*"),#("HLT_Ele25_eta2p1_WPLoose_Gsf_v*"),
                                    elPaths7     = cms.vstring("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*"),#("HLT_Ele25_eta2p1_WPTight_Gsf_v*"),#4
                                    elPaths8     = cms.vstring("HLT_Ele27_WPTight_Gsf_v*"),
                                    muPaths1     = cms.vstring("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v*"),#EXO-15-002
                                    muPaths2     = cms.vstring("HLT_Mu50_v*"), #B2G-15-005
                                    muPaths3     = cms.vstring("HLT_TkMu50_v*"), #B2G-15-005
                                    muPaths4     = cms.vstring("HLT_OldMu100_v*"), #MET
                                    muPaths5     = cms.vstring("HLT_TkMu100_v*"), #MET
                                    muPaths6     = cms.vstring("HLT_PFHT900_v*"),
                                    muPaths7     = cms.vstring("HLT_PFJet450_v*"),
                                    muPaths8     = cms.vstring("HLT_PFJet500_v*"),
                                    muPaths9     = cms.vstring("HLT_AK8PFJet450_v*"),
                                    muPaths10     = cms.vstring("HLT_AK8PFJet500_v*"),
                                    muPaths11     = cms.vstring("HLT_AK8PFJet360_TrimMass30_v*"),#e1,MET110,e7,m1,m5-m11,HLT_hadronic
                                    muPaths12     = cms.vstring("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*"),
                                    noiseFilter = cms.InputTag('TriggerResults','', hltFiltersProcessName),
                                    noiseFilterSelection_HBHENoiseFilter = cms.string('Flag_HBHENoiseFilter'),
                                    noiseFilterSelection_HBHENoiseIsoFilter = cms.string("Flag_HBHENoiseIsoFilter"),
                                    noiseFilterSelection_GlobalTightHaloFilter = cms.string('Flag_globalSuperTightHalo2016Filter'),
                                    noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter = cms.string('Flag_EcalDeadCellTriggerPrimitiveFilter'),
                                    noiseFilterSelection_goodVertices = cms.string('Flag_goodVertices'),
                                    noiseFilterSelection_eeBadScFilter = cms.string('Flag_eeBadScFilter'),
                                    noiseFilterSelection_badMuon = cms.InputTag('BadPFMuonFilter'),
                                    noiseFilterSelection_badChargedHadron = cms.InputTag('BadChargedCandidateFilter'),
                                    CoutOrNot  = cms.bool(False)
                                    )


if option=='GEN':
    process.treeDumper.metSrc = 'genMetTrue'
    process.treeDumper.isGen  = True
 
if runOnMC:
	process.analysis = cms.Path(process.JetUserDataak4 +
                            process.JetUserData +
			    process.leptonSequence +
                            process.HBHENoiseFilterResultProducer+
                            process.ApplyBaselineHBHENoiseFilter+
                            process.ApplyBaselineHBHEIsoNoiseFilter+
                            process.jetSequence +
                            process.metfilterSequence +
                            process.ecalBadCalibReducedMINIAODFilter*process.prefiringweight*process.treeDumper)

if not runOnMC:
        process.analysis = cms.Path(
                            process.leptonSequence +
                            process.HBHENoiseFilterResultProducer+
                            process.ApplyBaselineHBHENoiseFilter+
                            process.ApplyBaselineHBHEIsoNoiseFilter+
                            process.jetSequence +
                            process.metfilterSequence +
                            process.ecalBadCalibReducedMINIAODFilter*process.prefiringweight*process.treeDumper)



process.load("ExoDiBosonResonances.EDBRCommon.data.RSGravitonToWW_kMpl01_M_1000_Tune4C_13TeV_pythia8")
#process.source.inputCommands = ['keep *','drop *_isolatedTracks_*_*']
process.source.fileNames = [
#'/store/mc/RunIISummer16MiniAODv3/WkkToWRadionToWWW_M5000-R0-7-TuneCUETP8M1_13TeV-madgraph-pythia/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/280000/DCB4F5A0-753F-EA11-90C3-509A4C8339DE.root'
#'root://eosuser.cern.ch//eos/user/q/qiguo/aQGC/WWW_private/WWW_aQGC_miniAOD_v1/crab_WWW_aQGC_miniAOD_v1/200504_084734/0000/SUS-RunIISummer16MiniAODv3-00241_3.root, root://eosuser.cern.ch//eos/user/q/qiguo/aQGC/WWW_private/WWW_aQGC_miniAOD_v1/crab_WWW_aQGC_miniAOD_v1/200504_084734/0000/SUS-RunIISummer16MiniAODv3-00241_2.root, root://eosuser.cern.ch//eos/user/q/qiguo/aQGC/WWW_private/WWW_aQGC_miniAOD_v1/crab_WWW_aQGC_miniAOD_v1/200504_084734/0000/SUS-RunIISummer16MiniAODv3-00241_1.root'
#'root://eosuser.cern.ch//eos/user/q/qiguo/aQGC/WWW_private/WWW_aQGC_miniAOD_v1/crab_WWW_aQGC_miniAOD_v1/200504_084734/0000/SUS-RunIISummer16MiniAODv3-00241_3.root',' root://eosuser.cern.ch//eos/user/q/qiguo/aQGC/WWW_private/WWW_aQGC_miniAOD_v1/crab_WWW_aQGC_miniAOD_v1/200504_084734/0000/SUS-RunIISummer16MiniAODv3-00241_2.root', 'root://eosuser.cern.ch//eos/user/q/qiguo/aQGC/WWW_private/WWW_aQGC_miniAOD_v1/crab_WWW_aQGC_miniAOD_v1/200504_084734/0000/SUS-RunIISummer16MiniAODv3-00241_1.root'
#'file:/eos/user/q/qiguo/aQGC/condor/private_sample/WZZ/v1/miniaod_8428521_0.root'#,'file:/eos/user/q/qiguo/temp/miniaod_8350145_1.root','file:/eos/user/q/qiguo/temp/miniaod.root'
#'file:/eos/user/q/qiguo/aQGC/condor/private_sample/WWW/v2/miniaod_8764611_28.root'#,'file:/eos/user/q/qiguo/temp/miniaod_8350145_1.root','file:/eos/user/q/qiguo/temp/miniaod.root'
'/store/mc/RunIISummer16MiniAODv3/WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/20000/5473D0C8-A0F0-E811-9ABF-003048FFCC16.root'
#'/store/mc/RunIISummer16MiniAODv3/WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/20000/366CD306-71F0-E811-910D-0CC47A745298.root'
]

process.maxEvents.input = 20000
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.FwkReport.limit = 99999999
print "hh"
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("RStreeEDBR_pickup_NLO_2.root")
                                   )


