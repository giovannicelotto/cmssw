# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: --scenario pp --era Run3_2024 --customise Configuration/DataProcessing/Utils.addMonitoring --step NANO:@JME --conditions 150X_mcRun3_2024_realistic_v2 --datatier NANOAODSIM --eventcontent NANOEDMAODSIM1 --python_filename GEN-RunIII2024Summer24NanoAODv15-00385_1_cfg.py --fileout file:GEN-RunIII2024Summer24NanoAODv15-00385.root --filein dbs:/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM --number 1763 --number_out 1763 --no_exec --mc
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_2024_cff import Run3_2024
from FWCore.ParameterSet.VarParsing import VarParsing
import sys 

process = cms.Process('NANO',Run3_2024)
options = VarParsing('python')

#options.register('outputName', "",
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.string,
#    "output Name"
#)
options.register('collection', "candidate", VarParsing.multiplicity.singleton,VarParsing.varType.string,    "collection chosen")
options.register('SVLoose',     0.,         VarParsing.multiplicity.singleton,  VarParsing.varType.float,   "Threshold for SV score to enter in IVF")
options.register('SVSeed',      0.,         VarParsing.multiplicity.singleton,  VarParsing.varType.float,   "Threshold for SV score to be a seed for IVF")
options.register('SVEdge',      0.,         VarParsing.multiplicity.singleton,  VarParsing.varType.float,   "Threshold for edgeScore_ij of track j to be included in seed of track i")

options.register('f1',          0.7,         VarParsing.multiplicity.singleton,  VarParsing.varType.float,   "f1")
options.register('minS1',       2.,         VarParsing.multiplicity.singleton,  VarParsing.varType.float,   "minS1")
options.register('f2',          1.,         VarParsing.multiplicity.singleton,  VarParsing.varType.float,   "f2")
options.register('minS2',       0.,         VarParsing.multiplicity.singleton,  VarParsing.varType.float,   "minS2")

options.parseArguments()
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000),
    output = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( (
        '/store/mc/RunIII2024Summer24MiniAODv6/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/150X_mcRun3_2024_realistic_v2-v2/110000/01af4637-202e-4830-abe8-b129a6fc6a0c.root',
        #"/store/mc/RunIII2024Summer24MiniAODv6/QCD_Bin-PT-15to7000_Par-PT-Flat_TuneCH3_13p6TeV_herwig7/MINIAODSIM/NoPU_150X_mcRun3_2024_realistic_v2_ext2-v2/100000/003bdd8a-afd1-402d-ab83-c300d729b8d9.root"
        #'file:/work/gcelotto/btv_mini_rerun/CMSSW_15_0_6/src/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8.root'
     ) ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('--scenario nevts:1763'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

if options.outputFile=="output.root":
    if options.collection=="":
        processName = "outputBTV.root" 
    else:
        processName = options.collection+".root"
else:
    processName = options.outputFile
process.NANOEDMAODSIM1output = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:'+processName),
    outputCommands = cms.untracked.vstring(
        #'keep *',   
    'drop *',  # Drop everything by default
    "keep *_svTable_*_*",  # Keep event-level FlatTables
    "keep *_svCandidateTable_*_*",  # Keep event-level FlatTables
    'keep *_genVertexProducer_*_*',
    "keep *_genParticleTable_*_*",  # Keep the merged gen particles
    #"keep *_DummyTrackValueMap_*_*",
    "keep *_dummyValueMap_selectedTrackTable_*",
    "keep *_DummyTrackValueMap_selectedTrackTable_*",
    "keep *_gvProducer_*_*",
    "keep *_gvCentralProducer_*_*",
    "keep *_genCandidateVertexProducer_*_*",
    "keep nanoaodFlatTable_*Table*_*_*",  # Keep event-level FlatTables
    #"keep nanoaodFlatTable_svCandidateTable_*_*",  # for central
    #"keep nanoaodFlatTable_vertexTable_*_*",  # for central

    "keep nanoaodUniqueString_nanoMetadata_*_*",  # Keep basic metadata
    "keep nanoaodMergeableCounterTable_*_*_*",


)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '150X_mcRun3_2024_realistic_v2', '')

































sys.path.append("/work/gcelotto/btv_new/CMSSW_15_0_6/src/PhysicsTools/SimpleVertexTable/python")
process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")
from mergedGenParticles import mergedGenParticles, finalGenParticles, genParticleTable
process.mergedGenParticles = mergedGenParticles 
process.finalGenParticles = finalGenParticles 
process.genParticleTable = genParticleTable 



from genVertices import custom_GV_producer
from sv_candidate_cff import custom_sv_candidate #inclusiveCandidateVertexFinder, candidateVertexMerger, CandidateVertexArbitrator, myCandidateInclusiveSecondaryVertices, myvertexTable, svCandidateTable, slimmedSecondaryVertices
from sv_reco_cff import custom_sv_tracks
process = custom_GV_producer(process, collection=options.collection)

#print("\n===== GV PRODUCER DEBUG =====")
#print(process.gvProducer.dumpPython())
#print("================================\n")



print(options.collection)
if options.collection=="candidate":
    process = custom_sv_candidate(process)
    process.nanoAOD_step = cms.Path(
                                process.mergedGenParticles*
                                process.finalGenParticles*
                                process.genParticleTable+
                                process.sv_candidate*
                                process.genVertexProducer_sequence )

elif options.collection=="track":
    process = custom_sv_tracks(process, threshold_values=(options.SVLoose, options.SVSeed, options.SVEdge), f1=options.f1, minSig1=options.minS1, f2=options.f2, minSig2=options.minS2)
    process.nanoAOD_step = cms.Path(
                                process.mergedGenParticles*
                                process.finalGenParticles*
                                process.genParticleTable+
                                process.sv_track*
                                process.genVertexProducer_sequence )
elif options.collection=="central":
    #process.nanoTableTaskCommon = cms.Task(process.boostedTauTask, process.electronTask, process.jetForMETTask, process.jetPuppiForMETTask, process.jetPuppiTablesTask, process.jetPuppiTask, process.jetTablesTask, process.jetTask, process.linkedObjects, process.lowPtElectronTablesTask, process.lowPtElectronTask, process.metTablesTask, process.muonTablesTask, process.muonTask, process.nanoMetadata, process.photonTask, process.svCandidateTable, process.tauTablesTask, process.tauTask, process.vertexTable, process.vertexTask)
    #process.nanoSequenceCommon = cms.Sequence(process.nanoTableTaskCommon)
    #process.nanoTableTaskFS = cms.Task()
    process.nanoAOD_step = cms.Path(
                                process.mergedGenParticles*
                                process.finalGenParticles*
                                process.genParticleTable+
                                process.nanoSequenceMC*
                                process.genVertexProducer_sequence )
else:
    print("No optinos chosen")
















# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOEDMAODSIM1output_step = cms.EndPath(process.NANOEDMAODSIM1output)

# Schedule definition
process.schedule = cms.Schedule(
                                process.nanoAOD_step,    # central nanoAOD step already included in nanoSequenceMC
                                process.endjob_step,
                                process.NANOEDMAODSIM1output_step)







from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# Automatic addition of the customisation function from PhysicsTools.NanoAOD.nano_cff
from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeCommon 

#call to customisation function nanoAOD_customizeCommon imported from PhysicsTools.NanoAOD.nano_cff
process = nanoAOD_customizeCommon(process)

# Automatic addition of the customisation function from PhysicsTools.NanoAOD.custom_jme_cff
from PhysicsTools.NanoAOD.custom_jme_cff import PrepJMECustomNanoAOD 

#call to customisation function PrepJMECustomNanoAOD imported from PhysicsTools.NanoAOD.custom_jme_cff
process = PrepJMECustomNanoAOD(process)

# End of customisation functions


# Customisation from command line

process.source.delayReadingEventProducts = cms.untracked.bool(False)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
