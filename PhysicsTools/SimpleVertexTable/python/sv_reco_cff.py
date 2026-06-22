import FWCore.ParameterSet.Config as cms

#packedPFcandidates are filtered (if(c.hasTrackDetails() && c.charge() != 0 && c.numberOfHits()> 0))
unpackedTracksAndVertices = cms.EDProducer('PATTrackAndVertexUnpacker',
    slimmedVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    slimmedSecondaryVertices = cms.InputTag("slimmedSecondaryVertices"),
    additionalTracks = cms.InputTag("lostTracks"),
    packedCandidates = cms.InputTag("packedPFCandidates")
)

dummyValueMap = cms.EDProducer("DummyTrackValueMap",
    src = cms.InputTag("unpackedTracksAndVertices"),
    pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    model_path = cms.FileInPath("PhysicsTools/data/submod_out128_hyper_1802.onnx"),
    threshold = cms.double(-10.) # this needs to be < zero always for the moment
)


# IVF parameter : https://github.com/cms-sw/cmssw/blob/55251374c7e82ee5ee7626de6248007aec863e1c/RecoVertex/AdaptiveVertexFinder/python/inclusiveVertexFinder_cfi.py#L15C1-L16C49
inclusiveVertexFinder = cms.EDProducer('InclusiveVertexFinder',
  svScores = cms.InputTag("dummyValueMap", "SVscore"),
  edgeScores = cms.InputTag("dummyValueMap", "edgeScores"),
  edgeIndices = cms.InputTag("dummyValueMap", "edgeIndices"),
  svScoreThreshold = cms.double(0.1),   # minimal SV score to enter in IVF
  seedScoreThreshold = cms.double(0.8), # minimal SV score to be a seed for IVF
  edgeScoreThreshold = cms.double(0.7),  # minimal edgeScore_ij of track j to be included in seed of track i
  beamSpot = cms.InputTag('offlineBeamSpot'),
  clusterizer = cms.PSet(
    #clusterMaxDistance = cms.double(0.1), # 0.05 default | requirement on adding a track to the cluster
    #clusterMaxSignificance = cms.double(9.), # 4.5 default |requirement on adding a track to the cluster
    #clusterMinAngleCosine = cms.double(0.25), # -2 to disable. | 0.5 default requirement on adding a track to the cluster
    #distanceRatio = cms.double(10.),
    clusterMaxDistance = cms.double(0.1),
    clusterMaxSignificance = cms.double(9.0),
    clusterMinAngleCosine = cms.double(0.25),
    distanceRatio = cms.double(10.),
    maxTimeSignificance = cms.double(3.5),
    seedMax3DIPSignificance = cms.double(9999), #disabled
    seedMax3DIPValue = cms.double(9999), #disabled
    #seedMin3DIPSignificance = cms.double(0.6), # which tracks can start a cluster
    #seedMin3DIPValue = cms.double(0.0025),  # which tracks can start a cluster
    seedMin3DIPSignificance = cms.double(0.6),
    seedMin3DIPValue = cms.double(0.0025),
  ),
  fitterRatio = cms.double(0.25),
  fitterSigmacut = cms.double(3),
  fitterTini = cms.double(256),
  maxNTracks = cms.uint32(30),
  maximumLongitudinalImpactParameter = cms.double(0.3),
  maximumTimeSignificance = cms.double(3), # new?
  minHits = cms.uint32(4), #8
  minPt = cms.double(0.4),
  primaryVertices = cms.InputTag('unpackedTracksAndVertices'),
  #tracks = cms.InputTag('dummyValueMap', 'selectedTracks'),
  tracks = cms.InputTag('unpackedTracksAndVertices'),
  useDirectVertexFitter = cms.bool(True),
  useVertexReco = cms.bool(True),
  vertexMinAngleCosine = cms.double(0.95),
  vertexMinDLen2DSig = cms.double(2.5),
  vertexMinDLenSig = cms.double(0.5),
  vertexReco = cms.PSet(
    finder = cms.string('avr'),
    primcut = cms.double(1),
    seccut = cms.double(3),
    smoothing = cms.bool(True)
  ),
)

#Vertex Merger step1 https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_X/RecoVertex/AdaptiveVertexFinder/python/vertexMerger_cfi.py
vertexMerger = cms.EDProducer( "VertexMerger",
    secondaryVertices = cms.InputTag("inclusiveVertexFinder"),  
    maxFraction = cms.double(0.7), 
    minSignificance = cms.double(2.0)
)

#Arbitrator step
# https://github.com/cms-sw/cmssw/blob/55251374c7e82ee5ee7626de6248007aec863e1c/RecoVertex/AdaptiveVertexFinder/python/trackVertexArbitrator_cfi.py#L16
trackVertexArbitrator = cms.EDProducer("TrackVertexArbitrator",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    primaryVertices = cms.InputTag("unpackedTracksAndVertices"),
    tracks = cms.InputTag("unpackedTracksAndVertices"),
    secondaryVertices = cms.InputTag("vertexMerger"),
    dLenFraction = cms.double(0.333),
    dRCut = cms.double(0.4),
    distCut = cms.double(0.04),
    sigCut = cms.double(5),
    fitterSigmacut =  cms.double(3),
    fitterTini = cms.double(256),
    fitterRatio = cms.double(0.25),
    trackMinLayers = cms.int32(4),
    trackMinPt = cms.double(0.4),
    trackMinPixels = cms.int32(1)
    # plus any additional parameters it requires
)

#Vertex Merger step2 https://github.com/cms-sw/cmssw/blob/557f39bce1d5cba35316c2358a89e888901a07e5/RecoVertex/AdaptiveVertexFinder/python/inclusiveVertexing_cff.py#L7
myFinalInclusiveSecondaryVertices = vertexMerger.clone(
    secondaryVertices = "trackVertexArbitrator",
    maxFraction = cms.double(1.0), #0.2 default
    minSignificance = cms.double(0.) ) #10 default


svTable = cms.EDProducer("SVTableProducer", 
                        pvSrc=cms.InputTag("offlineSlimmedPrimaryVertices"),
                        src = cms.InputTag("myFinalInclusiveSecondaryVertices"),
                        dlenSigMin = cms.double(3.0))


# Missing cut in dlen and dlenSig
# Missing cut in dlen and dlenSig
# Missing cut in dlen and dlenSig



def custom_sv_tracks(process, threshold_values=(0.0, 0., 0.)):
  process.unpackedTracksAndVertices = unpackedTracksAndVertices
  process.inclusiveVertexFinder = inclusiveVertexFinder.clone(
    svScoreThreshold=cms.double(threshold_values[0]),
    seedScoreThreshold=cms.double(threshold_values[1]),
    edgeScoreThreshold=cms.double(threshold_values[2]))
  print(f"Using custom thresholds for SV reconstruction: svScoreThreshold={threshold_values[0]}, seedScoreThreshold={threshold_values[1]}, edgeScoreThreshold={threshold_values[2]}")
  process.vertexMerger = vertexMerger
  process.trackVertexArbitrator = trackVertexArbitrator
  process.myFinalInclusiveSecondaryVertices = myFinalInclusiveSecondaryVertices
  process.svTable = svTable
  process.dummyValueMap = dummyValueMap
  process.sv_track = cms.Sequence(    process.unpackedTracksAndVertices*
                                      process.dummyValueMap*
                                      process.inclusiveVertexFinder*
                                      process.vertexMerger*
                                      process.trackVertexArbitrator*
                                      process.myFinalInclusiveSecondaryVertices*
                                      process.svTable)
  return process