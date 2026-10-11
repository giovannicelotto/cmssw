#include <cmath>
#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/isFinite.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "RecoVertex/AdaptiveVertexFinder/interface/TTHelpers.h"

// Computes, per track in `tracks` (e.g. unpackedTracksAndVertices) w.r.t. the
// leading vertex in `primaryVertices`, the SAME two quantities that
// TemplatedInclusiveVertexFinder::produce() cuts on when building its `tts`:
//
//   dz       = track.dz(pv.position())                                         (== "maximumLongitudinalImpactParameter" quantity)
//   timeSig  = |track.timeExt() - pv.t()| / sqrt(track.dtErrorExt()^2 + pv.covariance(3,3))
//              only defined when timing info is available, exactly as in IVF   (== "maximumTimeSignificance" quantity)
//
// Also computes whether the track would pass TracksClusteringFromDisplacedSeed's SEED
// selection, i.e. the window cut in clusters():
//
//   ip = IPTools::absoluteImpactParameter3D(tt, pv);
//   ip.first && ip.value() in [seedMin3DIPValue, seedMax3DIPValue]
//            && ip.significance() in [seedMin3DIPSignificance, seedMax3DIPSignificance]
//
// Both the raw values and pass/fail flags (against configurable thresholds,
// defaulting to IVF's own defaults) are produced as edm::ValueMap<...>, so they
// plug into SimpleTrackFlatTableProducer's externalVariables/ExtVar exactly like
// genPartIdx/SVscore already do in your track table.
//
// NOTE: minHits/minPt do NOT need this producer -- they only depend on the track
// itself and can be written directly as `Var(...)` expressions in the flat table
class TrackVertexVars : public edm::stream::EDProducer<> {
public:
  explicit TrackVertexVars(const edm::ParameterSet &params);
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &event, const edm::EventSetup &es) override;

  edm::EDGetTokenT<std::vector<reco::Track>> token_tracks;
  edm::EDGetTokenT<reco::VertexCollection> token_primaryVertex;
  edm::EDGetTokenT<reco::BeamSpot> token_beamSpot;
  edm::EDGetTokenT<std::vector<pat::Muon>> token_muons;
  edm::EDGetTokenT<std::vector<pat::Electron>> token_electrons;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> token_trackBuilder;

  double maxLIP_;      // "maximumLongitudinalImpactParameter", same meaning/default as IVF
  double maxTimeSig_;  // "maximumTimeSignificance", same meaning/default as IVF

  // seed-selection window, same meaning/defaults as TracksClusteringFromDisplacedSeed's clusterizer PSet
  double seedMin3DIPValue_;
  double seedMax3DIPValue_;
  double seedMin3DIPSignificance_;
  double seedMax3DIPSignificance_;

  double dRMatchMax_;      // NEW: same style/meaning as trackGVFeaturesProducer's dRMatchMax_
  double relPtMatchMax_;   // NEW

  static constexpr float SENTINEL = -9999.f;  // used whenever a quantity can't be computed (no PV, no timing, invalid TT)
};

TrackVertexVars::TrackVertexVars(const edm::ParameterSet &params)
    : token_tracks(consumes<std::vector<reco::Track>>(params.getParameter<edm::InputTag>("tracks"))),
      token_primaryVertex(consumes<reco::VertexCollection>(params.getParameter<edm::InputTag>("primaryVertices"))),
      token_beamSpot(consumes<reco::BeamSpot>(params.getParameter<edm::InputTag>("beamSpot"))),
      token_muons(consumes<std::vector<pat::Muon>>(params.getParameter<edm::InputTag>("muons"))),           // NEW
      token_electrons(consumes<std::vector<pat::Electron>>(params.getParameter<edm::InputTag>("electrons"))), // NEW
      token_trackBuilder(
          esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
      maxLIP_(params.getParameter<double>("maximumLongitudinalImpactParameter")),
      maxTimeSig_(params.getParameter<double>("maximumTimeSignificance")),
      seedMin3DIPValue_(params.getParameter<double>("seedMin3DIPValue")),
      seedMax3DIPValue_(params.getParameter<double>("seedMax3DIPValue")),
      seedMin3DIPSignificance_(params.getParameter<double>("seedMin3DIPSignificance")),
      seedMax3DIPSignificance_(params.getParameter<double>("seedMax3DIPSignificance")),
      dRMatchMax_(params.getParameter<double>("dRMatchMax")),           // NEW
      relPtMatchMax_(params.getParameter<double>("relPtMatchMax")) {    // NEW
  produces<edm::ValueMap<float>>("dz");
  produces<edm::ValueMap<float>>("dzSignificance");
  produces<edm::ValueMap<float>>("ip2dValue");
  produces<edm::ValueMap<float>>("ip2dSignificance");
  produces<edm::ValueMap<float>>("timeSig");
  produces<edm::ValueMap<int>>("passLIP");
  produces<edm::ValueMap<int>>("passTimeSig");
  produces<edm::ValueMap<float>>("ip3dValue");
  produces<edm::ValueMap<float>>("ip3dSignificance");
  produces<edm::ValueMap<int>>("passSeed");
  produces<edm::ValueMap<int>>("isMuon");       // NEW
  produces<edm::ValueMap<int>>("isElectron");   // NEW
}

void TrackVertexVars::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("unpackedTracksAndVertices"));
  desc.add<edm::InputTag>("primaryVertices", edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpot"));
  desc.add<edm::InputTag>("muons", edm::InputTag("slimmedMuons"));           // NEW
  desc.add<edm::InputTag>("electrons", edm::InputTag("slimmedElectrons"));   // NEW
  desc.add<double>("maximumLongitudinalImpactParameter", 0.3);
  desc.add<double>("maximumTimeSignificance", 3.0);
  desc.add<double>("seedMin3DIPValue", 0.005);
  desc.add<double>("seedMax3DIPValue", 9999.0);
  desc.add<double>("seedMin3DIPSignificance", 1.2);
  desc.add<double>("seedMax3DIPSignificance", 9999.0);
  desc.add<double>("dRMatchMax", 0.01);      // NEW — tight by default, tracks are matched 1:1 to a lepton, not a jet-level cone
  desc.add<double>("relPtMatchMax", 0.1);    // NEW
  descriptions.addWithDefaultLabel(desc);
}

void TrackVertexVars::produce(edm::Event &event, const edm::EventSetup &es) {
  using namespace reco;

  edm::Handle<std::vector<Track>> tracks;
  event.getByToken(token_tracks, tracks);

  edm::Handle<VertexCollection> primaryVertices;
  event.getByToken(token_primaryVertex, primaryVertices);

  edm::Handle<BeamSpot> beamSpot;
  event.getByToken(token_beamSpot, beamSpot);


  edm::Handle<std::vector<pat::Muon>> muons;          // NEW
  event.getByToken(token_muons, muons);                // NEW
  edm::Handle<std::vector<pat::Electron>> electrons;    // NEW
  event.getByToken(token_electrons, electrons);         // NEW

  edm::ESHandle<TransientTrackBuilder> trackBuilder = es.getHandle(token_trackBuilder);

  size_t nTrk = tracks->size();
  std::vector<float> dz(nTrk, SENTINEL);
  std::vector<float> dzSig(nTrk, SENTINEL);
  std::vector<float> ip2dValue(nTrk, SENTINEL);
  std::vector<float> ip2dSignificance(nTrk, SENTINEL);
  std::vector<float> timeSig(nTrk, SENTINEL);
  std::vector<int> passLIP(nTrk, 0);
  std::vector<int> passTimeSig(nTrk, 1);  // IVF does not cut when timing is unavailable -> treat as "pass" by default, same as IVF's `continue` logic (no continue = pass)
  std::vector<float> ip3dValue(nTrk, SENTINEL);
  std::vector<float> ip3dSignificance(nTrk, SENTINEL);
  std::vector<int> passSeed(nTrk, 0);
  std::vector<int> isMuon(nTrk, 0);       // NEW
  std::vector<int> isElectron(nTrk, 0);   // NEW

  if (!primaryVertices->empty()) {
    const Vertex &pv = (*primaryVertices)[0];

    for (size_t i = 0; i < nTrk; ++i) {
      const Track &trk = (*tracks)[i];

      // dz: as computed in IVF 
      float trkDz = trk.dz(pv.position());
      if (!std::isfinite(trkDz)) {
        passLIP[i] = 0;
        continue;   // dz[i], dzSig[i] stay at SENTINEL
      }
      dz[i] = trkDz;
      passLIP[i] = (std::abs(trkDz) <= maxLIP_) ? 1 : 0;
      // dz significance
      float dzErr = std::sqrt(std::pow(trk.dzError(), 2) + pv.covariance(2, 2));
      if (dzErr > 0.f && std::isfinite(dzErr)) {
        dzSig[i] = trkDz / dzErr;
      }

      // needs a TransientTrack for timeExt()/dtErrorExt() and for the IP tools below,
      // built the same way IVF builds tts, then setBeamSpot the same way IVF does
      // before it hands tracks to the clusterizer's seed selection
      TransientTrack tt(tthelpers::buildTT(tracks, trackBuilder, i));
      if (!tt.isValid())
        continue;
      tt.setBeamSpot(*beamSpot);
      
      // IP2D
      std::pair<bool, Measurement1D> ip2d = IPTools::absoluteTransverseImpactParameter(tt, pv);
      if (ip2d.first) {
          float val = ip2d.second.value();
          float sig = ip2d.second.significance();
          if (edm::isFinite(val)) ip2dValue[i] = val;
          if (edm::isFinite(sig)) ip2dSignificance[i] = sig;
      }
      // IP3D
      std::pair<bool, Measurement1D> ip = IPTools::absoluteImpactParameter3D(tt, pv);
      if (ip.first) {
        float ipVal = ip.second.value();
        float ipSig = ip.second.significance();
        bool valOk = edm::isFinite(ipVal);
        bool sigOk = edm::isFinite(ipSig);
        if (valOk) ip3dValue[i] = ipVal;
        if (sigOk) ip3dSignificance[i] = ipSig;
        // only mark as passing/failing the seed window if BOTH quantities are well defined;
        // otherwise leave passSeed at its default (0), same as the ip.first==false case
        if (valOk && sigOk) {
          passSeed[i] = (ipVal >= seedMin3DIPValue_ && ipVal <= seedMax3DIPValue_ &&
                        ipSig >= seedMin3DIPSignificance_ && ipSig <= seedMax3DIPSignificance_)
                            ? 1
                            : 0;
          }
      }

      // timeSig
      if (edm::isFinite(tt.timeExt()) && pv.covariance(3, 3) > 0.) {
        auto tError = std::sqrt(std::pow(tt.dtErrorExt(), 2) + pv.covariance(3, 3));
        float sig = std::abs(tt.timeExt() - pv.t()) / tError;
        timeSig[i] = sig;
        passTimeSig[i] = (sig <= maxTimeSig_) ? 1 : 0;
      }
      // else: timeSig stays SENTINEL and passTimeSig stays 1 (no timing -> IVF does not reject on this cut)

      // else: ip3dValue/ip3dSignificance stay SENTINEL, passSeed stays 0 (same as IVF: ip.first==false -> not a seed)
    }
  }
  
  
  // ---- Track-lepton matching (isMuon / isElectron) ----
  // Same style as trackGVFeaturesProducer's daughter matching: best-dR match
  // within a relative-pT window, evaluated independently of PV availability.
  for (size_t i = 0; i < nTrk; ++i) {
    const Track &trk = (*tracks)[i];
    float trkEta = trk.eta();
    float trkPhi = trk.phi();
    float trkPt  = trk.pt();

    float bestDR_mu = std::numeric_limits<float>::infinity();
    for (const auto &mu : *muons) {
      float ptRatio = std::abs(trkPt - mu.pt()) / mu.pt();
      if (ptRatio > relPtMatchMax_) continue;
      float dr = reco::deltaR(trkEta, trkPhi, mu.eta(), mu.phi());
      if (dr >= dRMatchMax_ || dr >= bestDR_mu) continue;
      bestDR_mu = dr;
      isMuon[i] = 1;
    }
    //if (std::isfinite(bestDR_mu)) {
    //  std::cout << "[TrackVertexVars][muon] trk " << i
    //            << " closest dR=" << bestDR_mu
    //            << " matched=" << isMuon[i] << std::endl;
    //}

    float bestDR_e = std::numeric_limits<float>::infinity();
    for (const auto &ele : *electrons) {
      float ptRatio = std::abs(trkPt - ele.pt()) / ele.pt();
      if (ptRatio > relPtMatchMax_) continue;
      float dr = reco::deltaR(trkEta, trkPhi, ele.eta(), ele.phi());
      if (dr >= dRMatchMax_ || dr >= bestDR_e) continue;
      bestDR_e = dr;
      isElectron[i] = 1;
    }
    //if (std::isfinite(bestDR_e)) {
    //  std::cout << "[TrackVertexVars][electron] trk " << i
    //            << " closest dR=" << bestDR_e
    //            << " matched=" << isElectron[i] << std::endl;
    //}
  }
  auto dzMap = std::make_unique<edm::ValueMap<float>>();
  edm::ValueMap<float>::Filler dzFiller(*dzMap);
  dzFiller.insert(tracks, dz.begin(), dz.end());
  dzFiller.fill();
  event.put(std::move(dzMap), "dz");

  auto timeSigMap = std::make_unique<edm::ValueMap<float>>();
  edm::ValueMap<float>::Filler timeSigFiller(*timeSigMap);
  timeSigFiller.insert(tracks, timeSig.begin(), timeSig.end());
  timeSigFiller.fill();
  event.put(std::move(timeSigMap), "timeSig");

  auto passLIPMap = std::make_unique<edm::ValueMap<int>>();
  edm::ValueMap<int>::Filler passLIPFiller(*passLIPMap);
  passLIPFiller.insert(tracks, passLIP.begin(), passLIP.end());
  passLIPFiller.fill();
  event.put(std::move(passLIPMap), "passLIP");

  auto passTimeSigMap = std::make_unique<edm::ValueMap<int>>();
  edm::ValueMap<int>::Filler passTimeSigFiller(*passTimeSigMap);
  passTimeSigFiller.insert(tracks, passTimeSig.begin(), passTimeSig.end());
  passTimeSigFiller.fill();
  event.put(std::move(passTimeSigMap), "passTimeSig");

  auto ip3dValueMap = std::make_unique<edm::ValueMap<float>>();
  edm::ValueMap<float>::Filler ip3dValueFiller(*ip3dValueMap);
  ip3dValueFiller.insert(tracks, ip3dValue.begin(), ip3dValue.end());
  ip3dValueFiller.fill();
  event.put(std::move(ip3dValueMap), "ip3dValue");

  auto ip3dSigMap = std::make_unique<edm::ValueMap<float>>();
  edm::ValueMap<float>::Filler ip3dSigFiller(*ip3dSigMap);
  ip3dSigFiller.insert(tracks, ip3dSignificance.begin(), ip3dSignificance.end());
  ip3dSigFiller.fill();
  event.put(std::move(ip3dSigMap), "ip3dSignificance");

  auto passSeedMap = std::make_unique<edm::ValueMap<int>>();
  edm::ValueMap<int>::Filler passSeedFiller(*passSeedMap);
  passSeedFiller.insert(tracks, passSeed.begin(), passSeed.end());
  passSeedFiller.fill();
  event.put(std::move(passSeedMap), "passSeed");

  auto dzSigMap = std::make_unique<edm::ValueMap<float>>();
  edm::ValueMap<float>::Filler dzSigFiller(*dzSigMap);
  dzSigFiller.insert(tracks, dzSig.begin(), dzSig.end());
  dzSigFiller.fill();
  event.put(std::move(dzSigMap), "dzSignificance");

  auto ip2dValueMap = std::make_unique<edm::ValueMap<float>>();
  edm::ValueMap<float>::Filler ip2dValueFiller(*ip2dValueMap);
  ip2dValueFiller.insert(tracks, ip2dValue.begin(), ip2dValue.end());
  ip2dValueFiller.fill();
  event.put(std::move(ip2dValueMap), "ip2dValue");

  auto ip2dSigMap = std::make_unique<edm::ValueMap<float>>();
  edm::ValueMap<float>::Filler ip2dSigFiller(*ip2dSigMap);
  ip2dSigFiller.insert(tracks, ip2dSignificance.begin(), ip2dSignificance.end());
  ip2dSigFiller.fill();
  event.put(std::move(ip2dSigMap), "ip2dSignificance");

  auto isMuonMap = std::make_unique<edm::ValueMap<int>>();
  edm::ValueMap<int>::Filler isMuonFiller(*isMuonMap);
  isMuonFiller.insert(tracks, isMuon.begin(), isMuon.end());
  isMuonFiller.fill();
  event.put(std::move(isMuonMap), "isMuon");

  auto isElectronMap = std::make_unique<edm::ValueMap<int>>();
  edm::ValueMap<int>::Filler isElectronFiller(*isElectronMap);
  isElectronFiller.insert(tracks, isElectron.begin(), isElectron.end());
  isElectronFiller.fill();
  event.put(std::move(isElectronMap), "isElectron");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TrackVertexVars);