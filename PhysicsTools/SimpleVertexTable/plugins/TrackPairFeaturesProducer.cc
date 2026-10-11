// Physics-only track-pair (edge) feature producer.
//
// This is the truth-independent half of what used to be TrackGVFeaturesProducer's
// "trackPair" table: everything needed to compute dca/dcaSig/cptopv/pvtoPCA/dotprod/
// pairMom/pairInvMass/pt/eta/deltaR/deltaPhi/deltaEta, with NO dependency on gvTable/
// gvDaughtersTable/genPartIdx (those are generator-truth quantities that only exist
// in simulated MC events and would make this producer unusable on real data).
//
// The truth-dependent half (per-track hadidx/label, per-pair edgeLabel) lives in
// TrackPairLabelProducer, which consumes THIS producer's "trackPair" table (by
// trk1Idx/trk2Idx) to add an aligned edgeLabel column, and only needs to run in the
// training/ntuple-production sequence, never in general reconstruction.
//
// modelInference.cc reads this producer's "trackPair" table directly (by column name)
// instead of recomputing any of this inline, so there is exactly one implementation
// of the edge-feature math, shared between training-ntuple production and inference.

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include <cmath>
#include <algorithm>
#include <vector>

class TrackPairFeaturesProducer : public edm::stream::EDProducer<> {
public:
  explicit TrackPairFeaturesProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> pvToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;

  // Cut thresholds -- defaults match exactly what TrackGVFeaturesProducer used at
  // training time. Kept configurable rather than hardcoded so this stays the single
  // source of truth; if these ever change, everything downstream (labeling, training
  // ntuples, inference) picks it up automatically as long as they share one cfi.
  double maxDR_;
  double maxDEta_;
  double maxDPhi_;
  double maxPairEta_;
  double maxInvMass_;
  double maxDca_;

  static constexpr float PION_MASS = 0.13957018f;
};

TrackPairFeaturesProducer::TrackPairFeaturesProducer(const edm::ParameterSet& iConfig)
    : tracksToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
      pvToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
      ttbToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      maxDR_(iConfig.getParameter<double>("maxDR")),
      maxDEta_(iConfig.getParameter<double>("maxDEta")),
      maxDPhi_(iConfig.getParameter<double>("maxDPhi")),
      maxPairEta_(iConfig.getParameter<double>("maxPairEta")),
      maxInvMass_(iConfig.getParameter<double>("maxInvMass")),
      maxDca_(iConfig.getParameter<double>("maxDca")) {
  produces<nanoaod::FlatTable>("trackPair");
}

void TrackPairFeaturesProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("unpackedTracksAndVertices"));
  desc.add<edm::InputTag>("primaryVertices", edm::InputTag("unpackedTracksAndVertices"));
  desc.add<double>("maxDR", 0.8);
  desc.add<double>("maxDEta", 0.4);
  desc.add<double>("maxDPhi", 0.5);
  desc.add<double>("maxPairEta", 2.5);
  desc.add<double>("maxInvMass", 5.0);
  desc.add<double>("maxDca", 0.4);
  descriptions.addWithDefaultLabel(desc);
}

void TrackPairFeaturesProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<std::vector<reco::Track>> tracksH;
  iEvent.getByToken(tracksToken_, tracksH);
  const auto& tracks = *tracksH;

  edm::Handle<std::vector<reco::Vertex>> pvH;
  iEvent.getByToken(pvToken_, pvH);
  const reco::Vertex& pv = pvH->front();

  const auto& theB = iSetup.getData(ttbToken_);

  size_t nTrk = tracks.size();
  std::vector<reco::TransientTrack> t_trks(nTrk);
  for (size_t i = 0; i < nTrk; ++i) t_trks[i] = theB.build(tracks[i]);

  std::vector<int> trk1Idx, trk2Idx;
  std::vector<float> outDeltaR, outDeltaPhi, outDeltaEta, outDca, outDcaSig, outCptopv, outPvtoPCA1,
      outPvtoPCA2, outDotprod1, outDotprod2, outPairMom, outPairInvMass, outPt, outEta;

  VertexDistance3D vdist3D;

  // Full ordered pairs (i,j) with i != j -- NOT i<j. dotprod1/pvtoPCA1 treat i as
  // "seed" and dotprod2/pvtoPCA2 treat j as "track", so (i,j) and (j,i) are genuinely
  // different feature vectors, not a symmetric edge to deduplicate.
  for (size_t i = 0; i < nTrk; ++i) {
    if (!t_trks[i].isValid()) continue;
    for (size_t j = 0; j < nTrk; ++j) {
      if (i == j) continue;
      if (!t_trks[j].isValid()) continue;

      float dr = reco::deltaR(tracks[i].eta(), tracks[i].phi(), tracks[j].eta(), tracks[j].phi());
      float dEta = std::abs(tracks[i].eta() - tracks[j].eta());
      float dPhi = reco::deltaPhi(tracks[i].phi(), tracks[j].phi());
      if (dr > maxDR_) continue;
      if (dEta > maxDEta_) continue;
      if (std::abs(dPhi) > maxDPhi_) continue;

      float e1 = std::sqrt(tracks[i].p() * tracks[i].p() + PION_MASS * PION_MASS);
      float e2 = std::sqrt(tracks[j].p() * tracks[j].p() + PION_MASS * PION_MASS);
      float sumPx = tracks[i].px() + tracks[j].px();
      float sumPy = tracks[i].py() + tracks[j].py();
      float sumPz = tracks[i].pz() + tracks[j].pz();
      float pt_pair = std::sqrt(sumPx * sumPx + sumPy * sumPy);
      float compositeP = std::sqrt(sumPx * sumPx + sumPy * sumPy + sumPz * sumPz);
      float eta_pair;
      if (pt_pair < 1e-6f) {
        eta_pair = (sumPz > 0) ? 9999.f : -9999.f;
      } else {
        eta_pair = 0.5f * std::log((compositeP + sumPz) / (compositeP - sumPz));
      }
      if (std::abs(eta_pair) > maxPairEta_) continue;

      float sumE = e1 + e2;
      float invMass = std::sqrt(std::max(0.f, sumE * sumE - (sumPx * sumPx + sumPy * sumPy + sumPz * sumPz)));
      if (invMass > maxInvMass_) continue;

      float dca_val = -1.f, dcaSig_val = -1.f, cptopv_val = -1.f;
      float pvToPCA1_val = -1.f, pvToPCA2_val = -1.f, dot1_val = -999.f, dot2_val = -999.f, pairMom_val = -1.f;

      TwoTrackMinimumDistance minDist;
      if (minDist.calculate(t_trks[i].impactPointState(), t_trks[j].impactPointState())) {
        auto m = vdist3D.distance(
            VertexState(minDist.points().second, t_trks[i].impactPointState().cartesianError().position()),
            VertexState(minDist.points().first, t_trks[j].impactPointState().cartesianError().position()));
        dca_val = m.value();
        if (m.error() > 0) dcaSig_val = m.value() / m.error();

        GlobalPoint cp(minDist.crossingPoint());
        GlobalPoint pvp(pv.x(), pv.y(), pv.z());
        GlobalPoint seedPCA = minDist.points().second;
        GlobalPoint trackPCA = minDist.points().first;

        pvToPCA1_val = (seedPCA - pvp).mag();
        pvToPCA2_val = (trackPCA - pvp).mag();
        cptopv_val = (cp - pvp).mag();
        dot2_val = (trackPCA - pvp).unit().dot(t_trks[j].impactPointState().globalDirection().unit());
        dot1_val = (seedPCA - pvp).unit().dot(t_trks[i].impactPointState().globalDirection().unit());

        GlobalVector pairMomentum((Basic3DVector<float>)(t_trks[i].track().momentum() + t_trks[j].track().momentum()));
        pairMom_val = pairMomentum.mag();
      }

      if (dca_val > maxDca_) continue;

      trk1Idx.push_back(static_cast<int>(i));
      trk2Idx.push_back(static_cast<int>(j));
      outDeltaR.push_back(dr);
      outDeltaPhi.push_back(dPhi);
      outDeltaEta.push_back(dEta);
      outDca.push_back(dca_val);
      outDcaSig.push_back(dcaSig_val);
      outCptopv.push_back(cptopv_val);
      outPvtoPCA1.push_back(pvToPCA1_val);
      outPvtoPCA2.push_back(pvToPCA2_val);
      outDotprod1.push_back(dot1_val);
      outDotprod2.push_back(dot2_val);
      outPairMom.push_back(pairMom_val);
      outPairInvMass.push_back(invMass);
      outPt.push_back(pt_pair);
      outEta.push_back(eta_pair);
    }
  }

  auto pairTable = std::make_unique<nanoaod::FlatTable>(trk1Idx.size(), "trackPair", false);
  pairTable->addColumn<int>("trk1Idx", trk1Idx, "index of first track in the pair (into the source track collection)");
  pairTable->addColumn<int>("trk2Idx", trk2Idx, "index of second track in the pair");
  pairTable->addColumn<float>("deltaR", outDeltaR, "deltaR between the two tracks");
  pairTable->addColumn<float>("deltaPhi", outDeltaPhi, "deltaPhi between the two tracks");
  pairTable->addColumn<float>("deltaEta", outDeltaEta, "abs(eta1-eta2)");
  pairTable->addColumn<float>("dca", outDca, "distance of closest approach");
  pairTable->addColumn<float>("dcaSig", outDcaSig, "DCA significance");
  pairTable->addColumn<float>("cptopv", outCptopv, "distance from crossing point to PV");
  pairTable->addColumn<float>("pvtoPCA1", outPvtoPCA1, "distance PV to track1 PCA");
  pairTable->addColumn<float>("pvtoPCA2", outPvtoPCA2, "distance PV to track2 PCA");
  pairTable->addColumn<float>("dotprod1", outDotprod1, "direction dot product, track1 side");
  pairTable->addColumn<float>("dotprod2", outDotprod2, "direction dot product, track2 side");
  pairTable->addColumn<float>("pairMom", outPairMom, "magnitude of the pair momentum sum");
  pairTable->addColumn<float>("pairInvMass", outPairInvMass, "invariant mass, pion mass hypothesis");
  pairTable->addColumn<float>("pt", outPt, "transverse momentum of the pair");
  pairTable->addColumn<float>("eta", outEta, "eta of the pair");
  iEvent.put(std::move(pairTable), "trackPair");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TrackPairFeaturesProducer);