// Truth-dependent labeling producer. This is the OTHER half of what used to be
// TrackGVFeaturesProducer: everything that needs gvTable / gvDaughtersTable /
// genPartIdx (generator-truth quantities, MC-only) lives here. It never needs to
// run in a real-data reconstruction sequence -- only in training/ntuple production,
// alongside TrackPairFeaturesProducer.
//
// To keep per-pair labels aligned with TrackPairFeaturesProducer's row order (rather
// than recomputing/reordering pairs independently), this producer reads that table's
// trk1Idx/trk2Idx columns directly and emits a same-sized, same-order "trackPair"
// FlatTable extension containing just edgeLabel.

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <limits>
#include <vector>
#include <iostream>

class TrackPairLabelProducer : public edm::stream::EDProducer<> {
public:
  explicit TrackPairLabelProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;
  edm::EDGetTokenT<nanoaod::FlatTable> gvTableToken_;
  edm::EDGetTokenT<nanoaod::FlatTable> gvDauTableToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> genPartIdxToken_;
  edm::EDGetTokenT<nanoaod::FlatTable> trackPairToken_;  // TrackPairFeaturesProducer's output

  double relPtMatchMax_;
  double dRMatchMax_;
  double trkPtCut_;
};

TrackPairLabelProducer::TrackPairLabelProducer(const edm::ParameterSet& iConfig)
    : tracksToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
      gvTableToken_(consumes<nanoaod::FlatTable>(iConfig.getParameter<edm::InputTag>("gvTable"))),
      gvDauTableToken_(consumes<nanoaod::FlatTable>(iConfig.getParameter<edm::InputTag>("gvDaughtersTable"))),
      genPartIdxToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("genPartIdx"))),
      trackPairToken_(consumes<nanoaod::FlatTable>(iConfig.getParameter<edm::InputTag>("trackPairSrc"))),
      relPtMatchMax_(iConfig.getParameter<double>("relPtMatchMax")),
      dRMatchMax_(iConfig.getParameter<double>("dRMatchMax")),
      trkPtCut_(iConfig.getParameter<double>("trkPtCut")) {
  produces<nanoaod::FlatTable>("track");      // extension of trackTable: hadidx, label
  produces<nanoaod::FlatTable>("trackPair");  // extension of TrackPairFeaturesProducer's trackPair table: edgeLabel
}

void TrackPairLabelProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("unpackedTracksAndVertices"));
  desc.add<edm::InputTag>("gvTable", edm::InputTag("genVertexProducer", "GVTable"));
  desc.add<edm::InputTag>("gvDaughtersTable", edm::InputTag("genVertexProducer", "GVDaughtersTable"));
  desc.add<edm::InputTag>("genPartIdx", edm::InputTag("trackGenMatch", "genPartIdx"));
  desc.add<edm::InputTag>("trackPairSrc", edm::InputTag("trackPairFeatures", "trackPair"));
  desc.add<double>("relPtMatchMax", 0.2);
  desc.add<double>("dRMatchMax", 0.02);
  desc.add<double>("trkPtCut", 0.4);
  descriptions.addWithDefaultLabel(desc);
}

void TrackPairLabelProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<std::vector<reco::Track>> tracksH;
  iEvent.getByToken(tracksToken_, tracksH);
  const auto& tracks = *tracksH;

  edm::Handle<edm::ValueMap<int>> genPartIdxH;
  iEvent.getByToken(genPartIdxToken_, genPartIdxH);

  edm::Handle<nanoaod::FlatTable> gvTable;
  iEvent.getByToken(gvTableToken_, gvTable);
  edm::Handle<nanoaod::FlatTable> gvDauTable;
  iEvent.getByToken(gvDauTableToken_, gvDauTable);

  edm::Handle<nanoaod::FlatTable> pairTable;
  iEvent.getByToken(trackPairToken_, pairTable);

  auto GV_isB = gvTable->columnData<int>(gvTable->columnIndex("isB"));
  auto GV_isD = gvTable->columnData<int>(gvTable->columnIndex("isD"));
  auto GV_fromHF = gvTable->columnData<int>(gvTable->columnIndex("fromHF"));

  size_t nDau = gvDauTable->size();
  auto Dau_pt = gvDauTable->columnData<float>(gvDauTable->columnIndex("pt"));
  auto Dau_eta = gvDauTable->columnData<float>(gvDauTable->columnIndex("eta"));
  auto Dau_phi = gvDauTable->columnData<float>(gvDauTable->columnIndex("phi"));
  auto Dau_hadIdx = gvDauTable->columnData<int>(gvDauTable->columnIndex("hadronIndex"));

  // Track - Daughters Matching (same logic as TrackGVFeaturesProducer had)
  size_t nTrk = tracks.size();
  std::vector<int> trk_hadidx(nTrk, -1);
  std::vector<int> trk_label(nTrk, 0);

  for (size_t it = 0; it < nTrk; ++it) {
    if (tracks[it].pt() < trkPtCut_) continue;
    float bestDR = std::numeric_limits<float>::infinity();
    int bestDau = -1;
    int genPartIdx = genPartIdxH->get(tracksH.id(), it);
    // 0=PU if no genMatch, 1=fromPV if genMatched but not from SV hadron
    trk_label[it] = (genPartIdx >= 0) ? 1 : 0;

    for (size_t id = 0; id < nDau; ++id) {
      if (Dau_pt[id] <= 0.f) continue;
      float ptRatio = std::abs(tracks[it].pt() - Dau_pt[id]) / Dau_pt[id];
      if (ptRatio > relPtMatchMax_) continue;
      float dr = reco::deltaR(tracks[it].eta(), tracks[it].phi(), Dau_eta[id], Dau_phi[id]);
      if (dr >= dRMatchMax_ || dr >= bestDR) continue;
      bestDR = dr;
      bestDau = static_cast<int>(id);
    }

    if (bestDau >= 0) {
      trk_hadidx[it] = Dau_hadIdx[bestDau];
      // 0: PU, 1: fromPV, 2: fromB, 3: fromB/C(HF), 4: fromC, 5: other secondary
      if (trk_hadidx[it] == -1) {
        trk_label[it] = 0;
        std::cout << "  -> label=0 (hadidx == -1, daughter has no hadron)" << std::endl;
      } else if (GV_isB[trk_hadidx[it]] == 1) {
        trk_label[it] = 2;
      } else if (GV_isD[trk_hadidx[it]] == 1 && GV_fromHF[trk_hadidx[it]] == 1) {
        trk_label[it] = 3;
      } else if (GV_isD[trk_hadidx[it]] == 1 && GV_fromHF[trk_hadidx[it]] == 0) {
        trk_label[it] = 4;
      } else {
        trk_label[it] = 5;
      }
    }
  }

  auto trackTableOut = std::make_unique<nanoaod::FlatTable>(nTrk, "track", false, true);  // extension=true
  trackTableOut->addColumn<int>("hadidx", trk_hadidx, "matched GV daughter's hadron index, -1 if unmatched");
  trackTableOut->addColumn<int>("label", trk_label, "truth label (currently mirrors flav; adjust as needed)");
  iEvent.put(std::move(trackTableOut), "track");

  // ---- Per-pair edgeLabel, aligned to TrackPairFeaturesProducer's row order ----
  size_t nPairs = pairTable->size();
  auto trk1Idx = pairTable->columnData<int>(pairTable->columnIndex("trk1Idx"));
  auto trk2Idx = pairTable->columnData<int>(pairTable->columnIndex("trk2Idx"));

  std::vector<float> edgeLabel(nPairs, 0.f);
  for (size_t k = 0; k < nPairs; ++k) {
    int i = trk1Idx[k];
    int j = trk2Idx[k];
    bool sameHadron = trk_hadidx[i] >= 0 && trk_hadidx[i] == trk_hadidx[j];
    edgeLabel[k] = sameHadron ? 1.f : 0.f;
  }

  auto pairLabelOut = std::make_unique<nanoaod::FlatTable>(nPairs, "trackPair", false, true);  // extension=true
  pairLabelOut->addColumn<float>("edgeLabel", edgeLabel, "1 if both tracks truth-matched to the same GV, else 0");
  iEvent.put(std::move(pairLabelOut), "trackPair");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TrackPairLabelProducer);