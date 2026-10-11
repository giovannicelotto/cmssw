#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
using namespace cms::Ort;
#include <iostream>
#include <cmath>
#include <limits>
#include <array>
#include <string>
#include "DataFormats/Common/interface/Wrapper.h"
#include <vector>


class modelInference : public edm::stream::EDProducer<edm::GlobalCache<ONNXRuntime>> {
public:
    explicit modelInference(const edm::ParameterSet&, const ONNXRuntime*);
    void produce(edm::Event &iEvent,const edm::EventSetup &iSetup) override;
    static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet&);
    static void globalEndJob(const ONNXRuntime*);

private:
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> pvsToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muonsToken_;
  edm::EDGetTokenT<std::vector<pat::Electron>> electronsToken_;
  edm::EDGetTokenT<nanoaod::FlatTable> trackPairToken_;  // TrackPairFeaturesProducer's output
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBToken;

  // Same style/meaning as TrackVertexVars' lepton-match window. IMPORTANT: the *training*
  // config (trackVertexVars_cfi in the python fragment) sets these to 0.001 / 0.2, which
  // differ from TrackVertexVars.cc's own fillDescriptions() defaults (0.01 / 0.1). Match
  // whatever was actually used in the ntuple production for these values.
  double dRMatchMax_;
  double relPtMatchMax_;

  static constexpr float PION_MASS = 0.13957018f;
  static constexpr int N_NODE_CLASSES = 6;  // sv_pred's output width (GraphEventGNN.sv_pred -> 6-way)
  // Class indices summed for the combined "is SV" score. Update if the class
  // ordering from config["truths"] turns out different from classes 2-5.
  static constexpr std::array<int, 4> SV_CLASS_INDICES = {{2, 3, 4, 5}};
};

modelInference::modelInference(const edm::ParameterSet& iConfig, const ONNXRuntime *cache):
  tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("src"))),
  pvsToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
  beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronsToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
  trackPairToken_(consumes<nanoaod::FlatTable>(iConfig.getParameter<edm::InputTag>("trackPairSrc"))),
  theTTBToken(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
  dRMatchMax_(iConfig.getParameter<double>("dRMatchMax")),
  relPtMatchMax_(iConfig.getParameter<double>("relPtMatchMax"))
{
    // One scalar ValueMap<float> per output class of sv_pred (6-way node classifier).
    // Placeholder names -- rename to match config["truths"] class order once settled;
    // index 0..5 here corresponds directly to sv_logits_flat[node*N_NODE_CLASSES + k].
    for (int k = 0; k < N_NODE_CLASSES; ++k) {
        produces<edm::ValueMap<float>>("nodeScore" + std::to_string(k));
    }
    // Combined "is SV" score: (score2+score3+score4+score5) / (score0+...+score5).
    // Denominator is included explicitly rather than assumed == 1, so this stays
    // correct even if N_NODE_CLASSES or which indices count as "SV" changes later.
    produces<edm::ValueMap<float>>("nodeScoreSVCombined");
    produces<edm::ValueMap<std::vector<float>>>("edgeScores");
    produces<edm::ValueMap<std::vector<int>>>("edgeIndices");
    // Extension of TrackPairFeaturesProducer's "trackPair" table: the model's raw
    // sigmoid score per pair, same row order/size, so it merges onto that table.
    produces<nanoaod::FlatTable>("trackPair");
}
std::unique_ptr<ONNXRuntime> modelInference::initializeGlobalCache(const edm::ParameterSet &iConfig) 
{
    return std::make_unique<ONNXRuntime>(iConfig.getParameter<edm::FileInPath>("model_path").fullPath());
}

void modelInference::globalEndJob(const ONNXRuntime *cache) {}

void modelInference::produce(edm::Event& iEvent,const edm::EventSetup &iSetup)
{

    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(tracksToken_, tracks);
    edm::Handle<std::vector<reco::Vertex>> pvs;
    iEvent.getByToken(pvsToken_, pvs);
    edm::Handle<reco::BeamSpot> beamSpot;
    iEvent.getByToken(beamSpotToken_, beamSpot);
    edm::Handle<std::vector<pat::Muon>> muons;
    iEvent.getByToken(muonsToken_, muons);
    edm::Handle<std::vector<pat::Electron>> electrons;
    iEvent.getByToken(electronsToken_, electrons);
    edm::Handle<nanoaod::FlatTable> pairTable;
    iEvent.getByToken(trackPairToken_, pairTable);


    const auto& PV0 = pvs->front();
    const auto& ttBuilder = iSetup.getData(theTTBToken);
    

    // retrieve information of tracks. Everything here mirrors TrackVertexVars.cc exactly
    // (same PV, same setBeamSpot call, same IPTools entry points, same dz-significance
    // formula) so that inference features match what the model was trained on.
    std::vector<float> trk_pt, trk_eta, trk_phi;
    std::vector<float> trk_dz, trk_dzSignificance;
    std::vector<float> trk_ip2dValue, trk_ip2dSignificance, trk_ip3dValue, trk_ip3dSignificance;
    std::vector<float> trk_p, trk_charge, trk_energy, trk_chi2;
    std::vector<float> trk_numberOfValidHits, trk_numberOfValidPixelHits, trk_numberOfValidStripHits;
    std::vector<float> trk_px, trk_py, trk_pz, trk_qdotp, trk_qoverp;
    std::vector<float> trk_isElectron, trk_isMuon;

    std::vector<float> track_SVscore;

    int num_tracks = tracks->size();
    track_SVscore.reserve(num_tracks);

    std::vector<int> origToNode(num_tracks, -1); // length of tracks. -1 for non valid ttrk. otherwise counts 0,1, 2, -1, 3, ..
    std::vector<int> nodeToOrig;
    nodeToOrig.reserve(num_tracks);
    int node = 0;

    // fill the vectors for each valid transient track
    for (size_t i = 0; i < tracks->size(); ++i) {
        const auto& trk = (*tracks)[i];
        reco::TransientTrack ttrk = ttBuilder.build(trk);

        if (!(ttrk.isValid())) continue;

        origToNode[i] = node;
 	    nodeToOrig.push_back((int)i);

        trk_pt.push_back(trk.pt());
        trk_eta.push_back(trk.eta());
        trk_phi.push_back(trk.phi());
        trk_p.push_back(trk.p());
        trk_charge.push_back(trk.charge());
        trk_numberOfValidHits.push_back(trk.numberOfValidHits());
        trk_numberOfValidPixelHits.push_back(trk.hitPattern().numberOfValidPixelHits());
        trk_numberOfValidStripHits.push_back(trk.hitPattern().numberOfValidStripHits());
        trk_energy.push_back(std::sqrt(trk.p() * trk.p() + PION_MASS * PION_MASS));
        trk_chi2.push_back(trk.chi2());
        trk_px.push_back(trk.px());
        trk_py.push_back(trk.py());
        trk_pz.push_back(trk.pz());
        trk_qdotp.push_back(trk.charge() * trk.p());
        trk_qoverp.push_back(trk.qoverp());

        // dz / dz significance -- identical formula to TrackVertexVars.cc
        float dzVal = trk.dz(PV0.position());
        float dzErr = std::sqrt(std::pow(trk.dzError(), 2) + PV0.covariance(2, 2));
        trk_dz.push_back(dzVal);
        trk_dzSignificance.push_back((dzErr > 0.f && std::isfinite(dzErr)) ? dzVal / dzErr : -9999.f);

        // IP2D / IP3D -- same IPTools entry points as TrackVertexVars.cc, same
        // setBeamSpot call beforehand (this was missing before and changes the result)
        ttrk.setBeamSpot(*beamSpot);

        float ip2dVal = -9999.f, ip2dSig = -9999.f;
        std::pair<bool, Measurement1D> ip2d = IPTools::absoluteTransverseImpactParameter(ttrk, PV0);
        if (ip2d.first) {
            if (edm::isFinite(ip2d.second.value())) ip2dVal = ip2d.second.value();
            if (edm::isFinite(ip2d.second.significance())) ip2dSig = ip2d.second.significance();
        }
        trk_ip2dValue.push_back(ip2dVal);
        trk_ip2dSignificance.push_back(ip2dSig);

        float ip3dVal = -9999.f, ip3dSig = -9999.f;
        std::pair<bool, Measurement1D> ip3d = IPTools::absoluteImpactParameter3D(ttrk, PV0);
        if (ip3d.first) {
            if (edm::isFinite(ip3d.second.value())) ip3dVal = ip3d.second.value();
            if (edm::isFinite(ip3d.second.significance())) ip3dSig = ip3d.second.significance();
        }
        trk_ip3dValue.push_back(ip3dVal);
        trk_ip3dSignificance.push_back(ip3dSig);

        node++;
    }

    // Track-lepton matching (isMuon / isElectron) -- same best-dR-within-relPt-window
    // logic as TrackVertexVars.cc, using the SAME dRMatchMax/relPtMatchMax values as
    // the training-time trackVertexVars configuration (pass 0.001 / 0.2 in the cfi,
    // not TrackVertexVars.cc's own fillDescriptions defaults of 0.01 / 0.1).
    trk_isMuon.assign(node, 0.f);
    trk_isElectron.assign(node, 0.f);
    for (int ni = 0; ni < node; ++ni) {
        int oi = nodeToOrig[ni];
        const auto& trk = (*tracks)[oi];
        float trkEta = trk.eta();
        float trkPhi = trk.phi();
        float trkPt = trk.pt();

        float bestDR_mu = std::numeric_limits<float>::infinity();
        for (const auto& mu : *muons) {
            float ptRatio = std::abs(trkPt - mu.pt()) / mu.pt();
            if (ptRatio > relPtMatchMax_) continue;
            float dr = reco::deltaR(trkEta, trkPhi, mu.eta(), mu.phi());
            if (dr >= dRMatchMax_ || dr >= bestDR_mu) continue;
            bestDR_mu = dr;
            trk_isMuon[ni] = 1.f;
        }

        float bestDR_e = std::numeric_limits<float>::infinity();
        for (const auto& ele : *electrons) {
            float ptRatio = std::abs(trkPt - ele.pt()) / ele.pt();
            if (ptRatio > relPtMatchMax_) continue;
            float dr = reco::deltaR(trkEta, trkPhi, ele.eta(), ele.phi());
            if (dr >= dRMatchMax_ || dr >= bestDR_e) continue;
            bestDR_e = dr;
            trk_isElectron[ni] = 1.f;
        }
    }

    // Edge (track-pair) features now come straight from TrackPairFeaturesProducer's
    // "trackPair" FlatTable -- read by column name, so this is immune to whatever
    // internal column order that producer happens to use. That producer is the
    // single shared implementation of the pair loop/cuts/geometry, used both for
    // training-ntuple production and here, so there is no risk of the two drifting
    // apart the way the old inline recomputation could.
    size_t nPairs = pairTable->size();
    auto tp_trk1Idx = pairTable->columnData<int>(pairTable->columnIndex("trk1Idx"));
    auto tp_trk2Idx = pairTable->columnData<int>(pairTable->columnIndex("trk2Idx"));
    auto tp_deltaEta = pairTable->columnData<float>(pairTable->columnIndex("deltaEta"));
    auto tp_deltaPhi = pairTable->columnData<float>(pairTable->columnIndex("deltaPhi"));
    auto tp_pt = pairTable->columnData<float>(pairTable->columnIndex("pt"));
    auto tp_eta = pairTable->columnData<float>(pairTable->columnIndex("eta"));
    auto tp_dca = pairTable->columnData<float>(pairTable->columnIndex("dca"));
    auto tp_deltaR = pairTable->columnData<float>(pairTable->columnIndex("deltaR"));
    auto tp_dcaSig = pairTable->columnData<float>(pairTable->columnIndex("dcaSig"));
    auto tp_cptopv = pairTable->columnData<float>(pairTable->columnIndex("cptopv"));
    auto tp_pvtoPCA1 = pairTable->columnData<float>(pairTable->columnIndex("pvtoPCA1"));
    auto tp_pvtoPCA2 = pairTable->columnData<float>(pairTable->columnIndex("pvtoPCA2"));
    auto tp_dotprod1 = pairTable->columnData<float>(pairTable->columnIndex("dotprod1"));
    auto tp_dotprod2 = pairTable->columnData<float>(pairTable->columnIndex("dotprod2"));
    auto tp_pairMom = pairTable->columnData<float>(pairTable->columnIndex("pairMom"));
    auto tp_pairInvMass = pairTable->columnData<float>(pairTable->columnIndex("pairInvMass"));

    std::vector<int> trk_idxs, trk_jdxs;
    std::vector<size_t> kept_row_idx;  // maps each kept entry back to its row in pairTable
    std::vector<float> pair_deltaEta, pair_deltaPhi, pair_pt, pair_eta;
    std::vector<float> deltaR, dca, dca_sig, cptopv, pvtoPCA_j, pvtoPCA_i, dotprod_j, dotprod_i, pair_mom, pair_invmass;

    for (size_t k = 0; k < nPairs; ++k) {
        int oi = tp_trk1Idx[k];
        int oj = tp_trk2Idx[k];
        // Defensive: TrackPairFeaturesProducer uses the same isValid() criterion to
        // decide which tracks can appear in a pair, so these should always resolve,
        // but skip cleanly rather than index with -1 if that ever changes upstream.
        if (origToNode[oi] < 0 || origToNode[oj] < 0) continue;

        trk_idxs.push_back(oi);
        trk_jdxs.push_back(oj);
        kept_row_idx.push_back(k);
        pair_deltaEta.push_back(tp_deltaEta[k]);
        pair_deltaPhi.push_back(tp_deltaPhi[k]);
        pair_pt.push_back(tp_pt[k]);
        pair_eta.push_back(tp_eta[k]);
        deltaR.push_back(tp_deltaR[k]);
        dca.push_back(tp_dca[k]);
        dca_sig.push_back(tp_dcaSig[k]);
        cptopv.push_back(tp_cptopv[k]);
        pvtoPCA_i.push_back(tp_pvtoPCA1[k]);
        pvtoPCA_j.push_back(tp_pvtoPCA2[k]);
        dotprod_i.push_back(tp_dotprod1[k]);
        dotprod_j.push_back(tp_dotprod2[k]);
        pair_mom.push_back(tp_pairMom[k]);
        pair_invmass.push_back(tp_pairInvMass[k]);
    }
    

    // === Assemble the 22 node (track) features, in the EXACT order the model was
    // trained on. Keep this list and its order in lockstep with graph_evt_0p8.yml's
    // cpf_candidates.
    std::vector<std::vector<float>> track_features;
    static constexpr int N_NODE_FEATURES = 22;
    static constexpr int N_EDGE_FEATURES = 14;

    for (int ni = 0; ni < node; ++ni) {
        std::vector<float> features = {
            trk_pt[ni],
            trk_eta[ni],
            trk_phi[ni],
            trk_dz[ni],
            trk_dzSignificance[ni],
            trk_ip2dSignificance[ni],
            trk_ip2dValue[ni],
            trk_ip3dSignificance[ni],
            trk_ip3dValue[ni],
            trk_numberOfValidHits[ni],
            trk_charge[ni],
            trk_energy[ni],
            trk_chi2[ni],
            trk_numberOfValidStripHits[ni],
            trk_numberOfValidPixelHits[ni],
            trk_px[ni],
            trk_py[ni],
            trk_pz[ni],
            trk_qdotp[ni],
            trk_qoverp[ni],
            trk_isElectron[ni],
            trk_isMuon[ni],
        };

        bool has_nan = false;
        for (float val : features) {
            if (!std::isfinite(val)) {
                has_nan = true;
                break;
            }
        }

        if (has_nan) {
            // Whole-track invalid sentinel. This is the -999.0 that GraphEventGNN.forward()
            // checks on feature[0] (track_pt) to mask the node with the model's learned
            // dummy token -- confirm this matches whatever DatasetConstructorTask pads
            // invalid/missing tracks with in cpf_arr, since it's a DIFFERENT convention
            // from the per-quantity -9999.f sentinel used above for individual IP/dz
            // failures on an otherwise-valid track.
            features.assign(N_NODE_FEATURES, -999.0f);
        }

        track_features.push_back(features);
    }

    // Build edge_index / per-edge feature rows from the (already column-name-matched)
    // pair data assembled above.
    std::vector<int64_t> edge_i, edge_j;
    std::vector<std::vector<float>> edge_features;
   	for (size_t idx = 0; idx < trk_idxs.size(); ++idx) {
   	    int oi = trk_idxs[idx];
   	    int oj = trk_jdxs[idx];

	    int ni = origToNode[oi];
	    int nj = origToNode[oj];

   	    edge_i.push_back(ni);
   	    edge_j.push_back(nj);
   	    
   	     edge_features.push_back({
   	     pair_deltaEta[idx],
   	     pair_deltaPhi[idx],
   	     pair_pt[idx],
   	     pair_eta[idx],
   	     dca[idx],
   	     deltaR[idx],
   	     dca_sig[idx],
   	     cptopv[idx],
   	     pvtoPCA_i[idx],
   	     pvtoPCA_j[idx],
   	     dotprod_i[idx],
   	     dotprod_j[idx],
   	     pair_mom[idx],
   	     pair_invmass[idx]
   	     }); 
   	 }

   	std::vector<float> x_in_flat;
   	x_in_flat.reserve(track_features.size() * N_NODE_FEATURES);
   	for (const auto& feat : track_features)
   	    x_in_flat.insert(x_in_flat.end(), feat.begin(), feat.end());

   	std::vector<float> edge_index_flat_f;
   	edge_index_flat_f.reserve(edge_i.size() * 2);
   	for (size_t k = 0; k < edge_i.size(); ++k) {
   	    edge_index_flat_f.push_back(static_cast<float>(edge_i[k]));
   	}

   	for (size_t k = 0; k < edge_j.size(); ++k) {
   	    edge_index_flat_f.push_back(static_cast<float>(edge_j[k]));
   	}

   	std::vector<float> edge_attr_flat;
   	edge_attr_flat.reserve(edge_features.size() * N_EDGE_FEATURES);
   	for (const auto& feat : edge_features){
   	    edge_attr_flat.insert(edge_attr_flat.end(), feat.begin(), feat.end());
   	    }   

   	std::vector<std::string> input_names_ = {"x_in", "edge_index", "edge_attr"};

   	std::vector<std::vector<int64_t>> input_shapes_ = {
   	 {1, static_cast<int64_t>(track_features.size()), N_NODE_FEATURES},  // x_in
   	 {1, 2, static_cast<int64_t>(edge_i.size())},                        // edge_index
   	 {1, static_cast<int64_t>(edge_features.size()), N_EDGE_FEATURES}    // edge_attr
   	};
   	   
   	std::vector<std::vector<float>> data_ = {
   	 x_in_flat,
   	 edge_index_flat_f,
   	 edge_attr_flat
   	};


   	std::vector<std::vector<float>> output = globalCache()->run(input_names_, data_, input_shapes_);

   	const std::vector<float>& sv_logits_flat = output[0];
   	const std::vector<float>& edge_logits_flat  = output[1]; // [E]

    // Full 6-way softmax over sv_logits_flat, which is [node * N_NODE_CLASSES + k]
    // row-major (matches GraphEventGNN.forward()'s sv_logits shape [N, 6]).
    auto softmax_node = [&](int ni) -> std::array<float, N_NODE_CLASSES> {
        std::array<float, N_NODE_CLASSES> probs{};
        float m = -std::numeric_limits<float>::infinity();
        for (int k = 0; k < N_NODE_CLASSES; ++k) {
            m = std::max(m, sv_logits_flat[ni * N_NODE_CLASSES + k]);
        }
        float denom = 0.f;
        for (int k = 0; k < N_NODE_CLASSES; ++k) {
            probs[k] = std::exp(sv_logits_flat[ni * N_NODE_CLASSES + k] - m);
            denom += probs[k];
        }
        for (int k = 0; k < N_NODE_CLASSES; ++k) {
            probs[k] /= denom;
        }
        return probs;
    };
    auto sigmoid = [](float x) { return 1.0f / (1.0f + std::exp(-x)); };
    std::vector<std::vector<float>> score_map(num_tracks);
    std::vector<std::vector<int>> idx_map(num_tracks);
    
    // Each entry in edge_i/edge_j is already a distinct directed edge (both (i,j) and
    // (j,i) were enumerated separately in the pair loop above, each with its own
    // edge_attr and therefore its own model-predicted score) -- do NOT symmetrize
    // here by writing the same score into both directions, that would overwrite the
    // genuinely different (j,i) prediction with (i,j)'s score.
    for (size_t e = 0; e < edge_i.size(); ++e) {
        int oi = nodeToOrig[edge_i[e]];
        int oj = nodeToOrig[edge_j[e]];
        float score = sigmoid(edge_logits_flat[e]);
        idx_map[oi].push_back(oj);
        score_map[oi].push_back(score);
    }
    
    auto vmScores = std::make_unique<edm::ValueMap<std::vector<float>>>();
    edm::ValueMap<std::vector<float>>::Filler fillerS(*vmScores);
    fillerS.insert(tracks, score_map.begin(), score_map.end());
    fillerS.fill();
    iEvent.put(std::move(vmScores), "edgeScores");
    
    auto vmIndices = std::make_unique<edm::ValueMap<std::vector<int>>>();
    edm::ValueMap<std::vector<int>>::Filler fillerI(*vmIndices);
    fillerI.insert(tracks, idx_map.begin(), idx_map.end());
    fillerI.fill();
    iEvent.put(std::move(vmIndices), "edgeIndices");

    // Same score, but shaped to extend TrackPairFeaturesProducer's "trackPair" table:
    // one row per pair, same order/size as that table, -1.f for any row this producer
    // didn't score (only possible via the defensive skip above).
    std::vector<float> pair_score(nPairs, -1.0f);
    for (size_t e = 0; e < edge_i.size(); ++e) {
        pair_score[kept_row_idx[e]] = sigmoid(edge_logits_flat[e]);
    }
    auto pairScoreTable = std::make_unique<nanoaod::FlatTable>(nPairs, "trackPair", false, true);  // extension=true
    pairScoreTable->addColumn<float>("score", pair_score, "model-predicted probability that this pair belongs to the same secondary vertex");
    iEvent.put(std::move(pairScoreTable), "trackPair");

    // One vector<float> per class, each indexed by original track index i, holding
    // that class's probability (or -1.0f sentinel for tracks with no valid node/TT).
    std::vector<std::vector<float>> selected_scores(N_NODE_CLASSES, std::vector<float>(tracks->size(), -1.0f));
    std::vector<float> selected_sv_combined(tracks->size(), -1.0f);

    for (size_t i = 0; i < tracks->size(); ++i) {
        int ni = origToNode[i];
        if (ni < 0) continue;  // stays -1.0f in every class, same convention as before

        auto probs = softmax_node(ni);
        for (int k = 0; k < N_NODE_CLASSES; ++k) {
            selected_scores[k][i] = probs[k];
        }

        float numer = 0.f, denom = 0.f;
        for (int k = 0; k < N_NODE_CLASSES; ++k) {
            denom += probs[k];
        }
        for (int idx : SV_CLASS_INDICES) {
            numer += probs[idx];
        }
        selected_sv_combined[i] = (denom > 0.f) ? (numer / denom) : -1.0f;
    }

    for (int k = 0; k < N_NODE_CLASSES; ++k) {
        auto valMap = std::make_unique<edm::ValueMap<float>>();
        edm::ValueMap<float>::Filler filler(*valMap);
        filler.insert(tracks, selected_scores[k].begin(), selected_scores[k].end());
        filler.fill();
        iEvent.put(std::move(valMap), "nodeScore" + std::to_string(k));
    }

    auto svCombinedMap = std::make_unique<edm::ValueMap<float>>();
    edm::ValueMap<float>::Filler svCombinedFiller(*svCombinedMap);
    svCombinedFiller.insert(tracks, selected_sv_combined.begin(), selected_sv_combined.end());
    svCombinedFiller.fill();
    iEvent.put(std::move(svCombinedMap), "nodeScoreSVCombined");
}

DEFINE_FWK_MODULE(modelInference);