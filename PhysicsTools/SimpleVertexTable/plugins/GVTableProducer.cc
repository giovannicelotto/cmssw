#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include <vector>
#include <unordered_set>
#include <limits>
#include <tuple>
#include <cmath>
#include <algorithm>

#include "Math/SMatrix.h"
#include "Math/SVector.h"
//#include "Math/SMatrixFunctions.h"
//typedef ROOT::Math::SMatrix<float, 3, 3> Matrix3x3;
//typedef ROOT::Math::SVector<float, 3> Vector3;
typedef ROOT::Math::SVector<double,3> Vector3D;
typedef reco::Vertex::CovarianceMatrix CovMatrix;
class GenVertexProducer : public edm::stream::EDProducer<> {
public:
    explicit GenVertexProducer(const edm::ParameterSet&);
    void produce(edm::Event&, const edm::EventSetup&) override;

private:

    int checkPDG(int abs_pdg) const;

    std::optional<std::tuple<float, float, float>>isAncestor(const reco::Candidate* mother,const reco::Candidate* daughter) const;
    bool hasHFAncestor(const reco::Candidate* hadron) const;
    bool hasHFDescendant(const reco::Candidate* hadron) const;

    std::vector<std::vector<float>> computeDistanceMatrix(
                    const std::vector<float>& SV_x,const std::vector<float>& SV_y,const std::vector<float>& SV_z,
                    std::vector<CovMatrix> SV_cov,
                    const std::vector<float>& Hadron_GVx,const std::vector<float>& Hadron_GVy,const std::vector<float>& Hadron_GVz);
    void printDistanceMatrix(const std::vector<std::vector<float>>& distances);
    std::tuple<std::vector<int>, std::vector<float>, std::vector<float>, std::vector<int>, std::vector<int>, std::vector<int>>  matchHadronsToSV(
                                                                        std::vector<std::vector<float>> distances,
                                                                        const std::vector<float>& SVtrk_pt, const std::vector<float>& SVtrk_eta, const std::vector<float>& SVtrk_phi, const std::vector<int>& SVtrk_SVidx,
                                                                        const std::vector<float>& Daughters_pt,     //genparticles
                                                                        const std::vector<float>& Daughters_eta,  //genparticles
                                                                        const std::vector<float>& Daughters_phi,  //genparticles
                                                                        const std::vector<int>& Daughters_GVidx, // hadron index per daughter
                                                                        const std::vector<float>& SV_eta,
                                                                        const std::vector<float>& SV_phi,
                                                                        const std::vector<float>& GV_eta,
                                                                        const std::vector<float>& GV_phi,
                                                                        int n_Hadrons,
                                                                        int nRequiredCommonTracks,
                                                                        double dR_max,
                                                                        double relPt_max,
                                                                        bool doubleMatching,
                                                                        int doubleMatching_nRequiredCommonTracks,
                                                                        double doubleMatching_maxSignificance,
                                                                        double doubleMatching_dR_max,
                                                                        double doubleMatching_relPt_max
                                                                    );

    // NEW: matches GVDaughters (gen daughters) to reconstructed tracks.
    // Returns, indexed the same way as the Daughters_* vectors:
    //   trkIdx      : index into the `tracks` collection of the matched track, -1 if unmatched
    //   isMatched   : 1 if matched, 0 otherwise
    //   matchDeltaR : deltaR to the matched track, -1 if unmatched
    //   matchDPtRel : |pt_trk - pt_dau| / pt_dau of the matched track, -1 if unmatched
    std::tuple<std::vector<int>, std::vector<int>, std::vector<float>, std::vector<float>> matchDaughtersToTracks(
                                                                        const std::vector<float>& Daughters_pt,
                                                                        const std::vector<float>& Daughters_eta,
                                                                        const std::vector<float>& Daughters_phi,
                                                                        const std::vector<int>& Daughters_charge,
                                                                        const std::vector<reco::Track>& tracks,
                                                                        double maxDeltaR,
                                                                        double maxDPtRel,
                                                                        bool checkCharge,
                                                                        bool resolveAmbiguities) const;

    const edm::EDGetTokenT<std::vector<reco::Vertex>> pvs_;
    edm::EDGetTokenT<edm::View<reco::Candidate>> genToken_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> svToken_;
    edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;  // NEW: general tracks used for GVDaughters matching
    int nRequiredCommonTracks_;
    double dlenSigMin_;
    double dR_max_;
    double relPt_max_;
    bool doubleMatching_;
    int doubleMatching_nRequiredCommonTracks_;
    double doubleMatching_maxSignificance_;
    double doubleMatching_dR_max_;
    double doubleMatching_relPt_max_;

    // NEW: GVDaughters <-> tracks matching configuration
    double trkMaxDeltaR_;
    double trkMaxDPtRel_;
    bool trkCheckCharge_;
    bool trkResolveAmbiguities_;
};


GenVertexProducer::GenVertexProducer(const edm::ParameterSet& iConfig):
    pvs_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
    genToken_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
    svToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("secondaryVertices"))),
    tracksToken_(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"))),
    nRequiredCommonTracks_(iConfig.getParameter<int>("nRequiredCommonTracks")),
    dlenSigMin_(iConfig.getParameter<double>("dlenSigMin")),
    dR_max_(iConfig.getParameter<double>("dR_max")),
    relPt_max_(iConfig.getParameter<double>("relPt_max")),
    doubleMatching_(iConfig.getParameter<bool>("doubleMatching")),
    doubleMatching_nRequiredCommonTracks_(iConfig.getParameter<int>("doubleMatching_nRequiredCommonTracks")),
    doubleMatching_maxSignificance_(iConfig.getParameter<double>("doubleMatching_maxSignificance")),
    doubleMatching_dR_max_(iConfig.getParameter<double>("doubleMatching_dR_max")),
    doubleMatching_relPt_max_(iConfig.getParameter<double>("doubleMatching_relPt_max")),
    trkMaxDeltaR_(iConfig.getParameter<double>("trkMaxDeltaR")),
    trkMaxDPtRel_(iConfig.getParameter<double>("trkMaxDPtRel")),
    trkCheckCharge_(iConfig.getParameter<bool>("trkCheckCharge")),
    trkResolveAmbiguities_(iConfig.getParameter<bool>("trkResolveAmbiguities"))
{
    produces<nanoaod::FlatTable>("GVTable");
    produces<nanoaod::FlatTable>("rejectedGVTable");
    produces<nanoaod::FlatTable>("GVDaughtersTable");
    produces<nanoaod::FlatTable>("GVDirectDaughters");
    produces<nanoaod::FlatTable>("SVGVMatchTable");
    produces<nanoaod::FlatTable>("SVGVtrkMatchTable"); 

}


void GenVertexProducer::produce(edm::Event& iEvent,
             const edm::EventSetup&) 
    {
        edm::Handle<edm::View<reco::Candidate>> genHandle;
        iEvent.getByToken(genToken_, genHandle);
        edm::Handle<std::vector<reco::Vertex>> svHandle;
        iEvent.getByToken(svToken_, svHandle);
        auto pvsIn = iEvent.getHandle(pvs_);
        edm::Handle<std::vector<reco::Track>> tracksHandle;   // NEW
        iEvent.getByToken(tracksToken_, tracksHandle);        // NEW

        



        const auto& genParticles = genHandle;
        const auto& secondaryVertices = svHandle;
        const auto& tracks = *tracksHandle;                   // NEW

        // Output vectors
        std::vector<float> Hadron_pt, Hadron_eta, Hadron_phi;
        std::vector<float> SV_x, SV_y, SV_z, SV_eta, SV_phi;
        std::vector<CovMatrix> SV_cov;
        std::vector<float> Hadron_GVx, Hadron_GVy, Hadron_GVz;
        std::vector<float>  Hadron_GVx_i, Hadron_GVy_i, Hadron_GVz_i;
        std::vector<int> Hadron_pdgId;
        std::vector<int> Hadron_pdgClass, Hadron_isB, Hadron_isD;
        std::vector<int> Hadron_fromHF, Hadron_toHF;
        std::vector<int> GV_nDaughters; // NEW: total number of GVDaughters per GV (denominator for nDaughtersMatchedToTracks)
        std::vector<float> GV_maxDaughterPairDeltaR; // NEW: max pairwise deltaR among a GV's own GVDaughters (no tracks involved)
        std::vector<float> Daughters_pt, Daughters_eta, Daughters_phi;
        std::vector<float> Daughters_vx, Daughters_vy, Daughters_vz;
        std::vector<int> Daughters_charge, Daughters_GVidx, Daughters_pdgId;
        
        std::vector<float> allHadron_GVx, allHadron_GVy, allHadron_GVz;
        std::vector<float>  allHadron_GVx_i, allHadron_GVy_i, allHadron_GVz_i;
        std::vector<int> allHadron_pdgId;
        std::vector<float> allHadron_pt, allHadron_eta, allHadron_phi;

        std::vector<float> directDaughters_pt, directDaughters_eta, directDaughters_phi;
        std::vector<int> directDaughters_charge, directDaughters_GVidx, directDaughters_pdgId;
        VertexDistance3D vdist;
        const auto& PV0 = pvsIn->front();

        // save coordinates of SV (will be used for matching with GV)
        for (auto const& sv : *secondaryVertices) {
            Measurement1D dl = vdist.distance(PV0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
            if (dl.value() > 0 and dl.significance() > dlenSigMin_) {
                SV_x.push_back(sv.x());
                SV_y.push_back(sv.y());
                SV_z.push_back(sv.z());

                // Get Eta and Phi from tracks
                TLorentzVector p4s_SV = TLorentzVector(0,0,0,0);
                for (auto it = sv.tracks_begin(); it != sv.tracks_end(); ++it) {
                    const edm::RefToBase<reco::Track>& trkRef = *it;
                    TLorentzVector p4;
                    p4.SetPtEtaPhiM(trkRef->pt(),trkRef->eta(),trkRef->phi(),0.13957039);
                    p4s_SV += p4;
                }
                SV_eta.push_back(p4s_SV.Eta());
                SV_phi.push_back(p4s_SV.Phi());
                SV_cov.push_back(sv.covariance());
            }
        }

        // Filling Hadrons and Daughters vectors
        int ngv=0;
        int nRejectedGV=0;
        int ngv_b=0, ngv_d=0, ngv_s=0, ngv_tau=0;
        for(size_t i=0; i<genParticles->size(); ++i){
            const reco::Candidate* hadron = &(*genParticles)[i];
            //std::cout<<"Hadron "<<i<<" PDG ID: "<<hadron->pdgId()<<", pt: "<<hadron->pt()<<", eta: "<<hadron->eta()<<std::endl;
            if(!(hadron->pt()>10 && std::abs(hadron->eta())<2.5)) continue;

            int hadPDG = checkPDG(std::abs(hadron->pdgId())); // 1: Beauty, 2: Charmed, 3: Strange,  4: Tau,  0: Else
            if(hadPDG==0) continue;
            //     code here
            //    
            //    
            //    
            //    
            //  

                


            //  Collect stable charged daughters
            std::vector<float> temp_pt, temp_eta, temp_phi, temp_vx, temp_vy, temp_vz; // kinematics of gen daughters of the hadron in the loop
            std::vector<int> temp_charge, temp_GVidx, temp_flav, temp_pdgId;
            int nPack=0;
            float vx=std::numeric_limits<float>::quiet_NaN();
            float vy=std::numeric_limits<float>::quiet_NaN();
            float vz=std::numeric_limits<float>::quiet_NaN();

            for(size_t j=0; j<genParticles->size(); ++j){
                const reco::Candidate* dau = &(*genParticles)[j];
                if(dau==hadron) continue;
                if(!(dau->status()==1 && dau->charge()!=0 && dau->pt()>0.8 && std::abs(dau->eta())<2.5)) continue;

                auto GV = isAncestor(hadron,dau); //takes the x,y,z of the daughter (decay point of the hadron) if daughters otherwise return nan
                if(GV.has_value()){
                    std::tie(vx,vy,vz) = *GV;
                    if(!std::isnan(vx)){
                        nPack++;
                        temp_pt.push_back(dau->pt());
                        temp_eta.push_back(dau->eta());
                        temp_phi.push_back(dau->phi());
                        temp_vx.push_back(dau->vx());
                        temp_vy.push_back(dau->vy());
                        temp_vz.push_back(dau->vz());
                        temp_charge.push_back(dau->charge());
                        temp_pdgId.push_back(dau->pdgId());
                        temp_GVidx.push_back(ngv); // hadron index
                        //temp_flav.push_back(hadPDG);
                    }
                }
            }
            // If has more than 2 good daughters, the Hadron is Good, we found a GV:
            if(nPack>=2){
                // Save hadron info
                //std::cout<<"Found hadron "<<ngv<<" PDG ID: "<<hadron->pdgId()<<", pt: "<<hadron->pt()<<", eta: "<<hadron->eta()<<std::endl;
                Hadron_pt.push_back(hadron->pt());
                Hadron_eta.push_back(hadron->eta());
                Hadron_phi.push_back(hadron->phi());
                Hadron_pdgId.push_back(hadron->pdgId());
                Hadron_pdgClass.push_back(hadPDG);
                
                

                // Save GenVertex
                ngv++;
                if(hadPDG==1) {
                    ngv_b++;
                    Hadron_isB.push_back(1);
                    }
                else{
                    Hadron_isB.push_back(0);
                }
                if(hadPDG==2) {
                    ngv_d++;
                    Hadron_isD.push_back(1);
                    }
                else{
                    Hadron_isD.push_back(0);
                }
                if(hadPDG==3) ngv_s++;
                if(hadPDG==4) ngv_tau++;
                Hadron_fromHF.push_back(hasHFAncestor(hadron) ? 1 : 0);
                Hadron_toHF.push_back(hasHFDescendant(hadron) ? 1 : 0);
                GV_nDaughters.push_back(nPack);  // NEW: total number of GVDaughters belonging to this GV

                // NEW: max pairwise deltaR among this GV's own daughters (gen-level, no tracks involved)
                float maxPairDR = 0.f;
                for (size_t a = 0; a < temp_eta.size(); ++a) {
                    for (size_t b = a + 1; b < temp_eta.size(); ++b) {
                        float dR = deltaR(temp_eta[a], temp_phi[a], temp_eta[b], temp_phi[b]);
                        if (dR > maxPairDR) maxPairDR = dR;
                    }
                }
                GV_maxDaughterPairDeltaR.push_back(maxPairDR);
                
                Hadron_GVx.push_back(vx);               // point of decay of the hadron
                Hadron_GVy.push_back(vy);               // point of decay of the hadron
                Hadron_GVz.push_back(vz);               // point of decay of the hadron
                Hadron_GVx_i.push_back(hadron->vx());   // point of origin of the hadron
                Hadron_GVy_i.push_back(hadron->vy());   // point of origin of the hadron
                Hadron_GVz_i.push_back(hadron->vz());   // point of origin of the hadron

                // Save daughters
                Daughters_pt.insert(Daughters_pt.end(), temp_pt.begin(), temp_pt.end());
                Daughters_vx.insert(Daughters_vx.end(), temp_vx.begin(), temp_vx.end());
                Daughters_vy.insert(Daughters_vy.end(), temp_vy.begin(), temp_vy.end());
                Daughters_vz.insert(Daughters_vz.end(), temp_vz.begin(), temp_vz.end());
                Daughters_eta.insert(Daughters_eta.end(), temp_eta.begin(), temp_eta.end());
                Daughters_phi.insert(Daughters_phi.end(), temp_phi.begin(), temp_phi.end());
                Daughters_charge.insert(Daughters_charge.end(), temp_charge.begin(), temp_charge.end());
                Daughters_pdgId.insert(Daughters_pdgId.end(), temp_pdgId.begin(), temp_pdgId.end());
                Daughters_GVidx.insert(Daughters_GVidx.end(), temp_GVidx.begin(), temp_GVidx.end());
                for(size_t j=0; j<genParticles->size(); ++j){   
                    const reco::Candidate* dau = &(*genParticles)[j];
                    if(dau==hadron) continue;
                    if (dau->numberOfMothers() > 0){
                        const reco::Candidate* mother = dau->mother(0);
                        if (mother == hadron){
                            directDaughters_pt.push_back(dau->pt());
                            directDaughters_eta.push_back(dau->eta());
                            directDaughters_phi.push_back(dau->phi());
                            directDaughters_charge.push_back(dau->charge());
                            directDaughters_pdgId.push_back(dau->pdgId());
                            directDaughters_GVidx.push_back(ngv-1);
                        }
                    }
                }
            }
            else{
                //fallback to save the decay point of the hadron if it has no daughters, but is still a rejected GV
                if (std::isnan(vx) && hadron->numberOfDaughters() > 0) {
                    const reco::Candidate* directDau = hadron->daughter(0);
                    vx = directDau->vx();
                    vy = directDau->vy();
                    vz = directDau->vz();
                }
                nRejectedGV++;
                allHadron_pt.push_back(hadron->pt());
                allHadron_eta.push_back(hadron->eta());
                allHadron_phi.push_back(hadron->phi());
                allHadron_pdgId.push_back(hadron->pdgId());
                allHadron_GVx.push_back(vx);               // point of decay of the hadron
                allHadron_GVy.push_back(vy);               // point of decay of the hadron
                allHadron_GVz.push_back(vz);               // point of decay of the hadron
                allHadron_GVx_i.push_back(hadron->vx());   // point of origin of the hadron
                allHadron_GVy_i.push_back(hadron->vy());   // point of origin of the hadron
                allHadron_GVz_i.push_back(hadron->vz());   // point of origin of the hadron
            }
        }


        // Filling tracks from reco SV
        std::vector<float> SVtrk_pt, SVtrk_eta, SVtrk_phi;
        std::vector<int> SVtrk_SVidx;
        int SV_index=0;
        for (const auto &sv : *secondaryVertices) {
            Measurement1D dl = vdist.distance(PV0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
            if (dl.value() > 0 and dl.significance() > dlenSigMin_) {
            for (auto it = sv.tracks_begin(); it != sv.tracks_end(); ++it) {
                const edm::RefToBase<reco::Track>& trkRef = *it;
                if (trkRef.isNull()) continue;
                TLorentzVector p4;
                p4.SetPtEtaPhiM(trkRef->pt(),trkRef->eta(),trkRef->phi(),0.13957039);
                SVtrk_pt.push_back(trkRef->pt());
                SVtrk_eta.push_back(trkRef->eta());
                SVtrk_phi.push_back(trkRef->phi());
                SVtrk_SVidx.push_back(SV_index); 
            }
            SV_index++;
            }
        }

        
        // Compute matrix of distances between SV and GV
        auto distances = computeDistanceMatrix(SV_x, SV_y, SV_z, SV_cov,Hadron_GVx, Hadron_GVy, Hadron_GVz);
        //printDistanceMatrix(distances);
        
        std::vector<int> Hadron_SVIdx(ngv, -1); // 
        std::vector<float> Hadron_SVDistance(ngv, -1); // 
        std::vector<float> Hadron_minDistNotMatched(ngv, 999.f);

        // perform matching based on distance matrix and track-to-daughter matching
        auto result = matchHadronsToSV(distances,SVtrk_pt, SVtrk_eta, SVtrk_phi, SVtrk_SVidx,Daughters_pt, Daughters_eta, Daughters_phi, Daughters_GVidx,
                                        SV_eta,SV_phi,Hadron_eta,Hadron_phi,ngv,nRequiredCommonTracks_,dR_max_,relPt_max_,doubleMatching_ , doubleMatching_nRequiredCommonTracks_, doubleMatching_maxSignificance_ , doubleMatching_dR_max_, doubleMatching_relPt_max_  );
        
        Hadron_SVIdx             = std::get<0>(result);
        Hadron_SVDistance        = std::get<1>(result);
        Hadron_minDistNotMatched = std::get<2>(result);
        std::vector<int> SVtrk_isMatched = std::get<3>(result);
        std::vector<int> SVtrk_GVIdx     = std::get<4>(result);
        std::vector<int> SVtrk_daughterIdx = std::get<5>(result);  // NEW
        auto svTrkGVTable = std::make_unique<nanoaod::FlatTable>(SVtrk_pt.size(), "mySVtrks", false, true);
        svTrkGVTable->addColumn<int>("isMatched", SVtrk_isMatched, "1 if track is matched to a genParticle daughter of the matched GV");
        svTrkGVTable->addColumn<int>("GVIdx", SVtrk_GVIdx, "Index of matched GenVertex hadron, -1 if unmatched");
        svTrkGVTable->addColumn<int>("daughterIdx", SVtrk_daughterIdx, "Index into GVDaughters table of the matched gen daughter, -1 if unmatched"); // NEW

        iEvent.put(std::move(svTrkGVTable), "SVGVtrkMatchTable");

        // Build per-SV matched flag / GV index (same ordering as SV_x, since same dlenSig cut)
        std::vector<int> SV_GVIdx(SV_x.size(), -1);
        std::vector<int> SV_isMatched(SV_x.size(), 0);
        for (int had = 0; had < ngv; ++had) {
            int sv = Hadron_SVIdx[had];
            if (sv >= 0) {
                SV_GVIdx[sv] = had;
                SV_isMatched[sv] = 1;
            }
        }

        // extension=true, name must match SVTableProducer's table name ("mySV")
        auto svGVTable = std::make_unique<nanoaod::FlatTable>(SV_x.size(), "mySV", false, true);
        svGVTable->addColumn<int>("GVIdx", SV_GVIdx, "Index of matched GenVertex hadron, -1 if unmatched");
        svGVTable->addColumn<int>("isMatched", SV_isMatched, "1 if SV matched to a GV");
        iEvent.put(std::move(svGVTable), "SVGVMatchTable");
        //std::cout<<Hadron_minDistNotMatched<<" is the min distance of unmatched hadrons"<<std::endl;

        // NEW: match GVDaughters (gen daughters) to reconstructed tracks
        auto trkMatchResult = matchDaughtersToTracks(
            Daughters_pt, Daughters_eta, Daughters_phi, Daughters_charge, tracks,
            trkMaxDeltaR_, trkMaxDPtRel_, trkCheckCharge_, trkResolveAmbiguities_);

        std::vector<int>   Daughters_trkIdx       = std::get<0>(trkMatchResult);
        std::vector<int>   Daughters_isTrkMatched = std::get<1>(trkMatchResult);
        std::vector<float> Daughters_trkDeltaR     = std::get<2>(trkMatchResult);
        std::vector<float> Daughters_trkDPtRel     = std::get<3>(trkMatchResult);

        // NEW: per-GV count of daughters matched to a track
        std::vector<int> GV_nDaughtersMatchedToTracks(ngv, 0);
        for (size_t i = 0; i < Daughters_GVidx.size(); ++i) {
            if (!Daughters_isTrkMatched[i]) continue;
            int gvIdx = Daughters_GVidx[i];
            if (gvIdx >= 0 && gvIdx < ngv) {
                GV_nDaughtersMatchedToTracks[gvIdx]++;
            }
        }

        // NEW: per-GV max trkDeltaR among its (matched) daughters (GVDaughters, not GVDirectDaughters).
        // Stays at -1 for a GV that has no track-matched daughters.
        std::vector<float> GV_maxDaughterTrkDeltaR(ngv, -1.f);
        for (size_t i = 0; i < Daughters_GVidx.size(); ++i) {
            if (!Daughters_isTrkMatched[i]) continue;
            int gvIdx = Daughters_GVidx[i];
            if (gvIdx >= 0 && gvIdx < ngv) {
                if (Daughters_trkDeltaR[i] > GV_maxDaughterTrkDeltaR[gvIdx]) {
                    GV_maxDaughterTrkDeltaR[gvIdx] = Daughters_trkDeltaR[i];
                }
            }
        }

        //  Build FlatTables 
        auto rejectedGVTable = std::make_unique<nanoaod::FlatTable>(nRejectedGV,"RejectedGV",false);
        rejectedGVTable->addColumn<float>("pt",allHadron_pt,"Rejected Hadron pt");
        rejectedGVTable->addColumn<float>("eta",allHadron_eta,"Rejected Hadron eta");
        rejectedGVTable->addColumn<float>("phi",allHadron_phi,"Rejected Hadron phi");
        rejectedGVTable->addColumn<float>("x",allHadron_GVx,"Rejected GV x");
        rejectedGVTable->addColumn<float>("y",allHadron_GVy,"Rejected GV y");
        rejectedGVTable->addColumn<float>("z",allHadron_GVz,"Rejected GV z");
        rejectedGVTable->addColumn<float>("x_i",allHadron_GVx_i,"Born x coordinate of Rejected GV ");
        rejectedGVTable->addColumn<float>("y_i",allHadron_GVy_i,"Born y coordinate of Rejected GV ");
        rejectedGVTable->addColumn<float>("z_i",allHadron_GVz_i,"Born z coordinate of Rejected GV ");

        auto gvTable = std::make_unique<nanoaod::FlatTable>(ngv,"GV",false);
        gvTable->addColumn<float>("pt",Hadron_pt,"Hadron pt");
        gvTable->addColumn<float>("eta",Hadron_eta,"Hadron eta");
        gvTable->addColumn<float>("phi",Hadron_phi,"Hadron phi");
        gvTable->addColumn<float>("x",Hadron_GVx,"GV x");
        gvTable->addColumn<float>("y",Hadron_GVy,"GV y");
        gvTable->addColumn<float>("z",Hadron_GVz,"GV z");
        gvTable->addColumn<float>("x_i",Hadron_GVx_i,"Born x coordinate of GV ");
        gvTable->addColumn<float>("y_i",Hadron_GVy_i,"Born y coordinate of GV ");
        gvTable->addColumn<float>("z_i",Hadron_GVz_i,"Born z coordinate of GV ");
        gvTable->addColumn<int>("Hadron_SVIdx",Hadron_SVIdx,"SVIdx");
        gvTable->addColumn<int>("Hadron_pdgId",Hadron_pdgId,"Hadron_pdgId");
        gvTable->addColumn<float>("SV_distanceSig",Hadron_SVDistance,"SV_distanceSig");
        // new class
        gvTable->addColumn<int>("isB",Hadron_isB,"isB");
        gvTable->addColumn<int>("isD",Hadron_isD,"isD");
        gvTable->addColumn<int>("fromHF",Hadron_fromHF,"1 if an ancestor in the decay chain is itself an HF/long-lived hadron (B/D/S/Tau)");
        gvTable->addColumn<int>("toHF",Hadron_toHF,"1 if a descendant in the decay chain is itself an HF/long-lived hadron (B/D/S/Tau)");
        gvTable->addColumn<int>("pdgClass",Hadron_pdgId,"pdgClass");
        gvTable->addColumn<float>("minDistNotMatched",Hadron_minDistNotMatched,"Minimum distance to SV among unmatched hadrons");
        gvTable->addColumn<int>("nDaughters",GV_nDaughters,"Total number of GVDaughters belonging to this GV (denominator for nDaughtersMatchedToTracks)"); // NEW
        gvTable->addColumn<float>("maxDaughterPairDeltaR",GV_maxDaughterPairDeltaR,"Max pairwise deltaR among this GV's own GVDaughters (gen-level, no tracks involved)"); // NEW
        gvTable->addColumn<int>("nDaughtersMatchedToTracks",GV_nDaughtersMatchedToTracks,"Number of GVDaughters of this GV matched to a reconstructed track"); // NEW
        gvTable->addColumn<float>("maxDaughterTrkDeltaR",GV_maxDaughterTrkDeltaR,"Max trkDeltaR among this GV's track-matched GVDaughters (not GVDirectDaughters); -1 if none matched"); // NEW
        
        //

        auto dauTable = std::make_unique<nanoaod::FlatTable>(Daughters_pt.size(),"GVDaughters",false);
        dauTable->addColumn<float>("pt",Daughters_pt,"Daughter pt");
        dauTable->addColumn<float>("eta",Daughters_eta,"Daughter eta");
        dauTable->addColumn<float>("phi",Daughters_phi,"Daughter phi");
        dauTable->addColumn<float>("vx", Daughters_vx, "Daughters_vx");
        dauTable->addColumn<float>("vy", Daughters_vy, "Daughters_vy");
        dauTable->addColumn<float>("vz", Daughters_vz, "Daughters_vz");
        dauTable->addColumn<int>("charge",Daughters_charge,"Daughter charge");
        dauTable->addColumn<int>("pdgId",Daughters_pdgId,"Daughter pdgId");
        dauTable->addColumn<int>("hadronIndex",Daughters_GVidx,"Hadron index");
        dauTable->addColumn<int>("trkIdx",Daughters_trkIdx,"Index of matched track in the input track collection, -1 if unmatched");     // NEW
        dauTable->addColumn<int>("isTrkMatched",Daughters_isTrkMatched,"1 if daughter matched to a reconstructed track");                // NEW
        dauTable->addColumn<float>("trkDeltaR",Daughters_trkDeltaR,"deltaR to matched track, -1 if unmatched");                          // NEW
        dauTable->addColumn<float>("trkDPtRel",Daughters_trkDPtRel,"relative pT difference to matched track, -1 if unmatched");          // NEW

        auto directdauTable = std::make_unique<nanoaod::FlatTable>(directDaughters_pt.size(),"GVDirectDaughters",false);
        directdauTable->addColumn<float>("pt",directDaughters_pt,"Daughter pt");
        directdauTable->addColumn<float>("eta",directDaughters_eta,"Daughter eta");
        directdauTable->addColumn<float>("phi",directDaughters_phi,"Daughter phi");
        directdauTable->addColumn<int>("charge",directDaughters_charge,"Daughter charge");
        directdauTable->addColumn<int>("pdgId",directDaughters_pdgId,"Daughter pdgId");
        directdauTable->addColumn<int>("hadronIndex",directDaughters_GVidx,"Hadron index");


        //dauTable->addColumn<int>("hadronFlav",Daughters_flav,"Hadron flavor");


        //
        iEvent.put(std::move(gvTable),"GVTable");
        iEvent.put(std::move(rejectedGVTable),"rejectedGVTable");
        iEvent.put(std::move(dauTable),"GVDaughtersTable");
        iEvent.put(std::move(directdauTable),"GVDirectDaughters");

    }






//  checkPDG() 
int GenVertexProducer::checkPDG(int abs_pdg) const {
    std::vector<int> pdgList_B = {521,511,531,541,5122,5132,5232,5332,5142,5242,5342,5512,5532,5542,5554};
    std::vector<int> pdgList_D = {411,421,431,4122,4232,4132,4332,4412,4422,4432,4444};
    std::vector<int> pdgList_S = {3122,3222,3212,3312,3322,3334};
    std::vector<int> pdgList_Tau = {15};

    if(std::find(pdgList_B.begin(),pdgList_B.end(),abs_pdg)!=pdgList_B.end()) return 1;
    if(std::find(pdgList_D.begin(),pdgList_D.end(),abs_pdg)!=pdgList_D.end()) return 2;
    if(std::find(pdgList_S.begin(),pdgList_S.end(),abs_pdg)!=pdgList_S.end()) return 3;
    if(std::find(pdgList_Tau.begin(),pdgList_Tau.end(),abs_pdg)!=pdgList_Tau.end()) return 4;
    return 0;
}


//  hasHFAncestor() 
// Walks upstream from `hadron` (excluding itself). Returns true if any
// ancestor along the mother(0) chain is itself an HF/long-lived particle
// (B, D, S baryon, or Tau, per the same PDG lists as checkPDG()).
bool GenVertexProducer::hasHFAncestor(const reco::Candidate* hadron) const {
    const reco::Candidate* current = hadron;
    while (current != nullptr && current->numberOfMothers() > 0) {
        const reco::Candidate* mother = current->mother(0);
        if (checkPDG(std::abs(mother->pdgId())) != 0) {
            return true;
        }
        current = mother;
    }
    return false;
}

//  hasHFDescendant() 
// Walks the full daughter tree of `hadron` (all daughters, not just stable
// charged ones). Returns true if any descendant at any depth is itself an
// HF/long-lived particle (B, D, S baryon, or Tau).
bool GenVertexProducer::hasHFDescendant(const reco::Candidate* hadron) const {
    for (size_t i = 0; i < hadron->numberOfDaughters(); ++i) {
        const reco::Candidate* dau = hadron->daughter(i);
        if (checkPDG(std::abs(dau->pdgId())) != 0) {
            return true;
        }
        if (hasHFDescendant(dau)) {
            return true;
        }
    }
    return false;
}

//  isAncestor() 
std::optional<std::tuple<float, float, float>> GenVertexProducer::isAncestor(const reco::Candidate* ancestor, const reco::Candidate* particle) const
    {
    std::vector<int> pdgList_B = {521,511,531,541,5122,5132,5232,5332,5142,5242,5342,5512,5532,5542,5554};
    std::vector<int> pdgList_D = {411,421,431,4122,4232,4132,4332,4412,4422,4432,4444};
    std::vector<int> pdgList_S = {3122,3222,3212,3312,3322,3334};
    std::vector<int> pdgList_Tau = {15};
    std::unordered_set<int> pdgSet_D(pdgList_D.begin(), pdgList_D.end());
    std::unordered_set<int> pdgSet_B(pdgList_B.begin(), pdgList_B.end());
    std::unordered_set<int> pdgSet_S(pdgList_S.begin(), pdgList_S.end());
    std::unordered_set<int> pdgSet_Tau(pdgList_Tau.begin(), pdgList_Tau.end());
    const reco::Candidate* current = particle;
    //const reco::Candidate* child = nullptr;

    while (current != nullptr && current->numberOfMothers() > 0) {
        const reco::Candidate* mother = current->mother(0);
        if (mother == ancestor) {
            // Found the ancestor; return the vertex of the current particle (i.e., the direct daughter)
            return std::make_optional(std::make_tuple(current->vx(), current->vy(), current->vz()));
        }
        int mother_pdg = std::abs(mother->pdgId());
        if (pdgSet_B.count(mother_pdg) || pdgSet_D.count(mother_pdg) || pdgSet_S.count(mother_pdg) || pdgSet_Tau.count(mother_pdg)) break;
        current = mother;
    }

    // If we reached here, the ancestor was not found in the chain
    return std::nullopt;
}


std::vector<std::vector<float>> GenVertexProducer::computeDistanceMatrix(
                const std::vector<float>& SV_x,
                const std::vector<float>& SV_y,
                const std::vector<float>& SV_z,
                std::vector<CovMatrix> SV_cov,
                const std::vector<float>& Hadron_GVx,
                const std::vector<float>& Hadron_GVy,
                const std::vector<float>& Hadron_GVz) {
    // computeDistanceMatrix
    // Returns
    // distances = Matrix of Euclidean distances between SV and GV
    
    
    size_t nSV = SV_x.size();
    size_t nHadron = Hadron_GVx.size();
    
    // 2D vector initialized to 999 nSV x nHadron
    std::vector<std::vector<float>> distances(nSV, std::vector<float>(nHadron, 999.0));

    for (size_t i = 0; i < nSV; ++i) {
        //CovMatrix covInv = SV_cov[i].Inverse();
        CovMatrix covInv = SV_cov[i];
        covInv.Invert();
        for (size_t j = 0; j < nHadron; ++j) {
            float dx = SV_x[i] - Hadron_GVx[j];
            float dy = SV_y[i] - Hadron_GVy[j];
            float dz = SV_z[i] - Hadron_GVz[j];
            float chi2 =
              dx * (covInv(0,0) * dx +
                    covInv(0,1) * dy +
                    covInv(0,2) * dz)

            + dy * (covInv(1,0) * dx +
                    covInv(1,1) * dy +
                    covInv(1,2) * dz)

            + dz * (covInv(2,0) * dx +
                    covInv(2,1) * dy +
                    covInv(2,2) * dz);

            float dist = std::sqrt(chi2);
            //float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            distances[i][j] = dist;
        }
    }

    return distances;
}


void GenVertexProducer::printDistanceMatrix(
    const std::vector<std::vector<float>>& distances
) {
    size_t nSV = distances.size();
    if (nSV == 0) return;

    size_t nGV = distances[0].size();

    std::cout << "\n Distance Matrix (SV rows × GV cols):\n\n";

    // Print header row
    std::cout << std::setw(8) << "SV/GV";
    for (size_t j = 0; j < nGV; ++j) {
        std::cout << std::setw(10) << "GV[" + std::to_string(j) + "]";
    }
    std::cout << "\n";

    // Print matrix values
    for (size_t i = 0; i < nSV; ++i) {
        std::cout << std::setw(8) << "SV[" + std::to_string(i) + "]";
        for (size_t j = 0; j < nGV; ++j) {
            std::cout << std::setw(10) << std::fixed << std::setprecision(3) << distances[i][j];
        }
        std::cout << "\n";
    }
}
float computeDR_SV_Had(int bestSV,
                       int bestHad,
                       const std::vector<float>& SV_eta,
                       const std::vector<float>& SV_phi,
                       const std::vector<float>& GV_eta,
                       const std::vector<float>& GV_phi)
{
    float dEta = SV_eta[bestSV] - GV_eta[bestHad];

    float dPhi = SV_phi[bestSV] - GV_phi[bestHad];

    // wrap phi into [-pi, pi]
    while (dPhi > M_PI)  dPhi -= 2.0 * M_PI;
    while (dPhi < -M_PI) dPhi += 2.0 * M_PI;

    return std::sqrt(dEta * dEta + dPhi * dPhi);
}

std::tuple<std::vector<int>, std::vector<float>, std::vector<float>, std::vector<int>, std::vector<int>, std::vector<int>>  GenVertexProducer::matchHadronsToSV(
    std::vector<std::vector<float>> distances,
    const std::vector<float>& SVtrk_pt,
    const std::vector<float>& SVtrk_eta,
    const std::vector<float>& SVtrk_phi,
    const std::vector<int>& SVtrk_SVidx,
    const std::vector<float>& Daughters_pt,     //genparticles
    const std::vector<float>& Daughters_eta,  //genparticles
    const std::vector<float>& Daughters_phi,  //genparticles
    const std::vector<int>& Daughters_GVidx, // hadron index per daughter
    const std::vector<float>& SV_eta,
    const std::vector<float>& SV_phi,
    const std::vector<float>& GV_eta,
    const std::vector<float>& GV_phi,
    int n_Hadrons,
    int nRequiredCommonTracks,
    double dR_max,
    double relPt_max,
    bool doubleMatching,
    int doubleMatching_nRequiredCommonTracks,
    double doubleMatching_maxSignificance,
    double doubleMatching_dR_max,
    double doubleMatching_relPt_max
) {
    // Returns
    // Hadron_SVIdx = Array of lenght = Hadron_pt.size() with index of the SV
    size_t nSV = distances.size(); //distances is nSV x nHadron matrix (the first dimension is nSV)
    std::vector<int> Hadron_SVIdx(n_Hadrons, -1); // Output
    std::vector<float> Hadron_SVDistance(n_Hadrons, -1); // Output
    std::vector<int> SVtrk_isMatched(SVtrk_pt.size(), 0);   // NEW
    std::vector<int> SVtrk_GVIdx(SVtrk_pt.size(), -1);      // NEW
    std::vector<int> SVtrk_daughterIdx(SVtrk_pt.size(), -1);
    std::vector<size_t> svTrackIdxs_fromBestSV;
    std::vector<std::vector<float>> distancesOriginal = distances;
    while (true) {
        float minDist = 999.0;
        int bestSV = -1;
        int bestHad = -1;

        // Find minimum distance in current matrix
        // store in 
        // - minDist
        // - bestSV
        // - bestHad
        for (size_t sv = 0; sv < nSV; ++sv) {
            for (int had = 0; had < n_Hadrons; ++had) {
                if (distances[sv][had] < minDist) {
                    minDist = distances[sv][had];
                    bestSV = sv;
                    bestHad = had;
                }
            }
        }
        //std::cout<<"\n Considering best pair: SV["<<bestSV<<"] and Hadron["<<bestHad<<"] with distance "<<minDist<<std::endl;
        if (minDist >= 997.0) break;  // done

        // Select tracks from SV
        // svTrackIdxs_fromBestSV is initialized every time
        // it stores the index of the tracks which originate from SV candidate in this loop
        
        svTrackIdxs_fromBestSV.clear();
        for (size_t i = 0; i < SVtrk_SVidx.size(); ++i) {
            // among all tracks from all SV, select those from the candidate SV
            //std::cout<<" Track index "<<i<<" SVtrk_SVidx: "<<SVtrk_SVidx[i]<<" SVtrk_pt: "<<SVtrk_pt[i]<<" Best SV :"<<bestSV<<std::endl;
            if (SVtrk_SVidx[i] == bestSV && SVtrk_pt[i] > 0.8 && std::fabs(SVtrk_eta[i]) < 2.5) {
                svTrackIdxs_fromBestSV.push_back(i);
            }
        }

        // Select daughters of Hadron
        std::vector<size_t> GenDaughtersIdxs_fromBestHad;
        for (size_t i = 0; i < Daughters_GVidx.size(); ++i) {
            if (Daughters_GVidx[i] == bestHad) {
                GenDaughtersIdxs_fromBestHad.push_back(i);
            }
        }

        // Match logic: check for 1 (2) or more matched tracks by ΔR & dPt/pT
        //int nRequiredCommonTracks = 1;
        int common = 0;
        //std::vector<size_t> matchedTrackIdxs; // global SVtrk indices matched in this SV/Hadron trial
        std::vector<std::pair<size_t,size_t>> matchedTrackToDaughter; // 
        for (size_t iSV : svTrackIdxs_fromBestSV) {
            bool trackMatched = false;
            size_t matchedDauIdx = 0;
            //std::cout<<" Checking SVtrack index "<<iSV<<std::endl;
            for (size_t iHad : GenDaughtersIdxs_fromBestHad) {
                //std::cout<<" Checking GVDaughters index "<<iHad<<std::endl;
                float dR = deltaR(SVtrk_eta[iSV], SVtrk_phi[iSV], Daughters_eta[iHad], Daughters_phi[iHad]);
                float relPt = std::fabs(SVtrk_pt[iSV] - Daughters_pt[iHad]) / Daughters_pt[iHad];
                //std::cout<<"Comparing SV track (pt: "<<SVtrk_pt[iSV]<<", eta: "<<SVtrk_eta[iSV]<<", phi: "<<SVtrk_phi[iSV]<<") with Daughter (pt: "<<Daughters_pt[iHad]<<", eta: "<<Daughters_eta[iHad]<<", phi: "<<Daughters_phi[iHad]<<") => dR: "<<dR<<", relPt: "<<relPt<<std::endl;
                if (dR < dR_max && relPt < relPt_max) {
                    
                    trackMatched = true;
                    matchedDauIdx = iHad;
                    break;
                    //std::cout<<"  -> Matched! Common tracks: "<<common<<std::endl;
                    //if (common >= nRequiredCommonTracks) break; // break the iHad cycle
                }
            }
            if (trackMatched) {
                ++common;
                matchedTrackToDaughter.push_back({iSV, matchedDauIdx});
            }
        }

        if (common >= nRequiredCommonTracks) {
            Hadron_SVIdx[bestHad] = bestSV;
            Hadron_SVDistance[bestHad] = minDist;
            for (auto& [iSV, iHad] : matchedTrackToDaughter) {
                SVtrk_isMatched[iSV] = 1;
                SVtrk_GVIdx[iSV] = bestHad;
                SVtrk_daughterIdx[iSV] = iHad;
            }
            for (int h = 0; h < n_Hadrons; ++h) distances[bestSV][h] = 1000.0;
            for (size_t s = 0; s < nSV; ++s) distances[s][bestHad] = 1000.0;
        } else {
            distances[bestSV][bestHad] = 998.0;
        }
    }




    //printDistanceMatrix(distances);
    if (doubleMatching){

        while (true){
            float minDist = 999.0;
            int bestSV = -1;
            int bestHad = -1;
            for (size_t sv = 0; sv < nSV; ++sv) {
                for (int had = 0; had < n_Hadrons; ++had) {
                    if (distances[sv][had] < minDist) {
                        minDist = distances[sv][had];
                        bestSV = sv;
                        bestHad = had;
                    }
                }
            }

            if (minDist >= 998.5) break;  // done
            //check whether deltaR bewteen bestSV and bestHad is less than doublematching_dR_max
            float dR_SV_Had = computeDR_SV_Had(bestSV, bestHad, SV_eta, SV_phi, GV_eta, GV_phi);
            

            svTrackIdxs_fromBestSV.clear();
            for (size_t i = 0; i < SVtrk_SVidx.size(); ++i) {
                // among all tracks from all SV, select those from the candidate SV
                //std::cout<<" Track index "<<i<<" SVtrk_SVidx: "<<SVtrk_SVidx[i]<<" SVtrk_pt: "<<SVtrk_pt[i]<<" Best SV :"<<bestSV<<std::endl;
                if (SVtrk_SVidx[i] == bestSV && SVtrk_pt[i] > 0.8 && std::fabs(SVtrk_eta[i]) < 2.5) {
                    svTrackIdxs_fromBestSV.push_back(i);
                }
            }

            // Select daughters of Hadron
            std::vector<size_t> GenDaughtersIdxs_fromBestHad;
            for (size_t i = 0; i < Daughters_GVidx.size(); ++i) {
                if (Daughters_GVidx[i] == bestHad) {
                    GenDaughtersIdxs_fromBestHad.push_back(i);
                }
            }

            // Match logic: check for 1 (2) or more matched tracks by ΔR & dPt/pT
            //int nRequiredCommonTracks = 1;
            int common = 0;
            for (size_t iSV : svTrackIdxs_fromBestSV) {
                //std::cout<<" Checking SVtrack index "<<iSV<<std::endl;
                for (size_t iHad : GenDaughtersIdxs_fromBestHad) {
                    //std::cout<<" Checking GVDaughters index "<<iHad<<std::endl;
                    float dR = deltaR(SVtrk_eta[iSV], SVtrk_phi[iSV], Daughters_eta[iHad], Daughters_phi[iHad]);
                    float relPt = std::fabs(SVtrk_pt[iSV] - Daughters_pt[iHad]) / Daughters_pt[iHad];
                    //std::cout<<"Comparing SV track (pt: "<<SVtrk_pt[iSV]<<", eta: "<<SVtrk_eta[iSV]<<", phi: "<<SVtrk_phi[iSV]<<") with Daughter (pt: "<<Daughters_pt[iHad]<<", eta: "<<Daughters_eta[iHad]<<", phi: "<<Daughters_phi[iHad]<<") => dR: "<<dR<<", relPt: "<<relPt<<std::endl;
                    if (dR < doubleMatching_dR_max && relPt < doubleMatching_relPt_max) {
                        ++common;
                        //std::cout<<"  -> Matched! Common tracks: "<<common<<std::endl;
                        if (common >= 1) break; // break the iHad cycle
                    }
                }
                if (common >= doubleMatching_nRequiredCommonTracks) break; // break the iSV cycle
            }
            //std::cout<<"\n [DoubleMatching] Pair: SV["<<bestSV<<"] and Hadron["<<bestHad<<"] with distance "<<distancesOriginal[bestSV][bestHad]<<" and dR "<<dR_SV_Had<<" and common tracks "<<common<<std::endl;
            if (common >= doubleMatching_nRequiredCommonTracks && distancesOriginal[bestSV][bestHad] < doubleMatching_maxSignificance ) {
            //if (common >= 1 && distancesOriginal[bestSV][bestHad] < doubleMatching_maxSignificance && dR_SV_Had < doubleMatching_dR_max ) {

                Hadron_SVIdx[bestHad] = bestSV;
                Hadron_SVDistance[bestHad]= -distancesOriginal[bestSV][bestHad];
                for (int h = 0; h < n_Hadrons; ++h) distances[bestSV][h] = 1000.0;  // remove SV row (MATCHED)
                for (size_t s = 0; s < nSV; ++s) distances[s][bestHad] = 1000.0;   // remove Hadron column (MATCHED)
                //std::cout << "[V] Matched Hadron[" << bestHad << "] to SV[" << bestSV << "] (distance = " << distancesOriginal[bestSV][bestHad] << ", common tracks = " << common << ")\n";
            } else {
                distances[bestSV][bestHad] = 999.0;  // exclude this pair
            }



        }

    }

    // check the minDist not mathced using original matrix (not modified by matching procedure)
    std::vector<float> minDistNotMatched(n_Hadrons, 999.);
    for (int had = 0; had < n_Hadrons; ++had) {
        for (size_t sv = 0; sv < nSV; ++sv) {
            if (distancesOriginal[sv][had] < minDistNotMatched[had] && Hadron_SVIdx[had] == -1) {
                minDistNotMatched[had] = distancesOriginal[sv][had];
            }
        }
    }
    return std::make_tuple(Hadron_SVIdx, Hadron_SVDistance, minDistNotMatched, SVtrk_isMatched, SVtrk_GVIdx, SVtrk_daughterIdx);
    }


// NEW: matchDaughtersToTracks()
// Matches gen-level GVDaughters to reconstructed tracks by deltaR and relative pT,
// following the same TrackGenMatcher-style logic (optional charge check, optional
// one-to-one "resolveAmbiguities" assignment).
//
// Strategy: build the list of all (daughter, track) pairs that pass the deltaR/relPt
// (and, if requested, charge) cuts, sort them by ascending deltaR, then greedily assign
// pairs. When resolveAmbiguities is true, once a track is used it can't be reused by
// another daughter (one-to-one matching, analogous to TrackGenMatcher's resolveAmbiguities).
// When false, several daughters may be matched to the same track.
std::tuple<std::vector<int>, std::vector<int>, std::vector<float>, std::vector<float>>
GenVertexProducer::matchDaughtersToTracks(
    const std::vector<float>& Daughters_pt,
    const std::vector<float>& Daughters_eta,
    const std::vector<float>& Daughters_phi,
    const std::vector<int>& Daughters_charge,
    const std::vector<reco::Track>& tracks,
    double maxDeltaR,
    double maxDPtRel,
    bool checkCharge,
    bool resolveAmbiguities) const
{
    size_t nDau = Daughters_pt.size();
    size_t nTrk = tracks.size();

    std::vector<int>   trkIdx(nDau, -1);
    std::vector<int>   isMatched(nDau, 0);
    std::vector<float> matchDeltaR(nDau, -1.f);
    std::vector<float> matchDPtRel(nDau, -1.f);

    struct Pair { float dR; size_t dau; size_t trk; };
    std::vector<Pair> pairs;
    pairs.reserve(nDau);

    for (size_t i = 0; i < nDau; ++i) {
        for (size_t j = 0; j < nTrk; ++j) {
            if (checkCharge && tracks[j].charge() != Daughters_charge[i]) continue;

            float dR = deltaR(Daughters_eta[i], Daughters_phi[i], tracks[j].eta(), tracks[j].phi());
            if (dR >= maxDeltaR) continue;

            float dPtRel = std::fabs(tracks[j].pt() - Daughters_pt[i]) / Daughters_pt[i];
            if (dPtRel >= maxDPtRel) continue;

            pairs.push_back({dR, i, j});
        }
    }

    // best (smallest deltaR) candidates get assigned first
    std::sort(pairs.begin(), pairs.end(), [](const Pair& a, const Pair& b) { return a.dR < b.dR; });

    std::vector<bool> dauUsed(nDau, false);
    std::vector<bool> trkUsed(nTrk, false);
    for (const auto& p : pairs) {
        if (dauUsed[p.dau]) continue;                        // each daughter is matched at most once
        if (resolveAmbiguities && trkUsed[p.trk]) continue;   // one-to-one if requested

        trkIdx[p.dau]      = static_cast<int>(p.trk);
        isMatched[p.dau]   = 1;
        matchDeltaR[p.dau] = p.dR;
        matchDPtRel[p.dau] = std::fabs(tracks[p.trk].pt() - Daughters_pt[p.dau]) / Daughters_pt[p.dau];

        dauUsed[p.dau] = true;
        if (resolveAmbiguities) trkUsed[p.trk] = true;
    }

    return std::make_tuple(trkIdx, isMatched, matchDeltaR, matchDPtRel);
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenVertexProducer);