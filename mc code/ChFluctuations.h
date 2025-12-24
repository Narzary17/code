#ifndef Analyzers_Cumulants_ChFluctuations_h
#define Analyzers_Cumulants_ChFluctuations_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"

#include <vector>
#include <string>

class ChFluctuations : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ChFluctuations(const edm::ParameterSet&);
      ~ChFluctuations();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // Tokens
      edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trackTags_;
      edm::EDGetTokenT<edm::View<pat::PackedGenParticle>> trackTagsgen_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxTags_;
      edm::EDGetTokenT<reco::Centrality> centralityTags_;

      // Parameters
      int noffmin_;
      int noffmax_;
      double ptnoffmin_;
      double ptnoffmax_;
      double dzdzerrornoff_;
      double d0d0errornoff_;
      double pterrorptnoff_;

      double etamin_;
      double etamax_;
      double ptmin_;
      double ptmax_;

      edm::InputTag fname_;
      bool cweight_;
      bool fIsMC;

      std::vector<double> centbins_;
      std::vector<double> etaBins_;

      // Efficiency
      TFile* feff_;
      TH3F* hEff_3D;
      TH3F* hFak3D;
      TH3F* hMul_3D;
      TH3F* hSec_3D;

      // --- Histograms (FIXED) ---
      TH1F* hihf_pf;     // Added
      TH1F* hNoff_;      // Added
      TH1F* hNch;        // Added
      TH1F* hZBestVtx_;  // Renamed from hZBestVtx to match .cc logic

      TH1F* hcent_bin;
      TH1F* hcentbin;
      TH1F* hptbin;
      TH1F* hetabin;

      std::vector<TH1F*> hpt;
      std::vector<TH1F*> heta;
      std::vector<TH1F*> hphi;
      std::vector<TH1F*> hnHits;
      std::vector<TH1F*> hptreso;
      std::vector<TH1F*> hchi2;
      std::vector<TH1F*> hDCAZ;
      std::vector<TH1F*> hDCAXY;

      std::vector<TH1F*> hpt_eff;
      std::vector<TH1F*> heta_eff;
      std::vector<TH1F*> hphi_eff;
      std::vector<TH1F*> hnHits_eff;
      std::vector<TH1F*> hptreso_eff;
      std::vector<TH1F*> hchi2_eff;
      std::vector<TH1F*> hDCAZ_eff;
      std::vector<TH1F*> hDCAXY_eff;

      std::vector<TH1F*> hpt_gen;
      std::vector<TH1F*> heta_gen;
      std::vector<TH1F*> hphi_gen;

      // Global
      THnSparseI* fThnCentPosNegEta;
      THnSparseI* fThnCentPosNegEtaMC;
      
      TH1D* fhEvent;
      TH1D* fhCent;
      TH1D* fhEta;
      TH1D* fhVz;
      TH1D* fhPhi;
      TH1D* fhPtPlus;
      TH1D* fhPtMinus;
      TH1D* fheta;

      const int nCentTableBins = 200;
      float binBoundaries[201]; 
};

#endif
