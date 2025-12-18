// -*- Header -*-
//
// Package:       Analyzers/ChFluctuations in MiniAOD and AOD for PbPb data
// Class:         ChFluctuations
//
//
// Author:       Subash Behera
// Created on: 8th April 2025



#ifndef Analyzers_ChFluctuations_ChFluctuations_h
#define Analyzers_ChFluctuations_ChFluctuations_h

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "THnSparse.h"

// Centrality Includes
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include <vector>

class ChFluctuations : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ChFluctuations(const edm::ParameterSet&);
      ~ChFluctuations();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      void LoopCMWVertices(const edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------
      
      // Tokens
      edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trackTags_;
      edm::EDGetTokenT<edm::View<pat::PackedGenParticle>> trackTagsgen_;
      edm::EDGetTokenT<edm::ValueMap<float>> chi2Map_;
      edm::EDGetTokenT<reco::VertexCollection> vtxTags_;

      // Centrality & Tower Tokens
      edm::EDGetTokenT<reco::Centrality> centralityTag_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_; 
      std::vector<double> centralityMap_;

      // Config parameters
      edm::InputTag fGeneral_;
      edm::InputTag fGeneral2_;
      double zminVtx_;
      double zmaxVtx_;
      double xBestVtx_;
      double yBestVtx_;
      double zBestVtx_;
      bool fIsMC;
      double fPtMin;
      double fPtMax;
      double fEtaMax;
      
      std::vector<double> pTBins_;
      std::vector<double> etaBins_;
      std::vector<double> centbins_;
      std::vector<int> algoParameters_;

      // Histograms
      TH1F *hZBestVtx;
      TH1F *hihf_pf; // Renamed from hHFEnergy
      TH1F *hNch;    // Added Multiplicity Plot
      
      TH1F *hcent_bin;
      TH1F *hcentbin;
      TH1F *hptbin;
      TH1F *hetabin;

      std::vector<TH1F*> hpt_gen;
      std::vector<TH1F*> heta_gen;
      std::vector<TH1F*> hphi_gen;

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

      THnSparseI *fThnCentPosNegEta;
      THnSparseI *fThnCentPosNegEtaMC;
      TH1D *fhEvent;
      TH1D *fhCent;
      TH1D *fhEta;
      TH1D *fhVz;
      TH1D *fhPhi;
      TH1D *fhPtPlus;
      TH1D *fhPtMinus;
      TH1D *fheta;

      TFile *feff_;
      TFile *feff2_;
      TH3F *hEff3D;
      TH3F *hFak3D;
      TH3F *hSec3D;
      TH3F *hMul3D;
};
#endif
