// -*- C++ -*-
//
//         Author: Subash Behera (Updated for OO)
//         Class: ChFluctuations
//

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Analyzers/ChFluctuations/interface/ChFluctuations.h"
#include <string>
#include <vector>
#include <algorithm>

ChFluctuations::ChFluctuations(const edm::ParameterSet& iConfig) :
  // 1. Tokens
  trackTags_(consumes< edm::View< pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("tracks"))), 
  trackTagsgen_(consumes< edm::View< pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("tracksgen"))), 
  chi2Map_( consumes< edm::ValueMap< float > >( iConfig.getParameter< edm::InputTag >( "trackschi2" ) ) ), 
  vtxTags_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex"))),

  // 2. Centrality
  centralityTag_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralityRaw"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("srcTower"))), 
  centralityMap_(iConfig.getUntrackedParameter<std::vector<double>>("centralityMap")),

  // 3. Efficiency
  fGeneral_(iConfig.getUntrackedParameter<edm::InputTag>("fGeneral")),
  fGeneral2_(iConfig.getUntrackedParameter<edm::InputTag>("fGeneral2")),

  // 4. Cuts
  zminVtx_(iConfig.getUntrackedParameter<double>("zminVtx")),
  zmaxVtx_(iConfig.getUntrackedParameter<double>("zmaxVtx")),
  fIsMC( iConfig.getUntrackedParameter<bool>("IsMC") ), 
  fPtMin( iConfig.getUntrackedParameter<double>("ptMin") ),
  fPtMax( iConfig.getUntrackedParameter<double>("ptMax") ),
  fEtaMax( iConfig.getUntrackedParameter<double>("etaMax") ),

  // 5. Bins
  pTBins_(iConfig.getUntrackedParameter< std::vector < double > >("pTBins")),
  etaBins_(iConfig.getUntrackedParameter< std::vector < double > >("etaBins")),
  centbins_(iConfig.getUntrackedParameter< std::vector < double > >("centbins")),
  algoParameters_(iConfig.getParameter< std::vector < int > >("algoParameters"))
{
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  feff_ = 0x0;
  TString f_General(fGeneral_.label().c_str());
  if(!f_General.IsNull()) {
      edm::FileInPath f1 (Form("Analyzers/ChFluctuations/data/EFF/general_tracks/%s",f_General.Data()));
      feff_ = new TFile(f1.fullPath().c_str(),"READ");
  }

  feff2_ = 0x0;
  TString f_General2(fGeneral2_.label().c_str());
  if(!f_General2.IsNull()) {
      edm::FileInPath f2 (Form("Analyzers/ChFluctuations/data/EFF/general_tracks/%s",f_General2.Data()));
      feff2_ = new TFile(f2.fullPath().c_str(),"READ");
  }

  hEff3D = nullptr; hFak3D = nullptr; hSec3D = nullptr; hMul3D = nullptr;

  if(feff2_) hEff3D = (TH3F*)feff2_->Get("hEff_3D");
  if(feff_) {
      hFak3D = (TH3F*)feff_->Get("hFak_3D");
      hSec3D = (TH3F*)feff_->Get("hSec_3D");
      hMul3D = (TH3F*)feff_->Get("hMul_3D");
  }
  if(feff_) feff_->Close();
  if(feff2_) feff2_->Close();

  usesResource("TFileService");
  edm::Service<TFileService> fss;

  TFileDirectory fGlobalHist1 = fss->mkdir("QAplots");
  hZBestVtx    = fGlobalHist1.make<TH1F>("hZvtx", "", 600, -30.0, 30.0);
  
  // hihf_pf
  hihf_pf      = fGlobalHist1.make<TH1F>("hihf_pf", "", 1000, 0.0, 5000.0);
  
  // Nch
  hNch         = fGlobalHist1.make<TH1F>("hNch", "", 500, 0.0, 500.0);

  hcent_bin    = fGlobalHist1.make<TH1F>("hcent_bin", "", 200, 0.0, 200.0);
  hcentbin     = fGlobalHist1.make<TH1F>("hcentbin", "", centbins_.size()-1, &centbins_[0]);
  hptbin     = fGlobalHist1.make<TH1F>("hptbin", "", 500, 0., 50.);
  hetabin     = fGlobalHist1.make<TH1F>("hetabin", "", etaBins_.size()-1, &etaBins_[0]);

  hpt_gen.resize(centbins_.size());
  heta_gen.resize(centbins_.size());
  hphi_gen.resize(centbins_.size());
  hpt.resize(centbins_.size());
  heta.resize(centbins_.size());
  hphi.resize(centbins_.size());
  hnHits.resize(centbins_.size());
  hptreso.resize(centbins_.size());
  hchi2.resize(centbins_.size());
  hDCAZ.resize(centbins_.size());
  hDCAXY.resize(centbins_.size());
  hpt_eff.resize(centbins_.size());
  heta_eff.resize(centbins_.size());
  hphi_eff.resize(centbins_.size());
  hnHits_eff.resize(centbins_.size());
  hptreso_eff.resize(centbins_.size());
  hchi2_eff.resize(centbins_.size());
  hDCAZ_eff.resize(centbins_.size());
  hDCAXY_eff.resize(centbins_.size());

  for(unsigned int i = 0; i < centbins_.size(); i++) {
      hpt_gen[i]          = fGlobalHist1.make<TH1F>(Form("hpt_gen_%d", i), "", 500, 0., 50.);
      heta_gen[i]         = fGlobalHist1.make<TH1F>(Form("heta_gen_%d",i), "", etaBins_.size()-1, &etaBins_[0]);
      hphi_gen[i]         = fGlobalHist1.make<TH1F>(Form("hphi_gen_%d",i), "", 62, -TMath::Pi(), TMath::Pi());
      hpt[i]              = fGlobalHist1.make<TH1F>(Form("hpt_%d", i), "", 500, 0., 50.);
      heta[i]             = fGlobalHist1.make<TH1F>(Form("heta_%d", i), "", etaBins_.size()-1, &etaBins_[0]);
      hphi[i]             = fGlobalHist1.make<TH1F>(Form("hphi_%d", i), "", 62, -TMath::Pi(), TMath::Pi());
      hnHits[i]           = fGlobalHist1.make<TH1F>(Form("hnHits_%d", i), "", 50, 0., 50.);
      hptreso[i]          = fGlobalHist1.make<TH1F>(Form("hptreso_%d", i), "", 100, 0., 1.);
      hchi2[i]            = fGlobalHist1.make<TH1F>(Form("hchi2_%d", i), "", 100, 0, 10.);
      hDCAZ[i]            = fGlobalHist1.make<TH1F>(Form("hDCAZ_%d", i), "", 80, -4., 4.);
      hDCAXY[i]           = fGlobalHist1.make<TH1F>(Form("hDCAXY_%d", i), "", 80, -4., 4.);
      hpt_eff[i]          = fGlobalHist1.make<TH1F>(Form("hpt_eff_%d", i), "", 500, 0., 50.);
      heta_eff[i]         = fGlobalHist1.make<TH1F>(Form("heta_eff_%d", i), "", etaBins_.size()-1, &etaBins_[0]);
      hphi_eff[i]         = fGlobalHist1.make<TH1F>(Form("hphi_eff_%d", i), "", 62, -TMath::Pi(), TMath::Pi());
      hnHits_eff[i]       = fGlobalHist1.make<TH1F>(Form("hnHits_eff_%d", i), "", 50, 0., 50.);
      hptreso_eff[i]      = fGlobalHist1.make<TH1F>(Form("hptreso_eff_%d", i), "", 100, 0., 1.);
      hchi2_eff[i]        = fGlobalHist1.make<TH1F>(Form("hchi2_eff_%d", i), "", 100, 0, 10.);
      hDCAZ_eff[i]        = fGlobalHist1.make<TH1F>(Form("hDCAZ_eff_%d", i), "", 80, -4., 4.);
      hDCAXY_eff[i]       = fGlobalHist1.make<TH1F>(Form("hDCAXY_eff_%d", i), "", 80, -4., 4.);
  }

  edm::Service<TFileService> fs;
  TFileDirectory fGlobalHist  = fs->mkdir("Global");
  const int kdim = 25; 
  int bin[kdim];
  double min[kdim], max[kdim];
  bin[0] = 81;
  for(int ibin = 1; ibin < kdim; ibin++) bin[ibin] = 3000;
  for(int ibin = 0; ibin < kdim; ibin++) min[ibin] = -0.5;
  max[0] = 80.5;
  for(int ibin = 1; ibin < kdim; ibin++) max[ibin] = 2999.5;

  fThnCentPosNegEta = fGlobalHist.make<THnSparseI>("fThnCentPosNegEta", "Cent-Etabinwise-Nplus-Nminus", kdim, bin, min, max);
  fThnCentPosNegEtaMC = fGlobalHist.make<THnSparseI>("fThnCentPosNegEtaMC", "Cent-Etabinwise-Nplus-Nminus", kdim, bin, min, max);
  
  fhEvent = fGlobalHist.make<TH1D>("fhEvent", "Event counter", 8, -0.5, 7.5 );
  fhCent = fGlobalHist.make<TH1D>("fhCent", "centrality", 200, 0, 200); 
  fhEta = fGlobalHist.make<TH1D>("fhEta", "EtaBining", etaBins_.size()-1, &etaBins_[0]);
  fhVz = fGlobalHist.make<TH1D>("fhVz", "V_{z} distribution", 320, -16.0, 16.0);
  fhPhi = fGlobalHist.make<TH1D>("fhPhi", "#phi distributions", 640, -3.2, 3.2);
  fhPtPlus = fGlobalHist.make<TH1D>("fhPtPlus", "N^{+} p_{T} distributions", 120, 0., 6.0);
  fhPtMinus = fGlobalHist.make<TH1D>("fhPtMinus", "N^{-} p_{T} distributions", 120, 0., 6.0);
  fheta = fGlobalHist.make<TH1D>("fheta", "#eta distribution", 100, -2.5, 2.5);
}

ChFluctuations::~ChFluctuations() {}
void ChFluctuations::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) { LoopCMWVertices(iEvent, iSetup); }
void ChFluctuations::beginJob() {}
void ChFluctuations::endJob() {}
void ChFluctuations::fillDescriptions(edm::ConfigurationDescriptions&  descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void ChFluctuations::LoopCMWVertices( const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;

  auto trks = iEvent.getHandle( trackTags_ );
  auto chi2Map = iEvent.getHandle( chi2Map_ );
  fhEvent->Fill(0);
  auto pvs = iEvent.getHandle( vtxTags_ );

  double bestvzError, bestvxError, bestvyError; 
  math::XYZPoint bestvtx;
  if ( !pvs->empty() ) {
      const reco::Vertex& vtx = (*pvs)[0];
      
      // 1. Primary Vertex Filter 
      if (vtx.isFake()) return;

      bestvzError = vtx.zError();
      bestvxError = vtx.xError();
      bestvyError = vtx.yError();
      bestvtx = vtx.position();
  } else { return; }

  // 2. Pileup Filter 
  int nGoodVertices = 0;
  for(auto const& v : *pvs) {
      if( !v.isFake() && v.ndof() > 4 && fabs(v.z()) <= 25.0 && v.position().Rho() <= 2.0 ) {
          nGoodVertices++;
      }
  }
  if ( nGoodVertices > 1 ) return; // Reject Pileup

  xBestVtx_ = bestvtx.x();
  yBestVtx_ = bestvtx.y();
  zBestVtx_ = bestvtx.z();
  fhEvent->Fill(0);
  if ( zBestVtx_ < zminVtx_ || zBestVtx_ >= zmaxVtx_ ) return;

  // ----------------- MANUAL CENTRALITY -------------------------------
  fhEvent->Fill(2);

  double etHFtowerSumPlus = 0;
  double etHFtowerSumMinus = 0;
  float hf = 0.0;

  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByToken(pfToken_, pfs);

  for(unsigned int i = 0, n = pfs->size(); i < n; ++i) {
      const pat::PackedCandidate &pf = (*pfs)[i];
      // pdgId 1 and 2 for hadronic and em particles in HF
      if(pf.pdgId() == 1 || pf.pdgId() == 2){
          if(pf.eta() > 3.0 && pf.eta() < 5.2){
              etHFtowerSumPlus += pf.et();
          }       
          if(pf.eta() < -3.0 && pf.eta() > -5.2){
              etHFtowerSumMinus += pf.et();
          }
      } 
  } 
  hf = etHFtowerSumPlus + etHFtowerSumMinus;

  // *** FILL hihf_pf ***
  hihf_pf->Fill(hf);

  int cbinVal = -1;
  int nBins = centralityMap_.size() - 1; 

  for(size_t i = 0; i < centralityMap_.size() - 1; ++i) {
      if(hf >= centralityMap_[i] && hf < centralityMap_[i+1]) {
          cbinVal = nBins - 1 - i; 
          break;
      }
  }
  if(hf >= centralityMap_.back()) cbinVal = 0;
  float centralityBin = (float)cbinVal;
  
  if(centralityBin < 0) return; 

  int centbin = -1;
  for(unsigned int i=0; i < centbins_.size()-1; i++) {
      if(centralityBin >= centbins_[i] && centralityBin < centbins_[i+1]) {
          centbin = i;
          break;
      }
  }

  fhEvent->Fill(3);
  hZBestVtx -> Fill(zBestVtx_);
  hcent_bin -> Fill(centralityBin);
  hcentbin -> Fill(centralityBin); 
  fhCent-> Fill(centralityBin);

  if( centbin < 0 ) return; 

  // --- GEN LOOP ---
  const int kdim = 25;
  double NplusMC = 0, NminusMC = 0;
  double NplusNminusMC[24];
  double ThnContainerMC[kdim];
  for(int idx = 0; idx < 24; idx++) NplusNminusMC[idx] = 0.;
  for(int jdx = 0; jdx < kdim; jdx++) ThnContainerMC[jdx] = 0.;

  if(!fIsMC){ 
    edm::Handle<edm::View<pat::PackedGenParticle>> trksgen;
    iEvent.getByToken(trackTagsgen_, trksgen);

    if(trksgen.isValid()) {
        for (auto const& iter_tk_gn : *trksgen) {
          if(iter_tk_gn.status() != 1) continue;
          double gen_pt = iter_tk_gn.pt();
          double gen_eta = iter_tk_gn.eta();
          int gen_charge = iter_tk_gn.charge();
          double gen_phi = iter_tk_gn.phi();

          if( gen_charge == 0 ) continue;
          if( gen_pt <= 0.3 || gen_pt >= 3.0 ) continue;
          if( gen_eta <= -2.4 || gen_eta >= 2.4 ) continue;

          hpt_gen[centbin]->Fill(gen_pt);
          heta_gen[centbin]->Fill(gen_eta);
          hphi_gen[centbin]->Fill(gen_phi);

          if( gen_pt < fPtMin || gen_pt > fPtMax ) continue;
          if( TMath::Abs(gen_eta) > fEtaMax ) continue;

          Int_t etabin = fhEta->FindBin( gen_eta ) - 1;
          if( etabin < 0 || etabin >= 12 ) continue; 

          if (gen_charge > 0 ) { NplusMC +=1.; NplusNminusMC[etabin] += 1.; }
          else if(gen_charge < 0 ) { NminusMC +=1.; NplusNminusMC[etabin + 12] += 1.; }
        }
        ThnContainerMC[0] = centbin;
        for(int idx = 0; idx < 24; idx++) ThnContainerMC[idx+1] = NplusNminusMC[idx];
        fThnCentPosNegEtaMC->Fill( ThnContainerMC );
        if( NplusMC > 2000 || NminusMC > 2000) fhEvent->Fill(4);
    }
  }

  // --- RECO LOOP ---
  double Nplus = 0, Nminus = 0;
  double NplusNminus[24];
  double ThnContainer[kdim];
  for(int idx = 0; idx < 24; idx++) NplusNminus[idx] = 0.;
  for(int jdx = 0; jdx < kdim; jdx++) ThnContainer[jdx] = 0.;

  int nch_event = 0; // Nch counter
  int trkIndx = -1;

  for (auto const& trk : *trks) {
      trkIndx++;
      if ( !trk.hasTrackDetails() ) continue;
      auto iter_tk = trk.pseudoTrack();

      double pterror = iter_tk.ptError();
      double vzErr=bestvzError;
      double vxErr=bestvxError;
      double dxy = iter_tk.dxy(bestvtx);
      double dz = iter_tk.dz(bestvtx);
      double dxysigma = sqrt(iter_tk.d0Error()*iter_tk.d0Error()+vxErr*bestvyError); 
      double dzsigma = sqrt(iter_tk.dzError()*iter_tk.dzError()+vzErr*vzErr);

      double pt = iter_tk.pt();
      double eta = iter_tk.eta();
      int charge = iter_tk.charge();
      double phi = iter_tk.phi();
      auto hit_pattern = iter_tk.hitPattern();
      double chi2ndof = ( double ) ( *chi2Map )[ trks->ptrAt( trkIndx ) ]; 
      double dcaxy = (dxy/dxysigma);
      double dcaz = (dz/dzsigma);
      double ptreso = (fabs(pterror) / pt);
      int nHits = iter_tk.numberOfValidHits();
      double chi2n = ( chi2ndof / hit_pattern.trackerLayersWithMeasurement() );
      int algo  = iter_tk.algo();

      if( charge == 0 ) continue;
      if( iter_tk.quality(reco::TrackBase::qualityByName("highPurity")) != 1 ) continue;
      if(pt > 10) { if( ptreso > 0.1) continue; } else { if( ptreso < 0.) continue; }
      if( fabs(dcaz) > 3.0 ) continue;
      if( fabs(dcaxy) > 3.0 ) continue;
      if( chi2n < 0.) continue;
      if( nHits < 0. ) continue;

      int count = 0;
      for(unsigned i = 0; i < algoParameters_.size(); i++) { if( algo == algoParameters_[i] ) count++; }
      if( count == 0 ) continue;

      if( pt <= 0.3 || pt >= 3.0 ) continue;
      if( eta <= -2.4 || eta >= 2.4 ) continue;

      nch_event++;

      float eff = 1.0, fak = 0.0;
      if(hEff3D && hFak3D) {
          eff = hEff3D->GetBinContent(hEff3D->GetXaxis()->FindBin(eta), hEff3D->GetYaxis()->FindBin(pt), hEff3D->GetZaxis()->FindBin(centralityBin));
          fak = hFak3D->GetBinContent(hFak3D->GetXaxis()->FindBin(eta), hFak3D->GetYaxis()->FindBin(pt), hFak3D->GetZaxis()->FindBin(centralityBin));
      }
      if(eff == 0) eff = 1.0; 
      float weight = (1. - fak)/eff;
      //float weight = 1.0

      hpt[centbin]->Fill(pt);
      heta[centbin]->Fill(eta);
      hphi[centbin]->Fill(phi);
      hnHits[centbin]->Fill(nHits);
      hptreso[centbin]->Fill(ptreso);
      hchi2[centbin]->Fill(chi2n);
      hDCAZ[centbin]->Fill(dcaz);
      hDCAXY[centbin]->Fill(dcaxy);

      hpt_eff[centbin]->Fill(pt, weight);
      heta_eff[centbin]->Fill(eta,weight);
      hphi_eff[centbin]->Fill(phi, weight);
      hnHits_eff[centbin]->Fill(nHits, weight);
      hptreso_eff[centbin]->Fill(ptreso, weight);
      hchi2_eff[centbin]->Fill(chi2n, weight);
      hDCAZ_eff[centbin]->Fill(dcaz, weight);
      hDCAXY_eff[centbin]->Fill(dcaxy, weight);
      
      fhPhi->Fill(phi);
      fheta->Fill(eta,weight);
      if( charge > 0 ) fhPtPlus->Fill(pt,weight);
      else fhPtMinus->Fill(pt,weight);

      Int_t etabin = fhEta->FindBin( eta ) - 1;
      if( etabin < 0 || etabin >= 12 ) continue; 

      if (charge > 0 ) { Nplus +=weight; NplusNminus[etabin] += weight; }
      else if(charge < 0 ) { Nminus +=weight; NplusNminus[etabin + 12] += weight; }
  } 

  // Fill Nch 
  hNch->Fill(nch_event);

  ThnContainer[0] = centbin;
  for(int idx = 0; idx < 24; idx++) ThnContainer[idx+1] = NplusNminus[idx];
  fThnCentPosNegEta->Fill( ThnContainer );
  if( Nplus > 2000 || Nminus > 2000) fhEvent->Fill(5);
}

DEFINE_FWK_MODULE(ChFluctuations);
