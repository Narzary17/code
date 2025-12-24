// -*- C++ -*-
// Author: Subash Behera
// Class: ChFluctuations
// Updated for O-O

#include "Analyzers/Cumulants/interface/ChFluctuations.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace edm;
using namespace std;

ChFluctuations::ChFluctuations(const edm::ParameterSet& iConfig) :
  // Tokens
  trackTags_(consumes< edm::View< pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("tracks"))),
  trackTagsgen_(consumes< edm::View< pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("tracksgen"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("srcTower"))),
  vtxTags_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex"))),
  centralityTags_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"))),

  noffmin_(iConfig.getUntrackedParameter<int>("noffmin")),
  noffmax_(iConfig.getUntrackedParameter<int>("noffmax")),
  ptnoffmin_(iConfig.getUntrackedParameter<double>("ptnoffmin")),
  ptnoffmax_(iConfig.getUntrackedParameter<double>("ptnoffmax")),
  dzdzerrornoff_(iConfig.getUntrackedParameter<double>("dzdzerrornoff")),
  d0d0errornoff_(iConfig.getUntrackedParameter<double>("d0d0errornoff")),
  pterrorptnoff_(iConfig.getUntrackedParameter<double>("pterrorptnoff")),

  etamin_(iConfig.getUntrackedParameter<double>("etamin")),
  etamax_(iConfig.getUntrackedParameter<double>("etamax")),
  ptmin_(iConfig.getUntrackedParameter<double>("ptmin")),
  ptmax_(iConfig.getUntrackedParameter<double>("ptmax")),

  fname_(iConfig.getUntrackedParameter<edm::InputTag>("fname")),
  cweight_(iConfig.getUntrackedParameter<bool>("cweight")),
  fIsMC(iConfig.getUntrackedParameter<bool>("IsMC")),

  // Bins
  centbins_(iConfig.getUntrackedParameter< std::vector < double > >("centbins")),
  etaBins_(iConfig.getUntrackedParameter< std::vector < double > >("etaBins"))
{
   TString filename(fname_.label().c_str());
   feff_ = 0x0;

   if(cweight_ && !filename.IsNull()) {
      edm::FileInPath fip(Form("Analyzers/Cumulants/data/EFF/%s",filename.Data()));
      feff_ = new TFile(fip.fullPath().c_str(),"READ");
      if(feff_) {
        hEff_3D = (TH3F*)feff_->Get("hEff_3D");
        hFak3D = (TH3F*)feff_->Get("hFak_3D");
        hMul_3D = (TH3F*)feff_->Get("hMul_3D");
        hSec_3D = (TH3F*)feff_->Get("hSec_3D");
      }
   }

   // --- CENTRALITY BINNING SETUP ---
   std::vector<double> defaultBins = {
        0,0.717572,1.43514,2.15272,2.87029,3.58786,4.30543,5.023,5.74057,6.45815,
        7.17572,7.89329,8.61086,9.32843,10.046,10.5615,10.9358,11.3222,11.7114,12.1054,
        12.5071,12.9114,13.3242,13.7529,14.1895,14.6323,15.0788,15.5384,16.0017,16.4718,
        16.954,17.4529,17.965,18.4818,19.0081,19.5384,20.0798,20.6295,21.1904,21.7624,
        22.3333,22.9144,23.5082,24.108,24.7103,25.3327,25.9643,26.6021,27.2404,27.8889,
        28.551,29.2201,29.8838,30.5609,31.2582,31.959,32.6755,33.3996,34.127,34.8736,
        35.624,36.3835,37.1513,37.9462,38.7306,39.5312,40.3522,41.1784,41.9996,42.8382,
        43.7007,44.5735,45.4614,46.3438,47.2291,48.1364,49.0614,49.9993,50.9485,51.9135,
        52.8921,53.8808,54.8964,55.893,56.9136,57.9568,59.0177,60.0888,61.1833,62.2824,
        63.399,64.5337,65.6807,66.8337,68.0214,69.219,70.4243,71.651,72.9014,74.165,
        75.4546,76.7537,78.0703,79.4045,80.7569,82.1322,83.5298,84.9421,86.3744,87.8423,
        89.3304,90.8548,92.361,93.9159,95.4992,97.0799,98.6655,100.303,101.96,103.647,
        105.32,107.036,108.77,110.55,112.332,114.116,115.936,117.785,119.708,121.683,
        123.605,125.589,127.588,129.583,131.654,133.746,135.887,138.016,140.191,142.401,
        144.623,146.867,149.14,151.514,153.909,156.277,158.74,161.209,163.662,166.177,
        168.728,171.298,173.901,176.575,179.302,182.038,184.793,187.576,190.47,193.355,
        196.279,199.296,202.355,205.406,208.567,211.746,215.007,218.296,221.624,224.981,
        228.464,231.972,235.5,239.147,242.785,246.538,250.308,254.2,258.246,262.353,
        266.531,270.819,275.228,279.794,284.51,289.356,294.368,299.561,305.114,310.888,
        316.984,323.451,330.516,337.981,346.265,355.527,366.193,379.001,395.659,421.298,1000
   };

   std::vector<double> userBins = iConfig.getUntrackedParameter<std::vector<double>>("binTable", defaultBins);

   for(size_t i = 0; i < userBins.size() && i <= 200; i++) {
       binBoundaries[i] = (float)userBins[i];
   }

   usesResource("TFileService");
}

ChFluctuations::~ChFluctuations() {
    if(feff_) feff_->Close();
}

void ChFluctuations::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void ChFluctuations::beginJob()
{
   edm::Service<TFileService> fs;
   TH1::SetDefaultSumw2();

   // QA Plots Directory
   TFileDirectory fVtxHist  = fs->mkdir("QAplots");

   hihf_pf = fVtxHist.make<TH1F>("hihf_pf", "Raw HF Energy", 1000, 0, 5000);
   hcent_bin = fVtxHist.make<TH1F>("hcent_bin", "Raw Centrality Index", 200, 0.0, 200.0);
   hcentbin  = fVtxHist.make<TH1F>("hcentbin", "Analysis Bins", centbins_.size()-1, &centbins_[0]);
   hNoff_ = fVtxHist.make<TH1F>("hNoff", "Number of Offline Tracks", 500, 0.0, 500.0);
   hNch   = fVtxHist.make<TH1F>("hNch", "Charged Particle Multiplicity", 500, 0.0, 500.0);
   hZBestVtx_ = fVtxHist.make<TH1F>("hZvtx", "", 600, -30.0, 30.0);

   hptbin = fVtxHist.make<TH1F>("hptbin", "", 500, 0., 50.);
   hetabin = fVtxHist.make<TH1F>("hetabin", "", etaBins_.size()-1, &etaBins_[0]);

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

   hpt_gen.resize(centbins_.size());
   heta_gen.resize(centbins_.size());
   hphi_gen.resize(centbins_.size());

   for(unsigned int i = 0; i < centbins_.size(); i++) {
      // Raw
      hpt[i]  = fVtxHist.make<TH1F>(Form("hpt_%d", i), "", 500, 0., 50.);
      heta[i] = fVtxHist.make<TH1F>(Form("heta_%d", i), "", etaBins_.size()-1, &etaBins_[0]);
      hphi[i] = fVtxHist.make<TH1F>(Form("hphi_%d", i), "", 62, -3.2, 3.2);
      hnHits[i]   = fVtxHist.make<TH1F>(Form("hnHits_%d", i), "", 50, 0., 50.);
      hptreso[i]  = fVtxHist.make<TH1F>(Form("hptreso_%d", i), "", 100, 0., 1.);
      hchi2[i]    = fVtxHist.make<TH1F>(Form("hchi2_%d", i), "", 100, 0, 10.);
      hDCAZ[i]    = fVtxHist.make<TH1F>(Form("hDCAZ_%d", i), "", 80, -4., 4.);
      hDCAXY[i]   = fVtxHist.make<TH1F>(Form("hDCAXY_%d", i), "", 80, -4., 4.);

      // Corrected
      hpt_eff[i] = fVtxHist.make<TH1F>(Form("hpt_eff_%d", i), "", 500, 0., 50.);
      heta_eff[i] = fVtxHist.make<TH1F>(Form("heta_eff_%d", i), "", etaBins_.size()-1, &etaBins_[0]);
      hphi_eff[i] = fVtxHist.make<TH1F>(Form("hphi_eff_%d", i), "", 62, -3.2, 3.2);
      hnHits_eff[i]   = fVtxHist.make<TH1F>(Form("hnHits_eff_%d", i), "", 50, 0., 50.);
      hptreso_eff[i]  = fVtxHist.make<TH1F>(Form("hptreso_eff_%d", i), "", 100, 0., 1.);
      hchi2_eff[i]    = fVtxHist.make<TH1F>(Form("hchi2_eff_%d", i), "", 100, 0, 10.);
      hDCAZ_eff[i]    = fVtxHist.make<TH1F>(Form("hDCAZ_eff_%d", i), "", 80, -4., 4.);
      hDCAXY_eff[i]   = fVtxHist.make<TH1F>(Form("hDCAXY_eff_%d", i), "", 80, -4., 4.);

      // Gen
      hpt_gen[i]  = fVtxHist.make<TH1F>(Form("hpt_gen_%d", i), "", 500, 0., 50.);
      heta_gen[i] = fVtxHist.make<TH1F>(Form("heta_gen_%d", i), "", etaBins_.size()-1, &etaBins_[0]);
      hphi_gen[i] = fVtxHist.make<TH1F>(Form("hphi_gen_%d", i), "", 62, -3.2, 3.2);
   }

   // THnSparse Setup
   const int kdim = 25;
   int bin[kdim];
   double min[kdim], max[kdim];
   bin[0] = 81;
   min[0] = -0.5; max[0] = 80.5;
   for(int ibin = 1; ibin < kdim; ibin++) { bin[ibin] = 3000; min[ibin] = -0.5; max[ibin] = 2999.5; }

   // Global Directory Histograms
   TFileDirectory fGlobalHist  = fs->mkdir("Global");

   Double_t CentBinning[10] = { 0, 10, 20, 40, 60, 80, 100, 120, 140, 160 };
   Double_t EtaBinning[13] = { 0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4 };

   fhEvent = fGlobalHist.make<TH1D>("fhEvent", "Event counter", 8, -0.5, 7.5 );
   fhCent = fGlobalHist.make<TH1D>("fhCent", "centrality", 9, CentBinning );
   fhEta = fGlobalHist.make<TH1D>("fhEta", "EtaBining", 12,  EtaBinning);
   fhVz = fGlobalHist.make<TH1D>("fhVz", "V_{z} distribution", 320, -16.0, 16.0);
   fhPhi = fGlobalHist.make<TH1D>("fhPhi", "#phi distributions", 640, -3.2, 3.2);
   fhPtPlus = fGlobalHist.make<TH1D>("fhPtPlus", "N^{+} p_{T} distributions", 120, 0., 6.0);
   fhPtMinus = fGlobalHist.make<TH1D>("fhPtMinus", "N^{-} p_{T} distributions", 120, 0., 6.0);
   fheta = fGlobalHist.make<TH1D>("fheta", "#eta distribution", 100, -2.5, 2.5);
   
   fThnCentPosNegEta = fGlobalHist.make<THnSparseI>("fThnCentPosNegEta", "Reco Cent-Etabinwise-Nplus-Nminus", kdim, bin, min, max);
   fThnCentPosNegEtaMC = fGlobalHist.make<THnSparseI>("fThnCentPosNegEtaMC", "Gen Cent-Etabinwise-Nplus-Nminus", kdim, bin, min, max);
}

void ChFluctuations::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  fhEvent->Fill(0);

  // ----------------- MANUAL CENTRALITY CALCULATION  ----------------
  double etHFtowerSumPlus = 0;
  double etHFtowerSumMinus = 0;
  double etHFtowerSum = 0;

  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByToken(pfToken_, pfs);

  if(pfs.isValid()){
      for(unsigned int i = 0, n = pfs->size(); i < n; ++i) {
        const pat::PackedCandidate &pf = (*pfs)[i];
        if(pf.pdgId() == 1 || pf.pdgId() == 2){
          if(pf.eta()>3.0 && pf.eta()<5.2) etHFtowerSumPlus += pf.et();
          if(pf.eta()<-3.0 && pf.eta()>-5.2) etHFtowerSumMinus += pf.et();
        }
      }
  }
  etHFtowerSum = etHFtowerSumPlus + etHFtowerSumMinus;
  hihf_pf->Fill(etHFtowerSum);

  int cent = -1;
  for (int ic=0; ic < nCentTableBins; ic++){
    if(etHFtowerSum > binBoundaries[ic] && etHFtowerSum < binBoundaries[ic+1]){
      cent = nCentTableBins - 1 - ic;
    }
  }

  if(cent < 0 || cent >= 200) return;

  auto pvs = iEvent.getHandle( vtxTags_ );
  double bestvzError;
  math::XYZPoint bestvtx;
  math::Error<3>::type vtx_cov;

  if ( pvs.isValid() && !pvs->empty()) {
      const reco::Vertex& vtx = (*pvs)[0];
      bestvzError = vtx.zError();
      bestvtx = vtx.position();
      vtx_cov = vtx.covariance();
  } else {
      return;
  }

  hZBestVtx_->Fill(bestvtx.z());
  fhVz->Fill(bestvtx.z());
  fhEvent->Fill(2);

  // ----------------- Noff CALCULATION ----------------
  auto trks = iEvent.getHandle( trackTags_ );
  if(!trks.isValid()) return;

  int noff_ = 0;

  for (auto const& trk : *trks) {
       if ( !trk.hasTrackDetails() ) continue;
       auto iter_tk = trk.pseudoTrack();

       double dzvtx = iter_tk.dz( bestvtx );
       double dxyvtx = iter_tk.dxy( bestvtx );
       double dzerror = std::hypot( iter_tk.dzError(), bestvzError );
       double dxyerror = iter_tk.dxyError( bestvtx, vtx_cov );
       double pterror = iter_tk.ptError();

       double pt     = iter_tk.pt();
       double eta    = iter_tk.eta();
       double charge = iter_tk.charge();

       if( !iter_tk.quality(reco::TrackBase::highPurity) ) continue;
       if( fabs(pterror) / pt      > pterrorptnoff_ ) continue;
       if( fabs(dzvtx / dzerror)   > dzdzerrornoff_ ) continue;
       if( fabs(dxyvtx / dxyerror) > d0d0errornoff_ ) continue;
       if( pt < 0.0001 ) continue;
       if( charge == 0 ) continue;

       if(eta < -2.4 || eta > 2.4) continue;
       if(pt < ptnoffmin_) continue;

       noff_++;
  }
 // hNoff_->Fill(noff_);

  if( noff_ < noffmin_ || noff_ >= noffmax_ ) return;

  // ----------------- ANALYSIS LOOP ----------------
  int centbin = -1;
  for(unsigned int i=0; i < centbins_.size()-1; i++) {
      if((float)cent >= centbins_[i] && (float)cent < centbins_[i+1]) {
          centbin = i;
          break;
      }
  }

  fhEvent->Fill(3);
  hcent_bin->Fill(cent);
  hcentbin->Fill(cent);
  fhCent->Fill(cent);

  if( centbin < 0 ) return;

  // GEN LOOP
  if(!fIsMC) {
      auto trksgen = iEvent.getHandle( trackTagsgen_ );
      double NplusNminusMC[24];
      for(int k=0; k<24; k++) NplusNminusMC[k] = 0;

      if(trksgen.isValid()) {
         for (auto const& iter_tk_gn : *trksgen) {
             if(iter_tk_gn.status() != 1) continue;
             double gen_pt = iter_tk_gn.pt();
             double gen_eta = iter_tk_gn.eta();
             int gen_charge = iter_tk_gn.charge();

             if(gen_charge == 0) continue;
             if(gen_pt < ptmin_ || gen_pt > ptmax_) continue;
             if(gen_eta < etamin_ || gen_eta > etamax_) continue;

             hpt_gen[centbin]->Fill(gen_pt);
             heta_gen[centbin]->Fill(gen_eta);
             hphi_gen[centbin]->Fill(iter_tk_gn.phi());

             int etabin = -1;
             for(size_t ie = 0; ie < etaBins_.size()-1; ie++){
                 if(gen_eta >= etaBins_[ie] && gen_eta < etaBins_[ie+1]) { etabin = ie; break; }
             }

             if(etabin >= 0 && etabin < 12) {
                 if(gen_charge > 0) NplusNminusMC[etabin] += 1.0;
                 else               NplusNminusMC[etabin + 12] += 1.0;
             }
         }
         double thn_content_mc[25];
         thn_content_mc[0] = centbin;
         for(int k=0; k<24; k++) thn_content_mc[k+1] = NplusNminusMC[k];
         fThnCentPosNegEtaMC->Fill(thn_content_mc);
      }
  }

  // RECO LOOP
  double NplusNminus[24];
  for(int k=0; k<24; k++) NplusNminus[k] = 0;
  int nch_counter = 0;

  for (auto const& trk : *trks) {
       if ( !trk.hasTrackDetails() ) continue;
       auto iter_tk = trk.pseudoTrack();

       double dzvtx = iter_tk.dz( bestvtx );
       double dxyvtx = iter_tk.dxy( bestvtx );
       double dzerror = std::hypot( iter_tk.dzError(), bestvzError );
       double dxyerror = iter_tk.dxyError( bestvtx, vtx_cov );
       double pterror = iter_tk.ptError();

       double eta    = iter_tk.eta();
       double pt     = iter_tk.pt();
       double phi    = iter_tk.phi();
       int charge = iter_tk.charge();

       int nHits = iter_tk.numberOfValidHits();
       double ptreso = (fabs(pterror) / pt);
       double chi2n = iter_tk.normalizedChi2() / iter_tk.hitPattern().trackerLayersWithMeasurement();
       double dcaz = (dzerror > 0) ? dzvtx / dzerror : -99;
       double dcaxy = (dxyerror > 0) ? dxyvtx / dxyerror : -99;

       if( !iter_tk.quality(reco::TrackBase::highPurity) ) continue;
       if( fabs(dcaz)   > 3.0 ) continue;
       if( fabs(dcaxy) > 3.0 ) continue;
       if( charge == 0 ) continue;

       if(eta < etamin_ || eta > etamax_) continue;
       if(pt < ptmin_ || pt > ptmax_)     continue;

       nch_counter++;

       double weight = 1.0;
       if(hEff_3D && cweight_) {
          double eff = hEff_3D->GetBinContent(hEff_3D->GetXaxis()->FindBin(eta), hEff_3D->GetYaxis()->FindBin(pt), hEff_3D->GetZaxis()->FindBin(cent));
           double fake = hFak3D->GetBinContent(hFak3D->GetXaxis()->FindBin(eta), hFak3D->GetYaxis()->FindBin(pt), hFak3D->GetZaxis()->FindBin(cent));

           if(eff > 0.001) weight = (1 - fake)/eff;
       }

       // Fill Raw QA
       hpt[centbin]->Fill(pt);
       heta[centbin]->Fill(eta);
       hphi[centbin]->Fill(phi);
       hnHits[centbin]->Fill(nHits);
       hptreso[centbin]->Fill(ptreso);
       hchi2[centbin]->Fill(chi2n);
       hDCAZ[centbin]->Fill(dcaz);
       hDCAXY[centbin]->Fill(dcaxy);

       // Fill Corrected QA
       hpt_eff[centbin]->Fill(pt, weight);
       heta_eff[centbin]->Fill(eta, weight);
       hphi_eff[centbin]->Fill(phi, weight);
       hnHits_eff[centbin]->Fill(nHits, weight);
       hptreso_eff[centbin]->Fill(ptreso, weight);
       hchi2_eff[centbin]->Fill(chi2n, weight);
       hDCAZ_eff[centbin]->Fill(dcaz, weight);
       hDCAXY_eff[centbin]->Fill(dcaxy, weight);

       // Fill Global
       fhPhi->Fill(phi);
       fheta->Fill(eta, weight);
       if( charge > 0 ) fhPtPlus->Fill(pt, weight);
       else fhPtMinus->Fill(pt, weight);

       int etabin = -1;
       for(size_t ie = 0; ie < etaBins_.size()-1; ie++){
           if(eta >= etaBins_[ie] && eta < etaBins_[ie+1]) {
               etabin = ie;
               break;
           }
       }

       if(etabin >= 0 && etabin < 12) {
           if(charge > 0) NplusNminus[etabin] += weight;
           else           NplusNminus[etabin + 12] += weight;
       }
   }
   hNch->Fill(nch_counter);

   double thn_content[25];
   thn_content[0] = centbin;
   for(int k=0; k<24; k++) thn_content[k+1] = NplusNminus[k];
   fThnCentPosNegEta->Fill(thn_content);
}

void ChFluctuations::endJob() {}

DEFINE_FWK_MODULE(ChFluctuations);
