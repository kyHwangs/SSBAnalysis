#ifndef _ssb_analysis_

#define _ssb_analysis_

#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TEnv.h"
#include "TLorentzVector.h"
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>

#include "./../analysis/SSBTree.h"

// Textreader
#include "./../TextReader/TextReader.hpp"

#include "./../interface/BTagCalibrationStandalone.h"
#include "./../interface/EffTable.hpp"
#include "./../interface/LumiReWeighting.h"
#include "./../interface/ssb_cpviol.hpp"
#include "./../interface/ssb_eff.hpp"
#include "./../kinsol/TtFullLepKinSolver.hpp"

#include "./../RoccoR/RoccoR.h"

#include "./../KinSolv/KinematicReconstruction.h"
#include "./../KinSolv/KinematicReconstructionSolution.h"
#include "./../KinSolv/analysisUtils.h"

using namespace std;

class ssb_analysis : public SSBTree {
public:
  // declare functions
  ssb_analysis(TTree *tree = 0, 
               bool fIsData_ = 0,
               std::string fEra_ = 0,
               std::string fProcessName_ = 0, 
               bool fRoccor_enalbed_ = 0,
               bool fIDISO_enalbed_ = 0, 
               bool fTRIGG_enalbed_ = 0,
               float fMassCut_ = 0);

  virtual ~ssb_analysis();

  // basic frame
  virtual void Loop(char *logfile);
  void Start(int genLoopon);
  void End();
  double GetLHEmass();

  // user define functions
  void SetOutputFileName(char *outname);
  void DeclareHistos();

  bool GetIsData() { return fIsData; }
  std::string GetProcessName() { return fProcessName; }

  double GetNormalization();

  double GetPUweight(int npv) { return fPuReweighting->weight(npv); }
  double GetPUweight(float npv) { return fPuReweighting->weight(npv); }

  // double GetEventTriggerEfficiencyScaleFactor(double pt_1, double eta_1, double pt_2, double eta_2) const;

  // double GetIdisoSF(double pT, double eta);
  // double GetTriggerSF(double pT1, double eta1, double pt2, double eta2);


  // TBranch        *b_Muon_isPF;
  // TBranch        *b_Muon_isGlobal;
  // TBranch        *b_Muon_isTracker;

  // TBranch        *b_Muon_ValidTrackHitFraction;
  // TBranch        *b_Muon_GlobalTrackChi2;
  // TBranch        *b_Muon_TrackerSTAPositionMatch;
  // TBranch        *b_Muon_TrackKink;
  // TBranch        *b_Muon_TrackExist;
  // TBranch        *b_Muon_SegmentCompatibility;

  // TBranch        *b_Muon_GlobalMuonTrackChamberHit;
  // TBranch        *b_Muon_TunePGlobalMuonTrackChamberHit;
  // TBranch        *b_Muon_TunePpTError;
  // TBranch        *b_Muon_TunePpT;
  // TBranch        *b_Muon_StationsHasSegments;

  // TBranch        *b_Muon_dxy;
  // TBranch        *b_Muon_dB;
  // TBranch        *b_Muon_dz;
  // TBranch        *b_Muon_PixelHit;
  // TBranch        *b_Muon_TrackLayerWithHit;
  
  // ssbtree->Branch("Muon_isPF", &VectorBox_Bool["Muon_isPF"]);
  // ssbtree->Branch("Muon_isGlobal", &VectorBox_Bool["Muon_isGlobal"]);
  // ssbtree->Branch("Muon_isTracker", &VectorBox_Bool["Muon_isTracker"]);

  // ssbtree->Branch("Muon_ValidTrackHitFraction", &VectorBox_Double["Muon_ValidTrackHitFraction"]);
  // ssbtree->Branch("Muon_GlobalTrackChi2", &VectorBox_Double["Muon_GlobalTrackChi2"]);
  // ssbtree->Branch("Muon_TrackerSTAPositionMatch", &VectorBox_Double["Muon_TrackerSTAPositionMatch"]);
  // ssbtree->Branch("Muon_TrackKink", &VectorBox_Double["Muon_TrackKink"]);
  // ssbtree->Branch("Muon_TrackExist", &VectorBox_Bool["Muon_TrackExist"]);
  // ssbtree->Branch("Muon_SegmentCompatibility", &VectorBox_Double["Muon_SegmentCompatibility"]);

  // ssbtree->Branch("Muon_GlobalMuonTrackChamberHit", &VectorBox_Double["Muon_GlobalMuonTrackChamberHit"]);
  // ssbtree->Branch("Muon_TunePGlobalMuonTrackChamberHit", &VectorBox_Double["Muon_TunePGlobalMuonTrackChamberHit"]);
  // ssbtree->Branch("Muon_TunePpTError", &VectorBox_Double["Muon_TunePpTError"]);
  // ssbtree->Branch("Muon_TunePpT", &VectorBox_Double["Muon_TunePpT"]);
  // ssbtree->Branch("Muon_StationsHasSegments", &VectorBox_Double["Muon_StationsHasSegments"]);

  // ssbtree->Branch("Muon_dxy", &VectorBox_Double["Muon_dxy"]);
  // ssbtree->Branch("Muon_dB", &VectorBox_Double["Muon_dB"]);
  // ssbtree->Branch("Muon_dz", &VectorBox_Double["Muon_dz"]);
  // ssbtree->Branch("Muon_PixelHit", &VectorBox_Double["Muon_PixelHit"]);
  // ssbtree->Branch("Muon_TrackLayerWithHit", &VectorBox_Double["Muon_TrackLayerWithHit"]);

  struct stdMUON {
    TLorentzVector v;
    TLorentzVector _v;
    int charge;
    bool id;
    float iso;

    bool isPF;
    bool isGlobal;
    bool isTracker;

    double ValidTrackHitFraction;
    double GlobalTrackChi2;
    double TrackerSTAPositionMatch;
    double TrackKink;
    bool TrackExist;

    double SegmentCompatibility;
    double GlobalMuonTrackChamberHit;
    double TunePGlobalMuonTrackChamberHit;
    double TunePpTError;
    double TunePpT;

    double StationsHasSegments;
    double dxy;
    double dB;
    double dz;
    double PixelHit;
    double TrackLayerWithHit;

    int nTrackLayer;

    stdMUON(TLorentzVector v_, int charge_, bool id_, float iso_,
            bool isPF_,
            bool isGlobal_,
            bool isTracker_,
            double ValidTrackHitFraction_,
            double GlobalTrackChi2_,
            double TrackerSTAPositionMatch_,
            double TrackKink_,
            bool TrackExist_,
            double SegmentCompatibility_,
            double GlobalMuonTrackChamberHit_,
            double TunePGlobalMuonTrackChamberHit_,
            double TunePpTError_,
            double TunePpT_,
            double StationsHasSegments_,
            double dxy_,
            double dB_,
            double dz_,
            double PixelHit_,
            double TrackLayerWithHit_,
            int nTrackLayer_ = 0)
        : v(v_), _v(v_), charge(charge_), id(id_), iso(iso_),
          isPF(isPF_),
          isGlobal(isGlobal_),
          isTracker(isTracker_),
          ValidTrackHitFraction(ValidTrackHitFraction_),
          GlobalTrackChi2(GlobalTrackChi2_),
          TrackerSTAPositionMatch(TrackerSTAPositionMatch_),
          TrackKink(TrackKink_),
          TrackExist(TrackExist_),
          SegmentCompatibility(SegmentCompatibility_),
          GlobalMuonTrackChamberHit(GlobalMuonTrackChamberHit_),
          TunePGlobalMuonTrackChamberHit(TunePGlobalMuonTrackChamberHit_),
          TunePpTError(TunePpTError_),
          TunePpT(TunePpT_),
          StationsHasSegments(StationsHasSegments_),
          dxy(dxy_),
          dB(dB_),
          dz(dz_),
          PixelHit(PixelHit_),
          TrackLayerWithHit(TrackLayerWithHit_),
          nTrackLayer(nTrackLayer_)
          {}

    bool isPass() {
      return (v.Pt() > 10. && std::abs(v.Eta()) < 2.4 && id && iso < 0.15);
    }
  };

  struct stdELEC {
    TLorentzVector v;
    float sc_eta;
    bool id;
    float iso;

    stdELEC(TLorentzVector v_, float sc_eta_, bool id_, float iso_)
        : v(v_), sc_eta(sc_eta_), id(id_), iso(iso_) {}

    bool isPass() {
      return (v.Pt() > 10. && std::abs(v.Eta()) < 2.4 && id && iso < 0.0821 &&
              !(std::abs(sc_eta) > 1.4442 && std::abs(sc_eta) < 1.5660));
    }
  };

  struct stdJET {
    TLorentzVector v;
    int id;
    float bdisc;

    stdJET(TLorentzVector v_, int id_, float bdisc_)
        : v(v_), id(id_), bdisc(bdisc_) {}

    bool isPass() {
      return (v.Pt() > 30. && std::abs(v.Eta()) < 2.4 && id > 0);
    }

    bool isBJet() { return (bdisc > 0.5426); }

    bool isOverlap(std::vector<TLorentzVector> lep) {

      for (int j = 0; j < lep.size(); j++)
        if (v.DeltaR(lep.at(j)) < 0.4)
          return true;

      return false;
    }
  };

private:
  // put variables that you want
  char *outfile;
  TFile *fout;

  bool fIsData;
  std::string fEra;
  std::string fProcessName;
  bool fRoccor_enalbed;
  bool fIDISO_enalbed;
  bool fTRIGG_enalbed;

  double fLeadingMuPt;
  double fSubLeadingMuPt;
  float fMassCut;

  // vector for ChargeMisId
  edm::LumiReWeighting *fPuReweighting;
  RoccoR *fRoccoR;
  EffTable fID_SF;
  EffTable fISO_SF;
  EffTable fTRIG_DataEff;
  EffTable fTRIG_MCEff;

public:
  // declare histograms

  TH1D *h_LHEDimuonMass;
  TH1D *h_LHEnMuon;
  TH1D *h_GenWeight;

  TH1D *h_JetPt;
  TH1D *h_JetEta;
  TH1D *h_JetPhi;

  TH1D *h_BJetPt;
  TH1D *h_BJetEta;
  TH1D *h_BJetPhi;

  TH1D *h_nJet_before;
  TH1D *h_nJet_after;
  TH1D *h_nBJet_before;
  TH1D *h_nBJet_after;

  TH1D *h_LeadingMuonPt;
  TH1D *h_LeadingMuonEta;
  TH1D *h_LeadingMuonPhi;

  TH1D *h_SubleadingMuonPt;
  TH1D *h_SubleadingMuonEta;
  TH1D *h_SubleadingMuonPhi;

  TH1D *h_MuonPt;
  TH1D *h_MuonEta;
  TH1D *h_MuonPhi;

  TH1D *h_dimuonMass;
  TH1D *h_dimuonMass_wide;
  TH1D *h_dimuonPt;
  TH1D *h_dimuonRap;

  TH1D *h_PV_Count_before_corr;
  TH1D *h_PileUp_Count_Interaction_before_corr;
  TH1D *h_PileUp_Count_Intime_before_corr;

  TH1D *h_PV_Count_after_corr;
  TH1D *h_PileUp_Count_Interaction_after_corr;
  TH1D *h_PileUp_Count_Intime_after_corr;

  TH1D *h_isGlobal_BeforeSel;
  TH1D *h_Chi2GlobalMuonTrack_BeforeSel;
  TH1D *h_NChamberGlobalMuonTrack_BeforeSel;
  TH1D *h_NofMatchedStation_BeforeSel;
  TH1D *h_dxy_BeforeSel;
  TH1D *h_dB_BeforeSel;
  TH1D *h_dxyOrdB_BeforeSel;
  TH1D *h_dz_BeforeSel;
  TH1D *h_NPixelHit_BeforeSel;
  TH1D *h_ValidInnerTrackHit_BeforeSel;

  TH1D *h_isGlobal_AfterSel;
  TH1D *h_Chi2GlobalMuonTrack_AfterSel;
  TH1D *h_NChamberGlobalMuonTrack_AfterSel;
  TH1D *h_NofMatchedStation_AfterSel;
  TH1D *h_dxy_AfterSel;
  TH1D *h_dB_AfterSel;
  TH1D *h_dxyOrdB_AfterSel;
  TH1D *h_dz_AfterSel;
  TH1D *h_NPixelHit_AfterSel;
  TH1D *h_ValidInnerTrackHit_AfterSel;

  TH1D *h_JetPt_0J0B;
  TH1D *h_JetEta_0J0B;
  TH1D *h_JetPhi_0J0B;
  TH1D *h_BJetPt_0J0B;
  TH1D *h_BJetEta_0J0B;
  TH1D *h_BJetPhi_0J0B;

  TH1D *h_JetPt_0J;
  TH1D *h_JetEta_0J;
  TH1D *h_JetPhi_0J;
  TH1D *h_BJetPt_0J;
  TH1D *h_BJetEta_0J;
  TH1D *h_BJetPhi_0J;

  TH1D *h_JetPt_1J;
  TH1D *h_JetEta_1J;
  TH1D *h_JetPhi_1J;
  TH1D *h_BJetPt_1J;
  TH1D *h_BJetEta_1J;
  TH1D *h_BJetPhi_1J;

  TH1D *h_JetPt_mt1J;
  TH1D *h_JetEta_mt1J;
  TH1D *h_JetPhi_mt1J;
  TH1D *h_BJetPt_mt1J;
  TH1D *h_BJetEta_mt1J;
  TH1D *h_BJetPhi_mt1J;

  TH1D *h_JetPt_0BJ;
  TH1D *h_JetEta_0BJ;
  TH1D *h_JetPhi_0BJ;
  TH1D *h_BJetPt_0BJ;
  TH1D *h_BJetEta_0BJ;
  TH1D *h_BJetPhi_0BJ;

  TH1D *h_JetPt_1BJ;
  TH1D *h_JetEta_1BJ;
  TH1D *h_JetPhi_1BJ;
  TH1D *h_BJetPt_1BJ;
  TH1D *h_BJetEta_1BJ;
  TH1D *h_BJetPhi_1BJ;

  TH1D *h_JetPt_mt1BJ;
  TH1D *h_JetEta_mt1BJ;
  TH1D *h_JetPhi_mt1BJ;
  TH1D *h_BJetPt_mt1BJ;
  TH1D *h_BJetEta_mt1BJ;
  TH1D *h_BJetPhi_mt1BJ;

  TH1D *h_LeadingMuonPt_0J0B;
  TH1D *h_LeadingMuonEta_0J0B;
  TH1D *h_LeadingMuonPhi_0J0B;
  TH1D *h_SubleadingMuonPt_0J0B;
  TH1D *h_SubleadingMuonEta_0J0B;
  TH1D *h_SubleadingMuonPhi_0J0B;
  TH1D *h_MuonPt_0J0B;
  TH1D *h_MuonEta_0J0B;
  TH1D *h_MuonPhi_0J0B;
  TH1D *h_dimuonMass_0J0B;
  TH1D *h_dimuonMass_wide_0J0B;
  TH1D *h_dimuonPt_0J0B;
  TH1D *h_dimuonRap_0J0B;

  TH1D *h_LeadingMuonPt_0J;
  TH1D *h_LeadingMuonEta_0J;
  TH1D *h_LeadingMuonPhi_0J;
  TH1D *h_SubleadingMuonPt_0J;
  TH1D *h_SubleadingMuonEta_0J;
  TH1D *h_SubleadingMuonPhi_0J;
  TH1D *h_MuonPt_0J;
  TH1D *h_MuonEta_0J;
  TH1D *h_MuonPhi_0J;
  TH1D *h_dimuonMass_0J;
  TH1D *h_dimuonMass_wide_0J;
  TH1D *h_dimuonPt_0J;
  TH1D *h_dimuonRap_0J;

  TH1D *h_LeadingMuonPt_1J;
  TH1D *h_LeadingMuonEta_1J;
  TH1D *h_LeadingMuonPhi_1J;
  TH1D *h_SubleadingMuonPt_1J;
  TH1D *h_SubleadingMuonEta_1J;
  TH1D *h_SubleadingMuonPhi_1J;
  TH1D *h_MuonPt_1J;
  TH1D *h_MuonEta_1J;
  TH1D *h_MuonPhi_1J;
  TH1D *h_dimuonMass_1J;
  TH1D *h_dimuonMass_wide_1J;
  TH1D *h_dimuonPt_1J;
  TH1D *h_dimuonRap_1J;

  TH1D *h_LeadingMuonPt_mt1J;
  TH1D *h_LeadingMuonEta_mt1J;
  TH1D *h_LeadingMuonPhi_mt1J;
  TH1D *h_SubleadingMuonPt_mt1J;
  TH1D *h_SubleadingMuonEta_mt1J;
  TH1D *h_SubleadingMuonPhi_mt1J;
  TH1D *h_MuonPt_mt1J;
  TH1D *h_MuonEta_mt1J;
  TH1D *h_MuonPhi_mt1J;
  TH1D *h_dimuonMass_mt1J;
  TH1D *h_dimuonMass_wide_mt1J;
  TH1D *h_dimuonPt_mt1J;
  TH1D *h_dimuonRap_mt1J;

  TH1D *h_LeadingMuonPt_0BJ;
  TH1D *h_LeadingMuonEta_0BJ;
  TH1D *h_LeadingMuonPhi_0BJ;
  TH1D *h_SubleadingMuonPt_0BJ;
  TH1D *h_SubleadingMuonEta_0BJ;
  TH1D *h_SubleadingMuonPhi_0BJ;
  TH1D *h_MuonPt_0BJ;
  TH1D *h_MuonEta_0BJ;
  TH1D *h_MuonPhi_0BJ;
  TH1D *h_dimuonMass_0BJ;
  TH1D *h_dimuonMass_wide_0BJ;
  TH1D *h_dimuonPt_0BJ;
  TH1D *h_dimuonRap_0BJ;

  TH1D *h_LeadingMuonPt_1BJ;
  TH1D *h_LeadingMuonEta_1BJ;
  TH1D *h_LeadingMuonPhi_1BJ;
  TH1D *h_SubleadingMuonPt_1BJ;
  TH1D *h_SubleadingMuonEta_1BJ;
  TH1D *h_SubleadingMuonPhi_1BJ;
  TH1D *h_MuonPt_1BJ;
  TH1D *h_MuonEta_1BJ;
  TH1D *h_MuonPhi_1BJ;
  TH1D *h_dimuonMass_1BJ;
  TH1D *h_dimuonMass_wide_1BJ;
  TH1D *h_dimuonPt_1BJ;
  TH1D *h_dimuonRap_1BJ;

  TH1D *h_LeadingMuonPt_mt1BJ;
  TH1D *h_LeadingMuonEta_mt1BJ;
  TH1D *h_LeadingMuonPhi_mt1BJ;
  TH1D *h_SubleadingMuonPt_mt1BJ;
  TH1D *h_SubleadingMuonEta_mt1BJ;
  TH1D *h_SubleadingMuonPhi_mt1BJ;
  TH1D *h_MuonPt_mt1BJ;
  TH1D *h_MuonEta_mt1BJ;
  TH1D *h_MuonPhi_mt1BJ;
  TH1D *h_dimuonMass_mt1BJ;
  TH1D *h_dimuonMass_wide_mt1BJ;
  TH1D *h_dimuonPt_mt1BJ;
  TH1D *h_dimuonRap_mt1BJ;
  
};
#endif

#ifdef ssb_analysis_cxx

ssb_analysis::ssb_analysis(TTree *tree, 
                           bool fIsData_,
                           std::string fEra_,
                           std::string fProcessName_, 
                           bool fRoccor_enalbed_,
                           bool fIDISO_enalbed_, 
                           bool fTRIGG_enalbed_,
                           float fMassCut_)
    : fIsData(fIsData_), 
      fEra(fEra_),
      fProcessName(fProcessName_),
      fRoccor_enalbed(fRoccor_enalbed_), 
      fIDISO_enalbed(fIDISO_enalbed_),
      fTRIGG_enalbed(fTRIGG_enalbed_), 
      fMassCut(fMassCut_) 
{
  if (tree == 0) {
    printf("ERROR: Can't find any input tree.\n");
  
  }
  Init(tree);

  // initializing HERE !

  // EffTable fID_SF;
  // EffTable fISO_SF;
  // EffTable fTRIG_DataEff;
  // EffTable fTRIG_MCEff;

  if (fEra == "UL2016APV") {
    fRoccoR = new RoccoR("./RoccoR/RoccoR2016aUL.txt");
    fPuReweighting = new edm::LumiReWeighting("./pileuInfo/MC_2016.root", "./pileuInfo/PileupHistogram-goldenJSON-13tev-2016-preVFP-69200ub-99bins.root", "pileup", "pileup");
    std::cout << "MC PU : ./pileuInfo/MC_2016.root" << std::endl;
    std::cout << "DATA PU : ./pileuInfo/PileupHistogram-goldenJSON-13tev-2016-preVFP-69200ub-99bins.root" << std::endl;
    std::cout << "Roccor : ./RoccoR/RoccoR2016aUL.txt" << std::endl;

    fLeadingMuPt = 26.;
    fSubLeadingMuPt = 10.;
    fID_SF        = EffTable("./lepEff/eff_240728/Run2016_UL_HIPM_ID.txt");
    fISO_SF       = EffTable("./lepEff/eff_240728/Run2016_UL_HIPM_ISO.txt");
    fTRIG_DataEff = EffTable("./lepEff/eff_240728/Run2016_UL_HIPM_TRIG_DATAeff.txt");
    fTRIG_MCEff   = EffTable("./lepEff/eff_240728/Run2016_UL_HIPM_TRIG_MCeff.txt");

  } else if (fEra == "UL2016") {
    fRoccoR = new RoccoR("./RoccoR/RoccoR2016bUL.txt");
    fPuReweighting = new edm::LumiReWeighting("./pileuInfo/MC_2016.root", "./pileuInfo/PileupHistogram-goldenJSON-13tev-2016-postVFP-69200ub-99bins.root", "pileup", "pileup");
    std::cout << "MC PU : ./pileuInfo/MC_2016.root" << std::endl;
    std::cout << "DATA PU : ./pileuInfo/PileupHistogram-goldenJSON-13tev-2016-postVFP-69200ub-99bins.root" << std::endl;  
    std::cout << "Roccor : ./RoccoR/RoccoR2016bUL.txt" << std::endl;

    fLeadingMuPt = 26.;
    fSubLeadingMuPt = 10.;
    fID_SF        = EffTable("./lepEff/eff_240728/Run2016_UL_ID.txt");
    fISO_SF       = EffTable("./lepEff/eff_240728/Run2016_UL_ISO.txt");
    fTRIG_DataEff = EffTable("./lepEff/eff_240728/Run2016_UL_TRIG_DATAeff.txt");
    fTRIG_MCEff   = EffTable("./lepEff/eff_240728/Run2016_UL_TRIG_MCeff.txt");

  } else if (fEra == "UL2017") {
    fRoccoR = new RoccoR("./RoccoR/RoccoR2017UL.txt");    
    fPuReweighting = new edm::LumiReWeighting("./pileuInfo/MC_2017.root", "./pileuInfo/PileupHistogram-goldenJSON-13tev-2017-69200ub-99bins.root", "pileup", "pileup");
    std::cout << "MC PU : ./pileuInfo/MC_2017.root" << std::endl;
    std::cout << "DATA PU : ./pileuInfo/PileupHistogram-goldenJSON-13tev-2017-69200ub-99bins.root" << std::endl;  
    std::cout << "Roccor : ./RoccoR/RoccoR2017UL.txt" << std::endl;

    fLeadingMuPt = 29.;
    fSubLeadingMuPt = 10.;
    fID_SF        = EffTable("./lepEff/eff_240728/Run2017_UL_ID.txt");
    fISO_SF       = EffTable("./lepEff/eff_240728/Run2017_UL_ISO.txt");
    fTRIG_DataEff = EffTable("./lepEff/eff_240728/Run2017_UL_TRIG_DATAeff.txt");
    fTRIG_MCEff   = EffTable("./lepEff/eff_240728/Run2017_UL_TRIG_MCeff.txt");

  } else {
    fRoccoR = new RoccoR("./RoccoR/RoccoR2018UL.txt");
    fPuReweighting = new edm::LumiReWeighting("./pileuInfo/MC_2018.root", "./pileuInfo/PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root", "pileup", "pileup");
    std::cout << "MC PU : ./pileuInfo/MC_2018.root" << std::endl;
    std::cout << "DATA PU : ./pileuInfo/PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root" << std::endl;  
    std::cout << "Roccor : ./RoccoR/RoccoR2018UL.txt" << std::endl;
  
    fLeadingMuPt = 26.;
    fSubLeadingMuPt = 10.;
    fID_SF        = EffTable("./lepEff/eff_240728/Run2018_UL_ID.txt");
    fISO_SF       = EffTable("./lepEff/eff_240728/Run2018_UL_ISO.txt");
    fTRIG_DataEff = EffTable("./lepEff/eff_240728/Run2018_UL_TRIG_DATAeff.txt");
    fTRIG_MCEff   = EffTable("./lepEff/eff_240728/Run2018_UL_TRIG_MCeff.txt");
  }  
}

// double ssb_analysis::GetEventTriggerEfficiencyScaleFactor(double pt_1, double eta_1, double pt_2, double eta_2) const {

//   double mu_1_data = fTRIG_DataEff.getEfficiency(pt_1, eta_1);
//   double mu_2_data = fTRIG_DataEff.getEfficiency(pt_2, eta_2);

//   double mu_1_mc = fTRIG_MCEff.getEfficiency(pt_1, eta_1);
//   double mu_2_mc = fTRIG_MCEff.getEfficiency(pt_2, eta_2);

//   return ( 1 - (1 - mu_1_data) * (1 - mu_2_data) ) / ( 1 - (1 - mu_1_mc) * (1 - mu_2_mc) );

// }

double ssb_analysis::GetNormalization() {

  std::map<std::string, double> normSFmap;

  normSFmap.insert(std::make_pair("DY_10To50", 19.1343));
  normSFmap.insert(std::make_pair("DY_50", 1.46281));
  normSFmap.insert(std::make_pair("TTJets", 0.19271));
  normSFmap.insert(std::make_pair("WJets", 74.3554));
  normSFmap.insert(std::make_pair("WW", 0.533925));
  normSFmap.insert(std::make_pair("WZ", 0.59207));
  normSFmap.insert(std::make_pair("ZZ", 0.574227));
  normSFmap.insert(std::make_pair("tW_antitop", 0.184339));
  normSFmap.insert(std::make_pair("tW_top", 0.183816));
  normSFmap.insert(std::make_pair("Run2016B", 1.));
  normSFmap.insert(std::make_pair("Run2016C", 1.));
  normSFmap.insert(std::make_pair("Run2016D", 1.));
  normSFmap.insert(std::make_pair("Run2016E", 1.));
  normSFmap.insert(std::make_pair("Run2016F", 1.));
  normSFmap.insert(std::make_pair("Run2016G", 1.));
  normSFmap.insert(std::make_pair("Run2016HV2", 1.));
  normSFmap.insert(std::make_pair("Run2016HV3", 1.));

  return normSFmap.at(fProcessName);
}

ssb_analysis::~ssb_analysis() {
  if (!fChain)
    return;
  delete fChain->GetCurrentFile();
  delete fout;
  delete fPuReweighting;
  delete fRoccoR;
}

#endif
