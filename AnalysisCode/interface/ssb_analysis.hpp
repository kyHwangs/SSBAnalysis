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

#include "./../roccor.2016.v3/RoccoR.h"

#include "./../KinSolv/KinematicReconstruction.h"
#include "./../KinSolv/KinematicReconstructionSolution.h"
#include "./../KinSolv/analysisUtils.h"

using namespace std;

class ssb_analysis : public SSBTree {
public:
  // declare functions
  ssb_analysis(TTree *tree = 0, bool fIsData_ = 0,
               std::string fProcessName_ = 0, bool fRoccor_enalbed_ = 0,
               bool fIDISO_enalbed_ = 0, bool fTRIGG_enalbed_ = 0,
               float fMassCut_ = 0);
  virtual ~ssb_analysis();

  // basic frame
  virtual void Loop(char *logfile);
  void Start(int genLoopon);
  void End();

  // user define functions
  void SetOutputFileName(char *outname);
  void DeclareHistos();

  bool GetIsData() { return fIsData; }
  std::string GetProcessName() { return fProcessName; }

  double GetNormalization();

  double GetPUweight(int npv) { return fPuReweighting->weight(npv); }
  double GetPUweight(float npv) { return fPuReweighting->weight(npv); }

  // double GetIdisoSF(double pT, double eta);
  // double GetTriggerSF(double pT1, double eta1, double pt2, double eta2);

  struct stdMUON {
    TLorentzVector v;
    TLorentzVector _v;
    int charge;
    bool id;
    float iso;
    int nTrackLayer;

    stdMUON(TLorentzVector v_, int charge_, bool id_, float iso_,
            int nTrackLayer_ = 0)
        : v(v_), _v(v_), charge(charge_), id(id_), iso(iso_),
          nTrackLayer(nTrackLayer_) {}

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
  std::string fProcessName;
  bool fRoccor_enalbed;
  bool fIDISO_enalbed;
  bool fTRIGG_enalbed;

  float fMassCut;

  // vector for ChargeMisId
  edm::LumiReWeighting *fPuReweighting;
  RoccoR *fRoccoR;
  EffTable fIDISO_SF;
  EffTable fTRIG_SF;
  EffTable fTracking_SF;

public:
  // declare histograms
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

  TH1D *h_ValidInnerTrackHit_BeforeSel;
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

ssb_analysis::ssb_analysis(TTree *tree, bool fIsData_,
                           std::string fProcessName_, bool fRoccor_enalbed_,
                           bool fIDISO_enalbed_, bool fTRIGG_enalbed_,
                           float fMassCut_)
    : fIsData(fIsData_), fProcessName(fProcessName_),
      fRoccor_enalbed(fRoccor_enalbed_), fIDISO_enalbed(fIDISO_enalbed_),
      fTRIGG_enalbed(fTRIGG_enalbed_), fMassCut(fMassCut_) {
  if (tree == 0) {
    printf("ERROR: Can't find any input tree.\n");
  }
  Init(tree);

  // initializing HERE !
  fPuReweighting = new edm::LumiReWeighting(
      "./pileuInfo/MC_Moriond.root",
      "./pileuInfo/PU_2016_69p2_36000_XSecCentral.root", "pileup", "pileup");
  fRoccoR = new RoccoR("./roccor.2016.v3/rcdata.2016.v3");
  fIDISO_SF = EffTable("./lepEff/EMPP/Run2016UL_IDISO_eta_pt.txt");
  fTRIG_SF = EffTable("./lepEff/EMPP/Run2016UL_TRIG_eta_pt.txt");
  fTracking_SF = EffTable("./lepEff/EMPP/Run2016_Tracking_eta.txt");
}

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
