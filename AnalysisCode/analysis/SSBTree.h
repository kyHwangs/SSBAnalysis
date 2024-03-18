//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 18 22:12:59 2020 by ROOT version 6.06/01
// from TTree SSBTree/Tree for Physics Analyses at CMS
// found on file: dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/Run2016/MC/TT_TuneCUETP8M1_13TeV_powheg_80X_V2_forMoriond/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_20200517_160544/200517_071434/0000/SSBTree_100.root
//////////////////////////////////////////////////////////

#ifndef SSBTree_h
#define SSBTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "TClonesArray.h"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

using namespace std;

class SSBTree {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
  Int_t           Info_EventNumber;
  Int_t           Info_Luminosity;
  Int_t           Info_RunNumber;
  Char_t          Info_isData;
  Int_t           Channel_Idx;
  Int_t           Channel_Idx_Final;
  Int_t           Channel_Lepton_Count;
  Int_t           Channel_Lepton_Count_Final;
  Int_t           Channel_Jets;
  Int_t           Channel_Jets_Abs;
  Double_t        L1_PreFire_Central;
  Double_t        L1_PreFire_Up;
  Double_t        L1_PreFire_Down;
  vector<string>  *METFilter_Name;
  vector<bool>    *METFilter_isError;
  vector<bool>    *METFilter_isPass;
  vector<bool>    *METFilter_isRun;
  vector<string>  *METFilterAdd_Name;
  vector<bool>    *METFilterAdd_isPass;
  TClonesArray    *Elec;
  vector<int>     *Elec_Charge;
  vector<bool>    *Elec_ChargeId_GsfCtf;
  vector<bool>    *Elec_ChargeId_GsfCtfPx;
  vector<bool>    *Elec_ChargeId_GsfPx;
  vector<double>  *Elec_Charge_CtfTr;
  vector<double>  *Elec_Charge_GsfTr;
  vector<double>  *Elec_Charge_Px;
  vector<bool>    *Elec_Conversion;
  Int_t           Elec_Count;
  vector<int>     *Elec_Id_Loose;
  vector<int>     *Elec_Id_RobustHighEnergy;
  vector<int>     *Elec_Id_RobustLoose;
  vector<int>     *Elec_Id_RobustTight;
  vector<int>     *Elec_Id_Tight;
  vector<int>     *Elec_Inner_Hit;
  vector<bool>    *Elec_MVA_Medium;
  vector<bool>    *Elec_MVA_Tight;
  vector<bool>    *Elec_MVA_HZZ;
  vector<float>   *Elec_MVA_Values;
  vector<int>     *Elec_MVA_Categories;
  vector<float>   *Elec_MVA_HZZ_Values;
  vector<int>     *Elec_MVA_HZZ_Categories;
  vector<double>  *Elec_PFIsoRho03;
  vector<double>  *Elec_PFIsoRho04;
  vector<bool>    *Elec_PFIsoValid;
  vector<double>  *Elec_PFIsodBeta03;
  vector<double>  *Elec_PFIsodBeta04;
  vector<bool>    *Elec_SCB_Loose;
  vector<bool>    *Elec_SCB_Medium;
  vector<bool>    *Elec_SCB_Tight;
  vector<bool>    *Elec_SCB_Veto;
  vector<float>   *Elec_SCB_dEtaIn;
  vector<float>   *Elec_SCB_dPhiIn;
  vector<float>   *Elec_SCB_hOverE;
  vector<bool>    *Elec_SCB_HEEP;
  vector<float>   *Elec_SCB_ooEmooP;
  vector<float>   *Elec_SCB_sigmaIetaIeta;
  vector<double>  *Elec_ScaleUp;
  vector<double>  *Elec_ScaleDown;
  vector<double>  *Elec_SigmaUp;
  vector<double>  *Elec_SigmaDown;
  vector<double>  *Elec_ScSmUpUp;
  vector<double>  *Elec_ScSmUpDown;
  vector<double>  *Elec_ScSmDownUp;
  vector<double>  *Elec_ScSmDownDown;
  vector<double>  *Elec_Supercluster_Eta;
  vector<double>  *Elec_Track_CtfdXY;
  vector<double>  *Elec_Track_CtfdZ;
  vector<double>  *Elec_Track_GsfdXY;
  vector<double>  *Elec_Track_GsfdZ;
  vector<int>     *Elec_pdgId;
  vector<double>  *Elec_relIso03;
  vector<double>  *Elec_relIso04;
  vector<double>  *Elec_Rnd;
  Char_t          Filter_Greedy_Muon;
  Char_t          Filter_Inconsistent_MuonPt;
  Char_t          Filter_PFReco_Muon;
  vector<bool>    *Filter_PV;
  TClonesArray    *GenJet;
  Int_t           GenJet_Count;
  vector<float>   *GenJet_ECalEnergy;
  vector<float>   *GenJet_HCalEnergy;
  TClonesArray    *GenMET;
  Int_t           GenMET_Count;
  TClonesArray    *GenPar;
  Int_t           GenPar_Count;
  vector<int>     *GenPar_Dau1_Idx;
  vector<int>     *GenPar_Dau2_Idx;
  vector<int>     *GenPar_Dau_Counter;
  vector<int>     *GenPar_Idx;
  vector<int>     *GenPar_Mom1_Idx;
  vector<int>     *GenPar_Mom2_Idx;
  vector<int>     *GenPar_Mom_Counter;
  vector<int>     *GenPar_Status;
  vector<int>     *GenPar_pdgId;
  Double_t        Gen_EventWeight;
  TClonesArray    *GenTop;
  TClonesArray    *GenAnTop;
  TClonesArray    *GenBJet;
  Int_t           GenBJet_Count;
  TClonesArray    *GenBHad;
  Int_t           GenBHad_Count;
  vector<int>     *GenBHad_pdgId;
  vector<int>     *GenBHad_Flavour;
  vector<int>     *GenBHad_FromTopWeakDecay;
  TClonesArray    *Jet;
  vector<int>     *Jet_Charge;
  Int_t           Jet_Count;
  vector<double>  *Jet_EnShiftedDown;
  vector<double>  *Jet_EnShiftedUp;
  vector<int>     *Jet_PartonFlavour;
  vector<int>     *Jet_HadronFlavour;
  vector<int>     *Jet_PFId;
  vector<int>     *Jet_PFIdVeto;
  vector<float>   *Jet_PhiResolution_MC;
  vector<float>   *Jet_PhiResolution_DATA;
  vector<double>  *Jet_EnergyResolution_MC;
  vector<double>  *Jet_EnergyResolution_DATA;
  vector<double>  *Jet_EnergyResolution_SF;
  vector<double>  *Jet_EnergyResolution_SFDown;
  vector<double>  *Jet_EnergyResolution_SFUp;
  vector<int>     *Jet_PileUpId;
  vector<float>   *Jet_PileUpMVA;
  vector<float>   *Jet_bDisc;
  vector<float>   *Jet_bDisc_pfCombinedInclusiveSecondaryVertexV2BJetTags;
  vector<float>   *Jet_bDisc_softPFMuonBJetTags;
  vector<float>   *Jet_bDisc_softPFMuonByIP3dBJetTags;
  vector<float>   *Jet_bDisc_softPFElectronByPtBJetTags;
  vector<float>   *Jet_bDisc_softPFElectronBJetTags;
  vector<float>   *Jet_bDisc_softPFMuonByPtBJetTags;
  vector<float>   *Jet_bDisc_softPFElectronByIP3dBJetTags;
  vector<float>   *Jet_bDisc_softPFMuonByIP2dBJetTags;
  vector<float>   *Jet_bDisc_softPFElectronByIP2dBJetTags;
  vector<bool>    *Jet_isJet;
  Double_t        LHE_Central;
  vector<int>     *LHE_Id;
  vector<double>  *LHE_Weight;
  Double_t        Frag_Cen_Weight;
  Double_t        Frag_Up_Weight;
  Double_t        Frag_Down_Weight;
  Double_t        Frag_Peterson_Weight;
  Double_t        Semilep_BrUp_Weight;
  Double_t        Semilep_BrDown_Weight;
  TClonesArray    *MET;
  Double_t        MET_Significance;
  vector<double>  *MET_JetEnShiftedUp_PT;
  vector<double>  *MET_JetEnShiftedUp_Phi;
  vector<double>  *MET_JetEnShiftedDown_PT;
  vector<double>  *MET_JetEnShiftedDown_Phi;
  vector<double>  *MET_MuonEnShiftedUp_PT;
  vector<double>  *MET_MuonEnShiftedUp_Phi;
  vector<double>  *MET_MuonEnShiftedDown_PT;
  vector<double>  *MET_MuonEnShiftedDown_Phi;
  vector<double>  *MET_ElectronEnShiftedUp_PT;
  vector<double>  *MET_ElectronEnShiftedUp_Phi;
  vector<double>  *MET_ElectronEnShiftedDown_PT;
  vector<double>  *MET_ElectronEnShiftedDown_Phi;
  vector<double>  *MET_UnclusteredEnShiftedUp_PT;
  vector<double>  *MET_UnclusteredEnShiftedUp_Phi;
  vector<double>  *MET_UnclusteredEnShiftedDown_PT;
  vector<double>  *MET_UnclusteredEnShiftedDown_Phi;
  vector<double>  *MET_JetResShiftedUp_PT;
  vector<double>  *MET_JetResShiftedUp_Phi;
  vector<double>  *MET_JetResShiftedDown_PT;
  vector<double>  *MET_JetResShiftedDown_Phi;
  TClonesArray    *METEGClean;
  Double_t        METEGClean_Significance;
  vector<double>  *METEGClean_JetEnShiftedUp_PT;
  vector<double>  *METEGClean_JetEnShiftedUp_Phi;
  vector<double>  *METEGClean_JetEnShiftedDown_PT;
  vector<double>  *METEGClean_JetEnShiftedDown_Phi;
  vector<double>  *METEGClean_MuonEnShiftedUp_PT;
  vector<double>  *METEGClean_MuonEnShiftedUp_Phi;
  vector<double>  *METEGClean_MuonEnShiftedDown_PT;
  vector<double>  *METEGClean_MuonEnShiftedDown_Phi;
  vector<double>  *METEGClean_ElectronEnShiftedUp_PT;
  vector<double>  *METEGClean_ElectronEnShiftedUp_Phi;
  vector<double>  *METEGClean_ElectronEnShiftedDown_PT;
  vector<double>  *METEGClean_ElectronEnShiftedDown_Phi;
  vector<double>  *METEGClean_UnclusteredEnShiftedUp_PT;
  vector<double>  *METEGClean_UnclusteredEnShiftedUp_Phi;
  vector<double>  *METEGClean_UnclusteredEnShiftedDown_PT;
  vector<double>  *METEGClean_UnclusteredEnShiftedDown_Phi;
  vector<double>  *METEGClean_JetResShiftedUp_PT;
  vector<double>  *METEGClean_JetResShiftedUp_Phi;
  vector<double>  *METEGClean_JetResShiftedDown_PT;
  vector<double>  *METEGClean_JetResShiftedDown_Phi;
  TClonesArray    *METMUEGClean;
  Double_t        METMUEGClean_Significance;
  vector<double>  *METMUEGClean_JetEnShiftedUp_PT;
  vector<double>  *METMUEGClean_JetEnShiftedUp_Phi;
  vector<double>  *METMUEGClean_JetEnShiftedDown_PT;
  vector<double>  *METMUEGClean_JetEnShiftedDown_Phi;
  vector<double>  *METMUEGClean_MuonEnShiftedUp_PT;
  vector<double>  *METMUEGClean_MuonEnShiftedUp_Phi;
  vector<double>  *METMUEGClean_MuonEnShiftedDown_PT;
  vector<double>  *METMUEGClean_MuonEnShiftedDown_Phi;
  vector<double>  *METMUEGClean_ElectronEnShiftedUp_PT;
  vector<double>  *METMUEGClean_ElectronEnShiftedUp_Phi;
  vector<double>  *METMUEGClean_ElectronEnShiftedDown_PT;
  vector<double>  *METMUEGClean_ElectronEnShiftedDown_Phi;
  vector<double>  *METMUEGClean_UnclusteredEnShiftedUp_PT;
  vector<double>  *METMUEGClean_UnclusteredEnShiftedUp_Phi;
  vector<double>  *METMUEGClean_UnclusteredEnShiftedDown_PT;
  vector<double>  *METMUEGClean_UnclusteredEnShiftedDown_Phi;
  vector<double>  *METMUEGClean_JetResShiftedUp_PT;
  vector<double>  *METMUEGClean_JetResShiftedUp_Phi;
  vector<double>  *METMUEGClean_JetResShiftedDown_PT;
  vector<double>  *METMUEGClean_JetResShiftedDown_Phi;
  TClonesArray    *METMUEGCleanCor;
  Double_t        METMUEGCleanCor_Significance;
  vector<double>  *METMUEGCleanCor_JetEnShiftedUp_PT;
  vector<double>  *METMUEGCleanCor_JetEnShiftedUp_Phi;
  vector<double>  *METMUEGCleanCor_JetEnShiftedDown_PT;
  vector<double>  *METMUEGCleanCor_JetEnShiftedDown_Phi;
  vector<double>  *METMUEGCleanCor_MuonEnShiftedUp_PT;
  vector<double>  *METMUEGCleanCor_MuonEnShiftedUp_Phi;
  vector<double>  *METMUEGCleanCor_MuonEnShiftedDown_PT;
  vector<double>  *METMUEGCleanCor_MuonEnShiftedDown_Phi;
  vector<double>  *METMUEGCleanCor_ElectronEnShiftedUp_PT;
  vector<double>  *METMUEGCleanCor_ElectronEnShiftedUp_Phi;
  vector<double>  *METMUEGCleanCor_ElectronEnShiftedDown_PT;
  vector<double>  *METMUEGCleanCor_ElectronEnShiftedDown_Phi;
  vector<double>  *METMUEGCleanCor_UnclusteredEnShiftedUp_PT;
  vector<double>  *METMUEGCleanCor_UnclusteredEnShiftedUp_Phi;
  vector<double>  *METMUEGCleanCor_UnclusteredEnShiftedDown_PT;
  vector<double>  *METMUEGCleanCor_UnclusteredEnShiftedDown_Phi;
  vector<double>  *METMUEGCleanCor_JetResShiftedUp_PT;
  vector<double>  *METMUEGCleanCor_JetResShiftedUp_Phi;
  vector<double>  *METMUEGCleanCor_JetResShiftedDown_PT;
  vector<double>  *METMUEGCleanCor_JetResShiftedDown_Phi;
  TClonesArray    *METMUCleanCor;
  Double_t        METMUCleanCor_Significance;
  vector<double>  *METMUCleanCor_JetEnShiftedUp_PT;
  vector<double>  *METMUCleanCor_JetEnShiftedUp_Phi;
  vector<double>  *METMUCleanCor_JetEnShiftedDown_PT;
  vector<double>  *METMUCleanCor_JetEnShiftedDown_Phi;
  vector<double>  *METMUCleanCor_MuonEnShiftedUp_PT;
  vector<double>  *METMUCleanCor_MuonEnShiftedUp_Phi;
  vector<double>  *METMUCleanCor_MuonEnShiftedDown_PT;
  vector<double>  *METMUCleanCor_MuonEnShiftedDown_Phi;
  vector<double>  *METMUCleanCor_ElectronEnShiftedUp_PT;
  vector<double>  *METMUCleanCor_ElectronEnShiftedUp_Phi;
  vector<double>  *METMUCleanCor_ElectronEnShiftedDown_PT;
  vector<double>  *METMUCleanCor_ElectronEnShiftedDown_Phi;
  vector<double>  *METMUCleanCor_UnclusteredEnShiftedUp_PT;
  vector<double>  *METMUCleanCor_UnclusteredEnShiftedUp_Phi;
  vector<double>  *METMUCleanCor_UnclusteredEnShiftedDown_PT;
  vector<double>  *METMUCleanCor_UnclusteredEnShiftedDown_Phi;
  vector<double>  *METMUCleanCor_JetResShiftedUp_PT;
  vector<double>  *METMUCleanCor_JetResShiftedUp_Phi;
  vector<double>  *METMUCleanCor_JetResShiftedDown_PT;
  vector<double>  *METMUCleanCor_JetResShiftedDown_Phi;
  TClonesArray    *METUnCor;
  Double_t        METUnCor_Significance;
  vector<double>  *METUnCor_JetEnShiftedUp_PT;
  vector<double>  *METUnCor_JetEnShiftedUp_Phi;
  vector<double>  *METUnCor_JetEnShiftedDown_PT;
  vector<double>  *METUnCor_JetEnShiftedDown_Phi;
  vector<double>  *METUnCor_MuonEnShiftedUp_PT;
  vector<double>  *METUnCor_MuonEnShiftedUp_Phi;
  vector<double>  *METUnCor_MuonEnShiftedDown_PT;
  vector<double>  *METUnCor_MuonEnShiftedDown_Phi;
  vector<double>  *METUnCor_ElectronEnShiftedUp_PT;
  vector<double>  *METUnCor_ElectronEnShiftedUp_Phi;
  vector<double>  *METUnCor_ElectronEnShiftedDown_PT;
  vector<double>  *METUnCor_ElectronEnShiftedDown_Phi;
  vector<double>  *METUnCor_UnclusteredEnShiftedUp_PT;
  vector<double>  *METUnCor_UnclusteredEnShiftedUp_Phi;
  vector<double>  *METUnCor_UnclusteredEnShiftedDown_PT;
  vector<double>  *METUnCor_UnclusteredEnShiftedDown_Phi;
  vector<double>  *METUnCor_JetResShiftedUp_PT;
  vector<double>  *METUnCor_JetResShiftedUp_Phi;
  vector<double>  *METUnCor_JetResShiftedDown_PT;
  vector<double>  *METUnCor_JetResShiftedDown_Phi;

  TClonesArray    *Muon;
  TClonesArray    *GenMuon;
  vector<int>     *Muon_Charge;
  Int_t           Muon_Count;
  vector<double>  *Muon_PFIsodBeta03;
  vector<double>  *Muon_PFIsodBeta04;
  vector<bool>    *Muon_isHighPt;
  vector<bool>    *Muon_isLoose;
  vector<bool>    *Muon_isMedium;
  vector<bool>    *Muon_isMedium2016;
  vector<bool>    *Muon_isSoft;
  vector<bool>    *Muon_isTight;
  vector<int>     *Muon_pdgId;
  vector<double>  *Muon_rand1;
  vector<double>  *Muon_rand2;
  vector<double>  *Muon_relIso03;
  vector<double>  *Muon_relIso04;
  vector<int>     *Muon_trackerLayers;


  // ssbtree->Branch("Muon", "TClonesArray", &VariableBox_LorentzVector["Muon"], 32000, 0);
  // VariableBox_LorentzVector["Muon"]->BypassStreamer();
  // ssbtree->Branch("GenMuon", "TClonesArray", &VariableBox_LorentzVector["GenMuon"], 32000, 0);
  // VariableBox_LorentzVector["GenMuon"]->BypassStreamer();
  // ssbtree->Branch("Muon_Charge", &VectorBox_Int["Muon_Charge"]);
  // ssbtree->Branch("Muon_Count", &VariableBox_Int["Muon_Count"], "Muon_Count/I");
  // ssbtree->Branch("Muon_PFIsodBeta03", &VectorBox_Double["Muon_PFIsodBeta03"]);
  // ssbtree->Branch("Muon_PFIsodBeta04", &VectorBox_Double["Muon_PFIsodBeta04"]);
  // ssbtree->Branch("Muon_isHighPt", &VectorBox_Bool["Muon_isHighPt"]);
  // ssbtree->Branch("Muon_isHighTrkPt", &VectorBox_Bool["Muon_isHighTrkPt"]);
  // ssbtree->Branch("Muon_isLoose", &VectorBox_Bool["Muon_isLoose"]);
  // ssbtree->Branch("Muon_isMedium", &VectorBox_Bool["Muon_isMedium"]);
  // ssbtree->Branch("Muon_isMediumPrompt", &VectorBox_Bool["Muon_isMediumPrompt"]);
  // ssbtree->Branch("Muon_isSoft", &VectorBox_Bool["Muon_isSoft"]);
  // ssbtree->Branch("Muon_isTight", &VectorBox_Bool["Muon_isTight"]);

  // ssbtree->Branch("Muon_isHighPt_calc", &VectorBox_Bool["Muon_isHighPt_calc"]);
  // ssbtree->Branch("Muon_isHighTrkPt_calc", &VectorBox_Bool["Muon_isHighTrkPt_calc"]);
  // ssbtree->Branch("Muon_isLoose_calc", &VectorBox_Bool["Muon_isLoose_calc"]);
  // ssbtree->Branch("Muon_isMedium_calc", &VectorBox_Bool["Muon_isMedium_calc"]);
  // ssbtree->Branch("Muon_isTight_calc", &VectorBox_Bool["Muon_isTight_calc"]);

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

  // ssbtree->Branch("Muon_isPFIsoVeryLoose", &VectorBox_Bool["Muon_isPFIsoVeryLoose"]);
  // ssbtree->Branch("Muon_isPFIsoLoose", &VectorBox_Bool["Muon_isPFIsoLoose"]);
  // ssbtree->Branch("Muon_isPFIsoMedium", &VectorBox_Bool["Muon_isPFIsoMedium"]);
  // ssbtree->Branch("Muon_isPFIsoTight", &VectorBox_Bool["Muon_isPFIsoTight"]);
  // ssbtree->Branch("Muon_isPFIsoVeryTight", &VectorBox_Bool["Muon_isPFIsoVeryTight"]);
  // ssbtree->Branch("Muon_isPFIsoVeryVeryTight", &VectorBox_Bool["Muon_isPFIsoVeryVeryTight"]);

  // ssbtree->Branch("Muon_pdgId", &VectorBox_Int["Muon_pdgId"]);
  // ssbtree->Branch("Muon_rand1", &VectorBox_Double["Muon_rand1"]);
  // ssbtree->Branch("Muon_rand2", &VectorBox_Double["Muon_rand2"]);
  // ssbtree->Branch("Muon_relIso03", &VectorBox_Double["Muon_relIso03"]);
  // ssbtree->Branch("Muon_relIso04", &VectorBox_Double["Muon_relIso04"]);
  // ssbtree->Branch("Muon_trackerLayers", &VectorBox_Int["Muon_trackerLayers"]);

  // ssbtree->Branch("Muon_tuneP_Pt", &VectorBox_Double["Muon_tuneP_Pt"]);
  // ssbtree->Branch("Muon_tuneP_Eta", &VectorBox_Double["Muon_tuneP_Eta"]);
  // ssbtree->Branch("Muon_tuneP_Phi", &VectorBox_Double["Muon_tuneP_Phi"]);
  // ssbtree->Branch("Muon_tuneP_Charge", &VectorBox_Int["Muon_tuneP_Charge"]);


  vector<double>  *PDFWeight_BjorkenX1;
  vector<double>  *PDFWeight_BjorkenX2;
  vector<double>  *PDFWeight_Cent;
  vector<double>  *PDFWeight_Cent_Down;
  vector<double>  *PDFWeight_Cent_Up;
  vector<int>     *PDFWeight_Id1;
  vector<int>     *PDFWeight_Id2;
  vector<double>  *PDFWeight_PDF1;
  vector<double>  *PDFWeight_PDF2;
  vector<double>  *PDFWeight_Q;
  vector<double>  *PDFWeight_Var1_Down;
  vector<double>  *PDFWeight_Var1_Up;
  vector<double>  *PDFWeight_Var2_Down;
  vector<double>  *PDFWeight_Var2_Up;
  TClonesArray    *Photon;
  Int_t           Photon_Count;
  vector<bool>    *Photon_SCB_Loose;
  vector<bool>    *Photon_SCB_Medium;
  vector<bool>    *Photon_SCB_Tight;
  vector<float>   *Photon_ChaHadIso;
  vector<float>   *Photon_NeuHadIso;
  vector<float>   *Photon_PhoIso;
  vector<float>   *Photon_WorstChargedIso;
  vector<float>   *Photon_ChaHadEffArea;
  vector<float>   *Photon_NeuHadEffArea;
  vector<float>   *Photon_PhoHadEffArea;
  vector<double>  *Photon_R9;
  vector<double>  *Photon_HoverE;
  vector<double>  *Photon_SuperCluster_Eta;
  vector<double>  *Photon_SuperCluster_Phi;
  vector<double>  *Photon_SuperCluster_EtaWidth;
  vector<double>  *Photon_SuperCluster_PhiWidth;
  vector<double>  *Photon_SuperCluster_Energy;
  vector<double>  *Photon_SuperCluster_RawEnergy;
  vector<bool>    *Photon_ElectronVeto;
  vector<double>  *Photon_Full5x5_SigmaIetaIeta;
  vector<bool>    *Photon_MVANonTrig_Tight;
  Int_t           PV_Count;
  Int_t           PileUp_Count_Interaction;
  Float_t         PileUp_Count_Intime;
  vector<double>  *Rho;
  vector<string>  *Trigger_Name;
  vector<unsigned int> *Trigger_PreScale;
  vector<bool>    *Trigger_isError;
  vector<bool>    *Trigger_isPass;
  vector<bool>    *Trigger_isRun;
  vector<double>  *Vertex_SumPtSquare;
  vector<double>  *Vertex_X;
  vector<double>  *Vertex_X_Error;
  vector<double>  *Vertex_Y;
  vector<double>  *Vertex_Y_Error;
  vector<double>  *Vertex_Z;
  vector<double>  *Vertex_Z_Error;

  // List of branches
  TBranch        *b_Info_EventNumber;   //!
  TBranch        *b_Info_Luminosity;   //!
  TBranch        *b_Info_RunNumber;   //!
  TBranch        *b_Info_isData;   //!
  TBranch        *b_Channel_Idx;   //!
  TBranch        *b_Channel_Idx_Final;   //!
  TBranch        *b_Channel_Lepton_Count;   //!
  TBranch        *b_Channel_Lepton_Count_Final;   //!
  TBranch        *b_Channel_Jets;   //!
  TBranch        *b_Channel_Jets_Abs;   //!
  TBranch        *b_L1_PreFire_Central;   //!
  TBranch        *b_L1_PreFire_Up;   //!
  TBranch        *b_L1_PreFire_Down;   //!
  TBranch        *b_METFilter_Name;   //!
  TBranch        *b_METFilter_isError;   //!
  TBranch        *b_METFilter_isPass;   //!
  TBranch        *b_METFilter_isRun;   //!
  TBranch        *b_METFilterAdd_Name;   //!
  TBranch        *b_METFilterAdd_isPass;   //!
  TBranch        *b_Elec;   //!
  TBranch        *b_Elec_Charge;   //!
  TBranch        *b_Elec_ChargeId_GsfCtf;   //!
  TBranch        *b_Elec_ChargeId_GsfCtfPx;   //!
  TBranch        *b_Elec_ChargeId_GsfPx;   //!
  TBranch        *b_Elec_Charge_CtfTr;   //!
  TBranch        *b_Elec_Charge_GsfTr;   //!
  TBranch        *b_Elec_Charge_Px;   //!
  TBranch        *b_Elec_Conversion;   //!
  TBranch        *b_Elec_Count;   //!
  TBranch        *b_Elec_Id_Loose;   //!
  TBranch        *b_Elec_Id_RobustHighEnergy;   //!
  TBranch        *b_Elec_Id_RobustLoose;   //!
  TBranch        *b_Elec_Id_RobustTight;   //!
  TBranch        *b_Elec_Id_Tight;   //!
  TBranch        *b_Elec_Inner_Hit;   //!
  TBranch        *b_Elec_MVA_Medium;   //!
  TBranch        *b_Elec_MVA_Tight;   //!
  TBranch        *b_Elec_MVA_HZZ;   //!
  TBranch        *b_Elec_MVA_Values;   //!
  TBranch        *b_Elec_MVA_Categories;   //!
  TBranch        *b_Elec_MVA_HZZ_Values;   //!
  TBranch        *b_Elec_MVA_HZZ_Categories;   //!
  TBranch        *b_Elec_PFIsoRho03;   //!
  TBranch        *b_Elec_PFIsoRho04;   //!
  TBranch        *b_Elec_PFIsoValid;   //!
  TBranch        *b_Elec_PFIsodBeta03;   //!
  TBranch        *b_Elec_PFIsodBeta04;   //!
  TBranch        *b_Elec_SCB_Loose;   //!
  TBranch        *b_Elec_SCB_Medium;   //!
  TBranch        *b_Elec_SCB_Tight;   //!
  TBranch        *b_Elec_SCB_Veto;   //!
  TBranch        *b_Elec_SCB_dEtaIn;   //!
  TBranch        *b_Elec_SCB_dPhiIn;   //!
  TBranch        *b_Elec_SCB_hOverE;   //!
  TBranch        *b_Elec_SCB_HEEP;   //!
  TBranch        *b_Elec_SCB_ooEmooP;   //!
  TBranch        *b_Elec_SCB_sigmaIetaIeta;   //!
  TBranch        *b_Elec_ScaleUp;   //!
  TBranch        *b_Elec_ScaleDown;   //!
  TBranch        *b_Elec_SigmaUp;   //!
  TBranch        *b_Elec_SigmaDown;   //!
  TBranch        *b_Elec_ScSmUpUp;   //!
  TBranch        *b_Elec_ScSmUpDown;   //!
  TBranch        *b_Elec_ScSmDownUp;   //!
  TBranch        *b_Elec_ScSmDownDown;   //!
  TBranch        *b_Elec_Supercluster_Eta;   //!
  TBranch        *b_Elec_Track_CtfdXY;   //!
  TBranch        *b_Elec_Track_CtfdZ;   //!
  TBranch        *b_Elec_Track_GsfdXY;   //!
  TBranch        *b_Elec_Track_GsfdZ;   //!
  TBranch        *b_Elec_pdgId;   //!
  TBranch        *b_Elec_relIso03;   //!
  TBranch        *b_Elec_relIso04;   //!
  TBranch        *b_Elec_Rnd;   //!
  TBranch        *b_Filter_Greedy_Muon;   //!
  TBranch        *b_Filter_Inconsistent_MuonPt;   //!
  TBranch        *b_Filter_PFReco_Muon;   //!
  TBranch        *b_Filter_PV;   //!
  TBranch        *b_GenJet;   //!
  TBranch        *b_GenJet_Count;   //!
  TBranch        *b_GenJet_ECalEnergy;   //!
  TBranch        *b_GenJet_HCalEnergy;   //!
  TBranch        *b_GenMET;   //!
  TBranch        *b_GenMET_Count;   //!
  TBranch        *b_GenPar;   //!
  TBranch        *b_GenPar_Count;   //!
  TBranch        *b_GenPar_Dau1_Idx;   //!
  TBranch        *b_GenPar_Dau2_Idx;   //!
  TBranch        *b_GenPar_Dau_Counter;   //!
  TBranch        *b_GenPar_Idx;   //!
  TBranch        *b_GenPar_Mom1_Idx;   //!
  TBranch        *b_GenPar_Mom2_Idx;   //!
  TBranch        *b_GenPar_Mom_Counter;   //!
  TBranch        *b_GenPar_Status;   //!
  TBranch        *b_GenPar_pdgId;   //!
  TBranch        *b_Gen_EventWeight;   //!
  TBranch        *b_GenTop;   //!
  TBranch        *b_GenAnTop;   //!
  TBranch        *b_GenBJet;   //!
  TBranch        *b_GenBJet_Count;   //!
  TBranch        *b_GenBHad;   //!
  TBranch        *b_GenBHad_Count;   //!
  TBranch        *b_GenBHad_pdgId;   //!
  TBranch        *b_GenBHad_Flavour;   //!
  TBranch        *b_GenBHad_FromTopWeakDecay;   //!
  TBranch        *b_Jet;   //!
  TBranch        *b_Jet_Charge;   //!
  TBranch        *b_Jet_Count;   //!
  TBranch        *b_Jet_EnShiftedDown;   //!
  TBranch        *b_Jet_EnShiftedUp;   //!
  TBranch        *b_Jet_PartonFlavour;   //!
  TBranch        *b_Jet_HadronFlavour;   //!
  TBranch        *b_Jet_PFId;   //!
  TBranch        *b_Jet_PFIdVeto;   //!
  TBranch        *b_Jet_PhiResolution_MC;   //!
  TBranch        *b_Jet_PhiResolution_DATA;   //!
  TBranch        *b_Jet_EnergyResolution_MC;   //!
  TBranch        *b_Jet_EnergyResolution_DATA;   //!
  TBranch        *b_Jet_EnergyResolution_SF;   //!
  TBranch        *b_Jet_EnergyResolution_SFDown;   //!
  TBranch        *b_Jet_EnergyResolution_SFUp;   //!
  TBranch        *b_Jet_PileUpId;   //!
  TBranch        *b_Jet_PileUpMVA;   //!
  TBranch        *b_Jet_bDisc;   //!
  TBranch        *b_Jet_bDisc_pfCombinedInclusiveSecondaryVertexV2BJetTags;   //!
  TBranch        *b_Jet_bDisc_softPFMuonBJetTags;   //!
  TBranch        *b_Jet_bDisc_softPFMuonByIP3dBJetTags;   //!
  TBranch        *b_Jet_bDisc_softPFElectronByPtBJetTags;   //!
  TBranch        *b_Jet_bDisc_softPFElectronBJetTags;   //!
  TBranch        *b_Jet_bDisc_softPFMuonByPtBJetTags;   //!
  TBranch        *b_Jet_bDisc_softPFElectronByIP3dBJetTags;   //!
  TBranch        *b_Jet_bDisc_softPFMuonByIP2dBJetTags;   //!
  TBranch        *b_Jet_bDisc_softPFElectronByIP2dBJetTags;   //!
  TBranch        *b_Jet_isJet;   //!
  TBranch        *b_LHE_Central;   //!
  TBranch        *b_LHE_Id;   //!
  TBranch        *b_LHE_Weight;   //!
  TBranch        *b_Frag_Cen_Weight;   //!
  TBranch        *b_Frag_Up_Weight;   //!
  TBranch        *b_Frag_Down_Weight;   //!
  TBranch        *b_Frag_Peterson_Weight;   //!
  TBranch        *b_Semilep_BrUp_Weight;   //!
  TBranch        *b_Semilep_BrDown_Weight;   //!
  TBranch        *b_MET;   //!
  TBranch        *b_MET_Significance;   //!
  TBranch        *b_MET_JetEnShiftedUp_PT;   //!
  TBranch        *b_MET_JetEnShiftedUp_Phi;   //!
  TBranch        *b_MET_JetEnShiftedDown_PT;   //!
  TBranch        *b_MET_JetEnShiftedDown_Phi;   //!
  TBranch        *b_MET_MuonEnShiftedUp_PT;   //!
  TBranch        *b_MET_MuonEnShiftedUp_Phi;   //!
  TBranch        *b_MET_MuonEnShiftedDown_PT;   //!
  TBranch        *b_MET_MuonEnShiftedDown_Phi;   //!
  TBranch        *b_MET_ElectronEnShiftedUp_PT;   //!
  TBranch        *b_MET_ElectronEnShiftedUp_Phi;   //!
  TBranch        *b_MET_ElectronEnShiftedDown_PT;   //!
  TBranch        *b_MET_ElectronEnShiftedDown_Phi;   //!
  TBranch        *b_MET_UnclusteredEnShiftedUp_PT;   //!
  TBranch        *b_MET_UnclusteredEnShiftedUp_Phi;   //!
  TBranch        *b_MET_UnclusteredEnShiftedDown_PT;   //!
  TBranch        *b_MET_UnclusteredEnShiftedDown_Phi;   //!
  TBranch        *b_MET_JetResShiftedUp_PT;   //!
  TBranch        *b_MET_JetResShiftedUp_Phi;   //!
  TBranch        *b_MET_JetResShiftedDown_PT;   //!
  TBranch        *b_MET_JetResShiftedDown_Phi;   //!
  TBranch        *b_METEGClean;   //!
  TBranch        *b_METEGClean_Significance;   //!
  TBranch        *b_METEGClean_JetEnShiftedUp_PT;   //!
  TBranch        *b_METEGClean_JetEnShiftedUp_Phi;   //!
  TBranch        *b_METEGClean_JetEnShiftedDown_PT;   //!
  TBranch        *b_METEGClean_JetEnShiftedDown_Phi;   //!
  TBranch        *b_METEGClean_MuonEnShiftedUp_PT;   //!
  TBranch        *b_METEGClean_MuonEnShiftedUp_Phi;   //!
  TBranch        *b_METEGClean_MuonEnShiftedDown_PT;   //!
  TBranch        *b_METEGClean_MuonEnShiftedDown_Phi;   //!
  TBranch        *b_METEGClean_ElectronEnShiftedUp_PT;   //!
  TBranch        *b_METEGClean_ElectronEnShiftedUp_Phi;   //!
  TBranch        *b_METEGClean_ElectronEnShiftedDown_PT;   //!
  TBranch        *b_METEGClean_ElectronEnShiftedDown_Phi;   //!
  TBranch        *b_METEGClean_UnclusteredEnShiftedUp_PT;   //!
  TBranch        *b_METEGClean_UnclusteredEnShiftedUp_Phi;   //!
  TBranch        *b_METEGClean_UnclusteredEnShiftedDown_PT;   //!
  TBranch        *b_METEGClean_UnclusteredEnShiftedDown_Phi;   //!
  TBranch        *b_METEGClean_JetResShiftedUp_PT;   //!
  TBranch        *b_METEGClean_JetResShiftedUp_Phi;   //!
  TBranch        *b_METEGClean_JetResShiftedDown_PT;   //!
  TBranch        *b_METEGClean_JetResShiftedDown_Phi;   //!
  TBranch        *b_METMUEGClean;   //!
  TBranch        *b_METMUEGClean_Significance;   //!
  TBranch        *b_METMUEGClean_JetEnShiftedUp_PT;   //!
  TBranch        *b_METMUEGClean_JetEnShiftedUp_Phi;   //!
  TBranch        *b_METMUEGClean_JetEnShiftedDown_PT;   //!
  TBranch        *b_METMUEGClean_JetEnShiftedDown_Phi;   //!
  TBranch        *b_METMUEGClean_MuonEnShiftedUp_PT;   //!
  TBranch        *b_METMUEGClean_MuonEnShiftedUp_Phi;   //!
  TBranch        *b_METMUEGClean_MuonEnShiftedDown_PT;   //!
  TBranch        *b_METMUEGClean_MuonEnShiftedDown_Phi;   //!
  TBranch        *b_METMUEGClean_ElectronEnShiftedUp_PT;   //!
  TBranch        *b_METMUEGClean_ElectronEnShiftedUp_Phi;   //!
  TBranch        *b_METMUEGClean_ElectronEnShiftedDown_PT;   //!
  TBranch        *b_METMUEGClean_ElectronEnShiftedDown_Phi;   //!
  TBranch        *b_METMUEGClean_UnclusteredEnShiftedUp_PT;   //!
  TBranch        *b_METMUEGClean_UnclusteredEnShiftedUp_Phi;   //!
  TBranch        *b_METMUEGClean_UnclusteredEnShiftedDown_PT;   //!
  TBranch        *b_METMUEGClean_UnclusteredEnShiftedDown_Phi;   //!
  TBranch        *b_METMUEGClean_JetResShiftedUp_PT;   //!
  TBranch        *b_METMUEGClean_JetResShiftedUp_Phi;   //!
  TBranch        *b_METMUEGClean_JetResShiftedDown_PT;   //!
  TBranch        *b_METMUEGClean_JetResShiftedDown_Phi;   //!
  TBranch        *b_METMUEGCleanCor;   //!
  TBranch        *b_METMUEGCleanCor_Significance;   //!
  TBranch        *b_METMUEGCleanCor_JetEnShiftedUp_PT;   //!
  TBranch        *b_METMUEGCleanCor_JetEnShiftedUp_Phi;   //!
  TBranch        *b_METMUEGCleanCor_JetEnShiftedDown_PT;   //!
  TBranch        *b_METMUEGCleanCor_JetEnShiftedDown_Phi;   //!
  TBranch        *b_METMUEGCleanCor_MuonEnShiftedUp_PT;   //!
  TBranch        *b_METMUEGCleanCor_MuonEnShiftedUp_Phi;   //!
  TBranch        *b_METMUEGCleanCor_MuonEnShiftedDown_PT;   //!
  TBranch        *b_METMUEGCleanCor_MuonEnShiftedDown_Phi;   //!
  TBranch        *b_METMUEGCleanCor_ElectronEnShiftedUp_PT;   //!
  TBranch        *b_METMUEGCleanCor_ElectronEnShiftedUp_Phi;   //!
  TBranch        *b_METMUEGCleanCor_ElectronEnShiftedDown_PT;   //!
  TBranch        *b_METMUEGCleanCor_ElectronEnShiftedDown_Phi;   //!
  TBranch        *b_METMUEGCleanCor_UnclusteredEnShiftedUp_PT;   //!
  TBranch        *b_METMUEGCleanCor_UnclusteredEnShiftedUp_Phi;   //!
  TBranch        *b_METMUEGCleanCor_UnclusteredEnShiftedDown_PT;   //!
  TBranch        *b_METMUEGCleanCor_UnclusteredEnShiftedDown_Phi;   //!
  TBranch        *b_METMUEGCleanCor_JetResShiftedUp_PT;   //!
  TBranch        *b_METMUEGCleanCor_JetResShiftedUp_Phi;   //!
  TBranch        *b_METMUEGCleanCor_JetResShiftedDown_PT;   //!
  TBranch        *b_METMUEGCleanCor_JetResShiftedDown_Phi;   //!
  TBranch        *b_METMUCleanCor;   //!
  TBranch        *b_METMUCleanCor_Significance;   //!
  TBranch        *b_METMUCleanCor_JetEnShiftedUp_PT;   //!
  TBranch        *b_METMUCleanCor_JetEnShiftedUp_Phi;   //!
  TBranch        *b_METMUCleanCor_JetEnShiftedDown_PT;   //!
  TBranch        *b_METMUCleanCor_JetEnShiftedDown_Phi;   //!
  TBranch        *b_METMUCleanCor_MuonEnShiftedUp_PT;   //!
  TBranch        *b_METMUCleanCor_MuonEnShiftedUp_Phi;   //!
  TBranch        *b_METMUCleanCor_MuonEnShiftedDown_PT;   //!
  TBranch        *b_METMUCleanCor_MuonEnShiftedDown_Phi;   //!
  TBranch        *b_METMUCleanCor_ElectronEnShiftedUp_PT;   //!
  TBranch        *b_METMUCleanCor_ElectronEnShiftedUp_Phi;   //!
  TBranch        *b_METMUCleanCor_ElectronEnShiftedDown_PT;   //!
  TBranch        *b_METMUCleanCor_ElectronEnShiftedDown_Phi;   //!
  TBranch        *b_METMUCleanCor_UnclusteredEnShiftedUp_PT;   //!
  TBranch        *b_METMUCleanCor_UnclusteredEnShiftedUp_Phi;   //!
  TBranch        *b_METMUCleanCor_UnclusteredEnShiftedDown_PT;   //!
  TBranch        *b_METMUCleanCor_UnclusteredEnShiftedDown_Phi;   //!
  TBranch        *b_METMUCleanCor_JetResShiftedUp_PT;   //!
  TBranch        *b_METMUCleanCor_JetResShiftedUp_Phi;   //!
  TBranch        *b_METMUCleanCor_JetResShiftedDown_PT;   //!
  TBranch        *b_METMUCleanCor_JetResShiftedDown_Phi;   //!
  TBranch        *b_METUnCor;   //!
  TBranch        *b_METUnCor_Significance;   //!
  TBranch        *b_METUnCor_JetEnShiftedUp_PT;   //!
  TBranch        *b_METUnCor_JetEnShiftedUp_Phi;   //!
  TBranch        *b_METUnCor_JetEnShiftedDown_PT;   //!
  TBranch        *b_METUnCor_JetEnShiftedDown_Phi;   //!
  TBranch        *b_METUnCor_MuonEnShiftedUp_PT;   //!
  TBranch        *b_METUnCor_MuonEnShiftedUp_Phi;   //!
  TBranch        *b_METUnCor_MuonEnShiftedDown_PT;   //!
  TBranch        *b_METUnCor_MuonEnShiftedDown_Phi;   //!
  TBranch        *b_METUnCor_ElectronEnShiftedUp_PT;   //!
  TBranch        *b_METUnCor_ElectronEnShiftedUp_Phi;   //!
  TBranch        *b_METUnCor_ElectronEnShiftedDown_PT;   //!
  TBranch        *b_METUnCor_ElectronEnShiftedDown_Phi;   //!
  TBranch        *b_METUnCor_UnclusteredEnShiftedUp_PT;   //!
  TBranch        *b_METUnCor_UnclusteredEnShiftedUp_Phi;   //!
  TBranch        *b_METUnCor_UnclusteredEnShiftedDown_PT;   //!
  TBranch        *b_METUnCor_UnclusteredEnShiftedDown_Phi;   //!
  TBranch        *b_METUnCor_JetResShiftedUp_PT;   //!
  TBranch        *b_METUnCor_JetResShiftedUp_Phi;   //!
  TBranch        *b_METUnCor_JetResShiftedDown_PT;   //!
  TBranch        *b_METUnCor_JetResShiftedDown_Phi;   //!
  TBranch        *b_Muon;   //!
  TBranch        *b_GenMuon;   //!
  TBranch        *b_Muon_Charge;   //!
  TBranch        *b_Muon_Count;   //!
  TBranch        *b_Muon_PFIsodBeta03;   //!
  TBranch        *b_Muon_PFIsodBeta04;   //!
  TBranch        *b_Muon_isHighPt;   //!
  TBranch        *b_Muon_isLoose;   //!
  TBranch        *b_Muon_isMedium;   //!
  TBranch        *b_Muon_isMedium2016;   //!
  TBranch        *b_Muon_isSoft;   //!
  TBranch        *b_Muon_isTight;   //!
  TBranch        *b_Muon_pdgId;   //!
  TBranch        *b_Muon_rand1;   //!
  TBranch        *b_Muon_rand2;   //!
  TBranch        *b_Muon_relIso03;   //!
  TBranch        *b_Muon_relIso04;   //!
  TBranch        *b_Muon_trackerLayers;   //!
  TBranch        *b_PDFWeight_BjorkenX1;   //!
  TBranch        *b_PDFWeight_BjorkenX2;   //!
  TBranch        *b_PDFWeight_Cent;   //!
  TBranch        *b_PDFWeight_Cent_Down;   //!
  TBranch        *b_PDFWeight_Cent_Up;   //!
  TBranch        *b_PDFWeight_Id1;   //!
  TBranch        *b_PDFWeight_Id2;   //!
  TBranch        *b_PDFWeight_PDF1;   //!
  TBranch        *b_PDFWeight_PDF2;   //!
  TBranch        *b_PDFWeight_Q;   //!
  TBranch        *b_PDFWeight_Var1_Down;   //!
  TBranch        *b_PDFWeight_Var1_Up;   //!
  TBranch        *b_PDFWeight_Var2_Down;   //!
  TBranch        *b_PDFWeight_Var2_Up;   //!
  TBranch        *b_Photon;   //!
  TBranch        *b_Photon_Count;   //!
  TBranch        *b_Photon_SCB_Loose;   //!
  TBranch        *b_Photon_SCB_Medium;   //!
  TBranch        *b_Photon_SCB_Tight;   //!
  TBranch        *b_Photon_ChaHadIso;   //!
  TBranch        *b_Photon_NeuHadIso;   //!
  TBranch        *b_Photon_PhoIso;   //!
  TBranch        *b_Photon_WorstChargedIso;   //!
  TBranch        *b_Photon_ChaHadEffArea;   //!
  TBranch        *b_Photon_NeuHadEffArea;   //!
  TBranch        *b_Photon_PhoHadEffArea;   //!
  TBranch        *b_Photon_R9;   //!
  TBranch        *b_Photon_HoverE;   //!
  TBranch        *b_Photon_SuperCluster_Eta;   //!
  TBranch        *b_Photon_SuperCluster_Phi;   //!
  TBranch        *b_Photon_SuperCluster_EtaWidth;   //!
  TBranch        *b_Photon_SuperCluster_PhiWidth;   //!
  TBranch        *b_Photon_SuperCluster_Energy;   //!
  TBranch        *b_Photon_SuperCluster_RawEnergy;   //!
  TBranch        *b_Photon_ElectronVeto;   //!
  TBranch        *b_Photon_Full5x5_SigmaIetaIeta;   //!
  TBranch        *b_Photon_MVANonTrig_Tight;   //!
  TBranch        *b_PV_Count;   //!
  TBranch        *b_PileUp_Count_Interaction;   //!
  TBranch        *b_PileUp_Count_Intime;   //!
  TBranch        *b_Rho;   //!
  TBranch        *b_Trigger_Name;   //!
  TBranch        *b_Trigger_PreScale;   //!
  TBranch        *b_Trigger_isError;   //!
  TBranch        *b_Trigger_isPass;   //!
  TBranch        *b_Trigger_isRun;   //!
  TBranch        *b_Vertex_SumPtSquare;   //!
  TBranch        *b_Vertex_X;   //!
  TBranch        *b_Vertex_X_Error;   //!
  TBranch        *b_Vertex_Y;   //!
  TBranch        *b_Vertex_Y_Error;   //!
  TBranch        *b_Vertex_Z;   //!
  TBranch        *b_Vertex_Z_Error;   //!

  SSBTree(TTree *tree=0);
  virtual ~SSBTree();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef SSBTree_cxx
SSBTree::SSBTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {/*
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/Run2016/MC/TT_TuneCUETP8M1_13TeV_powheg_80X_V2_forMoriond/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_20200517_160544/200517_071434/0000/SSBTree_100.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/Run2016/MC/TT_TuneCUETP8M1_13TeV_powheg_80X_V2_forMoriond/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_20200517_160544/200517_071434/0000/SSBTree_100.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/Run2016/MC/TT_TuneCUETP8M1_13TeV_powheg_80X_V2_forMoriond/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/crab_20200517_160544/200517_071434/0000/SSBTree_100.root:/ssbanalyzer");
      dir->GetObject("SSBTree",tree);
*/
   }
   Init(tree);
}

SSBTree::~SSBTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SSBTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SSBTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void SSBTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   METFilter_Name = 0;
   METFilter_isError = 0;
   METFilter_isPass = 0;
   METFilter_isRun = 0;
   METFilterAdd_Name = 0;
   METFilterAdd_isPass = 0;
   Elec = 0;
   Elec_Charge = 0;
   Elec_ChargeId_GsfCtf = 0;
   Elec_ChargeId_GsfCtfPx = 0;
   Elec_ChargeId_GsfPx = 0;
   Elec_Charge_CtfTr = 0;
   Elec_Charge_GsfTr = 0;
   Elec_Charge_Px = 0;
   Elec_Conversion = 0;
   Elec_Id_Loose = 0;
   Elec_Id_RobustHighEnergy = 0;
   Elec_Id_RobustLoose = 0;
   Elec_Id_RobustTight = 0;
   Elec_Id_Tight = 0;
   Elec_Inner_Hit = 0;
   Elec_MVA_Medium = 0;
   Elec_MVA_Tight = 0;
   Elec_MVA_HZZ = 0;
   Elec_MVA_Values = 0;
   Elec_MVA_Categories = 0;
   Elec_MVA_HZZ_Values = 0;
   Elec_MVA_HZZ_Categories = 0;
   Elec_PFIsoRho03 = 0;
   Elec_PFIsoRho04 = 0;
   Elec_PFIsoValid = 0;
   Elec_PFIsodBeta03 = 0;
   Elec_PFIsodBeta04 = 0;
   Elec_SCB_Loose = 0;
   Elec_SCB_Medium = 0;
   Elec_SCB_Tight = 0;
   Elec_SCB_Veto = 0;
   Elec_SCB_dEtaIn = 0;
   Elec_SCB_dPhiIn = 0;
   Elec_SCB_hOverE = 0;
   Elec_SCB_HEEP = 0;
   Elec_SCB_ooEmooP = 0;
   Elec_SCB_sigmaIetaIeta = 0;
   Elec_ScaleUp = 0;
   Elec_ScaleDown = 0;
   Elec_SigmaUp = 0;
   Elec_SigmaDown = 0;
   Elec_ScSmUpUp = 0;
   Elec_ScSmUpDown = 0;
   Elec_ScSmDownUp = 0;
   Elec_ScSmDownDown = 0;
   Elec_Supercluster_Eta = 0;
   Elec_Track_CtfdXY = 0;
   Elec_Track_CtfdZ = 0;
   Elec_Track_GsfdXY = 0;
   Elec_Track_GsfdZ = 0;
   Elec_pdgId = 0;
   Elec_relIso03 = 0;
   Elec_relIso04 = 0;
   Elec_Rnd = 0;
   Filter_PV = 0;
   GenJet = 0;
   GenJet_ECalEnergy = 0;
   GenJet_HCalEnergy = 0;
   GenMET = 0;
   GenPar = 0;
   GenPar_Dau1_Idx = 0;
   GenPar_Dau2_Idx = 0;
   GenPar_Dau_Counter = 0;
   GenPar_Idx = 0;
   GenPar_Mom1_Idx = 0;
   GenPar_Mom2_Idx = 0;
   GenPar_Mom_Counter = 0;
   GenPar_Status = 0;
   GenPar_pdgId = 0;
   GenTop = 0;
   GenAnTop = 0;
   GenBJet = 0;
   GenBHad = 0;
   GenBHad_pdgId = 0;
   GenBHad_Flavour = 0;
   GenBHad_FromTopWeakDecay = 0;
   Jet = 0;
   Jet_Charge = 0;
   Jet_EnShiftedDown = 0;
   Jet_EnShiftedUp = 0;
   Jet_PartonFlavour = 0;
   Jet_HadronFlavour = 0;
   Jet_PFId = 0;
   Jet_PFIdVeto = 0;
   Jet_PhiResolution_MC = 0;
   Jet_PhiResolution_DATA = 0;
   Jet_EnergyResolution_MC = 0;
   Jet_EnergyResolution_DATA = 0;
   Jet_EnergyResolution_SF = 0;
   Jet_EnergyResolution_SFDown = 0;
   Jet_EnergyResolution_SFUp = 0;
   Jet_PileUpId = 0;
   Jet_PileUpMVA = 0;
   Jet_bDisc = 0;
   Jet_bDisc_pfCombinedInclusiveSecondaryVertexV2BJetTags = 0;
   Jet_bDisc_softPFMuonBJetTags = 0;
   Jet_bDisc_softPFMuonByIP3dBJetTags = 0;
   Jet_bDisc_softPFElectronByPtBJetTags = 0;
   Jet_bDisc_softPFElectronBJetTags = 0;
   Jet_bDisc_softPFMuonByPtBJetTags = 0;
   Jet_bDisc_softPFElectronByIP3dBJetTags = 0;
   Jet_bDisc_softPFMuonByIP2dBJetTags = 0;
   Jet_bDisc_softPFElectronByIP2dBJetTags = 0;
   Jet_isJet = 0;
   LHE_Id = 0;
   LHE_Weight = 0;
   MET = 0;
   MET_JetEnShiftedUp_PT = 0;
   MET_JetEnShiftedUp_Phi = 0;
   MET_JetEnShiftedDown_PT = 0;
   MET_JetEnShiftedDown_Phi = 0;
   MET_MuonEnShiftedUp_PT = 0;
   MET_MuonEnShiftedUp_Phi = 0;
   MET_MuonEnShiftedDown_PT = 0;
   MET_MuonEnShiftedDown_Phi = 0;
   MET_ElectronEnShiftedUp_PT = 0;
   MET_ElectronEnShiftedUp_Phi = 0;
   MET_ElectronEnShiftedDown_PT = 0;
   MET_ElectronEnShiftedDown_Phi = 0;
   MET_UnclusteredEnShiftedUp_PT = 0;
   MET_UnclusteredEnShiftedUp_Phi = 0;
   MET_UnclusteredEnShiftedDown_PT = 0;
   MET_UnclusteredEnShiftedDown_Phi = 0;
   MET_JetResShiftedUp_PT = 0;
   MET_JetResShiftedUp_Phi = 0;
   MET_JetResShiftedDown_PT = 0;
   MET_JetResShiftedDown_Phi = 0;
   METEGClean = 0;
   METEGClean_JetEnShiftedUp_PT = 0;
   METEGClean_JetEnShiftedUp_Phi = 0;
   METEGClean_JetEnShiftedDown_PT = 0;
   METEGClean_JetEnShiftedDown_Phi = 0;
   METEGClean_MuonEnShiftedUp_PT = 0;
   METEGClean_MuonEnShiftedUp_Phi = 0;
   METEGClean_MuonEnShiftedDown_PT = 0;
   METEGClean_MuonEnShiftedDown_Phi = 0;
   METEGClean_ElectronEnShiftedUp_PT = 0;
   METEGClean_ElectronEnShiftedUp_Phi = 0;
   METEGClean_ElectronEnShiftedDown_PT = 0;
   METEGClean_ElectronEnShiftedDown_Phi = 0;
   METEGClean_UnclusteredEnShiftedUp_PT = 0;
   METEGClean_UnclusteredEnShiftedUp_Phi = 0;
   METEGClean_UnclusteredEnShiftedDown_PT = 0;
   METEGClean_UnclusteredEnShiftedDown_Phi = 0;
   METEGClean_JetResShiftedUp_PT = 0;
   METEGClean_JetResShiftedUp_Phi = 0;
   METEGClean_JetResShiftedDown_PT = 0;
   METEGClean_JetResShiftedDown_Phi = 0;
   METMUEGClean = 0;
   METMUEGClean_JetEnShiftedUp_PT = 0;
   METMUEGClean_JetEnShiftedUp_Phi = 0;
   METMUEGClean_JetEnShiftedDown_PT = 0;
   METMUEGClean_JetEnShiftedDown_Phi = 0;
   METMUEGClean_MuonEnShiftedUp_PT = 0;
   METMUEGClean_MuonEnShiftedUp_Phi = 0;
   METMUEGClean_MuonEnShiftedDown_PT = 0;
   METMUEGClean_MuonEnShiftedDown_Phi = 0;
   METMUEGClean_ElectronEnShiftedUp_PT = 0;
   METMUEGClean_ElectronEnShiftedUp_Phi = 0;
   METMUEGClean_ElectronEnShiftedDown_PT = 0;
   METMUEGClean_ElectronEnShiftedDown_Phi = 0;
   METMUEGClean_UnclusteredEnShiftedUp_PT = 0;
   METMUEGClean_UnclusteredEnShiftedUp_Phi = 0;
   METMUEGClean_UnclusteredEnShiftedDown_PT = 0;
   METMUEGClean_UnclusteredEnShiftedDown_Phi = 0;
   METMUEGClean_JetResShiftedUp_PT = 0;
   METMUEGClean_JetResShiftedUp_Phi = 0;
   METMUEGClean_JetResShiftedDown_PT = 0;
   METMUEGClean_JetResShiftedDown_Phi = 0;
   METMUEGCleanCor = 0;
   METMUEGCleanCor_JetEnShiftedUp_PT = 0;
   METMUEGCleanCor_JetEnShiftedUp_Phi = 0;
   METMUEGCleanCor_JetEnShiftedDown_PT = 0;
   METMUEGCleanCor_JetEnShiftedDown_Phi = 0;
   METMUEGCleanCor_MuonEnShiftedUp_PT = 0;
   METMUEGCleanCor_MuonEnShiftedUp_Phi = 0;
   METMUEGCleanCor_MuonEnShiftedDown_PT = 0;
   METMUEGCleanCor_MuonEnShiftedDown_Phi = 0;
   METMUEGCleanCor_ElectronEnShiftedUp_PT = 0;
   METMUEGCleanCor_ElectronEnShiftedUp_Phi = 0;
   METMUEGCleanCor_ElectronEnShiftedDown_PT = 0;
   METMUEGCleanCor_ElectronEnShiftedDown_Phi = 0;
   METMUEGCleanCor_UnclusteredEnShiftedUp_PT = 0;
   METMUEGCleanCor_UnclusteredEnShiftedUp_Phi = 0;
   METMUEGCleanCor_UnclusteredEnShiftedDown_PT = 0;
   METMUEGCleanCor_UnclusteredEnShiftedDown_Phi = 0;
   METMUEGCleanCor_JetResShiftedUp_PT = 0;
   METMUEGCleanCor_JetResShiftedUp_Phi = 0;
   METMUEGCleanCor_JetResShiftedDown_PT = 0;
   METMUEGCleanCor_JetResShiftedDown_Phi = 0;
   METMUCleanCor = 0;
   METMUCleanCor_JetEnShiftedUp_PT = 0;
   METMUCleanCor_JetEnShiftedUp_Phi = 0;
   METMUCleanCor_JetEnShiftedDown_PT = 0;
   METMUCleanCor_JetEnShiftedDown_Phi = 0;
   METMUCleanCor_MuonEnShiftedUp_PT = 0;
   METMUCleanCor_MuonEnShiftedUp_Phi = 0;
   METMUCleanCor_MuonEnShiftedDown_PT = 0;
   METMUCleanCor_MuonEnShiftedDown_Phi = 0;
   METMUCleanCor_ElectronEnShiftedUp_PT = 0;
   METMUCleanCor_ElectronEnShiftedUp_Phi = 0;
   METMUCleanCor_ElectronEnShiftedDown_PT = 0;
   METMUCleanCor_ElectronEnShiftedDown_Phi = 0;
   METMUCleanCor_UnclusteredEnShiftedUp_PT = 0;
   METMUCleanCor_UnclusteredEnShiftedUp_Phi = 0;
   METMUCleanCor_UnclusteredEnShiftedDown_PT = 0;
   METMUCleanCor_UnclusteredEnShiftedDown_Phi = 0;
   METMUCleanCor_JetResShiftedUp_PT = 0;
   METMUCleanCor_JetResShiftedUp_Phi = 0;
   METMUCleanCor_JetResShiftedDown_PT = 0;
   METMUCleanCor_JetResShiftedDown_Phi = 0;
   METUnCor = 0;
   METUnCor_JetEnShiftedUp_PT = 0;
   METUnCor_JetEnShiftedUp_Phi = 0;
   METUnCor_JetEnShiftedDown_PT = 0;
   METUnCor_JetEnShiftedDown_Phi = 0;
   METUnCor_MuonEnShiftedUp_PT = 0;
   METUnCor_MuonEnShiftedUp_Phi = 0;
   METUnCor_MuonEnShiftedDown_PT = 0;
   METUnCor_MuonEnShiftedDown_Phi = 0;
   METUnCor_ElectronEnShiftedUp_PT = 0;
   METUnCor_ElectronEnShiftedUp_Phi = 0;
   METUnCor_ElectronEnShiftedDown_PT = 0;
   METUnCor_ElectronEnShiftedDown_Phi = 0;
   METUnCor_UnclusteredEnShiftedUp_PT = 0;
   METUnCor_UnclusteredEnShiftedUp_Phi = 0;
   METUnCor_UnclusteredEnShiftedDown_PT = 0;
   METUnCor_UnclusteredEnShiftedDown_Phi = 0;
   METUnCor_JetResShiftedUp_PT = 0;
   METUnCor_JetResShiftedUp_Phi = 0;
   METUnCor_JetResShiftedDown_PT = 0;
   METUnCor_JetResShiftedDown_Phi = 0;
   Muon = 0;
   GenMuon = 0;
   Muon_Charge = 0;
   Muon_PFIsodBeta03 = 0;
   Muon_PFIsodBeta04 = 0;
   Muon_isHighPt = 0;
   Muon_isLoose = 0;
   Muon_isMedium = 0;
   Muon_isMedium2016 = 0;
   Muon_isSoft = 0;
   Muon_isTight = 0;
   Muon_pdgId = 0;
   Muon_rand1 = 0;
   Muon_rand2 = 0;
   Muon_relIso03 = 0;
   Muon_relIso04 = 0;
   Muon_trackerLayers = 0;
   PDFWeight_BjorkenX1 = 0;
   PDFWeight_BjorkenX2 = 0;
   PDFWeight_Cent = 0;
   PDFWeight_Cent_Down = 0;
   PDFWeight_Cent_Up = 0;
   PDFWeight_Id1 = 0;
   PDFWeight_Id2 = 0;
   PDFWeight_PDF1 = 0;
   PDFWeight_PDF2 = 0;
   PDFWeight_Q = 0;
   PDFWeight_Var1_Down = 0;
   PDFWeight_Var1_Up = 0;
   PDFWeight_Var2_Down = 0;
   PDFWeight_Var2_Up = 0;
   Photon = 0;
   Photon_SCB_Loose = 0;
   Photon_SCB_Medium = 0;
   Photon_SCB_Tight = 0;
   Photon_ChaHadIso = 0;
   Photon_NeuHadIso = 0;
   Photon_PhoIso = 0;
   Photon_WorstChargedIso = 0;
   Photon_ChaHadEffArea = 0;
   Photon_NeuHadEffArea = 0;
   Photon_PhoHadEffArea = 0;
   Photon_R9 = 0;
   Photon_HoverE = 0;
   Photon_SuperCluster_Eta = 0;
   Photon_SuperCluster_Phi = 0;
   Photon_SuperCluster_EtaWidth = 0;
   Photon_SuperCluster_PhiWidth = 0;
   Photon_SuperCluster_Energy = 0;
   Photon_SuperCluster_RawEnergy = 0;
   Photon_ElectronVeto = 0;
   Photon_Full5x5_SigmaIetaIeta = 0;
   Photon_MVANonTrig_Tight = 0;
   Rho = 0;
   Trigger_Name = 0;
   Trigger_PreScale = 0;
   Trigger_isError = 0;
   Trigger_isPass = 0;
   Trigger_isRun = 0;
   Vertex_SumPtSquare = 0;
   Vertex_X = 0;
   Vertex_X_Error = 0;
   Vertex_Y = 0;
   Vertex_Y_Error = 0;
   Vertex_Z = 0;
   Vertex_Z_Error = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Info_EventNumber", &Info_EventNumber, &b_Info_EventNumber);
   fChain->SetBranchAddress("Info_Luminosity", &Info_Luminosity, &b_Info_Luminosity);
   fChain->SetBranchAddress("Info_RunNumber", &Info_RunNumber, &b_Info_RunNumber);
   fChain->SetBranchAddress("Info_isData", &Info_isData, &b_Info_isData);
   fChain->SetBranchAddress("Channel_Idx", &Channel_Idx, &b_Channel_Idx);
   fChain->SetBranchAddress("Channel_Idx_Final", &Channel_Idx_Final, &b_Channel_Idx_Final);
   fChain->SetBranchAddress("Channel_Lepton_Count", &Channel_Lepton_Count, &b_Channel_Lepton_Count);
   fChain->SetBranchAddress("Channel_Lepton_Count_Final", &Channel_Lepton_Count_Final, &b_Channel_Lepton_Count_Final);
   fChain->SetBranchAddress("Channel_Jets", &Channel_Jets, &b_Channel_Jets);
   fChain->SetBranchAddress("Channel_Jets_Abs", &Channel_Jets_Abs, &b_Channel_Jets_Abs);
   fChain->SetBranchAddress("L1_PreFire_Central", &L1_PreFire_Central, &b_L1_PreFire_Central);
   fChain->SetBranchAddress("L1_PreFire_Up", &L1_PreFire_Up, &b_L1_PreFire_Up);
   fChain->SetBranchAddress("L1_PreFire_Down", &L1_PreFire_Down, &b_L1_PreFire_Down);
   fChain->SetBranchAddress("METFilter_Name", &METFilter_Name, &b_METFilter_Name);
   fChain->SetBranchAddress("METFilter_isError", &METFilter_isError, &b_METFilter_isError);
   fChain->SetBranchAddress("METFilter_isPass", &METFilter_isPass, &b_METFilter_isPass);
   fChain->SetBranchAddress("METFilter_isRun", &METFilter_isRun, &b_METFilter_isRun);
   fChain->SetBranchAddress("METFilterAdd_Name", &METFilterAdd_Name, &b_METFilterAdd_Name);
   fChain->SetBranchAddress("METFilterAdd_isPass", &METFilterAdd_isPass, &b_METFilterAdd_isPass);
   fChain->SetBranchAddress("Elec", &Elec, &b_Elec);
   fChain->SetBranchAddress("Elec_Charge", &Elec_Charge, &b_Elec_Charge);
   fChain->SetBranchAddress("Elec_ChargeId_GsfCtf", &Elec_ChargeId_GsfCtf, &b_Elec_ChargeId_GsfCtf);
   fChain->SetBranchAddress("Elec_ChargeId_GsfCtfPx", &Elec_ChargeId_GsfCtfPx, &b_Elec_ChargeId_GsfCtfPx);
   fChain->SetBranchAddress("Elec_ChargeId_GsfPx", &Elec_ChargeId_GsfPx, &b_Elec_ChargeId_GsfPx);
   fChain->SetBranchAddress("Elec_Charge_CtfTr", &Elec_Charge_CtfTr, &b_Elec_Charge_CtfTr);
   fChain->SetBranchAddress("Elec_Charge_GsfTr", &Elec_Charge_GsfTr, &b_Elec_Charge_GsfTr);
   fChain->SetBranchAddress("Elec_Charge_Px", &Elec_Charge_Px, &b_Elec_Charge_Px);
   fChain->SetBranchAddress("Elec_Conversion", &Elec_Conversion, &b_Elec_Conversion);
   fChain->SetBranchAddress("Elec_Count", &Elec_Count, &b_Elec_Count);
   fChain->SetBranchAddress("Elec_Id_Loose", &Elec_Id_Loose, &b_Elec_Id_Loose);
   fChain->SetBranchAddress("Elec_Id_RobustHighEnergy", &Elec_Id_RobustHighEnergy, &b_Elec_Id_RobustHighEnergy);
   fChain->SetBranchAddress("Elec_Id_RobustLoose", &Elec_Id_RobustLoose, &b_Elec_Id_RobustLoose);
   fChain->SetBranchAddress("Elec_Id_RobustTight", &Elec_Id_RobustTight, &b_Elec_Id_RobustTight);
   fChain->SetBranchAddress("Elec_Id_Tight", &Elec_Id_Tight, &b_Elec_Id_Tight);
   fChain->SetBranchAddress("Elec_Inner_Hit", &Elec_Inner_Hit, &b_Elec_Inner_Hit);
   fChain->SetBranchAddress("Elec_MVA_Medium", &Elec_MVA_Medium, &b_Elec_MVA_Medium);
   fChain->SetBranchAddress("Elec_MVA_Tight", &Elec_MVA_Tight, &b_Elec_MVA_Tight);
   fChain->SetBranchAddress("Elec_MVA_HZZ", &Elec_MVA_HZZ, &b_Elec_MVA_HZZ);
   fChain->SetBranchAddress("Elec_MVA_Values", &Elec_MVA_Values, &b_Elec_MVA_Values);
   fChain->SetBranchAddress("Elec_MVA_Categories", &Elec_MVA_Categories, &b_Elec_MVA_Categories);
   fChain->SetBranchAddress("Elec_MVA_HZZ_Values", &Elec_MVA_HZZ_Values, &b_Elec_MVA_HZZ_Values);
   fChain->SetBranchAddress("Elec_MVA_HZZ_Categories", &Elec_MVA_HZZ_Categories, &b_Elec_MVA_HZZ_Categories);
   fChain->SetBranchAddress("Elec_PFIsoRho03", &Elec_PFIsoRho03, &b_Elec_PFIsoRho03);
   fChain->SetBranchAddress("Elec_PFIsoRho04", &Elec_PFIsoRho04, &b_Elec_PFIsoRho04);
   fChain->SetBranchAddress("Elec_PFIsoValid", &Elec_PFIsoValid, &b_Elec_PFIsoValid);
   fChain->SetBranchAddress("Elec_PFIsodBeta03", &Elec_PFIsodBeta03, &b_Elec_PFIsodBeta03);
   fChain->SetBranchAddress("Elec_PFIsodBeta04", &Elec_PFIsodBeta04, &b_Elec_PFIsodBeta04);
   fChain->SetBranchAddress("Elec_SCB_Loose", &Elec_SCB_Loose, &b_Elec_SCB_Loose);
   fChain->SetBranchAddress("Elec_SCB_Medium", &Elec_SCB_Medium, &b_Elec_SCB_Medium);
   fChain->SetBranchAddress("Elec_SCB_Tight", &Elec_SCB_Tight, &b_Elec_SCB_Tight);
   fChain->SetBranchAddress("Elec_SCB_Veto", &Elec_SCB_Veto, &b_Elec_SCB_Veto);
   fChain->SetBranchAddress("Elec_SCB_dEtaIn", &Elec_SCB_dEtaIn, &b_Elec_SCB_dEtaIn);
   fChain->SetBranchAddress("Elec_SCB_dPhiIn", &Elec_SCB_dPhiIn, &b_Elec_SCB_dPhiIn);
   fChain->SetBranchAddress("Elec_SCB_hOverE", &Elec_SCB_hOverE, &b_Elec_SCB_hOverE);
   fChain->SetBranchAddress("Elec_SCB_HEEP", &Elec_SCB_HEEP, &b_Elec_SCB_HEEP);
   fChain->SetBranchAddress("Elec_SCB_ooEmooP", &Elec_SCB_ooEmooP, &b_Elec_SCB_ooEmooP);
   fChain->SetBranchAddress("Elec_SCB_sigmaIetaIeta", &Elec_SCB_sigmaIetaIeta, &b_Elec_SCB_sigmaIetaIeta);
   fChain->SetBranchAddress("Elec_ScaleUp", &Elec_ScaleUp, &b_Elec_ScaleUp);
   fChain->SetBranchAddress("Elec_ScaleDown", &Elec_ScaleDown, &b_Elec_ScaleDown);
   fChain->SetBranchAddress("Elec_SigmaUp", &Elec_SigmaUp, &b_Elec_SigmaUp);
   fChain->SetBranchAddress("Elec_SigmaDown", &Elec_SigmaDown, &b_Elec_SigmaDown);
   fChain->SetBranchAddress("Elec_ScSmUpUp", &Elec_ScSmUpUp, &b_Elec_ScSmUpUp);
   fChain->SetBranchAddress("Elec_ScSmUpDown", &Elec_ScSmUpDown, &b_Elec_ScSmUpDown);
   fChain->SetBranchAddress("Elec_ScSmDownUp", &Elec_ScSmDownUp, &b_Elec_ScSmDownUp);
   fChain->SetBranchAddress("Elec_ScSmDownDown", &Elec_ScSmDownDown, &b_Elec_ScSmDownDown);
   fChain->SetBranchAddress("Elec_Supercluster_Eta", &Elec_Supercluster_Eta, &b_Elec_Supercluster_Eta);
   fChain->SetBranchAddress("Elec_Track_CtfdXY", &Elec_Track_CtfdXY, &b_Elec_Track_CtfdXY);
   fChain->SetBranchAddress("Elec_Track_CtfdZ", &Elec_Track_CtfdZ, &b_Elec_Track_CtfdZ);
   fChain->SetBranchAddress("Elec_Track_GsfdXY", &Elec_Track_GsfdXY, &b_Elec_Track_GsfdXY);
   fChain->SetBranchAddress("Elec_Track_GsfdZ", &Elec_Track_GsfdZ, &b_Elec_Track_GsfdZ);
   fChain->SetBranchAddress("Elec_pdgId", &Elec_pdgId, &b_Elec_pdgId);
   fChain->SetBranchAddress("Elec_relIso03", &Elec_relIso03, &b_Elec_relIso03);
   fChain->SetBranchAddress("Elec_relIso04", &Elec_relIso04, &b_Elec_relIso04);
   fChain->SetBranchAddress("Elec_Rnd", &Elec_Rnd, &b_Elec_Rnd);
   fChain->SetBranchAddress("Filter_Greedy_Muon", &Filter_Greedy_Muon, &b_Filter_Greedy_Muon);
   fChain->SetBranchAddress("Filter_Inconsistent_MuonPt", &Filter_Inconsistent_MuonPt, &b_Filter_Inconsistent_MuonPt);
   fChain->SetBranchAddress("Filter_PFReco_Muon", &Filter_PFReco_Muon, &b_Filter_PFReco_Muon);
   fChain->SetBranchAddress("Filter_PV", &Filter_PV, &b_Filter_PV);
   fChain->SetBranchAddress("GenJet", &GenJet, &b_GenJet);
   fChain->SetBranchAddress("GenJet_Count", &GenJet_Count, &b_GenJet_Count);
   fChain->SetBranchAddress("GenJet_ECalEnergy", &GenJet_ECalEnergy, &b_GenJet_ECalEnergy);
   fChain->SetBranchAddress("GenJet_HCalEnergy", &GenJet_HCalEnergy, &b_GenJet_HCalEnergy);
   fChain->SetBranchAddress("GenMET", &GenMET, &b_GenMET);
   fChain->SetBranchAddress("GenMET_Count", &GenMET_Count, &b_GenMET_Count);
   fChain->SetBranchAddress("GenPar", &GenPar, &b_GenPar);
   fChain->SetBranchAddress("GenPar_Count", &GenPar_Count, &b_GenPar_Count);
   fChain->SetBranchAddress("GenPar_Dau1_Idx", &GenPar_Dau1_Idx, &b_GenPar_Dau1_Idx);
   fChain->SetBranchAddress("GenPar_Dau2_Idx", &GenPar_Dau2_Idx, &b_GenPar_Dau2_Idx);
   fChain->SetBranchAddress("GenPar_Dau_Counter", &GenPar_Dau_Counter, &b_GenPar_Dau_Counter);
   fChain->SetBranchAddress("GenPar_Idx", &GenPar_Idx, &b_GenPar_Idx);
   fChain->SetBranchAddress("GenPar_Mom1_Idx", &GenPar_Mom1_Idx, &b_GenPar_Mom1_Idx);
   fChain->SetBranchAddress("GenPar_Mom2_Idx", &GenPar_Mom2_Idx, &b_GenPar_Mom2_Idx);
   fChain->SetBranchAddress("GenPar_Mom_Counter", &GenPar_Mom_Counter, &b_GenPar_Mom_Counter);
   fChain->SetBranchAddress("GenPar_Status", &GenPar_Status, &b_GenPar_Status);
   fChain->SetBranchAddress("GenPar_pdgId", &GenPar_pdgId, &b_GenPar_pdgId);
   fChain->SetBranchAddress("Gen_EventWeight", &Gen_EventWeight, &b_Gen_EventWeight);
   fChain->SetBranchAddress("GenTop", &GenTop, &b_GenTop);
   fChain->SetBranchAddress("GenAnTop", &GenAnTop, &b_GenAnTop);
   fChain->SetBranchAddress("GenBJet", &GenBJet, &b_GenBJet);
   fChain->SetBranchAddress("GenBJet_Count", &GenBJet_Count, &b_GenBJet_Count);
   fChain->SetBranchAddress("GenBHad", &GenBHad, &b_GenBHad);
   fChain->SetBranchAddress("GenBHad_Count", &GenBHad_Count, &b_GenBHad_Count);
   fChain->SetBranchAddress("GenBHad_pdgId", &GenBHad_pdgId, &b_GenBHad_pdgId);
   fChain->SetBranchAddress("GenBHad_Flavour", &GenBHad_Flavour, &b_GenBHad_Flavour);
   fChain->SetBranchAddress("GenBHad_FromTopWeakDecay", &GenBHad_FromTopWeakDecay, &b_GenBHad_FromTopWeakDecay);
   fChain->SetBranchAddress("Jet", &Jet, &b_Jet);
   fChain->SetBranchAddress("Jet_Charge", &Jet_Charge, &b_Jet_Charge);
   fChain->SetBranchAddress("Jet_Count", &Jet_Count, &b_Jet_Count);
   fChain->SetBranchAddress("Jet_EnShiftedDown", &Jet_EnShiftedDown, &b_Jet_EnShiftedDown);
   fChain->SetBranchAddress("Jet_EnShiftedUp", &Jet_EnShiftedUp, &b_Jet_EnShiftedUp);
   fChain->SetBranchAddress("Jet_PartonFlavour", &Jet_PartonFlavour, &b_Jet_PartonFlavour);
   fChain->SetBranchAddress("Jet_HadronFlavour", &Jet_HadronFlavour, &b_Jet_HadronFlavour);
   fChain->SetBranchAddress("Jet_PFId", &Jet_PFId, &b_Jet_PFId);
   fChain->SetBranchAddress("Jet_PFIdVeto", &Jet_PFIdVeto, &b_Jet_PFIdVeto);
   fChain->SetBranchAddress("Jet_PhiResolution_MC", &Jet_PhiResolution_MC, &b_Jet_PhiResolution_MC);
   fChain->SetBranchAddress("Jet_PhiResolution_DATA", &Jet_PhiResolution_DATA, &b_Jet_PhiResolution_DATA);
   fChain->SetBranchAddress("Jet_EnergyResolution_MC", &Jet_EnergyResolution_MC, &b_Jet_EnergyResolution_MC);
   fChain->SetBranchAddress("Jet_EnergyResolution_DATA", &Jet_EnergyResolution_DATA, &b_Jet_EnergyResolution_DATA);
   fChain->SetBranchAddress("Jet_EnergyResolution_SF", &Jet_EnergyResolution_SF, &b_Jet_EnergyResolution_SF);
   fChain->SetBranchAddress("Jet_EnergyResolution_SFDown", &Jet_EnergyResolution_SFDown, &b_Jet_EnergyResolution_SFDown);
   fChain->SetBranchAddress("Jet_EnergyResolution_SFUp", &Jet_EnergyResolution_SFUp, &b_Jet_EnergyResolution_SFUp);
   fChain->SetBranchAddress("Jet_PileUpId", &Jet_PileUpId, &b_Jet_PileUpId);
   fChain->SetBranchAddress("Jet_PileUpMVA", &Jet_PileUpMVA, &b_Jet_PileUpMVA);
   fChain->SetBranchAddress("Jet_bDisc", &Jet_bDisc, &b_Jet_bDisc);
   fChain->SetBranchAddress("Jet_bDisc_pfCombinedInclusiveSecondaryVertexV2BJetTags", &Jet_bDisc_pfCombinedInclusiveSecondaryVertexV2BJetTags, &b_Jet_bDisc_pfCombinedInclusiveSecondaryVertexV2BJetTags);
   fChain->SetBranchAddress("Jet_bDisc_softPFMuonBJetTags", &Jet_bDisc_softPFMuonBJetTags, &b_Jet_bDisc_softPFMuonBJetTags);
   fChain->SetBranchAddress("Jet_bDisc_softPFMuonByIP3dBJetTags", &Jet_bDisc_softPFMuonByIP3dBJetTags, &b_Jet_bDisc_softPFMuonByIP3dBJetTags);
   fChain->SetBranchAddress("Jet_bDisc_softPFElectronByPtBJetTags", &Jet_bDisc_softPFElectronByPtBJetTags, &b_Jet_bDisc_softPFElectronByPtBJetTags);
   fChain->SetBranchAddress("Jet_bDisc_softPFElectronBJetTags", &Jet_bDisc_softPFElectronBJetTags, &b_Jet_bDisc_softPFElectronBJetTags);
   fChain->SetBranchAddress("Jet_bDisc_softPFMuonByPtBJetTags", &Jet_bDisc_softPFMuonByPtBJetTags, &b_Jet_bDisc_softPFMuonByPtBJetTags);
   fChain->SetBranchAddress("Jet_bDisc_softPFElectronByIP3dBJetTags", &Jet_bDisc_softPFElectronByIP3dBJetTags, &b_Jet_bDisc_softPFElectronByIP3dBJetTags);
   fChain->SetBranchAddress("Jet_bDisc_softPFMuonByIP2dBJetTags", &Jet_bDisc_softPFMuonByIP2dBJetTags, &b_Jet_bDisc_softPFMuonByIP2dBJetTags);
   fChain->SetBranchAddress("Jet_bDisc_softPFElectronByIP2dBJetTags", &Jet_bDisc_softPFElectronByIP2dBJetTags, &b_Jet_bDisc_softPFElectronByIP2dBJetTags);
   fChain->SetBranchAddress("Jet_isJet", &Jet_isJet, &b_Jet_isJet);
   fChain->SetBranchAddress("LHE_Central", &LHE_Central, &b_LHE_Central);
   fChain->SetBranchAddress("LHE_Id", &LHE_Id, &b_LHE_Id);
   fChain->SetBranchAddress("LHE_Weight", &LHE_Weight, &b_LHE_Weight);
   fChain->SetBranchAddress("Frag_Cen_Weight", &Frag_Cen_Weight, &b_Frag_Cen_Weight);
   fChain->SetBranchAddress("Frag_Up_Weight", &Frag_Up_Weight, &b_Frag_Up_Weight);
   fChain->SetBranchAddress("Frag_Down_Weight", &Frag_Down_Weight, &b_Frag_Down_Weight);
   fChain->SetBranchAddress("Frag_Peterson_Weight", &Frag_Peterson_Weight, &b_Frag_Peterson_Weight);
   fChain->SetBranchAddress("Semilep_BrUp_Weight", &Semilep_BrUp_Weight, &b_Semilep_BrUp_Weight);
   fChain->SetBranchAddress("Semilep_BrDown_Weight", &Semilep_BrDown_Weight, &b_Semilep_BrDown_Weight);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("MET_Significance", &MET_Significance, &b_MET_Significance);
   fChain->SetBranchAddress("MET_JetEnShiftedUp_PT", &MET_JetEnShiftedUp_PT, &b_MET_JetEnShiftedUp_PT);
   fChain->SetBranchAddress("MET_JetEnShiftedUp_Phi", &MET_JetEnShiftedUp_Phi, &b_MET_JetEnShiftedUp_Phi);
   fChain->SetBranchAddress("MET_JetEnShiftedDown_PT", &MET_JetEnShiftedDown_PT, &b_MET_JetEnShiftedDown_PT);
   fChain->SetBranchAddress("MET_JetEnShiftedDown_Phi", &MET_JetEnShiftedDown_Phi, &b_MET_JetEnShiftedDown_Phi);
   fChain->SetBranchAddress("MET_MuonEnShiftedUp_PT", &MET_MuonEnShiftedUp_PT, &b_MET_MuonEnShiftedUp_PT);
   fChain->SetBranchAddress("MET_MuonEnShiftedUp_Phi", &MET_MuonEnShiftedUp_Phi, &b_MET_MuonEnShiftedUp_Phi);
   fChain->SetBranchAddress("MET_MuonEnShiftedDown_PT", &MET_MuonEnShiftedDown_PT, &b_MET_MuonEnShiftedDown_PT);
   fChain->SetBranchAddress("MET_MuonEnShiftedDown_Phi", &MET_MuonEnShiftedDown_Phi, &b_MET_MuonEnShiftedDown_Phi);
   fChain->SetBranchAddress("MET_ElectronEnShiftedUp_PT", &MET_ElectronEnShiftedUp_PT, &b_MET_ElectronEnShiftedUp_PT);
   fChain->SetBranchAddress("MET_ElectronEnShiftedUp_Phi", &MET_ElectronEnShiftedUp_Phi, &b_MET_ElectronEnShiftedUp_Phi);
   fChain->SetBranchAddress("MET_ElectronEnShiftedDown_PT", &MET_ElectronEnShiftedDown_PT, &b_MET_ElectronEnShiftedDown_PT);
   fChain->SetBranchAddress("MET_ElectronEnShiftedDown_Phi", &MET_ElectronEnShiftedDown_Phi, &b_MET_ElectronEnShiftedDown_Phi);
   fChain->SetBranchAddress("MET_UnclusteredEnShiftedUp_PT", &MET_UnclusteredEnShiftedUp_PT, &b_MET_UnclusteredEnShiftedUp_PT);
   fChain->SetBranchAddress("MET_UnclusteredEnShiftedUp_Phi", &MET_UnclusteredEnShiftedUp_Phi, &b_MET_UnclusteredEnShiftedUp_Phi);
   fChain->SetBranchAddress("MET_UnclusteredEnShiftedDown_PT", &MET_UnclusteredEnShiftedDown_PT, &b_MET_UnclusteredEnShiftedDown_PT);
   fChain->SetBranchAddress("MET_UnclusteredEnShiftedDown_Phi", &MET_UnclusteredEnShiftedDown_Phi, &b_MET_UnclusteredEnShiftedDown_Phi);
   fChain->SetBranchAddress("MET_JetResShiftedUp_PT", &MET_JetResShiftedUp_PT, &b_MET_JetResShiftedUp_PT);
   fChain->SetBranchAddress("MET_JetResShiftedUp_Phi", &MET_JetResShiftedUp_Phi, &b_MET_JetResShiftedUp_Phi);
   fChain->SetBranchAddress("MET_JetResShiftedDown_PT", &MET_JetResShiftedDown_PT, &b_MET_JetResShiftedDown_PT);
   fChain->SetBranchAddress("MET_JetResShiftedDown_Phi", &MET_JetResShiftedDown_Phi, &b_MET_JetResShiftedDown_Phi);
   fChain->SetBranchAddress("METEGClean", &METEGClean, &b_METEGClean);
   fChain->SetBranchAddress("METEGClean_Significance", &METEGClean_Significance, &b_METEGClean_Significance);
   fChain->SetBranchAddress("METEGClean_JetEnShiftedUp_PT", &METEGClean_JetEnShiftedUp_PT, &b_METEGClean_JetEnShiftedUp_PT);
   fChain->SetBranchAddress("METEGClean_JetEnShiftedUp_Phi", &METEGClean_JetEnShiftedUp_Phi, &b_METEGClean_JetEnShiftedUp_Phi);
   fChain->SetBranchAddress("METEGClean_JetEnShiftedDown_PT", &METEGClean_JetEnShiftedDown_PT, &b_METEGClean_JetEnShiftedDown_PT);
   fChain->SetBranchAddress("METEGClean_JetEnShiftedDown_Phi", &METEGClean_JetEnShiftedDown_Phi, &b_METEGClean_JetEnShiftedDown_Phi);
   fChain->SetBranchAddress("METEGClean_MuonEnShiftedUp_PT", &METEGClean_MuonEnShiftedUp_PT, &b_METEGClean_MuonEnShiftedUp_PT);
   fChain->SetBranchAddress("METEGClean_MuonEnShiftedUp_Phi", &METEGClean_MuonEnShiftedUp_Phi, &b_METEGClean_MuonEnShiftedUp_Phi);
   fChain->SetBranchAddress("METEGClean_MuonEnShiftedDown_PT", &METEGClean_MuonEnShiftedDown_PT, &b_METEGClean_MuonEnShiftedDown_PT);
   fChain->SetBranchAddress("METEGClean_MuonEnShiftedDown_Phi", &METEGClean_MuonEnShiftedDown_Phi, &b_METEGClean_MuonEnShiftedDown_Phi);
   fChain->SetBranchAddress("METEGClean_ElectronEnShiftedUp_PT", &METEGClean_ElectronEnShiftedUp_PT, &b_METEGClean_ElectronEnShiftedUp_PT);
   fChain->SetBranchAddress("METEGClean_ElectronEnShiftedUp_Phi", &METEGClean_ElectronEnShiftedUp_Phi, &b_METEGClean_ElectronEnShiftedUp_Phi);
   fChain->SetBranchAddress("METEGClean_ElectronEnShiftedDown_PT", &METEGClean_ElectronEnShiftedDown_PT, &b_METEGClean_ElectronEnShiftedDown_PT);
   fChain->SetBranchAddress("METEGClean_ElectronEnShiftedDown_Phi", &METEGClean_ElectronEnShiftedDown_Phi, &b_METEGClean_ElectronEnShiftedDown_Phi);
   fChain->SetBranchAddress("METEGClean_UnclusteredEnShiftedUp_PT", &METEGClean_UnclusteredEnShiftedUp_PT, &b_METEGClean_UnclusteredEnShiftedUp_PT);
   fChain->SetBranchAddress("METEGClean_UnclusteredEnShiftedUp_Phi", &METEGClean_UnclusteredEnShiftedUp_Phi, &b_METEGClean_UnclusteredEnShiftedUp_Phi);
   fChain->SetBranchAddress("METEGClean_UnclusteredEnShiftedDown_PT", &METEGClean_UnclusteredEnShiftedDown_PT, &b_METEGClean_UnclusteredEnShiftedDown_PT);
   fChain->SetBranchAddress("METEGClean_UnclusteredEnShiftedDown_Phi", &METEGClean_UnclusteredEnShiftedDown_Phi, &b_METEGClean_UnclusteredEnShiftedDown_Phi);
   fChain->SetBranchAddress("METEGClean_JetResShiftedUp_PT", &METEGClean_JetResShiftedUp_PT, &b_METEGClean_JetResShiftedUp_PT);
   fChain->SetBranchAddress("METEGClean_JetResShiftedUp_Phi", &METEGClean_JetResShiftedUp_Phi, &b_METEGClean_JetResShiftedUp_Phi);
   fChain->SetBranchAddress("METEGClean_JetResShiftedDown_PT", &METEGClean_JetResShiftedDown_PT, &b_METEGClean_JetResShiftedDown_PT);
   fChain->SetBranchAddress("METEGClean_JetResShiftedDown_Phi", &METEGClean_JetResShiftedDown_Phi, &b_METEGClean_JetResShiftedDown_Phi);
   fChain->SetBranchAddress("METMUEGClean", &METMUEGClean, &b_METMUEGClean);
   fChain->SetBranchAddress("METMUEGClean_Significance", &METMUEGClean_Significance, &b_METMUEGClean_Significance);
   fChain->SetBranchAddress("METMUEGClean_JetEnShiftedUp_PT", &METMUEGClean_JetEnShiftedUp_PT, &b_METMUEGClean_JetEnShiftedUp_PT);
   fChain->SetBranchAddress("METMUEGClean_JetEnShiftedUp_Phi", &METMUEGClean_JetEnShiftedUp_Phi, &b_METMUEGClean_JetEnShiftedUp_Phi);
   fChain->SetBranchAddress("METMUEGClean_JetEnShiftedDown_PT", &METMUEGClean_JetEnShiftedDown_PT, &b_METMUEGClean_JetEnShiftedDown_PT);
   fChain->SetBranchAddress("METMUEGClean_JetEnShiftedDown_Phi", &METMUEGClean_JetEnShiftedDown_Phi, &b_METMUEGClean_JetEnShiftedDown_Phi);
   fChain->SetBranchAddress("METMUEGClean_MuonEnShiftedUp_PT", &METMUEGClean_MuonEnShiftedUp_PT, &b_METMUEGClean_MuonEnShiftedUp_PT);
   fChain->SetBranchAddress("METMUEGClean_MuonEnShiftedUp_Phi", &METMUEGClean_MuonEnShiftedUp_Phi, &b_METMUEGClean_MuonEnShiftedUp_Phi);
   fChain->SetBranchAddress("METMUEGClean_MuonEnShiftedDown_PT", &METMUEGClean_MuonEnShiftedDown_PT, &b_METMUEGClean_MuonEnShiftedDown_PT);
   fChain->SetBranchAddress("METMUEGClean_MuonEnShiftedDown_Phi", &METMUEGClean_MuonEnShiftedDown_Phi, &b_METMUEGClean_MuonEnShiftedDown_Phi);
   fChain->SetBranchAddress("METMUEGClean_ElectronEnShiftedUp_PT", &METMUEGClean_ElectronEnShiftedUp_PT, &b_METMUEGClean_ElectronEnShiftedUp_PT);
   fChain->SetBranchAddress("METMUEGClean_ElectronEnShiftedUp_Phi", &METMUEGClean_ElectronEnShiftedUp_Phi, &b_METMUEGClean_ElectronEnShiftedUp_Phi);
   fChain->SetBranchAddress("METMUEGClean_ElectronEnShiftedDown_PT", &METMUEGClean_ElectronEnShiftedDown_PT, &b_METMUEGClean_ElectronEnShiftedDown_PT);
   fChain->SetBranchAddress("METMUEGClean_ElectronEnShiftedDown_Phi", &METMUEGClean_ElectronEnShiftedDown_Phi, &b_METMUEGClean_ElectronEnShiftedDown_Phi);
   fChain->SetBranchAddress("METMUEGClean_UnclusteredEnShiftedUp_PT", &METMUEGClean_UnclusteredEnShiftedUp_PT, &b_METMUEGClean_UnclusteredEnShiftedUp_PT);
   fChain->SetBranchAddress("METMUEGClean_UnclusteredEnShiftedUp_Phi", &METMUEGClean_UnclusteredEnShiftedUp_Phi, &b_METMUEGClean_UnclusteredEnShiftedUp_Phi);
   fChain->SetBranchAddress("METMUEGClean_UnclusteredEnShiftedDown_PT", &METMUEGClean_UnclusteredEnShiftedDown_PT, &b_METMUEGClean_UnclusteredEnShiftedDown_PT);
   fChain->SetBranchAddress("METMUEGClean_UnclusteredEnShiftedDown_Phi", &METMUEGClean_UnclusteredEnShiftedDown_Phi, &b_METMUEGClean_UnclusteredEnShiftedDown_Phi);
   fChain->SetBranchAddress("METMUEGClean_JetResShiftedUp_PT", &METMUEGClean_JetResShiftedUp_PT, &b_METMUEGClean_JetResShiftedUp_PT);
   fChain->SetBranchAddress("METMUEGClean_JetResShiftedUp_Phi", &METMUEGClean_JetResShiftedUp_Phi, &b_METMUEGClean_JetResShiftedUp_Phi);
   fChain->SetBranchAddress("METMUEGClean_JetResShiftedDown_PT", &METMUEGClean_JetResShiftedDown_PT, &b_METMUEGClean_JetResShiftedDown_PT);
   fChain->SetBranchAddress("METMUEGClean_JetResShiftedDown_Phi", &METMUEGClean_JetResShiftedDown_Phi, &b_METMUEGClean_JetResShiftedDown_Phi);
   fChain->SetBranchAddress("METMUEGCleanCor", &METMUEGCleanCor, &b_METMUEGCleanCor);
   fChain->SetBranchAddress("METMUEGCleanCor_Significance", &METMUEGCleanCor_Significance, &b_METMUEGCleanCor_Significance);
   fChain->SetBranchAddress("METMUEGCleanCor_JetEnShiftedUp_PT", &METMUEGCleanCor_JetEnShiftedUp_PT, &b_METMUEGCleanCor_JetEnShiftedUp_PT);
   fChain->SetBranchAddress("METMUEGCleanCor_JetEnShiftedUp_Phi", &METMUEGCleanCor_JetEnShiftedUp_Phi, &b_METMUEGCleanCor_JetEnShiftedUp_Phi);
   fChain->SetBranchAddress("METMUEGCleanCor_JetEnShiftedDown_PT", &METMUEGCleanCor_JetEnShiftedDown_PT, &b_METMUEGCleanCor_JetEnShiftedDown_PT);
   fChain->SetBranchAddress("METMUEGCleanCor_JetEnShiftedDown_Phi", &METMUEGCleanCor_JetEnShiftedDown_Phi, &b_METMUEGCleanCor_JetEnShiftedDown_Phi);
   fChain->SetBranchAddress("METMUEGCleanCor_MuonEnShiftedUp_PT", &METMUEGCleanCor_MuonEnShiftedUp_PT, &b_METMUEGCleanCor_MuonEnShiftedUp_PT);
   fChain->SetBranchAddress("METMUEGCleanCor_MuonEnShiftedUp_Phi", &METMUEGCleanCor_MuonEnShiftedUp_Phi, &b_METMUEGCleanCor_MuonEnShiftedUp_Phi);
   fChain->SetBranchAddress("METMUEGCleanCor_MuonEnShiftedDown_PT", &METMUEGCleanCor_MuonEnShiftedDown_PT, &b_METMUEGCleanCor_MuonEnShiftedDown_PT);
   fChain->SetBranchAddress("METMUEGCleanCor_MuonEnShiftedDown_Phi", &METMUEGCleanCor_MuonEnShiftedDown_Phi, &b_METMUEGCleanCor_MuonEnShiftedDown_Phi);
   fChain->SetBranchAddress("METMUEGCleanCor_ElectronEnShiftedUp_PT", &METMUEGCleanCor_ElectronEnShiftedUp_PT, &b_METMUEGCleanCor_ElectronEnShiftedUp_PT);
   fChain->SetBranchAddress("METMUEGCleanCor_ElectronEnShiftedUp_Phi", &METMUEGCleanCor_ElectronEnShiftedUp_Phi, &b_METMUEGCleanCor_ElectronEnShiftedUp_Phi);
   fChain->SetBranchAddress("METMUEGCleanCor_ElectronEnShiftedDown_PT", &METMUEGCleanCor_ElectronEnShiftedDown_PT, &b_METMUEGCleanCor_ElectronEnShiftedDown_PT);
   fChain->SetBranchAddress("METMUEGCleanCor_ElectronEnShiftedDown_Phi", &METMUEGCleanCor_ElectronEnShiftedDown_Phi, &b_METMUEGCleanCor_ElectronEnShiftedDown_Phi);
   fChain->SetBranchAddress("METMUEGCleanCor_UnclusteredEnShiftedUp_PT", &METMUEGCleanCor_UnclusteredEnShiftedUp_PT, &b_METMUEGCleanCor_UnclusteredEnShiftedUp_PT);
   fChain->SetBranchAddress("METMUEGCleanCor_UnclusteredEnShiftedUp_Phi", &METMUEGCleanCor_UnclusteredEnShiftedUp_Phi, &b_METMUEGCleanCor_UnclusteredEnShiftedUp_Phi);
   fChain->SetBranchAddress("METMUEGCleanCor_UnclusteredEnShiftedDown_PT", &METMUEGCleanCor_UnclusteredEnShiftedDown_PT, &b_METMUEGCleanCor_UnclusteredEnShiftedDown_PT);
   fChain->SetBranchAddress("METMUEGCleanCor_UnclusteredEnShiftedDown_Phi", &METMUEGCleanCor_UnclusteredEnShiftedDown_Phi, &b_METMUEGCleanCor_UnclusteredEnShiftedDown_Phi);
   fChain->SetBranchAddress("METMUEGCleanCor_JetResShiftedUp_PT", &METMUEGCleanCor_JetResShiftedUp_PT, &b_METMUEGCleanCor_JetResShiftedUp_PT);
   fChain->SetBranchAddress("METMUEGCleanCor_JetResShiftedUp_Phi", &METMUEGCleanCor_JetResShiftedUp_Phi, &b_METMUEGCleanCor_JetResShiftedUp_Phi);
   fChain->SetBranchAddress("METMUEGCleanCor_JetResShiftedDown_PT", &METMUEGCleanCor_JetResShiftedDown_PT, &b_METMUEGCleanCor_JetResShiftedDown_PT);
   fChain->SetBranchAddress("METMUEGCleanCor_JetResShiftedDown_Phi", &METMUEGCleanCor_JetResShiftedDown_Phi, &b_METMUEGCleanCor_JetResShiftedDown_Phi);
   fChain->SetBranchAddress("METMUCleanCor", &METMUCleanCor, &b_METMUCleanCor);
   fChain->SetBranchAddress("METMUCleanCor_Significance", &METMUCleanCor_Significance, &b_METMUCleanCor_Significance);
   fChain->SetBranchAddress("METMUCleanCor_JetEnShiftedUp_PT", &METMUCleanCor_JetEnShiftedUp_PT, &b_METMUCleanCor_JetEnShiftedUp_PT);
   fChain->SetBranchAddress("METMUCleanCor_JetEnShiftedUp_Phi", &METMUCleanCor_JetEnShiftedUp_Phi, &b_METMUCleanCor_JetEnShiftedUp_Phi);
   fChain->SetBranchAddress("METMUCleanCor_JetEnShiftedDown_PT", &METMUCleanCor_JetEnShiftedDown_PT, &b_METMUCleanCor_JetEnShiftedDown_PT);
   fChain->SetBranchAddress("METMUCleanCor_JetEnShiftedDown_Phi", &METMUCleanCor_JetEnShiftedDown_Phi, &b_METMUCleanCor_JetEnShiftedDown_Phi);
   fChain->SetBranchAddress("METMUCleanCor_MuonEnShiftedUp_PT", &METMUCleanCor_MuonEnShiftedUp_PT, &b_METMUCleanCor_MuonEnShiftedUp_PT);
   fChain->SetBranchAddress("METMUCleanCor_MuonEnShiftedUp_Phi", &METMUCleanCor_MuonEnShiftedUp_Phi, &b_METMUCleanCor_MuonEnShiftedUp_Phi);
   fChain->SetBranchAddress("METMUCleanCor_MuonEnShiftedDown_PT", &METMUCleanCor_MuonEnShiftedDown_PT, &b_METMUCleanCor_MuonEnShiftedDown_PT);
   fChain->SetBranchAddress("METMUCleanCor_MuonEnShiftedDown_Phi", &METMUCleanCor_MuonEnShiftedDown_Phi, &b_METMUCleanCor_MuonEnShiftedDown_Phi);
   fChain->SetBranchAddress("METMUCleanCor_ElectronEnShiftedUp_PT", &METMUCleanCor_ElectronEnShiftedUp_PT, &b_METMUCleanCor_ElectronEnShiftedUp_PT);
   fChain->SetBranchAddress("METMUCleanCor_ElectronEnShiftedUp_Phi", &METMUCleanCor_ElectronEnShiftedUp_Phi, &b_METMUCleanCor_ElectronEnShiftedUp_Phi);
   fChain->SetBranchAddress("METMUCleanCor_ElectronEnShiftedDown_PT", &METMUCleanCor_ElectronEnShiftedDown_PT, &b_METMUCleanCor_ElectronEnShiftedDown_PT);
   fChain->SetBranchAddress("METMUCleanCor_ElectronEnShiftedDown_Phi", &METMUCleanCor_ElectronEnShiftedDown_Phi, &b_METMUCleanCor_ElectronEnShiftedDown_Phi);
   fChain->SetBranchAddress("METMUCleanCor_UnclusteredEnShiftedUp_PT", &METMUCleanCor_UnclusteredEnShiftedUp_PT, &b_METMUCleanCor_UnclusteredEnShiftedUp_PT);
   fChain->SetBranchAddress("METMUCleanCor_UnclusteredEnShiftedUp_Phi", &METMUCleanCor_UnclusteredEnShiftedUp_Phi, &b_METMUCleanCor_UnclusteredEnShiftedUp_Phi);
   fChain->SetBranchAddress("METMUCleanCor_UnclusteredEnShiftedDown_PT", &METMUCleanCor_UnclusteredEnShiftedDown_PT, &b_METMUCleanCor_UnclusteredEnShiftedDown_PT);
   fChain->SetBranchAddress("METMUCleanCor_UnclusteredEnShiftedDown_Phi", &METMUCleanCor_UnclusteredEnShiftedDown_Phi, &b_METMUCleanCor_UnclusteredEnShiftedDown_Phi);
   fChain->SetBranchAddress("METMUCleanCor_JetResShiftedUp_PT", &METMUCleanCor_JetResShiftedUp_PT, &b_METMUCleanCor_JetResShiftedUp_PT);
   fChain->SetBranchAddress("METMUCleanCor_JetResShiftedUp_Phi", &METMUCleanCor_JetResShiftedUp_Phi, &b_METMUCleanCor_JetResShiftedUp_Phi);
   fChain->SetBranchAddress("METMUCleanCor_JetResShiftedDown_PT", &METMUCleanCor_JetResShiftedDown_PT, &b_METMUCleanCor_JetResShiftedDown_PT);
   fChain->SetBranchAddress("METMUCleanCor_JetResShiftedDown_Phi", &METMUCleanCor_JetResShiftedDown_Phi, &b_METMUCleanCor_JetResShiftedDown_Phi);
   fChain->SetBranchAddress("METUnCor", &METUnCor, &b_METUnCor);
   fChain->SetBranchAddress("METUnCor_Significance", &METUnCor_Significance, &b_METUnCor_Significance);
   fChain->SetBranchAddress("METUnCor_JetEnShiftedUp_PT", &METUnCor_JetEnShiftedUp_PT, &b_METUnCor_JetEnShiftedUp_PT);
   fChain->SetBranchAddress("METUnCor_JetEnShiftedUp_Phi", &METUnCor_JetEnShiftedUp_Phi, &b_METUnCor_JetEnShiftedUp_Phi);
   fChain->SetBranchAddress("METUnCor_JetEnShiftedDown_PT", &METUnCor_JetEnShiftedDown_PT, &b_METUnCor_JetEnShiftedDown_PT);
   fChain->SetBranchAddress("METUnCor_JetEnShiftedDown_Phi", &METUnCor_JetEnShiftedDown_Phi, &b_METUnCor_JetEnShiftedDown_Phi);
   fChain->SetBranchAddress("METUnCor_MuonEnShiftedUp_PT", &METUnCor_MuonEnShiftedUp_PT, &b_METUnCor_MuonEnShiftedUp_PT);
   fChain->SetBranchAddress("METUnCor_MuonEnShiftedUp_Phi", &METUnCor_MuonEnShiftedUp_Phi, &b_METUnCor_MuonEnShiftedUp_Phi);
   fChain->SetBranchAddress("METUnCor_MuonEnShiftedDown_PT", &METUnCor_MuonEnShiftedDown_PT, &b_METUnCor_MuonEnShiftedDown_PT);
   fChain->SetBranchAddress("METUnCor_MuonEnShiftedDown_Phi", &METUnCor_MuonEnShiftedDown_Phi, &b_METUnCor_MuonEnShiftedDown_Phi);
   fChain->SetBranchAddress("METUnCor_ElectronEnShiftedUp_PT", &METUnCor_ElectronEnShiftedUp_PT, &b_METUnCor_ElectronEnShiftedUp_PT);
   fChain->SetBranchAddress("METUnCor_ElectronEnShiftedUp_Phi", &METUnCor_ElectronEnShiftedUp_Phi, &b_METUnCor_ElectronEnShiftedUp_Phi);
   fChain->SetBranchAddress("METUnCor_ElectronEnShiftedDown_PT", &METUnCor_ElectronEnShiftedDown_PT, &b_METUnCor_ElectronEnShiftedDown_PT);
   fChain->SetBranchAddress("METUnCor_ElectronEnShiftedDown_Phi", &METUnCor_ElectronEnShiftedDown_Phi, &b_METUnCor_ElectronEnShiftedDown_Phi);
   fChain->SetBranchAddress("METUnCor_UnclusteredEnShiftedUp_PT", &METUnCor_UnclusteredEnShiftedUp_PT, &b_METUnCor_UnclusteredEnShiftedUp_PT);
   fChain->SetBranchAddress("METUnCor_UnclusteredEnShiftedUp_Phi", &METUnCor_UnclusteredEnShiftedUp_Phi, &b_METUnCor_UnclusteredEnShiftedUp_Phi);
   fChain->SetBranchAddress("METUnCor_UnclusteredEnShiftedDown_PT", &METUnCor_UnclusteredEnShiftedDown_PT, &b_METUnCor_UnclusteredEnShiftedDown_PT);
   fChain->SetBranchAddress("METUnCor_UnclusteredEnShiftedDown_Phi", &METUnCor_UnclusteredEnShiftedDown_Phi, &b_METUnCor_UnclusteredEnShiftedDown_Phi);
   fChain->SetBranchAddress("METUnCor_JetResShiftedUp_PT", &METUnCor_JetResShiftedUp_PT, &b_METUnCor_JetResShiftedUp_PT);
   fChain->SetBranchAddress("METUnCor_JetResShiftedUp_Phi", &METUnCor_JetResShiftedUp_Phi, &b_METUnCor_JetResShiftedUp_Phi);
   fChain->SetBranchAddress("METUnCor_JetResShiftedDown_PT", &METUnCor_JetResShiftedDown_PT, &b_METUnCor_JetResShiftedDown_PT);
   fChain->SetBranchAddress("METUnCor_JetResShiftedDown_Phi", &METUnCor_JetResShiftedDown_Phi, &b_METUnCor_JetResShiftedDown_Phi);
   fChain->SetBranchAddress("Muon", &Muon, &b_Muon);
   fChain->SetBranchAddress("GenMuon", &GenMuon, &b_GenMuon);
   fChain->SetBranchAddress("Muon_Charge", &Muon_Charge, &b_Muon_Charge);
   fChain->SetBranchAddress("Muon_Count", &Muon_Count, &b_Muon_Count);
   fChain->SetBranchAddress("Muon_PFIsodBeta03", &Muon_PFIsodBeta03, &b_Muon_PFIsodBeta03);
   fChain->SetBranchAddress("Muon_PFIsodBeta04", &Muon_PFIsodBeta04, &b_Muon_PFIsodBeta04);
   fChain->SetBranchAddress("Muon_isHighPt", &Muon_isHighPt, &b_Muon_isHighPt);
   fChain->SetBranchAddress("Muon_isLoose", &Muon_isLoose, &b_Muon_isLoose);
   fChain->SetBranchAddress("Muon_isMedium", &Muon_isMedium, &b_Muon_isMedium);
   fChain->SetBranchAddress("Muon_isMedium2016", &Muon_isMedium2016, &b_Muon_isMedium2016);
   fChain->SetBranchAddress("Muon_isSoft", &Muon_isSoft, &b_Muon_isSoft);
   fChain->SetBranchAddress("Muon_isTight", &Muon_isTight, &b_Muon_isTight);
   fChain->SetBranchAddress("Muon_pdgId", &Muon_pdgId, &b_Muon_pdgId);
   fChain->SetBranchAddress("Muon_rand1", &Muon_rand1, &b_Muon_rand1);
   fChain->SetBranchAddress("Muon_rand2", &Muon_rand2, &b_Muon_rand2);
   fChain->SetBranchAddress("Muon_relIso03", &Muon_relIso03, &b_Muon_relIso03);
   fChain->SetBranchAddress("Muon_relIso04", &Muon_relIso04, &b_Muon_relIso04);
   fChain->SetBranchAddress("Muon_trackerLayers", &Muon_trackerLayers, &b_Muon_trackerLayers);
   fChain->SetBranchAddress("PDFWeight_BjorkenX1", &PDFWeight_BjorkenX1, &b_PDFWeight_BjorkenX1);
   fChain->SetBranchAddress("PDFWeight_BjorkenX2", &PDFWeight_BjorkenX2, &b_PDFWeight_BjorkenX2);
   fChain->SetBranchAddress("PDFWeight_Cent", &PDFWeight_Cent, &b_PDFWeight_Cent);
   fChain->SetBranchAddress("PDFWeight_Cent_Down", &PDFWeight_Cent_Down, &b_PDFWeight_Cent_Down);
   fChain->SetBranchAddress("PDFWeight_Cent_Up", &PDFWeight_Cent_Up, &b_PDFWeight_Cent_Up);
   fChain->SetBranchAddress("PDFWeight_Id1", &PDFWeight_Id1, &b_PDFWeight_Id1);
   fChain->SetBranchAddress("PDFWeight_Id2", &PDFWeight_Id2, &b_PDFWeight_Id2);
   fChain->SetBranchAddress("PDFWeight_PDF1", &PDFWeight_PDF1, &b_PDFWeight_PDF1);
   fChain->SetBranchAddress("PDFWeight_PDF2", &PDFWeight_PDF2, &b_PDFWeight_PDF2);
   fChain->SetBranchAddress("PDFWeight_Q", &PDFWeight_Q, &b_PDFWeight_Q);
   fChain->SetBranchAddress("PDFWeight_Var1_Down", &PDFWeight_Var1_Down, &b_PDFWeight_Var1_Down);
   fChain->SetBranchAddress("PDFWeight_Var1_Up", &PDFWeight_Var1_Up, &b_PDFWeight_Var1_Up);
   fChain->SetBranchAddress("PDFWeight_Var2_Down", &PDFWeight_Var2_Down, &b_PDFWeight_Var2_Down);
   fChain->SetBranchAddress("PDFWeight_Var2_Up", &PDFWeight_Var2_Up, &b_PDFWeight_Var2_Up);
   fChain->SetBranchAddress("Photon", &Photon, &b_Photon);
   fChain->SetBranchAddress("Photon_Count", &Photon_Count, &b_Photon_Count);
   fChain->SetBranchAddress("Photon_SCB_Loose", &Photon_SCB_Loose, &b_Photon_SCB_Loose);
   fChain->SetBranchAddress("Photon_SCB_Medium", &Photon_SCB_Medium, &b_Photon_SCB_Medium);
   fChain->SetBranchAddress("Photon_SCB_Tight", &Photon_SCB_Tight, &b_Photon_SCB_Tight);
   fChain->SetBranchAddress("Photon_ChaHadIso", &Photon_ChaHadIso, &b_Photon_ChaHadIso);
   fChain->SetBranchAddress("Photon_NeuHadIso", &Photon_NeuHadIso, &b_Photon_NeuHadIso);
   fChain->SetBranchAddress("Photon_PhoIso", &Photon_PhoIso, &b_Photon_PhoIso);
   fChain->SetBranchAddress("Photon_WorstChargedIso", &Photon_WorstChargedIso, &b_Photon_WorstChargedIso);
   fChain->SetBranchAddress("Photon_ChaHadEffArea", &Photon_ChaHadEffArea, &b_Photon_ChaHadEffArea);
   fChain->SetBranchAddress("Photon_NeuHadEffArea", &Photon_NeuHadEffArea, &b_Photon_NeuHadEffArea);
   fChain->SetBranchAddress("Photon_PhoHadEffArea", &Photon_PhoHadEffArea, &b_Photon_PhoHadEffArea);
   fChain->SetBranchAddress("Photon_R9", &Photon_R9, &b_Photon_R9);
   fChain->SetBranchAddress("Photon_HoverE", &Photon_HoverE, &b_Photon_HoverE);
   fChain->SetBranchAddress("Photon_SuperCluster_Eta", &Photon_SuperCluster_Eta, &b_Photon_SuperCluster_Eta);
   fChain->SetBranchAddress("Photon_SuperCluster_Phi", &Photon_SuperCluster_Phi, &b_Photon_SuperCluster_Phi);
   fChain->SetBranchAddress("Photon_SuperCluster_EtaWidth", &Photon_SuperCluster_EtaWidth, &b_Photon_SuperCluster_EtaWidth);
   fChain->SetBranchAddress("Photon_SuperCluster_PhiWidth", &Photon_SuperCluster_PhiWidth, &b_Photon_SuperCluster_PhiWidth);
   fChain->SetBranchAddress("Photon_SuperCluster_Energy", &Photon_SuperCluster_Energy, &b_Photon_SuperCluster_Energy);
   fChain->SetBranchAddress("Photon_SuperCluster_RawEnergy", &Photon_SuperCluster_RawEnergy, &b_Photon_SuperCluster_RawEnergy);
   fChain->SetBranchAddress("Photon_ElectronVeto", &Photon_ElectronVeto, &b_Photon_ElectronVeto);
   fChain->SetBranchAddress("Photon_Full5x5_SigmaIetaIeta", &Photon_Full5x5_SigmaIetaIeta, &b_Photon_Full5x5_SigmaIetaIeta);
   fChain->SetBranchAddress("Photon_MVANonTrig_Tight", &Photon_MVANonTrig_Tight, &b_Photon_MVANonTrig_Tight);
   fChain->SetBranchAddress("PV_Count", &PV_Count, &b_PV_Count);
   fChain->SetBranchAddress("PileUp_Count_Interaction", &PileUp_Count_Interaction, &b_PileUp_Count_Interaction);
   fChain->SetBranchAddress("PileUp_Count_Intime", &PileUp_Count_Intime, &b_PileUp_Count_Intime);
   fChain->SetBranchAddress("Rho", &Rho, &b_Rho);
   fChain->SetBranchAddress("Trigger_Name", &Trigger_Name, &b_Trigger_Name);
   fChain->SetBranchAddress("Trigger_PreScale", &Trigger_PreScale, &b_Trigger_PreScale);
   fChain->SetBranchAddress("Trigger_isError", &Trigger_isError, &b_Trigger_isError);
   fChain->SetBranchAddress("Trigger_isPass", &Trigger_isPass, &b_Trigger_isPass);
   fChain->SetBranchAddress("Trigger_isRun", &Trigger_isRun, &b_Trigger_isRun);
   fChain->SetBranchAddress("Vertex_SumPtSquare", &Vertex_SumPtSquare, &b_Vertex_SumPtSquare);
   fChain->SetBranchAddress("Vertex_X", &Vertex_X, &b_Vertex_X);
   fChain->SetBranchAddress("Vertex_X_Error", &Vertex_X_Error, &b_Vertex_X_Error);
   fChain->SetBranchAddress("Vertex_Y", &Vertex_Y, &b_Vertex_Y);
   fChain->SetBranchAddress("Vertex_Y_Error", &Vertex_Y_Error, &b_Vertex_Y_Error);
   fChain->SetBranchAddress("Vertex_Z", &Vertex_Z, &b_Vertex_Z);
   fChain->SetBranchAddress("Vertex_Z_Error", &Vertex_Z_Error, &b_Vertex_Z_Error);
   Notify();
}

Bool_t SSBTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SSBTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SSBTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SSBTree_cxx
