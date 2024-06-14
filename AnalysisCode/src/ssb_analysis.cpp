#define ssb_analysis_cxx

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "TCanvas.h"
#include "TH2.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TStyle.h"

#include "./../CommonTools.hpp"
#include "./../interface/ssb_analysis.hpp"

using namespace std;

void ssb_analysis::Loop(char *logfile) {

  //////////
  if (fChain == 0)
    return;
  //////////

  //////////
  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  //////////

  /// My variables
  Long64_t __tot_evt = 0;
  double totalGenWeight = 0;

  /// Check Total Event

  ////////////////////////
  /// start event loop ///
  ////////////////////////

  double globalWeight = 1;
  double normSF = 1;

  // if (!fIsData) {
  //   normSF = GetNormalization();
  //   globalWeight *= normSF;
  // }

  // std::cout << fProcessName << " "
  //           << fIsData      << " "
  //           << normSF       << " "
  //           << globalWeight << " "
  //           << std::endl;

  std::cout << " DimuonMassCut : " << fMassCut << std::endl;

  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) {
      printf("ERROR: Could not load tree!!!\n");
      break;
    }

    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    __tot_evt++;
    // std::cout << jentry << " " << Gen_EventWeight << std::endl;

    if (jentry % 10000 == 0)
      printf("Event %lld\n", jentry); //%lld supports Long64_t

    std::vector<TLorentzVector> LHE_leptons;
    for (int i = 0; i < LHE_particleID->size(); i++) {
      // std::cout << jentry << " events  - " << i << " " << LHE_particleID->at(i) << " " << LHE_Status->at(i) << std::endl;

      if (std::fabs(LHE_particleID->at(i)) == 13) {
        TLorentzVector l;
        l.SetPxPyPzE(LHE_Px->at(i), LHE_Py->at(i), LHE_Pz->at(i), LHE_E->at(i));
        LHE_leptons.push_back(l);
      }
    }

    int LHE_leptons_size = LHE_leptons.size();
    double LHE_mass = 0;
    
    if ( LHE_leptons_size == 2 ) {
      auto LHE_lepton_leading = LHE_leptons.at(0);
      auto LHE_lepton_subleading = LHE_leptons.at(1);
      auto LHE_dilepton = LHE_lepton_leading + LHE_lepton_subleading;

      LHE_mass = LHE_dilepton.M();
    }

    if (fProcessName.find("NNLO_inc") != std::string::npos && LHE_mass > 100.)
      continue;

    // totalGenWeight += Gen_EventWeight;
    // globalWeight = Gen_EventWeight;

    if fProcessName.find("NNLO") != std::string::npos() {
      std::cout << "INCLUDING NNLO : gen weight set to +- 1" << std::endl;

      if (Gen_EventWeight < 0) globalWeight = -1;
      else                     globalWeight = +1;

    } else {

      globalWeight = Gen_EventWeight;

    }

    // if (LHE_mass <= 10.) 
    //   std::cout << "LHE dimuon mass: " << jentry << " " << globalWeight << " " << LHE_mass << std::endl;

    totalGenWeight += globalWeight;

    FillHisto(h_GenWeight, globalWeight, 1.);

    // globalWeight *= normSF;

    FillHisto(h_LHEnMuon, LHE_leptons_size, globalWeight);

    if (LHE_leptons.size() != 2) {
      FillHisto(h_LHEDimuonMass, 0., globalWeight);
    } else {
      FillHisto(h_LHEDimuonMass, LHE_mass, globalWeight);
    }

    double puReweightFactor = GetPUweight(PileUp_Count_Intime);
    if (!fIsData)
      globalWeight *= puReweightFactor;


    ////////////////////////////////////////
    /// start Main Loop Function and Cut ///
    ////////////////////////////////////////

    int index_TrigIsoMu24 = 0;
    int index_TrigIsoTkMu24 = 0;

    for (int i = 0; i < Trigger_Name->size(); i++) {
      if (Trigger_Name->at(i).find("IsoMu24") != std::string::npos)
        index_TrigIsoMu24 = i;

      if (Trigger_Name->at(i).find("IsoTkMu24") != std::string::npos)
        index_TrigIsoTkMu24 = i;
    }

    if (!(Trigger_isPass->at(index_TrigIsoMu24) ||
          Trigger_isPass->at(index_TrigIsoMu24)))
      continue;

    int index_HBHENoiseFilter = 0;
    int index_HBHENoiseIsoFilter = 0;
    int index_CSCTightHaloFilter = 0;
    int index_EcalDeadCellTriggerPrimitiveFilter = 0;
    int index_goodVertices = 0;
    int index_chargedHadronTrackResolutionFilter = 0;
    int index_muonBadTrackFilter = 0;
    int index_eeBadScFilter = 0;

    for (int idx = 0; idx < METFilter_Name->size(); idx++) {
      if (METFilter_Name->at(idx).find("Flag_HBHENoiseFilter") != string::npos)
        index_HBHENoiseFilter = idx;
      if (METFilter_Name->at(idx).find("Flag_HBHENoiseIsoFilter") !=
          string::npos)
        index_HBHENoiseIsoFilter = idx;
      if (METFilter_Name->at(idx).find("Flag_CSCTightHaloFilter") !=
          string::npos)
        index_CSCTightHaloFilter = idx;
      if (METFilter_Name->at(idx).find(
              "Flag_EcalDeadCellTriggerPrimitiveFilter") != string::npos)
        index_EcalDeadCellTriggerPrimitiveFilter = idx;
      if (METFilter_Name->at(idx).find("Flag_goodVertices") != string::npos)
        index_goodVertices = idx;
      if (METFilter_Name->at(idx).find(
              "Flag_chargedHadronTrackResolutionFilter") != string::npos)
        index_chargedHadronTrackResolutionFilter = idx;
      if (METFilter_Name->at(idx).find("Flag_muonBadTrackFilter") !=
          string::npos)
        index_muonBadTrackFilter = idx;
      if (METFilter_Name->at(idx).find("Flag_eeBadScFilter") != string::npos)
        index_eeBadScFilter = idx;
    }

    bool passed_filter =
        (METFilter_isPass->at(index_HBHENoiseFilter) &&
         METFilter_isPass->at(index_HBHENoiseIsoFilter) &&
         METFilter_isPass->at(index_CSCTightHaloFilter) &&
         METFilter_isPass->at(index_EcalDeadCellTriggerPrimitiveFilter) &&
         METFilter_isPass->at(index_goodVertices) &&
         METFilter_isPass->at(index_chargedHadronTrackResolutionFilter) &&
         METFilter_isPass->at(index_muonBadTrackFilter));

    if (fIsData)
      passed_filter =
          (passed_filter && METFilter_isPass->at(index_eeBadScFilter));

    if (!passed_filter)
      continue;

    // vector<bool>    *Muon_isPF;
    // vector<bool>    *Muon_isGlobal;
    // vector<bool>    *Muon_isTracker;

    // vector<double>  *Muon_ValidTrackHitFraction;
    // vector<double>  *Muon_GlobalTrackChi2;
    // vector<double>  *Muon_TrackerSTAPositionMatch;
    // vector<double>  *Muon_TrackKink;
    // vector<bool>    *Muon_TrackExist;
    // vector<double>  *Muon_SegmentCompatibility;

    // vector<double>  *Muon_GlobalMuonTrackChamberHit;
    // vector<double>  *Muon_TunePGlobalMuonTrackChamberHit;
    // vector<double>  *Muon_TunePpTError;
    // vector<double>  *Muon_TunePpT;
    // vector<double>  *Muon_StationsHasSegments;

    // vector<double>  *Muon_dxy;
    // vector<double>  *Muon_dB;
    // vector<double>  *Muon_dz;
    // vector<double>  *Muon_PixelHit;
    // vector<double>  *Muon_TrackLayerWithHit;


    std::vector<stdMUON> recoMuon;
    recoMuon.reserve(Muon->GetEntries());
    for (int i = 0; i < Muon->GetEntries(); ++i)
      recoMuon.emplace_back(stdMUON(*((TLorentzVector *)Muon->At(i)),
                                    Muon_Charge->at(i), Muon_isTight->at(i),
                                    Muon_PFIsodBeta04->at(i),
                                    Muon_isPF->at(i),
                                    Muon_isGlobal->at(i),
                                    Muon_isTracker->at(i),
                                    Muon_ValidTrackHitFraction->at(i),
                                    Muon_GlobalTrackChi2->at(i),
                                    Muon_TrackerSTAPositionMatch->at(i),
                                    Muon_TrackKink->at(i),
                                    Muon_TrackExist->at(i),
                                    Muon_SegmentCompatibility->at(i),
                                    Muon_GlobalMuonTrackChamberHit->at(i),
                                    Muon_TunePGlobalMuonTrackChamberHit->at(i),
                                    Muon_TunePpTError->at(i),
                                    Muon_TunePpT->at(i),
                                    Muon_StationsHasSegments->at(i),
                                    Muon_dxy->at(i),
                                    Muon_dB->at(i),
                                    Muon_dz->at(i),
                                    Muon_PixelHit->at(i),
                                    Muon_TrackLayerWithHit->at(i),
                                    Muon_trackerLayers->at(i)));

    for (int i = 0; i < recoMuon.size(); i++) {
      FillHisto(h_isGlobal_BeforeSel,                recoMuon.at(i).isGlobal, globalWeight);
      FillHisto(h_Chi2GlobalMuonTrack_BeforeSel,     recoMuon.at(i).GlobalTrackChi2, globalWeight);
      FillHisto(h_NChamberGlobalMuonTrack_BeforeSel, recoMuon.at(i).GlobalMuonTrackChamberHit, globalWeight);
      FillHisto(h_NofMatchedStation_BeforeSel,       recoMuon.at(i).StationsHasSegments, globalWeight);
      FillHisto(h_dxy_BeforeSel,                     recoMuon.at(i).dxy, globalWeight);
      FillHisto(h_dB_AfterSel,                       recoMuon.at(i).dB, globalWeight);
      FillHisto(h_dxyOrdB_AfterSel,                  (recoMuon.at(i).dxy < 0.2 || recoMuon.at(i).dB < 0.2), globalWeight);
      FillHisto(h_dz_BeforeSel,                      recoMuon.at(i).dz, globalWeight);
      FillHisto(h_NPixelHit_BeforeSel,               recoMuon.at(i).PixelHit, globalWeight);
      FillHisto(h_ValidInnerTrackHit_BeforeSel,      recoMuon.at(i).TrackLayerWithHit, globalWeight);
    }

    std::vector<TLorentzVector> genStdMuon;
    genStdMuon.reserve(GenMuon->GetEntries());
    for (int i = 0; i < GenMuon->GetEntries(); ++i)
      genStdMuon.emplace_back(*((TLorentzVector *)GenMuon->At(i)));

    if (fRoccor_enalbed) {
      for (int i = 0; i < recoMuon.size(); i++) {
        if (fIsData) {
          recoMuon.at(i).v *= fRoccoR->kScaleDT(
              recoMuon.at(i).charge, recoMuon.at(i).v.Pt(),
              recoMuon.at(i).v.Eta(), recoMuon.at(i).v.Phi(), 0, 0);
        } else {

          double drmin = 999.;
          bool match = false;
          int genMuonIdx = 0;

          for (int j = 0; j < genStdMuon.size(); j++) {
            if (recoMuon.at(i).v.DeltaR(genStdMuon.at(j)) < 0.1 &&
                recoMuon.at(i).v.DeltaR(genStdMuon.at(j)) < drmin) {
              match = true;
              genMuonIdx = j;
              drmin = recoMuon.at(i).v.DeltaR(genStdMuon.at(j));
            }
          }

          if (match) {
            recoMuon.at(i).v *= fRoccoR->kScaleFromGenMC(
                recoMuon.at(i).charge, recoMuon.at(i).v.Pt(),
                recoMuon.at(i).v.Eta(), recoMuon.at(i).v.Phi(),
                Muon_trackerLayers->at(i), genStdMuon.at(genMuonIdx).Pt(),
                Muon_rand2->at(i), 0, 0);
          } else {
            recoMuon.at(i).v *= fRoccoR->kScaleAndSmearMC(
                recoMuon.at(i).charge, recoMuon.at(i).v.Pt(),
                recoMuon.at(i).v.Eta(), recoMuon.at(i).v.Phi(),
                Muon_trackerLayers->at(i), Muon_rand1->at(i), Muon_rand2->at(i),
                0, 0);
          }
        }
      }
    }

    std::sort(recoMuon.begin(), recoMuon.end(),
              [](const stdMUON &lhs, const stdMUON &rhs) {
                return lhs.v.Pt() > rhs.v.Pt();
              });

    std::vector<stdMUON> setMuonsPassingCut;
    for (int i = 0; i < recoMuon.size(); i++)
      if (recoMuon.at(i).isPass())
        setMuonsPassingCut.push_back(recoMuon.at(i));

    // int pos = 0;
    // int neg = 0;

    // for (int i = 0; i < setMuonsPassingCut.size(); i++) {
    //   if (setMuonsPassingCut.at(i).charge > 0.)
    //     pos++;
    //   else
    //     neg++;
    // }

    // if (pos == 0 || neg == 0)
    //   continue;

    if (!(setMuonsPassingCut.at(0).v.Pt() > 25.))
      continue;

    int idx_subleading = -1;
    for (int i = 1; i < setMuonsPassingCut.size(); i++) {
      if (setMuonsPassingCut.at(0).charge * setMuonsPassingCut.at(i).charge < 0) {
        idx_subleading = i;
        break;
      }
    }

    if (idx_subleading == -1)
      continue;

    TLorentzVector leadingMuon = setMuonsPassingCut.at(0).v;
    TLorentzVector subleadingMuon = setMuonsPassingCut.at(idx_subleading).v;

    if (!((leadingMuon + subleadingMuon).M() > fMassCut))
      continue;

    if (!fIsData) {

      if (fIDISO_enalbed) {
        globalWeight *= fTracking_SF.getEfficiency(setMuonsPassingCut.at(0)._v.Pt(), setMuonsPassingCut.at(0)._v.Eta());
        globalWeight *= fTracking_SF.getEfficiency(setMuonsPassingCut.at(idx_subleading)._v.Pt(), setMuonsPassingCut.at(idx_subleading)._v.Eta());

        globalWeight *= fIDISO_SF.getEfficiency(setMuonsPassingCut.at(0)._v.Pt(), setMuonsPassingCut.at(0)._v.Eta());
        globalWeight *= fIDISO_SF.getEfficiency(setMuonsPassingCut.at(idx_subleading)._v.Pt(), setMuonsPassingCut.at(idx_subleading)._v.Eta());

        // std::cout << normSF << " "
        //           << puReweightFactor << " "
        //           << fIDISO_SF.getEfficiency(leadingMuon.Pt(),
        //           leadingMuon.Eta()) << " "
        //           << fIDISO_SF.getEfficiency(subleadingMuon.Pt(),
        //           subleadingMuon.Eta()) << " "
        //           << globalWeight << " "
        //           << std::endl;
      }

      if (fTRIGG_enalbed) {
        globalWeight *= fTRIG_SF.getTriggerEfficiency(setMuonsPassingCut.at(0)._v.Pt(), setMuonsPassingCut.at(0)._v.Eta(), setMuonsPassingCut.at(idx_subleading)._v.Pt(), setMuonsPassingCut.at(idx_subleading)._v.Eta());

        // std::cout << normSF << " "
        //           << puReweightFactor << " "
        //           << fIDISO_SF.getEfficiency(leadingMuon.Pt(),
        //           leadingMuon.Eta()) << " "
        //           << fIDISO_SF.getEfficiency(subleadingMuon.Pt(),
        //           subleadingMuon.Eta()) << " "
        //           << fTRIG_SF.getTriggerEfficiency(leadingMuon.Pt(),
        //           leadingMuon.Eta(), subleadingMuon.Pt(),
        //           subleadingMuon.Eta()) << " "
        //           << globalWeight << " "
        //           << std::endl;
      }
    }


    std::vector<int> muonIndex = {0, idx_subleading};
    for (int i = 0; i < 2; i++) {
      if ( !(setMuonsPassingCut.at(muonIndex.at(i)).dxy < 0.2 || setMuonsPassingCut.at(muonIndex.at(i)).dB < 0.2) || setMuonsPassingCut.at(muonIndex.at(i)).dz > 0.5 )
        std::cout << "MUON passing ID but suspecious : " 
                  << setMuonsPassingCut.at(muonIndex.at(i)).v.Pt() << " "
                  << setMuonsPassingCut.at(muonIndex.at(i)).v.Eta() << " "
                  << setMuonsPassingCut.at(muonIndex.at(i)).v.Phi() << " "
                  << (leadingMuon + subleadingMuon).M() << " "
                  << setMuonsPassingCut.at(muonIndex.at(i)).dxy << " "
                  << setMuonsPassingCut.at(muonIndex.at(i)).dB << " "
                  << setMuonsPassingCut.at(muonIndex.at(i)).dz << std::endl;
    }


    FillHisto(h_isGlobal_AfterSel,                setMuonsPassingCut.at(0).isGlobal, globalWeight);
    FillHisto(h_Chi2GlobalMuonTrack_AfterSel,     setMuonsPassingCut.at(0).GlobalTrackChi2, globalWeight);
    FillHisto(h_NChamberGlobalMuonTrack_AfterSel, setMuonsPassingCut.at(0).GlobalMuonTrackChamberHit, globalWeight);
    FillHisto(h_NofMatchedStation_AfterSel,       setMuonsPassingCut.at(0).StationsHasSegments, globalWeight);
    FillHisto(h_dxy_AfterSel,                     setMuonsPassingCut.at(0).dxy, globalWeight);
    FillHisto(h_dB_AfterSel,                      setMuonsPassingCut.at(0).dB, globalWeight);
    FillHisto(h_dxyOrdB_AfterSel,                 (setMuonsPassingCut.at(0).dxy < 0.2 || setMuonsPassingCut.at(0).dB < 0.2), globalWeight);
    FillHisto(h_dz_AfterSel,                      setMuonsPassingCut.at(0).dz, globalWeight);
    FillHisto(h_NPixelHit_AfterSel,               setMuonsPassingCut.at(0).PixelHit, globalWeight);
    FillHisto(h_ValidInnerTrackHit_AfterSel,      setMuonsPassingCut.at(0).TrackLayerWithHit, globalWeight);

    FillHisto(h_isGlobal_AfterSel,                setMuonsPassingCut.at(idx_subleading).isGlobal, globalWeight);
    FillHisto(h_Chi2GlobalMuonTrack_AfterSel,     setMuonsPassingCut.at(idx_subleading).GlobalTrackChi2, globalWeight);
    FillHisto(h_NChamberGlobalMuonTrack_AfterSel, setMuonsPassingCut.at(idx_subleading).GlobalMuonTrackChamberHit, globalWeight);
    FillHisto(h_NofMatchedStation_AfterSel,       setMuonsPassingCut.at(idx_subleading).StationsHasSegments, globalWeight);
    FillHisto(h_dxy_AfterSel,                     setMuonsPassingCut.at(idx_subleading).dxy, globalWeight);
    FillHisto(h_dB_AfterSel,                      setMuonsPassingCut.at(idx_subleading).dB, globalWeight);
    FillHisto(h_dxyOrdB_AfterSel,                 (setMuonsPassingCut.at(idx_subleading).dxy < 0.2 || setMuonsPassingCut.at(idx_subleading).dB < 0.2), globalWeight);    
    FillHisto(h_dz_AfterSel,                      setMuonsPassingCut.at(idx_subleading).dz, globalWeight);
    FillHisto(h_NPixelHit_AfterSel,               setMuonsPassingCut.at(idx_subleading).PixelHit, globalWeight);
    FillHisto(h_ValidInnerTrackHit_AfterSel,      setMuonsPassingCut.at(idx_subleading).TrackLayerWithHit, globalWeight);

    std::vector<stdELEC> recoElec;
    recoElec.reserve(Elec->GetEntries());
    for (int i = 0; i < Elec->GetEntries(); ++i)
      recoElec.emplace_back(stdELEC(
          *((TLorentzVector *)Elec->At(i)), Elec_Supercluster_Eta->at(i),
          Elec_SCB_Tight->at(i), Elec_PFIsoRho03->at(i)));

    std::vector<stdJET> recoJet;
    recoJet.reserve(Jet->GetEntries());
    for (int i = 0; i < Jet->GetEntries(); ++i)
      recoJet.emplace_back(stdJET(*((TLorentzVector *)Jet->At(i)),
                                  Jet_PFId->at(i), Jet_bDisc->at(i)));

    std::vector<TLorentzVector> idTightLep;

    for (int i = 0; i < recoMuon.size(); i++)
      if (recoMuon.at(i).isPass())
        idTightLep.push_back(recoMuon.at(i).v);

    for (int i = 0; i < recoElec.size(); i++)
      if (recoElec.at(i).isPass())
        idTightLep.push_back(recoElec.at(i).v);

    std::vector<TLorentzVector> selStdJet;
    std::vector<TLorentzVector> selStdBJet;
    for (int i = 0; i < recoJet.size(); i++) {

      if (!(recoJet.at(i).isPass() && !recoJet.at(i).isOverlap(idTightLep)))
        continue;

      selStdJet.push_back(recoJet.at(i).v);

      if (recoJet.at(i).isBJet())
        selStdBJet.push_back(recoJet.at(i).v);
    }

    int nJet_selected = selStdJet.size();
    int nBJet_selected = selStdBJet.size();

    for (int i = 0; i < nJet_selected; i++) {
      FillHisto(h_JetPt, selStdJet.at(i).Pt(), globalWeight);
      FillHisto(h_JetEta, selStdJet.at(i).Eta(), globalWeight);
      FillHisto(h_JetPhi, selStdJet.at(i).Phi(), globalWeight);
    }

    for (int i = 0; i < nBJet_selected; i++) {
      FillHisto(h_BJetPt, selStdBJet.at(i).Pt(), globalWeight);
      FillHisto(h_BJetEta, selStdBJet.at(i).Eta(), globalWeight);
      FillHisto(h_BJetPhi, selStdBJet.at(i).Phi(), globalWeight);
    }

    FillHisto(h_nJet_before, nJet_selected, globalWeight);
    FillHisto(h_nJet_after, nJet_selected, globalWeight);

    FillHisto(h_nBJet_before, nBJet_selected, globalWeight);
    FillHisto(h_nBJet_after, nBJet_selected, globalWeight);

    TLorentzVector diMuon = leadingMuon + subleadingMuon;

    FillHisto(h_LeadingMuonPt, leadingMuon.Pt(), globalWeight);
    FillHisto(h_LeadingMuonEta, leadingMuon.Eta(), globalWeight);
    FillHisto(h_LeadingMuonPhi, leadingMuon.Phi(), globalWeight);

    FillHisto(h_MuonPt, leadingMuon.Pt(), globalWeight);
    FillHisto(h_MuonEta, leadingMuon.Eta(), globalWeight);
    FillHisto(h_MuonPhi, leadingMuon.Phi(), globalWeight);

    FillHisto(h_SubleadingMuonPt, subleadingMuon.Pt(), globalWeight);
    FillHisto(h_SubleadingMuonEta, subleadingMuon.Eta(), globalWeight);
    FillHisto(h_SubleadingMuonPhi, subleadingMuon.Phi(), globalWeight);

    FillHisto(h_MuonPt, subleadingMuon.Pt(), globalWeight);
    FillHisto(h_MuonEta, subleadingMuon.Eta(), globalWeight);
    FillHisto(h_MuonPhi, subleadingMuon.Phi(), globalWeight);

    FillHisto(h_dimuonMass, diMuon.M(), globalWeight);
    FillHisto(h_dimuonMass_wide, diMuon.M(), globalWeight);
    FillHisto(h_dimuonPt, diMuon.Pt(), globalWeight);
    FillHisto(h_dimuonRap, diMuon.Rapidity(), globalWeight);

    FillHisto(h_PV_Count_before_corr, PV_Count, globalWeight / puReweightFactor);
    FillHisto(h_PileUp_Count_Interaction_before_corr, PileUp_Count_Interaction, globalWeight / puReweightFactor);
    FillHisto(h_PileUp_Count_Intime_before_corr, PileUp_Count_Intime, globalWeight / puReweightFactor);

    FillHisto(h_PV_Count_after_corr, PV_Count, globalWeight);
    FillHisto(h_PileUp_Count_Interaction_after_corr, PileUp_Count_Interaction, globalWeight);
    FillHisto(h_PileUp_Count_Intime_after_corr, PileUp_Count_Intime, globalWeight);

    if (nJet_selected == 0 && nBJet_selected == 0) {
      FillHisto(h_LeadingMuonPt_0J0B, leadingMuon.Pt(), globalWeight);
      FillHisto(h_LeadingMuonEta_0J0B, leadingMuon.Eta(), globalWeight);
      FillHisto(h_LeadingMuonPhi_0J0B, leadingMuon.Phi(), globalWeight);
      FillHisto(h_MuonPt_0J0B, leadingMuon.Pt(), globalWeight);
      FillHisto(h_MuonEta_0J0B, leadingMuon.Eta(), globalWeight);
      FillHisto(h_MuonPhi_0J0B, leadingMuon.Phi(), globalWeight);
      FillHisto(h_SubleadingMuonPt_0J0B, subleadingMuon.Pt(), globalWeight);
      FillHisto(h_SubleadingMuonEta_0J0B, subleadingMuon.Eta(), globalWeight);
      FillHisto(h_SubleadingMuonPhi_0J0B, subleadingMuon.Phi(), globalWeight);
      FillHisto(h_MuonPt_0J0B, subleadingMuon.Pt(), globalWeight);
      FillHisto(h_MuonEta_0J0B, subleadingMuon.Eta(), globalWeight);
      FillHisto(h_MuonPhi_0J0B, subleadingMuon.Phi(), globalWeight);
      FillHisto(h_dimuonMass_0J0B, diMuon.M(), globalWeight);
      FillHisto(h_dimuonMass_wide_0J0B, diMuon.M(), globalWeight);
      FillHisto(h_dimuonPt_0J0B, diMuon.Pt(), globalWeight);
      FillHisto(h_dimuonRap_0J0B, diMuon.Rapidity(), globalWeight);

      for (int i = 0; i < nJet_selected; i++) {
        FillHisto(h_JetPt_0J0B, selStdJet.at(i).Pt(), globalWeight);
        FillHisto(h_JetEta_0J0B, selStdJet.at(i).Eta(), globalWeight);
        FillHisto(h_JetPhi_0J0B, selStdJet.at(i).Phi(), globalWeight);
      }

      for (int i = 0; i < nBJet_selected; i++) {
        FillHisto(h_BJetPt_0J0B, selStdBJet.at(i).Pt(), globalWeight);
        FillHisto(h_BJetEta_0J0B, selStdBJet.at(i).Eta(), globalWeight);
        FillHisto(h_BJetPhi_0J0B, selStdBJet.at(i).Phi(), globalWeight);
      }
    }

    if (nJet_selected == 0) {
      FillHisto(h_LeadingMuonPt_0J, leadingMuon.Pt(), globalWeight);
      FillHisto(h_LeadingMuonEta_0J, leadingMuon.Eta(), globalWeight);
      FillHisto(h_LeadingMuonPhi_0J, leadingMuon.Phi(), globalWeight);
      FillHisto(h_MuonPt_0J, leadingMuon.Pt(), globalWeight);
      FillHisto(h_MuonEta_0J, leadingMuon.Eta(), globalWeight);
      FillHisto(h_MuonPhi_0J, leadingMuon.Phi(), globalWeight);
      FillHisto(h_SubleadingMuonPt_0J, subleadingMuon.Pt(), globalWeight);
      FillHisto(h_SubleadingMuonEta_0J, subleadingMuon.Eta(), globalWeight);
      FillHisto(h_SubleadingMuonPhi_0J, subleadingMuon.Phi(), globalWeight);
      FillHisto(h_MuonPt_0J, subleadingMuon.Pt(), globalWeight);
      FillHisto(h_MuonEta_0J, subleadingMuon.Eta(), globalWeight);
      FillHisto(h_MuonPhi_0J, subleadingMuon.Phi(), globalWeight);
      FillHisto(h_dimuonMass_0J, diMuon.M(), globalWeight);
      FillHisto(h_dimuonMass_wide_0J, diMuon.M(), globalWeight);
      FillHisto(h_dimuonPt_0J, diMuon.Pt(), globalWeight);
      FillHisto(h_dimuonRap_0J, diMuon.Rapidity(), globalWeight);

      for (int i = 0; i < nJet_selected; i++) {
        FillHisto(h_JetPt_0J, selStdJet.at(i).Pt(), globalWeight);
        FillHisto(h_JetEta_0J, selStdJet.at(i).Eta(), globalWeight);
        FillHisto(h_JetPhi_0J, selStdJet.at(i).Phi(), globalWeight);
      }

      for (int i = 0; i < nBJet_selected; i++) {
        FillHisto(h_BJetPt_0J, selStdBJet.at(i).Pt(), globalWeight);
        FillHisto(h_BJetEta_0J, selStdBJet.at(i).Eta(), globalWeight);
        FillHisto(h_BJetPhi_0J, selStdBJet.at(i).Phi(), globalWeight);
      }
    }

    if (nJet_selected == 1) {
      FillHisto(h_LeadingMuonPt_1J, leadingMuon.Pt(), globalWeight);
      FillHisto(h_LeadingMuonEta_1J, leadingMuon.Eta(), globalWeight);
      FillHisto(h_LeadingMuonPhi_1J, leadingMuon.Phi(), globalWeight);
      FillHisto(h_MuonPt_1J, leadingMuon.Pt(), globalWeight);
      FillHisto(h_MuonEta_1J, leadingMuon.Eta(), globalWeight);
      FillHisto(h_MuonPhi_1J, leadingMuon.Phi(), globalWeight);
      FillHisto(h_SubleadingMuonPt_1J, subleadingMuon.Pt(), globalWeight);
      FillHisto(h_SubleadingMuonEta_1J, subleadingMuon.Eta(), globalWeight);
      FillHisto(h_SubleadingMuonPhi_1J, subleadingMuon.Phi(), globalWeight);
      FillHisto(h_MuonPt_1J, subleadingMuon.Pt(), globalWeight);
      FillHisto(h_MuonEta_1J, subleadingMuon.Eta(), globalWeight);
      FillHisto(h_MuonPhi_1J, subleadingMuon.Phi(), globalWeight);
      FillHisto(h_dimuonMass_1J, diMuon.M(), globalWeight);
      FillHisto(h_dimuonMass_wide_1J, diMuon.M(), globalWeight);
      FillHisto(h_dimuonPt_1J, diMuon.Pt(), globalWeight);
      FillHisto(h_dimuonRap_1J, diMuon.Rapidity(), globalWeight);

      for (int i = 0; i < nJet_selected; i++) {
        FillHisto(h_JetPt_1J, selStdJet.at(i).Pt(), globalWeight);
        FillHisto(h_JetEta_1J, selStdJet.at(i).Eta(), globalWeight);
        FillHisto(h_JetPhi_1J, selStdJet.at(i).Phi(), globalWeight);
      }

      for (int i = 0; i < nBJet_selected; i++) {
        FillHisto(h_BJetPt_1J, selStdBJet.at(i).Pt(), globalWeight);
        FillHisto(h_BJetEta_1J, selStdBJet.at(i).Eta(), globalWeight);
        FillHisto(h_BJetPhi_1J, selStdBJet.at(i).Phi(), globalWeight);
      }
    }

    if (nJet_selected > 1) {
      FillHisto(h_LeadingMuonPt_mt1J, leadingMuon.Pt(), globalWeight);
      FillHisto(h_LeadingMuonEta_mt1J, leadingMuon.Eta(), globalWeight);
      FillHisto(h_LeadingMuonPhi_mt1J, leadingMuon.Phi(), globalWeight);
      FillHisto(h_MuonPt_mt1J, leadingMuon.Pt(), globalWeight);
      FillHisto(h_MuonEta_mt1J, leadingMuon.Eta(), globalWeight);
      FillHisto(h_MuonPhi_mt1J, leadingMuon.Phi(), globalWeight);
      FillHisto(h_SubleadingMuonPt_mt1J, subleadingMuon.Pt(), globalWeight);
      FillHisto(h_SubleadingMuonEta_mt1J, subleadingMuon.Eta(), globalWeight);
      FillHisto(h_SubleadingMuonPhi_mt1J, subleadingMuon.Phi(), globalWeight);
      FillHisto(h_MuonPt_mt1J, subleadingMuon.Pt(), globalWeight);
      FillHisto(h_MuonEta_mt1J, subleadingMuon.Eta(), globalWeight);
      FillHisto(h_MuonPhi_mt1J, subleadingMuon.Phi(), globalWeight);
      FillHisto(h_dimuonMass_mt1J, diMuon.M(), globalWeight);
      FillHisto(h_dimuonMass_wide_mt1J, diMuon.M(), globalWeight);
      FillHisto(h_dimuonPt_mt1J, diMuon.Pt(), globalWeight);
      FillHisto(h_dimuonRap_mt1J, diMuon.Rapidity(), globalWeight);

      for (int i = 0; i < nJet_selected; i++) {
        FillHisto(h_JetPt_mt1J, selStdJet.at(i).Pt(), globalWeight);
        FillHisto(h_JetEta_mt1J, selStdJet.at(i).Eta(), globalWeight);
        FillHisto(h_JetPhi_mt1J, selStdJet.at(i).Phi(), globalWeight);
      }

      for (int i = 0; i < nBJet_selected; i++) {
        FillHisto(h_BJetPt_mt1J, selStdBJet.at(i).Pt(), globalWeight);
        FillHisto(h_BJetEta_mt1J, selStdBJet.at(i).Eta(), globalWeight);
        FillHisto(h_BJetPhi_mt1J, selStdBJet.at(i).Phi(), globalWeight);
      }
    }

    if (nBJet_selected == 0) {
      FillHisto(h_LeadingMuonPt_0BJ, leadingMuon.Pt(), globalWeight);
      FillHisto(h_LeadingMuonEta_0BJ, leadingMuon.Eta(), globalWeight);
      FillHisto(h_LeadingMuonPhi_0BJ, leadingMuon.Phi(), globalWeight);
      FillHisto(h_MuonPt_0BJ, leadingMuon.Pt(), globalWeight);
      FillHisto(h_MuonEta_0BJ, leadingMuon.Eta(), globalWeight);
      FillHisto(h_MuonPhi_0BJ, leadingMuon.Phi(), globalWeight);
      FillHisto(h_SubleadingMuonPt_0BJ, subleadingMuon.Pt(), globalWeight);
      FillHisto(h_SubleadingMuonEta_0BJ, subleadingMuon.Eta(), globalWeight);
      FillHisto(h_SubleadingMuonPhi_0BJ, subleadingMuon.Phi(), globalWeight);
      FillHisto(h_MuonPt_0BJ, subleadingMuon.Pt(), globalWeight);
      FillHisto(h_MuonEta_0BJ, subleadingMuon.Eta(), globalWeight);
      FillHisto(h_MuonPhi_0BJ, subleadingMuon.Phi(), globalWeight);
      FillHisto(h_dimuonMass_0BJ, diMuon.M(), globalWeight);
      FillHisto(h_dimuonMass_wide_0BJ, diMuon.M(), globalWeight);
      FillHisto(h_dimuonPt_0BJ, diMuon.Pt(), globalWeight);
      FillHisto(h_dimuonRap_0BJ, diMuon.Rapidity(), globalWeight);

      for (int i = 0; i < nJet_selected; i++) {
        FillHisto(h_JetPt_0BJ, selStdJet.at(i).Pt(), globalWeight);
        FillHisto(h_JetEta_0BJ, selStdJet.at(i).Eta(), globalWeight);
        FillHisto(h_JetPhi_0BJ, selStdJet.at(i).Phi(), globalWeight);
      }

      for (int i = 0; i < nBJet_selected; i++) {
        FillHisto(h_BJetPt_0BJ, selStdBJet.at(i).Pt(), globalWeight);
        FillHisto(h_BJetEta_0BJ, selStdBJet.at(i).Eta(), globalWeight);
        FillHisto(h_BJetPhi_0BJ, selStdBJet.at(i).Phi(), globalWeight);
      }
    }

    if (nBJet_selected == 1) {
      FillHisto(h_LeadingMuonPt_1BJ, leadingMuon.Pt(), globalWeight);
      FillHisto(h_LeadingMuonEta_1BJ, leadingMuon.Eta(), globalWeight);
      FillHisto(h_LeadingMuonPhi_1BJ, leadingMuon.Phi(), globalWeight);
      FillHisto(h_MuonPt_1BJ, leadingMuon.Pt(), globalWeight);
      FillHisto(h_MuonEta_1BJ, leadingMuon.Eta(), globalWeight);
      FillHisto(h_MuonPhi_1BJ, leadingMuon.Phi(), globalWeight);
      FillHisto(h_SubleadingMuonPt_1BJ, subleadingMuon.Pt(), globalWeight);
      FillHisto(h_SubleadingMuonEta_1BJ, subleadingMuon.Eta(), globalWeight);
      FillHisto(h_SubleadingMuonPhi_1BJ, subleadingMuon.Phi(), globalWeight);
      FillHisto(h_MuonPt_1BJ, subleadingMuon.Pt(), globalWeight);
      FillHisto(h_MuonEta_1BJ, subleadingMuon.Eta(), globalWeight);
      FillHisto(h_MuonPhi_1BJ, subleadingMuon.Phi(), globalWeight);
      FillHisto(h_dimuonMass_1BJ, diMuon.M(), globalWeight);
      FillHisto(h_dimuonMass_wide_1BJ, diMuon.M(), globalWeight);
      FillHisto(h_dimuonPt_1BJ, diMuon.Pt(), globalWeight);
      FillHisto(h_dimuonRap_1BJ, diMuon.Rapidity(), globalWeight);

      for (int i = 0; i < nJet_selected; i++) {
        FillHisto(h_JetPt_1BJ, selStdJet.at(i).Pt(), globalWeight);
        FillHisto(h_JetEta_1BJ, selStdJet.at(i).Eta(), globalWeight);
        FillHisto(h_JetPhi_1BJ, selStdJet.at(i).Phi(), globalWeight);
      }

      for (int i = 0; i < nBJet_selected; i++) {
        FillHisto(h_BJetPt_1BJ, selStdBJet.at(i).Pt(), globalWeight);
        FillHisto(h_BJetEta_1BJ, selStdBJet.at(i).Eta(), globalWeight);
        FillHisto(h_BJetPhi_1BJ, selStdBJet.at(i).Phi(), globalWeight);
      }
    }

    if (nBJet_selected > 1) {
      FillHisto(h_LeadingMuonPt_mt1BJ, leadingMuon.Pt(), globalWeight);
      FillHisto(h_LeadingMuonEta_mt1BJ, leadingMuon.Eta(), globalWeight);
      FillHisto(h_LeadingMuonPhi_mt1BJ, leadingMuon.Phi(), globalWeight);
      FillHisto(h_MuonPt_mt1BJ, leadingMuon.Pt(), globalWeight);
      FillHisto(h_MuonEta_mt1BJ, leadingMuon.Eta(), globalWeight);
      FillHisto(h_MuonPhi_mt1BJ, leadingMuon.Phi(), globalWeight);
      FillHisto(h_SubleadingMuonPt_mt1BJ, subleadingMuon.Pt(), globalWeight);
      FillHisto(h_SubleadingMuonEta_mt1BJ, subleadingMuon.Eta(), globalWeight);
      FillHisto(h_SubleadingMuonPhi_mt1BJ, subleadingMuon.Phi(), globalWeight);
      FillHisto(h_MuonPt_mt1BJ, subleadingMuon.Pt(), globalWeight);
      FillHisto(h_MuonEta_mt1BJ, subleadingMuon.Eta(), globalWeight);
      FillHisto(h_MuonPhi_mt1BJ, subleadingMuon.Phi(), globalWeight);
      FillHisto(h_dimuonMass_mt1BJ, diMuon.M(), globalWeight);
      FillHisto(h_dimuonMass_wide_mt1BJ, diMuon.M(), globalWeight);
      FillHisto(h_dimuonPt_mt1BJ, diMuon.Pt(), globalWeight);
      FillHisto(h_dimuonRap_mt1BJ, diMuon.Rapidity(), globalWeight);

      for (int i = 0; i < nJet_selected; i++) {
        FillHisto(h_JetPt_mt1BJ, selStdJet.at(i).Pt(), globalWeight);
        FillHisto(h_JetEta_mt1BJ, selStdJet.at(i).Eta(), globalWeight);
        FillHisto(h_JetPhi_mt1BJ, selStdJet.at(i).Phi(), globalWeight);
      }

      for (int i = 0; i < nBJet_selected; i++) {
        FillHisto(h_BJetPt_mt1BJ, selStdBJet.at(i).Pt(), globalWeight);
        FillHisto(h_BJetEta_mt1BJ, selStdBJet.at(i).Eta(), globalWeight);
        FillHisto(h_BJetPhi_mt1BJ, selStdBJet.at(i).Phi(), globalWeight);
      }
    }



  } // event loop

  printf("Total processed number of events: %lld\n", __tot_evt);
  std::cout << "Total Gen Weigh : " << totalGenWeight << " out of " << __tot_evt
            << std::endl;

} // end Loop function

void ssb_analysis::Start(int genLoopon) {
  if (genLoopon == 0) {
    fout = new TFile(Form("output/%s", outfile), "RECREATE");
  } else if (genLoopon == 1) {
    fout = new TFile(Form("output/%s", outfile), "UPDATE");
  } else {
    cout << "genLoopon error" << endl;
  }
  fout->cd("");

  std::cout << "is Data?     : " << GetIsData() << std::endl;
  std::cout << "Process Name : " << GetProcessName() << std::endl;

  TDirectory *dir = gDirectory;
  dir->cd();

  DeclareHistos();
}

double ssb_analysis::GetLHEmass() {

  // for (int i = 0; i < LHE_nParticle; i++) {
  //   std::cout << "  - " << i << " " << LHE_particleID->at(i) << " " << LHE_Status->at(i) << std::endl;
  // }

  // std::cout << " " << std::endl;
  // std::cout << " " << std::endl;

  return 1;
}

void ssb_analysis::DeclareHistos() {

  // Test For Systematic All-in-One Code

  h_LHEDimuonMass = new TH1D(Form("h_LHEDimuonMass"), Form("h_LHEDimuonMass"), 6000, 0., 6000.);
  h_GenWeight = new TH1D(Form("h_GenWeight"), Form("h_GenWeight"), 20000, -10000., 10000.);
  h_LHEnMuon = new TH1D(Form("h_LHEnMuon"), Form("h_LHEnMuon"), 10, 0., 10.);

  h_nJet_before =
      new TH1D(Form("h_nJet_before"), Form("nJet_before"), 20, 0., 20.);
  h_nJet_before->Sumw2();
  h_nJet_after =
      new TH1D(Form("h_nJet_after"), Form("nJet_after"), 20, 0., 20.);
  h_nJet_after->Sumw2();

  h_nBJet_before =
      new TH1D(Form("h_nBJet_before"), Form("nJet_before"), 20, 0., 20.);
  h_nBJet_before->Sumw2();
  h_nBJet_after =
      new TH1D(Form("h_nBJet_after"), Form("nJet_after"), 20, 0., 20.);
  h_nBJet_after->Sumw2();

  h_JetPt = new TH1D(Form("h_JetPt"), Form("Jet_pT"), 1000, 0, 1000);
  h_JetPt->Sumw2();
  h_JetEta = new TH1D(Form("h_JetEta"), Form("Jet_Eta"), 60, -3., 3.);
  h_JetEta->Sumw2();
  h_JetPhi =
      new TH1D(Form("h_JetPhi"), Form("Jet_Phi"), 60, -3.141594, 3.141594);
  h_JetPhi->Sumw2();
  h_BJetPt = new TH1D(Form("h_BJetPt"), Form("BJet_pT"), 1000, 0, 1000);
  h_BJetPt->Sumw2();
  h_BJetEta = new TH1D(Form("h_BJetEta"), Form("BJet_Eta"), 60, -3., 3.);
  h_BJetEta->Sumw2();
  h_BJetPhi =
      new TH1D(Form("h_BJetPhi"), Form("BJet_Phi"), 60, -3.141594, 3.141594);
  h_BJetPhi->Sumw2();

  h_LeadingMuonPt =
      new TH1D(Form("h_LeadingMuonPt"), Form("Muon_pT"), 1000, 0, 1000);
  h_LeadingMuonPt->Sumw2();
  h_LeadingMuonEta =
      new TH1D(Form("h_LeadingMuonEta"), Form("Muon_Eta"), 60, -3., 3.);
  h_LeadingMuonEta->Sumw2();
  h_LeadingMuonPhi = new TH1D(Form("h_LeadingMuonPhi"), Form("Muon_Phi"), 60,
                              -3.141594, 3.141594);
  h_LeadingMuonPhi->Sumw2();

  h_SubleadingMuonPt =
      new TH1D(Form("h_SubleadingMuonPt"), Form("Muon_pT"), 1000, 0, 1000);
  h_SubleadingMuonPt->Sumw2();
  h_SubleadingMuonEta =
      new TH1D(Form("h_SubleadingMuonEta"), Form("Muon_Eta"), 60, -3., 3.);
  h_SubleadingMuonEta->Sumw2();
  h_SubleadingMuonPhi = new TH1D(Form("h_SubleadingMuonPhi"), Form("Muon_Phi"),
                                 60, -3.141594, 3.141594);
  h_SubleadingMuonPhi->Sumw2();

  h_MuonPt = new TH1D(Form("h_MuonPt"), Form("Muon_pT"), 1000, 0, 1000);
  h_MuonPt->Sumw2();
  h_MuonEta = new TH1D(Form("h_MuonEta"), Form("Muon_Eta"), 60, -3., 3.);
  h_MuonEta->Sumw2();
  h_MuonPhi =
      new TH1D(Form("h_MuonPhi"), Form("Muon_Phi"), 60, -3.141594, 3.141594);
  h_MuonPhi->Sumw2();

  std::vector<float> xbins = {
      0,   10,  15,  20,  25,  30,  35,  40,   45,   50,   55,  60,
      64,  68,  72,  76,  81,  86,  91,  96,   101,  106,  110, 115,
      120, 126, 133, 141, 150, 160, 171, 185,  200,  220,  243, 273,
      320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000, 3010};

  h_dimuonMass =
      new TH1D(Form("h_dimuonMass"), Form("inv_Mass"), 100, 41., 141.);
  h_dimuonMass->Sumw2();
  h_dimuonMass_wide = new TH1D(Form("h_dimuonMass_wide"), Form("inv_Mass"),
                               xbins.size() - 1, &(xbins[0]));
  h_dimuonMass_wide->Sumw2();

  h_dimuonPt = new TH1D(Form("h_dimuonPt"), Form("dimuon_pT"), 1000, 0., 1000.);
  h_dimuonPt->Sumw2();
  h_dimuonRap = new TH1D(Form("h_dimuonRap"), Form("dimuon_rap"), 60, -3., 3.);
  h_dimuonRap->Sumw2();

  h_PV_Count_before_corr = new TH1D(Form("h_PV_Count_before_PUcorr"),
                                    Form("PV_Count"), 100, 0., 100.);
  h_PV_Count_before_corr->Sumw2();
  h_PileUp_Count_Interaction_before_corr =
      new TH1D(Form("h_PileUp_Count_Interaction_before_PUcorr"),
               Form("PileUp_Count_Interaction"), 1000, 0., 100.);
  h_PileUp_Count_Interaction_before_corr->Sumw2();
  h_PileUp_Count_Intime_before_corr =
      new TH1D(Form("h_PileUp_Count_Intime_before_PUcorr"),
               Form("PileUp_Count_Intime"), 1000, 0., 100.);
  h_PileUp_Count_Intime_before_corr->Sumw2();

  h_PV_Count_after_corr = new TH1D(Form("h_PV_Count_after_PUcorr"),
                                   Form("PV_Count"), 100, 0., 100.);
  h_PV_Count_after_corr->Sumw2();
  h_PileUp_Count_Interaction_after_corr =
      new TH1D(Form("h_PileUp_Count_Interaction_after_PUcorr"),
               Form("PileUp_Count_Interaction"), 1000, 0., 100.);
  h_PileUp_Count_Interaction_after_corr->Sumw2();
  h_PileUp_Count_Intime_after_corr =
      new TH1D(Form("h_PileUp_Count_Intime_after_PUcorr"),
               Form("PileUp_Count_Intime"), 1000, 0., 100.);
  h_PileUp_Count_Intime_after_corr->Sumw2();


  h_isGlobal_BeforeSel = new TH1D(Form("h_isGlobal_BeforeSel"), Form("h_isGlobal_BeforeSel"), 2, 0., 2.); h_isGlobal_BeforeSel->Sumw2();
  h_Chi2GlobalMuonTrack_BeforeSel = new TH1D(Form("h_Chi2GlobalMuonTrack_BeforeSel"), Form("h_Chi2GlobalMuonTrack_BeforeSel"), 200, 0., 100.); h_Chi2GlobalMuonTrack_BeforeSel->Sumw2();
  h_NChamberGlobalMuonTrack_BeforeSel = new TH1D(Form("h_NChamberGlobalMuonTrack_BeforeSel"), Form("h_NChamberGlobalMuonTrack_BeforeSel"), 100, 0., 100.); h_NChamberGlobalMuonTrack_BeforeSel->Sumw2();
  h_NofMatchedStation_BeforeSel = new TH1D(Form("h_NofMatchedStation_BeforeSel"), Form("h_NofMatchedStation_BeforeSel"), 100, 0., 100.); h_NofMatchedStation_BeforeSel->Sumw2();
  h_dxy_BeforeSel = new TH1D(Form("h_dxy_BeforeSel"), Form("h_dxy_BeforeSel"), 210, 0., 21.); h_dxy_BeforeSel->Sumw2();
  h_dB_BeforeSel = new TH1D(Form("h_dB_BeforeSel"), Form("h_dB_BeforeSel"), 210, 0., 21.); h_dB_BeforeSel->Sumw2();
  h_dxyOrdB_BeforeSel = new TH1D(Form("h_dxyOrdB_BeforeSel"), Form("h_dxyOrdB_BeforeSel"), 2, 0., 2.); h_dxyOrdB_BeforeSel->Sumw2();
  h_dz_BeforeSel = new TH1D(Form("h_dz_BeforeSel"), Form("h_dz_BeforeSel"), 210, 0., 21.); h_dz_BeforeSel->Sumw2();
  h_NPixelHit_BeforeSel = new TH1D(Form("h_NPixelHit_BeforeSel"), Form("h_NPixelHit_BeforeSel"), 100, 0., 100.); h_NPixelHit_BeforeSel->Sumw2();
  h_ValidInnerTrackHit_BeforeSel = new TH1D(Form("h_ValidInnerTrackHit_BeforeSel"), Form("h_ValidInnerTrackHit_BeforeSel"), 100, 0., 100.); h_ValidInnerTrackHit_BeforeSel->Sumw2();

  h_isGlobal_AfterSel = new TH1D(Form("h_isGlobal_AfterSel"), Form("h_isGlobal_AfterSel"), 2, 0., 2.); h_isGlobal_AfterSel->Sumw2();
  h_Chi2GlobalMuonTrack_AfterSel = new TH1D(Form("h_Chi2GlobalMuonTrack_AfterSel"), Form("h_Chi2GlobalMuonTrack_AfterSel"), 200, 0., 100.); h_Chi2GlobalMuonTrack_AfterSel->Sumw2();
  h_NChamberGlobalMuonTrack_AfterSel = new TH1D(Form("h_NChamberGlobalMuonTrack_AfterSel"), Form("h_NChamberGlobalMuonTrack_AfterSel"), 100, 0., 100.); h_NChamberGlobalMuonTrack_AfterSel->Sumw2();
  h_NofMatchedStation_AfterSel = new TH1D(Form("h_NofMatchedStation_AfterSel"), Form("h_NofMatchedStation_AfterSel"), 100, 0., 100.); h_NofMatchedStation_AfterSel->Sumw2();
  h_dxy_AfterSel = new TH1D(Form("h_dxy_AfterSel"), Form("h_dxy_AfterSel"), 210, 0., 21.); h_dxy_AfterSel->Sumw2();
  h_dB_AfterSel = new TH1D(Form("h_dB_AfterSel"), Form("h_dB_AfterSel"), 210, 0., 21.); h_dB_AfterSel->Sumw2();
  h_dxyOrdB_AfterSel = new TH1D(Form("h_dxyOrdB_AfterSel"), Form("h_dxyOrdB_AfterSel"), 2, 0., 2.); h_dxyOrdB_AfterSel->Sumw2();  
  h_dz_AfterSel = new TH1D(Form("h_dz_AfterSel"), Form("h_dz_AfterSel"), 210, 0., 21.); h_dz_AfterSel->Sumw2();
  h_NPixelHit_AfterSel = new TH1D(Form("h_NPixelHit_AfterSel"), Form("h_NPixelHit_AfterSel"), 100, 0., 100.); h_NPixelHit_AfterSel->Sumw2();
  h_ValidInnerTrackHit_AfterSel = new TH1D(Form("h_ValidInnerTrackHit_AfterSel"), Form("h_ValidInnerTrackHit_AfterSel"), 100, 0., 100.); h_ValidInnerTrackHit_AfterSel->Sumw2();




  h_JetPt_0J0B = new TH1D(Form("h_JetPt_0J0B"), Form("Jet_pT"), 1000, 0, 1000);
  h_JetPt_0J0B->Sumw2();
  h_JetEta_0J0B = new TH1D(Form("h_JetEta_0J0B"), Form("Jet_Eta"), 60, -3., 3.);
  h_JetEta_0J0B->Sumw2();
  h_JetPhi_0J0B =
      new TH1D(Form("h_JetPhi_0J0B"), Form("Jet_Phi"), 60, -3.141594, 3.141594);
  h_JetPhi_0J0B->Sumw2();
  h_BJetPt_0J0B =
      new TH1D(Form("h_BJetPt_0J0B"), Form("BJet_pT"), 1000, 0, 1000);
  h_BJetPt_0J0B->Sumw2();
  h_BJetEta_0J0B =
      new TH1D(Form("h_BJetEta_0J0B"), Form("BJet_Eta"), 60, -3., 3.);
  h_BJetEta_0J0B->Sumw2();
  h_BJetPhi_0J0B = new TH1D(Form("h_BJetPhi_0J0B"), Form("BJet_Phi"), 60,
                            -3.141594, 3.141594);
  h_BJetPhi_0J0B->Sumw2();

  h_JetPt_0J = new TH1D(Form("h_JetPt_0J"), Form("Jet_pT"), 1000, 0, 1000);
  h_JetPt_0J->Sumw2();
  h_JetEta_0J = new TH1D(Form("h_JetEta_0J"), Form("Jet_Eta"), 60, -3., 3.);
  h_JetEta_0J->Sumw2();
  h_JetPhi_0J =
      new TH1D(Form("h_JetPhi_0J"), Form("Jet_Phi"), 60, -3.141594, 3.141594);
  h_JetPhi_0J->Sumw2();
  h_BJetPt_0J = new TH1D(Form("h_BJetPt_0J"), Form("BJet_pT"), 1000, 0, 1000);
  h_BJetPt_0J->Sumw2();
  h_BJetEta_0J = new TH1D(Form("h_BJetEta_0J"), Form("BJet_Eta"), 60, -3., 3.);
  h_BJetEta_0J->Sumw2();
  h_BJetPhi_0J =
      new TH1D(Form("h_BJetPhi_0J"), Form("BJet_Phi"), 60, -3.141594, 3.141594);
  h_BJetPhi_0J->Sumw2();

  h_JetPt_1J = new TH1D(Form("h_JetPt_1J"), Form("Jet_pT"), 1000, 0, 1000);
  h_JetPt_1J->Sumw2();
  h_JetEta_1J = new TH1D(Form("h_JetEta_1J"), Form("Jet_Eta"), 60, -3., 3.);
  h_JetEta_1J->Sumw2();
  h_JetPhi_1J =
      new TH1D(Form("h_JetPhi_1J"), Form("Jet_Phi"), 60, -3.141594, 3.141594);
  h_JetPhi_1J->Sumw2();
  h_BJetPt_1J = new TH1D(Form("h_BJetPt_1J"), Form("BJet_pT"), 1000, 0, 1000);
  h_BJetPt_1J->Sumw2();
  h_BJetEta_1J = new TH1D(Form("h_BJetEta_1J"), Form("BJet_Eta"), 60, -3., 3.);
  h_BJetEta_1J->Sumw2();
  h_BJetPhi_1J =
      new TH1D(Form("h_BJetPhi_1J"), Form("BJet_Phi"), 60, -3.141594, 3.141594);
  h_BJetPhi_1J->Sumw2();

  h_JetPt_mt1J = new TH1D(Form("h_JetPt_mt1J"), Form("Jet_pT"), 1000, 0, 1000);
  h_JetPt_mt1J->Sumw2();
  h_JetEta_mt1J = new TH1D(Form("h_JetEta_mt1J"), Form("Jet_Eta"), 60, -3., 3.);
  h_JetEta_mt1J->Sumw2();
  h_JetPhi_mt1J =
      new TH1D(Form("h_JetPhi_mt1J"), Form("Jet_Phi"), 60, -3.141594, 3.141594);
  h_JetPhi_mt1J->Sumw2();
  h_BJetPt_mt1J =
      new TH1D(Form("h_BJetPt_mt1J"), Form("BJet_pT"), 1000, 0, 1000);
  h_BJetPt_mt1J->Sumw2();
  h_BJetEta_mt1J =
      new TH1D(Form("h_BJetEta_mt1J"), Form("BJet_Eta"), 60, -3., 3.);
  h_BJetEta_mt1J->Sumw2();
  h_BJetPhi_mt1J = new TH1D(Form("h_BJetPhi_mt1J"), Form("BJet_Phi"), 60,
                            -3.141594, 3.141594);
  h_BJetPhi_mt1J->Sumw2();

  h_JetPt_0BJ = new TH1D(Form("h_JetPt_0BJ"), Form("Jet_pT"), 1000, 0, 1000);
  h_JetPt_0BJ->Sumw2();
  h_JetEta_0BJ = new TH1D(Form("h_JetEta_0BJ"), Form("Jet_Eta"), 60, -3., 3.);
  h_JetEta_0BJ->Sumw2();
  h_JetPhi_0BJ =
      new TH1D(Form("h_JetPhi_0BJ"), Form("Jet_Phi"), 60, -3.141594, 3.141594);
  h_JetPhi_0BJ->Sumw2();
  h_BJetPt_0BJ = new TH1D(Form("h_BJetPt_0BJ"), Form("BJet_pT"), 1000, 0, 1000);
  h_BJetPt_0BJ->Sumw2();
  h_BJetEta_0BJ =
      new TH1D(Form("h_BJetEta_0BJ"), Form("BJet_Eta"), 60, -3., 3.);
  h_BJetEta_0BJ->Sumw2();
  h_BJetPhi_0BJ = new TH1D(Form("h_BJetPhi_0BJ"), Form("BJet_Phi"), 60,
                           -3.141594, 3.141594);
  h_BJetPhi_0BJ->Sumw2();

  h_JetPt_1BJ = new TH1D(Form("h_JetPt_1BJ"), Form("Jet_pT"), 1000, 0, 1000);
  h_JetPt_1BJ->Sumw2();
  h_JetEta_1BJ = new TH1D(Form("h_JetEta_1BJ"), Form("Jet_Eta"), 60, -3., 3.);
  h_JetEta_1BJ->Sumw2();
  h_JetPhi_1BJ =
      new TH1D(Form("h_JetPhi_1BJ"), Form("Jet_Phi"), 60, -3.141594, 3.141594);
  h_JetPhi_1BJ->Sumw2();
  h_BJetPt_1BJ = new TH1D(Form("h_BJetPt_1BJ"), Form("BJet_pT"), 1000, 0, 1000);
  h_BJetPt_1BJ->Sumw2();
  h_BJetEta_1BJ =
      new TH1D(Form("h_BJetEta_1BJ"), Form("BJet_Eta"), 60, -3., 3.);
  h_BJetEta_1BJ->Sumw2();
  h_BJetPhi_1BJ = new TH1D(Form("h_BJetPhi_1BJ"), Form("BJet_Phi"), 60,
                           -3.141594, 3.141594);
  h_BJetPhi_1BJ->Sumw2();

  h_JetPt_mt1BJ =
      new TH1D(Form("h_JetPt_mt1BJ"), Form("Jet_pT"), 1000, 0, 1000);
  h_JetPt_mt1BJ->Sumw2();
  h_JetEta_mt1BJ =
      new TH1D(Form("h_JetEta_mt1BJ"), Form("Jet_Eta"), 60, -3., 3.);
  h_JetEta_mt1BJ->Sumw2();
  h_JetPhi_mt1BJ = new TH1D(Form("h_JetPhi_mt1BJ"), Form("Jet_Phi"), 60,
                            -3.141594, 3.141594);
  h_JetPhi_mt1BJ->Sumw2();
  h_BJetPt_mt1BJ =
      new TH1D(Form("h_BJetPt_mt1BJ"), Form("BJet_pT"), 1000, 0, 1000);
  h_BJetPt_mt1BJ->Sumw2();
  h_BJetEta_mt1BJ =
      new TH1D(Form("h_BJetEta_mt1BJ"), Form("BJet_Eta"), 60, -3., 3.);
  h_BJetEta_mt1BJ->Sumw2();
  h_BJetPhi_mt1BJ = new TH1D(Form("h_BJetPhi_mt1BJ"), Form("BJet_Phi"), 60,
                             -3.141594, 3.141594);
  h_BJetPhi_mt1BJ->Sumw2();

  h_LeadingMuonPt_0J0B =
      new TH1D(Form("h_LeadingMuonPt_0J0B"), Form("Muon_pT"), 1000, 0, 1000);
  h_LeadingMuonPt_0J0B->Sumw2();
  h_LeadingMuonEta_0J0B =
      new TH1D(Form("h_LeadingMuonEta_0J0B"), Form("Muon_Eta"), 60, -3., 3.);
  h_LeadingMuonEta_0J0B->Sumw2();
  h_LeadingMuonPhi_0J0B = new TH1D(Form("h_LeadingMuonPhi_0J0B"),
                                   Form("Muon_Phi"), 60, -3.141594, 3.141594);
  h_LeadingMuonPhi_0J0B->Sumw2();
  h_SubleadingMuonPt_0J0B =
      new TH1D(Form("h_SubleadingMuonPt_0J0B"), Form("Muon_pT"), 1000, 0, 1000);
  h_SubleadingMuonPt_0J0B->Sumw2();
  h_SubleadingMuonEta_0J0B =
      new TH1D(Form("h_SubleadingMuonEta_0J0B"), Form("Muon_Eta"), 60, -3., 3.);
  h_SubleadingMuonEta_0J0B->Sumw2();
  h_SubleadingMuonPhi_0J0B =
      new TH1D(Form("h_SubleadingMuonPhi_0J0B"), Form("Muon_Phi"), 60,
               -3.141594, 3.141594);
  h_SubleadingMuonPhi_0J0B->Sumw2();
  h_MuonPt_0J0B =
      new TH1D(Form("h_MuonPt_0J0B"), Form("Muon_pT"), 1000, 0, 1000);
  h_MuonPt_0J0B->Sumw2();
  h_MuonEta_0J0B =
      new TH1D(Form("h_MuonEta_0J0B"), Form("Muon_Eta"), 60, -3., 3.);
  h_MuonEta_0J0B->Sumw2();
  h_MuonPhi_0J0B = new TH1D(Form("h_MuonPhi_0J0B"), Form("Muon_Phi"), 60,
                            -3.141594, 3.141594);
  h_MuonPhi_0J0B->Sumw2();
  h_dimuonMass_0J0B =
      new TH1D(Form("h_dimuonMass_0J0B"), Form("inv_Mass"), 100, 41., 141.);
  h_dimuonMass_0J0B->Sumw2();
  h_dimuonMass_wide_0J0B =
      new TH1D(Form("h_dimuonMass_wide_0J0B"), Form("inv_Mass"),
               xbins.size() - 1, &(xbins[0]));
  h_dimuonMass_wide_0J0B->Sumw2();
  h_dimuonPt_0J0B =
      new TH1D(Form("h_dimuonPt_0J0B"), Form("dimuon_pT"), 1000, 0., 1000.);
  h_dimuonPt_0J0B->Sumw2();
  h_dimuonRap_0J0B =
      new TH1D(Form("h_dimuonRap_0J0B"), Form("dimuon_rap"), 60, -3., 3.);
  h_dimuonRap_0J0B->Sumw2();

  h_LeadingMuonPt_0J =
      new TH1D(Form("h_LeadingMuonPt_0J"), Form("Muon_pT"), 1000, 0, 1000);
  h_LeadingMuonPt_0J->Sumw2();
  h_LeadingMuonEta_0J =
      new TH1D(Form("h_LeadingMuonEta_0J"), Form("Muon_Eta"), 60, -3., 3.);
  h_LeadingMuonEta_0J->Sumw2();
  h_LeadingMuonPhi_0J = new TH1D(Form("h_LeadingMuonPhi_0J"), Form("Muon_Phi"),
                                 60, -3.141594, 3.141594);
  h_LeadingMuonPhi_0J->Sumw2();
  h_SubleadingMuonPt_0J =
      new TH1D(Form("h_SubleadingMuonPt_0J"), Form("Muon_pT"), 1000, 0, 1000);
  h_SubleadingMuonPt_0J->Sumw2();
  h_SubleadingMuonEta_0J =
      new TH1D(Form("h_SubleadingMuonEta_0J"), Form("Muon_Eta"), 60, -3., 3.);
  h_SubleadingMuonEta_0J->Sumw2();
  h_SubleadingMuonPhi_0J = new TH1D(Form("h_SubleadingMuonPhi_0J"),
                                    Form("Muon_Phi"), 60, -3.141594, 3.141594);
  h_SubleadingMuonPhi_0J->Sumw2();
  h_MuonPt_0J = new TH1D(Form("h_MuonPt_0J"), Form("Muon_pT"), 1000, 0, 1000);
  h_MuonPt_0J->Sumw2();
  h_MuonEta_0J = new TH1D(Form("h_MuonEta_0J"), Form("Muon_Eta"), 60, -3., 3.);
  h_MuonEta_0J->Sumw2();
  h_MuonPhi_0J =
      new TH1D(Form("h_MuonPhi_0J"), Form("Muon_Phi"), 60, -3.141594, 3.141594);
  h_MuonPhi_0J->Sumw2();
  h_dimuonMass_0J =
      new TH1D(Form("h_dimuonMass_0J"), Form("inv_Mass"), 100, 41., 141.);
  h_dimuonMass_0J->Sumw2();
  h_dimuonMass_wide_0J =
      new TH1D(Form("h_dimuonMass_wide_0J"), Form("inv_Mass"), xbins.size() - 1,
               &(xbins[0]));
  h_dimuonMass_wide_0J->Sumw2();
  h_dimuonPt_0J =
      new TH1D(Form("h_dimuonPt_0J"), Form("dimuon_pT"), 1000, 0., 1000.);
  h_dimuonPt_0J->Sumw2();
  h_dimuonRap_0J =
      new TH1D(Form("h_dimuonRap_0J"), Form("dimuon_rap"), 60, -3., 3.);
  h_dimuonRap_0J->Sumw2();

  h_LeadingMuonPt_1J =
      new TH1D(Form("h_LeadingMuonPt_1J"), Form("Muon_pT"), 1000, 0, 1000);
  h_LeadingMuonPt_1J->Sumw2();
  h_LeadingMuonEta_1J =
      new TH1D(Form("h_LeadingMuonEta_1J"), Form("Muon_Eta"), 60, -3., 3.);
  h_LeadingMuonEta_1J->Sumw2();
  h_LeadingMuonPhi_1J = new TH1D(Form("h_LeadingMuonPhi_1J"), Form("Muon_Phi"),
                                 60, -3.141594, 3.141594);
  h_LeadingMuonPhi_1J->Sumw2();
  h_SubleadingMuonPt_1J =
      new TH1D(Form("h_SubleadingMuonPt_1J"), Form("Muon_pT"), 1000, 0, 1000);
  h_SubleadingMuonPt_1J->Sumw2();
  h_SubleadingMuonEta_1J =
      new TH1D(Form("h_SubleadingMuonEta_1J"), Form("Muon_Eta"), 60, -3., 3.);
  h_SubleadingMuonEta_1J->Sumw2();
  h_SubleadingMuonPhi_1J = new TH1D(Form("h_SubleadingMuonPhi_1J"),
                                    Form("Muon_Phi"), 60, -3.141594, 3.141594);
  h_SubleadingMuonPhi_1J->Sumw2();
  h_MuonPt_1J = new TH1D(Form("h_MuonPt_1J"), Form("Muon_pT"), 1000, 0, 1000);
  h_MuonPt_1J->Sumw2();
  h_MuonEta_1J = new TH1D(Form("h_MuonEta_1J"), Form("Muon_Eta"), 60, -3., 3.);
  h_MuonEta_1J->Sumw2();
  h_MuonPhi_1J =
      new TH1D(Form("h_MuonPhi_1J"), Form("Muon_Phi"), 60, -3.141594, 3.141594);
  h_MuonPhi_1J->Sumw2();
  h_dimuonMass_1J =
      new TH1D(Form("h_dimuonMass_1J"), Form("inv_Mass"), 100, 41., 141.);
  h_dimuonMass_1J->Sumw2();
  h_dimuonMass_wide_1J =
      new TH1D(Form("h_dimuonMass_wide_1J"), Form("inv_Mass"), xbins.size() - 1,
               &(xbins[0]));
  h_dimuonMass_wide_1J->Sumw2();
  h_dimuonPt_1J =
      new TH1D(Form("h_dimuonPt_1J"), Form("dimuon_pT"), 1000, 0., 1000.);
  h_dimuonPt_1J->Sumw2();
  h_dimuonRap_1J =
      new TH1D(Form("h_dimuonRap_1J"), Form("dimuon_rap"), 60, -3., 3.);
  h_dimuonRap_1J->Sumw2();

  h_LeadingMuonPt_mt1J =
      new TH1D(Form("h_LeadingMuonPt_mt1J"), Form("Muon_pT"), 1000, 0, 1000);
  h_LeadingMuonPt_mt1J->Sumw2();
  h_LeadingMuonEta_mt1J =
      new TH1D(Form("h_LeadingMuonEta_mt1J"), Form("Muon_Eta"), 60, -3., 3.);
  h_LeadingMuonEta_mt1J->Sumw2();
  h_LeadingMuonPhi_mt1J = new TH1D(Form("h_LeadingMuonPhi_mt1J"),
                                   Form("Muon_Phi"), 60, -3.141594, 3.141594);
  h_LeadingMuonPhi_mt1J->Sumw2();
  h_SubleadingMuonPt_mt1J =
      new TH1D(Form("h_SubleadingMuonPt_mt1J"), Form("Muon_pT"), 1000, 0, 1000);
  h_SubleadingMuonPt_mt1J->Sumw2();
  h_SubleadingMuonEta_mt1J =
      new TH1D(Form("h_SubleadingMuonEta_mt1J"), Form("Muon_Eta"), 60, -3., 3.);
  h_SubleadingMuonEta_mt1J->Sumw2();
  h_SubleadingMuonPhi_mt1J =
      new TH1D(Form("h_SubleadingMuonPhi_mt1J"), Form("Muon_Phi"), 60,
               -3.141594, 3.141594);
  h_SubleadingMuonPhi_mt1J->Sumw2();
  h_MuonPt_mt1J =
      new TH1D(Form("h_MuonPt_mt1J"), Form("Muon_pT"), 1000, 0, 1000);
  h_MuonPt_mt1J->Sumw2();
  h_MuonEta_mt1J =
      new TH1D(Form("h_MuonEta_mt1J"), Form("Muon_Eta"), 60, -3., 3.);
  h_MuonEta_mt1J->Sumw2();
  h_MuonPhi_mt1J = new TH1D(Form("h_MuonPhi_mt1J"), Form("Muon_Phi"), 60,
                            -3.141594, 3.141594);
  h_MuonPhi_mt1J->Sumw2();
  h_dimuonMass_mt1J =
      new TH1D(Form("h_dimuonMass_mt1J"), Form("inv_Mass"), 100, 41., 141.);
  h_dimuonMass_mt1J->Sumw2();
  h_dimuonMass_wide_mt1J =
      new TH1D(Form("h_dimuonMass_wide_mt1J"), Form("inv_Mass"),
               xbins.size() - 1, &(xbins[0]));
  h_dimuonMass_wide_mt1J->Sumw2();
  h_dimuonPt_mt1J =
      new TH1D(Form("h_dimuonPt_mt1J"), Form("dimuon_pT"), 1000, 0., 1000.);
  h_dimuonPt_mt1J->Sumw2();
  h_dimuonRap_mt1J =
      new TH1D(Form("h_dimuonRap_mt1J"), Form("dimuon_rap"), 60, -3., 3.);
  h_dimuonRap_mt1J->Sumw2();

  h_LeadingMuonPt_0BJ =
      new TH1D(Form("h_LeadingMuonPt_0BJ"), Form("Muon_pT"), 1000, 0, 1000);
  h_LeadingMuonPt_0BJ->Sumw2();
  h_LeadingMuonEta_0BJ =
      new TH1D(Form("h_LeadingMuonEta_0BJ"), Form("Muon_Eta"), 60, -3., 3.);
  h_LeadingMuonEta_0BJ->Sumw2();
  h_LeadingMuonPhi_0BJ = new TH1D(Form("h_LeadingMuonPhi_0BJ"),
                                  Form("Muon_Phi"), 60, -3.141594, 3.141594);
  h_LeadingMuonPhi_0BJ->Sumw2();
  h_SubleadingMuonPt_0BJ =
      new TH1D(Form("h_SubleadingMuonPt_0BJ"), Form("Muon_pT"), 1000, 0, 1000);
  h_SubleadingMuonPt_0BJ->Sumw2();
  h_SubleadingMuonEta_0BJ =
      new TH1D(Form("h_SubleadingMuonEta_0BJ"), Form("Muon_Eta"), 60, -3., 3.);
  h_SubleadingMuonEta_0BJ->Sumw2();
  h_SubleadingMuonPhi_0BJ = new TH1D(Form("h_SubleadingMuonPhi_0BJ"),
                                     Form("Muon_Phi"), 60, -3.141594, 3.141594);
  h_SubleadingMuonPhi_0BJ->Sumw2();
  h_MuonPt_0BJ = new TH1D(Form("h_MuonPt_0BJ"), Form("Muon_pT"), 1000, 0, 1000);
  h_MuonPt_0BJ->Sumw2();
  h_MuonEta_0BJ =
      new TH1D(Form("h_MuonEta_0BJ"), Form("Muon_Eta"), 60, -3., 3.);
  h_MuonEta_0BJ->Sumw2();
  h_MuonPhi_0BJ = new TH1D(Form("h_MuonPhi_0BJ"), Form("Muon_Phi"), 60,
                           -3.141594, 3.141594);
  h_MuonPhi_0BJ->Sumw2();
  h_dimuonMass_0BJ =
      new TH1D(Form("h_dimuonMass_0BJ"), Form("inv_Mass"), 100, 41., 141.);
  h_dimuonMass_0BJ->Sumw2();
  h_dimuonMass_wide_0BJ =
      new TH1D(Form("h_dimuonMass_wide_0BJ"), Form("inv_Mass"),
               xbins.size() - 1, &(xbins[0]));
  h_dimuonMass_wide_0BJ->Sumw2();
  h_dimuonPt_0BJ =
      new TH1D(Form("h_dimuonPt_0BJ"), Form("dimuon_pT"), 1000, 0., 1000.);
  h_dimuonPt_0BJ->Sumw2();
  h_dimuonRap_0BJ =
      new TH1D(Form("h_dimuonRap_0BJ"), Form("dimuon_rap"), 60, -3., 3.);
  h_dimuonRap_0BJ->Sumw2();

  h_LeadingMuonPt_1BJ =
      new TH1D(Form("h_LeadingMuonPt_1BJ"), Form("Muon_pT"), 1000, 0, 1000);
  h_LeadingMuonPt_1BJ->Sumw2();
  h_LeadingMuonEta_1BJ =
      new TH1D(Form("h_LeadingMuonEta_1BJ"), Form("Muon_Eta"), 60, -3., 3.);
  h_LeadingMuonEta_1BJ->Sumw2();
  h_LeadingMuonPhi_1BJ = new TH1D(Form("h_LeadingMuonPhi_1BJ"),
                                  Form("Muon_Phi"), 60, -3.141594, 3.141594);
  h_LeadingMuonPhi_1BJ->Sumw2();
  h_SubleadingMuonPt_1BJ =
      new TH1D(Form("h_SubleadingMuonPt_1BJ"), Form("Muon_pT"), 1000, 0, 1000);
  h_SubleadingMuonPt_1BJ->Sumw2();
  h_SubleadingMuonEta_1BJ =
      new TH1D(Form("h_SubleadingMuonEta_1BJ"), Form("Muon_Eta"), 60, -3., 3.);
  h_SubleadingMuonEta_1BJ->Sumw2();
  h_SubleadingMuonPhi_1BJ = new TH1D(Form("h_SubleadingMuonPhi_1BJ"),
                                     Form("Muon_Phi"), 60, -3.141594, 3.141594);
  h_SubleadingMuonPhi_1BJ->Sumw2();
  h_MuonPt_1BJ = new TH1D(Form("h_MuonPt_1BJ"), Form("Muon_pT"), 1000, 0, 1000);
  h_MuonPt_1BJ->Sumw2();
  h_MuonEta_1BJ =
      new TH1D(Form("h_MuonEta_1BJ"), Form("Muon_Eta"), 60, -3., 3.);
  h_MuonEta_1BJ->Sumw2();
  h_MuonPhi_1BJ = new TH1D(Form("h_MuonPhi_1BJ"), Form("Muon_Phi"), 60,
                           -3.141594, 3.141594);
  h_MuonPhi_1BJ->Sumw2();
  h_dimuonMass_1BJ =
      new TH1D(Form("h_dimuonMass_1BJ"), Form("inv_Mass"), 100, 41., 141.);
  h_dimuonMass_1BJ->Sumw2();
  h_dimuonMass_wide_1BJ =
      new TH1D(Form("h_dimuonMass_wide_1BJ"), Form("inv_Mass"),
               xbins.size() - 1, &(xbins[0]));
  h_dimuonMass_wide_1BJ->Sumw2();
  h_dimuonPt_1BJ =
      new TH1D(Form("h_dimuonPt_1BJ"), Form("dimuon_pT"), 1000, 0., 1000.);
  h_dimuonPt_1BJ->Sumw2();
  h_dimuonRap_1BJ =
      new TH1D(Form("h_dimuonRap_1BJ"), Form("dimuon_rap"), 60, -3., 3.);
  h_dimuonRap_1BJ->Sumw2();

  h_LeadingMuonPt_mt1BJ =
      new TH1D(Form("h_LeadingMuonPt_mt1BJ"), Form("Muon_pT"), 1000, 0, 1000);
  h_LeadingMuonPt_mt1BJ->Sumw2();
  h_LeadingMuonEta_mt1BJ =
      new TH1D(Form("h_LeadingMuonEta_mt1BJ"), Form("Muon_Eta"), 60, -3., 3.);
  h_LeadingMuonEta_mt1BJ->Sumw2();
  h_LeadingMuonPhi_mt1BJ = new TH1D(Form("h_LeadingMuonPhi_mt1BJ"),
                                    Form("Muon_Phi"), 60, -3.141594, 3.141594);
  h_LeadingMuonPhi_mt1BJ->Sumw2();
  h_SubleadingMuonPt_mt1BJ = new TH1D(Form("h_SubleadingMuonPt_mt1BJ"),
                                      Form("Muon_pT"), 1000, 0, 1000);
  h_SubleadingMuonPt_mt1BJ->Sumw2();
  h_SubleadingMuonEta_mt1BJ = new TH1D(Form("h_SubleadingMuonEta_mt1BJ"),
                                       Form("Muon_Eta"), 60, -3., 3.);
  h_SubleadingMuonEta_mt1BJ->Sumw2();
  h_SubleadingMuonPhi_mt1BJ =
      new TH1D(Form("h_SubleadingMuonPhi_mt1BJ"), Form("Muon_Phi"), 60,
               -3.141594, 3.141594);
  h_SubleadingMuonPhi_mt1BJ->Sumw2();
  h_MuonPt_mt1BJ =
      new TH1D(Form("h_MuonPt_mt1BJ"), Form("Muon_pT"), 1000, 0, 1000);
  h_MuonPt_mt1BJ->Sumw2();
  h_MuonEta_mt1BJ =
      new TH1D(Form("h_MuonEta_mt1BJ"), Form("Muon_Eta"), 60, -3., 3.);
  h_MuonEta_mt1BJ->Sumw2();
  h_MuonPhi_mt1BJ = new TH1D(Form("h_MuonPhi_mt1BJ"), Form("Muon_Phi"), 60,
                             -3.141594, 3.141594);
  h_MuonPhi_mt1BJ->Sumw2();
  h_dimuonMass_mt1BJ =
      new TH1D(Form("h_dimuonMass_mt1BJ"), Form("inv_Mass"), 100, 41., 141.);
  h_dimuonMass_mt1BJ->Sumw2();
  h_dimuonMass_wide_mt1BJ =
      new TH1D(Form("h_dimuonMass_wide_mt1BJ"), Form("inv_Mass"),
               xbins.size() - 1, &(xbins[0]));
  h_dimuonMass_wide_mt1BJ->Sumw2();
  h_dimuonPt_mt1BJ =
      new TH1D(Form("h_dimuonPt_mt1BJ"), Form("dimuon_pT"), 1000, 0., 1000.);
  h_dimuonPt_mt1BJ->Sumw2();
  h_dimuonRap_mt1BJ =
      new TH1D(Form("h_dimuonRap_mt1BJ"), Form("dimuon_rap"), 60, -3., 3.);
  h_dimuonRap_mt1BJ->Sumw2();

  // sample

  // std::vector<float> xbins = { 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64,
  // 68, 72, 76, 81, 86,
  //                              91, 96, 101, 106, 110, 115, 120, 126, 133,
  //                              141, 150, 160, 171, 185, 200, 220, 243, 273,
  //                              320, 380, 440, 510, 600, 700, 830, 1000, 1500,
  //                              3000};

  // qcut_pT           = new TH1D(Form("qcut_pT"),Form("Muon_pT"), 1000, 0,
  // 1000); qcut_pT->Sumw2(); qcut_eta          = new
  // TH1D(Form("qcut_eta"),Form("Muon_eta"), 60, -3., 3.); qcut_eta->Sumw2();
  // qcut_phi          = new TH1D(Form("qcut_phi"),Form("Muon_phi"), 60,
  // -3.141593, 3.141593); qcut_phi->Sumw2();

  // muon_mass           = new TH1D(Form("muon_mass"),Form("Muon_mass"),
  // 200, 41., 141.); muon_mass->Sumw2(); muon_mass_wide          = new
  // TH1D(Form("muon_mass_wide"),Form("Muon_mass"), xbins.size()-1,
  // &(xbins[0])); muon_mass_wide->Sumw2();
}

void ssb_analysis::End() {
  fout->Write();
  fout->Close();
}

void ssb_analysis::SetOutputFileName(char *outname) { outfile = outname; }
