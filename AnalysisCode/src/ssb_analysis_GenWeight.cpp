#define ssb_analysis_cxx

#include <iostream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "TMath.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"

#include "./../interface/ssb_analysis.hpp"
#include "./../CommonTools.hpp"

using namespace std;


void ssb_analysis::Loop( char *logfile )
{

   //////////
   if (fChain == 0) return;
   //////////

   //////////
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   //////////

   ///My variables
   Long64_t __tot_evt = 0;

   /// Check Total Event
 
   ////////////////////////
   /// start event loop ///
   ////////////////////////

   // float globalWeight = 1;
   // float normSF = 1;

   // if ( !fIsData ) {
   //    normSF = GetNormalization();
   //    globalWeight *= normSF;
   // }

   // // std::cout << fProcessName << " "
   // //           << fIsData      << " "
   // //           << normSF       << " "
   // //           << globalWeight << " "
   // //           << std::endl;

   double totalGenWeight = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
      {
         printf("ERROR: Could not load tree!!!\n");
         break;
      }

      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
      if (jentry % 10000 == 0) printf("Event %lld\n", jentry); //%lld supports Long64_t

      __tot_evt++;

      if (Gen_EventWeight != 1.) 
         std::cout << jentry << " " << Gen_EventWeight << std::endl;

      totalGenWeight += Gen_EventWeight;


   //    globalWeight = 1;
   //    globalWeight *= normSF;

   //    if ( !fIsData )
   //       globalWeight *= GetPUweight(PileUp_Count_Intime);

   //    // std::cout << jentry << " "
   //    //           << normSF << " "
   //    //           << GetPUweight(PileUp_Count_Intime) << " "
   //    //           << globalWeight << " "
   //    //           << std::endl;

   //    ////////////////////////////////////////
   //    /// start Main Loop Function and Cut ///
   //    ////////////////////////////////////////

   //    int index_TrigIsoMu24 = 0;
   //    int index_TrigIsoTkMu24 = 0;

   //    for (int i = 0; i < Trigger_Name->size(); i++) {
   //       if (Trigger_Name->at(i).find("IsoMu24") != std::string::npos)
   //          index_TrigIsoMu24 = i;

   //       if (Trigger_Name->at(i).find("IsoTkMu24") != std::string::npos)
   //          index_TrigIsoTkMu24 = i;
   //    }

   //    // std::cout << index_TrigIsoMu24 << " " << Trigger_Name->at(index_TrigIsoMu24) << " " 
   //    //           << index_TrigIsoTkMu24 << " " << Trigger_Name->at(index_TrigIsoTkMu24) << std::endl;

   //    if ( !(Trigger_isPass->at(index_TrigIsoMu24) || Trigger_isPass->at(index_TrigIsoMu24)) )
   //       continue;


   //    std::vector<TLorentzVector*> recoMuon;
   //    recoMuon.reserve(Muon->GetEntries());
   //    for ( int i = 0; i < Muon->GetEntries(); ++i ) {
   //       auto aMuon = (TLorentzVector*)Muon->At(i);
   //       recoMuon.emplace_back(aMuon);
   //    }

   //    std::vector<TLorentzVector*> genlMuon;
   //    genlMuon.reserve(GenMuon->GetEntries());
   //    for ( int i = 0; i < GenMuon->GetEntries(); ++i ) {
   //       auto aMuon = (TLorentzVector*)GenMuon->At(i);
   //       genlMuon.emplace_back(aMuon);
   //    }
      
   //    std::vector<TLorentzVector> recoStdMuon;
   //    recoStdMuon.reserve(recoMuon.size());
   //    for ( int i = 0; i < recoMuon.size(); i++) 
   //       recoStdMuon.emplace_back(*(recoMuon.data()[i]));

   //    std::vector<TLorentzVector> genStdMuon;
   //    genStdMuon.reserve(genlMuon.size());
   //    for ( int i = 0; i < genlMuon.size(); i++) 
   //       genStdMuon.emplace_back(*(genlMuon.data()[i]));


   //    if (fRoccor_enalbed) {
   //       for (int i = 0; i < recoStdMuon.size(); i++) {
   //          if ( fIsData ) {
   //             recoStdMuon.at(i) *= fRoccoR->kScaleDT(Muon_Charge->at(i), 
   //                                                    recoStdMuon.at(i).Pt(), 
   //                                                    recoStdMuon.at(i).Eta(), 
   //                                                    recoStdMuon.at(i).Phi(), 
   //                                                    0, 
   //                                                    0);
   //          } else {

   //             double drmin=999.;
   //             bool match =false;
   //             int genMuonIdx = 0;

   //             for ( int j = 0; j < genStdMuon.size(); j++ ) {
   //                if( recoStdMuon.at(i).DeltaR(genStdMuon.at(j)) < 0.1 && recoStdMuon.at(i).DeltaR(genStdMuon.at(j)) < drmin ) {
   //                   match = true;
   //                   genMuonIdx = j;
   //                   drmin = recoStdMuon.at(i).DeltaR(genStdMuon.at(j));
   //                }
   //             }

   //             if ( match ) {
   //                recoStdMuon.at(i) *= fRoccoR->kScaleFromGenMC(Muon_Charge->at(i), 
   //                                                              recoStdMuon.at(i).Pt(), 
   //                                                              recoStdMuon.at(i).Eta(), 
   //                                                              recoStdMuon.at(i).Phi(), 
   //                                                              Muon_trackerLayers->at(i), 
   //                                                              genStdMuon.at(genMuonIdx).Pt(), 
   //                                                              Muon_rand2->at(i), 
   //                                                              0, 
   //                                                              0);
   //             } else {
   //                recoStdMuon.at(i) *= fRoccoR->kScaleAndSmearMC(Muon_Charge->at(i), 
   //                                                               recoStdMuon.at(i).Pt(),
   //                                                               recoStdMuon.at(i).Eta(),
   //                                                               recoStdMuon.at(i).Phi(),
   //                                                               Muon_trackerLayers->at(i),
   //                                                               Muon_rand1->at(i),
   //                                                               Muon_rand2->at(i),
   //                                                               0,
   //                                                               0);
   //             }
   //          }
   //       }
   //    }

   //    std::vector<TLorentzVector> setMuonsPassingCut;
   //    std::vector<int> setMuonsChargePassingCut;

   //    TLorentzVector aMuon;
   //    for( int i = 0; i < recoStdMuon.size(); i++ )
   //    {
   //       aMuon = recoStdMuon.at(i);

   //       if ( !(aMuon.Pt() > 20.) )
   //          continue;

   //       if ( !( aMuon.Eta() < 2.4 && aMuon.Eta() > -2.4) )
   //          continue;

   //       if ( !(Muon_isTight->at(i) == 1.) )
   //          continue;

   //       if ( !(Muon_PFIsodBeta04->at(i) < 0.15) )
   //          continue;

   //       setMuonsPassingCut.push_back(aMuon);
   //       setMuonsChargePassingCut.push_back(Muon_Charge->at(i));
   //    }

   

   //    int pos = 0;
   //    int neg = 0;

   //    for ( int i = 0; i < setMuonsChargePassingCut.size(); i++ ) {
   //       if ( setMuonsChargePassingCut.at(i) > 0. ) {
   //          pos++;
   //       } else {
   //          neg++;
   //       }
   //    }

   //    if ( pos == 0 || neg == 0 )
   //       continue;

   //    if ( !(setMuonsPassingCut.at(0).Pt() > 26.) )
   //       continue;

   //    int idx_subleading = 0;
   //    for ( int i = 1; i < setMuonsChargePassingCut.size(); i++ ) {
   //       if ( !(setMuonsChargePassingCut.at(0) == setMuonsChargePassingCut.at(i)) ) {
   //          idx_subleading = i;
   //          break;
   //       }
   //    }

   //    TLorentzVector leadingMuon    = setMuonsPassingCut.at(0);
   //    TLorentzVector subleadingMuon = setMuonsPassingCut.at(idx_subleading);

   //    // if (fIDISO_enalbed) {

   //    // }

   //    // if (fTRIGG_enalbed) {
         
   //    // }

   //    TLorentzVector diMuon         = leadingMuon + subleadingMuon;

   //    FillHisto(h_LeadingMuonPt     , leadingMuon.Pt()    , globalWeight);
   //    FillHisto(h_MuonPt            , leadingMuon.Pt()    , globalWeight);
   //    FillHisto(h_MuonEta           , leadingMuon.Eta()    , globalWeight);
   //    FillHisto(h_MuonPhi           , leadingMuon.Phi()    , globalWeight);

   //    FillHisto(h_SubleadingMuonPt  , subleadingMuon.Pt() , globalWeight);
   //    FillHisto(h_MuonPt            , subleadingMuon.Pt() , globalWeight);
   //    FillHisto(h_MuonEta           , subleadingMuon.Eta() , globalWeight);
   //    FillHisto(h_MuonPhi           , subleadingMuon.Phi() , globalWeight);

   //    FillHisto(h_dimuonMass        , diMuon.M()          , globalWeight);
   //    FillHisto(h_dimuonMass_wide   , diMuon.M()          , globalWeight);
   //    FillHisto(h_dimuonPt          , diMuon.Pt()         , globalWeight);
   //    FillHisto(h_dimuonRap         , diMuon.Rapidity()   , globalWeight);

   //    FillHisto(h_PV_Count_before_corr                 , PV_Count                 , normSF);
   //    FillHisto(h_PileUp_Count_Interaction_before_corr , PileUp_Count_Interaction , normSF);
   //    FillHisto(h_PileUp_Count_Intime_before_corr      , PileUp_Count_Intime      , normSF);

   //    FillHisto(h_PV_Count_after_corr                  , PV_Count                 , globalWeight);
   //    FillHisto(h_PileUp_Count_Interaction_after_corr  , PileUp_Count_Interaction , globalWeight);
   //    FillHisto(h_PileUp_Count_Intime_after_corr       , PileUp_Count_Intime      , globalWeight);

   }//event loop
  
   // printf("Total processed number of events: %lld out of %f \n", __tot_evt, totalGenWeight);
   // printf("Total Gen Weigh : %f out of %lld nentries\n", totalGenWeight, nentries);
   // printf("Total Gen Weigh : %f out of %lld __tot_evt\n", totalGenWeight, __tot_evt);
   std::cout << "Total Gen Weigh : " << totalGenWeight << " out of " << __tot_evt << std::endl;


}//end Loop function

void ssb_analysis::Start( int genLoopon )
{
   if      ( genLoopon == 0 ){ fout = new TFile(Form("output/%s",outfile),"RECREATE");}
   else if      ( genLoopon == 1 ){ fout = new TFile(Form("output/%s",outfile),"UPDATE");}
   else {cout << "genLoopon error" << endl;}
   fout->cd("");

   std::cout << "is Data?     : " << GetIsData() << std::endl;
   std::cout << "Process Name : " << GetProcessName() << std::endl;

   TDirectory *dir = gDirectory;
   dir->cd();

   DeclareHistos();

}

void ssb_analysis::DeclareHistos()
{

   // Test For Systematic All-in-One Code

   // h_LeadingMuonPt        = new TH1D(Form("h_LeadingMuonPt"),Form("Muon_pT"), 1000, 0, 1000); h_LeadingMuonPt->Sumw2();
   // h_SubleadingMuonPt     = new TH1D(Form("h_SubleadingMuonPt"),Form("Muon_pT"), 1000, 0, 1000); h_SubleadingMuonPt->Sumw2();
   // h_MuonPt               = new TH1D(Form("h_MuonPt"),Form("Muon_pT"), 1000, 0, 1000); h_MuonPt->Sumw2();
   // h_MuonEta              = new TH1D(Form("h_MuonEta"),Form("Muon_Eta"), 60, -3., 3.); h_MuonEta->Sumw2();
   // h_MuonPhi              = new TH1D(Form("h_MuonPhi"),Form("Muon_Phi"), 60, -3.141593, 3.141593); h_MuonPhi->Sumw2();

   // std::vector<float> xbins = { 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 
   //                              91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 
   //                              185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

   // h_dimuonMass           = new TH1D(Form("h_dimuonMass"),Form("inv_Mass"), 100, 41., 141.); h_dimuonMass->Sumw2();
   // h_dimuonMass_wide      = new TH1D(Form("h_dimuonMass_wide"),Form("inv_Mass"), xbins.size()-1, &(xbins[0])); h_dimuonMass_wide->Sumw2();

   // h_dimuonPt             = new TH1D(Form("h_dimuonPt"),Form("dimuon_pT"), 1000, 0., 1000.); h_dimuonPt->Sumw2();
   // h_dimuonRap            = new TH1D(Form("h_dimuonRap"),Form("dimuon_rap"), 60, -3., 3.); h_dimuonRap->Sumw2();

   // h_PV_Count_before_corr                  = new TH1D(Form("h_PV_Count_before_PUcorr"), Form("PV_Count"), 100, 0., 100.); h_PV_Count_before_corr->Sumw2();
   // h_PileUp_Count_Interaction_before_corr  = new TH1D(Form("h_PileUp_Count_Interaction_before_PUcorr"), Form("PileUp_Count_Interaction"), 100, 0., 100.); h_PileUp_Count_Interaction_before_corr->Sumw2();
   // h_PileUp_Count_Intime_before_corr       = new TH1D(Form("h_PileUp_Count_Intime_before_PUcorr"), Form("PileUp_Count_Intime"), 100, 0., 100.); h_PileUp_Count_Intime_before_corr->Sumw2();

   // h_PV_Count_after_corr                  = new TH1D(Form("h_PV_Count_after_PUcorr"), Form("PV_Count"), 100, 0., 100.); h_PV_Count_after_corr->Sumw2();
   // h_PileUp_Count_Interaction_after_corr  = new TH1D(Form("h_PileUp_Count_Interaction_after_PUcorr"), Form("PileUp_Count_Interaction"), 100, 0., 100.); h_PileUp_Count_Interaction_after_corr->Sumw2();
   // h_PileUp_Count_Intime_after_corr       = new TH1D(Form("h_PileUp_Count_Intime_after_PUcorr"), Form("PileUp_Count_Intime"), 100, 0., 100.); h_PileUp_Count_Intime_after_corr->Sumw2();




   // sample

   // std::vector<float> xbins = { 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 
   //                              91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 
   //                              185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

   // qcut_pT           = new TH1D(Form("qcut_pT"),Form("Muon_pT"), 1000, 0, 1000); qcut_pT->Sumw2();
   // qcut_eta          = new TH1D(Form("qcut_eta"),Form("Muon_eta"), 60, -3., 3.); qcut_eta->Sumw2();
   // qcut_phi          = new TH1D(Form("qcut_phi"),Form("Muon_phi"), 60, -3.141593, 3.141593); qcut_phi->Sumw2();

   // muon_mass           = new TH1D(Form("muon_mass"),Form("Muon_mass"), 200, 41., 141.); muon_mass->Sumw2();
   // muon_mass_wide          = new TH1D(Form("muon_mass_wide"),Form("Muon_mass"), xbins.size()-1, &(xbins[0])); muon_mass_wide->Sumw2();

}


void ssb_analysis::End()
{
   fout->Write();
   fout->Close();
}

void ssb_analysis::SetOutputFileName(char *outname)
{   
   outfile = outname;
}


