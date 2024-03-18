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

// TClonesArray    *Muon;
// TClonesArray    *GenMuon;
// vector<int>     *Muon_Charge;
// Int_t           Muon_Count;
// vector<double>  *Muon_PFIsodBeta03;
// vector<double>  *Muon_PFIsodBeta04;
// vector<bool>    *Muon_isHighPt;
// vector<bool>    *Muon_isLoose;
// vector<bool>    *Muon_isMedium;
// vector<bool>    *Muon_isMedium2016;
// vector<bool>    *Muon_isSoft;
// vector<bool>    *Muon_isTight;
// vector<int>     *Muon_pdgId;
// vector<double>  *Muon_rand1;
// vector<double>  *Muon_rand2;
// vector<double>  *Muon_relIso03;
// vector<double>  *Muon_relIso04;
// vector<int>     *Muon_trackerLayers;

// Muon 
// Iso :  Muon_PFIsodBeta04 < 0.15 
// Eta : |eta| < 2.4 
// pT : pT > 20 GeV
// iD : Muon_isTight == true 


// <MC>
// 8 HLT_IsoMu24_v4 1 0 1 1
// 10 HLT_IsoTkMu24_v4 1 0 1 1
// 11 HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v6 1 0 1 1
// 12 HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7 1 0 1 1
// 13 HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5 1 0 1 1
// 14 HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6 1 0 1 1

// <DATA>
// 7 HLT_IsoMu24_v1 1 0 0 1
// 8 HLT_IsoTkMu24_v1 1 0 0 1
// 9 HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v2 1 0 0 1
// 10 HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2 1 0 0 1
// 11 HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v2 1 0 0 1
// 12 HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2 1 0 0 1

bool inAcceptance(const TLorentzVector* aMuon, float pTcut, float etaCut) {

   return ( aMuon->Pt() > pTcut && std::abs(aMuon->Eta()) < etaCut );
} 

bool isTightID(const bool isTight) {

   return isTight;
}

bool cutPFRelIso(const double pfIso, double cutPF) {


   return pfIso < cutPF;
}

std::vector<TLorentzVector*> sortByPt(const std::vector<TLorentzVector*>& vec) {

   if (vec.size() <= 1) 
      return vec;

   std::vector<TLorentzVector*> retrunVec;
   std::vector<TLorentzVector*> vecIter = vec;

   for ( int i = 0; i < vec.size(); i++ ) {
      TLorentzVector* aVec = vecIter.at(0);
      int delIdx = 0;

      for (int j = 0; j < vecIter.size(); j++ ) {

         if ( aVec->Pt() < vecIter.at(j)->Pt() ) {
            aVec = vecIter.at(j);
            delIdx = j;
         }

      }

      retrunVec.push_back(aVec);
      vecIter.erase(vecIter.begin() + delIdx);
   }

   return retrunVec;
}



void ssb_analysis::Loop( char *logfile )
{

   //////////
   if (fChain == 0) return;
   //////////

   //////////
   Long64_t nentries = fChain->GetEntriesFast();
   // std::cout << nentries << std::endl; 
   // return;   

   Long64_t nbytes = 0, nb = 0;
   //////////

   ///My variables
   Long64_t __tot_evt = 0;

   /// Check Total Event
 
   ////////////////////////
   /// start event loop ///
   ////////////////////////

   bool isMC = true;

   int index_IsoMu24 = 0;
   int index_IsoTkMu24 = 0;

   if ( isMC ) {
      index_IsoMu24   = 8;
      index_IsoTkMu24 = 10;
   } else {
      index_IsoMu24   = 7;
      index_IsoTkMu24 = 8;
   }

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

      // if (__tot_evt == 1000)
      //    return;


      ////////////////////////////////////////
      /// start Main Loop Function and Cut ///
      ////////////////////////////////////////

      double eventWeight = 1;

      if ( !(Trigger_isPass->at(index_IsoMu24) || Trigger_isPass->at(index_IsoTkMu24)) )
         continue;

      std::vector<TLorentzVector*> passMuon;
      std::vector<int> chargeMuon;

      TLorentzVector* aMuon;
      for(int i = 0; i < Muon->GetEntries(); ++i ) {

         aMuon = (TLorentzVector*)Muon->At(i);

         FillHisto(h_cf_Muon_pT, aMuon->Pt());
         FillHisto(h_cf_Muon_Eta, aMuon->Eta());
         FillHisto(h_cf_Muon_Phi, aMuon->Phi());

         if ( !inAcceptance(aMuon, 12., 2.4) )
            continue;

         if ( !isTightID(Muon_isTight->at(i)) )
            continue;

         if ( !cutPFRelIso(Muon_PFIsodBeta04->at(i), 0.15) )
            continue;

         chargeMuon.push_back(Muon_Charge->at(i));
         passMuon.push_back(aMuon);

      }

      int nPos = 0;
      int nNeg = 0;

      if ( chargeMuon.size() == 0 || passMuon.size() == 0 )
         continue;

      if ( !(passMuon.at(0)->Pt() > 26.) )
         continue;

      for (int i = 0; i < chargeMuon.size(); i++) {
         if ( chargeMuon.at(i) < 0 )
            nNeg++;

         if ( chargeMuon.at(i) > 0 )
            nPos++;
      }

      if (nPos < 1 || nNeg < 1)
         continue;

      TLorentzVector leadingMuon = *( passMuon.data()[0] );
      TLorentzVector subleadingMuon; // = *( passMuon.data()[1] );

      for ( int i = 1; i < chargeMuon.size(); i++ ) {
         if (chargeMuon.at(i) * chargeMuon.at(0) < 0) {
            subleadingMuon = *( passMuon.data()[i] );
            break;
         }
      }



      FillHisto(h_cf_Muon_pT_pair, leadingMuon.Pt());
      FillHisto(h_cf_Muon_Eta_pair, leadingMuon.Eta());
      FillHisto(h_cf_Muon_Phi_pair, leadingMuon.Phi());
      FillHisto(h_cf_Muon_pT_pair, subleadingMuon.Pt());
      FillHisto(h_cf_Muon_Eta_pair, subleadingMuon.Eta());
      FillHisto(h_cf_Muon_Phi_pair, subleadingMuon.Phi());

      FillHisto(h_cf_Muon_pT_leading, leadingMuon.Pt());
      FillHisto(h_cf_Muon_Eta_leading, leadingMuon.Eta());
      FillHisto(h_cf_Muon_Phi_leading, leadingMuon.Phi());

      FillHisto(h_cf_Muon_pT_subleading, subleadingMuon.Pt());
      FillHisto(h_cf_Muon_Eta_subleading, subleadingMuon.Eta());
      FillHisto(h_cf_Muon_Phi_subleading, subleadingMuon.Phi());


      FillHisto(h_cf_Muon_Mass_pair, (leadingMuon + subleadingMuon).M() );
      FillHisto(h_cf_Muon_Mass_wide_pair, (leadingMuon + subleadingMuon).M() );
   }//event loop
  
   printf("Total processed number of events: %lld\n", __tot_evt);


}//end Loop function

void ssb_analysis::Start( int genLoopon )
{
   if      ( genLoopon == 0 ){ fout = new TFile(Form("output/%s",outfile),"RECREATE");}
   else if      ( genLoopon == 1 ){ fout = new TFile(Form("output/%s",outfile),"UPDATE");}
   else {cout << "genLoopon error" << endl;}
   fout->cd("");

   TDirectory *dir = gDirectory;
   dir->cd();

   DeclareHistos();

}

void ssb_analysis::DeclareHistos()
{

   std::vector<float> xbins = { 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 
                                91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 
                                185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1500, 3000};

   /// Test For Systematic All-in-One Code ///
   h_cf_Muon_pT             = new TH1D(Form("_h_cf_Muon_pT_"),Form("Muon pT"), 1000, 0, 1000); h_cf_Muon_pT->Sumw2();
   h_cf_Muon_Eta            = new TH1D(Form("_h_cf_Muon_Eta_"),Form("Muon Eta"), 60, -3., 3.); h_cf_Muon_Eta->Sumw2();
   h_cf_Muon_Phi            = new TH1D(Form("_h_cf_Muon_Phi_"),Form("Muon Phi"), 60, -3.141593, 3.141593); h_cf_Muon_Phi->Sumw2();

   h_cf_Muon_pT_inAccep     = new TH1D(Form("_h_cf_Muon_pT_inAccep"),Form("Muon pT"), 1000, 0, 1000); h_cf_Muon_pT_inAccep->Sumw2();
   h_cf_Muon_Eta_inAccep    = new TH1D(Form("_h_cf_Muon_Eta_inAccep"),Form("Muon Eta"), 60, -3., 3.); h_cf_Muon_Eta_inAccep->Sumw2();
   h_cf_Muon_Phi_inAccep    = new TH1D(Form("_h_cf_Muon_Phi_inAccep"),Form("Muon Phi"), 60, -3.141593, 3.141593); h_cf_Muon_Phi_inAccep->Sumw2();

   h_cf_Muon_pT_passID      = new TH1D(Form("_h_cf_Muon_pT_passID"),Form("Muon pT"), 1000, 0, 1000); h_cf_Muon_pT_passID->Sumw2();
   h_cf_Muon_Eta_passID     = new TH1D(Form("_h_cf_Muon_Eta_passID"),Form("Muon Eta"), 60, -3., 3.); h_cf_Muon_Eta_passID->Sumw2();
   h_cf_Muon_Phi_passID     = new TH1D(Form("_h_cf_Muon_Phi_passID"),Form("Muon Phi"), 60, -3.141593, 3.141593); h_cf_Muon_Phi_passID->Sumw2();

   h_cf_Muon_pT_passISO     = new TH1D(Form("_h_cf_Muon_pT_passISO"),Form("Muon pT"), 1000, 0, 1000); h_cf_Muon_pT_passISO->Sumw2();
   h_cf_Muon_Eta_passISO    = new TH1D(Form("_h_cf_Muon_Eta_passISO"),Form("Muon Eta"), 60, -3., 3.); h_cf_Muon_Eta_passISO->Sumw2();
   h_cf_Muon_Phi_passISO    = new TH1D(Form("_h_cf_Muon_Phi_passISO"),Form("Muon Phi"), 60, -3.141593, 3.141593); h_cf_Muon_Phi_passISO->Sumw2();

   h_cf_Muon_pT_pair        = new TH1D(Form("_h_cf_Muon_pT_pair"),Form("Muon pT"), 1000, 0, 1000); h_cf_Muon_pT_pair->Sumw2();
   h_cf_Muon_Eta_pair       = new TH1D(Form("_h_cf_Muon_Eta_pair"),Form("Muon Eta"), 60, -3., 3.); h_cf_Muon_Eta_pair->Sumw2();
   h_cf_Muon_Phi_pair       = new TH1D(Form("_h_cf_Muon_Phi_pair"),Form("Muon Phi"), 60, -3.141593, 3.141593); h_cf_Muon_Phi_pair->Sumw2();

   h_cf_Muon_Mass_pair      = new TH1D(Form("_h_cf_Muon_Mass_pair"),Form("diMuon inv. Mass"), 200, 41, 141); h_cf_Muon_Mass_pair->Sumw2();
   h_cf_Muon_Mass_wide_pair = new TH1D(Form("_h_cf_Muon_Mass_wide_pair"),Form("diMuon inv. Mass"), xbins.size()-1, &(xbins[0])); h_cf_Muon_Mass_wide_pair->Sumw2();

   h_cf_Muon_pT_leading     = new TH1D(Form("_h_cf_Muon_pT_leading"),Form("Muon pT"), 1000, 0, 1000); h_cf_Muon_pT_leading->Sumw2();
   h_cf_Muon_Eta_leading    = new TH1D(Form("_h_cf_Muon_Eta_leading"),Form("Muon Eta"), 60, -3., 3.); h_cf_Muon_Eta_leading->Sumw2();
   h_cf_Muon_Phi_leading    = new TH1D(Form("_h_cf_Muon_Phi_leading"),Form("Muon Phi"), 60, -3.141593, 3.141593); h_cf_Muon_Phi_leading->Sumw2();

   h_cf_Muon_pT_subleading  = new TH1D(Form("_h_cf_Muon_pT_subleading"),Form("Muon pT"), 1000, 0, 1000); h_cf_Muon_pT_subleading->Sumw2();
   h_cf_Muon_Eta_subleading = new TH1D(Form("_h_cf_Muon_Eta_subleading"),Form("Muon Eta"), 60, -3., 3.); h_cf_Muon_Eta_subleading->Sumw2();
   h_cf_Muon_Phi_subleading = new TH1D(Form("_h_cf_Muon_Phi_subleading"),Form("Muon Phi"), 60, -3.141593, 3.141593); h_cf_Muon_Phi_subleading->Sumw2();



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


