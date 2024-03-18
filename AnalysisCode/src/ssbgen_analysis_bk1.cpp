#define ssbgen_analysis_cxx

#include <iostream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "TMath.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"

#include "./../interface/ssbgen_analysis.hpp"
#include "./../CommonTools.hpp"

using namespace std;

void ssbgen_analysis::GetVariables()
{

}
void ssbgen_analysis::TLVInitial()
{

}


void ssbgen_analysis::Loop()
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

   num_dilep =0;
   num_semilep =0;
   num_had =0;
   num_non =0;
   num_dimu =0;
   num_diel =0;
   num_muel =0;
   num_dimu_v1 =0;
   num_diel_v1 =0;
   num_muel_v1 =0;
   num_mutau_v1 =0;
   num_eltau_v1 =0;
   num_ditau_v1 =0;
   ////////////////////////
   /// start event loop ///
   ////////////////////////
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
      {
         printf("ERROR: Could not load tree!!!\n");
         break;
      }

      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
      if (jentry % 100000 == 0) printf("Event %lld\n", jentry); //%lld supports Long64_t

      __tot_evt++;

      ////////////////////////////////////////
      /// start Main Loop Function and Cut ///
      ////////////////////////////////////////

      TLVInitial();

      /////////////////////////////
      /// Find Dilepton Channel ///
      /////////////////////////////
      num_lep = 0;

      if ( Channel_Idx_Final == 26 ) {num_dimu++;}
      if ( Channel_Idx_Final == 22 ) {num_diel++;}
      if ( Channel_Idx_Final == 24 ) {num_muel++;}

      if ( Channel_Idx == 26 ) {num_dimu_v1++;}  // not final mumu
      if ( Channel_Idx == 22 ) {num_diel_v1++;}  // not final mumu
      if ( Channel_Idx == 24 ) {num_muel_v1++;}  // not final mumu

      if ( Channel_Idx == -2 ) {num_mutau_v1++;}  // not final mumu
      if ( Channel_Idx == -4 ) {num_eltau_v1++;}  // not final mumu
      if ( Channel_Idx == -30 ) {num_ditau_v1++;}  // not final mumu

      if ( Channel_Lepton_Count == 2 ) {num_dilep++;}
      if ( Channel_Lepton_Count == 1 ) {num_semilep++;}
      if ( Channel_Lepton_Count == 0 ) {num_had++;}
/*      if ( Channel_Lepton_Count == 2 ) {num_dilep++;

//         if ( Info_EventNumber == 81928690 ){
            cout << "Info_EventNumber  ? " <<  Info_EventNumber << endl;
            for ( int i = 0; i < GenPar_Status->size(); ++i )
            {
               cout <<"Idx? "<< GenPar_Idx->at(i) << " status? " << GenPar_Status->at(i) << " pdgId ? " << GenPar_pdgId->at(i)
                    << " Mom1Idx? " << GenPar_Mom1_Idx->at(i) << " Mom2Idx? " << GenPar_Mom2_Idx->at(i) << endl;
            }
            cout << "-- -- -- -- -- -- " << endl;
//         }
      }
      else if ( Channel_Lepton_Count == 1 ) {num_semilep++;}
      else if ( Channel_Lepton_Count == 0 ) {num_had++;}
      else {cout << "somthing wrong"<< endl; cout << "Channel_Lepton_Count_Final ? " << Channel_Lepton_Count << endl;}*/

      //////////////////////
      /// Fill Histogram ///
      //////////////////////

   }//event loop

   cout << "num_had ?  " << num_had << endl;
   cout << "num_semilep ?  " << num_semilep << endl;
   cout << "num_dilep ?  " << num_dilep << endl;
   cout << "num_non ?  " << num_non << endl;

   cout << "num_dimu ?  " << num_dimu << endl;
   cout << "num_diel ?  " << num_diel << endl;
   cout << "num_muel ?  " << num_muel << endl;

   cout << "num_dimu_v1 ?  " << num_dimu_v1 << endl;
   cout << "num_diel_v1 ?  " << num_diel_v1 << endl;
   cout << "num_muel_v1 ?  " << num_muel_v1 << endl;

   cout << "num_mutau_v1 ?  " << num_mutau_v1 << endl;
   cout << "num_eltau_v1 ?  " << num_eltau_v1 << endl;
   cout << "num_ditau_v1 ?  " << num_ditau_v1 << endl;

   cout << " MC sf v1 ? " << 831.76*0.011*1000/num_dimu_v1 << endl;
   cout << " MC sf ?  "   << 831.76*0.011*1000/num_dimu    << endl;

   printf("before Met Filter Total processed number of events: %lld\n", __tot_evt);
   
}//end Loop function


void ssbgen_analysis::GetNtupleTotalEvent( unsigned int totevent )// basic frame
{
   NtupletotalEvent = totevent;
}

void ssbgen_analysis::Start()
{
   fout = new TFile(Form("output/%s",outfile),"RECREATE");
   fout->cd("");

   TDirectory *dir = gDirectory;
   dir->cd();

   DeclareHistos();

}

void ssbgen_analysis::DeclareHistos()
{
   double  pi = TMath::Pi();

   Channel_Had = new TH1F(Form("Channel_Had" ), Form(" Num. leptons "), 10, 0.0, 10); 
   Channel_Semi = new TH1F(Form("Channel_Semi" ), Form(" Num. leptons "), 10, 0.0, 10); 
   Channel_Dilep = new TH1F(Form("Channel_Dilep" ), Form(" Num. leptons "), 10, 0.0, 10); 
   Channel_non = new TH1F(Form("Channel_non" ), Form(" Num. leptons "), 10, 0.0, 10); 

}

void ssbgen_analysis::End()
{
   fout->Write();
   fout->Close();
}

void ssbgen_analysis::SetOutputFileName(char *outname)
{
   outfile = outname;
}

