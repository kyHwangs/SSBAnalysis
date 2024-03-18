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
      ////////////////////////////////////////
      /// start Main Loop Function and Cut ///
      ////////////////////////////////////////
      TLorentzVector* ele;
      for(int i = 0; i < Elec->GetEntries(); ++i )
      {
         ele = (TLorentzVector*)Elec->At(i);
         //cout << "Electron pT " << ele->Pt() << endl;
         FillHisto(h_cf_elept, ele->Pt());
      }
      TLorentzVector* met;
      for(int i = 0; i < MET->GetEntries(); ++i )
      {
         met = (TLorentzVector*)MET->At(i);
         //cout << "METtron pT " << met->Pt() << endl;
         FillHisto(h_cf_metpt, met->Pt());
      }
  
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

   /// Test For Systematic All-in-One Code ///
   h_cf_elept          = new TH1D(Form("_h_cf_elept_"),Form("Electron pT"), 1000, 0, 1000); h_cf_elept->Sumw2();
   h_cf_metpt          = new TH1D(Form("_h_cf_metpt_"),Form("MET"), 1000, 0, 1000); h_cf_metpt->Sumw2();


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


