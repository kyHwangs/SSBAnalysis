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

void ssb_analysis::SetConfig()
{
 
   // Lep Infor ( Di Muon or Di Electron )
   lepisotype = SSBConfReader->GetText( "Iso_type" );
   lepId      = SSBConfReader->GetText( "Lep_ID"   );
   n_elep_Id  = SSBConfReader->Size( "Lep_ID_cut"  );
  
   // Muon Infor for Only MuEl channel
   Muonisotype = SSBConfReader->GetText( "MuonIso_type" );
   Muonid      = SSBConfReader->GetText( "Muon_ID"      );

   // Electron Infor Only for MuEl channel
   Elecisotype = SSBConfReader->GetText( "ElecIso_type" );
   Elecid      = SSBConfReader->GetText( "Elec_ID"      );
   n_elec_Id   = SSBConfReader->Size(    "Elec_ID_cut"  );
   // Jet Infor
   jetId      = SSBConfReader->GetText("Jet_ID");
   jetbtag    = SSBConfReader->GetText("Jet_btag");


   // Kinematic variables for Object //(not MuEl)
   lep_pt     = SSBConfReader->GetNumber( "Lep_pt"  );
   lep_eta    = SSBConfReader->GetNumber( "Lep_eta" );
   lepisocut  = SSBConfReader->GetNumber( "Iso_cut" );
   jet_pt     = SSBConfReader->GetNumber( "Jet_pt"  );
   jet_eta    = SSBConfReader->GetNumber( "Jet_eta" );
   met_cut    = SSBConfReader->GetNumber( "MET_cut" );

   // Kinematic variables for Object (MuEl)
   muon_pt     = SSBConfReader->GetNumber( "Muon_pt"     );
   muon_eta    = SSBConfReader->GetNumber( "Muon_eta"    );
   muonisocut  = SSBConfReader->GetNumber( "MuonIso_cut" );
   elec_pt     = SSBConfReader->GetNumber( "Elec_pt"     );
   elec_eta    = SSBConfReader->GetNumber( "Elec_eta"    );
   elecisocut  = SSBConfReader->GetNumber( "ElecIso_cut" );
   

   // variables for Jet Cleaning 
   muiso_type = SSBConfReader->GetText( "Muiso_for_Jetcleanning" );
   muId       = SSBConfReader->GetText( "Muid_for_Jetcleanning"  );
   eliso_type = SSBConfReader->GetText( "Eliso_for_Jetcleanning" );
   elId       = SSBConfReader->GetText( "Elid_for_Jetcleanning"  );

   // Kinematic variables for Jet Cleaning 
   mu_pt_jetcl   = SSBConfReader->GetNumber( "Mupt_for_Jetcleanning"     );
   mu_eta_jetcl  = SSBConfReader->GetNumber( "Mueta_for_Jetcleanning"    );
   mu_iso_jetcl  = SSBConfReader->GetNumber( "Muisocut_for_Jetcleanning" );
   el_pt_jetcl   = SSBConfReader->GetNumber( "Elpt_for_Jetcleanning"     );
   el_eta_jetcl  = SSBConfReader->GetNumber( "Eleta_for_Jetcleanning"    );
   el_iso_jetcl  = SSBConfReader->GetNumber( "Elisocut_for_Jetcleanning" );
   n_el_id_jetcl = SSBConfReader->Size(      "Elidcut_for_Jetcleanning"  );

   if (n_el_id_jetcl == 1)
   {
      el_id_jetcl_1 = SSBConfReader->GetNumber("Elidcut_for_Jetcleanning");
      el_id_jetcl_2 = -999;
   }
   else if (n_el_id_jetcl == 2 )
   {
      el_id_jetcl_1 = SSBConfReader->GetNumber("Elidcut_for_Jetcleanning", 1);
      el_id_jetcl_2 = SSBConfReader->GetNumber("Elidcut_for_Jetcleanning", 2);
   }
   else {cout << " Electron Jetcleaning error" <<endl; }


   /////////////////////////////
   // For KinematicSolver !!! //
   /////////////////////////////

//   double nupair[] = { 30.7137, 56.2880, 23.0744, 59.1015,24.9145 }; // for 8TeV
   double nupair[] = { 1/4.53731e+06,  5.81649e+01,  2.31911e+01,  1.19038e+04, 
                       -9.78134e+01 , -2.76894e+01, -2.01114e-02, -5.07004e+00,
                         5.76608e+01,  2.34419e+01,  1.18119e+04,  5.26711e+03,
                         1.64080e+03,  5.86524e+00, -8.07173e+01, -5.44618e+00 };
   int nupair_size = sizeof(nupair)/sizeof(nupair[0]);
   
   for (int i = 0; i < nupair_size; ++i ){ nupars_.push_back(nupair[i]); }

   topmass_begin = 50.0;
   topmass_end   = 500.0;
   topmass_step_ = 1.0;


}
void ssb_analysis::GetVariables()
{

   ///////////////////////
   /// Case of Di-Muon ///
   /////////////////////// 
   if ( TString(Decaymode).Contains( "dimuon" ) ) 
   { 

      //////////////////////
      /// Muon condition ///
      //////////////////////

      /// Muon Iso type
      if ( TString(lepisotype).Contains( "PFIsodbeta03" ) )       { v_lepton_iso = Muon_PFIsodBeta03; }
      else if ( TString(lepisotype).Contains( "PFIsodbeta04" ) )  { v_lepton_iso = Muon_PFIsodBeta04; }
      else { cout << "Muon Condition Iso Error" << endl; } 

      /// Muon ID
      if      (TString(lepId).Contains( "Loose"  ) ) { v_lepton_Id = Muon_isLoose;  }
      else if (TString(lepId).Contains( "Tight"  ) ) { v_lepton_Id = Muon_isTight;  }
      else { cout << "Muon ID Error" << endl; } 

   }
   ///////////////////////////
   /// Case of Di-Electron ///
   ///////////////////////////
   else if (TString(Decaymode).Contains( "dielec" ) )
   {
      /// Electron iso type
      if      ( TString( lepisotype ).Contains( "PFIsoRho03"   ) ) { v_lepton_iso = Elec_PFIsoRho03;   }
      else if ( TString( lepisotype ).Contains( "PFIsoRho04"   ) ) { v_lepton_iso = Elec_PFIsoRho04;   }
      else { cout << "Electron Iso type Error at Case of Di-Electron" << endl; }
   
      // Electron ID

      /// simple cut based electron ID
      if      (TString(lepId).Contains( "Loose"  ) ) { v_lepton_Id = Elec_SCB_Loose;  }
      else if (TString(lepId).Contains( "Medium" ) ) { v_lepton_Id = Elec_SCB_Medium; }
      else if (TString(lepId).Contains( "Tight"  ) ) { v_lepton_Id = Elec_SCB_Tight;  }
      else if (TString(lepId).Contains( "Veto"   ) ) { v_lepton_Id = Elec_SCB_Veto;   }
      else { cout << "Electron ID Error" << endl; }

   } 
   ////////////////////
   /// Case of MuEl ///
   ////////////////////
   else if (TString(Decaymode).Contains( "muel" ) )
   {
      /// Muon Iso type
      if      ( TString(Muonisotype).Contains( "PFIsodbeta03" ) )  { v_mlep_iso = Muon_PFIsodBeta03; }
      else if ( TString(Muonisotype).Contains( "PFIsodbeta04" ) )  { v_mlep_iso = Muon_PFIsodBeta04; }
      else { cout << "Muon Condition Iso Error at Case of MuEl" << endl; } 

      /// Muon ID
      if      ( TString(Muonid).Contains( "Loose"  ) ) { v_muon_Id = Muon_isLoose;  }
      else if ( TString(Muonid).Contains( "Tight"  ) ) { v_muon_Id = Muon_isTight;  }
      else { cout << "Muon Condition ID Error" << endl; } 

      /// Electron iso type
      if      ( TString( Elecisotype ).Contains( "PFIsoRho03" ) )   { v_elec_iso = Elec_PFIsoRho03;   }
      else if ( TString( Elecisotype ).Contains( "PFIsoRho04" ) )   { v_elec_iso = Elec_PFIsoRho04;   }
      else { cout << "Electron Iso type Error at Case of MuEl" << endl; }
   
      // Electron ID SCB
     
      if      (TString(Elecid).Contains( "Loose"  ) ) { v_electron_Id = Elec_SCB_Loose;  }
      else if (TString(Elecid).Contains( "Medium" ) ) { v_electron_Id = Elec_SCB_Medium; }
      else if (TString(Elecid).Contains( "Tight"  ) ) { v_electron_Id = Elec_SCB_Tight;  }
      else if (TString(Elecid).Contains( "Veto"   ) ) { v_electron_Id = Elec_SCB_Veto;   }
      else { cout << "Electron ID Error" << endl; } 

   }
   else { cout << "Error Decaymode !! " << endl; }

   /////////////////////////////
   /// *** Jet condition *** ///
   /////////////////////////////

   if      ( TString(jetId).Contains( "PFLoose" ) ) { v_jet_Id = Jet_PFId; jet_id = 1; }
   else if ( TString(jetId).Contains( "PFTight" ) ) { v_jet_Id = Jet_PFId; jet_id = 2; }
   else {cout << "Jet condition error" << endl;}

   if      ( TString(jetbtag).Contains( "CSVL"  ) )    { bdisccut = 0.244; }
   else if ( TString(jetbtag).Contains( "CSVM"  ) )    { bdisccut = 0.679; }
   else if ( TString(jetbtag).Contains( "CSVT"  ) )    { bdisccut = 0.898; }
   else if ( TString(jetbtag).Contains( "CISVL" ) )    { bdisccut = 0.423; }
   else if ( TString(jetbtag).Contains( "CISVM" ) )    { bdisccut = 0.814; }
   else if ( TString(jetbtag).Contains( "CISVT" ) )    { bdisccut = 0.941; }
   else if ( TString(jetbtag).Contains( "pfCSVV2L" ) ) { bdisccut = 0.5426; }
   else if ( TString(jetbtag).Contains( "pfCSVV2M" ) ) { bdisccut = 0.8484; }
   else if ( TString(jetbtag).Contains( "pfCSVV2T" ) ) { bdisccut = 0.9535; }
   else { cout << "bscriminator error !!" << endl; }
   //////////////////////////
   // *** JetCleanning ***///
   //////////////////////////

   ///////////////////////////
   // Muon for Jet Cleaning //
   ///////////////////////////

   /// Muon iso type
   if      ( TString( muiso_type ).Contains( "PFIsodbeta03" ) ) { v_mu_iso_jcl = Muon_PFIsodBeta03; }
   else if ( TString( muiso_type ).Contains( "PFIsodbeta04" ) ) { v_mu_iso_jcl = Muon_PFIsodBeta04; }
   else { cout << "Muon Iso type error for Jet cleanning" << endl;}

   /// Muon ID
   if      ( TString( muId ).Contains( "Loose" ) ) { v_mu_Id_jcl = Muon_isLoose; }
   else if ( TString( muId ).Contains( "Tight" ) ) { v_mu_Id_jcl = Muon_isTight; }
   else { cout << "Muon ID error for Jet cleanning" << endl;}
 
   ///////////////////////////////
   // Electron for Jet Cleaning //
   ///////////////////////////////

   /// Electron iso type
   if      ( TString( eliso_type ).Contains( "PFIsoRho03"   ) ) { v_el_iso_jcl = Elec_PFIsoRho03;   }
   else if ( TString( eliso_type ).Contains( "PFIsoRho04"   ) ) { v_el_iso_jcl = Elec_PFIsoRho04;   }
   else { cout << "Electron Iso type Error for Jet Cleaning" << endl;}

   // Electron ID SCB 

   if      (TString(elId).Contains( "Loose"  ) ) { v_el_Id_jcl = Elec_SCB_Loose;  }
   else if (TString(elId).Contains( "Medium" ) ) { v_el_Id_jcl = Elec_SCB_Medium; }
   else if (TString(elId).Contains( "Tight"  ) ) { v_el_Id_jcl = Elec_SCB_Tight;  }
   else if (TString(elId).Contains( "Veto"   ) ) { v_el_Id_jcl = Elec_SCB_Veto;   }
   else { cout << "Electron ID Error" << endl; }

}

void ssb_analysis::TLVInitial()
{
      Lep1      = new TLorentzVector(); // Leading Lepton 
      Lep2      = new TLorentzVector(); // Second Leading Lepton

      TMuon     = new TLorentzVector(); // Leading Lepton 
      TElectron = new TLorentzVector(); // Second Leading Lepton

      Jet1      = new TLorentzVector(); // Leading Jet 
      Jet2      = new TLorentzVector(); // Second Leading Jet

      Jet1Up      = new TLorentzVector(); // Leading Jet 
      Jet2Up      = new TLorentzVector(); // Second Leading Jet

      Jet1Dn      = new TLorentzVector(); // Leading Jet 
      Jet2Dn      = new TLorentzVector(); // Second Leading Jet

      bJet1     = new TLorentzVector(); // Leading Jet 
      bJet2     = new TLorentzVector(); // Second Leading Jet
      bJet1Up     = new TLorentzVector(); // Leading Jet 
      bJet2Up     = new TLorentzVector(); // Second Leading Jet
      bJet1Dn     = new TLorentzVector(); // Leading Jet 
      bJet2Dn     = new TLorentzVector(); // Second Leading Jet

      Met       = new TLorentzVector(); // MET
      MetJESUp  = new TLorentzVector(); // MET
      MetJESDn  = new TLorentzVector(); // MET
  
      Nu1       = new TLorentzVector(); // Nutrino 1 ( For Lep1 )
      Nu2       = new TLorentzVector(); // Nutrino 2 ( For Lep2 )

      W1        = new TLorentzVector(); // W boson 1 ( For Lep1 )  
      W2        = new TLorentzVector(); // W boson 2 ( For Lep2 )

      Top       = new TLorentzVector(); // Top
      AnTop     = new TLorentzVector(); // Anti-top

      Top1      = new TLorentzVector(); // Top Leading Pt
      Top2      = new TLorentzVector(); // Top Second Leading Pt 

      bJet      = new TLorentzVector(); // bJet 
      AnbJet    = new TLorentzVector(); // Anti-bJet

      Lep       = new TLorentzVector(); // Lepton 
      AnLep     = new TLorentzVector(); // Anti-Lepton

      Nu       = new TLorentzVector(); // Lepton 
      AnNu     = new TLorentzVector(); // Anti-Lepton

      // Lepton selection
      Mu_lep_sel  = new TLorentzVector();
      Ele_lep_sel = new TLorentzVector();

      // Veto muon and electron
      vetoMu  = new TLorentzVector();
      vetoEle = new TLorentzVector();

      // Jet
      Mu_jet_clean           = new TLorentzVector();     
      Ele_jet_clean          = new TLorentzVector();
      selected_Ele_jet_clean = new TLorentzVector();
      selected_Mu_jet_clean  = new TLorentzVector();
      TJet                   = new TLorentzVector();
}


void ssb_analysis::Loop( char *logfile )
{

   GetTotalEvent();
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
   cout << "Ntuple Total Event Check !! " << NtupletotalEvent << endl;

   MCSF();
   cout << "MC_SF Check !! " << mc_sf_ << endl;

   Np_eventw_2 = 0.0;
   Nm_eventw_2 = 0.0;
   Num_cpo3_p  = 0; 
   Num_cpo3_m  = 0; 

   int count_emu = 0;
  
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

      // Make TL to use SSBMiniTree //
      MakeVecforTL();

      // initailizing TLorentzVector
      TLVInitial();

      evt_weight_ = 1;

      GetVariables();
      // Apply MC Scale Factor //
      MCSFApply();
      // Apply GenWeight //
      GenWeightApply();
      // Apply PileUpReWeight

      PDFWeightApply();
      PileUpReWeightApply();
      RMDuplEvt();
      //TopPtReweightApply(); /// Apply Top pT Reweight ///
      ////////////////////////////////////////
      /// ** Di-Lepton Channel Analysis ** ///
      ////////////////////////////////////////
    
      if ( TString(Decaymode).Contains( "dielec" ) ||
           TString(Decaymode).Contains( "dimuon" ) || 
           TString(Decaymode).Contains( "muel" )     )
      {
//         if ( EveRun() == false ) { continue;}
      
         /////////////////////////////////////////////
         /// Finding Di-Lep. Channel at Gen.Level. ///
         /////////////////////////////////////////////

         if ( ChannelIndex() == false) {continue;}
         // Num. primary vertex counter
         NumPVCount();
         // Lepton Define //
         LeptonSelector();
         LeptonOrder();
         // Jet Define //
         JetSelector();
         // Met Define //
         METDefiner();
         BDsicApply();
         BJetDefiner(); 

         /// To check up the vertex distribution before pre-selection ///
         FillHisto(h_Num_PV_BeforePreSel, num_pv, evt_weight_);

         if (METFilterAPP() == true){
            FillHisto(h_Num_PV_AfterMetFilter, num_pv, evt_weight_);
         }
         if (Trigger() == true){
            FillHisto(h_Num_PV_AfterTrigger, num_pv, evt_weight_);
         }
         ////////////////////
         /// Event Filter ///
         ////////////////////
         if ( METFilterAPP() == false ) {continue;}
         if ( Filter_PV->at(0) == false ) {continue;}

         //////////////////////////////////
         /// trigger requirement step 0 ///
         //////////////////////////////////
         if ( Trigger() == false ) {continue;}
         // Apply Trigger SF //

         FillHisto( h_EventWeight[0] , evt_weight_  );
         FillHisto( h_cf_NLeptons[0], Muon_Count, evt_weight_);
         FillHisto( h_cf_NJets[0], v_jet_idx.size(), evt_weight_ );
         FillHisto( h_cf_NPV[0]    , num_pv, evt_weight_ );
         FillHisto( h_Num_PV[0]    , num_pv, evt_weight_ );
      
         if ( NumIsoLeptons() == false ){continue;}
         TriggerSFApply();
         ///////////////////////
         /// 3rd Lepton Veto ///
         ///////////////////////
         if ( ThirdLeptonVeto() == false ){continue;}
         //if (LeptonsPtAddtional() == false ) {continue;}
         // TLorentzVector define

         ////////////////////////////////
         /// DiLeptonMassCut() step 1 ///
         ////////////////////////////////
         if ( DiLeptonMassCut() == false ) {continue;}
         LeptonSFApply();

         FillHisto( h_EventWeight[1] , evt_weight_  );
         FillHisto( h_cf_NPV[1]    , num_pv        , evt_weight_ );
         FillHisto( h_cf_NLeptons[1], v_lepton_idx.size(), evt_weight_ );
         FillHisto( h_cf_Lep1pt[1] , (*Lep1).Pt()  , evt_weight_ );
         FillHisto( h_cf_Lep1eta[1], (*Lep1).Eta() , evt_weight_ );
         FillHisto( h_cf_Lep1phi[1], (*Lep1).Phi() , evt_weight_ );
         FillHisto( h_cf_Lep2pt[1] , (*Lep2).Pt()  , evt_weight_ );
         FillHisto( h_cf_Lep2eta[1], (*Lep2).Eta() , evt_weight_ );
         FillHisto( h_cf_Lep2phi[1], (*Lep2).Phi() , evt_weight_ );
         FillHisto( h_cf_dilep_inv_mass[1], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
         FillHisto( h_cf_NJets[1]  , v_jet_idx.size(), evt_weight_ );
         if(v_jet_idx.size() > 0 )
         {
            FillHisto( h_cf_Jet1pt[1] , (*Jet1).Pt()  , evt_weight_ );
            FillHisto( h_cf_Jet1eta[1], (*Jet1).Eta() , evt_weight_ );
            FillHisto( h_cf_Jet1phi[1], (*Jet1).Phi() , evt_weight_ );
            if (v_jet_idx.size() > 1)
            {
               FillHisto( h_cf_Jet2pt[1] , (*Jet2).Pt()  , evt_weight_ );
               FillHisto( h_cf_Jet2eta[1], (*Jet2).Eta() , evt_weight_ );
               FillHisto( h_cf_Jet2phi[1], (*Jet2).Phi() , evt_weight_ );
            }
         }
         FillHisto( h_cf_metpt[1] ,Met->Pt()  , evt_weight_ );
         FillHisto( h_cf_metphi[1],Met->Phi() , evt_weight_ );
         FillHisto( h_DiLepMass[1], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
         FillHisto( h_Num_PV[1], num_pv, evt_weight_ );
         FillHisto( h_Lep1pt[1] , (*Lep1).Pt()  , evt_weight_ );
         FillHisto( h_Lep1eta[1], (*Lep1).Eta() , evt_weight_ );
         FillHisto( h_Lep1phi[1], (*Lep1).Phi() , evt_weight_ );
         FillHisto( h_Lep2pt[1] , (*Lep2).Pt()  , evt_weight_ );
         FillHisto( h_Lep2eta[1], (*Lep2).Eta() , evt_weight_ );
         FillHisto( h_Lep2phi[1], (*Lep2).Phi() , evt_weight_ );
         FillHisto( h_METpt[1]   , Met->Pt()  , evt_weight_ );
         FillHisto( h_METphi[1]  , Met->Phi()  , evt_weight_ );
         if(v_jet_idx.size() > 0 )
         {
            FillHisto( h_Jet1pt[1] , (*Jet1).Pt()  , evt_weight_ );
            FillHisto( h_Jet1eta[1], (*Jet1).Eta() , evt_weight_ );
            FillHisto( h_Jet1phi[1], (*Jet1).Phi() , evt_weight_ );
            if (v_jet_idx.size() > 1)
            {
               FillHisto( h_Jet2pt[1] , (*Jet2).Pt()  , evt_weight_ );
               FillHisto( h_Jet2eta[1], (*Jet2).Eta() , evt_weight_ );
               FillHisto( h_Jet2phi[1], (*Jet2).Phi() , evt_weight_ );
            }
         }
         if (isAllSyst == true){
            for ( int i = 0; i < v_SystFullName.size(); ++i )
            {  
               if (TString(v_SystFullName[i]).Contains("JetEn")){continue;}
               FillHisto( h_cf_sys_NLeptons[i][1], v_lepton_idx.size(), v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep1pt[i][1] , (*Lep1).Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep1eta[i][1], (*Lep1).Eta() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep1phi[i][1], (*Lep1).Phi() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep2pt[i][1] , (*Lep2).Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep2eta[i][1], (*Lep2).Eta() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep2phi[i][1], (*Lep2).Phi() , v_SystEvt[i] );
               FillHisto( h_cf_sys_NPV[i][1]    , num_pv        , v_SystEvt[i] );
               FillHisto( h_cf_sys_NJets[i][1]  , v_jet_idx.size(), v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet1pt[i][1] , (*Jet1).Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet1eta[i][1], (*Jet1).Eta() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet1phi[i][1], (*Jet1).Phi() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet2pt[i][1] , (*Jet2).Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet2eta[i][1], (*Jet2).Eta() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet2phi[i][1], (*Jet2).Phi() , v_SystEvt[i] );
               FillHisto( h_cf_sys_metpt[i][1] ,Met->Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_metphi[i][1],Met->Phi() , v_SystEvt[i] );
               
               FillHisto( h_cf_sys_dilep_inv_mass[i][1], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
               FillHisto( h_cf_sys_metpt[i][1] , Met->Pt() , v_SystEvt[i] );
               FillHisto( h_cf_sys_metphi[i][1], Met->Phi(), v_SystEvt[i] );
               
               FillHisto( h_sys_Num_PV[i][1], num_pv, v_SystEvt[i] );
               FillHisto( h_sys_DiLepMass[i][1], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
               FillHisto( h_sys_Lep1pt[i][1] , (*Lep1).Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_Lep1eta[i][1], (*Lep1).Eta() , v_SystEvt[i] );
               FillHisto( h_sys_Lep1phi[i][1], (*Lep1).Phi() , v_SystEvt[i] );
               FillHisto( h_sys_Lep2pt[i][1] , (*Lep2).Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_Lep2eta[i][1], (*Lep2).Eta() , v_SystEvt[i] );
               FillHisto( h_sys_Lep2phi[i][1], (*Lep2).Phi() , v_SystEvt[i] );
               FillHisto( h_sys_Jet1pt[i][1] , (*Jet1).Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_Jet1eta[i][1], (*Jet1).Eta() , v_SystEvt[i] );
               FillHisto( h_sys_Jet1phi[i][1], (*Jet1).Phi() , v_SystEvt[i] );
               FillHisto( h_sys_Jet2pt[i][1] , (*Jet2).Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_Jet2eta[i][1], (*Jet2).Eta() , v_SystEvt[i] );
               FillHisto( h_sys_Jet2phi[i][1], (*Jet2).Phi() , v_SystEvt[i] );
               FillHisto( h_sys_Num_Jets[i][1], v_jet_idx.size(), v_SystEvt[i] );
               FillHisto( h_sys_METpt[i][1]   , Met->Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_METphi[i][1]  , Met->Phi()  , v_SystEvt[i] );
            }     
         }
      
         // ZVetocut step 2 
         if ( ZVetoCut() == false ) {continue;}
         FillHisto( h_EventWeight[2], evt_weight_  );
         FillHisto( h_cf_NLeptons[2], v_lepton_idx.size(), evt_weight_ );
         FillHisto( h_cf_Lep1pt[2] , (*Lep1).Pt()  , evt_weight_ );
         FillHisto( h_cf_Lep1eta[2], (*Lep1).Eta() , evt_weight_ );
         FillHisto( h_cf_Lep1phi[2], (*Lep1).Phi() , evt_weight_ );
         FillHisto( h_cf_Lep2pt[2] , (*Lep2).Pt()  , evt_weight_ );
         FillHisto( h_cf_Lep2eta[2], (*Lep2).Eta() , evt_weight_ );
         FillHisto( h_cf_Lep2phi[2], (*Lep2).Phi() , evt_weight_ );
         FillHisto( h_cf_NPV[2]    , num_pv        , evt_weight_ );
         FillHisto( h_cf_NJets[2]  , v_jet_idx.size(), evt_weight_ );
         FillHisto( h_cf_dilep_inv_mass[2], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
         FillHisto( h_cf_metpt[2] , Met->Pt() , evt_weight_ );
         FillHisto( h_cf_metphi[2], Met->Phi(), evt_weight_ );
         if(v_jet_idx.size() > 0 )
         {
            FillHisto( h_cf_Jet1pt[2] , (*Jet1).Pt()  , evt_weight_ );
            FillHisto( h_cf_Jet1eta[2], (*Jet1).Eta() , evt_weight_ );
            FillHisto( h_cf_Jet1phi[2], (*Jet1).Phi() , evt_weight_ );
            if (v_jet_idx.size() > 1)
            {
               FillHisto( h_cf_Jet2pt[2] , (*Jet2).Pt()  , evt_weight_ );
               FillHisto( h_cf_Jet2eta[2], (*Jet2).Eta() , evt_weight_ );
               FillHisto( h_cf_Jet2phi[2], (*Jet2).Phi() , evt_weight_ );
            }
         }

         FillHisto( h_Num_PV[2], num_pv, evt_weight_ );
         FillHisto( h_DiLepMass[2], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
         FillHisto( h_Lep1pt[2] , (*Lep1).Pt()  , evt_weight_ );
         FillHisto( h_Lep1eta[2], (*Lep1).Eta() , evt_weight_ );
         FillHisto( h_Lep1phi[2], (*Lep1).Phi() , evt_weight_ );
         FillHisto( h_Lep2pt[2] , (*Lep2).Pt()  , evt_weight_ );
         FillHisto( h_Lep2eta[2], (*Lep2).Eta() , evt_weight_ );
         FillHisto( h_Lep2phi[2], (*Lep2).Phi() , evt_weight_ );
         FillHisto( h_Num_Jets[2]  , v_jet_idx.size(), evt_weight_ );
         FillHisto( h_METpt[2]   , Met->Pt()  , evt_weight_ );
         FillHisto( h_METphi[2]  , Met->Phi()  , evt_weight_ );
         if(v_jet_idx.size() > 0 )
         {
            FillHisto( h_Jet1pt[2] , (*Jet1).Pt()  , evt_weight_ );
            FillHisto( h_Jet1eta[2], (*Jet1).Eta() , evt_weight_ );
            FillHisto( h_Jet1phi[2], (*Jet1).Phi() , evt_weight_ );
            if (v_jet_idx.size() > 1)
            {
               FillHisto( h_Jet2pt[2] , (*Jet2).Pt()  , evt_weight_ );
               FillHisto( h_Jet2eta[2], (*Jet2).Eta() , evt_weight_ );
               FillHisto( h_Jet2phi[2], (*Jet2).Phi() , evt_weight_ );
            }
         }
         //FillHisto( h_cf_met[2], , evt_weight_ );
         //FillHisto( h_cf_Nbjets[2], , evt_weight_ );
         if (isAllSyst == true){
            for ( int i = 0; i < v_SystFullName.size(); ++i )
            {  
            //   if (TString(v_SystFullName[i]).Contains("JetEn")){continue;}
               FillHisto( h_cf_sys_NLeptons[i][2], v_lepton_idx.size(), v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep1pt[i][2] , (*Lep1).Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep1eta[i][2], (*Lep1).Eta() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep1phi[i][2], (*Lep1).Phi() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep2pt[i][2] , (*Lep2).Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep2eta[i][2], (*Lep2).Eta() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Lep2phi[i][2], (*Lep2).Phi() , v_SystEvt[i] );
               FillHisto( h_cf_sys_NPV[i][2]    , num_pv        , v_SystEvt[i] );
               FillHisto( h_cf_sys_NJets[i][2]  , v_jet_idx.size(), v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet1pt[i][2] , (*Jet1).Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet1eta[i][2], (*Jet1).Eta() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet1phi[i][2], (*Jet1).Phi() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet2pt[i][2] , (*Jet2).Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet2eta[i][2], (*Jet2).Eta() , v_SystEvt[i] );
               FillHisto( h_cf_sys_Jet2phi[i][2], (*Jet2).Phi() , v_SystEvt[i] );
               FillHisto( h_cf_sys_metpt[i][2] ,Met->Pt()  , v_SystEvt[i] );
               FillHisto( h_cf_sys_metphi[i][2],Met->Phi() , v_SystEvt[i] );
               
               FillHisto( h_cf_sys_dilep_inv_mass[i][2], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
               FillHisto( h_cf_sys_metpt[i][2] , Met->Pt() , v_SystEvt[i] );
               FillHisto( h_cf_sys_metphi[i][2], Met->Phi(), v_SystEvt[i] );
               
               FillHisto( h_sys_Num_PV[i][2], num_pv, v_SystEvt[i] );
               FillHisto( h_sys_DiLepMass[i][2], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
               FillHisto( h_sys_Lep1pt[i][2] , (*Lep1).Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_Lep1eta[i][2], (*Lep1).Eta() , v_SystEvt[i] );
               FillHisto( h_sys_Lep1phi[i][2], (*Lep1).Phi() , v_SystEvt[i] );
               FillHisto( h_sys_Lep2pt[i][2] , (*Lep2).Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_Lep2eta[i][2], (*Lep2).Eta() , v_SystEvt[i] );
               FillHisto( h_sys_Lep2phi[i][2], (*Lep2).Phi() , v_SystEvt[i] );
               FillHisto( h_sys_Jet1pt[i][2] , (*Jet1).Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_Jet1eta[i][2], (*Jet1).Eta() , v_SystEvt[i] );
               FillHisto( h_sys_Jet1phi[i][2], (*Jet1).Phi() , v_SystEvt[i] );
               FillHisto( h_sys_Jet2pt[i][2] , (*Jet2).Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_Jet2eta[i][2], (*Jet2).Eta() , v_SystEvt[i] );
               FillHisto( h_sys_Jet2phi[i][2], (*Jet2).Phi() , v_SystEvt[i] );
               FillHisto( h_sys_Num_Jets[i][2], v_jet_idx.size(), v_SystEvt[i] );
               FillHisto( h_sys_METpt[i][2]   , Met->Pt()  , v_SystEvt[i] );
               FillHisto( h_sys_METphi[i][2]  , Met->Phi()  , v_SystEvt[i] );
            }
         }
         /////////////////////////      
         // Num. Jet cut step 3 //
         /////////////////////////
         JetDefiner();
         //if ( NumJetCut(v_jet_idx) == false ) {continue;}
         if ( NumJetCut(v_jet_idx) == true ) 
         {

            // Define Leading Jet and Second Leading Jet 
            FillHisto( h_EventWeight[3] , evt_weight_  );
            FillHisto( h_cf_NLeptons[3], v_lepton_idx.size(), evt_weight_ );
            FillHisto( h_cf_Lep1pt[3] , (*Lep1).Pt()  , evt_weight_ );
            FillHisto( h_cf_Lep1eta[3], (*Lep1).Eta() , evt_weight_ );
            FillHisto( h_cf_Lep1phi[3], (*Lep1).Phi() , evt_weight_ );
            FillHisto( h_cf_Lep2pt[3] , (*Lep2).Pt()  , evt_weight_ );
            FillHisto( h_cf_Lep2eta[3], (*Lep2).Eta() , evt_weight_ );
            FillHisto( h_cf_Lep2phi[3], (*Lep2).Phi() , evt_weight_ );
            FillHisto( h_cf_NPV[3]    , num_pv        , evt_weight_ );
            FillHisto( h_cf_NJets[3]  , v_jet_idx.size(), evt_weight_ );
            FillHisto( h_cf_Jet1pt[3] , (*Jet1).Pt()  , evt_weight_ );
            FillHisto( h_cf_Jet1eta[3], (*Jet1).Eta() , evt_weight_ );
            FillHisto( h_cf_Jet1phi[3], (*Jet1).Phi() , evt_weight_ );
            FillHisto( h_cf_Jet2pt[3] , (*Jet2).Pt()  , evt_weight_ );
            FillHisto( h_cf_Jet2eta[3], (*Jet2).Eta() , evt_weight_ );
            FillHisto( h_cf_Jet2phi[3], (*Jet2).Phi() , evt_weight_ );
            FillHisto( h_cf_metpt[3] ,Met->Pt()  , evt_weight_ );
            FillHisto( h_cf_metphi[3],Met->Phi() , evt_weight_ );
            
            FillHisto( h_cf_dilep_inv_mass[3], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
            FillHisto( h_cf_metpt[3] , Met->Pt() , evt_weight_ );
            FillHisto( h_cf_metphi[3], Met->Phi(), evt_weight_ );
            
            FillHisto( h_Num_PV[3], num_pv, evt_weight_ );
            FillHisto( h_DiLepMass[3], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
            FillHisto( h_Lep1pt[3] , (*Lep1).Pt()  , evt_weight_ );
            FillHisto( h_Lep1eta[3], (*Lep1).Eta() , evt_weight_ );
            FillHisto( h_Lep1phi[3], (*Lep1).Phi() , evt_weight_ );
            FillHisto( h_Lep2pt[3] , (*Lep2).Pt()  , evt_weight_ );
            FillHisto( h_Lep2eta[3], (*Lep2).Eta() , evt_weight_ );
            FillHisto( h_Lep2phi[3], (*Lep2).Phi() , evt_weight_ );
            FillHisto( h_Jet1pt[3] , (*Jet1).Pt()  , evt_weight_ );
            FillHisto( h_Jet1eta[3], (*Jet1).Eta() , evt_weight_ );
            FillHisto( h_Jet1phi[3], (*Jet1).Phi() , evt_weight_ );
            FillHisto( h_Jet2pt[3] , (*Jet2).Pt()  , evt_weight_ );
            FillHisto( h_Jet2eta[3], (*Jet2).Eta() , evt_weight_ );
            FillHisto( h_Jet2phi[3], (*Jet2).Phi() , evt_weight_ );
            FillHisto( h_Num_Jets[3], v_jet_idx.size(), evt_weight_ );
            FillHisto( h_METpt[3]   , Met->Pt()  , evt_weight_ );
            FillHisto( h_METphi[3]  , Met->Phi()  , evt_weight_ );

            if (isAllSyst == true){
               for ( int i = 0; i < v_SystFullName.size(); ++i )
               {  
                  if (TString(v_SystFullName[i]).Contains("JetEn")){continue;}
                  FillHisto( h_cf_sys_NLeptons[i][3], v_lepton_idx.size(), v_SystEvt[i] );
                  FillHisto( h_cf_sys_Lep1pt[i][3] , (*Lep1).Pt()  , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Lep1eta[i][3], (*Lep1).Eta() , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Lep1phi[i][3], (*Lep1).Phi() , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Lep2pt[i][3] , (*Lep2).Pt()  , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Lep2eta[i][3], (*Lep2).Eta() , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Lep2phi[i][3], (*Lep2).Phi() , v_SystEvt[i] );
                  FillHisto( h_cf_sys_NPV[i][3]    , num_pv        , v_SystEvt[i] );
                  FillHisto( h_cf_sys_NJets[i][3]  , v_jet_idx.size(), v_SystEvt[i] );
                  FillHisto( h_cf_sys_Jet1pt[i][3] , (*Jet1).Pt()  , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Jet1eta[i][3], (*Jet1).Eta() , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Jet1phi[i][3], (*Jet1).Phi() , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Jet2pt[i][3] , (*Jet2).Pt()  , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Jet2eta[i][3], (*Jet2).Eta() , v_SystEvt[i] );
                  FillHisto( h_cf_sys_Jet2phi[i][3], (*Jet2).Phi() , v_SystEvt[i] );
                  FillHisto( h_cf_sys_metpt[i][3] ,Met->Pt()  , v_SystEvt[i] );
                  FillHisto( h_cf_sys_metphi[i][3],Met->Phi() , v_SystEvt[i] );
                  
                  FillHisto( h_cf_sys_dilep_inv_mass[i][3], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                  FillHisto( h_cf_sys_metpt[i][3] , Met->Pt() , v_SystEvt[i] );
                  FillHisto( h_cf_sys_metphi[i][3], Met->Phi(), v_SystEvt[i] );
                  
                  FillHisto( h_sys_Num_PV[i][3], num_pv, v_SystEvt[i] );
                  FillHisto( h_sys_DiLepMass[i][3], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                  FillHisto( h_sys_Lep1pt[i][3] , (*Lep1).Pt()  , v_SystEvt[i] );
                  FillHisto( h_sys_Lep1eta[i][3], (*Lep1).Eta() , v_SystEvt[i] );
                  FillHisto( h_sys_Lep1phi[i][3], (*Lep1).Phi() , v_SystEvt[i] );
                  FillHisto( h_sys_Lep2pt[i][3] , (*Lep2).Pt()  , v_SystEvt[i] );
                  FillHisto( h_sys_Lep2eta[i][3], (*Lep2).Eta() , v_SystEvt[i] );
                  FillHisto( h_sys_Lep2phi[i][3], (*Lep2).Phi() , v_SystEvt[i] );
                  FillHisto( h_sys_Jet1pt[i][3] , (*Jet1).Pt()  , v_SystEvt[i] );
                  FillHisto( h_sys_Jet1eta[i][3], (*Jet1).Eta() , v_SystEvt[i] );
                  FillHisto( h_sys_Jet1phi[i][3], (*Jet1).Phi() , v_SystEvt[i] );
                  FillHisto( h_sys_Jet2pt[i][3] , (*Jet2).Pt()  , v_SystEvt[i] );
                  FillHisto( h_sys_Jet2eta[i][3], (*Jet2).Eta() , v_SystEvt[i] );
                  FillHisto( h_sys_Jet2phi[i][3], (*Jet2).Phi() , v_SystEvt[i] );
                  FillHisto( h_sys_Num_Jets[i][3], v_jet_idx.size(), v_SystEvt[i] );
                  FillHisto( h_sys_METpt[i][3]   , Met->Pt()  , v_SystEvt[i] );
                  FillHisto( h_sys_METphi[i][3]  , Met->Phi()  , v_SystEvt[i] );
               }
            }
            //////////////////// 
            // MET cut step 4 //
            ////////////////////
            //if ( METCut(Met) == false ) {continue;}
            if ( METCut(Met) == true ) {            
               FillHisto( h_EventWeight[4] , evt_weight_  );
               FillHisto( h_cf_NLeptons[4], v_lepton_idx.size(), evt_weight_ );
               
               FillHisto( h_cf_Lep1pt[4] , (*Lep1).Pt()  , evt_weight_ );
               FillHisto( h_cf_Lep1eta[4], (*Lep1).Eta() , evt_weight_ );
               FillHisto( h_cf_Lep1phi[4], (*Lep1).Phi() , evt_weight_ );
               FillHisto( h_cf_Lep2pt[4] , (*Lep2).Pt()  , evt_weight_ );
               FillHisto( h_cf_Lep2eta[4], (*Lep2).Eta() , evt_weight_ );
               FillHisto( h_cf_Lep2phi[4], (*Lep2).Phi() , evt_weight_ );
               FillHisto( h_cf_NPV[4]    , num_pv        , evt_weight_ );
               FillHisto( h_cf_NJets[4]  , v_jet_idx.size(), evt_weight_ );
               FillHisto( h_cf_Jet1pt[4] , (*Jet1).Pt()  , evt_weight_ );
               FillHisto( h_cf_Jet1eta[4], (*Jet1).Eta() , evt_weight_ );
               FillHisto( h_cf_Jet1phi[4], (*Jet1).Phi() , evt_weight_ );
               FillHisto( h_cf_Jet2pt[4] , (*Jet2).Pt()  , evt_weight_ );
               FillHisto( h_cf_Jet2eta[4], (*Jet2).Eta() , evt_weight_ );
               FillHisto( h_cf_Jet2phi[4], (*Jet2).Phi() , evt_weight_ );
               
               FillHisto( h_cf_dilep_inv_mass[4], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
               FillHisto( h_cf_metpt[4] , Met->Pt() , evt_weight_ );
               FillHisto( h_cf_metphi[4], Met->Phi(), evt_weight_ );
               
               FillHisto( h_Lep1pt[4]  , Lep1->Pt() , evt_weight_ );
               FillHisto( h_Lep2pt[4]  , Lep2->Pt() , evt_weight_ );
               FillHisto( h_Lep1eta[4] , Lep1->Eta(), evt_weight_ );
               FillHisto( h_Lep2eta[4] , Lep2->Eta(), evt_weight_ );
               FillHisto( h_Lep1phi[4] , Lep1->Phi(), evt_weight_ );
               FillHisto( h_Lep2phi[4] , Lep2->Phi(), evt_weight_ );
               
               FillHisto( h_Jet1pt[4]  , Jet1->Pt() , evt_weight_ );
               FillHisto( h_Jet2pt[4]  , Jet2->Pt() , evt_weight_ );
               FillHisto( h_Jet1eta[4] , Jet1->Eta(), evt_weight_ );
               FillHisto( h_Jet2eta[4] , Jet2->Eta(), evt_weight_ );
               FillHisto( h_Jet1phi[4] , Jet1->Phi(), evt_weight_ );
               FillHisto( h_Jet2phi[4] , Jet2->Phi(), evt_weight_ );
               FillHisto( h_METpt[4]   , Met->Pt()  , evt_weight_ );
               FillHisto( h_METphi[4]  , Met->Phi() , evt_weight_ );
               
               FillHisto( h_DiLepMass[4], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
               FillHisto( h_Num_PV[4], num_pv, evt_weight_ );
               FillHisto( h_Num_Jets[4], v_jet_idx.size(), evt_weight_ );
               if (isAllSyst == true)
               {
                  for ( int i = 0; i < v_SystFullName.size(); ++i )
                  {  
                     if (TString(v_SystFullName[i]).Contains("JetEn")){continue;}
                     FillHisto( h_cf_sys_NLeptons[i][4], v_lepton_idx.size(), v_SystEvt[i] );
                     FillHisto( h_cf_sys_Lep1pt[i][4] , (*Lep1).Pt()  , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Lep1eta[i][4], (*Lep1).Eta() , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Lep1phi[i][4], (*Lep1).Phi() , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Lep2pt[i][4] , (*Lep2).Pt()  , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Lep2eta[i][4], (*Lep2).Eta() , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Lep2phi[i][4], (*Lep2).Phi() , v_SystEvt[i] );
                     FillHisto( h_cf_sys_NPV[i][4]    , num_pv        , v_SystEvt[i] );
                     FillHisto( h_cf_sys_NJets[i][4]  , v_jet_idx.size(), v_SystEvt[i] );
                     FillHisto( h_cf_sys_Jet1pt[i][4] , (*Jet1).Pt()  , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Jet1eta[i][4], (*Jet1).Eta() , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Jet1phi[i][4], (*Jet1).Phi() , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Jet2pt[i][4] , (*Jet2).Pt()  , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Jet2eta[i][4], (*Jet2).Eta() , v_SystEvt[i] );
                     FillHisto( h_cf_sys_Jet2phi[i][4], (*Jet2).Phi() , v_SystEvt[i] );
                     FillHisto( h_cf_sys_metpt[i][4] ,Met->Pt()  , v_SystEvt[i] );
                     FillHisto( h_cf_sys_metphi[i][4],Met->Phi() , v_SystEvt[i] );
                     
                     FillHisto( h_cf_sys_dilep_inv_mass[i][4], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                     FillHisto( h_cf_sys_metpt[i][4] , Met->Pt() , v_SystEvt[i] );
                     FillHisto( h_cf_sys_metphi[i][4], Met->Phi(), v_SystEvt[i] );
                     
                     FillHisto( h_sys_Num_PV[i][4], num_pv, v_SystEvt[i] );
                     FillHisto( h_sys_DiLepMass[i][4], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                     FillHisto( h_sys_Lep1pt[i][4] , (*Lep1).Pt()  , v_SystEvt[i] );
                     FillHisto( h_sys_Lep1eta[i][4], (*Lep1).Eta() , v_SystEvt[i] );
                     FillHisto( h_sys_Lep1phi[i][4], (*Lep1).Phi() , v_SystEvt[i] );
                     FillHisto( h_sys_Lep2pt[i][4] , (*Lep2).Pt()  , v_SystEvt[i] );
                     FillHisto( h_sys_Lep2eta[i][4], (*Lep2).Eta() , v_SystEvt[i] );
                     FillHisto( h_sys_Lep2phi[i][4], (*Lep2).Phi() , v_SystEvt[i] );
                     FillHisto( h_sys_Jet1pt[i][4] , (*Jet1).Pt()  , v_SystEvt[i] );
                     FillHisto( h_sys_Jet1eta[i][4], (*Jet1).Eta() , v_SystEvt[i] );
                     FillHisto( h_sys_Jet1phi[i][4], (*Jet1).Phi() , v_SystEvt[i] );
                     FillHisto( h_sys_Jet2pt[i][4] , (*Jet2).Pt()  , v_SystEvt[i] );
                     FillHisto( h_sys_Jet2eta[i][4], (*Jet2).Eta() , v_SystEvt[i] );
                     FillHisto( h_sys_Jet2phi[i][4], (*Jet2).Phi() , v_SystEvt[i] );
                     FillHisto( h_sys_Num_Jets[i][4], v_jet_idx.size(), v_SystEvt[i] );
                     FillHisto( h_sys_METpt[i][4]   , Met->Pt()  , v_SystEvt[i] );
                     FillHisto( h_sys_METphi[i][4]  , Met->Phi()  , v_SystEvt[i] );
                  }
               }
               double AllJetpt = 0;
               for ( int ijet = 0; ijet < v_jet_idx.size(); ++ijet )
               { 
//                  TLorentzVector *htJet = (TLorentzVector*)Jet->At( v_jet_idx[ijet] );
                  TLorentzVector *htJet = v_jet_TL[ijet];
                  AllJetpt += htJet->Pt(); 
               }
               
               FillHisto( h_HT[4], AllJetpt, evt_weight_);
               
               /////////////////////////////////////////
               /// One or more b-Tagging Requirement ///
               /////////////////////////////////////////
               
               BTaggigSFApply();
               //if ( BJetCut(v_bjet_idx) == false ) {continue;} // one or more b-tagging 
               if ( BJetCut(v_bjet_idx) == true ) // one or more b-tagging 
               {  
               //   cout << " btag sf : "<< evt_weight_/evt_weight_beforeBtag_<< endl;
                  FillHisto(h_bTagWeight,evt_weight_/evt_weight_beforeBtag_); 
                  FillHisto( h_EventWeight[5] , evt_weight_  );
                  FillHisto( h_cf_NLeptons[5], v_lepton_idx.size(), evt_weight_ );
                  
                  FillHisto( h_cf_Lep1pt[5] , (*Lep1).Pt()  , evt_weight_ );
                  FillHisto( h_cf_Lep1eta[5], (*Lep1).Eta() , evt_weight_ );
                  FillHisto( h_cf_Lep1phi[5], (*Lep1).Phi() , evt_weight_ );
                  FillHisto( h_cf_Lep2pt[5] , (*Lep2).Pt()  , evt_weight_ );
                  FillHisto( h_cf_Lep2eta[5], (*Lep2).Eta() , evt_weight_ );
                  FillHisto( h_cf_Lep2phi[5], (*Lep2).Phi() , evt_weight_ );
                  FillHisto( h_cf_NPV[5]    , num_pv        , evt_weight_ );
                  FillHisto( h_cf_NJets[5]  , v_jet_idx.size(), evt_weight_ );
                  FillHisto( h_cf_Jet1pt[5] , (*Jet1).Pt()  , evt_weight_ );
                  FillHisto( h_cf_Jet1eta[5], (*Jet1).Eta() , evt_weight_ );
                  FillHisto( h_cf_Jet1phi[5], (*Jet1).Phi() , evt_weight_ );
                  FillHisto( h_cf_Jet2pt[5] , (*Jet2).Pt()  , evt_weight_ );
                  FillHisto( h_cf_Jet2eta[5], (*Jet2).Eta() , evt_weight_ );
                  FillHisto( h_cf_Jet2phi[5], (*Jet2).Phi() , evt_weight_ );
                  
                  FillHisto( h_cf_dilep_inv_mass[5], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
                  FillHisto( h_cf_metpt[5] , Met->Pt() , evt_weight_ );
                  FillHisto( h_cf_metphi[5], Met->Phi(), evt_weight_ );
                  
                  
                  FillHisto( h_Lep1pt[5]  , Lep1->Pt() , evt_weight_ );
                  FillHisto( h_Lep2pt[5]  , Lep2->Pt() , evt_weight_ );
                  FillHisto( h_Lep1eta[5] , Lep1->Eta(), evt_weight_ );
                  FillHisto( h_Lep2eta[5] , Lep2->Eta(), evt_weight_ );
                  FillHisto( h_Lep1phi[5] , Lep1->Phi(), evt_weight_ );
                  FillHisto( h_Lep2phi[5] , Lep2->Phi(), evt_weight_ );
                  
                  if (TString(Decaymode).Contains("muel"))
                  { 
                     FillHisto( h_Muonpt[5]  , TMuon->Pt()      , evt_weight_ );
                     FillHisto( h_Elecpt[5]  , TElectron->Pt()  , evt_weight_ );
                     FillHisto( h_Muoneta[5] , TMuon->Eta()     , evt_weight_ );
                     FillHisto( h_Eleceta[5] , TElectron->Eta() , evt_weight_ );
                     FillHisto( h_Muonphi[5] , TMuon->Phi()     , evt_weight_ );
                     FillHisto( h_Elecphi[5] , TElectron->Phi() , evt_weight_ );
                  }
                  
                  FillHisto( h_Jet1pt[5]  , Jet1->Pt() , evt_weight_ );
                  FillHisto( h_Jet2pt[5]  , Jet2->Pt() , evt_weight_ );
                  FillHisto( h_Jet1eta[5] , Jet1->Eta(), evt_weight_ );
                  FillHisto( h_Jet2eta[5] , Jet2->Eta(), evt_weight_ );
                  FillHisto( h_Jet1phi[5] , Jet1->Phi(), evt_weight_ );
                  FillHisto( h_Jet2phi[5] , Jet2->Phi(), evt_weight_ );
                  FillHisto( h_METpt[5]   , Met->Pt()  , evt_weight_ );
                  FillHisto( h_METphi[5]  , Met->Phi() , evt_weight_ );
                  
                  
                  FillHisto( h_DiLepMass[5], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
                  
                  FillHisto( h_Num_PV[5], num_pv, evt_weight_ );
                  FillHisto( h_Num_Jets[5], v_jet_idx.size(), evt_weight_ );
                  FillHisto( h_Num_bJets[5], nbtagged, evt_weight_ );
                  
                  FillHisto( h_HT[5], AllJetpt, evt_weight_);
                  if (isAllSyst == true) {
                     for ( int i = 0; i < v_SystFullName.size(); ++i )
                     {  
                        if (TString(v_SystFullName[i]).Contains("JetEn")){continue;}
                        FillHisto( h_cf_sys_NLeptons[i][5], v_lepton_idx.size(), v_SystEvt[i] );
                        FillHisto( h_cf_sys_Lep1pt[i][5] , (*Lep1).Pt()  , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Lep1eta[i][5], (*Lep1).Eta() , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Lep1phi[i][5], (*Lep1).Phi() , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Lep2pt[i][5] , (*Lep2).Pt()  , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Lep2eta[i][5], (*Lep2).Eta() , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Lep2phi[i][5], (*Lep2).Phi() , v_SystEvt[i] );
                        FillHisto( h_cf_sys_NPV[i][5]    , num_pv        , v_SystEvt[i] );
                        FillHisto( h_cf_sys_NJets[i][5]  , v_jet_idx.size(), v_SystEvt[i] );
                        FillHisto( h_cf_sys_Jet1pt[i][5] , (*Jet1).Pt()  , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Jet1eta[i][5], (*Jet1).Eta() , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Jet1phi[i][5], (*Jet1).Phi() , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Jet2pt[i][5] , (*Jet2).Pt()  , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Jet2eta[i][5], (*Jet2).Eta() , v_SystEvt[i] );
                        FillHisto( h_cf_sys_Jet2phi[i][5], (*Jet2).Phi() , v_SystEvt[i] );
                        FillHisto( h_cf_sys_metpt[i][5] ,Met->Pt()  , v_SystEvt[i] );
                        FillHisto( h_cf_sys_metphi[i][5],Met->Phi() , v_SystEvt[i] );
                        
                        FillHisto( h_cf_sys_dilep_inv_mass[i][5], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                        FillHisto( h_cf_sys_metpt[i][5] , Met->Pt() , v_SystEvt[i] );
                        FillHisto( h_cf_sys_metphi[i][5], Met->Phi(), v_SystEvt[i] );
                        
                        FillHisto( h_sys_Num_PV[i][5], num_pv, v_SystEvt[i] );
                        FillHisto( h_sys_DiLepMass[i][5], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                        FillHisto( h_sys_Lep1pt[i][5] , (*Lep1).Pt()  , v_SystEvt[i] );
                        FillHisto( h_sys_Lep1eta[i][5], (*Lep1).Eta() , v_SystEvt[i] );
                        FillHisto( h_sys_Lep1phi[i][5], (*Lep1).Phi() , v_SystEvt[i] );
                        FillHisto( h_sys_Lep2pt[i][5] , (*Lep2).Pt()  , v_SystEvt[i] );
                        FillHisto( h_sys_Lep2eta[i][5], (*Lep2).Eta() , v_SystEvt[i] );
                        FillHisto( h_sys_Lep2phi[i][5], (*Lep2).Phi() , v_SystEvt[i] );
                        FillHisto( h_sys_Jet1pt[i][5] , (*Jet1).Pt()  , v_SystEvt[i] );
                        FillHisto( h_sys_Jet1eta[i][5], (*Jet1).Eta() , v_SystEvt[i] );
                        FillHisto( h_sys_Jet1phi[i][5], (*Jet1).Phi() , v_SystEvt[i] );
                        FillHisto( h_sys_Jet2pt[i][5] , (*Jet2).Pt()  , v_SystEvt[i] );
                        FillHisto( h_sys_Jet2eta[i][5], (*Jet2).Eta() , v_SystEvt[i] );
                        FillHisto( h_sys_Jet2phi[i][5], (*Jet2).Phi() , v_SystEvt[i] );
                        FillHisto( h_sys_Num_Jets[i][5], v_jet_idx.size(), v_SystEvt[i] );
                        FillHisto( h_sys_METpt[i][5]   , Met->Pt()  , v_SystEvt[i] );
                        FillHisto( h_sys_METphi[i][5]  , Met->Phi()  , v_SystEvt[i] );
                     }                 
                  }
                  TVector3 lepvec( Lep1->Px() + Lep2->Px(), Lep1->Py() + Lep2->Py(), Lep1->Pz() + Lep2->Pz() );
                  TMatrixD matrix(3,3);
                  
                  double  vec[3];
                  double  lep2 = lepvec.Mag2();
                  double  _apla = 0;
                  double  _sphe = 0;
                  double  _plan = 0;
                  
                  for ( int i=0; i<3; i++)
                  {
                     for ( int j=0; j<3; j++ )
                     {
                        matrix(i,j) = lepvec(i)*lepvec(j);
                        double norm = lep2;
                        for ( int k=0; k<v_jet_idx.size(); k++ )
                        {
                           //TLorentzVector *topoJet = (TLorentzVector*)Jet->At( v_jet_idx[k] );
                           TLorentzVector *topoJet =  v_jet_TL[k];
                           vec[0] = topoJet->Px();
                           vec[1] = topoJet->Py();
                           vec[2] = topoJet->Pz();
                           matrix(i,j) += vec[i]*vec[j];
                           norm        += vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
                        }
                        if (norm > 0) matrix(i,j) /= norm;
                     }
                  }
                  if ( v_jet_idx.size() > 0 )
                  {
                     TMatrixDEigen evmatrix(matrix);
                     TVectorD eigenv=evmatrix.GetEigenValuesRe();
                     matrix.EigenVectors(eigenv);
                  
                     _apla  = 1.5*eigenv(2);
                     _sphe  = 1.5*(eigenv(2) + eigenv(1));
                     _plan  = eigenv(1) - eigenv(2);
                  }
                  
                  FillHisto(h_Topo_Apla[0] , _apla, evt_weight_);
                  FillHisto(h_Topo_Sphe[0] , _sphe, evt_weight_);
                  FillHisto(h_Topo_Plan[0] , _plan, evt_weight_);
                  
                  ///////////////////////////////////////
                  /// 2 or More B-Tagging Requriement ///
                  /////////////////////////////////////// 
                  
                  if ( DoubleBtag(v_bjet_idx) == true ) 
                  {
                     FillHisto( h_EventWeight[6], evt_weight_  );
                     FillHisto( h_Lep1pt[6] , Lep1->Pt() , evt_weight_ );
                     FillHisto( h_Lep2pt[6] , Lep2->Pt() , evt_weight_ );
                     FillHisto( h_Lep1eta[6], Lep1->Eta(), evt_weight_ );
                     FillHisto( h_Lep2eta[6], Lep2->Eta(), evt_weight_ );
                     FillHisto( h_Lep1phi[6], Lep1->Phi(), evt_weight_ );
                     FillHisto( h_Lep2phi[6], Lep2->Phi(), evt_weight_ );
                     
                     if (TString(Decaymode).Contains("muel"))
                     { 
                        FillHisto( h_Muonpt[6]  , TMuon->Pt()      , evt_weight_ );
                        FillHisto( h_Elecpt[6]  , TElectron->Pt()  , evt_weight_ );
                        FillHisto( h_Muoneta[6] , TMuon->Eta()     , evt_weight_ );
                        FillHisto( h_Eleceta[6] , TElectron->Eta() , evt_weight_ );
                        FillHisto( h_Muonphi[6] , TMuon->Phi()     , evt_weight_ );
                        FillHisto( h_Elecphi[6] , TElectron->Phi() , evt_weight_ );
                     }
                     FillHisto( h_Jet1pt[6] , Jet1->Pt() , evt_weight_ );
                     FillHisto( h_Jet2pt[6] , Jet2->Pt() , evt_weight_ );
                     FillHisto( h_Jet1eta[6], Jet1->Eta(), evt_weight_ );
                     FillHisto( h_Jet2eta[6], Jet2->Eta(), evt_weight_ );
                     FillHisto( h_Jet1phi[6], Jet1->Phi(), evt_weight_ );
                     FillHisto( h_Jet2phi[6], Jet2->Phi(), evt_weight_ );
                     FillHisto( h_METpt[6]  , Met->Pt()  , evt_weight_ );
                     FillHisto( h_METphi[6] , Met->Phi() , evt_weight_ );
                     FillHisto( h_HT[6] , AllJetpt, evt_weight_);
                   
                     FillHisto( h_DiLepMass[6], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
                     FillHisto( h_Num_PV[6]   , num_pv          , evt_weight_ );
                     FillHisto( h_Num_Jets[6] , v_jet_idx.size(), evt_weight_ );
                     FillHisto( h_Num_bJets[6], nbtagged        , evt_weight_ );
                   
                     FillHisto( h_cf_NLeptons[6], v_lepton_idx.size(), evt_weight_ );
                     FillHisto( h_cf_Lep1pt[6] , (*Lep1).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Lep1eta[6], (*Lep1).Eta() , evt_weight_ );
                     FillHisto( h_cf_Lep1phi[6], (*Lep1).Phi() , evt_weight_ );
                     FillHisto( h_cf_Lep2pt[6] , (*Lep2).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Lep2eta[6], (*Lep2).Eta() , evt_weight_ );
                     FillHisto( h_cf_Lep2phi[6], (*Lep2).Phi() , evt_weight_ );
                     FillHisto( h_cf_NPV[6]    , num_pv        , evt_weight_ );
                     FillHisto( h_cf_NJets[6]  , v_jet_idx.size(), evt_weight_ );
                     FillHisto( h_cf_Jet1pt[6] , (*Jet1).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Jet1eta[6], (*Jet1).Eta() , evt_weight_ );
                     FillHisto( h_cf_Jet1phi[6], (*Jet1).Phi() , evt_weight_ );
                     FillHisto( h_cf_Jet2pt[6] , (*Jet2).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Jet2eta[6], (*Jet2).Eta() , evt_weight_ );
                     FillHisto( h_cf_Jet2phi[6], (*Jet2).Phi() , evt_weight_ );
                   
                     FillHisto( h_cf_dilep_inv_mass[6], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
                     FillHisto( h_cf_metpt[6] , Met->Pt() , evt_weight_ );
                     FillHisto( h_cf_metphi[6], Met->Phi(), evt_weight_ );
                   
                   
                     FillHisto( h_Topo_Apla[1] , _apla, evt_weight_);
                     FillHisto( h_Topo_Sphe[1] , _sphe, evt_weight_);
                     FillHisto( h_Topo_Plan[1] , _plan, evt_weight_);
                     if (isAllSyst == true) {
                        for ( int i = 0; i < v_SystFullName.size(); ++i )
                        {  
                           if (TString(v_SystFullName[i]).Contains("JetEn")){continue;}
//                           if ( TString(v_SystFullName[i]).Contains("PileUpUp") ) { cout << "evt_weight_ : " << v_SystEvt[i] << " Jet1 Pt : " << Jet1->Pt() << endl;}
                           FillHisto( h_cf_sys_NLeptons[i][6], v_lepton_idx.size(), v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep1pt[i][6] , (*Lep1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep1eta[i][6], (*Lep1).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep1phi[i][6], (*Lep1).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep2pt[i][6] , (*Lep2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep2eta[i][6], (*Lep2).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep2phi[i][6], (*Lep2).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_NPV[i][6]    , num_pv        , v_SystEvt[i] );
                           FillHisto( h_cf_sys_NJets[i][6]  , v_jet_idx.size(), v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet1pt[i][6] , (*Jet1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet1eta[i][6], (*Jet1).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet1phi[i][6], (*Jet1).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet2pt[i][6] , (*Jet2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet2eta[i][6], (*Jet2).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet2phi[i][6], (*Jet2).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_metpt[i][6] ,Met->Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_metphi[i][6],Met->Phi() , v_SystEvt[i] );
                           
                           FillHisto( h_cf_sys_dilep_inv_mass[i][6], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                           FillHisto( h_cf_sys_metpt[i][6] , Met->Pt() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_metphi[i][6], Met->Phi(), v_SystEvt[i] );
                           
                           FillHisto( h_sys_Num_PV[i][6], num_pv, v_SystEvt[i] );
                           FillHisto( h_sys_DiLepMass[i][6], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                           FillHisto( h_sys_Lep1pt[i][6] , (*Lep1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Lep1eta[i][6], (*Lep1).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Lep1phi[i][6], (*Lep1).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Lep2pt[i][6] , (*Lep2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Lep2eta[i][6], (*Lep2).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Lep2phi[i][6], (*Lep2).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet1pt[i][6] , (*Jet1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Jet1eta[i][6], (*Jet1).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet1phi[i][6], (*Jet1).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet2pt[i][6] , (*Jet2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Jet2eta[i][6], (*Jet2).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet2phi[i][6], (*Jet2).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Num_Jets[i][6], v_jet_idx.size(), v_SystEvt[i] );
                           FillHisto( h_sys_METpt[i][6]   , Met->Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_METphi[i][6]  , Met->Phi()  , v_SystEvt[i] );
                        }
                     }
                  } // One or More b-tagging // Step 5 
                   
                  if (nbtagged ==2)// exactly 2b-tagging ..
                  {
//                     ApplyBTagWeight(-2);
                  
                     FillHisto( h_EventWeight[7], evt_weight_  );
                     FillHisto( h_cf_NLeptons[7], v_lepton_idx.size(), evt_weight_ );
                  
                     FillHisto( h_cf_Lep1pt[7] , (*Lep1).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Lep1eta[7], (*Lep1).Eta() , evt_weight_ );
                     FillHisto( h_cf_Lep1phi[7], (*Lep1).Phi() , evt_weight_ );
                     FillHisto( h_cf_Lep2pt[7] , (*Lep2).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Lep2eta[7], (*Lep2).Eta() , evt_weight_ );
                     FillHisto( h_cf_Lep2phi[7], (*Lep2).Phi() , evt_weight_ );
                     FillHisto( h_cf_NPV[7]    , num_pv        , evt_weight_ );
                     FillHisto( h_cf_NJets[7]  , v_jet_idx.size(), evt_weight_ );
                     FillHisto( h_cf_Jet1pt[7] , (*Jet1).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Jet1eta[7], (*Jet1).Eta() , evt_weight_ );
                     FillHisto( h_cf_Jet1phi[7], (*Jet1).Phi() , evt_weight_ );
                     FillHisto( h_cf_Jet2pt[7] , (*Jet2).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Jet2eta[7], (*Jet2).Eta() , evt_weight_ );
                     FillHisto( h_cf_Jet2phi[7], (*Jet2).Phi() , evt_weight_ );
                  
                     FillHisto( h_cf_dilep_inv_mass[7], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
                     FillHisto( h_cf_metpt[7] , Met->Pt() , evt_weight_ );
                     FillHisto( h_cf_metphi[7], Met->Phi(), evt_weight_ );
                  
                     FillHisto( h_Lep1pt[7] , Lep1->Pt() , evt_weight_ );
                     FillHisto( h_Lep2pt[7] , Lep2->Pt() , evt_weight_ );
                     FillHisto( h_Lep1eta[7], Lep1->Eta(), evt_weight_ );
                     FillHisto( h_Lep2eta[7], Lep2->Eta(), evt_weight_ );
                     FillHisto( h_Lep1phi[7], Lep1->Phi(), evt_weight_ );
                     FillHisto( h_Lep2phi[7], Lep2->Phi(), evt_weight_ );
                    
                     if (TString(Decaymode).Contains("muel"))
                     { 
                        FillHisto( h_Muonpt[7]  , TMuon->Pt()      , evt_weight_ );
                        FillHisto( h_Elecpt[7]  , TElectron->Pt()  , evt_weight_ );
                        FillHisto( h_Muoneta[7] , TMuon->Eta()     , evt_weight_ );
                        FillHisto( h_Eleceta[7] , TElectron->Eta() , evt_weight_ );
                        FillHisto( h_Muonphi[7] , TMuon->Phi()     , evt_weight_ );
                        FillHisto( h_Elecphi[7] , TElectron->Phi() , evt_weight_ );
                     }
                     FillHisto( h_Jet1pt[7] , Jet1->Pt() , evt_weight_ );
                     FillHisto( h_Jet2pt[7] , Jet2->Pt() , evt_weight_ );
                     FillHisto( h_Jet1eta[7], Jet1->Eta(), evt_weight_ );
                     FillHisto( h_Jet2eta[7], Jet2->Eta(), evt_weight_ );
                     FillHisto( h_Jet1phi[7], Jet1->Phi(), evt_weight_ );
                     FillHisto( h_Jet2phi[7], Jet2->Phi(), evt_weight_ );
                     FillHisto( h_METpt[7]  , Met->Pt()  , evt_weight_ );
                     FillHisto( h_METphi[7] , Met->Phi() , evt_weight_ );
                     FillHisto( h_HT[7] , AllJetpt, evt_weight_);
                  
                     FillHisto( h_DiLepMass[7], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
                  
                     FillHisto( h_Num_PV[7]   , num_pv          , evt_weight_ );
                     FillHisto( h_Num_Jets[7] , v_jet_idx.size(), evt_weight_ );
                     FillHisto( h_Num_bJets[7], nbtagged        , evt_weight_ );
                  
                     FillHisto( h_Topo_Apla[2] , _apla, evt_weight_);
                     FillHisto( h_Topo_Sphe[2] , _sphe, evt_weight_);
                     FillHisto( h_Topo_Plan[2] , _plan, evt_weight_);
                     if (isAllSyst == true) {
                        for ( int i = 0; i < v_SystFullName.size(); ++i )
                        {  
                           if (TString(v_SystFullName[i]).Contains("JetEn")){continue;}
                           FillHisto( h_cf_sys_NLeptons[i][7], v_lepton_idx.size(), v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep1pt[i][7] , (*Lep1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep1eta[i][7], (*Lep1).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep1phi[i][7], (*Lep1).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep2pt[i][7] , (*Lep2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep2eta[i][7], (*Lep2).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep2phi[i][7], (*Lep2).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_NPV[i][7]    , num_pv        , v_SystEvt[i] );
                           FillHisto( h_cf_sys_NJets[i][7]  , v_jet_idx.size(), v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet1pt[i][7] , (*Jet1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet1eta[i][7], (*Jet1).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet1phi[i][7], (*Jet1).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet2pt[i][7] , (*Jet2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet2eta[i][7], (*Jet2).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet2phi[i][7], (*Jet2).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_metpt[i][7] ,Met->Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_metphi[i][7],Met->Phi() , v_SystEvt[i] );
                           
                           FillHisto( h_cf_sys_dilep_inv_mass[i][7], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                           FillHisto( h_cf_sys_metpt[i][7] , Met->Pt() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_metphi[i][7], Met->Phi(), v_SystEvt[i] );
                           
                           FillHisto( h_sys_Num_PV[i][7], num_pv, v_SystEvt[i] );
                           FillHisto( h_sys_DiLepMass[i][7], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                           FillHisto( h_sys_Lep1pt[i][7] , (*Lep1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Lep1eta[i][7], (*Lep1).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Lep1phi[i][7], (*Lep1).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Lep2pt[i][7] , (*Lep2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Lep2eta[i][7], (*Lep2).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Lep2phi[i][7], (*Lep2).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet1pt[i][7] , (*Jet1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Jet1eta[i][7], (*Jet1).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet1phi[i][7], (*Jet1).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet2pt[i][7] , (*Jet2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Jet2eta[i][7], (*Jet2).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet2phi[i][7], (*Jet2).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Num_Jets[i][7], v_jet_idx.size(), v_SystEvt[i] );
                           FillHisto( h_sys_METpt[i][7]   , Met->Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_METphi[i][7]  , Met->Phi()  , v_SystEvt[i] );
                        }
                     }
                  }         
                  // Reconstruction of Top with KinSol.
                  KinSol();   
                  if ( ksolweight_ != -1 )
                  {
                     FillHisto( h_EventWeight[8], evt_weight_  );
                     FillHisto( h_cf_NLeptons[8], v_lepton_idx.size(), evt_weight_ );
                  
                     FillHisto( h_cf_Lep1pt[8] , (*Lep1).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Lep1eta[8], (*Lep1).Eta() , evt_weight_ );
                     FillHisto( h_cf_Lep1phi[8], (*Lep1).Phi() , evt_weight_ );
                     FillHisto( h_cf_Lep2pt[8] , (*Lep2).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Lep2eta[8], (*Lep2).Eta() , evt_weight_ );
                     FillHisto( h_cf_Lep2phi[8], (*Lep2).Phi() , evt_weight_ );
                     FillHisto( h_cf_NPV[8]    , num_pv        , evt_weight_ );
                     FillHisto( h_cf_NJets[8]  , v_jet_idx.size(), evt_weight_ );
                     FillHisto( h_cf_Jet1pt[8] , (*Jet1).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Jet1eta[8], (*Jet1).Eta() , evt_weight_ );
                     FillHisto( h_cf_Jet1phi[8], (*Jet1).Phi() , evt_weight_ );
                     FillHisto( h_cf_Jet2pt[8] , (*Jet2).Pt()  , evt_weight_ );
                     FillHisto( h_cf_Jet2eta[8], (*Jet2).Eta() , evt_weight_ );
                     FillHisto( h_cf_Jet2phi[8], (*Jet2).Phi() , evt_weight_ );
                  
                     FillHisto( h_cf_dilep_inv_mass[8], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
                     FillHisto( h_cf_metpt[8] , Met->Pt() , evt_weight_ );
                     FillHisto( h_cf_metphi[8], Met->Phi(), evt_weight_ );
                  
                     FillHisto( h_Lep1pt[8] , Lep1->Pt() , evt_weight_ );
                     FillHisto( h_Lep2pt[8] , Lep2->Pt() , evt_weight_ );
                     FillHisto( h_Lep1eta[8], Lep1->Eta(), evt_weight_ );
                     FillHisto( h_Lep2eta[8], Lep2->Eta(), evt_weight_ );
                     FillHisto( h_Lep1phi[8], Lep1->Phi(), evt_weight_ );
                     FillHisto( h_Lep2phi[8], Lep2->Phi(), evt_weight_ );
                    
                     if (TString(Decaymode).Contains("muel"))
                     { 
                        FillHisto( h_Muonpt[8]  , TMuon->Pt()      , evt_weight_ );
                        FillHisto( h_Elecpt[8]  , TElectron->Pt()  , evt_weight_ );
                        FillHisto( h_Muoneta[8] , TMuon->Eta()     , evt_weight_ );
                        FillHisto( h_Eleceta[8] , TElectron->Eta() , evt_weight_ );
                        FillHisto( h_Muonphi[8] , TMuon->Phi()     , evt_weight_ );
                        FillHisto( h_Elecphi[8] , TElectron->Phi() , evt_weight_ );
                     }
                     FillHisto( h_Jet1pt[8] , Jet1->Pt() , evt_weight_ );
                     FillHisto( h_Jet2pt[8] , Jet2->Pt() , evt_weight_ );
                     FillHisto( h_Jet1eta[8], Jet1->Eta(), evt_weight_ );
                     FillHisto( h_Jet2eta[8], Jet2->Eta(), evt_weight_ );
                     FillHisto( h_Jet1phi[8], Jet1->Phi(), evt_weight_ );
                     FillHisto( h_Jet2phi[8], Jet2->Phi(), evt_weight_ );
                     FillHisto( h_METpt[8]  , Met->Pt()  , evt_weight_ );
                     FillHisto( h_METphi[8] , Met->Phi() , evt_weight_ );
                     FillHisto( h_HT[8]     , AllJetpt   , evt_weight_);
                  
                     FillHisto( h_DiLepMass[8], ( (*Lep1)+(*Lep2) ).M(), evt_weight_ );
                  
                     FillHisto( h_Num_PV[8]   , num_pv          , evt_weight_ );
                     FillHisto( h_Num_Jets[8] , v_jet_idx.size(), evt_weight_ );
                     FillHisto( h_Num_bJets[8], nbtagged        , evt_weight_ );
                  
                     FillHisto( h_Topo_Apla[3] , _apla, evt_weight_);
                     FillHisto( h_Topo_Sphe[3] , _sphe, evt_weight_);
                     FillHisto( h_Topo_Plan[3] , _plan, evt_weight_);
                     if ( Top->Pt() > AnTop->Pt() ) { (*Top1) = (*Top); (*Top2) = (*AnTop); }
                     else { (*Top1) = (*AnTop); (*Top2) = (*Top); }

                     if (isAllSyst == true){
                        for ( int i = 0; i < v_SystFullName.size(); ++i )
                        {  
                           if (TString(v_SystFullName[i]).Contains("JetEn")){continue;}
                           FillHisto( h_cf_sys_NLeptons[i][8], v_lepton_idx.size(), v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep1pt[i][8] , (*Lep1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep1eta[i][8], (*Lep1).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep1phi[i][8], (*Lep1).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep2pt[i][8] , (*Lep2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep2eta[i][8], (*Lep2).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Lep2phi[i][8], (*Lep2).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_NPV[i][8]    , num_pv        , v_SystEvt[i] );
                           FillHisto( h_cf_sys_NJets[i][8]  , v_jet_idx.size(), v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet1pt[i][8] , (*Jet1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet1eta[i][8], (*Jet1).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet1phi[i][8], (*Jet1).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet2pt[i][8] , (*Jet2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet2eta[i][8], (*Jet2).Eta() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_Jet2phi[i][8], (*Jet2).Phi() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_metpt[i][8] ,Met->Pt()  , v_SystEvt[i] );
                           FillHisto( h_cf_sys_metphi[i][8],Met->Phi() , v_SystEvt[i] );
                           
                           FillHisto( h_cf_sys_dilep_inv_mass[i][8], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                           FillHisto( h_cf_sys_metpt[i][8] , Met->Pt() , v_SystEvt[i] );
                           FillHisto( h_cf_sys_metphi[i][8], Met->Phi(), v_SystEvt[i] );
                           
                           FillHisto( h_sys_Num_PV[i][8], num_pv, v_SystEvt[i] );
                           FillHisto( h_sys_DiLepMass[i][8], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[i] );
                           FillHisto( h_sys_Lep1pt[i][8] , (*Lep1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Lep1eta[i][8], (*Lep1).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Lep1phi[i][8], (*Lep1).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Lep2pt[i][8] , (*Lep2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Lep2eta[i][8], (*Lep2).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Lep2phi[i][8], (*Lep2).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet1pt[i][8] , (*Jet1).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Jet1eta[i][8], (*Jet1).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet1phi[i][8], (*Jet1).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet2pt[i][8] , (*Jet2).Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_Jet2eta[i][8], (*Jet2).Eta() , v_SystEvt[i] );
                           FillHisto( h_sys_Jet2phi[i][8], (*Jet2).Phi() , v_SystEvt[i] );
                           FillHisto( h_sys_Num_Jets[i][8], v_jet_idx.size(), v_SystEvt[i] );
                           FillHisto( h_sys_METpt[i][8]   , Met->Pt()  , v_SystEvt[i] );
                           FillHisto( h_sys_METphi[i][8]  , Met->Phi()  , v_SystEvt[i] );

                           FillHisto( h_sys_Top1Mass_[i]    , Top1->M()        , v_SystEvt[i] );
                           FillHisto( h_sys_Top1pt_[i]      , Top1->Pt()       , v_SystEvt[i] );
                           FillHisto( h_sys_Top1phi_[i]     , Top1->Phi()      , v_SystEvt[i] );
                           FillHisto( h_sys_Top1Rapidity_[i], Top1->Rapidity() , v_SystEvt[i] );
                           FillHisto( h_sys_Top1Energy_[i]  , Top1->Energy()   , v_SystEvt[i] );
                           FillHisto( h_sys_Top2Mass_[i]    , Top2->M()        , v_SystEvt[i] );
                           FillHisto( h_sys_Top2pt_[i]      , Top2->Pt()       , v_SystEvt[i] );
                           FillHisto( h_sys_Top2phi_[i]     , Top2->Phi()      , v_SystEvt[i] );
                           FillHisto( h_sys_Top2Rapidity_[i], Top2->Rapidity() , v_SystEvt[i] );
                           FillHisto( h_sys_Top2Energy_[i]  , Top2->Energy()   , v_SystEvt[i] );

                           FillHisto( h_sys_TopMass_[i]      , Top->M()         , v_SystEvt[i] );
                           FillHisto( h_sys_Toppt_[i]        , Top->Pt()        , v_SystEvt[i] );
                           FillHisto( h_sys_Topphi_[i]       , Top->Phi()       , v_SystEvt[i] );
                           FillHisto( h_sys_TopRapidity_[i]  , Top->Rapidity()  , v_SystEvt[i] );
                           FillHisto( h_sys_TopEnergy_[i]    , Top->Energy()    , v_SystEvt[i] );
                           FillHisto( h_sys_AnTopMass_[i]    , AnTop->M()       , v_SystEvt[i] );
                           FillHisto( h_sys_AnToppt_[i]      , AnTop->Pt()      , v_SystEvt[i] );
                           FillHisto( h_sys_AnTopphi_[i]     , AnTop->Phi()     , v_SystEvt[i] );
                           FillHisto( h_sys_AnTopRapidity_[i], AnTop->Rapidity(), v_SystEvt[i] );
                           FillHisto( h_sys_AnTopEnergy_[i]  , AnTop->Energy()  , v_SystEvt[i] );

                        }
                     }
                  
                     FillHisto( h_TopMass      , Top->M()         , evt_weight_ );
                     FillHisto( h_Toppt        , Top->Pt()        , evt_weight_ );
                     FillHisto( h_Topphi       , Top->Phi()       , evt_weight_ );
                     FillHisto( h_TopRapidity  , Top->Rapidity()  , evt_weight_ );
                     FillHisto( h_TopEnergy    , Top->Energy()    , evt_weight_ );
                     FillHisto( h_AnTopMass    , AnTop->M()       , evt_weight_ );
                     FillHisto( h_AnToppt      , AnTop->Pt()      , evt_weight_ );
                     FillHisto( h_AnTopphi     , AnTop->Phi()     , evt_weight_ );
                     FillHisto( h_AnTopRapidity, AnTop->Rapidity(), evt_weight_ );
                     FillHisto( h_AnTopEnergy  , AnTop->Energy()  , evt_weight_ );
                  
                  
                     FillHisto( h_W1Mass , W1->M()  , evt_weight_ );
                     FillHisto( h_W2Mass , W2->M()  , evt_weight_ );
                  
                     FillHisto( h_W1Mt , W1->Mt()  , evt_weight_ );
                     FillHisto( h_W2Mt , W2->Mt()  , evt_weight_ );
                  
                     FillHisto( h_bJet1Energy , bJet1->Energy()  , evt_weight_ );
                     FillHisto( h_bJet2Energy , bJet2->Energy()  , evt_weight_ );
                  
                     FillHisto( h_bJetEnergy   , bJet->Energy()   , evt_weight_ );
                     FillHisto( h_AnbJetEnergy , AnbJet->Energy() , evt_weight_ );
                     FillHisto( h_bJetPt       , bJet->Pt()   , evt_weight_ );
                     FillHisto( h_AnbJetPt     , AnbJet->Pt() , evt_weight_ );
                     FillHisto( h_LepEnergy    , Lep->Energy()    , evt_weight_ );
                     FillHisto( h_AnLepEnergy  , AnLep->Energy()  , evt_weight_ );
                     FillHisto( h_NuEnergy     , Nu->Energy()     , evt_weight_ );
                     FillHisto( h_AnNuEnergy   , AnNu->Energy()   , evt_weight_ );
                  
                     /// Kin Solver Purity ... ///
                     FillHisto( h2_TopMassVsLepBMass    , Top->M() ,((*Lep)+(*AnbJet)).M()   , evt_weight_ );
                     FillHisto( h2_AnTopMassVsLepBMass  , AnTop->M() ,((*AnLep)+(*bJet)).M()   , evt_weight_ );
                     FillHisto( h2_AnLepBMassVsLepBMass , ((*Lep)+(*AnbJet)).M() ,((*AnLep)+(*bJet)).M()   , evt_weight_ );
                  
                     FillHisto( h_LepbJetMass   , ((*Lep)+(*AnbJet)).M() , evt_weight_ );
                     FillHisto( h_AnLepbJetMass , ((*AnLep)+(*bJet)).M() , evt_weight_ );
                  
                     LepAnLepMisCharge();
                     
                  
                     // Calculating O3 variable
                     double cpO3 = ssbcpviol->getO3Vari( bJet , AnbJet, AnLep, Lep );
                     double cpO3_jprup = ssbcpviol->getO3VariJPRUp( bJet ,AnbJet , AnLep ,Lep );
                     double cpO3_jprdown = ssbcpviol->getO3VariJPRDown( bJet ,AnbJet , AnLep ,Lep );
                  
                     double cpOb = ssbcpviol->getObVari( bJet , AnbJet, AnLep, Lep );
                     double cpO5 = ssbcpviol->getO5Vari( bJet , AnbJet, AnLep, Lep );
                     if (cpO3 == 0.0 ){cout << "--00000000000000--" << cpO3 << endl;}
                     // Calculating O3 error. 
                     if (cpO3 > 0.0) {Num_cpo3_p++;
      //                  cout << "evt_weight_ ^2 " << evt_weight_*evt_weight_ << endl;
                        FillHisto( h_CPO3_evtweight,0.5,evt_weight_);
                        FillHisto( h_CPO3_evtweight,6.5,evt_weight_*evt_weight_);
                        FillHisto( h_CPO3_Plus,cpO3,evt_weight_);
                     }
                     else           { Num_cpo3_m++; Nm_eventw_2 += (evt_weight_*evt_weight_);
                        FillHisto( h_CPO3_evtweight,2.5,evt_weight_);
                        FillHisto( h_CPO3_evtweight,8.5,evt_weight_*evt_weight_);
                        laecO3 += 50;
                        FillHisto( h_CPO3_Minus,cpO3,evt_weight_);
                     }
                  
                     if (cpOb > 0.0) {
                        FillHisto( h_CPOb_evtweight,0.5,evt_weight_);
                        FillHisto( h_CPOb_evtweight,6.5,evt_weight_*evt_weight_);
                     }
                     else           { 
                        FillHisto( h_CPOb_evtweight,2.5,evt_weight_);
                        FillHisto( h_CPOb_evtweight,8.5,evt_weight_*evt_weight_);
                        laecOb += 50;
                     }
                  
                     if (cpO5 > 0.0) {
                        FillHisto( h_CPO5_evtweight,0.5,evt_weight_);
                        FillHisto( h_CPO5_evtweight,6.5,evt_weight_*evt_weight_);
                     }
                     else           { 
                        FillHisto( h_CPO5_evtweight,2.5,evt_weight_);
                        FillHisto( h_CPO5_evtweight,8.5,evt_weight_*evt_weight_);
                        laecO5 += 50;
                     }
                  
                     /// Jet pT Variation Study For Dilution
                     if (cpO3 > 0.0 && cpO3_jprup > 0.0 )FillHisto( h_CPO3_JPRUp_Plus_Plus ,  cpO3_jprup, evt_weight_ );
                     if (cpO3 > 0.0 && cpO3_jprup < 0.0 )FillHisto( h_CPO3_JPRUp_Plus_Minus,  cpO3_jprup, evt_weight_ );
                     if (cpO3 < 0.0 && cpO3_jprup < 0.0 )FillHisto( h_CPO3_JPRUp_Minus_Minus, cpO3_jprup, evt_weight_ );
                     if (cpO3 < 0.0 && cpO3_jprup > 0.0 )FillHisto( h_CPO3_JPRUp_Minus_Plus,  cpO3_jprup, evt_weight_ );
                  
                     if (cpO3 > 0.0 && cpO3_jprdown > 0.0 )FillHisto( h_CPO3_JPRDown_Plus_Plus ,  cpO3_jprdown, evt_weight_ );
                     if (cpO3 > 0.0 && cpO3_jprdown < 0.0 )FillHisto( h_CPO3_JPRDown_Plus_Minus,  cpO3_jprdown, evt_weight_ );
                     if (cpO3 < 0.0 && cpO3_jprdown < 0.0 )FillHisto( h_CPO3_JPRDown_Minus_Minus, cpO3_jprdown, evt_weight_ );
                     if (cpO3 < 0.0 && cpO3_jprdown > 0.0 )FillHisto( h_CPO3_JPRDown_Minus_Plus,  cpO3_jprdown, evt_weight_ );
                  
                  
//                     cout << "laec ? " << laec << endl;
                     
                     // Calculating Bjorken variables
                     Bjorken( Lep , AnLep , bJet, AnbJet , Nu , AnNu ); 
                  
                     FillHisto( h_CPO3_reco        , cpO3         , evt_weight_ );
                     FillHisto( h_CPO3_reco_JPRUp  , cpO3_jprup   , evt_weight_ );
                     FillHisto( h_CPO3_reco_JPRDown, cpO3_jprdown , evt_weight_ );
                     FillHisto( h_CPOb_reco, ssbcpviol->getObVari( bJet , AnbJet, AnLep, Lep ) , evt_weight_ );
                     FillHisto( h_CPO5_reco, ssbcpviol->getO5Vari( bJet , AnbJet, AnLep, Lep ) , evt_weight_ );
                     FillHisto( h_BjorkenX1, x1_bj, evt_weight_ );
                     FillHisto( h_BjorkenX2, x2_bj, evt_weight_ );
                     FillHisto( h_BjorkenX3, x3_bj, evt_weight_ );
                     
                     // Get CP-Violation Variables //
                     v_recocp_O.clear();
                     v_recocp_O.push_back( ssbcpviol->getO1Vari( Top, AnTop, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO2Vari( Top, AnTop, bJet, AnbJet ) );
                     v_recocp_O.push_back( ssbcpviol->getO3Vari( bJet, AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO4Vari( AnbJet, bJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO5Vari( bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO6Vari( bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO7Vari( Top , AnTop, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO8Vari( Top, AnTop, bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO9Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO10Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO11Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO12Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO13Vari( bJet , AnbJet, AnLep, Lep )  );
                     for (int i = 0; i < v_recocp_O.size(); ++ i)
                     {
                        FillHisto( h_Reco_CPO_[i], v_recocp_O[i] , evt_weight_ );
                        FillHisto( h_Reco_CPO_ReRange_[i], v_recocp_O[i] , evt_weight_ );
                     }
                     if (isAllSyst == true){
                        for (int i = 0; i <  v_SystEvt.size(); ++ i)
                        {
                           if (TString(v_SystFullName[i]).Contains("JetEn")){continue;}
                           for (int j = 0; j < v_recocp_O.size(); ++ j )
                           {
                              FillHisto( h_sys_Reco_CPO_[i][j], v_recocp_O[j] , v_SystEvt[i]   );
                              FillHisto( h_sys_Reco_CPO_ReRange_[i][j], v_recocp_O[j] , v_SystEvt[i]  );
                           }
                        }
                     } 
                     FillHisto( h_Top1Mass    , Top1->M()        , evt_weight_ );
                     FillHisto( h_Top1pt      , Top1->Pt()       , evt_weight_ );
                     FillHisto( h_Top1phi     , Top1->Phi()      , evt_weight_ );
                     FillHisto( h_Top1Rapidity, Top1->Rapidity() , evt_weight_ );
                     FillHisto( h_Top1Energy  , Top1->Energy()   , evt_weight_ );
                  
                     FillHisto( h_Top2Mass    , Top2->M()        , evt_weight_ );
                     FillHisto( h_Top2pt      , Top2->Pt()       , evt_weight_ );
                     FillHisto( h_Top2phi     , Top2->Phi()      , evt_weight_ );
                     FillHisto( h_Top2Rapidity, Top2->Rapidity() , evt_weight_ );
                     FillHisto( h_Top2Energy  , Top2->Energy()   , evt_weight_ );
                     // BJet - BJet //
                     FillHisto( h_Diff_BBar_Pt    , fabs( bJet->Pt() - AnbJet->Pt() )         , evt_weight_ );
                     FillHisto( h_Diff_BBar_Energy, fabs( bJet->Energy() - AnbJet->Energy() ) , evt_weight_ );
                     FillHisto( h_Diff_BBar_P     , fabs( bJet->P() - AnbJet->P() )           , evt_weight_ );
                  
                     FillHisto( h_Diff_LepAnLep_Pt    , fabs( Lep->Pt() - AnLep->Pt() )         , evt_weight_ );
                     FillHisto( h_Diff_LepAnLep_Energy, fabs( Lep->Energy() - AnLep->Energy() ) , evt_weight_ );
                     FillHisto( h_Diff_LepAnLep_P     , fabs( Lep->P() - AnLep->P() )           , evt_weight_ );
                  
                     FillHisto( h_LepAnLepEngCheckO3    , laecO3 , evt_weight_ );
                     FillHisto( h_LepAnLepEngCheckOb    , laecOb , evt_weight_ );
                     FillHisto( h_LepAnLepEngCheckO5    , laecO5 , evt_weight_ );
                  
                     ////////////////////////////////
                     /// JetAngular Study Area ...///
                     ////////////////////////////////
                     double recocpO3_etavariup   = ssbcpviol->getO3Vari( ssbcpviol->JetAngEta(bJet,"up"  ),ssbcpviol->JetAngEta(AnbJet,"up"  ), AnLep, Lep );
                     double recocpO3_etavaridown = ssbcpviol->getO3Vari( ssbcpviol->JetAngEta(bJet,"down"),ssbcpviol->JetAngEta(AnbJet,"down"), AnLep, Lep );
                     double recocpO3_phivariup   = ssbcpviol->getO3Vari( ssbcpviol->JetAngPhi(bJet,"up"  ),ssbcpviol->JetAngPhi(AnbJet,"up"  ), AnLep, Lep );
                     double recocpO3_phivaridown = ssbcpviol->getO3Vari( ssbcpviol->JetAngPhi(bJet,"down"),ssbcpviol->JetAngPhi(AnbJet,"down"), AnLep, Lep );
                  
                     FillHisto( h_CPO3_EtaVariUp,  recocpO3_etavariup,  evt_weight_);
                     FillHisto( h_CPO3_EtaVariDown,recocpO3_etavaridown,evt_weight_);
                     FillHisto( h_CPO3_PhiVariUp,  recocpO3_phivariup,  evt_weight_);
                     FillHisto( h_CPO3_PhiVariDown,recocpO3_phivaridown,evt_weight_);
                     if (cpO3 > 0.0) {
                        FillHisto( h_CPO3reco_evtweight,0.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight,6.5,evt_weight_*evt_weight_);
                     }
                     else{
                        FillHisto( h_CPO3reco_evtweight,2.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight,8.5,evt_weight_*evt_weight_);
                     }
                  
                     //////////////////////////////////////////////////////////// 
                     /// CPO3 with Eta Variation (Up) at Reconstruction Level ///
                     //////////////////////////////////////////////////////////// 
                     if ( recocpO3_etavariup > 0.0) {
                        FillHisto( h_CPO3reco_evtweight_EtaVariUp,0.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_EtaVariUp,6.5,evt_weight_*evt_weight_);
                     }
                     else{
                        FillHisto( h_CPO3reco_evtweight_EtaVariUp,2.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_EtaVariUp,8.5,evt_weight_*evt_weight_);
                     }
                  
                     if ( cpO3 > 0.0  && recocpO3_etavariup > 0.0 )
                     {
                        FillHisto( h_CPO3_EtaVariUp_Plus_Plus,recocpO3_etavariup,evt_weight_);
                     }
                     if ( cpO3 > 0.0  && recocpO3_etavariup < 0.0 )
                     {
                        FillHisto( h_CPO3_EtaVariUp_Plus_Minus,recocpO3_etavariup,evt_weight_);
                     }
                  
                     if ( cpO3 < 0.0  && recocpO3_etavariup > 0.0 )
                     {
                        FillHisto( h_CPO3_EtaVariUp_Minus_Plus,recocpO3_etavariup,evt_weight_);
                     }
                     if ( cpO3 < 0.0  && recocpO3_etavariup < 0.0 )
                     {
                        FillHisto( h_CPO3_EtaVariUp_Minus_Minus,recocpO3_etavariup,evt_weight_);
                     }
                  
                     ////////////////////////////////////////////////////////////// 
                     /// CPO3 with Eta Variation (Down) at Reconstruction Level ///
                     //////////////////////////////////////////////////////////////
                     if (recocpO3_etavaridown > 0.0) {
                        FillHisto( h_CPO3reco_evtweight_EtaVariDown,0.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_EtaVariDown,6.5,evt_weight_*evt_weight_);
                     }
                     else{
                        FillHisto( h_CPO3reco_evtweight_EtaVariDown,2.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_EtaVariDown,8.5,evt_weight_*evt_weight_);
                     }
                     if ( cpO3 > 0.0  && recocpO3_etavaridown > 0.0 )
                     {
                        FillHisto( h_CPO3_EtaVariDown_Plus_Plus,recocpO3_etavaridown,evt_weight_);
                     }
                     if ( cpO3 > 0.0  && recocpO3_etavaridown < 0.0 )
                     {
                        FillHisto( h_CPO3_EtaVariDown_Plus_Minus,recocpO3_etavaridown,evt_weight_);
                     }
                  
                     if ( cpO3 < 0.0  && recocpO3_etavaridown > 0.0 )
                     {
                        FillHisto( h_CPO3_EtaVariDown_Minus_Plus,recocpO3_etavaridown,evt_weight_);
                     }
                     if ( cpO3 < 0.0  && recocpO3_etavaridown < 0.0 )
                     {
                        FillHisto( h_CPO3_EtaVariDown_Minus_Minus,recocpO3_etavaridown,evt_weight_);
                     }
                  
                     //////////////////////////////////////////////////////////// 
                     /// CPO3 with Phi Variation (Up) at Reconstruction Level ///
                     //////////////////////////////////////////////////////////// 
                     if (recocpO3_phivariup > 0.0) {
                        FillHisto( h_CPO3reco_evtweight_PhiVariUp,0.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_PhiVariUp,6.5,evt_weight_*evt_weight_);
                     }
                     else{
                        FillHisto( h_CPO3reco_evtweight_PhiVariUp,2.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_PhiVariUp,8.5,evt_weight_*evt_weight_);
                     }
                     if ( recocpO3_phivariup > 0.0) {
                        FillHisto( h_CPO3reco_evtweight_PhiVariUp,0.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_PhiVariUp,6.5,evt_weight_*evt_weight_);
                     }
                     else{
                        FillHisto( h_CPO3reco_evtweight_PhiVariUp,2.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_PhiVariUp,8.5,evt_weight_*evt_weight_);
                     }
                  
                     if ( cpO3 > 0.0  && recocpO3_phivariup > 0.0 )
                     {
                        FillHisto( h_CPO3_PhiVariUp_Plus_Plus,recocpO3_phivariup,evt_weight_);
                     }
                     if ( cpO3 > 0.0  && recocpO3_phivariup < 0.0 )
                     {
                        FillHisto( h_CPO3_PhiVariUp_Plus_Minus,recocpO3_phivariup,evt_weight_);
                     }
                  
                     if ( cpO3 < 0.0  && recocpO3_phivariup > 0.0 )
                     {
                        FillHisto( h_CPO3_PhiVariUp_Minus_Plus,recocpO3_phivariup,evt_weight_);
                     }
                     if ( cpO3 < 0.0  && recocpO3_phivariup < 0.0 )
                     {
                        FillHisto( h_CPO3_PhiVariUp_Minus_Minus,recocpO3_phivariup,evt_weight_);
                     }
                     ////////////////////////////////////////////////////////////// 
                     /// CPO3 with Phi Variation (Down) at Reconstruction Level ///
                     //////////////////////////////////////////////////////////////
                  
                     if (recocpO3_phivaridown > 0.0) {
                        FillHisto( h_CPO3reco_evtweight_PhiVariDown,0.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_PhiVariDown,6.5,evt_weight_*evt_weight_);
                     }
                     else{
                        FillHisto( h_CPO3reco_evtweight_PhiVariDown,2.5,evt_weight_);
                        FillHisto( h_CPO3reco_evtweight_PhiVariDown,8.5,evt_weight_*evt_weight_);
                     }
                     if ( cpO3 > 0.0  && recocpO3_phivaridown > 0.0 )
                     {
                        FillHisto( h_CPO3_PhiVariDown_Plus_Plus,recocpO3_phivaridown,evt_weight_);
                     }
                     if ( cpO3 > 0.0  && recocpO3_phivaridown < 0.0 )
                     {
                        FillHisto( h_CPO3_PhiVariDown_Plus_Minus,recocpO3_phivaridown,evt_weight_);
                     }
                  
                     if ( cpO3 < 0.0  && recocpO3_phivaridown > 0.0 )
                     {
                        FillHisto( h_CPO3_PhiVariDown_Minus_Plus,recocpO3_phivaridown,evt_weight_);
                     }
                     if ( cpO3 < 0.0  && recocpO3_phivaridown < 0.0 )
                     {
                        FillHisto( h_CPO3_PhiVariDown_Minus_Minus,recocpO3_phivaridown,evt_weight_);
                     }
                  } // KinSolver // for central !!
               } // bJetcut 1- or More //
            }  // MET Cut 
         } // Num Jet Cut 

         ///////////////////////////////
         /// JES Up Systematic Study ///
         ///////////////////////////////
         if ( NumJetCut(v_jetup_idx) == true && isAllSyst == true ) 
         {
            FillHisto( h_cf_sys_NLeptons[idx_jecup][3], v_lepton_idx.size(), v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Lep1pt[idx_jecup][3] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Lep1eta[idx_jecup][3], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Lep1phi[idx_jecup][3], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Lep2pt[idx_jecup][3] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Lep2eta[idx_jecup][3], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Lep2phi[idx_jecup][3], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_NPV[idx_jecup][3]    , num_pv        , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_NJets[idx_jecup][3]  , v_jetup_idx.size(), v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Jet1pt[idx_jecup][3] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Jet1eta[idx_jecup][3], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Jet1phi[idx_jecup][3], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Jet2pt[idx_jecup][3] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Jet2eta[idx_jecup][3], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_Jet2phi[idx_jecup][3], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_metpt[idx_jecup][3] ,MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_metphi[idx_jecup][3],MetJESUp->Phi() , v_SystEvt[idx_jecup] );
            
            FillHisto( h_cf_sys_dilep_inv_mass[idx_jecup][3], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_metpt[idx_jecup][3] , MetJESUp->Pt() , v_SystEvt[idx_jecup] );
            FillHisto( h_cf_sys_metphi[idx_jecup][3], MetJESUp->Phi(), v_SystEvt[idx_jecup] );
            
            FillHisto( h_sys_Num_PV[idx_jecup][3], num_pv, v_SystEvt[idx_jecup] );
            FillHisto( h_sys_DiLepMass[idx_jecup][3], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Lep1pt[idx_jecup][3] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Lep1eta[idx_jecup][3], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Lep1phi[idx_jecup][3], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Lep2pt[idx_jecup][3] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Lep2eta[idx_jecup][3], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Lep2phi[idx_jecup][3], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Jet1pt[idx_jecup][3] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Jet1eta[idx_jecup][3], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Jet1phi[idx_jecup][3], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Jet2pt[idx_jecup][3] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Jet2eta[idx_jecup][3], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Jet2phi[idx_jecup][3], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_Num_Jets[idx_jecup][3], v_jetup_idx.size(), v_SystEvt[idx_jecup] );
            FillHisto( h_sys_METpt[idx_jecup][3]   , MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
            FillHisto( h_sys_METphi[idx_jecup][3]  , MetJESUp->Phi()  , v_SystEvt[idx_jecup] );
            if ( METCut(MetJESUp) == true ) 
            {
               FillHisto( h_cf_sys_NLeptons[idx_jecup][4], v_lepton_idx.size(), v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Lep1pt[idx_jecup][4] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Lep1eta[idx_jecup][4], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Lep1phi[idx_jecup][4], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Lep2pt[idx_jecup][4] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Lep2eta[idx_jecup][4], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Lep2phi[idx_jecup][4], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_NPV[idx_jecup][4]    , num_pv        , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_NJets[idx_jecup][4]  , v_jet_idx.size(), v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Jet1pt[idx_jecup][4] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Jet1eta[idx_jecup][4], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Jet1phi[idx_jecup][4], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Jet2pt[idx_jecup][4] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Jet2eta[idx_jecup][4], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_Jet2phi[idx_jecup][4], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_metpt[idx_jecup][4] ,MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_metphi[idx_jecup][4],MetJESUp->Phi() , v_SystEvt[idx_jecup] );
               
               FillHisto( h_cf_sys_dilep_inv_mass[idx_jecup][4], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_metpt[idx_jecup][4] , Met->Pt() , v_SystEvt[idx_jecup] );
               FillHisto( h_cf_sys_metphi[idx_jecup][4], Met->Phi(), v_SystEvt[idx_jecup] );
               
               FillHisto( h_sys_Num_PV[idx_jecup][4], num_pv, v_SystEvt[idx_jecup] );
               FillHisto( h_sys_DiLepMass[idx_jecup][4], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Lep1pt[idx_jecup][4] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Lep1eta[idx_jecup][4], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Lep1phi[idx_jecup][4], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Lep2pt[idx_jecup][4] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Lep2eta[idx_jecup][4], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Lep2phi[idx_jecup][4], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Jet1pt[idx_jecup][4] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Jet1eta[idx_jecup][4], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Jet1phi[idx_jecup][4], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Jet2pt[idx_jecup][4] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Jet2eta[idx_jecup][4], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Jet2phi[idx_jecup][4], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_Num_Jets[idx_jecup][4], v_jetup_idx.size(), v_SystEvt[idx_jecup] );
               FillHisto( h_sys_METpt[idx_jecup][4]   , MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
               FillHisto( h_sys_METphi[idx_jecup][4]  , MetJESUp->Phi()  , v_SystEvt[idx_jecup] );
                  /*cout << "Info_EventNumber : " << Info_EventNumber << 
                           " Jet1Up : " << Jet1Up->Pt() << " Jet2Up : " << Jet2Up->Pt() << 
                           " Jet1 : " << Jet1->Pt() << " Jet2 : " << Jet2->Pt() << 
                            " evt_: " << v_SystEvt[i] << endl;*/

               /// Apply BTagging Scale Factor for JES Up
               BTaggigSFApplyJES("JetEnUp");
               if ( BJetCut(v_bjetup_idx) == true ) // one or more b-tagging 
               {
                  FillHisto( h_cf_sys_NLeptons[idx_jecup][5], v_lepton_idx.size(), v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Lep1pt[idx_jecup][5] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Lep1eta[idx_jecup][5], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Lep1phi[idx_jecup][5], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Lep2pt[idx_jecup][5] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Lep2eta[idx_jecup][5], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Lep2phi[idx_jecup][5], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_NPV[idx_jecup][5]    , num_pv        , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_NJets[idx_jecup][5]  , v_jetup_idx.size(), v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Jet1pt[idx_jecup][5] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Jet1eta[idx_jecup][5], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Jet1phi[idx_jecup][5], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Jet2pt[idx_jecup][5] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Jet2eta[idx_jecup][5], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_Jet2phi[idx_jecup][5], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_metpt[idx_jecup][5] ,MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_metphi[idx_jecup][5],MetJESUp->Phi() , v_SystEvt[idx_jecup] );
                  
                  FillHisto( h_cf_sys_dilep_inv_mass[idx_jecup][5], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_metpt[idx_jecup][5] , MetJESUp->Pt() , v_SystEvt[idx_jecup] );
                  FillHisto( h_cf_sys_metphi[idx_jecup][5], MetJESUp->Phi(), v_SystEvt[idx_jecup] );
                  
                  FillHisto( h_sys_Num_PV[idx_jecup][5], num_pv, v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_DiLepMass[idx_jecup][5], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Lep1pt[idx_jecup][5] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Lep1eta[idx_jecup][5], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Lep1phi[idx_jecup][5], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Lep2pt[idx_jecup][5] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Lep2eta[idx_jecup][5], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Lep2phi[idx_jecup][5], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Jet1pt[idx_jecup][5] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Jet1eta[idx_jecup][5], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Jet1phi[idx_jecup][5], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Jet2pt[idx_jecup][5] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Jet2eta[idx_jecup][5], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Jet2phi[idx_jecup][5], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_Num_Jets[idx_jecup][5], v_jetup_idx.size(), v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_METpt[idx_jecup][5]   , MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
                  FillHisto( h_sys_METphi[idx_jecup][5]  , MetJESUp->Phi()  , v_SystEvt[idx_jecup] );

                  ///////////////////////////////////////
                  /// 2 or More B-Tagging Requriement ///
                  /////////////////////////////////////// 
                  if ( DoubleBtag(v_bjetup_idx) == true )
                  {
                     FillHisto( h_cf_sys_NLeptons[idx_jecup][6], v_lepton_idx.size(), v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep1pt[idx_jecup][6] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep1eta[idx_jecup][6], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep1phi[idx_jecup][6], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep2pt[idx_jecup][6] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep2eta[idx_jecup][6], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep2phi[idx_jecup][6], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_NPV[idx_jecup][6]    , num_pv        , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_NJets[idx_jecup][6]  , v_jet_idx.size(), v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet1pt[idx_jecup][6] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet1eta[idx_jecup][6], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet1phi[idx_jecup][6], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet2pt[idx_jecup][6] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet2eta[idx_jecup][6], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet2phi[idx_jecup][6], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metpt[idx_jecup][6] ,MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metphi[idx_jecup][6],MetJESUp->Phi() , v_SystEvt[idx_jecup] );
                     
                     FillHisto( h_cf_sys_dilep_inv_mass[idx_jecup][6], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metpt[idx_jecup][6] , MetJESUp->Pt() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metphi[idx_jecup][6], MetJESUp->Phi(), v_SystEvt[idx_jecup] );
                     
                     FillHisto( h_sys_Num_PV[idx_jecup][6], num_pv, v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_DiLepMass[idx_jecup][6], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep1pt[idx_jecup][6] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep1eta[idx_jecup][6], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep1phi[idx_jecup][6], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep2pt[idx_jecup][6] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep2eta[idx_jecup][6], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep2phi[idx_jecup][6], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet1pt[idx_jecup][6] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet1eta[idx_jecup][6], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet1phi[idx_jecup][6], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet2pt[idx_jecup][6] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet2eta[idx_jecup][6], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet2phi[idx_jecup][6], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Num_Jets[idx_jecup][6], v_jetup_idx.size(), v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_METpt[idx_jecup][6]   , MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_METphi[idx_jecup][6]  , MetJESUp->Phi()  , v_SystEvt[idx_jecup] );

                  }//// JES Up  2 or more b-tagging  
                  ////
                  
                  if (v_bjetup_idx.size() ==2)// exactly 2b-tagging ..
                  {
                     FillHisto( h_cf_sys_NLeptons[idx_jecup][7], v_lepton_idx.size(), v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep1pt[idx_jecup][7] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep1eta[idx_jecup][7], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep1phi[idx_jecup][7], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep2pt[idx_jecup][7] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep2eta[idx_jecup][7], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep2phi[idx_jecup][7], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_NPV[idx_jecup][7]    , num_pv        , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_NJets[idx_jecup][7]  , v_jet_idx.size(), v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet1pt[idx_jecup][7] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet1eta[idx_jecup][7], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet1phi[idx_jecup][7], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet2pt[idx_jecup][7] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet2eta[idx_jecup][7], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet2phi[idx_jecup][7], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metpt[idx_jecup][7] ,MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metphi[idx_jecup][7],MetJESUp->Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_dilep_inv_mass[idx_jecup][7], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metpt[idx_jecup][7] , MetJESUp->Pt() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metphi[idx_jecup][7], MetJESUp->Phi(), v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Num_PV[idx_jecup][7], num_pv, v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_DiLepMass[idx_jecup][7], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep1pt[idx_jecup][7] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep1eta[idx_jecup][7], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep1phi[idx_jecup][7], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep2pt[idx_jecup][7] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep2eta[idx_jecup][7], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep2phi[idx_jecup][7], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet1pt[idx_jecup][7] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet1eta[idx_jecup][7], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet1phi[idx_jecup][7], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet2pt[idx_jecup][7] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet2eta[idx_jecup][7], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet2phi[idx_jecup][7], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Num_Jets[idx_jecup][7], v_jetup_idx.size(), v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_METpt[idx_jecup][7]   , MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_METphi[idx_jecup][7]  , MetJESUp->Phi()  , v_SystEvt[idx_jecup] );
                  }// JES Up  exactly 2b-tagging ..
                  /// Top Reconstruction for JES Up 
                  KinSolSys(Lep,AnLep,bJet1Up,bJet2Up,MetJESUp);
                  if ( ksolweight_ != -1 )
                  {
                     if ( Top->Pt() > AnTop->Pt() ) { (*Top1) = (*Top); (*Top2) = (*AnTop); }
                     else { (*Top1) = (*AnTop); (*Top2) = (*Top); }
                     FillHisto( h_cf_sys_NLeptons[idx_jecup][8], v_lepton_idx.size(), v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep1pt[idx_jecup][8] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep1eta[idx_jecup][8], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep1phi[idx_jecup][8], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep2pt[idx_jecup][8] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep2eta[idx_jecup][8], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Lep2phi[idx_jecup][8], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_NPV[idx_jecup][8]    , num_pv        , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_NJets[idx_jecup][8]  , v_jet_idx.size(), v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet1pt[idx_jecup][8] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet1eta[idx_jecup][8], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet1phi[idx_jecup][8], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet2pt[idx_jecup][8] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet2eta[idx_jecup][8], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_Jet2phi[idx_jecup][8], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metpt[idx_jecup][8] ,MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metphi[idx_jecup][8],MetJESUp->Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_dilep_inv_mass[idx_jecup][8], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metpt[idx_jecup][8] , MetJESUp->Pt() , v_SystEvt[idx_jecup] );
                     FillHisto( h_cf_sys_metphi[idx_jecup][8], MetJESUp->Phi(), v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Num_PV[idx_jecup][8], num_pv, v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_DiLepMass[idx_jecup][8], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep1pt[idx_jecup][8] , (*Lep1).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep1eta[idx_jecup][8], (*Lep1).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep1phi[idx_jecup][8], (*Lep1).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep2pt[idx_jecup][8] , (*Lep2).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep2eta[idx_jecup][8], (*Lep2).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Lep2phi[idx_jecup][8], (*Lep2).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet1pt[idx_jecup][8] , (*Jet1Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet1eta[idx_jecup][8], (*Jet1Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet1phi[idx_jecup][8], (*Jet1Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet2pt[idx_jecup][8] , (*Jet2Up).Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet2eta[idx_jecup][8], (*Jet2Up).Eta() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Jet2phi[idx_jecup][8], (*Jet2Up).Phi() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Num_Jets[idx_jecup][8], v_jetup_idx.size(), v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_METpt[idx_jecup][8]   , MetJESUp->Pt()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_METphi[idx_jecup][8]  , MetJESUp->Phi()  , v_SystEvt[idx_jecup] );

                     FillHisto( h_sys_Top1Mass_[idx_jecup]    , Top1->M()        , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Top1pt_[idx_jecup]      , Top1->Pt()       , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Top1phi_[idx_jecup]     , Top1->Phi()      , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Top1Rapidity_[idx_jecup], Top1->Rapidity() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Top1Energy_[idx_jecup]  , Top1->Energy()   , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Top2Mass_[idx_jecup]    , Top2->M()        , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Top2pt_[idx_jecup]      , Top2->Pt()       , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Top2phi_[idx_jecup]     , Top2->Phi()      , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Top2Rapidity_[idx_jecup], Top2->Rapidity() , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Top2Energy_[idx_jecup]  , Top2->Energy()   , v_SystEvt[idx_jecup] );

                     FillHisto( h_sys_TopMass_[idx_jecup]      , Top->M()         , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Toppt_[idx_jecup]        , Top->Pt()        , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_Topphi_[idx_jecup]       , Top->Phi()       , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_TopRapidity_[idx_jecup]  , Top->Rapidity()  , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_TopEnergy_[idx_jecup]    , Top->Energy()    , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_AnTopMass_[idx_jecup]    , AnTop->M()       , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_AnToppt_[idx_jecup]      , AnTop->Pt()      , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_AnTopphi_[idx_jecup]     , AnTop->Phi()     , v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_AnTopRapidity_[idx_jecup], AnTop->Rapidity(), v_SystEvt[idx_jecup] );
                     FillHisto( h_sys_AnTopEnergy_[idx_jecup]  , AnTop->Energy()  , v_SystEvt[idx_jecup] );



                     v_recocp_O.clear();
                     v_recocp_O.push_back( ssbcpviol->getO1Vari( Top, AnTop, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO2Vari( Top, AnTop, bJet, AnbJet ) );
                     v_recocp_O.push_back( ssbcpviol->getO3Vari( bJet, AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO4Vari( AnbJet, bJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO5Vari( bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO6Vari( bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO7Vari( Top , AnTop, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO8Vari( Top, AnTop, bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO9Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO10Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO11Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO12Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO13Vari( bJet , AnbJet, AnLep, Lep )  );
                     
                     for (int j = 0; j < v_recocp_O.size(); ++ j )
                     {
                        FillHisto( h_sys_Reco_CPO_[idx_jecup][j], v_recocp_O[j] , v_SystEvt[idx_jecup]   );
                        FillHisto( h_sys_Reco_CPO_ReRange_[idx_jecup][j], v_recocp_O[j] , v_SystEvt[idx_jecup]  );
                     }
                  } // Kinmatic Solver //
               }// JES Up  one or more b-tagging 
            }// JES Up MET Cut   
         } // JES Up Num Jet Cut
         if ( NumJetCut(v_jetdn_idx) == true && isAllSyst == true ) 
         {
            FillHisto( h_cf_sys_NLeptons[idx_jecdn][3], v_lepton_idx.size(), v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Lep1pt[idx_jecdn][3] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Lep1eta[idx_jecdn][3], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Lep1phi[idx_jecdn][3], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Lep2pt[idx_jecdn][3] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Lep2eta[idx_jecdn][3], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Lep2phi[idx_jecdn][3], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_NPV[idx_jecdn][3]    , num_pv        , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_NJets[idx_jecdn][3]  , v_jetup_idx.size(), v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Jet1pt[idx_jecdn][3] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Jet1eta[idx_jecdn][3], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Jet1phi[idx_jecdn][3], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Jet2pt[idx_jecdn][3] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Jet2eta[idx_jecdn][3], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_Jet2phi[idx_jecdn][3], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_metpt[idx_jecdn][3] ,MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_metphi[idx_jecdn][3],MetJESDn->Phi() , v_SystEvt[idx_jecdn] );
            
            FillHisto( h_cf_sys_dilep_inv_mass[idx_jecdn][3], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_metpt[idx_jecdn][3] , MetJESDn->Pt() , v_SystEvt[idx_jecdn] );
            FillHisto( h_cf_sys_metphi[idx_jecdn][3], MetJESDn->Phi(), v_SystEvt[idx_jecdn] );
            
            FillHisto( h_sys_Num_PV[idx_jecdn][3], num_pv, v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_DiLepMass[idx_jecdn][3], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Lep1pt[idx_jecdn][3] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Lep1eta[idx_jecdn][3], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Lep1phi[idx_jecdn][3], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Lep2pt[idx_jecdn][3] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Lep2eta[idx_jecdn][3], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Lep2phi[idx_jecdn][3], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Jet1pt[idx_jecdn][3] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Jet1eta[idx_jecdn][3], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Jet1phi[idx_jecdn][3], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Jet2pt[idx_jecdn][3] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Jet2eta[idx_jecdn][3], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Jet2phi[idx_jecdn][3], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_Num_Jets[idx_jecdn][3], v_jetup_idx.size(), v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_METpt[idx_jecdn][3]   , MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
            FillHisto( h_sys_METphi[idx_jecdn][3]  , MetJESDn->Phi()  , v_SystEvt[idx_jecdn] );

            if ( METCut(MetJESDn) == true )
            {
               FillHisto( h_cf_sys_NLeptons[idx_jecdn][4], v_lepton_idx.size(), v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Lep1pt[idx_jecdn][4] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Lep1eta[idx_jecdn][4], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Lep1phi[idx_jecdn][4], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Lep2pt[idx_jecdn][4] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Lep2eta[idx_jecdn][4], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Lep2phi[idx_jecdn][4], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_NPV[idx_jecdn][4]    , num_pv        , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_NJets[idx_jecdn][4]  , v_jet_idx.size(), v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Jet1pt[idx_jecdn][4] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Jet1eta[idx_jecdn][4], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Jet1phi[idx_jecdn][4], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Jet2pt[idx_jecdn][4] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Jet2eta[idx_jecdn][4], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_Jet2phi[idx_jecdn][4], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_metpt[idx_jecdn][4] ,MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_metphi[idx_jecdn][4],MetJESDn->Phi() , v_SystEvt[idx_jecdn] );
               
               FillHisto( h_cf_sys_dilep_inv_mass[idx_jecdn][4], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_metpt[idx_jecdn][4] , Met->Pt() , v_SystEvt[idx_jecdn] );
               FillHisto( h_cf_sys_metphi[idx_jecdn][4], Met->Phi(), v_SystEvt[idx_jecdn] );
               
               FillHisto( h_sys_Num_PV[idx_jecdn][4], num_pv, v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_DiLepMass[idx_jecdn][4], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Lep1pt[idx_jecdn][4] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Lep1eta[idx_jecdn][4], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Lep1phi[idx_jecdn][4], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Lep2pt[idx_jecdn][4] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Lep2eta[idx_jecdn][4], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Lep2phi[idx_jecdn][4], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Jet1pt[idx_jecdn][4] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Jet1eta[idx_jecdn][4], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Jet1phi[idx_jecdn][4], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Jet2pt[idx_jecdn][4] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Jet2eta[idx_jecdn][4], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Jet2phi[idx_jecdn][4], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_Num_Jets[idx_jecdn][4], v_jetup_idx.size(), v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_METpt[idx_jecdn][4]   , MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
               FillHisto( h_sys_METphi[idx_jecdn][4]  , MetJESDn->Phi()  , v_SystEvt[idx_jecdn] );
               /// Apply BTagging Scale Factor for JES Dn
               BTaggigSFApplyJES("JetEnDown");
               if ( BJetCut(v_bjetdn_idx) == true ) // one or more b-tagging 
               {
                  FillHisto( h_cf_sys_NLeptons[idx_jecdn][5], v_lepton_idx.size(), v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Lep1pt[idx_jecdn][5] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Lep1eta[idx_jecdn][5], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Lep1phi[idx_jecdn][5], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Lep2pt[idx_jecdn][5] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Lep2eta[idx_jecdn][5], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Lep2phi[idx_jecdn][5], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_NPV[idx_jecdn][5]    , num_pv        , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_NJets[idx_jecdn][5]  , v_jetup_idx.size(), v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Jet1pt[idx_jecdn][5] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Jet1eta[idx_jecdn][5], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Jet1phi[idx_jecdn][5], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Jet2pt[idx_jecdn][5] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Jet2eta[idx_jecdn][5], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_Jet2phi[idx_jecdn][5], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_metpt[idx_jecdn][5] ,MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_metphi[idx_jecdn][5],MetJESDn->Phi() , v_SystEvt[idx_jecdn] );
                  
                  FillHisto( h_cf_sys_dilep_inv_mass[idx_jecdn][5], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_metpt[idx_jecdn][5] , MetJESDn->Pt() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_cf_sys_metphi[idx_jecdn][5], MetJESDn->Phi(), v_SystEvt[idx_jecdn] );
                  
                  FillHisto( h_sys_Num_PV[idx_jecdn][5], num_pv, v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_DiLepMass[idx_jecdn][5], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Lep1pt[idx_jecdn][5] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Lep1eta[idx_jecdn][5], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Lep1phi[idx_jecdn][5], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Lep2pt[idx_jecdn][5] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Lep2eta[idx_jecdn][5], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Lep2phi[idx_jecdn][5], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Jet1pt[idx_jecdn][5] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Jet1eta[idx_jecdn][5], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Jet1phi[idx_jecdn][5], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Jet2pt[idx_jecdn][5] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Jet2eta[idx_jecdn][5], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Jet2phi[idx_jecdn][5], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_Num_Jets[idx_jecdn][5], v_jetdn_idx.size(), v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_METpt[idx_jecdn][5]   , MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
                  FillHisto( h_sys_METphi[idx_jecdn][5]  , MetJESDn->Phi()  , v_SystEvt[idx_jecdn] );
                  ///////////////////////////////////////
                  /// 2 or More B-Tagging Requriement ///
                  /////////////////////////////////////// 
                  if ( DoubleBtag(v_bjetdn_idx) == true )
                  {
                     FillHisto( h_cf_sys_NLeptons[idx_jecdn][6], v_lepton_idx.size(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep1pt[idx_jecdn][6] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep1eta[idx_jecdn][6], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep1phi[idx_jecdn][6], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep2pt[idx_jecdn][6] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep2eta[idx_jecdn][6], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep2phi[idx_jecdn][6], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_NPV[idx_jecdn][6]    , num_pv        , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_NJets[idx_jecdn][6]  , v_jet_idx.size(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet1pt[idx_jecdn][6] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet1eta[idx_jecdn][6], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet1phi[idx_jecdn][6], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet2pt[idx_jecdn][6] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet2eta[idx_jecdn][6], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet2phi[idx_jecdn][6], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metpt[idx_jecdn][6] ,MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metphi[idx_jecdn][6],MetJESDn->Phi() , v_SystEvt[idx_jecdn] );
                     
                     FillHisto( h_cf_sys_dilep_inv_mass[idx_jecdn][6], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metpt[idx_jecdn][6] , MetJESDn->Pt() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metphi[idx_jecdn][6], MetJESDn->Phi(), v_SystEvt[idx_jecdn] );
                     
                     FillHisto( h_sys_Num_PV[idx_jecdn][6], num_pv, v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_DiLepMass[idx_jecdn][6], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep1pt[idx_jecdn][6] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep1eta[idx_jecdn][6], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep1phi[idx_jecdn][6], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep2pt[idx_jecdn][6] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep2eta[idx_jecdn][6], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep2phi[idx_jecdn][6], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet1pt[idx_jecdn][6] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet1eta[idx_jecdn][6], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet1phi[idx_jecdn][6], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet2pt[idx_jecdn][6] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet2eta[idx_jecdn][6], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet2phi[idx_jecdn][6], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Num_Jets[idx_jecdn][6], v_jetdn_idx.size(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_METpt[idx_jecdn][6]   , MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_METphi[idx_jecdn][6]  , MetJESDn->Phi()  , v_SystEvt[idx_jecdn] );

                  }//// JES Down  2 or more b-tagging  
                  if (v_bjetdn_idx.size() ==2)// exactly 2b-tagging ..
                  {
                     FillHisto( h_cf_sys_NLeptons[idx_jecdn][7], v_lepton_idx.size(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep1pt[idx_jecdn][7] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep1eta[idx_jecdn][7], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep1phi[idx_jecdn][7], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep2pt[idx_jecdn][7] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep2eta[idx_jecdn][7], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep2phi[idx_jecdn][7], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_NPV[idx_jecdn][7]    , num_pv        , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_NJets[idx_jecdn][7]  , v_jet_idx.size(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet1pt[idx_jecdn][7] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet1eta[idx_jecdn][7], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet1phi[idx_jecdn][7], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet2pt[idx_jecdn][7] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet2eta[idx_jecdn][7], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet2phi[idx_jecdn][7], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metpt[idx_jecdn][7] ,MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metphi[idx_jecdn][7],MetJESDn->Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_dilep_inv_mass[idx_jecdn][7], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metpt[idx_jecdn][7] , MetJESDn->Pt() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metphi[idx_jecdn][7], MetJESDn->Phi(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Num_PV[idx_jecdn][7], num_pv, v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_DiLepMass[idx_jecdn][7], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep1pt[idx_jecdn][7] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep1eta[idx_jecdn][7], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep1phi[idx_jecdn][7], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep2pt[idx_jecdn][7] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep2eta[idx_jecdn][7], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep2phi[idx_jecdn][7], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet1pt[idx_jecdn][7] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet1eta[idx_jecdn][7], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet1phi[idx_jecdn][7], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet2pt[idx_jecdn][7] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet2eta[idx_jecdn][7], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet2phi[idx_jecdn][7], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Num_Jets[idx_jecdn][7], v_jetdn_idx.size(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_METpt[idx_jecdn][7]   , MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_METphi[idx_jecdn][7]  , MetJESDn->Phi()  , v_SystEvt[idx_jecdn] );
                  }// JES Dn  exactly 2b-tagging ..
               
                  /// Top Reconstruction for JES Down 
                  KinSolSys(Lep,AnLep,bJet1Dn,bJet2Dn,MetJESDn);
                  if ( ksolweight_ != -1 )
                  {
                     if ( Top->Pt() > AnTop->Pt() ) { (*Top1) = (*Top); (*Top2) = (*AnTop); }
                     else { (*Top1) = (*AnTop); (*Top2) = (*Top); }
                     FillHisto( h_cf_sys_NLeptons[idx_jecdn][8], v_lepton_idx.size(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep1pt[idx_jecdn][8] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep1eta[idx_jecdn][8], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep1phi[idx_jecdn][8], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep2pt[idx_jecdn][8] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep2eta[idx_jecdn][8], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Lep2phi[idx_jecdn][8], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_NPV[idx_jecdn][8]    , num_pv        , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_NJets[idx_jecdn][8]  , v_jet_idx.size(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet1pt[idx_jecdn][8] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet1eta[idx_jecdn][8], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet1phi[idx_jecdn][8], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet2pt[idx_jecdn][8] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet2eta[idx_jecdn][8], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_Jet2phi[idx_jecdn][8], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metpt[idx_jecdn][8] ,MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metphi[idx_jecdn][8],MetJESDn->Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_dilep_inv_mass[idx_jecdn][8], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metpt[idx_jecdn][8] , MetJESDn->Pt() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_cf_sys_metphi[idx_jecdn][8], MetJESDn->Phi(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Num_PV[idx_jecdn][8], num_pv, v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_DiLepMass[idx_jecdn][8], ( (*Lep1)+(*Lep2) ).M(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep1pt[idx_jecdn][8] , (*Lep1).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep1eta[idx_jecdn][8], (*Lep1).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep1phi[idx_jecdn][8], (*Lep1).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep2pt[idx_jecdn][8] , (*Lep2).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep2eta[idx_jecdn][8], (*Lep2).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Lep2phi[idx_jecdn][8], (*Lep2).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet1pt[idx_jecdn][8] , (*Jet1Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet1eta[idx_jecdn][8], (*Jet1Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet1phi[idx_jecdn][8], (*Jet1Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet2pt[idx_jecdn][8] , (*Jet2Dn).Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet2eta[idx_jecdn][8], (*Jet2Dn).Eta() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Jet2phi[idx_jecdn][8], (*Jet2Dn).Phi() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Num_Jets[idx_jecdn][8], v_jetup_idx.size(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_METpt[idx_jecdn][8]   , MetJESDn->Pt()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_METphi[idx_jecdn][8]  , MetJESDn->Phi()  , v_SystEvt[idx_jecdn] );

                     FillHisto( h_sys_Top1Mass_[idx_jecdn]    , Top1->M()        , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Top1pt_[idx_jecdn]      , Top1->Pt()       , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Top1phi_[idx_jecdn]     , Top1->Phi()      , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Top1Rapidity_[idx_jecdn], Top1->Rapidity() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Top1Energy_[idx_jecdn]  , Top1->Energy()   , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Top2Mass_[idx_jecdn]    , Top2->M()        , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Top2pt_[idx_jecdn]      , Top2->Pt()       , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Top2phi_[idx_jecdn]     , Top2->Phi()      , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Top2Rapidity_[idx_jecdn], Top2->Rapidity() , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Top2Energy_[idx_jecdn]  , Top2->Energy()   , v_SystEvt[idx_jecdn] );

                     FillHisto( h_sys_TopMass_[idx_jecdn]      , Top->M()         , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Toppt_[idx_jecdn]        , Top->Pt()        , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_Topphi_[idx_jecdn]       , Top->Phi()       , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_TopRapidity_[idx_jecdn]  , Top->Rapidity()  , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_TopEnergy_[idx_jecdn]    , Top->Energy()    , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_AnTopMass_[idx_jecdn]    , AnTop->M()       , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_AnToppt_[idx_jecdn]      , AnTop->Pt()      , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_AnTopphi_[idx_jecdn]     , AnTop->Phi()     , v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_AnTopRapidity_[idx_jecdn], AnTop->Rapidity(), v_SystEvt[idx_jecdn] );
                     FillHisto( h_sys_AnTopEnergy_[idx_jecdn]  , AnTop->Energy()  , v_SystEvt[idx_jecdn] );


                     v_recocp_O.clear();
                     v_recocp_O.push_back( ssbcpviol->getO1Vari( Top, AnTop, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO2Vari( Top, AnTop, bJet, AnbJet ) );
                     v_recocp_O.push_back( ssbcpviol->getO3Vari( bJet, AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO4Vari( AnbJet, bJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO5Vari( bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO6Vari( bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO7Vari( Top , AnTop, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO8Vari( Top, AnTop, bJet , AnbJet, AnLep, Lep ) );
                     v_recocp_O.push_back( ssbcpviol->getO9Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO10Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO11Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO12Vari( bJet , AnbJet, AnLep, Lep )  );
                     v_recocp_O.push_back( ssbcpviol->getO13Vari( bJet , AnbJet, AnLep, Lep )  );

                     for (int j = 0; j < v_recocp_O.size(); ++ j )
                     {
                        FillHisto( h_sys_Reco_CPO_[idx_jecdn][j], v_recocp_O[j] , v_SystEvt[idx_jecdn]   );
                        FillHisto( h_sys_Reco_CPO_ReRange_[idx_jecdn][j], v_recocp_O[j] , v_SystEvt[idx_jecdn]  );
                     }
                  } // Kinmatic Solver //
               }// Num bJet Cut // one or more b-tagging //
            } // JES Down MET Cut // 
         } // JES Down Num Jet Cut 
         //////////////////////
         /// Fill Histogram ///
         //////////////////////
         //ClearVectors();*/
     
      }//Di-Lepton Analysis//

      ///////////////////////////////////
      /// ** Lepton + Jet Analysis ** ///
      ///////////////////////////////////
      else if ( TString(Decaymode).Contains( "muonJet" )  ) 
      {
         cout << " -- We don't have Muon + Jet code ---" << endl;   
      }
      else {cout << "Decaymode Error : " << Decaymode << endl; }
   }//event loop
  
   printf("Total processed number of events: %lld\n", __tot_evt);


}//end Loop function

void ssb_analysis::GetNtupleTotalEvent( unsigned int totevent )// Not Used Function
{
   NtupletotalEvent = totevent;
}

void ssb_analysis::Start( int genLoopon )
{
   if      ( genLoopon == 0 ){ fout = new TFile(Form("output/%s",outfile),"RECREATE");}
   else if      ( genLoopon == 1 ){ fout = new TFile(Form("output/%s",outfile),"UPDATE");}
   //if      ( genLoopon == 0 ){ fout = new TFile(Form("gsidcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/SSB_CPviolation/output/%s",outfile),"RECREATE");}
   //else if ( genLoopon == 1 ){ fout = new TFile(Form("gsidcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sha/SSB_CPviolation/output/%s",outfile),"UPDATE"  );}
   else {cout << "genLoopon error" << endl;}
   fout->cd("");

   TDirectory *dir = gDirectory;
   dir->cd();

   DeclareHistos();

   /////////////////////////////////
   /// For All-in-One Systematic ///
   /////////////////////////////////
   if (isAllSyst == true)
   {
      for (int i = 0; i < v_SystFullName.size(); ++i )
      {
         gSystem->mkdir(Form("output/%s",v_SystFullName[i].Data()));
         a_fout[i] = new TFile(Form("output/%s/%s",v_SystFullName[i].Data(),v_outName[i].Data()),"RECREATE");
         a_fout[i]->cd();
         DeclareHistosSyst(i);
      }
   }
   ReadDuplList(); 
}

void ssb_analysis::DeclareHistos()
{

   /// Test For Systematic All-in-One Code ///

   for (int i =0 ; i < 10 ; i++)
   {
      h_cf_NLeptons[i]       = new TH1F(Form("_h_cf_NLeptons_%d_"      , i),Form("Num. Lepton %s"                ,cutflowName[i].Data()), 20  , 0    , 20  ); h_cf_NLeptons[i]->Sumw2();         
      h_cf_Lep1pt[i]         = new TH1F(Form("_h_cf_Lep1pt_%d_"        , i),Form("Leading Lepton pT %s"          ,cutflowName[i].Data()), 1000, 0    , 1000); h_cf_Lep1pt[i]->Sumw2();           
      h_cf_Lep1phi[i]        = new TH1F(Form("_h_cf_Lep1phi_%d_"       , i),Form("Leading Lepton phi %s"        ,cutflowName[i].Data()), 24  , -1*pi, 1*pi); h_cf_Lep1phi[i]->Sumw2();          
      h_cf_Lep1eta[i]        = new TH1F(Form("_h_cf_Lep1eta_%d_"       , i),Form("Leading Lepton eta %s"        ,cutflowName[i].Data()), 50  , -2.5 , 2.5 ); h_cf_Lep1eta[i]->Sumw2();          
      h_cf_Lep2pt[i]         = new TH1F(Form("_h_cf_Lep2pt_%d_"        , i),Form("Second Leading Lepton pT %s"   ,cutflowName[i].Data()), 1000, 0    , 1000); h_cf_Lep2pt[i]->Sumw2();           
      h_cf_Lep2phi[i]        = new TH1F(Form("_h_cf_Lep2phi_%d_"       , i),Form("Second Leading Lepton phi %s"  ,cutflowName[i].Data()), 24  , -1*pi, 1*pi); h_cf_Lep2phi[i]->Sumw2();      
      h_cf_Lep2eta[i]        = new TH1F(Form("_h_cf_Lep2eta_%d_"       , i),Form("Second Leading Lepton eta %s"  ,cutflowName[i].Data()), 50  , -2.5 , 2.5 ); h_cf_Lep2eta[i]->Sumw2();       
      h_cf_dilep_inv_mass[i] = new TH1F(Form("_h_cf_dilep_inv_mass_%d_", i),Form("Invariant Mass of Dilepton %s" ,cutflowName[i].Data()), 1000, 0    , 1000); h_cf_dilep_inv_mass[i]->Sumw2();
      h_cf_Jet1pt[i]         = new TH1F(Form("_h_cf_Jet1pt_%d_"        , i),Form("Leading Jet pT %s"             ,cutflowName[i].Data()), 1000, 0    , 1000); h_cf_Jet1pt[i]->Sumw2();       
      h_cf_Jet1phi[i]        = new TH1F(Form("_h_cf_Jet1phi_%d_"       , i),Form("Leading Jet phi %s"            ,cutflowName[i].Data()), 24  , -1*pi, 1*pi); h_cf_Jet1phi[i]->Sumw2();       
      h_cf_Jet1eta[i]        = new TH1F(Form("_h_cf_Jet1eta_%d_"       , i),Form("Leading Jet eta %s"            ,cutflowName[i].Data()), 50  , -2.5 , 2.5 ); h_cf_Jet1eta[i]->Sumw2();       
      h_cf_Jet2pt[i]         = new TH1F(Form("_h_cf_Jet2pt_%d_"        , i),Form("Second Leading Jet pT %s"      ,cutflowName[i].Data()), 1000, 0    , 1000); h_cf_Jet2pt[i]->Sumw2();        
      h_cf_Jet2phi[i]        = new TH1F(Form("_h_cf_Jet2phi_%d_"       , i),Form("Second Leading Jet phi %s"     ,cutflowName[i].Data()), 24  , -1*pi, 1*pi); h_cf_Jet2phi[i]->Sumw2();       
      h_cf_Jet2eta[i]        = new TH1F(Form("_h_cf_Jet2eta_%d_"       , i),Form("Second Leading Jet eta %s"     ,cutflowName[i].Data()), 50  , -2.5 , 2.5 ); h_cf_Jet2eta[i]->Sumw2();      
      h_cf_NJets[i]          = new TH1F(Form("_h_cf_NJets_%d_"         , i),Form("Num. of Jets %s"               ,cutflowName[i].Data()), 20  , 0    , 20  ); h_cf_NJets[i]->Sumw2();        
      h_cf_metpt[i]          = new TH1F(Form("_h_cf_metpt_%d_"         , i),Form("MET %s"                        ,cutflowName[i].Data()), 1000, 0    , 1000); h_cf_metpt[i]->Sumw2();        
      h_cf_metphi[i]         = new TH1F(Form("_h_cf_metphi_%d_"        , i),Form("MET Phi %s"                    ,cutflowName[i].Data()), 24  , -1*pi, 1*pi); h_cf_metphi[i]->Sumw2();       
      h_cf_Nbjets[i]         = new TH1F(Form("_h_cf_Nbjets_%d_"        , i),Form("Num. of b-jets %s"             ,cutflowName[i].Data()), 20  , 0    , 20  ); h_cf_Nbjets[i]->Sumw2();        
      h_cf_NPV[i]            = new TH1F(Form("_h_cf_NPV_%d_"           , i),Form("Num. of Primary Vertex %s"     ,cutflowName[i].Data()), 100 , 0    , 100 ); h_cf_NPV[i]->Sumw2();          


      h_EventWeight[i] = new TH1F(Form("h_EventWeight_%d",i), Form("Event Weight %s",cutflowName[i].Data()), 1000, -5, 5); h_EventWeight[i]->Sumw2();

      h_Lep1pt[i]  = new TH1F(Form("h_Lep1pt_%d" ,i), Form("Leading Lepton pT %s"        ,cutflowName[i].Data()), 250, 0.0, 250); h_Lep1pt[i]->Sumw2(); 
      h_Lep2pt[i]  = new TH1F(Form("h_Lep2pt_%d" ,i), Form("Second Leading Lepton pT %s" ,cutflowName[i].Data()), 250, 0.0, 250); h_Lep2pt[i]->Sumw2();
      h_Lep1eta[i] = new TH1F(Form("h_Lep1eta_%d",i), Form("Leading Lepton Eta    %s"    ,cutflowName[i].Data()), 50, -2.5, 2.5); h_Lep1eta[i]->Sumw2();
      h_Lep2eta[i] = new TH1F(Form("h_Lep2eta_%d",i), Form("Second Leading Lepton Eta %s",cutflowName[i].Data()), 50, -2.5, 2.5); h_Lep2eta[i]->Sumw2();
      h_Lep1phi[i] = new TH1F(Form("h_Lep1phi_%d",i), Form("Leading Lepton Phi %s"       ,cutflowName[i].Data()), 24, -1*pi, pi); h_Lep1phi[i]->Sumw2();
      h_Lep2phi[i] = new TH1F(Form("h_Lep2phi_%d",i), Form("Second Leading Lepton Phi %s",cutflowName[i].Data()), 24, -1*pi, pi); h_Lep2phi[i]->Sumw2();
   
      h_Muonpt[i]  = new TH1F(Form("h_Muonpt_%d",i ), Form("Muon pT %s"         ,cutflowName[i].Data()), 250, 0.0, 250); h_Muonpt[i]->Sumw2();
      h_Elecpt[i]  = new TH1F(Form("h_Elecpt_%d",i ), Form("Electron pT %s"     ,cutflowName[i].Data()), 250, 0.0, 250); h_Elecpt[i]->Sumw2();
      h_Muoneta[i] = new TH1F(Form("h_Muoneta_%d",i), Form("Muon Eta %s"        ,cutflowName[i].Data()), 50, -2.5, 2.5); h_Muoneta[i]->Sumw2();
      h_Eleceta[i] = new TH1F(Form("h_Eleceta_%d",i), Form("Electron Eta %s"    ,cutflowName[i].Data()), 50, -2.5, 2.5); h_Eleceta[i]->Sumw2();
      h_Muonphi[i] = new TH1F(Form("h_Muonphi_%d",i), Form("Muon Phi %s"        ,cutflowName[i].Data()), 24, -1*pi, pi); h_Muonphi[i]->Sumw2();
      h_Elecphi[i] = new TH1F(Form("h_Elecphi_%d",i), Form("Electron Phi %s"    ,cutflowName[i].Data()), 24, -1*pi, pi); h_Elecphi[i]->Sumw2();
   
      h_Jet1pt[i]  = new TH1F(Form("h_Jet1pt_%d", i), Form("Leading Jet pT %s"        ,cutflowName[i].Data()), 250, 0.0, 250); h_Jet1pt[i]->Sumw2();
      h_Jet2pt[i]  = new TH1F(Form("h_Jet2pt_%d", i), Form("Second Leading Jet pT %s" ,cutflowName[i].Data()), 250, 0.0, 250); h_Jet2pt[i]->Sumw2();
      h_Jet1eta[i] = new TH1F(Form("h_Jet1eta_%d",i), Form("Leading Jet Eta %s"       ,cutflowName[i].Data()), 50, -2.5, 2.5); h_Jet1eta[i]->Sumw2();
      h_Jet2eta[i] = new TH1F(Form("h_Jet2eta_%d",i), Form("Second Leading Jet Eta %s",cutflowName[i].Data()), 50, -2.5, 2.5); h_Jet2eta[i]->Sumw2();
      h_Jet1phi[i] = new TH1F(Form("h_Jet1phi_%d",i), Form("Leading Jet Phi %s"       ,cutflowName[i].Data()), 24, -1*pi, pi); h_Jet1phi[i]->Sumw2();
      h_Jet2phi[i] = new TH1F(Form("h_Jet2phi_%d",i), Form("Second Leading Jet Phi %s",cutflowName[i].Data()), 24, -1*pi, pi); h_Jet2phi[i]->Sumw2();
    
     
      h_METpt[i]  = new TH1F(Form("h_METpt_%d",i), Form("MET pT %s" ,cutflowName[i].Data()), 200, 0.0, 200); h_METpt[i]->Sumw2();
      h_METphi[i] = new TH1F(Form("h_METphi_%d",i), Form("MET Phi %s",cutflowName[i].Data()), 24, -1*pi, pi); h_METphi[i]->Sumw2();

      h_DiLepMass[i] = new TH1F(Form("h_DiLepMass_%d",i),Form("Di-Lepton Invariant Mass %s",cutflowName[i].Data()), 300, 0.0, 300); h_DiLepMass[i]->Sumw2();
      h_Num_PV[i]    = new TH1F(Form("h_Num_PV_%d",i),     Form("Num of Primary Vertex after %s",cutflowName[i].Data()), 100, 0.0, 100); h_Num_PV[i]->Sumw2();
      h_Num_Jets[i]  = new TH1F(Form("h_Num_Jets_%d",i), Form("Num. of Jets after %s",cutflowName[i].Data()), 20, 0.0, 20); h_Num_Jets[i]->Sumw2();
      h_Num_bJets[i] = new TH1F(Form("h_Num_bJets_%d",i),Form("Num. of b Jets after %s",cutflowName[i].Data()), 20, 0.0, 20); h_Num_bJets[i]->Sumw2();
      
      if(i>3){h_HT[i]      = new TH1F(Form("h_HT_%d",i), Form("Sum . pT of All Jet %s",cutflowName[i].Data()), 1000, 0.0, 1000); h_HT[i]->Sumw2();}
   }
   for (int i =0; i < 5; ++i)
   {
      h_Topo_Apla[i] = new TH1F(Form("h_Topo_Apla_%d",i), Form("Topologycal Vari. Aplanarity after %s",cutflowName[i+5].Data()), 300, -1, 2); h_Topo_Apla[i]->Sumw2();
      h_Topo_Sphe[i] = new TH1F(Form("h_Topo_Sphe_%d",i), Form("Topologycal Vari. Sphericity after %s",cutflowName[i+5].Data()), 300, -1, 2); h_Topo_Sphe[i]->Sumw2();
      h_Topo_Plan[i] = new TH1F(Form("h_Topo_Plan_%d",i), Form("Topologycal Vari. Planarity after %s",cutflowName[i+5].Data()), 300, -1, 2); h_Topo_Plan[i]->Sumw2();
      h_bTagEff[i]  = new TH1F(Form("h_bTagEff_%d",i), Form("BTagging Weight %s",cutflowName[i+5].Data()), 1000, -5, 5); h_bTagEff[i]->Sumw2();
   }
   h_PileUp      = new TH1F(Form("h_PileUp"), Form("hPileUp"), 1000, 0.0, 100); h_PileUp->Sumw2();// For Closure test // 
   h_PileUp_Up   = new TH1F(Form("h_PileUp_Up"), Form("hPileUp_Up"), 1000, 0.0, 100); h_PileUp_Up->Sumw2();// For Closure test // 
   h_PileUp_Down = new TH1F(Form("h_PileUp_Down"), Form("hPileUp_Down"), 1000, 0.0, 100); h_PileUp_Down->Sumw2();// For Closure test // 

   h_NumEl  = new TH1F(Form("h_NumEl" ), Form("hNumEl" ), 10, 0.0, 10); h_NumEl->Sumw2(); 
   h_NumMu  = new TH1F(Form("h_NumMu" ), Form("hNumMu" ), 10, 0.0, 10); h_NumMu->Sumw2(); 

   h_Num_PV_BeforePreSel = new TH1F(Form("h_Num_PV_BeforePreSel"), Form("Num. of Primary Vertex before Pre-Selection"), 100, 0.0, 100); h_Num_PV_BeforePreSel->Sumw2();// For Closure test // 
   h_Num_PV_AfterMetFilter = new TH1F(Form("h_Num_PV_AfterMetFilter"), Form("Num. of Primary Vertex After MetFilter"), 100, 0.0, 100); h_Num_PV_AfterMetFilter->Sumw2();// For Closure test // 
   h_Num_PV_AfterTrigger = new TH1F(Form("h_Num_PV_AfterTrigger"), Form("Num. of Primary Vertex After Trigger"), 100, 0.0, 100); h_Num_PV_AfterTrigger->Sumw2();// For Closure test // 

   h_Top1Mass     = new TH1F(Form("h_Top1Mass"   ), Form("Top1 Mass"   ), 1000, 0.0, 1000); h_Top1Mass->Sumw2(); 
   h_Top1pt       = new TH1F(Form("h_Top1pt"   ), Form("Top1 pt"   ), 1000, 0.0, 1000); h_Top1pt->Sumw2(); 
   h_Top1Rapidity = new TH1F(Form("h_Top1Rapidity"   ), Form("Top1 Rapidity"   ), 100, -5, 5); h_Top1Rapidity->Sumw2(); 
   h_Top1phi      = new TH1F(Form("h_Top1phi"   ), Form("Top1 phi"   ), 24, -1*pi, pi); h_Top1phi->Sumw2(); 
   h_Top1Energy   = new TH1F(Form("h_Top1Energy"   ), Form("Top1 Energy"   ), 1000, 0.0, 1000); h_Top1Energy->Sumw2(); 
   h_Top2Mass     = new TH1F(Form("h_Top2Mass" ), Form("Top2 Mass" ), 1000, 0.0, 1000); h_Top2Mass->Sumw2();
   h_Top2pt       = new TH1F(Form("h_Top2pt"   ), Form("Top2 pt"   ), 1000, 0.0, 1000); h_Top2pt->Sumw2(); 
   h_Top2Rapidity = new TH1F(Form("h_Top2Rapidity"   ), Form("Top2 Rapidity"   ), 100, -5, 5); h_Top2Rapidity->Sumw2(); 
   h_Top2phi      = new TH1F(Form("h_Top2phi"   ), Form("Top2 phi"   ), 24, -1*pi, pi); h_Top2phi->Sumw2(); 
   h_Top2Energy = new TH1F(Form("h_Top2Energy"   ), Form("Top2 Energy"   ), 1000, 0.0, 1000); h_Top2Energy->Sumw2(); 

   h_TopMass       = new TH1F(Form("h_TopMass"   ), Form("Top Mass"   ), 1000, 0.0, 1000); h_TopMass->Sumw2(); 
   h_Toppt         = new TH1F(Form("h_Toppt"   ), Form("Top pt"   ), 1000, 0.0, 1000); h_Toppt->Sumw2(); 
   h_TopRapidity   = new TH1F(Form("h_TopRapidity"   ), Form("Top Rapidity"   ), 100, -5, 5); h_TopRapidity->Sumw2(); 
   h_Topphi        = new TH1F(Form("h_Topphi"   ), Form("Top phi"   ), 24, -1*pi, pi); h_Topphi->Sumw2(); 
   h_TopEnergy     = new TH1F(Form("h_TopEnergy"   ), Form("Top Energy"   ), 1000, 0.0, 1000); h_TopEnergy->Sumw2(); 
   h_AnTopMass     = new TH1F(Form("h_AnTopMass" ), Form("AnTop Mass" ), 1000, 0.0, 1000); h_AnTopMass->Sumw2();
   h_AnToppt       = new TH1F(Form("h_AnToppt"   ), Form("AnTop pt"   ), 1000, 0.0, 1000); h_AnToppt->Sumw2(); 
   h_AnTopRapidity = new TH1F(Form("h_AnTopRapidity"   ), Form("AnTop Rapidity"   ), 100, -5, 5); h_AnTopRapidity->Sumw2(); 
   h_AnTopphi      = new TH1F(Form("h_AnTopphi"   ), Form("AnTop phi"   ), 24, -1*pi, pi); h_AnTopphi->Sumw2(); 
   h_AnTopEnergy   = new TH1F(Form("h_AnTopEnergy"   ), Form("AnTop Energy"   ), 1000, 0.0, 1000); h_AnTopEnergy->Sumw2(); 


   h_W1Mass     = new TH1F(Form("h_W1Mass"  ), Form("W1 Mass" ), 300, 0.0, 300); h_W1Mass->Sumw2();
   h_W2Mass     = new TH1F(Form("h_W2Mass"  ), Form("W2 Mass" ), 300, 0.0, 300); h_W2Mass->Sumw2();

   h_W1Mt     = new TH1F(Form("h_W1Mt"  ), Form("W1 Transverse Mass" ), 300, 0.0, 300); h_W1Mt->Sumw2();
   h_W2Mt     = new TH1F(Form("h_W2Mt"  ), Form("W2 Transverse Mass" ), 300, 0.0, 300); h_W2Mt->Sumw2();

   h_bJet1Energy = new TH1F(Form("h_bJet1Energy" ), Form("Leading bJet Energy"   ), 500, 0.0, 500); h_bJet1Energy->Sumw2(); 
   h_bJet2Energy = new TH1F(Form("h_bJet2Energy" ), Form("Second Leading bJet Energy" ), 500, 0.0, 500); h_bJet2Energy->Sumw2();

   h_bJetEnergy   = new TH1F(Form("h_bJetEnergy" ), Form("bJet Energy"   ), 1000, 0.0, 1000); h_bJetEnergy->Sumw2(); 
   h_AnbJetEnergy = new TH1F(Form("h_AnbJetEnergy" ), Form("b-barJet Energy" ), 1000, 0.0, 1000); h_AnbJetEnergy->Sumw2(); 
   h_bJetPt   = new TH1F(Form("h_bJetPt" ), Form("bJet Pt"   ), 1000, 0.0, 1000); h_bJetPt->Sumw2(); 
   h_AnbJetPt = new TH1F(Form("h_AnbJetPt" ), Form("b-barJet Pt" ), 1000, 0.0, 1000); h_AnbJetPt->Sumw2();

   h_Lep1Energy = new TH1F(Form("h_Lep1Energy" ), Form("Leading Lepton Energy"   ), 400, 0.0, 400); h_Lep1Energy->Sumw2(); 
   h_Lep2Energy = new TH1F(Form("h_Lep2Energy" ), Form("Second Leading Lepton Energy" ), 400, 0.0, 400); h_Lep2Energy->Sumw2();

   h_LepEnergy   = new TH1F(Form("h_LepEnergy" ), Form("Lepton Energy"   ), 400, 0.0, 400); h_LepEnergy->Sumw2(); 
   h_AnLepEnergy = new TH1F(Form("h_AnLepEnergy" ), Form("Anti-Lepton Energy" ), 400, 0.0, 400); h_AnLepEnergy->Sumw2();

   h_Nu1Energy = new TH1F(Form("h_Nu1Energy" ), Form("Leading Nuetrino Energy"   ), 400, 0.0, 400); h_Nu1Energy->Sumw2(); 
   h_Nu2Energy = new TH1F(Form("h_Nu2Energy" ), Form("Second Leading Nuetrino Energy" ), 400, 0.0, 400); h_Nu2Energy->Sumw2(); 

   h_NuEnergy   = new TH1F(Form("h_NuEnergy" ), Form("Nuetrino Energy"   ), 400, 0.0, 400); h_NuEnergy->Sumw2(); 
   h_AnNuEnergy = new TH1F(Form("h_AnNuEnergy" ), Form("anti-Nuetrino Energy" ), 400, 0.0, 400); h_AnNuEnergy->Sumw2();

   h_GenNuEnergy   = new TH1F(Form("h_GenNuEnergy" ), Form("Nuetrino Energy At Genrator Level"   ), 2000, 0.0, 2000); h_GenNuEnergy->Sumw2();
   h_GenAnNuEnergy = new TH1F(Form("h_GenAnNuEnergy" ), Form("anti-Nuetrino Energy At Genrator Level" ), 2000, 0.0, 2000); h_GenAnNuEnergy->Sumw2();

   h_CPO3_reco         = new TH1F(Form("h_CPO3_reco"   ), Form("CPO3_reco"   ), 200, -10, 10); h_CPO3_reco->Sumw2();
   h_CPO3_reco_JPRUp   = new TH1F(Form("h_CPO3_reco_JPRUp"     ), Form("CPO3_reco_JPRUp"     ), 200, -10, 10); h_CPO3_reco_JPRUp->Sumw2();
   h_CPO3_reco_JPRDown = new TH1F(Form("h_CPO3_reco_JPRDown"   ), Form("CPO3_reco_JPRDown"   ), 200, -10, 10); h_CPO3_reco_JPRDown->Sumw2();
   h_CPO3_reco_TopRes  = new TH1F(Form("h_CPO3_reco_TopRes"   ), Form("CPO3_reco in the Top Mass Window"   ), 200, -10, 10); h_CPO3_reco_TopRes->Sumw2();
   h_CPOb_reco         = new TH1F(Form("h_CPOb_reco"   ), Form("CPOb_reco"   ), 200, -10, 10); h_CPOb_reco->Sumw2();
   h_CPO5_reco         = new TH1F(Form("h_CPO5_reco"   ), Form("CPO5_reco"   ), 200, -10, 10); h_CPO5_reco->Sumw2();

   h_CPO3_evtweight  = new TH1F(Form("h_CPO3_evtweight"   ), Form("CPO3_evtweight"   ), 10, 0, 10); h_CPO3_evtweight->Sumw2();
   h_CPOb_evtweight  = new TH1F(Form("h_CPOb_evtweight"   ), Form("CPOb_evtweight"   ), 10, 0, 10); h_CPOb_evtweight->Sumw2();
   h_CPO5_evtweight  = new TH1F(Form("h_CPO5_evtweight"   ), Form("CPO5_evtweight"   ), 10, 0, 10); h_CPO5_evtweight->Sumw2();

   h_BjorkenX1 = new TH1F(Form("h_BjorkenX1" ), Form("Bjorken X1" ), 120, -0.1  , 1.1); h_BjorkenX1->Sumw2();  
   h_BjorkenX2 = new TH1F(Form("h_BjorkenX2" ), Form("Bjorken X2" ), 120, -0.1, 1.1); h_BjorkenX2->Sumw2(); 
   h_BjorkenX3 = new TH1F(Form("h_BjorkenX3" ), Form("Bjorken X3" ), 2000, 0.0  , 2000); h_BjorkenX3->Sumw2();

   /////////////////////////////
   /// To Check Evetn Weight ///
   /////////////////////////////

   h_LepbJetMass          = new TH1F(Form("h_LepbJetMass"), Form("Lep + bJet Invariant Mass"), 1000, 0, 1000); h_LepbJetMass->Sumw2();
   h_AnLepbJetMass        = new TH1F(Form("h_AnLepbJetMass"), Form("AnLep + bJet Invariant Mass"), 1000, 0, 1000); h_AnLepbJetMass->Sumw2();
   h2_TopMassVsLepBMass    = new TH2F(Form("TopMassVsLepBMass"), Form("Top Mass Vs (Lep+B) Mass"), 1000, 0, 1000, 1000,0,1000); h2_TopMassVsLepBMass->Sumw2();
   h2_AnTopMassVsLepBMass  = new TH2F(Form("AnTopMassVsLepBMass"), Form("Anti-Top Mass Vs (Lep+B) Mass"), 1000, 0, 1000, 1000,0,1000); h2_AnTopMassVsLepBMass->Sumw2();
   h2_AnLepBMassVsLepBMass     = new TH2F(Form("AnLepBMassVsLepBMass"), Form("(AntiLep+B) Mass Vs (Lep+B) Mass"), 1000, 0, 1000, 1000,0,1000); h2_AnLepBMassVsLepBMass->Sumw2();

   h_Diff_BBar_Pt     = new TH1F(Form("h_Diff_BBar_Pt"    ), Form("PT Difference between B and B bar "     ), 500, 0, 500); h_Diff_BBar_Pt->Sumw2();
   h_Diff_BBar_Energy = new TH1F(Form("h_Diff_BBar_Energy"), Form("Energy Difference between B and B bar"  ), 500, 0, 500); h_Diff_BBar_Energy->Sumw2();
   h_Diff_BBar_P      = new TH1F(Form("h_Diff_BBar_P"     ), Form("Momentum Difference between B and B bar"), 500, 0, 500); h_Diff_BBar_P->Sumw2();

   h_Diff_LepAnLep_Pt     = new TH1F(Form("h_Diff_LepAnLep_Pt"    ), Form("PT Difference between Lep. And AntiLep."      ), 500, 0, 500); h_Diff_LepAnLep_Pt->Sumw2();
   h_Diff_LepAnLep_Energy = new TH1F(Form("h_Diff_LepAnLep_Energy"), Form("Energy Difference between Lep. And AntiLep."  ), 500, 0, 500); h_Diff_LepAnLep_Energy->Sumw2();
   h_Diff_LepAnLep_P      = new TH1F(Form("h_Diff_LepAnLep_P"     ), Form("Momentum Difference between Lep. And AntiLep."), 500, 0, 500); h_Diff_LepAnLep_P->Sumw2();

   h_LepAnLepEngCheckO3 = new TH1F(Form("h_LepAnLepEngCheckO3"), Form("Lep1 and Lep2 pt check wheather they over the 300GeV for O3"), 100, 0, 100); h_LepAnLepEngCheckO3->Sumw2();
   h_LepAnLepEngCheckOb = new TH1F(Form("h_LepAnLepEngCheckOb"), Form("Lep1 and Lep2 pt check wheather they over the 300GeV for Ob"), 100, 0, 100); h_LepAnLepEngCheckOb->Sumw2();
   h_LepAnLepEngCheckO5 = new TH1F(Form("h_LepAnLepEngCheckO5"), Form("Lep1 and Lep2 pt check wheather they over the 300GeV for O5"), 100, 0, 100); h_LepAnLepEngCheckO5->Sumw2();

   h_PileUpCheck  = new TH1F(Form("h_PileUpCheck"),  Form("NewPileUpMethod - OriginMethod"), 4000, -2, 2); h_PileUpCheck->Sumw2();
   //////////////////////////////////
   /// For CP-Violation Variables ///
   //////////////////////////////////
   for (int i =0; i < 13; ++i)
   {
      h_Reco_CPO_[i] = new TH1F(Form("h_Reco_CPO%d",i+1 ), Form("CPO%d",i+1   ), 200, -10, 10); h_Reco_CPO_[i]->Sumw2();
      h_Reco_CPO_ReRange_[i] = new TH1F(Form("h_Reco_CPO%d_ReRange",i+1 ), Form("CPO%d",i+1   ), 40, -2, 2); h_Reco_CPO_ReRange_[i]->Sumw2();
   }
   ////////////////////////////////
   /// To Check O3 dilution ... ///
   ////////////////////////////////
   h_CPO3_Plus              = new TH1F(Form("h_CPO3_Plus"),              Form("CPO3 Positive"), 200, -10, 10); h_CPO3_Plus->Sumw2(); 
   h_CPO3_Minus             = new TH1F(Form("h_CPO3_Minus"),             Form("CPO3 Negative"), 200, -10, 10); h_CPO3_Minus->Sumw2(); 
   h_CPO3_JPRUp_Plus_Plus   = new TH1F(Form("h_CPO3_JPRUp_Plus_Plus"),   Form("CPO3_JPRUp Plus -> Plus"),   200, -10, 10); h_CPO3_JPRUp_Plus_Plus->Sumw2();                   
   h_CPO3_JPRUp_Plus_Minus  = new TH1F(Form("h_CPO3_JPRUp_Plus_Minus"),  Form("CPO3_JPRUp Plus -> Minus"),  200, -10, 10); h_CPO3_JPRUp_Plus_Minus->Sumw2();                   
   h_CPO3_JPRUp_Minus_Minus = new TH1F(Form("h_CPO3_JPRUp_Minus_Minus"), Form("CPO3_JPRUp Minus -> Minus"), 200, -10, 10); h_CPO3_JPRUp_Minus_Minus->Sumw2();                   
   h_CPO3_JPRUp_Minus_Plus  = new TH1F(Form("h_CPO3_JPRUp_Minus_Plus"),  Form("CPO3_JPRUp Minus -> Minus"), 200, -10, 10); h_CPO3_JPRUp_Minus_Plus->Sumw2();                   

   h_CPO3_JPRDown_Plus_Plus   = new TH1F(Form("h_CPO3_JPRDown_Plus_Plus"),   Form("CPO3_JPRDown Plus-> Plus"),   200, -10, 10); h_CPO3_JPRDown_Plus_Plus->Sumw2();
   h_CPO3_JPRDown_Plus_Minus  = new TH1F(Form("h_CPO3_JPRDown_Plus_Minus"),  Form("CPO3_JPRDown Plus-> Minus"),  200, -10, 10); h_CPO3_JPRDown_Plus_Minus->Sumw2();
   h_CPO3_JPRDown_Minus_Minus = new TH1F(Form("h_CPO3_JPRDown_Minus_Minus"), Form("CPO3_JPRDown Minus-> Minus"), 200, -10, 10); h_CPO3_JPRDown_Minus_Minus->Sumw2();
   h_CPO3_JPRDown_Minus_Plus  = new TH1F(Form("h_CPO3_JPRDown_Minus_Plus"),  Form("CPO3_JPRDown Minus-> Minus"), 200, -10, 10); h_CPO3_JPRDown_Minus_Plus->Sumw2();


   h_CPO3reco_evtweight             = new TH1D(Form("h_CPO3reco_evtweight"),             Form("Event Weight for CPO3 at Generator Level"),   10, 0, 10); h_CPO3reco_evtweight->Sumw2();
   h_CPO3reco_evtweight_EtaVariUp   = new TH1D(Form("h_CPO3reco_evtweight_EtaVariUp"),   Form("Event Weight for CPO3 at Generator Level"),   10, 0, 10); h_CPO3reco_evtweight_EtaVariUp->Sumw2();
   h_CPO3reco_evtweight_EtaVariDown = new TH1D(Form("h_CPO3reco_evtweight_EtaVariDown"), Form("Event Weight for CPO3 at Generator Level"),   10, 0, 10); h_CPO3reco_evtweight_EtaVariDown->Sumw2();
   h_CPO3reco_evtweight_PhiVariUp   = new TH1D(Form("h_CPO3reco_evtweight_PhiVariUp"),   Form("Event Weight for CPO3 at Generator Level"),   10, 0, 10); h_CPO3reco_evtweight_PhiVariUp->Sumw2();
   h_CPO3reco_evtweight_PhiVariDown = new TH1D(Form("h_CPO3reco_evtweight_PhiVariDown"), Form("Event Weight for CPO3 at Generator Level"),   10, 0, 10); h_CPO3reco_evtweight_PhiVariDown->Sumw2();

   h_CPO3_EtaVariUp   = new TH1F(Form("h_CPO3_EtaVariUp"),   Form("CPO3_EtaVariUp"),   200, -10, 10); h_CPO3_EtaVariUp->Sumw2(); 
   h_CPO3_EtaVariUp_Plus_Plus   = new TH1F(Form("h_CPO3_EtaVariUp_Plus_Plus"),   Form("CPO3_EtaVariUp Plus-> Plus"),   200, -10, 10); h_CPO3_EtaVariUp_Plus_Plus->Sumw2();
   h_CPO3_EtaVariUp_Plus_Minus  = new TH1F(Form("h_CPO3_EtaVariUp_Plus_Minus"),  Form("CPO3_EtaVariUp Plus-> Minus"),  200, -10, 10); h_CPO3_EtaVariUp_Plus_Minus->Sumw2();
   h_CPO3_EtaVariUp_Minus_Minus = new TH1F(Form("h_CPO3_EtaVariUp_Minus_Minus"), Form("CPO3_EtaVariUp Minus-> Minus"), 200, -10, 10); h_CPO3_EtaVariUp_Minus_Minus->Sumw2();
   h_CPO3_EtaVariUp_Minus_Plus  = new TH1F(Form("h_CPO3_EtaVariUp_Minus_Plus"),  Form("CPO3_EtaVariUp Minus-> Minus"), 200, -10, 10); h_CPO3_EtaVariUp_Minus_Plus->Sumw2();

   h_CPO3_EtaVariDown   = new TH1F(Form("h_CPO3_EtaVariDown"),   Form("CPO3_EtaVariDown"),   200, -10, 10); h_CPO3_EtaVariDown->Sumw2(); 
   h_CPO3_EtaVariDown_Plus_Plus   = new TH1F(Form("h_CPO3_EtaVariDown_Plus_Plus"),   Form("CPO3_EtaVariDown Plus-> Plus"),   200, -10, 10); h_CPO3_EtaVariDown_Plus_Plus->Sumw2();
   h_CPO3_EtaVariDown_Plus_Minus  = new TH1F(Form("h_CPO3_EtaVariDown_Plus_Minus"),  Form("CPO3_EtaVariDown Plus-> Minus"),  200, -10, 10); h_CPO3_EtaVariDown_Plus_Minus->Sumw2();
   h_CPO3_EtaVariDown_Minus_Minus = new TH1F(Form("h_CPO3_EtaVariDown_Minus_Minus"), Form("CPO3_EtaVariDown Minus-> Minus"), 200, -10, 10); h_CPO3_EtaVariDown_Minus_Minus->Sumw2();
   h_CPO3_EtaVariDown_Minus_Plus  = new TH1F(Form("h_CPO3_EtaVariDown_Minus_Plus"),  Form("CPO3_EtaVariDown Minus-> Minus"), 200, -10, 10); h_CPO3_EtaVariDown_Minus_Plus->Sumw2();

   h_CPO3_PhiVariUp   = new TH1F(Form("h_CPO3_PhiVariUp"),   Form("CPO3_PhiVariUp"),   200, -10, 10); h_CPO3_PhiVariUp->Sumw2(); 
   h_CPO3_PhiVariUp_Plus_Plus   = new TH1F(Form("h_CPO3_PhiVariUp_Plus_Plus"),   Form("CPO3_PhiVariUp Plus-> Plus"),   200, -10, 10); h_CPO3_PhiVariUp_Plus_Plus->Sumw2();
   h_CPO3_PhiVariUp_Plus_Minus  = new TH1F(Form("h_CPO3_PhiVariUp_Plus_Minus"),  Form("CPO3_PhiVariUp Plus-> Minus"),  200, -10, 10); h_CPO3_PhiVariUp_Plus_Minus->Sumw2();
   h_CPO3_PhiVariUp_Minus_Minus = new TH1F(Form("h_CPO3_PhiVariUp_Minus_Minus"), Form("CPO3_PhiVariUp Minus-> Minus"), 200, -10, 10); h_CPO3_PhiVariUp_Minus_Minus->Sumw2();
   h_CPO3_PhiVariUp_Minus_Plus  = new TH1F(Form("h_CPO3_PhiVariUp_Minus_Plus"),  Form("CPO3_PhiVariUp Minus-> Minus"), 200, -10, 10); h_CPO3_PhiVariUp_Minus_Plus->Sumw2();

   h_CPO3_PhiVariDown   = new TH1F(Form("h_CPO3_PhiVariDown"),   Form("CPO3_PhiVariDown"),   200, -10, 10); h_CPO3_PhiVariDown->Sumw2();
   h_CPO3_PhiVariDown_Plus_Plus   = new TH1F(Form("h_CPO3_PhiVariDown_Plus_Plus"),   Form("CPO3_PhiVariDown Plus-> Plus"),   200, -10, 10); h_CPO3_PhiVariDown_Plus_Plus->Sumw2();
   h_CPO3_PhiVariDown_Plus_Minus  = new TH1F(Form("h_CPO3_PhiVariDown_Plus_Minus"),  Form("CPO3_PhiVariDown Plus-> Minus"),  200, -10, 10); h_CPO3_PhiVariDown_Plus_Minus->Sumw2();
   h_CPO3_PhiVariDown_Minus_Minus = new TH1F(Form("h_CPO3_PhiVariDown_Minus_Minus"), Form("CPO3_PhiVariDown Minus-> Minus"), 200, -10, 10); h_CPO3_PhiVariDown_Minus_Minus->Sumw2();
   h_CPO3_PhiVariDown_Minus_Plus  = new TH1F(Form("h_CPO3_PhiVariDown_Minus_Plus"),  Form("CPO3_PhiVariDown Minus-> Minus"), 200, -10, 10); h_CPO3_PhiVariDown_Minus_Plus->Sumw2();

   h_bTagWeight   = new TH1F(Form("h_bTagWeight"),   Form("bTag Weight "),   200, -2, 2); h_bTagWeight->Sumw2();
}

void ssb_analysis::DeclareHistosSyst(int index_)
{
   for (int j =0; j <10; j++) // Cutflow Loop //
   {
      // Cut Flow histogram //
      h_cf_sys_NLeptons[index_][j]       = new TH1F(Form("_h_cf_sys_NLeptons_%s_%d_"      ,v_SystFullName[index_].Data(), j),Form("Num. Lepton %s (%s)"                ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 20  , 0    , 20  ); h_cf_sys_NLeptons[index_][j]->Sumw2();         
      h_cf_sys_Lep1pt[index_][j]         = new TH1F(Form("_h_cf_sys_Lep1pt_%s_%d_"        ,v_SystFullName[index_].Data(), j),Form("Leading Lepton pT %s (%s)"          ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 1000, 0    , 1000); h_cf_sys_Lep1pt[index_][j]->Sumw2();           
      h_cf_sys_Lep1phi[index_][j]        = new TH1F(Form("_h_cf_sys_Lep1phi_%s_%d_"       ,v_SystFullName[index_].Data(), j),Form("Leading Lepton phi %s (%s)"        ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 24  , -1*pi, 1*pi); h_cf_sys_Lep1phi[index_][j]->Sumw2();          
      h_cf_sys_Lep1eta[index_][j]        = new TH1F(Form("_h_cf_sys_Lep1eta_%s_%d_"       ,v_SystFullName[index_].Data(), j),Form("Leading Lepton eta %s (%s)"        ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 50  , -2.5 , 2.5 ); h_cf_sys_Lep1eta[index_][j]->Sumw2();          
      h_cf_sys_Lep2pt[index_][j]         = new TH1F(Form("_h_cf_sys_Lep2pt_%s_%d_"        ,v_SystFullName[index_].Data(), j),Form("Second Leading Lepton pT %s (%s)"   ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 1000, 0    , 1000); h_cf_sys_Lep2pt[index_][j]->Sumw2();           
      h_cf_sys_Lep2phi[index_][j]        = new TH1F(Form("_h_cf_sys_Lep2phi_%s_%d_"       ,v_SystFullName[index_].Data(), j),Form("Second Leading Lepton phi %s (%s)"  ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 24  , -1*pi, 1*pi); h_cf_sys_Lep2phi[index_][j]->Sumw2();      
      h_cf_sys_Lep2eta[index_][j]        = new TH1F(Form("_h_cf_sys_Lep2eta_%s_%d_"       ,v_SystFullName[index_].Data(), j),Form("Second Leading Lepton eta %s (%s)"  ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 50  , -2.5 , 2.5 ); h_cf_sys_Lep2eta[index_][j]->Sumw2();       
      h_cf_sys_dilep_inv_mass[index_][j] = new TH1F(Form("_h_cf_sys_dilep_inv_mass_%s_%d_",v_SystFullName[index_].Data(), j),Form("Invariant Mass of Dilepton %s (%s)" ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 1000, 0    , 1000); h_cf_sys_dilep_inv_mass[index_][j]->Sumw2();
      h_cf_sys_Jet1pt[index_][j]         = new TH1F(Form("_h_cf_sys_Jet1pt_%s_%d_"        ,v_SystFullName[index_].Data(), j),Form("Leading Jet pT %s (%s)"             ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 1000, 0    , 1000); h_cf_sys_Jet1pt[index_][j]->Sumw2();       
      h_cf_sys_Jet1phi[index_][j]        = new TH1F(Form("_h_cf_sys_Jet1phi_%s_%d_"       ,v_SystFullName[index_].Data(), j),Form("Leading Jet phi %s (%s)"            ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 24  , -1*pi, 1*pi); h_cf_sys_Jet1phi[index_][j]->Sumw2();       
      h_cf_sys_Jet1eta[index_][j]        = new TH1F(Form("_h_cf_sys_Jet1eta_%s_%d_"       ,v_SystFullName[index_].Data(), j),Form("Leading Jet eta %s (%s)"            ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 50  , -2.5 , 2.5 ); h_cf_sys_Jet1eta[index_][j]->Sumw2();       
      h_cf_sys_Jet2pt[index_][j]         = new TH1F(Form("_h_cf_sys_Jet2pt_%s_%d_"        ,v_SystFullName[index_].Data(), j),Form("Second Leading Jet pT %s (%s)"      ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 1000, 0    , 1000); h_cf_sys_Jet2pt[index_][j]->Sumw2();        
      h_cf_sys_Jet2phi[index_][j]        = new TH1F(Form("_h_cf_sys_Jet2phi_%s_%d_"       ,v_SystFullName[index_].Data(), j),Form("Second Leading Jet phi %s (%s)"     ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 24  , -1*pi, 1*pi); h_cf_sys_Jet2phi[index_][j]->Sumw2();       
      h_cf_sys_Jet2eta[index_][j]        = new TH1F(Form("_h_cf_sys_Jet2eta_%s_%d_"       ,v_SystFullName[index_].Data(), j),Form("Second Leading Jet eta %s (%s)"     ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 50  , -2.5 , 2.5 ); h_cf_sys_Jet2eta[index_][j]->Sumw2();      
      h_cf_sys_NJets[index_][j]          = new TH1F(Form("_h_cf_sys_NJets_%s_%d_"         ,v_SystFullName[index_].Data(), j),Form("Num. of Jets %s (%s)"               ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 20  , 0    , 20  ); h_cf_sys_NJets[index_][j]->Sumw2();        
      h_cf_sys_metpt[index_][j]          = new TH1F(Form("_h_cf_sys_metpt_%s_%d_"         ,v_SystFullName[index_].Data(), j),Form("MET %s (%s)"                        ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 1000, 0    , 1000); h_cf_sys_metpt[index_][j]->Sumw2();        
      h_cf_sys_metphi[index_][j]         = new TH1F(Form("_h_cf_sys_metphi_%s_%d_"        ,v_SystFullName[index_].Data(), j),Form("MET Phi %s (%s)"                    ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 24  , -1*pi, 1*pi); h_cf_sys_metphi[index_][j]->Sumw2();       
      h_cf_sys_Nbjets[index_][j]         = new TH1F(Form("_h_cf_sys_Nbjets_%s_%d_"        ,v_SystFullName[index_].Data(), j),Form("Num. of b-jets %s (%s)"             ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 20  , 0    , 20  ); h_cf_sys_Nbjets[index_][j]->Sumw2();        
      h_cf_sys_NPV[index_][j]            = new TH1F(Form("_h_cf_sys_NPV_%s_%d_"           ,v_SystFullName[index_].Data(), j),Form("Num. of Primary Vertex %s (%s)"     ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 100 , 0    , 100 ); h_cf_sys_NPV[index_][j]->Sumw2();

      h_sys_Lep1pt[index_][j]  = new TH1F(Form("h_sys_Lep1pt_%s_%d"   ,v_SystFullName[index_].Data(), j), Form("Leading Lepton pT %s (%s)"        ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 250, 0.0, 250); h_sys_Lep1pt[index_][j]->Sumw2(); 
      h_sys_Lep2pt[index_][j]  = new TH1F(Form("h_sys_Lep2pt_%s_%d_"  ,v_SystFullName[index_].Data(), j), Form("Second Leading Lepton pT %s (%s)" ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 250, 0.0, 250); h_sys_Lep2pt[index_][j]->Sumw2();
      h_sys_Lep1eta[index_][j] = new TH1F(Form("h_sys_Lep1eta_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Leading Lepton Eta    %s (%s)"    ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 50, -2.5, 2.5); h_sys_Lep1eta[index_][j]->Sumw2();
      h_sys_Lep2eta[index_][j] = new TH1F(Form("h_sys_Lep2eta_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Second Leading Lepton Eta %s (%s)",cutflowName[j].Data(),v_SystFullName[index_].Data()), 50, -2.5, 2.5); h_sys_Lep2eta[index_][j]->Sumw2();
      h_sys_Lep1phi[index_][j] = new TH1F(Form("h_sys_Lep1phi_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Leading Lepton Phi %s (%s)"       ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 24, -1*pi, pi); h_sys_Lep1phi[index_][j]->Sumw2();
      h_sys_Lep2phi[index_][j] = new TH1F(Form("h_sys_Lep2phi_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Second Leading Lepton Phi %s (%s)",cutflowName[j].Data(),v_SystFullName[index_].Data()), 24, -1*pi, pi); h_sys_Lep2phi[index_][j]->Sumw2();
   
      h_sys_Muonpt[index_][j]  = new TH1F(Form("h_sys_Muonpt_%s_%d_"  ,v_SystFullName[index_].Data(), j), Form("Muon pT %s (%s)"         ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 250, 0.0, 250); h_sys_Muonpt[index_][j]->Sumw2();
      h_sys_Elecpt[index_][j]  = new TH1F(Form("h_sys_Elecpt_%s_%d_"  ,v_SystFullName[index_].Data(), j), Form("Electron pT %s (%s)"     ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 250, 0.0, 250); h_sys_Elecpt[index_][j]->Sumw2();
      h_sys_Muoneta[index_][j] = new TH1F(Form("h_sys_Muoneta_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Muon Eta %s (%s)"        ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 50, -2.5, 2.5); h_sys_Muoneta[index_][j]->Sumw2();
      h_sys_Eleceta[index_][j] = new TH1F(Form("h_sys_Eleceta_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Electron Eta %s (%s)"    ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 50, -2.5, 2.5); h_sys_Eleceta[index_][j]->Sumw2();
      h_sys_Muonphi[index_][j] = new TH1F(Form("h_sys_Muonphi_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Muon Phi %s (%s)"        ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 24, -1*pi, pi); h_sys_Muonphi[index_][j]->Sumw2();
      h_sys_Elecphi[index_][j] = new TH1F(Form("h_sys_Elecphi_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Electron Phi %s (%s)"    ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 24, -1*pi, pi); h_sys_Elecphi[index_][j]->Sumw2();
   
      h_sys_Jet1pt[index_][j]  = new TH1F(Form("h_sys_Jet1pt_%s_%d_",  v_SystFullName[index_].Data(), j), Form("Leading Jet pT %s (%s)"        ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 250, 0.0, 250); h_sys_Jet1pt[index_][j]->Sumw2();
      h_sys_Jet2pt[index_][j]  = new TH1F(Form("h_sys_Jet2pt_%s_%d_",  v_SystFullName[index_].Data(), j), Form("Second Leading Jet pT %s (%s)" ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 250, 0.0, 250); h_sys_Jet2pt[index_][j]->Sumw2();
      h_sys_Jet1eta[index_][j] = new TH1F(Form("h_sys_Jet1eta_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Leading Jet Eta %s (%s)"       ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 50, -2.5, 2.5); h_sys_Jet1eta[index_][j]->Sumw2();
      h_sys_Jet2eta[index_][j] = new TH1F(Form("h_sys_Jet2eta_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Second Leading Jet Eta %s (%s)",cutflowName[j].Data(),v_SystFullName[index_].Data()), 50, -2.5, 2.5); h_sys_Jet2eta[index_][j]->Sumw2();
      h_sys_Jet1phi[index_][j] = new TH1F(Form("h_sys_Jet1phi_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Leading Jet Phi %s (%s)"       ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 24, -1*pi, pi); h_sys_Jet1phi[index_][j]->Sumw2();
      h_sys_Jet2phi[index_][j] = new TH1F(Form("h_sys_Jet2phi_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("Second Leading Jet Phi %s (%s)",cutflowName[j].Data(),v_SystFullName[index_].Data()), 24, -1*pi, pi); h_sys_Jet2phi[index_][j]->Sumw2();
    
     
      h_sys_METpt[index_][j]  = new TH1F(Form("h_sys_METpt_%s_%d_"  ,v_SystFullName[index_].Data(), j), Form("MET pT %s (%s)" ,cutflowName[j].Data(),v_SystFullName[index_].Data()), 200, 0.0, 200); h_sys_METpt[index_][j]->Sumw2();
      h_sys_METphi[index_][j] = new TH1F(Form("h_sys_METphi_%s_%d_" ,v_SystFullName[index_].Data(), j), Form("MET Phi %s (%s)",cutflowName[j].Data(),v_SystFullName[index_].Data()), 24, -1*pi, pi); h_sys_METphi[index_][j]->Sumw2();
   
      h_sys_DiLepMass[index_][j] = new TH1F(Form("h_sys_DiLepMass_%s_%d_"     ,v_SystFullName[index_].Data(), j),Form("Di-Lepton Invariant Mass %s (%s)",cutflowName[j].Data(),v_SystFullName[index_].Data()), 300, 0.0, 300); h_sys_DiLepMass[index_][j]->Sumw2();
      h_sys_Num_PV[index_][j]    = new TH1F(Form("h_sys_Num_PV_%s_%d_"        ,v_SystFullName[index_].Data(), j),     Form("Num of Primary Vertex after %s (%s)",cutflowName[j].Data(),v_SystFullName[index_].Data()), 100, 0.0, 100); h_sys_Num_PV[index_][j]->Sumw2();
      h_sys_Num_Jets[index_][j]  = new TH1F(Form("h_sys_Num_Jets_%s_%d_"      ,v_SystFullName[index_].Data(), j), Form("Num. of Jets after %s (%s)",cutflowName[j].Data(),v_SystFullName[index_].Data()), 20, 0.0, 20); h_sys_Num_Jets[index_][j]->Sumw2();
      h_sys_Num_bJets[index_][j] = new TH1F(Form("h_sys_Num_bJets_%s_%d_"     ,v_SystFullName[index_].Data(), j),Form("Num. of b Jets after %s (%s)",cutflowName[j].Data(),v_SystFullName[index_].Data()), 20, 0.0, 20); h_sys_Num_bJets[index_][j]->Sumw2();
   } // Cut Flow LOOP
   for (int j =0; j < 13; ++j)
   {
      h_sys_Reco_CPO_[index_][j]         = new TH1F(Form("h_sys_Reco_CPO%d_%s"        ,j+1 ,v_SystFullName[index_].Data() ), Form("CPO%d (%s)",j+1,v_SystFullName[index_].Data()   ), 200, -10, 10); h_sys_Reco_CPO_[index_][j]->Sumw2();
      h_sys_Reco_CPO_ReRange_[index_][j] = new TH1F(Form("h_sys_Reco_CPO%d_%s_ReRange",j+1 ,v_SystFullName[index_].Data() ), Form("CPO%d (%s)",j+1,v_SystFullName[index_].Data()   ), 40, -2, 2); h_sys_Reco_CPO_ReRange_[index_][j]->Sumw2();
   }

   h_sys_Top1Mass_[index_]     = new TH1F(Form("h_sys_Top1Mass_%s"     ,v_SystFullName[index_].Data() ), Form("Top1 Mass (%s)"     ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_Top1Mass_[index_]->Sumw2(); 
   h_sys_Top1pt_[index_]       = new TH1F(Form("h_sys_Top1pt_%s"       ,v_SystFullName[index_].Data() ), Form("Top1 pt (%s)"       ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_Top1pt_[index_]->Sumw2(); 
   h_sys_Top1Rapidity_[index_] = new TH1F(Form("h_sys_Top1Rapidity_%s" ,v_SystFullName[index_].Data() ), Form("Top1 Rapidity (%s)" ,v_SystFullName[index_].Data() ), 100, -5, 5); h_sys_Top1Rapidity_[index_]->Sumw2(); 
   h_sys_Top1phi_[index_]      = new TH1F(Form("h_sys_Top1phi_%s"      ,v_SystFullName[index_].Data() ), Form("Top1 phi (%s)"      ,v_SystFullName[index_].Data() ), 24, -1*pi, pi); h_sys_Top1phi_[index_]->Sumw2(); 
   h_sys_Top1Energy_[index_]   = new TH1F(Form("h_sys_Top1Energy_%s"   ,v_SystFullName[index_].Data() ), Form("Top1 Energy (%s)"   ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_Top1Energy_[index_]->Sumw2(); 
   h_sys_Top2Mass_[index_]     = new TH1F(Form("h_sys_Top2Mass_%s"     ,v_SystFullName[index_].Data() ), Form("Top2 Mass (%s)"     ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_Top2Mass_[index_]->Sumw2();
   h_sys_Top2pt_[index_]       = new TH1F(Form("h_sys_Top2pt_%s"       ,v_SystFullName[index_].Data() ), Form("Top2 pt (%s)"       ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_Top2pt_[index_]->Sumw2(); 
   h_sys_Top2Rapidity_[index_] = new TH1F(Form("h_sys_Top2Rapidity_%s" ,v_SystFullName[index_].Data() ), Form("Top2 Rapidity (%s)" ,v_SystFullName[index_].Data() ), 100, -5, 5); h_sys_Top2Rapidity_[index_]->Sumw2(); 
   h_sys_Top2phi_[index_]      = new TH1F(Form("h_sys_Top2phi_%s"      ,v_SystFullName[index_].Data() ), Form("Top2 phi (%s)"      ,v_SystFullName[index_].Data() ), 24, -1*pi, pi); h_sys_Top2phi_[index_]->Sumw2(); 
   h_sys_Top2Energy_[index_]   = new TH1F(Form("h_sys_Top2Energy_%s"   ,v_SystFullName[index_].Data() ), Form("Top2 Energy (%s)"   ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_Top2Energy_[index_]->Sumw2(); 

   h_sys_TopMass_[index_]       = new TH1F(Form("h_sys_TopMass_%s"       ,v_SystFullName[index_].Data() ), Form("Top Mass (%s)"       ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_TopMass_[index_]->Sumw2(); 
   h_sys_Toppt_[index_]         = new TH1F(Form("h_sys_Toppt_%s"         ,v_SystFullName[index_].Data() ), Form("Top pt (%s)"         ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_Toppt_[index_]->Sumw2(); 
   h_sys_TopRapidity_[index_]   = new TH1F(Form("h_sys_TopRapidity_%s"   ,v_SystFullName[index_].Data() ), Form("Top Rapidity (%s)"   ,v_SystFullName[index_].Data() ), 100, -5, 5); h_sys_TopRapidity_[index_]->Sumw2(); 
   h_sys_Topphi_[index_]        = new TH1F(Form("h_sys_Topphi_%s"        ,v_SystFullName[index_].Data() ), Form("Top phi (%s)"        ,v_SystFullName[index_].Data() ), 24, -1*pi, pi); h_sys_Topphi_[index_]->Sumw2(); 
   h_sys_TopEnergy_[index_]     = new TH1F(Form("h_sys_TopEnergy_%s"     ,v_SystFullName[index_].Data() ), Form("Top Energy (%s)"     ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_TopEnergy_[index_]->Sumw2(); 
   h_sys_AnTopMass_[index_]     = new TH1F(Form("h_sys_AnTopMass_%s"     ,v_SystFullName[index_].Data() ), Form("AnTop Mass (%s)"     ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_AnTopMass_[index_]->Sumw2();
   h_sys_AnToppt_[index_]       = new TH1F(Form("h_sys_AnToppt_%s"       ,v_SystFullName[index_].Data() ), Form("AnTop pt (%s)"       ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_AnToppt_[index_]->Sumw2(); 
   h_sys_AnTopRapidity_[index_] = new TH1F(Form("h_sys_AnTopRapidity_%s" ,v_SystFullName[index_].Data() ), Form("AnTop Rapidity (%s)" ,v_SystFullName[index_].Data() ), 100, -5, 5); h_sys_AnTopRapidity_[index_]->Sumw2(); 
   h_sys_AnTopphi_[index_]      = new TH1F(Form("h_sys_AnTopphi_%s"      ,v_SystFullName[index_].Data() ), Form("AnTop phi (%s)"      ,v_SystFullName[index_].Data() ), 24, -1*pi, pi); h_sys_AnTopphi_[index_]->Sumw2(); 
   h_sys_AnTopEnergy_[index_]   = new TH1F(Form("h_sys_AnTopEnergy_%s"   ,v_SystFullName[index_].Data() ), Form("AnTop Energy (%s)"   ,v_SystFullName[index_].Data() ), 1000, 0.0, 1000); h_sys_AnTopEnergy_[index_]->Sumw2(); 

}

void ssb_analysis::End()
{
   fout->Write();
   fout->Close();
   /////////////////////////////////
   /// For All-in-One Systematic ///
   /////////////////////////////////
   if (isAllSyst == true)
   {
      for (int i = 0; i < v_SystFullName.size(); ++i )
      {
         a_fout[i]->Write();
         a_fout[i]->Close();
      }
   }
}

void ssb_analysis::SetInputFileName( char *inname )
{
   FileName_ = inname;
}
void ssb_analysis::SetOutputFileName(char *outname)
{   
   outfile = outname;
   TString rm_rootout = TString(outfile).ReplaceAll(".root","");
   v_outName.clear();
   if (isAllSyst == true)
   {
      for (int i =0; i < v_SystFullName.size(); ++i)
      { TString on = Form("%s_%s.root",rm_rootout.Data(),v_SystFullName[i].Data());  v_outName.push_back(on); }
   }
}

void ssb_analysis::GetTotalEvent()
{
   totalEvent = 0;
   cout << " FileName_ ? at Get TotalEvent () "<< FileName_ << endl;
   if ( TString( CenOfE ).Contains( "8TeV" ) )
   {
      if      ( TString(FileName_).Contains( "TTJets_FullLept" ) )     { totalEvent = SSBConfReader->GetNumber("TotalEvent", 1 );}
      else if ( TString(FileName_).Contains( "DYJetsToLL_M_10To50" ) ) { totalEvent = SSBConfReader->GetNumber("TotalEvent", 2 );}
      else if ( TString(FileName_).Contains( "DYJetsToLL_M_50" ) )     { totalEvent = SSBConfReader->GetNumber("TotalEvent", 3 );}
      else if ( TString(FileName_).Contains( "TTJets_HadronicMG" ) )   { totalEvent = SSBConfReader->GetNumber("TotalEvent", 4 );}
      else if ( TString(FileName_).Contains( "TTJets_SemiLeptMG" ) )   { totalEvent = SSBConfReader->GetNumber("TotalEvent", 5 );}
      else if ( TString(FileName_).Contains( "T_tW-channel" ) )        { totalEvent = SSBConfReader->GetNumber("TotalEvent", 6 );}
      else if ( TString(FileName_).Contains( "Tbar_tW-channel" ) )     { totalEvent = SSBConfReader->GetNumber("TotalEvent", 7 );}
      else if ( TString(FileName_).Contains( "WW" ) )                  { totalEvent = SSBConfReader->GetNumber("TotalEvent", 8 );}
      else if ( TString(FileName_).Contains( "WZ" ) )                  { totalEvent = SSBConfReader->GetNumber("TotalEvent", 9 );}
      else if ( TString(FileName_).Contains( "ZZ" ) )                  { totalEvent = SSBConfReader->GetNumber("TotalEvent", 10);}
      else if ( TString(FileName_).Contains( "WJetsToLNu" ) )          { totalEvent = SSBConfReader->GetNumber("TotalEvent", 11);}
      else if ( TString(FileName_).Contains( "Data" ) )                { totalEvent = SSBConfReader->GetNumber("TotalEvent", 12);}
      else { cout << " File Name Error !! at GetTotalEvent()  " << endl; }
   }
   else if ( TString( CenOfE ).Contains( "13TeV" ) )
   {
      cout << "CENOFE works well !! " << endl;
      if      ( TString(FileName_).Contains( "TTJets" ) )              { totalEvent = SSBConfReader->GetNumber("TotalEvent", 1 );}
      else if ( TString(FileName_).Contains( "WJetsToLNu" ) )          { totalEvent = SSBConfReader->GetNumber("TotalEvent", 2 );}
      else if ( TString(FileName_).Contains( "DYJetsToLL_M_10To50" ) ) { totalEvent = SSBConfReader->GetNumber("TotalEvent", 3 );}
      else if ( TString(FileName_).Contains( "DYJetsToLL_M_50" ) )     { totalEvent = SSBConfReader->GetNumber("TotalEvent", 4 );}
      else if ( TString(FileName_).Contains( "ST_tW_top" ) )           { totalEvent = SSBConfReader->GetNumber("TotalEvent", 5 );}
      else if ( TString(FileName_).Contains( "ST_tW_antitop" ) )       { totalEvent = SSBConfReader->GetNumber("TotalEvent", 6 );}
      else if ( TString(FileName_).Contains( "WW" ) )                  { totalEvent = SSBConfReader->GetNumber("TotalEvent", 7 );}
      else if ( TString(FileName_).Contains( "WZ" ) )                  { totalEvent = SSBConfReader->GetNumber("TotalEvent", 8 );}
      else if ( TString(FileName_).Contains( "ZZ" ) )                  { totalEvent = SSBConfReader->GetNumber("TotalEvent", 9 );}
      else if ( TString(FileName_).Contains( "Data" ) )                { totalEvent = SSBConfReader->GetNumber("TotalEvent", 10);}
      else { cout << " File Name Error !! at GetTotalEvent()  " << endl; }
   }
   else { cout << "CenOfE Error at GetTotalEvent() " << endl; }
   cout << "totalEvent ? " << totalEvent << endl;
}

// MC scale factor function
void ssb_analysis::MCSF()
{

//   double lumi = 19.6*1000;
//   double lumi = Lumi*1000;
//   double lumi = 40.24;
//   double lumi = 2.11*1000;
   double lumi = Lumi/1000000;
   double br   = 0.0159;
  // cout << "Lumi ? " << Lumi/1000000 << endl;
   cout << "lumi ? " << lumi<< endl;
   if ( TString(UsingTotEnv).Contains( "False" ) || TString(UsingTotEnv).Contains( "false" ) ) 
   {
      
      if ( TString( CenOfE ).Contains( "8TeV" ) ) // for 8 TeV
      {
         if      ( TString(FileName_).Contains( "TTJets_FullLept"     ) ) { mc_sf_ = 25.3*    lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "DYJetsToLL_M_10To50" ) ) { mc_sf_ = 860.5*   lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "DYJetsToLL_M_50"     ) ) { mc_sf_ = 3532.8*  lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "TTJets_HadronicMG"   ) ) { mc_sf_ = 106.9*   lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "TTJets_SemiLeptMG"   ) ) { mc_sf_ = 103*     lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "T_tW-channel"        ) ) { mc_sf_ = 11.2*    lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "Tbar_tW-channel"     ) ) { mc_sf_ = 11.2*    lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "WW"                  ) ) { mc_sf_ = 5.8*     lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "WZ"                  ) ) { mc_sf_ = 22.4*    lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "ZZ"                  ) ) { mc_sf_ = 9.0*     lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "WJetsToLNu"          ) ) { mc_sf_ = 37509.0* lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "Data"                ) ) { mc_sf_ = 1;                                }
         else { cout << " File Name Error !! at MCSF() " << endl; }
      }
      else if ( TString( CenOfE ).Contains( "13TeV" ) ) // for 13 TeV
      {
         if      ( TString(FileName_).Contains( "TTJets"              ) ) { mc_sf_ = 831.76*   lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "WJetsToLNu"          ) ) { mc_sf_ = 61526*    lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "DYJetsToLL_M_10To50" ) ) { mc_sf_ = 18610.0*    lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "DYJetsToLL_M_50"     ) ) { mc_sf_ = 6025*     lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "ST_tW_top"           ) ) { mc_sf_ = 35.6*     lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "ST_tW_antitop"       ) ) { mc_sf_ = 35.6*     lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "WW"                  ) ) { mc_sf_ = 118.7*    lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "WZ"                  ) ) { mc_sf_ = 65.9*     lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "ZZ"                  ) ) { mc_sf_ = 31.8*     lumi / NtupletotalEvent; }
         else if ( TString(FileName_).Contains( "Data"                ) ) { mc_sf_ = 1;                                 }
         else { cout << " File Name Error !! at MCSF()" << endl;}
      }
      else { cout << "Center Of Energy Error in MCSF()" << endl;  }
   }
   else if ( TString(UsingTotEnv).Contains( "True" ) || TString(UsingTotEnv).Contains( "true" ) )
   {
      if ( TString( CenOfE ).Contains( "8TeV" )  )
      {
         if      ( TString(FileName_).Contains( "TTJets_FullLept"     ) ) { mc_sf_ = 25.3*    lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "DYJetsToLL_M_10To50" ) ) { mc_sf_ = 860.5*   lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "DYJetsToLL_M_50"     ) ) { mc_sf_ = 3532.8*  lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "TTJets_HadronicMG"   ) ) { mc_sf_ = 106.9*   lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "TTJets_SemiLeptMG"   ) ) { mc_sf_ = 103*     lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "T_tW-channel"        ) ) { mc_sf_ = 11.2*    lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "Tbar_tW-channel"     ) ) { mc_sf_ = 11.2*    lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "WW"                  ) ) { mc_sf_ = 5.8*     lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "WZ"                  ) ) { mc_sf_ = 22.4*    lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "ZZ"                  ) ) { mc_sf_ = 9.0*     lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "WJetsToLNu"          ) ) { mc_sf_ = 37509.0* lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "Data"                ) ) { mc_sf_ = 1;                          }
         else { cout << " File Name Error !! " << endl; }
      }
      else if ( TString(CenOfE).Contains("13TeV") )
      {
         if      ( TString(FileName_).Contains( "TTJets"              ) ) { mc_sf_ = 831.76*   lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "WJetsToLNu"          ) ) { mc_sf_ = 61526*    lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "DYJetsToLL_M_10To50" ) ) { mc_sf_ = 18610.0*    lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "DYJetsToLL_M_50"     ) ) { mc_sf_ = 6025*     lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "ST_tW_top"           ) ) { mc_sf_ = 35.6*     lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "ST_tW_antitop"       ) ) { mc_sf_ = 35.6*     lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "WW"                  ) ) { mc_sf_ = 118.7*    lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "WZ"                  ) ) { mc_sf_ = 65.9*     lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "ZZ"                  ) ) { mc_sf_ = 31.8*     lumi / totalEvent; }
         else if ( TString(FileName_).Contains( "Data"                ) ) { mc_sf_ = 1;                                 }
         else { cout << " File Name Error !! at MCSF()" << endl;}

      }
      else { cout << "Center Of Energy Error in MCSF()" << endl; }
   } 
   else {cout << "MCSF error " << endl; }
}

// Apply MC SF To Event //
void ssb_analysis::MCSFApply()
{
   evt_weight_beforemcsf_ =1; // Initailize evt_weight_beforemcsf_ //
   evt_weight_beforemcsf_ = evt_weight_; // keep event weight // 
   if ( !TString(FileName_).Contains( "Data") ){ evt_weight_ = evt_weight_*mc_sf_; } // apply MC scale factor // 
   else {evt_weight_ = 1;}
   if ( isAllSyst == true )
   {
      for( int i =0; i < v_SystFullName.size(); ++i )
      { 
         if ( !TString(FileName_).Contains( "Data") )
         {
            v_SystEvt[i]=evt_weight_beforemcsf_*mc_sf_;
         }
         else {
            v_SystEvt[i]=1;
         }
         m_Syst_EvtW[ v_SystFullName[i] ] = v_SystEvt[i];
      } // apply MC scale factor //
   }
}
void ssb_analysis::GenWeightApply()
{
   double genweight = 1.0;
   evt_weight_beforegenweight_ =1; // Initailize evt_weight_beforegenweight_ //
   evt_weight_beforegenweight_ = evt_weight_; // keep event weight //
   if ( !TString(FileName_).Contains( "Data") ){
      if (Gen_EventWeight > 0.0){genweight =1;}
      else {genweight =-1;} 
      evt_weight_ = evt_weight_*genweight;
   }
   else {evt_weight_ = 1;}
   if ( isAllSyst == true )
   {
      for( int i =0; i < v_SystFullName.size(); ++i )
      { 
         if ( !TString(FileName_).Contains( "Data") )
         {
            v_SystEvt[i]=evt_weight_beforegenweight_*genweight;
         }
         else {
            v_SystEvt[i]=1;
         }
         m_Syst_EvtW[ v_SystFullName[i] ] = v_SystEvt[i];
      } // apply MC scale factor //
   }
}
void ssb_analysis::PDFWeightApply()
{
/*   double pdfweight = 1.0;
   evt_weight_beforepdfweight_ = evt_weight_;
   if ( !TString(FileName_).Contains("Data") )
   {
      if (PDFSys == -1 ){ pdfweight =1; }
      else if ( PDFSys == 0 ) { pdfweight = PDFWeight_Cent->at(1); cout << "--- " << pdfweight << endl;}
      else if ( PDFSys > 0 )  
      {
         if      ( PDFSys%2 == 1 ) { pdfweight = PDFWeight_Var1_Up->at(PDFSys/2); } //odd number is PDFUp
         else if ( PDFSys%2 == 0 ) { pdfweight = PDFWeight_Var1_Down->at(PDFSys/2 -1); } //even number is PDFDown
         else {cout << "Somethig Wrong in ther PDF sys ... " << endl;}
      }
      evt_weight_ = evt_weight_*pdfweight;
   }
   else { evt_weight_ = 1; }*/
}
// Apply Trigger SF To Event //
void ssb_analysis::TriggerSFApply()
{
   double triggersf_ =1;
   evt_weight_beforeTrigger_ =1; // Initailize evt_weight_beforeTrigger_ //
   evt_weight_beforeTrigger_ = evt_weight_; // keep event weight // 
   if ( !TString(FileName_).Contains( "Data") )
    {  //evt_weight_ = evt_weight_*triggersf_;
/*      if ( TString(Decaymode).Contains( "dielec" ) ) { triggersf_ = 0.966; }
      else if (TString(Decaymode).Contains( "muel" )  )  { triggersf_ = 0.962;}
      else if (TString(Decaymode).Contains( "dimu" )  )  { triggersf_ = 0.980;}*/
      if ( TString(Decaymode).Contains( "dielec" ) ) { triggersf_ =SSBEffcal->TrigDiElec_Eff(Lep1,Lep2); }
      else if (TString(Decaymode).Contains( "muel" )  )  { triggersf_ = SSBEffcal->TrigMuElec_Eff(TElectron,TMuon);
//      cout << "muon " << TMuon->Pt() << " elec " << TElectron->Pt() << endl;
      }
      else if (TString(Decaymode).Contains( "dimu" )  )  { triggersf_ = SSBEffcal->TrigDiMuon_Eff(Lep1,Lep2);}
      else { cout << "CHECK OUT TRIGGER EFF"<< endl;}
      evt_weight_ = evt_weight_*triggersf_;
   }// apply Trigger scale factor //
   else { evt_weight_ = evt_weight_; } 
}
void ssb_analysis::PileUpReWeightApply()
{
   evt_weight_beforePileup_ = 1;
   evt_weight_beforePileup_ = evt_weight_; // keep event weight // 
   double puweight_ = 1.;
   double pu_weight_central = puweight->weight( PileUp_Count_Intime );
   double pu_weight_up = puweightup->weight( PileUp_Count_Intime );
   double pu_weight_dn = puweightdn->weight( PileUp_Count_Intime );

   if ( !TString(FileName_).Contains( "Data") )
   {
//      puweight_ = puweight->weight( PileUp_Count_Intime );
      if (TString(PileUpSys).Contains("central") ) { puweight_   = pu_weight_central;}
      else if (TString(PileUpSys).Contains("up") ) { puweight_   = pu_weight_up; }
      else if (TString(PileUpSys).Contains("down") ) { puweight_ = pu_weight_dn; }
      else {  
      cout << "PileUp sys Error ... Defalut is Weight_PileUp ... : " << PileUpSys << endl;
   }
   evt_weight_ = evt_weight_*puweight_; } // apply PileUpReweight //
   else {evt_weight_ = 1;}
   if ( isAllSyst == true )
   {
      for( int i =0; i < v_SystFullName.size(); ++i )
      {
         if( !TString(FileName_).Contains("Data") )
         {
            if ( v_SystFullName[i].Contains("PileUpUp"))
            {
               v_SystEvt[i] = v_SystEvt[i]*pu_weight_up;
            }
            else if ( v_SystFullName[i].Contains("PileUpDown"))
            {
               v_SystEvt[i] = v_SystEvt[i]*pu_weight_dn;
            }
            else {
               v_SystEvt[i]=v_SystEvt[i]*pu_weight_central;
            }
         }
         else { // Not for Data // 
            v_SystEvt[i]=1;
         }
         m_Syst_EvtW[ v_SystFullName[i] ] = v_SystEvt[i];
      }
   }
}
void ssb_analysis::Weight()
{
   double puweight_ = 1;
   double k_fac_;
   double genweight = 1;
   double pdfweight = 1;

/*   /////////////////////////
   /// To Get PDF Weight ///
   /////////////////////////
   if ( !TString(FileName_).Contains("Data") )
   {
      if (PDFSys == -1 ){ pdfweight =1; }
      else if ( PDFSys == 0 ) { pdfweight = PDFWeight_Cent->at(1); }
      else if ( PDFSys > 0 )  
      {
         if      ( PDFSys%2 == 1 ) { pdfweight = PDFWeight_Var1_Up->at(PDFSys/2); } //odd number is PDFUp
         else if ( PDFSys%2 == 0 ) { pdfweight = PDFWeight_Var1_Down->at(PDFSys/2 -1); } //even number is PDFDown
         else {cout << "Somethig Wrong in ther PDF sys ... " << endl;}
      }
   }*/
   evt_weight_ = 1;
   
   ///////////////////////
   // PileUp Systematic //
   ///////////////////////
   puweight_ = puweight->weight( PileUp_Count_Intime );
   if (TString(PileUpSys).Contains("central") ) { puweight_   = puweight->weight( PileUp_Count_Intime ); }
   else if (TString(PileUpSys).Contains("up") ) { puweight_   = puweightup->weight( PileUp_Count_Intime ); }
   else if (TString(PileUpSys).Contains("down") ) { puweight_ = puweightdn->weight( PileUp_Count_Intime ); }
   else {  
   cout << "PileUp sys Error ... Defalut is Weight_PileUp ... : " << PileUpSys << endl;
   }
   if ( TString(FileName_).Contains( "TTJets"   ) )
   {  
      k_fac_ = 1;
      evt_weight_ = k_fac_*mc_sf_*puweight_*lep_eff*pdfweight;
   }
//   { k_fac_ = 252.89/235.2; evt_weight_ = k_fac_*mc_sf_; }
   else if ( TString(FileName_).Contains( "DYJetsToLL" ) || 
             TString(FileName_).Contains( "WW"         ) ||
             TString(FileName_).Contains( "WZ"         ) ||
             TString(FileName_).Contains( "ZZ"         ) ||
             TString(FileName_).Contains( "WJetsToLNu" ) ||
             TString(FileName_).Contains( "tW" )   )
   {
      k_fac_ =1; 
      if ( TString(FileName_).Contains( "DYJetsToLL" ) || TString(FileName_).Contains( "WJetsToLNu" ) ) 
      {
         if (Gen_EventWeight > 0){genweight =1;}
         else {genweight =-1;}
      }
      evt_weight_ = k_fac_*mc_sf_*puweight_*lep_eff*genweight*pdfweight; 
   }
   else if ( TString(FileName_).Contains( "Data" ) )
   {
      k_fac_ =1; evt_weight_ = 1; puweight_=1; 
   }
   else {cout << "Event-Weight error !!" << endl; evt_weight_ =0; }
}
void ssb_analysis::NumPVCount()
{
   num_pv = 0;
   //cout << " start NumPVCount "<< endl;
   for (int i =0; i < Filter_PV->size(); ++i)
   {
      //cout <<  " Filter_PV->at(i) : " <<  Filter_PV->at(i) << endl;
      if ( Filter_PV->at(i) == true ) {num_pv++;}
   }
}

bool ssb_analysis::EveRun()
{
   bool everun = false;
   if (TString(FileName_).Contains( "Data" ) ) // Only Data
   {
      if ( TString(FileName_).Contains( "Jul" ) ) //Jul-Data
      {
         if ( TString(Decaymode).Contains( "dimuon" ) || TString(Decaymode).Contains( "dielec" ) || TString(Decaymode).Contains( "muel" )  )
         {
            if (Info_RunNumber >= 246908 && Info_RunNumber < 251604 ) { everun = true; }
         }
         else { cout << "YOU DON'T USE THIS FUNCTION , CHECH DECAY MODE" << endl;}
      }
      else 
      {
         if ( TString(Decaymode).Contains( "dimuon" ) || TString(Decaymode).Contains( "dielec" ) || TString(Decaymode).Contains( "muel" ) )
         {
            if (Info_RunNumber >= 251604 && Info_RunNumber <= 251883 ) { everun = true; }
         }
         else { cout << "YOU DON'T USE THIS FUNCTION , CHECH DECAY MODE" << endl;}
      }//prompt
   }
   else { everun = true;} //MC
   return everun;
}
bool ssb_analysis::ChannelIndex()
{
   bool chan_cut = true;
   int decay_idx = 0;
   if ( TString(Decaymode).Contains( "dielec" ) ) decay_idx = 22;
   if ( TString(Decaymode).Contains( "muel" ) )   decay_idx = 24;
   if ( TString(Decaymode).Contains( "dimuon" ) ) decay_idx = 26;

   if ( Channel_Idx != decay_idx && ( TString(FileName_).Contains("other") == false ) && TString(FileName_).Contains("TTJets") ) {chan_cut=false;} // only TTbar -> Di-Lepton Signal.
   if ( Channel_Idx == decay_idx && TString(FileName_).Contains("other") ) {chan_cut = false;} // only TTbar -> Di-Lepton Signal.
   return chan_cut;
}
bool ssb_analysis::METFilterAPP()
{
   bool metfilt_ = true;
   for ( int i = 0; i < METFilter_Name->size(); ++i )
   { 
      TString METFiltName = METFilter_Name->at(i);
      //cout << "METFiltName : " << METFilter_Name->at(i) << endl;
      if ( !(TString(FileName_).Contains("Data")) && 
            TString(METFiltName).Contains("Flag_eeBadScFilter")  ) {continue;}
      for( int j =0; j < v_METFilterName.size(); ++j )
      {
         if ( TString( METFiltName ).Contains(v_METFilterName.at(j)) )
         {
            if ( !(METFilter_isPass->at(i)) ){ metfilt_ = false; }
         } 
      }
   }
   for ( int i = 0; i < METFilterAdd_Name->size(); ++i )
   {
      TString METFiltName = METFilterAdd_Name->at(i);
      if ( !(TString(FileName_).Contains("Data")) && 
            TString(METFiltName).Contains("Flag_eeBadScFilter")  ) {continue;}
      for( int j =0; j < v_METFilterName.size(); ++j )
      {
         if ( TString( METFiltName ).Contains(v_METFilterName.at(j)) )
         {
            if ( !( METFilterAdd_isPass->at(i) ) ){ metfilt_ = false; }
         } 
      }
   }
   return metfilt_;
   //return true;
}

// Trigger Requirement Function
bool ssb_analysis::Trigger()
{
   /// Variable for Trigger Function
   int ptrigindex;
   bool trigpass;

   ptrigindex =0;
   trigpass = false;
   TString trgName = ""; 
//   if (!TString(FileName_).Contains( "Data" ) ) { trigpass = true; return trigpass;}
   for (int i =0; i < Trigger_Name->size(); i++)
   {
//      cout << "Ntuple Triggers : " <<  Trigger_Name->at(i) << endl; 
      for (int j = 0; j < trigName.size(); j++)
      {
         trgName = trigName[j];
         if (TString(Decaymode).Contains( "dimuon" ) ) {
            if ( TString(FileName_).Contains( "Run2016H" ) ) 
            { 
               if (trgName.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v") ||
                   trgName.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v") ) {continue;}
            }
            else {
               if (trgName.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ||
                   trgName.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") ) {continue;}
            }
         }
         if (TString(Decaymode).Contains( "muel" ) ) {
            if ( TString(FileName_).Contains( "Run2016H" ) ) 
            { 
               if (trgName.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") ||
                   trgName.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ) {continue;}
            }
            else {
               if (trgName.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") ||
                   trgName.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) {continue;}
            }
         }

         if ( TString( Trigger_Name->at(i) ).Contains( trigName.at(j) ) )//TString clone 
         {
//            cout << "trigName.at(j) ?" << trigName.at(j) << endl;
            if ( ( Trigger_isPass->at(i)  ) && 
                !( Trigger_isError->at(i) ) && 
                 ( Trigger_isRun->at(i) )      ) 
            { ptrigindex = ptrigindex+1; }
         }
      }
   }
   if ( ptrigindex > 0 ) { trigpass = true; }
   return trigpass;
}
// Function for Lepton Selection
void ssb_analysis::LeptonSelector()
{
   v_lepton_idx.clear();
   v_muon_idx.clear();
   v_electron_idx.clear();

   // dummy vector //
   v_lep_idx_temp.clear();
   v_muon_idx_temp.clear();
   v_electron_idx_temp.clear();

   // Case of Dimuon //
   if (TString(Decaymode).Contains("dimuon")) 
   { 
      Int_t nmu = Muon->GetEntriesFast();
      for (int i =0; i < nmu; ++i )
      {
         Mu_lep_sel = (TLorentzVector*)Muon->At(i);
         if ( ( v_lepton_iso->at(i) > lepisocut )   ||
              (*Mu_lep_sel).Pt()    < lep_pt        || 
              fabs( (*Mu_lep_sel).Eta() ) > lep_eta || 
              ( v_lepton_Id->at(i) == false )          ){continue;}         
         v_lep_idx_temp.push_back( i );

         if (v_lep_idx_temp.size()==1 && v_lepton_idx.size() == 0 ){  v_lepton_idx.push_back(i);}
         if (v_lep_idx_temp.size()>0)
         {
            if ( Muon_Charge->at( v_lep_idx_temp.at(0) ) != Muon_Charge->at( i ) ) { v_lepton_idx.push_back(i); }
         }
      }
      
   }
   // Case of Dielectron //
   else if (TString(Decaymode).Contains("dielec")) 
   { 
      Int_t nel = Elec->GetEntriesFast();

      for (int i =0; i < nel; ++i )
      {
         Ele_lep_sel = (TLorentzVector*)Elec->At(i);
         if ( (*Ele_lep_sel).Pt() != (*Ele_lep_sel).Pt() ){continue;}
         if ( (*Ele_lep_sel).Pt() < lep_pt )          {continue;}
         if ( fabs( (*Ele_lep_sel).Eta() ) > lep_eta ){continue;}   

         // Requirement of Electron Isolation with Electron Eta 
         //if (Elec_SCB_Medium->at(i) == false) {continue;}    
         //if (v_electron_Id->at(i) == false) {continue;}    
         if (v_lepton_Id->at(i) == false) {continue;}    
         if ( fabs( Elec_Supercluster_Eta->at(i) ) > 1.4442 && fabs( Elec_Supercluster_Eta->at(i) ) < 1.566 ) {continue;} 
         if ( Elec_ChargeId_GsfPx->at(i) == false ){continue;}//
         if ( fabs( Elec_Supercluster_Eta->at(i) ) <= 1.479 ) 
         { 
            if (v_lepton_iso->at(i) > 0.0695 ) {continue;}
            if ( fabs( Elec_Track_GsfdZ->at(i) ) > 0.10 || fabs( Elec_Track_GsfdXY->at(i) ) > 0.05 )  {continue;} 
         }
         if ( fabs( Elec_Supercluster_Eta->at(i) ) > 1.479 )  
         {  
            if (v_lepton_iso->at(i) > 0.0821 ) {continue;}
            if (fabs(Elec_Track_GsfdZ->at(i)) > 0.20 || fabs(Elec_Track_GsfdXY->at(i)) > 0.10 )  {continue;} 
         }

         v_lep_idx_temp.push_back( i );

         if ( v_lep_idx_temp.size()==1 && v_lepton_idx.size() == 0 ){  v_lepton_idx.push_back(i); }
         if ( v_lep_idx_temp.size() > 0 )
         {
            if ( Elec_Charge->at( v_lep_idx_temp.at(0) ) != Elec_Charge->at( i ) ) { v_lepton_idx.push_back(i); }
         }

      }
   }
   // Case of MuEl //
   else if (TString(Decaymode).Contains("muel"))
   {

      Int_t nmu = Muon->GetEntriesFast();
      for (int i =0; i < nmu; ++i )
      {
         Mu_lep_sel = (TLorentzVector*)Muon->At(i);
         if ( ( v_mlep_iso->at(i) > muonisocut ) || (*Mu_lep_sel).Pt() < muon_pt || fabs( (*Mu_lep_sel).Eta() ) > muon_eta || v_muon_Id->at(i) == false  ){continue;}         
         //if ( ( v_mlep_iso->at(i) > muonisocut ) || (*Mu_lep_sel).Pt() < 35 || fabs( (*Mu_lep_sel).Eta() ) > muon_eta || v_muon_Id->at(i) == false  ){continue;}         
         v_muon_idx.push_back( i );
      }

      Int_t nel = Elec->GetEntriesFast();
      for (int i =0; i < nel; ++i )
      {
         Ele_lep_sel = (TLorentzVector*)Elec->At(i);

         ////
         if ( (*Ele_lep_sel).Pt() < elec_pt ){continue;}
         //if ( (*Ele_lep_sel).Pt() < 35 ){continue;}

         if ( fabs( (*Ele_lep_sel).Eta() ) > elec_eta ){continue;}

         //if (Elec_SCB_Medium->at(i) == false) {continue;}    
         if (v_electron_Id->at(i) == false) {continue;}    
         if ( Elec_Conversion->at(i) == false ){continue;}
         if ( fabs( Elec_Supercluster_Eta->at(i) ) > 1.4442 && fabs( Elec_Supercluster_Eta->at(i) ) < 1.566 ) {continue;}
         if ( Elec_ChargeId_GsfPx->at(i) == false ){continue;}        
         if ( fabs( Elec_Supercluster_Eta->at(i) ) <= 1.479 ) 
         { 
            if ( fabs( Elec_Track_GsfdZ->at(i) ) > 0.10 || fabs( Elec_Track_GsfdXY->at(i) ) > 0.05 )  {continue;} 
         }
         if ( fabs( Elec_Supercluster_Eta->at(i) ) > 1.479 )  
         {  
            if (fabs(Elec_Track_GsfdZ->at(i)) > 0.20 || fabs(Elec_Track_GsfdXY->at(i)) > 0.10 )  {continue;} 
         }
         v_electron_idx_temp.push_back( i );

      }

      if (v_muon_idx.size() > 0 && v_electron_idx_temp.size() > 0 )
      {  
         for ( int i = 0; i < v_electron_idx_temp.size(); ++i )
         {
            if ( Muon_Charge->at( v_muon_idx.at(0) ) != Elec_Charge->at( v_electron_idx_temp.at(i) ) ) { v_electron_idx.push_back(v_electron_idx_temp.at(i) ); }
         }
      }
      
      if ( v_electron_idx.size() > 0 && v_muon_idx.size() > 0 ) 
      { 
         if ( Muon_Charge->at( v_muon_idx.at(0) ) == Elec_Charge->at( v_electron_idx.at(0) ) ) { cout << "MuEl seclection wrong !! " << endl; } 
      }

   }
   else { cout << "Lepton Selection error" << endl;}
}
bool ssb_analysis::NumIsoLeptons()//YOU SHOULD REQUIRE THIS FUNCTION AFTER LEPTONSELETOR //
{
   bool numLeptons = true;
   if ( TString(Decaymode).Contains( "dimuon" ) ||
        TString(Decaymode).Contains( "dielec" )    ){  if ( v_lepton_idx.size() <= 1 ) {numLeptons=false;} } //Requiring  Two isolated Lepton .

   else if ( TString(Decaymode).Contains( "muel" ) ){ if (v_muon_idx.size() < 1 || v_electron_idx.size() < 1 ) {numLeptons=false; } else {}}
   else if ( TString(Decaymode).Contains( "muonJet" ) ){if( v_lepton_idx.size() != 1){numLeptons=false;} }
   else { cout << "?? something wrong " << endl; }
   return numLeptons;
}
bool ssb_analysis::LeptonsPtAddtional()//YOU SHOULD REQUIRE THIS FUNCTION AFTER NumIsoLeptons //
{
   bool lepptadd = false;
   if ( TString(Decaymode).Contains( "dimuon" ) ||
        TString(Decaymode).Contains( "dielec" ) ||
        TString(Decaymode).Contains( "muel" )      ){  if ( Lep1->Pt() > 25 && Lep2->Pt() > 20 ) {lepptadd=true;} } //Pt of Leading Lepton should be over than 25 GeV and Second Leading Lepton Pt should be over thand 20 GeV.

   else if ( TString(Decaymode).Contains( "muonJet" ) ){if( v_lepton_idx.size() == 1){lepptadd=true;} }
   else { cout << "?? something wrong " << endl; }  
   return lepptadd;
}
bool ssb_analysis::MuVeto()
{
   bool muveto;
   muveto = false;
   int vetolep = 0;

   if ( TString(Decaymode).Contains( "dimuon" ) )
   {
      for ( int ilumuon = 0; ilumuon < Muon_Count; ++ilumuon )
      {
         if ( ilumuon != v_lepton_idx.at(0) && ilumuon != v_lepton_idx.at(1) )
         {
            if ( Muon_isLoose->at( ilumuon ) == true && Muon_PFIsodBeta04->at( ilumuon ) < 0.25 )
            {   
               vetoMu = (TLorentzVector*)Muon->At(ilumuon);
               if ( vetoMu->Pt() > 20 && fabs(vetoMu->Eta()) < 2.4 ) { vetolep++; }
            }
         }
      }

   }
   else if ( TString(Decaymode).Contains( "dielec" ) )
   {
      // 3rd muon veto
      TLorentzVector *vetoMu = new TLorentzVector();
      for ( int ilumuon = 0; ilumuon < Muon_Count; ++ilumuon )
      {
         if ( Muon_isLoose->at( ilumuon ) == true && Muon_PFIsodBeta04->at( ilumuon ) < 0.25 )
         {   
            vetoMu = (TLorentzVector*)Muon->At(ilumuon);
            if ( vetoMu->Pt() > 20 && fabs(vetoMu->Eta()) < 2.4 ) { vetolep++; }
         }
      }
   }
   else if ( TString(Decaymode).Contains( "muel" ) )
   {
      for ( int ilumuon = 0; ilumuon < Muon_Count; ++ilumuon )
      {        
         if ( ilumuon != v_muon_idx.at(0) )
         {             
            if ( Muon_isLoose->at( ilumuon ) == true && Muon_PFIsodBeta04->at( ilumuon ) < 0.25 )
            {   //vetolep++;
               vetoMu = (TLorentzVector*)Muon->At(ilumuon);
               if ( vetoMu->Pt() > 20 && fabs(vetoMu->Eta()) < 2.4 ) { vetolep++; }
            }  
         }        
      }
   }
   else { cout << "Something Wrong At MuVeto. Check Decaymode-" << endl;}

   if ( vetolep == 0 ){ muveto = true; }
   return muveto;
}
bool ssb_analysis::ElVeto()
{
   bool elveto;
   elveto = false;
   int vetolep = 0;
   if ( TString(Decaymode).Contains( "dimuon" ) )// case of dimuon ...
   {
      for ( int iele = 0; iele < Elec_Count; ++iele )
      {
         vetoEle = (TLorentzVector*)Elec->At(iele);
         if(  vetoEle->Pt() != vetoEle->Pt()) {continue;}
         if (
             Elec_SCB_Veto->at(iele) == true &&
             vetoEle->Pt() > 20                 &&
             fabs( vetoEle->Eta() ) < 2.5       &&
             Elec_Conversion->at(iele) == true     )
         {
          vetolep++; 
         }

      }
   } // dimuon requriement...
   else if( TString(Decaymode).Contains( "dielec" ) ) 
   {
      // 3rd Electron veto ..
      for ( int iele = 0; iele < Elec_Count; ++iele )
      {
         if ( iele != v_lepton_idx.at(0) && iele != v_lepton_idx.at(1) )
         {
            vetoEle = (TLorentzVector*)Elec->At(iele);
            if (
                 Elec_SCB_Veto->at(iele) == true &&
                 vetoEle->Pt()            > 20      &&
                 fabs( vetoEle->Eta() )   < 2.5     &&
                 Elec_Conversion->at(iele) == true     )
            { 
               if ( fabs( Elec_Supercluster_Eta->at(iele) ) <= 1.479 )
               {
                  if (
                       v_el_iso_jcl->at( iele ) < 0.175 &&
                       fabs( Elec_Track_GsfdXY->at(iele) ) < 0.05 &&
                       fabs( Elec_Track_GsfdZ->at(iele) ) < 0.10 
                     ){ vetolep++; }
               
               } 
               else 
               {
                  if (
                       v_el_iso_jcl->at( iele ) < 0.159 &&           
                       fabs( Elec_Track_GsfdXY->at(iele) ) < 0.10 &&
                       fabs( Elec_Track_GsfdXY->at(iele) ) < 0.20 
                     ){ vetolep++; }
               }
            }
         }
      }

   } // case of dielectron...

   else if ( TString(Decaymode).Contains( "muel" ) ) 
   {
      for ( int iele = 0; iele < Elec_Count; ++iele )
      {
         if ( iele != v_electron_idx.at(0) )
         {
            vetoEle = (TLorentzVector*)Elec->At(iele);
            if(  vetoEle->Pt() != vetoEle->Pt()) {continue;}

            if (
                 Elec_SCB_Veto->at(iele) == true &&
                 vetoEle->Pt()            > 20   &&
                 fabs( vetoEle->Eta() )   < 2.5  &&
                 Elec_Conversion->at(iele) == true 
               )
            {
               if ( fabs( Elec_Supercluster_Eta->at(iele) ) <= 1.479 )
               {
                  if (
                       v_el_iso_jcl->at( iele ) < 0.175 &&
                       fabs( Elec_Track_GsfdXY->at(iele) ) < 0.05 &&
                       fabs( Elec_Track_GsfdZ->at(iele) ) < 0.10 
                     )
                  { vetolep++; }
               
               } 
               else 
               {
                  if (
                       v_el_iso_jcl->at( iele ) < 0.159 &&           
                       fabs( Elec_Track_GsfdXY->at(iele) ) < 0.10 &&
                       fabs( Elec_Track_GsfdXY->at(iele) ) < 0.20 
                     )
                  { vetolep++; }
               }              
            }
         }
      }
   }
   else {cout << "EleVto check !!! " << endl;}
   if ( vetolep == 0 ){ elveto = true; }

   return elveto;
}
bool ssb_analysis::ThirdLeptonVeto()
{
   bool third_veto = true;
   if ( TString(Decaymode).Contains( "dimuon" ) || TString(Decaymode).Contains( "dielec" )    )
   { 
      if ( v_lepton_idx.size() <= 1 ) {third_veto = false;;}
      //if ( v_lepton_idx.size() != 2 ) {third_veto = false;;}

      if( TString(Decaymode).Contains( "dimuon" ) ) 
      { 
         if ( Muon_Charge->at( v_lepton_idx.at(0) ) == Muon_Charge->at( v_lepton_idx.at(1) ) ) {third_veto = false;;}
         if ( MuVeto() == false ) {third_veto = false;;}
         if ( ElVeto() == false ) {third_veto = false;;}
      } 
      else if( TString(Decaymode).Contains( "dielec" ) ) 
      {
         if ( Elec_Charge->at( v_lepton_idx.at(0) ) == Elec_Charge->at( v_lepton_idx.at(1) ) ) {third_veto = false;;}
         if ( MuVeto() == false ) {third_veto = false;;}
         if ( ElVeto() == false ) {third_veto = false;;}
      }
   } //Requiring  Two isolated Lepton .

   else if ( TString(Decaymode).Contains( "muel" ) )
   { 
      if (v_muon_idx.size() < 1 || v_electron_idx.size() < 1 ) { third_veto = false;; } 
      if ( Muon_Charge->at( v_muon_idx.at(0) ) == Elec_Charge->at( v_electron_idx.at(0) ) ) { third_veto = false;; }
      if ( MuVeto() == false ) {third_veto = false;;}
      if ( ElVeto() == false ) {third_veto = false;;}
   }

   else { cout << "something wrong at ThirdLeptonVeto veto Function !!!!" << endl; }
   return third_veto;
}
void ssb_analysis::LeptonOrder()
{
   if ( TString(Decaymode).Contains( "dimuon" ) )
   {
      if (v_lepton_idx.size() > 1)
      {
         Lep1 = (TLorentzVector*)Muon->At( v_lepton_idx.at(0) );
         Lep2 = (TLorentzVector*)Muon->At( v_lepton_idx.at(1) );
         if (Muon_Charge->at( v_lepton_idx.at(0) ) < 0 )
         {
            Lep   = (TLorentzVector*)Muon->At( v_lepton_idx.at(0) );
            AnLep = (TLorentzVector*)Muon->At( v_lepton_idx.at(1) );
         }
         else
         {
            Lep   = (TLorentzVector*)Muon->At( v_lepton_idx.at(1) );
            AnLep = (TLorentzVector*)Muon->At( v_lepton_idx.at(0) );
         }
      }
   }
 
   else if ( TString(Decaymode).Contains( "dielec" ) )
   {
      if (v_lepton_idx.size() > 1)
      {
         Lep1 = (TLorentzVector*)Elec->At( v_lepton_idx.at(0) );
         Lep2 = (TLorentzVector*)Elec->At( v_lepton_idx.at(1) );
         if (Elec_Charge->at( v_lepton_idx.at(0) ) < 0 )
         {
            Lep   = (TLorentzVector*)Elec->At( v_lepton_idx.at(0) );
            AnLep = (TLorentzVector*)Elec->At( v_lepton_idx.at(1) );
         }
         else
         {
            Lep   = (TLorentzVector*)Elec->At( v_lepton_idx.at(1) );
            AnLep = (TLorentzVector*)Elec->At( v_lepton_idx.at(0) );
         }
      }
   }
   else if ( TString(Decaymode).Contains( "muel" ) )
   {
      if ( v_muon_idx.size() > 0 && v_electron_idx.size() > 0 )
      {
         TMuon = (TLorentzVector*)Muon->At( v_muon_idx.at(0) );
         TElectron = (TLorentzVector*)Elec->At( v_electron_idx.at(0) );
         if ((*TMuon).Pt() > (*TElectron).Pt())
         {
            Lep1 = (TLorentzVector*)Muon->At( v_muon_idx.at(0) );
            Lep2 = (TLorentzVector*)Elec->At( v_electron_idx.at(0) );
         }
         else
         {
            Lep1 = (TLorentzVector*)Elec->At( v_electron_idx.at(0) );
            Lep2 = (TLorentzVector*)Muon->At( v_muon_idx.at(0) );
         }
         if ( Muon_Charge->at( v_muon_idx.at(0) ) < 0 )
         {
            Lep   = (TLorentzVector*)Muon->At( v_muon_idx.at(0) );
            AnLep = (TLorentzVector*)Elec->At( v_electron_idx.at(0) );
         }
         else 
         {
            AnLep = (TLorentzVector*)Muon->At( v_muon_idx.at(0) );
            Lep   = (TLorentzVector*)Elec->At( v_electron_idx.at(0) );
         }
      }
   }
   else {cout << "Lepton TLorentzVector Error" << endl;}
}

void ssb_analysis::LeptonGetSF()
{

   lep_eff = 1.0;
   if ( TString(Decaymode).Contains( "dimuon" ) )
   {
      lep_eff = SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,LepIdSFSys,LepIsoSFSys,LepTrackSFSys); 
   }
   else if ( TString(Decaymode).Contains( "dielec" ) )
   {
      //lep_eff = SSBEffcal->DoubleElec_EffROOT(Lep1,Lep2,Elec_Supercluster_Eta->at(v_lepton_idx[0]) ,Elec_Supercluster_Eta->at(v_lepton_idx[1]) ,LepIdSFSys,LepIsoSFSys);
      lep_eff = SSBEffcal->DoubleElec_EffROOT(Lep1,Lep2,Elec_Supercluster_Eta->at(v_lepton_idx[0]) ,Elec_Supercluster_Eta->at(v_lepton_idx[1]) ,LepIdSFSys,LepRecoSFSys);// LepRecoSFSys is for Lepton //
   }
   else if ( TString(Decaymode).Contains( "muel" ) )
   {
      //lep_eff = 1.0;
      lep_eff = SSBEffcal->MuonElec_EffROOT(TMuon,TElectron,Elec_Supercluster_Eta->at(v_electron_idx[0]) ,LepIdSFSys,LepIsoSFSys,LepTrackSFSys,LepRecoSFSys);//LepRecoSFSys is for Electron
   }
   else { lep_eff = 1.0; cout << "LeptonGetSF error !!!!" << endl; }
}
void ssb_analysis::LeptonSFApply()
{

   double lep_sf_cent = 0.0;

   evt_weight_beforeLepsf_ = 1.0;
   evt_weight_beforeLepsf_ = evt_weight_; // keep event weight 
   if ( !TString(FileName_).Contains( "Data") )
   {
      LeptonGetSF(); // Getting Lepton SF //
      evt_weight_ = evt_weight_*lep_eff;
   }
   else {evt_weight_ = 1;}

   /// For All-in-One Systematic ///
   if ( isAllSyst == true )
   {
      /// Get Central Lepton SF for each channel
      if ( !TString(FileName_).Contains( "Data") )
      {
         if ( TString(Decaymode).Contains( "dimuon" ) )
         {
            lep_sf_cent = SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,"central","central","central"); // lep1 lep2 idsys isosys tracksys 
         }
         else if ( TString(Decaymode).Contains( "dielec" ) )
         {
            lep_sf_cent = SSBEffcal->DoubleElec_EffROOT(Lep1,Lep2,Elec_Supercluster_Eta->at(v_lepton_idx[0]) ,Elec_Supercluster_Eta->at(v_lepton_idx[1]) ,"central","central");// LepRecoSFSys is for Lepton //
         }
         else if ( TString(Decaymode).Contains( "muel" ) )
         {
            lep_sf_cent = SSBEffcal->MuonElec_EffROOT(TMuon,TElectron,Elec_Supercluster_Eta->at(v_electron_idx[0]) ,"central","central","central","central");//LepID, LepIso, Tracksys, ReconSys
         }
         else {lep_sf_cent = 0.0;cout << "Error In LeptonSFApply" << endl;}
      }
      for( int i =0; i < v_SystFullName.size(); ++i )
      {
         if ( !TString(FileName_).Contains( "Data") ){
            if ( v_SystFullName[i].Contains("LepIDUp"))
            {
               if ( TString(Decaymode).Contains( "dimuon" ) ) 
               {
                  //v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,"up","central"); // lep1 lep2 idsys isosys 
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,"up","central","central"); // lep1 lep2 idsys isosys tracksys 
               }
               else if ( TString(Decaymode).Contains( "dielec" ) )
               {
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleElec_EffROOT(Lep1,Lep2,Elec_Supercluster_Eta->at(v_lepton_idx[0]) ,Elec_Supercluster_Eta->at(v_lepton_idx[1]) ,"up","central");// LepRecoSFSys is for Lepton //
               }
               else if ( TString(Decaymode).Contains( "muel" ) )
               {
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->MuonElec_EffROOT(TMuon,TElectron,Elec_Supercluster_Eta->at(v_electron_idx[0]) ,"up","central","central","central");//LepID, LepIso, Tracksys, ReconSys
               }
               else {cout << "LeptonSFApply error for All-in-One Syst.!!!!" << endl;}
            }
            else if ( v_SystFullName[i].Contains("LepIDDown"))
            {
               if ( TString(Decaymode).Contains( "dimuon" ) ) 
               {
                  //v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,"down","central"); // lep1 lep2 idsys isosys 
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,"down","central","central"); // lep1 lep2 idsys isosys tracksys 
               }
               else if ( TString(Decaymode).Contains( "dielec" ) )
               {
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleElec_EffROOT(Lep1,Lep2,Elec_Supercluster_Eta->at(v_lepton_idx[0]) ,Elec_Supercluster_Eta->at(v_lepton_idx[1]) ,"down","central");// LepRecoSFSys is for Lepton //
               }
               else if ( TString(Decaymode).Contains( "muel" ) )
               {
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->MuonElec_EffROOT(TMuon,TElectron,Elec_Supercluster_Eta->at(v_electron_idx[0]) ,"down","central","central","central");//LepID, LepIso, Tracksys, ReconSys
               }
               else {cout << "LeptonSFApply error for All-in-One Syst.!!!!" << endl;}

            }
            else if ( v_SystFullName[i].Contains("LepIsoUp"))
            {
               if ( TString(Decaymode).Contains( "dimuon" ) ) 
               {
                  //v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,"central","up"); // lep1 lep2 idsys isosys 
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,"central","up","central"); // lep1 lep2 idsys isosys tracksys
               }
               else if ( TString(Decaymode).Contains( "dielec" ) )
               {
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleElec_EffROOT(Lep1,Lep2,Elec_Supercluster_Eta->at(v_lepton_idx[0]) ,Elec_Supercluster_Eta->at(v_lepton_idx[1]) ,"central","central");// LepRecoSFSys is for Lepton //
               }
               else if ( TString(Decaymode).Contains( "muel" ) )
               {
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->MuonElec_EffROOT(TMuon,TElectron,Elec_Supercluster_Eta->at(v_electron_idx[0]) ,"central","up","central","central");//LepID, LepIso, Tracksys, ReconSys
               }
               else {cout << "LeptonSFApply error for All-in-One Syst.!!!!" << endl;}
            }
            else if ( v_SystFullName[i].Contains("LepIsoDown"))
            {
               if ( TString(Decaymode).Contains( "dimuon" ) ) 
               {
                  //v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,"central","down"); // lep1 lep2 idsys isosys 
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleMuon_EffROOT(Lep1,Lep2,"central","down","central"); // lep1 lep2 idsys isosys tracksys
               }
               else if ( TString(Decaymode).Contains( "dielec" ) )
               {
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->DoubleElec_EffROOT(Lep1,Lep2,Elec_Supercluster_Eta->at(v_lepton_idx[0]) ,Elec_Supercluster_Eta->at(v_lepton_idx[1]) ,"central","central");// LepRecoSFSys is for Lepton //
               }
               else if ( TString(Decaymode).Contains( "muel" ) )
               {
                  v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->MuonElec_EffROOT(TMuon,TElectron,Elec_Supercluster_Eta->at(v_electron_idx[0]) ,"central","down","central","central");//LepID, LepIso, Tracksys, ReconSys
               }
               else {cout << "LeptonSFApply error for All-in-One Syst.!!!!" << endl;}

            }

            else {
               /// At here, you will use leptons sf as central value ///
               v_SystEvt[i] = v_SystEvt[i]*lep_sf_cent;
            }
         }
         else { // Data 
            v_SystEvt[i]=1;
         }
         m_Syst_EvtW[ v_SystFullName[i] ] = v_SystEvt[i]; 
//         cout << "v_SystFullName[" << i << "] : " << v_SystFullName[i] << " : " << v_SystEvt[i] << endl;
      }
   }
}
// Jet Selection and Jet cleaning 
void ssb_analysis::JetSelector()
{
 
   v_jet_idx.clear();
   v_jet_idx_cand.clear();
   v_jetup_idx_cand.clear();
   v_jetdn_idx_cand.clear();
   v_jet_TL_cand.clear();
   v_jetup_TL_cand.clear();
   v_jetdn_TL_cand.clear();

   v_jet_TL.clear();
   v_jetup_TL.clear();
   v_jetdn_TL.clear();

   v_ele_idx_cand.clear();
   v_mu_idx_cand.clear();
   v_jetup_idx.clear();
   v_jetdn_idx.clear();

   Int_t injet = Jet->GetEntriesFast();
   Int_t inele = Elec->GetEntriesFast();
   Int_t nmu   = Muon->GetEntriesFast();

   TJet = new TLorentzVector();

   for (int i =0; i < injet; i++)
   {
      TLorentzVector *Jetpre = (TLorentzVector*)Jet->At(i);

      /////////////////////
      /// Smearing Jet  ///
      /////////////////////
      Jetpre = JERSmearing(Jetpre,i);

      //////////////////////////////
      /// JES Systematic Method1 ///
      //////////////////////////////
      if      ( TString(JetEnSys).Contains("JetEnNorm"     ) ) { TJet=Jetpre; }
      else if ( TString(JetEnSys).Contains("JetEnShiftedUp") )
      { 
         TJet->SetPtEtaPhiE(Jetpre->Pt()*Jet_EnShiftedUp->at(i),Jetpre->Eta(),Jetpre->Phi(),Jetpre->Energy()*Jet_EnShiftedUp->at(i));
      }
      else if ( TString(JetEnSys).Contains("JetEnShiftedDown") )
      {
         TJet->SetPtEtaPhiE(Jetpre->Pt()*Jet_EnShiftedDown->at(i),Jetpre->Eta(),Jetpre->Phi(),Jetpre->Energy()*Jet_EnShiftedDown->at(i));
      }
      else { cout << "Check out Jet Systematic Config ... Default is Norminal Jet ..."<< endl; TJet=Jetpre; }

      

      if ( (*TJet).Pt() > jet_pt && fabs( (*TJet).Eta() ) < jet_eta && v_jet_Id->at(i) >= jet_id )
      {
         v_jet_idx_cand.push_back( i ); 
         v_jet_TL_cand.push_back( TJet ); 
         if ( (*TJet).Pt() < 30. ){cout << "Something wrong " << endl;}
      }

      ////////////////////////////////////////
      /// For JES Systematic with Method 2 ///
      ////////////////////////////////////////
      TLorentzVector *JetUp =new TLorentzVector();
      JetUp->SetPtEtaPhiE(Jetpre->Pt()*Jet_EnShiftedUp->at(i),Jetpre->Eta(),Jetpre->Phi(),Jetpre->Energy()*Jet_EnShiftedUp->at(i));
      TLorentzVector *JetDown = new TLorentzVector();
      JetDown->SetPtEtaPhiE(Jetpre->Pt()*Jet_EnShiftedDown->at(i),Jetpre->Eta(),Jetpre->Phi(),Jetpre->Energy()*Jet_EnShiftedDown->at(i));
      if ( (*JetUp).Pt() > jet_pt && fabs( (*JetUp).Eta() ) < jet_eta && v_jet_Id->at(i) >= jet_id )
      {
         v_jetup_idx_cand.push_back( i ); 
         v_jetup_TL_cand.push_back( JetUp ); 
         if ( (*JetUp).Pt() < 30. ){cout << "Something wrong " << endl;}
      }
      if ( (*JetDown).Pt() > jet_pt && fabs( (*JetDown).Eta() ) < jet_eta && v_jet_Id->at(i) >= jet_id )
      {
         v_jetdn_idx_cand.push_back( i ); 
         v_jetdn_TL_cand.push_back( JetDown ); 
         if ( (*JetDown).Pt() < 30. ){cout << "Something wrong " << endl;}
      }
   }

   ///////////////////////////////////////////
   /// Electron Selection For Jet Cleaning /// 
   ///////////////////////////////////////////
   for (int ine = 0; ine <inele; ine++ )
   {  
      Ele_jet_clean = (TLorentzVector*)Elec->At(ine); 

      if (
           Elec_SCB_Veto->at(ine) == true &&
           Ele_jet_clean->Pt() > el_pt_jetcl &&
           fabs( Ele_jet_clean->Eta() ) < el_eta_jetcl &&
           Elec_Conversion->at(ine) == true
         )
      {
//         v_ele_idx_cand.push_back( ine );
         if ( fabs( Elec_Supercluster_Eta->at(ine) ) <= 1.479 )
         {
            if (
                 v_el_iso_jcl->at( ine ) < 0.175 &&
                 fabs( Elec_Track_GsfdXY->at(ine) ) < 0.05 &&
                 fabs( Elec_Track_GsfdZ->at(ine) ) < 0.10 
               )
            { v_ele_idx_cand.push_back( ine ); }

         } 
         else 
         {
            if (
                 v_el_iso_jcl->at( ine ) < 0.159 &&           
                 fabs( Elec_Track_GsfdXY->at(ine) ) < 0.10 &&
                 fabs( Elec_Track_GsfdXY->at(ine) ) < 0.20 
               )
            { v_ele_idx_cand.push_back( ine ); }
        
         }
      }

   }
     
   // Muon selectron for JetCleaning
   for ( int inmu =0; inmu < nmu; ++inmu )
   {
      Mu_jet_clean = (TLorentzVector*)Muon->At(inmu);
      if ( 
           (*Mu_jet_clean).Pt() > mu_pt_jetcl         && //Pt cut 
           fabs((*Mu_jet_clean).Eta()) < mu_eta_jetcl && //Eta cut
           v_mu_iso_jcl->at(inmu) < mu_iso_jetcl      && //Isolation cut...
           v_mu_Id_jcl->at(inmu) == true                 // Id cut ...
         )
      { v_mu_idx_cand.push_back( inmu ); }
   }

   // Jet Cleaning 
   for ( int i =0; i < v_jet_idx_cand.size(); ++i )
   {
      TJet = new TLorentzVector();
      l_jet = 0;
      TJet = v_jet_TL_cand[i];
      // electron
      if ( v_ele_idx_cand.size() > 0 )
      {
         for ( int j = 0; j < v_ele_idx_cand.size(); ++j )
         {
            selected_Ele_jet_clean = (TLorentzVector*)Elec->At( v_ele_idx_cand.at(j) );
            if ( (*selected_Ele_jet_clean).DeltaR( (*TJet) ) < 0.4){ l_jet++; }
         }
      }

      // muon
      if ( v_mu_idx_cand.size() > 0 )
      {  
         for ( int k = 0; k < v_mu_idx_cand.size(); ++k )
         {
            selected_Mu_jet_clean = (TLorentzVector*)Muon->At( v_mu_idx_cand.at(k) );
            if ( (*selected_Mu_jet_clean).DeltaR( (*TJet) ) < 0.4){ l_jet++; }
         }
      }

      if ( l_jet == 0 ) 
      {
         v_jet_idx.push_back( v_jet_idx_cand.at(i) );
         v_jet_TL.push_back( v_jet_TL_cand.at(i) );
      }
   }
   // Jet Cleaning for All-in-One Syst. //
   for ( int i =0; i < v_jetup_idx_cand.size(); ++i ) // JES Up 
   {
      l_jet = 0;
      TJet =  v_jetup_TL_cand.at(i);
      // electron
      if ( v_ele_idx_cand.size() > 0 )
      {
         for ( int j = 0; j < v_ele_idx_cand.size(); ++j )
         {
            selected_Ele_jet_clean = (TLorentzVector*)Elec->At( v_ele_idx_cand.at(j) );
            if ( (*selected_Ele_jet_clean).DeltaR( (*TJet) ) < 0.4){ l_jet++; }
         }
      }
      // muon
      if ( v_mu_idx_cand.size() > 0 )
      {  
         for ( int k = 0; k < v_mu_idx_cand.size(); ++k )
         {
            selected_Mu_jet_clean = (TLorentzVector*)Muon->At( v_mu_idx_cand.at(k) );
            if ( (*selected_Mu_jet_clean).DeltaR( (*TJet) ) < 0.4){ l_jet++; }
         }
      }
      if ( l_jet == 0 ) 
      {
         v_jetup_idx.push_back( v_jetup_idx_cand.at(i) );
         v_jetup_TL.push_back( v_jetup_TL_cand.at(i) );
      }
   }

   for ( int i =0; i < v_jetdn_idx_cand.size(); ++i ) // Jet JES Down
   {
      l_jet = 0;
      //TJet = (TLorentzVector*)Jet->At( v_jetdn_idx_cand.at(i) );
      TJet = v_jetdn_TL_cand.at(i);
      // electron
      if ( v_ele_idx_cand.size() > 0 )
      {
         for ( int j = 0; j < v_ele_idx_cand.size(); ++j )
         {
            selected_Ele_jet_clean = (TLorentzVector*)Elec->At( v_ele_idx_cand.at(j) );
            if ( (*selected_Ele_jet_clean).DeltaR( (*TJet) ) < 0.4){ l_jet++; }
         }
      }

      // muon
      if ( v_mu_idx_cand.size() > 0 )
      {  
         for ( int k = 0; k < v_mu_idx_cand.size(); ++k )
         {
            selected_Mu_jet_clean = (TLorentzVector*)Muon->At( v_mu_idx_cand.at(k) );
            if ( (*selected_Mu_jet_clean).DeltaR( (*TJet) ) < 0.4){ l_jet++; }
         }
      }

      if ( l_jet == 0 ) 
      {
         v_jetdn_idx.push_back( v_jetdn_idx_cand.at(i) );
         v_jetdn_TL.push_back( v_jetdn_TL_cand.at(i) );
      }
   }
   if ( TString(FileName_).Contains( "Data" ) )
   {
      v_jetdn_idx.clear();
      v_jetup_idx.clear();
      v_jetdn_TL.clear();
      v_jetup_TL.clear();
      v_jetdn_idx = v_jet_idx;
      v_jetup_idx = v_jet_idx;
      v_jetdn_TL = v_jet_TL;
      v_jetup_TL = v_jet_TL;
   }
}
void ssb_analysis::JetDefiner()
{
   Jet1 = new TLorentzVector();
   Jet2 = new TLorentzVector();
   Jet1Up = new TLorentzVector();
   Jet2Up = new TLorentzVector();
   Jet1Dn = new TLorentzVector();
   Jet2Dn = new TLorentzVector();

   if ( v_jet_TL.size() == 1 ) { Jet1 = v_jet_TL[0]; }
   if ( v_jet_TL.size() > 1 ) { Jet1 = v_jet_TL[0]; Jet2 = v_jet_TL[1]; }
   if ( v_jetup_TL.size() == 1 ) { Jet1Up = v_jetup_TL[0]; }
   if ( v_jetup_TL.size() > 1 ) { Jet1Up = v_jetup_TL[0]; Jet2Up = v_jetup_TL[1]; }
   if ( v_jetdn_TL.size() == 1 ) { Jet1Dn = v_jetdn_TL[0]; }
   if ( v_jetdn_TL.size() > 1 ) { Jet1Dn = v_jetdn_TL[0]; Jet2Dn = v_jetdn_TL[1]; }

   if ( TString(JetEnSys).Contains("JetEnNorm"     ) )
   { 
      Jet1 = Jet1;
      Jet2 = Jet2;
      Jet1Up = Jet1Up;
      Jet2Up = Jet2Up;
      Jet1Dn = Jet1Dn;
      Jet2Dn = Jet2Dn;
   }
   else if ( TString(JetEnSys).Contains("JetEnShiftedUp") )
   {
      Jet1 = Jet1Up;
      Jet2 = Jet2Up;
   }
   else if ( TString(JetEnSys).Contains("JetEnShiftedDown") )
   {
      Jet1 = Jet1Dn;
      Jet2 = Jet2Dn;
   }
   else
   {
      cout << "Check out Jet Systematic Config ... Default is Norminal Jet ..."<< endl;
   }
}
void ssb_analysis::JetDefinerTL( std::vector<TLorentzVector*> v_jet, std::vector<TLorentzVector*> v_jetup, std::vector<TLorentzVector*> v_jetdn)
{
   Jet1 = new TLorentzVector();
   Jet2 = new TLorentzVector();
   Jet1Up = new TLorentzVector();
   Jet2Up = new TLorentzVector();
   Jet1Dn = new TLorentzVector();
   Jet2Dn = new TLorentzVector();

   if ( v_jet.size() == 1 ) { Jet1 = v_jet[0]; }
   if ( v_jet.size() > 1 ) { Jet1 = v_jet[0]; Jet2 = v_jet[1]; }
   if ( v_jetup.size() == 1 ) { Jet1Up = v_jetup[0]; }
   if ( v_jetup.size() > 1 ) { Jet1Up = v_jetup[0]; Jet2Up = v_jetup[1]; }
   if ( v_jetdn.size() == 1 ) { Jet1Dn = v_jetdn[0]; }
   if ( v_jetdn.size() > 1 ) { Jet1Dn = v_jetdn[0]; Jet2Dn = v_jetdn[1]; }

   if ( TString(JetEnSys).Contains("JetEnNorm"     ) )
   { 
/*      Jet1 = Jet1;
      Jet2 = Jet2;
      Jet1Up = Jet1Up;
      Jet2Up = Jet2Up;
      Jet1Dn = Jet1Dn;
      Jet2Dn = Jet2Dn;*/
   }
   else if ( TString(JetEnSys).Contains("JetEnShiftedUp") )
   {
      Jet1 = Jet1Up;
      Jet2 = Jet2Up;
   }
   else if ( TString(JetEnSys).Contains("JetEnShiftedDown") )
   {
      Jet1 = Jet1Dn;
      Jet2 = Jet2Dn;
   }
   else
   {
      cout << "Check out Jet Systematic Config ... Default is Norminal Jet ..."<< endl;
   }
}
void ssb_analysis::METDefiner()
{
 
   Met = new TLorentzVector();
   MetJESUp = new TLorentzVector();
   MetJESDn = new TLorentzVector();

   /// For Data, We don't need to apply JES systematic ///   
   if (TString(FileName_).Contains( "Data" ))
   {
      Met =  (TLorentzVector*)METMUEGCleanCor->At(0); 
      MetJESDn = Met; MetJESUp = Met;
   }
   else {
      Met =  (TLorentzVector*)METMUCleanCor->At(0);
      MetJESUp->SetPtEtaPhiE(METMUCleanCor_JetEnShiftedUp_PT->at(0) ,0, METMUCleanCor_JetEnShiftedUp_Phi->at(0),0);
      MetJESDn->SetPtEtaPhiE(METMUCleanCor_JetEnShiftedDown_PT->at(0) ,0, METMUCleanCor_JetEnShiftedDown_Phi->at(0),0);
      if ( TString(JetEnSys).Contains("JetEnNorm"     ) )
      {
         Met =  (TLorentzVector*)METMUCleanCor->At(0);
      }
      else if ( TString(JetEnSys).Contains("JetEnShiftedUp") )
      {
         Met = MetJESUp;
      }
      else if ( TString(JetEnSys).Contains("JetEnShiftedDown") )
      {
         Met = MetJESDn;
      }
      else 
      {
         Met = (TLorentzVector*)MET->At(0);
         cout << "Check out Jet Systematic Config ... Default is Norminal MET ..."<< endl;      
      }
   }
}
// DiMuon Mass Cut : step 1
bool ssb_analysis::DiLeptonMassCut()
{
   dimu_masscut = false;

   if ( ((*Lep1)+(*Lep2)).M() > 20 ){ dimu_masscut = true; }

   return dimu_masscut;
}

// ZVeto Cut : step 2
bool ssb_analysis::ZVetoCut()
{
   zvetocut = false;

   if ( TString(Decaymode).Contains( "dimuon" ) || TString(Decaymode).Contains( "dielec" ) ) 
   {
      if ( ((*Lep1)+(*Lep2)).M() <= 76 || ((*Lep1)+(*Lep2)).M() >= 106 ){ zvetocut = true; }
   }
   else if ( TString(Decaymode).Contains( "muel" ) ){ zvetocut = true; }
   else {cout << "ZVeto Error !!" << endl;}

   return zvetocut;
}

// Num.Jet Cut : Step 3
bool ssb_analysis::NumJetCut(std::vector<int> v_jets)
{
   numjetcut = false;
   //if ( v_jets.size() == 2 ){ numjetcut = true; }
   if ( v_jets.size() >= 2 ){ numjetcut = true; }
   return numjetcut;
}

//MET Cut : step 4
bool ssb_analysis::METCut(TLorentzVector* met)
{
   metcut = false;
   if ( TString(Decaymode).Contains( "dimuon" ) || TString(Decaymode).Contains( "dielec" ) )
   {
      if (met->Pt() > 40) { metcut =true; }
   }
   else if ( TString(Decaymode).Contains( "muel" ) ){ metcut =true; }
   else {cout << "METCut Error !!" << endl;}
   return metcut;
}
void ssb_analysis::BDsicApply()
{
   v_bjet_idx.clear();
   v_bjetup_idx.clear();
   v_bjetdn_idx.clear();
   v_bjet_TL.clear();
   v_bjetup_TL.clear();
   v_bjetdn_TL.clear();
   nbtagged =0;
   for (int ije = 0; ije < v_jet_idx.size(); ije++)
   {
       if (Jet_bDisc->at( v_jet_idx.at(ije) ) > bdisccut ) { nbtagged++;v_bjet_idx.push_back( v_jet_idx.at(ije) ); v_bjet_TL.push_back( v_jet_TL.at(ije) ); }
   }
   /// For JES Up 
   for (int ije = 0; ije < v_jetup_idx.size(); ije++)
   {
       if (Jet_bDisc->at( v_jetup_idx.at(ije) ) > bdisccut ) { v_bjetup_idx.push_back( v_jetup_idx.at(ije) ); v_bjetup_TL.push_back( v_jetup_TL.at(ije) ); }
   }
   /// For JES Down 
   for (int ije = 0; ije < v_jetdn_idx.size(); ije++)
   {
       if (Jet_bDisc->at( v_jetdn_idx.at(ije) ) > bdisccut ) { v_bjetdn_idx.push_back( v_jetdn_idx.at(ije) ); v_bjetdn_TL.push_back( v_jetdn_TL.at(ije) ); }
   }
}
//B Jet Cut : step 5
bool ssb_analysis::BJetCut( std::vector<int>v_bjets )
{
   bjetcut = false;
   if (v_bjets.size() >= 1) { bjetcut = true; }
   return bjetcut;
}
void ssb_analysis::BTaggigSFApply()
{
   evt_weight_beforeBtag_ = 1;
   evt_weight_beforeBtag_ = evt_weight_;
   if( TString(FileName_).Contains( "Data" ) ) {return;}
   std::vector<double> v_Jet_pT;
   std::vector<double> v_Jet_eta;
   std::vector<double> v_Jet_bDisc;
   std::vector<int> v_Jet_flav;
   v_Jet_pT.clear();
   v_Jet_eta.clear();
   v_Jet_bDisc.clear();
   v_Jet_flav.clear();
   TLorentzVector* jet_ = new TLorentzVector();
   for ( int ijet = 0; ijet < v_jet_idx.size(); ++ijet )
   {
      int idx_ = v_jet_idx[ijet];
//      TLorentzVector* jet_ = (TLorentzVector*)Jet->At(idx_);
//      jet_ = (TLorentzVector*)v_jet_TL[idx_];
      jet_ = (TLorentzVector*)v_jet_TL[ijet];
      v_Jet_pT.push_back(jet_->Pt()); 
      v_Jet_eta.push_back(jet_->Eta()); 
      v_Jet_bDisc.push_back(Jet_bDisc->at(idx_)); 
      v_Jet_flav.push_back(abs(Jet_PartonFlavour->at(idx_))); 
   } 

   if( !TString(FileName_).Contains( "Data" ) ) {evt_weight_ = SSBEffcal->Btagging_EvenWeight( v_Jet_pT, v_Jet_eta, v_Jet_bDisc, bdisccut,v_Jet_flav )*evt_weight_;}
   else { evt_weight_ = 1; }
   if ( isAllSyst == true )
   {
      double b_sf_cent = 1.0;
      b_sf_cent = SSBEffcal->Btagging_EvenWeightSys(v_Jet_pT,v_Jet_eta,v_Jet_bDisc,bdisccut,"Central");
      for( int i =0; i < v_SystFullName.size(); ++i )
      {
         if ( !TString(FileName_).Contains( "Data") ){
            if ( v_SystFullName[i].Contains("JetEnUp") || v_SystFullName[i].Contains("JetEnDown") ){continue;}
            else if ( v_SystFullName[i].Contains("BTagSF") || v_SystFullName[i].Contains("BTagEff") )
            {
               v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->Btagging_EvenWeightSys(v_Jet_pT, v_Jet_eta, v_Jet_bDisc, bdisccut,v_SystFullName[i]);
            }
            else {
               v_SystEvt[i] = v_SystEvt[i]*b_sf_cent;
            }
         }
         else { // Data 
            v_SystEvt[i]=1;
         }
         m_Syst_EvtW[ v_SystFullName[i] ] = v_SystEvt[i]; 
      }
   }
}
void ssb_analysis::BTaggigSFApplyJES(TString sys_)
{
   if( TString(FileName_).Contains( "Data" ) ) {return;}
   std::vector<double> v_JetUp_pT;
   std::vector<double> v_JetUp_eta;
   std::vector<double> v_JetUp_bDisc;
   v_JetUp_pT.clear();
   v_JetUp_eta.clear();
   v_JetUp_bDisc.clear();

   std::vector<double> v_JetDn_pT;
   std::vector<double> v_JetDn_eta;
   std::vector<double> v_JetDn_bDisc;
   v_JetDn_pT.clear();
   v_JetDn_eta.clear();
   v_JetDn_bDisc.clear();

   for ( int ijetup = 0; ijetup < v_jetup_idx.size(); ++ijetup )
   {
      int idx_ = v_jetup_idx[ijetup];
      TLorentzVector* jetup_ =v_jetup_TL[ijetup] ;
      v_JetUp_pT.push_back(jetup_->Pt()); 
      v_JetUp_eta.push_back(jetup_->Eta()); 
      v_JetUp_bDisc.push_back(Jet_bDisc->at(idx_)); 
   }
   for ( int ijetdn = 0; ijetdn < v_jetdn_idx.size(); ++ijetdn )
   {
      int idx_ = v_jetdn_idx[ijetdn];
      TLorentzVector* jetdn_ =v_jetdn_TL[ijetdn] ;
      v_JetDn_pT.push_back(jetdn_->Pt()); 
      v_JetDn_eta.push_back(jetdn_->Eta()); 
      v_JetDn_bDisc.push_back(Jet_bDisc->at(idx_)); 
   } 
   if ( isAllSyst == true )
   {
      for( int i =0; i < v_SystFullName.size(); ++i )
      {
         if ( !TString(FileName_).Contains( "Data") ){
            if ( v_SystFullName[i].Contains("JetEnUp") && sys_ ){
               v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->Btagging_EvenWeightSys(v_JetUp_pT,v_JetUp_eta,v_JetUp_bDisc,bdisccut,"Central"); 
               m_Syst_EvtW[ v_SystFullName[i] ] = v_SystEvt[i]; 
            }
            else if ( v_SystFullName[i].Contains("JetEnDown") && sys_ ){
               v_SystEvt[i] = v_SystEvt[i]*SSBEffcal->Btagging_EvenWeightSys(v_JetDn_pT,v_JetDn_eta,v_JetDn_bDisc,bdisccut,"Central"); 
               m_Syst_EvtW[ v_SystFullName[i] ] = v_SystEvt[i]; 
            }
            else { continue; }
         }
         else { // Data 
            v_SystEvt[i]=1;
         }
      }
   }
   else {} 
}
void ssb_analysis::BJetDefiner(std::vector<int> v_jets, std::vector<int> v_bjets)
{
   //////////////////
   // Define BJets //
   //////////////////
   bcandi = -9;
   if( v_bjet_TL.size() > 1)
   {
      bJet1 = v_bjet_TL.at(0); 
      bJet2 = v_bjet_TL.at(1);
   }

   else if (  v_bjet_idx.size() > 0 )
   {
      bJet1 = v_bjet_TL.at(0);
      for ( int k = 0; k < v_jet_idx.size(); ++k)
      {
         if ( v_jet_idx.at(k) != v_bjet_idx.at(0) ){ bJet2 =v_jet_TL[k]; break; }
      }
   }
 
   // Bjet TLorentzVector
   else if ( v_bjet_idx.size() == 0 )
   {
      if ( v_jet_idx.size() > 1 )
      {
         bJet1 = v_jet_TL[0];
         bJet2 = v_jet_TL[1];
      }
   }
   else { cout << "you should check out nbjet cut!! " << endl; }
}
void ssb_analysis::BJetDefiner()
{
   //////////////////
   // Define BJets //
   //////////////////
   /// For Central Jet ///
   if( v_bjet_TL.size() > 1)
   {
      bJet1 = v_bjet_TL.at(0); 
      bJet2 = v_bjet_TL.at(1);
   }
   else if ( v_bjet_TL.size() > 0 )
   {  
      bJet1 = v_bjet_TL.at(0);
      for ( int i = 0; i < v_jet_idx.size(); ++i )
      {
         if ( v_jet_idx[i] != v_bjet_idx[0] ){ bJet2 = v_jet_TL[i]; break; }
      }
   }
   else if ( v_bjet_TL.size() == 0 )
   {
      bJet1 = Jet1;
      bJet2 = Jet2;
   }
   else { cout << "you should check out nbjet cut!! Error about v_bjet_TL !! " << endl; }
   /// For JES ShiftedUp Jet ///
   if( v_bjetup_TL.size() > 1)
   {
      bJet1Up = v_bjetup_TL.at(0); 
      bJet2Up = v_bjetup_TL.at(1);
   }
   else if ( v_bjetup_TL.size() > 0 )
   {  
      bJet1Up = v_bjetup_TL.at(0);
      for ( int i = 0; i < v_jet_idx.size(); ++i )
      {
         if ( v_jetup_idx[i] != v_bjetup_idx[0] ){ bJet2Up = v_jetup_TL[i]; break; }
      }
   }
   else if ( v_bjetup_TL.size() == 0 )
   {
      bJet1Up = Jet1Up;
      bJet2Up = Jet2Up;
   }
   else { cout << "you should check out nbjet cut!! Error about v_bjetup_TL !! " << endl; }
   /// For JES ShiftedDown Jet ///
   if( v_bjetdn_TL.size() > 1)
   {
      bJet1Dn = v_bjetdn_TL.at(0); 
      bJet2Dn = v_bjetdn_TL.at(1);
   }
   else if ( v_bjetdn_TL.size() > 0 )
   {  
      bJet1Dn = v_bjetdn_TL.at(0);
      for ( int i = 0; i < v_jet_idx.size(); ++i )
      {
         if ( v_jet_idx[i] != v_bjet_idx[0] ){ bJet2Dn = v_jetdn_TL[i]; break; }
      }
   }
   else if ( v_bjetdn_TL.size() == 0 )
   {
      bJet1Dn = Jet1Dn;
      bJet2Dn = Jet2Dn;
   }

   else { cout << "you should check out nbjet cut!! Error about v_bjetdn_TL !! " << endl; }
}
/////////////////////
//Double batagging //
/////////////////////
bool ssb_analysis::DoubleBtag( std::vector<int>v_bjets )
{
   doublebtag = false;
   if ( v_bjets.size() >=2 ) {doublebtag = true;}
   return doublebtag;
}
//Kinematic Solver to reconstruct top
void ssb_analysis::KinSol()
{

   xconstraint = 0; yconstraint  = 0;
 
   ksolweight_ = -1.0;
   /// initialize array
   //cout << "At KinSol -- " << "Lep->Px() : " << Lep->Px() << " AnLep->Px() : " << AnLep->Px() << " bJet1->Px() : " << bJet1->Px() << " bJet2->Px() : " << bJet2->Px() << " Met->Px() : " << Met->Px()<< endl;
   xconstraint = Lep->Px() + AnLep->Px() + bJet1->Px() + bJet2->Px() + Met->Px();
   yconstraint = Lep->Py() + AnLep->Py() + bJet1->Py() + bJet2->Py() + Met->Py();

   ssbflsolver->SetConstraints( xconstraint , yconstraint );  

   TtFullLepKinSolver::NeutrinoSolution nuSol1  = ssbflsolver->getNuSolution( (*Lep),   (*AnLep), (*bJet1), (*bJet2) );
   TtFullLepKinSolver::NeutrinoSolution nuSol2  = ssbflsolver->getNuSolution( (*Lep),   (*AnLep), (*bJet2), (*bJet1) );


   if ( nuSol1.weight >nuSol2.weight )
   {
      ksolweight_ = nuSol1.weight;
      (*AnNu)     = nuSol1.neutrino;
      (*Nu)       = nuSol1.neutrinoBar;
      (*bJet)     = (*bJet2);
      (*AnbJet)   = (*bJet1);
   }
   else 
   {  
      ksolweight_ = nuSol2.weight;
      (*AnNu)     = nuSol2.neutrino; 
      (*Nu)       = nuSol2.neutrinoBar;
      (*bJet)     = (*bJet1);
      (*AnbJet)   = (*bJet2);
   }
   //(*Top)      = (*Lep)   + (*AnbJet) + (*AnNu);
   //(*AnTop)    = (*AnLep) + (*bJet)   + (*Nu);
   (*AnTop)     = (*Lep)   + (*AnbJet) + (*AnNu);
   (*Top)       = (*AnLep) + (*bJet)   + (*Nu);
   (*W1)        = (*Lep)   + (*AnNu);
   (*W2)        = (*AnLep) + (*Nu);

}
void ssb_analysis::KinSolSys( TLorentzVector* lep, TLorentzVector* alep, TLorentzVector* bj1, TLorentzVector* bj2, TLorentzVector* met )
{

   xconstraint = 0; yconstraint  = 0;
 
   ksolweight_ = -1.0;
   /// initialize array

   xconstraint = lep->Px() + alep->Px() + bj1->Px() + bJet2->Px() + met->Px();
   yconstraint = lep->Py() + alep->Py() + bj1->Py() + bJet2->Py() + met->Py();

   ssbflsolver->SetConstraints( xconstraint , yconstraint );  

   TtFullLepKinSolver::NeutrinoSolution nuSol1  = ssbflsolver->getNuSolution( (*lep),   (*alep), (*bj1), (*bj2) );
   TtFullLepKinSolver::NeutrinoSolution nuSol2  = ssbflsolver->getNuSolution( (*lep),   (*alep), (*bj2), (*bj1) );


   if ( nuSol1.weight >nuSol2.weight )
   {
      ksolweight_ = nuSol1.weight;
      (*AnNu)     = nuSol1.neutrino;
      (*Nu)       = nuSol1.neutrinoBar;
      (*bJet)     = (*bj2);
      (*AnbJet)   = (*bj1);
   }
   else 
   {  
      ksolweight_ = nuSol2.weight;
      (*AnNu)     = nuSol2.neutrino; 
      (*Nu)       = nuSol2.neutrinoBar;
      (*bJet)     = (*bj1);
      (*AnbJet)   = (*bj2);
   }
   (*AnTop)     = (*lep)  + (*AnbJet) + (*AnNu);
   (*Top)       = (*alep) + (*bJet)   + (*Nu);
   (*W1)        = (*lep)  + (*AnNu);
   (*W2)        = (*alep) + (*Nu);

}
void ssb_analysis::Bjorken( TLorentzVector *lep, TLorentzVector *alep,TLorentzVector *bj,TLorentzVector *abj,TLorentzVector *nu,TLorentzVector *anu )
{
   x1_bj = -999;
   x2_bj = -999;
   x3_bj = -999;
   x1_bj = ( lep->Pz() + alep->Pz() + bj->Pz() + abj->Pz() + nu->Pz() + anu->Pz() 
           + lep->E()  + alep->E()  + bj->E()  + abj->E()  + nu->E()  + anu->E()  )/(2*4*1000);
   x2_bj = ( lep->E()  + alep->E()  + bj->E()  + abj->E()  + nu->E()  + anu->E() 
           - lep->Pz() - alep->Pz() - bj->Pz() - abj->Pz() - nu->Pz() - anu->Pz() )/(2*4*1000);
   x3_bj = ( x1_bj + x2_bj )*(4*1000);
   
   
}

void ssb_analysis::LepAnLepMisCharge()
{
   laecO3 = 100; laecOb = 100; laecO5 = 100;
   if ( TString(Decaymode).Contains( "dimuon" ) )
   {
      if ( Lep->Pt() < 300 )
      {
         if ( AnLep->Pt() < 300 ){ laecO3 = 1.5; laecOb = 1.5; laecO5 = 1.5; }
         else                    { laecO3 = 3.5; laecOb = 3.5; laecO5 = 3.5; }
      }
      else 
      {
         if ( AnLep->Pt() < 300 ){ laecO3 = 5.5; laecOb = 5.5; laecO5 = 5.5; }
         else                    { laecO3 = 7.5; laecOb = 7.5; laecO5 = 7.5; }
      }
   }
   else if ( TString(Decaymode).Contains( "dielec" ) )
   {
      if( Lep->Pt() >20 && Lep->Pt() < 50 )
      {
         if( fabs( Lep->Eta() ) < 1.479  ) 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  ){
               if( fabs( AnLep->Eta() ) < 1.479 ) { laecO3 = 1.5; laecOb = 1.5; laecO5 = 1.5; }
               else                               { laecO3 = 3.5; laecOb = 3.5; laecO5 = 3.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { laecO3 = 5.5; laecOb = 5.5; laecO5 = 5.5; }
               else                                { laecO3 = 7.5; laecOb = 7.5; laecO5 = 7.5; }
            } //
         } // Lep eta < 1.479
         else // Lep 1.479  < eta < 2.5 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  )
            {
               if( fabs( AnLep->Eta() ) < 1.479  ) { laecO3 = 11.5; laecOb = 11.5; laecO5 = 11.5; }
               else                                { laecO3 = 13.5; laecOb = 13.5; laecO5 = 13.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { laecO3 = 15.5; laecOb = 15.5; laecO5 = 15.5; }
               else                                { laecO3 = 17.5; laecOb = 17.5; laecO5 = 17.5; }

            } //
         }
      }
      else // Lep >= 50GeV 
      {
         if( fabs( Lep->Eta() ) < 1.479  ) 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  ){
               if( fabs( AnLep->Eta() ) < 1.479 ) { laecO3 = 21.5; laecOb = 21.5; laecO5 = 21.5; }
               else                               { laecO3 = 23.5; laecOb = 23.5; laecO5 = 23.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { laecO3 = 25.5; laecOb = 25.5; laecO5 = 25.5; }
               else                                { laecO3 = 27.5; laecOb = 27.5; laecO5 = 27.5; }
            } //
         } // Lep eta < 1.479
         else // Lep 1.479  < eta < 2.5 
         {
            if( AnLep->Pt()  > 20 &&  AnLep->Pt()  < 50  )
            {
               if( fabs( AnLep->Eta() ) < 1.479  ) { laecO3 = 31.5; laecOb = 31.5; laecO5 = 31.5; }
               else                                { laecO3 = 33.5; laecOb = 33.5; laecO5 = 33.5; }
            }
            else { // AnLep pT >= 50 GeV
               if( fabs( AnLep->Eta() ) < 1.479  ) { laecO3 = 35.5; laecOb = 35.5; laecO5 = 35.5; }
               else                                { laecO3 = 37.5; laecOb = 37.5; laecO5 = 37.5; }

            } //
         }
      }
   }
   else if ( TString(Decaymode).Contains( "muel" ) )
   {
      if (TMuon->Pt() < 300 )
      {
         if( TElectron->Pt() > 20 && TElectron->Pt() < 50 )
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { laecO3 = 1.5; laecOb = 1.5; laecO5 = 1.5; } // Barrel
            else                                    { laecO3 = 3.5; laecOb = 3.5; laecO5 = 3.5; } // EndCap
         }
         else // Electron Pt > 50 GeV 
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { laecO3 = 5.5; laecOb = 5.5; laecO5 = 5.5; } // Barrel
            else                                    { laecO3 = 7.5; laecOb = 7.5; laecO5 = 7.5; } // Endcap
         }
      } 
      else 
      {
         if( TElectron->Pt() > 20 && TElectron->Pt() < 50 )
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { laecO3 = 11.5; laecOb = 11.5; laecO5 = 11.5; } // Barrel
            else                                    { laecO3 = 13.5; laecOb = 13.5; laecO5 = 13.5; } // EndCap
         }
         else // Electron Pt > 50 GeV 
         {
            if ( fabs( TElectron->Eta() ) < 1.479 ) { laecO3 = 15.5; laecOb = 15.5; laecO5 = 15.5; } // Barrel
            else                                    { laecO3 = 17.5; laecOb = 17.5; laecO5 = 17.5; } // Endcap
         }

      }

   }

}
double ssb_analysis::JetPtResol(TLorentzVector *jet)
{
   double jpr_;
   double jpt_;
   double jeta_;
   jpt_ = jet->Pt();
   jeta_ = jet->Eta();
   jpr_ = ssbcpviol->getJPR(jeta_,jpt_);
   return jpr_;
}

void ssb_analysis::TopPtReweightApply()
{
   if ( TString(FileName_).Contains( "TTJets" ) )
   {
//      TLorentzVector* genTop = v_GenTop[0];
//      TLorentzVector* genAnTop = v_GenAnTop[0];
      TLorentzVector* genTop   = (TLorentzVector*)GenTop->At(0);
      TLorentzVector* genAnTop = (TLorentzVector*)GenAnTop->At(0);
      double topptreweight;
      double a = 0.0615;
      double b = -0.0005;
      double toppt = genTop->Pt();
      double antoppt = genAnTop->Pt();
      double topsf = exp(a+b*toppt);
      double antopsf = exp(a+b*antoppt);
      topptreweight = sqrt(topsf*antopsf);
      evt_weight_ = evt_weight_*topptreweight;
   }
   else {}
   
   
}
TLorentzVector* ssb_analysis::JERSmearing(TLorentzVector *jet, int idx_)
{
   TLorentzVector* jetsmerd = new TLorentzVector();
   if ( dojer&!TString(FileName_).Contains( "Data" ) )
   {
//      cout << "Do the JER Smearing !!!" << endl;
      jetsmerd->SetPtEtaPhiE(jet->Pt()*Jet_EnergyResolution_SF->at(idx_),
                             jet->Eta(),jet->Phi(),
                            jet->Energy()*Jet_EnergyResolution_SF->at(idx_));
   }
   else {
      jetsmerd = jet;
   }
   //cout << "before smeared jet pt : " << jet->Pt()  << " Eta : " << jet->Eta() << " phi " << jet->Phi() << " Energy : " << jet->Energy() << endl;
   //cout << "after smeared jetsmerd pt : " << jetsmerd->Pt()  << " Eta : " << jetsmerd->Eta() << " phi " << jetsmerd->Phi() << " Energy : " << jetsmerd->Energy() << endl;
   return  jetsmerd;
}
void ssb_analysis::MakeVecforTL()
{
}
void ssb_analysis::ClearVectors()
{
   v_lepton_Id->clear();
   v_lepton_iso->clear();
   v_jet_Id->clear();
   v_mlep_iso->clear();
   v_elec_iso->clear();
   v_muon_Id->clear();  // for MuEl channel
   v_electron_Id->clear(); // for ElEl and MuEl channels

   /// User variables for JetCleanning cut
   v_mu_iso_jcl->clear();
   v_mu_Id_jcl->clear();
   v_el_iso_jcl->clear();
   v_el_Id_jcl->clear();
}
void ssb_analysis::ReadDuplList()
{
   map_Dilepevt.clear();
   /// Duplicate Event ///
   std::map<int,std::vector<int> > m_test;
   TString evtfileName = FileName_;
   char bar = '_';
   evtfileName.ReplaceAll("Single","Double");
   evtfileName.Remove(evtfileName.Last(bar));
   if (FileName_.Contains("Run2016B")) 
   {
      FILE *XSectionInfo;
      char LumiNum[1000];
      char EvtNum[1000];
      char RunNum[1000];
      char col1[1000];
      int lumi_number;
      int evt_number;
      int run_number;
      std::vector<int> v_EvtNums;
      v_EvtNums.clear();
      typedef std::map<int,std::vector<int> > map_lumi_evt;
      map_lumi_evt  m_lumi_evt;
      std::map<int,map_lumi_evt> m_run_lumievt;
      //map<TString,int> map_Dilepevt;
      for (int i = 1; i < 1+1; ++i)
      {
         ostringstream ostr;
         ostr.str("");
         ostr << i;
         string ostr_ = "_" + ostr.str() + ".txt";
         string evt_file = evtfileName.Data()  + ostr_;
         string slash = "/";
         string tttt= evtfileName.Data();
         string evt_filePath = "../EventList/" + tttt + slash + evt_file;
         cout << "evt_filePath : " << evt_filePath << endl;
         XSectionInfo = fopen(evt_filePath.c_str(),"r");                                                                               
         int evtnumber = 0;
         int runnumber = -1;
         int luminumber = -1;
         int idx_evt = 0;
         int idx_lumi = 0;
         int idx_run = 0;
         map<int, map<int, map<int, int> > > dilepevt;
         if (XSectionInfo !=NULL){
            int lines_ = 0;
            v_EvtNums.clear();
            m_lumi_evt.clear();
            while (fscanf(XSectionInfo, "%s %s %d %s %d %s %d \n",RunNum,col1,&run_number, LumiNum,&lumi_number, EvtNum,& evt_number ) != EOF)
            {
               lines_++;
               if ( evtnumber == 0 || runnumber == -1 || luminumber ==-1) { cout << "First Line  run_number : " << run_number << " lumi_number : " << lumi_number << " evt_number : " << evt_number << endl; }// First line 
               else {

               }

/*               else {
                  if (runnumber != run_number )
                  {  
                     if ( m_run_lumievt.count(runnumber) == 0 ){  m_run_lumievt[ runnumber ] = m_lumi_evt; }
                     else {
                        if ( m_run_lumievt[runnumber].count(luminumber) == 0 ){ // there is no duplicate lumi number same as this lumi number
                           cout << "runnumber : " << runnumber << endl;
                           m_run_lumievt[runnumber][luminumber] = v_EvtNums;
                           cout << "  m_run_lumievt[runnumber][luminumber] size : " <<  m_run_lumievt[runnumber][luminumber].size() << " v_EvtNums.size () " << v_EvtNums.size() << endl;
                        }
                        else {
                           for (int isubv = 0;isubv < v_EvtNums.size(); ++isubv )
                           {
                              m_run_lumievt[runnumber][luminumber].push_back( v_EvtNums[isubv] ); 
                           }
                        }
                     }
                     cout << "runnumber : " << runnumber << " luminumber : " << luminumber << " m_run_lumievt[runnumber][luminumber] size: " << m_run_lumievt[runnumber][luminumber].size() << endl;
                     m_lumi_evt.clear();
                     v_EvtNums.clear();
                  }
                  else {
                     if ( lumi_number != luminumber ) 
                     {
                        if (m_lumi_evt.count(luminumber) == 0 ){m_lumi_evt[ luminumber ] = v_EvtNums;}
                        else {
                           for (int isubv =0; isubv < v_EvtNums.size(); ++isubv){
                              m_lumi_evt[ luminumber ].push_back( v_EvtNums[isubv] );
                           }
                        }
                        v_EvtNums.clear();
                     }
                  }/// Run Num Check 
               }
               evtnumber = evt_number;
               runnumber = run_number;
               luminumber = lumi_number;
               v_EvtNums.push_back(evtnumber);
            */      
            }// End of ReadFile
            cout << " Last : runnumber " << runnumber << "Last : luminumber" << luminumber << endl;
            // For Last Section //
/*            if (m_lumi_evt.count(luminumber) == 0 ){m_lumi_evt[ luminumber ] = v_EvtNums;}
            else {
               for (int isubv =0; isubv < v_EvtNums.size(); ++isubv){
                  m_lumi_evt[ luminumber ].push_back( v_EvtNums[isubv] );
               }
            }
            if ( m_run_lumievt.count(runnumber) == 0 ){  m_run_lumievt[ runnumber ] = m_lumi_evt; }
            else { 
               if ( m_run_lumievt[runnumber].count(luminumber) == 0 ){ // there is no duplicate lumi number same as this lumi number
                  m_run_lumievt[runnumber][luminumber] = v_EvtNums;
               }
               else {
                  for (int isubv = 0;isubv < v_EvtNums.size(); ++isubv )
                  {
                     m_run_lumievt[runnumber][luminumber].push_back( v_EvtNums[isubv] ); 
                  }
               }
            }*/
            fclose(XSectionInfo);
            cout << "lines_ "<< lines_ << endl;
         }
         else {
            delete XSectionInfo;
         }
      }// End of loop of fileList //
      // Check out Dilepton List 
      int ilines = 0;
      int txtilines = 0;
      for (std::map<int,map_lumi_evt>::iterator iter = m_run_lumievt.begin(); iter != m_run_lumievt.end(); ++iter )
      {  //cout << "iter first : " << iter->first << endl;
         for (map_lumi_evt::iterator itsub = iter->second.begin(); itsub != iter->second.end(); ++itsub)
         {  
            //cout << "itsub first : " << itsub->first << endl;
            txtilines +=itsub->second.size();
            for (std::vector<int>::iterator ivec = itsub->second.begin(); ivec != itsub->second.end(); ++ivec)
            {
              ilines++; 
            }
         }
      }
      cout << "ilines ? " << ilines  << endl;
      cout << "txtilines ? " << txtilines  << endl;
   } // FileContains 
   cout << "End of Read Dilep Evt" << endl;
}
bool ssb_analysis::RMDuplEvt()
{
   /// Duplicate Event ///
//   cout << " map_Dilepevt : " << map_Dilepevt.size() <<endl;
//   cout << "Remove Duplicate Event " << endl;
   return true; 
}
