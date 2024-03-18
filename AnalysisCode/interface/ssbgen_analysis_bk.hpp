#ifndef _ssbgen_analysis_

#define _ssbgen_analysis_
  
#include <set>
#include <string>
#include <fstream>
#include <cassert>
#include <map>
#include <vector>
#include <string>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TGraphErrors.h>
#include "TLorentzVector.h"
#include "TEnv.h"
  
#include "./../analysis/SSBGenTree.h"
#include "./../interface/ssb_cpviol.hpp"

// Textreader
#include "./../TextReader/TextReader.hpp"

using namespace std;


class ssbgen_analysis : public SSBGenTree 
{
   public:
      //declare functions
      ssbgen_analysis(TTree *tree=0);
      virtual ~ssbgen_analysis();

      //basic frame
      virtual void Loop();
      void GetNtupleTotalEvent( unsigned int totevent );
      void Start();
      void End();

      //user define functions
      void SetOutputFileName( char *outname );
      void DeclareHistos();
      int MomFinder( int idx );
      // Read config files function
      void GetVariables();
      void TLVInitial();

      // Cut step Function 



   private:
      //put variables that you want
      char *outfile;
      TFile *fout;
      unsigned int NtupletotalEvent;
      unsigned int totalEvent;


      //TextReader from Jaehoon.
      TextReader *SSBConfReader;

      //CPViolation
      SSBCPViol *ssbcpviol; 

      //Map for Parton Finder
      std::map<int,int> idx_status;
      std::map<int,int> idx_pdgId;
      std::map<int,int> idx_mother1;
      std::map<int,int> idx_mother2;
      std::map<int,int> idx_doughter1;
      std::map<int,int> idx_doughter2;

      int num_lep;
      int num_mu;

      int num_dilep;
      int num_semilep;
      int num_had;
      int num_non;
      int num_dimu;
      int num_diel;
      int num_muel;
      int num_dimu_v1;
      int num_diel_v1;
      int num_muel_v1;
      int num_mutau_v1;
      int num_eltau_v1;
      int num_ditau_v1;

   public:

      //declare histograms
      TH1F *Channel_Had;
      TH1F *Channel_Semi;
      TH1F *Channel_Dilep;
      TH1F *Channel_non;

};
#endif

#ifdef ssbgen_analysis_cxx

ssbgen_analysis::ssbgen_analysis(TTree *tree)
{
   if (tree == 0)
   {
      printf("ERROR: Can't find any input tree.\n");
   }
   Init(tree);
   

   GetVariables();

   // Text Reader from Jaehoon.
   SSBConfReader = new TextReader();
   SSBConfReader->ReadFile("./configs/analysis_config.config");
   SSBConfReader->ReadVariables();
   // cpviolation 
   ssbcpviol   = new SSBCPViol();

}

ssbgen_analysis::~ssbgen_analysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

#endif
