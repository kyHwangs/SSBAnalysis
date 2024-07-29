////////////////////////
///                  ///  
///   Common Tools   ///
///                  ///
////////////////////////

////////////////////////////////////////////////
// Collection of tools for multiple purposes  //
////////////////////////////////////////////////

#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "THStack.h"
#include "TSystem.h"
#include "TEnv.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include <sstream>
#include "TStyle.h"
#include "TLorentzVector.h"

class read_my_params
{
   public:
      read_my_params (TEnv *env, std::vector<std::string> &fnames, EEnvLevel level=kEnvUser)
      : _env (env), _fnames (fnames), _level(level), _found_one (false){ }
      void SetupFile(std::string &dirname)
      {
         for(Int_t numfiles=0;numfiles<_fnames.size();numfiles++)
         {
            const std::string _fname=_fnames[numfiles];
            std::string fname (dirname + "/" + _fname);
            TString temp (fname.c_str());
            TString final_filename =temp; // (d0root_find_file (temp));

            if (final_filename != "") 
            {
               _env->ReadFile(final_filename.Data(), _level);
               _found_one = true;
            }
         }
      }
      bool GotOne (void) const {return _found_one;}

   private:
      TEnv *_env;
      std::vector<std::string> _fnames;
      EEnvLevel _level;
      bool _found_one;
};

void FillHisto(TH1 *hist, Double_t val, Double_t weight = 1.0);
void FillHisto(TH1D* hist, Double_t val, Double_t weight);
void FillHisto(TH2 *hist, Double_t valx, Double_t valy, Double_t weight = 1.0);
void FillHisto(TH3 *hist, Double_t valx, Double_t valy, Double_t valz, Double_t weight = 1.0);

