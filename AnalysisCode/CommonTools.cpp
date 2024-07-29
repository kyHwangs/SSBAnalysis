#include "CommonTools.hpp"
#include <stdlib.h>

void FillHisto(TH1 *hist, Double_t val, Double_t weight)
{
   Int_t nbins = hist->GetNbinsX();
   Double_t minval = hist->GetXaxis()->GetBinLowEdge(1);  
   Double_t maxval = hist->GetXaxis()->GetBinLowEdge(nbins) + hist->GetXaxis()->GetBinWidth(nbins);

   Double_t mincen = hist->GetXaxis()->GetBinCenter(1); 
   Double_t maxcen = hist->GetXaxis()->GetBinCenter(nbins); 

   // std::cout << hist->GetName() << " " << val << " " << weight << " " << minval << " " << maxval << std::endl;

   if (val < minval) {
      hist->Fill(mincen, weight);
      // std::cout << "Underflow: " << hist->GetEntries() << " " << hist->GetMean() << std::endl;
   } else if (val > maxval) {
      hist->Fill(maxcen, weight);
      // std::cout << "Overflow: " << hist->GetEntries() << " " << hist->GetMean() << std::endl;
   } else {
      hist->Fill(val, weight);
      // std::cout << "NORMAL: " << hist->GetEntries() << " " << hist->GetMean() << std::endl;
   }
}

void FillHisto(TH1D* hist, Double_t val, Double_t weight)
{
   Int_t nbins = hist->GetNbinsX();
   Double_t minval = hist->GetXaxis()->GetBinLowEdge(1);  
   Double_t maxval = hist->GetXaxis()->GetBinLowEdge(nbins) + hist->GetXaxis()->GetBinWidth(nbins);

   Double_t mincen = hist->GetXaxis()->GetBinCenter(1); 
   Double_t maxcen = hist->GetXaxis()->GetBinCenter(nbins); 

   // std::cout << hist->GetName() << " " << val << " " << weight << " " << minval << " " << maxval << std::endl;

   if (val < minval) {
      hist->Fill(mincen, weight);
      // std::cout << "Underflow: " << hist->GetEntries() << " " << hist->GetMean() << std::endl;
   } else if (val > maxval) {
      hist->Fill(maxcen, weight);
      // std::cout << "Overflow: " << hist->GetEntries() << " " << hist->GetMean() << std::endl;
   } else {
      hist->Fill(val, weight);
      // std::cout << "NORMAL: " << hist->GetEntries() << " " << hist->GetMean() << std::endl;
   }
}

void FillHisto(TH2 *hist, Double_t valx, Double_t valy, Double_t weight)
{
   Int_t nbinsx=hist->GetNbinsX();  
   Double_t minvalx=hist->GetXaxis()->GetBinCenter(1);  
   Double_t maxvalx=hist->GetXaxis()->GetBinCenter(nbinsx);        

   Int_t nbinsy=hist->GetNbinsY();
   Double_t minvaly=hist->GetYaxis()->GetBinCenter(1);  
   Double_t maxvaly=hist->GetYaxis()->GetBinCenter(nbinsy);  

   Double_t newvalx=valx;  
   Double_t newvaly=valy;

   if(valx< minvalx)    
      newvalx=minvalx;  
   else if(valx>maxvalx)
      newvalx=maxvalx;

   if(valy< minvaly)
      newvaly=minvaly;
   else if(valy>maxvaly)
      newvaly=maxvaly;

   hist->Fill(newvalx, newvaly, weight);
}

void FillHisto(TH3 *hist, Double_t valx, Double_t valy, Double_t valz, Double_t weight)
{
   Int_t nbinsx=hist->GetNbinsX();
   Double_t minvalx=hist->GetXaxis()->GetBinCenter(1);
   Double_t maxvalx=hist->GetXaxis()->GetBinCenter(nbinsx);

   Int_t nbinsy=hist->GetNbinsY();
   Double_t minvaly=hist->GetYaxis()->GetBinCenter(1);
   Double_t maxvaly=hist->GetYaxis()->GetBinCenter(nbinsy);

   Int_t nbinsz=hist->GetNbinsZ();
   Double_t minvalz=hist->GetZaxis()->GetBinCenter(1);
   Double_t maxvalz=hist->GetZaxis()->GetBinCenter(nbinsz);

   Double_t newvalx=valx;
   Double_t newvaly=valy;
   Double_t newvalz=valz;

   if(valx< minvalx)
      newvalx=minvalx;
   else if(valx>maxvalx)
      newvalx=maxvalx;

   if(valy< minvaly)
      newvaly=minvaly;
   else if(valy>maxvaly)
      newvaly=maxvaly;

   if(valz< minvalz)
      newvalz=minvalz;
   else if(valz>maxvalz)
      newvalz=maxvalz;

   hist->Fill(newvalx, newvaly, newvalz, weight);
}

