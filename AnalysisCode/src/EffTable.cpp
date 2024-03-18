#include "./../interface/EffTable.hpp"

#include <cmath>
#include <fstream>
#include <sstream>

void EffTable::init() {

  std::cout << " Loading " << fFilename << std::endl;
  std::ifstream file(fFilename);

  if (!file) {
    std::cout << " File " << fFilename << " does not exist" << std::endl;
    throw std::invalid_argument("File " + fFilename + " doesn't exist");
  }


  while (file) {
    std::vector<double> data;
    std::string oneLine;

    std::getline(file, oneLine);
    std::stringstream ss(oneLine);

    std::string value;

    while(std::getline(ss, value, ' ')) {

        double tmpDouble = 0;
        std::stringstream ssDouble(value);
        ssDouble >> tmpDouble;
        data.push_back(tmpDouble);
    }

    if (data.size() == 7)  _recd.push_back(record{data[0], data[1], data[2], data[3], data[4], data[5], data[6]});
    if (data.size() == 11) _recd_11.push_back(record_11{data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10]});
  }

  std::cout << "Eff table low : " << _recd.size() << std::endl;
  std::cout << "Eff table11 low : " << _recd_11.size() << std::endl;

}

bool EffTable::record::belongTo(double pt, double eta) const
{
    return (pt < ptHi && pt >= ptLow) && (eta < etaHi && eta >= etaLow);
}

bool EffTable::record_11::belongTo_forTrigger(double mu1pt, double mu1eta, double mu2pt, double mu2eta) const
{
    return (mu1pt < mu1ptHi && mu1pt >= mu1ptLow) && (mu1eta < mu1etaHi && mu1eta >= mu1etaLow) && 
           (mu2pt < mu2ptHi && mu2pt >= mu2ptLow) && (mu2eta < mu2etaHi && mu2eta >= mu2etaLow);
}

double EffTable::getEfficiency(double pt, double eta) const
{
  double hiPtBin = 0;
  
  for (unsigned int i = 0; i != _recd.size(); i++)
    if ((_recd[i]).belongTo(pt, eta))
      return _recd[i].effi;

  for (unsigned int i = 0; i != _recd.size(); i++)
    if ((_recd[i]).belongTo(0.5 * (_recd[i].ptHi + _recd[i].ptLow), eta))
      hiPtBin = _recd[i].effi;

  return hiPtBin;
}

double EffTable::getTriggerEfficiency(double mu1pt, double mu1eta, double mu2pt, double mu2eta) const
{
  double binEff = 0;

  for (unsigned int i = 0; i != _recd_11.size(); i++)
    if ((_recd_11[i]).belongTo_forTrigger(mu1pt, mu1eta, mu2pt, mu2eta))
      return _recd_11[i].effi;

  for (unsigned int i = 0; i != _recd_11.size(); i++)
    if ((_recd_11[i]).belongTo_forTrigger(0.5 * (_recd_11[i].mu1ptHi + _recd_11[i].mu1ptLow), mu1eta, 0.5 * (_recd_11[i].mu2ptHi + _recd_11[i].mu2ptLow), mu2eta))
      binEff = _recd_11[i].effi;

  return binEff;
}
