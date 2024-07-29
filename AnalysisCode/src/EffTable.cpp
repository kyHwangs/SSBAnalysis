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

    if (data.size() == 7) _recd.push_back(record{data[0], data[1], data[2], data[3], data[4], data[5], data[6]});
  }

  std::cout << "Eff table low : " << _recd.size() << std::endl;

}

bool EffTable::record::belongTo(double pt, double eta)
{
    return (pt < ptHi && pt >= ptLow) && (eta < etaHi && eta >= etaLow);
}


double EffTable::getEfficiency(double pt, double eta)
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

