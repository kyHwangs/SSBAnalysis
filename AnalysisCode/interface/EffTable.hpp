#include <set>
#include <string>
#include <fstream>
#include <iostream>
#include <cassert>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <TH1.h>
#include <TH2.h>

class EffTable
{

struct record
{
    double ptLow, ptHi, etaLow, etaHi, effi, effSystDown, effSystUp;

    bool belongTo(double, double) const;
};

std::vector<record> _recd;

public:

	explicit EffTable() = default;

	explicit EffTable(const std::string &fFilename_) 
	:fFilename(fFilename_) {
		std::cout << "loading efficiency SF table of " << fFilename_ << std::endl;
		init();
	};

	double getEfficiency(double pt, double eta) const;

private:

	void init();

	std::string fFilename;

};
