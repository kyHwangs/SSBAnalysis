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
    double etaLow, etaHi, ptLow, ptHi, effi, effiPre, effiPost;

    bool belongTo(double, double) const;
};

struct record_11
{
    double mu1etaLow, mu1etaHi, mu1ptLow, mu1ptHi, mu2etaLow, mu2etaHi, mu2ptLow, mu2ptHi, effi, effiPre, effiPost;

    bool belongTo_forTrigger(double, double, double, double) const;
};

std::vector<record> _recd;
std::vector<record_11> _recd_11;

public:

	explicit EffTable() = default;

	explicit EffTable(const std::string &fFilename_) 
	:fFilename(fFilename_) {
		std::cout << "loading efficiency SF table of " << fFilename_ << std::endl;
		init();
	};

	double getEfficiency(double pt, double eta) const;
	double getTriggerEfficiency(double mu1pt, double mu1eta, double mu2pt, double mu2eta) const;

private:

	void init();

	std::string fFilename;

};
