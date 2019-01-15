/** \class GEMQC8ConfRcdReader
 *
 *  An EDAnalyzer to print GEMQC8ConfRcd values
 *
 *  \author A. Juodagalvis - Vilnius University
 *  Oct 2018
 */

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/GEMObjects/interface/GEMQC8Conf.h"
#include "CondFormats/DataRecord/interface/GEMQC8ConfRcd.h"

#include <iostream>
#include <fstream>

//using namespace std;
//using namespace edm;

// class declaration
class GEMQC8ConfRcdReader : public edm::EDAnalyzer {
public:
  explicit GEMQC8ConfRcdReader( const edm::ParameterSet& );
  ~GEMQC8ConfRcdReader();
  void beginJob();
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
  void endJob();

private:
  std::string dumpFileName_;
  std::ofstream dumpFout_;
  //bool m_flag;
};

GEMQC8ConfRcdReader::GEMQC8ConfRcdReader( const edm::ParameterSet& iConfig ) :
  dumpFileName_(iConfig.getUntrackedParameter<std::string>("dumpFileName","")),
  dumpFout_()
{
}

GEMQC8ConfRcdReader::~GEMQC8ConfRcdReader(){}


void GEMQC8ConfRcdReader::beginJob()
{
  if (dumpFileName_.size()) {
    dumpFout_.open(dumpFileName_.c_str());
  }
}


void GEMQC8ConfRcdReader::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  std::cout << "====== GEMQC8ConfRcdReader" << std::endl;

  edm::ESHandle<GEMQC8Conf> hInfo;
  iSetup.get<GEMQC8ConfRcd>().get(hInfo);
  if (!hInfo.isValid()) std::cout << "invalid handle hInfo" << std::endl;
  const GEMQC8Conf* info=hInfo.product();
  if(!info) { std::cout << "info is null" << std::endl; return; }

  info->print(std::cout);
  std::cout << std::endl;

  if (dumpFout_.is_open()) info->print(dumpFout_);
}

void GEMQC8ConfRcdReader::endJob()
{
  if (dumpFout_.is_open()) {
    dumpFout_ << "GEMQC8ConfRcdReader: dump done\n";
    dumpFout_.close();
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(GEMQC8ConfRcdReader);
