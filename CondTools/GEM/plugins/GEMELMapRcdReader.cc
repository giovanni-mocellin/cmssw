/** \class GEMELMapRcdReader
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

#include "CondFormats/GEMObjects/interface/GEMELMap.h"
#include "CondFormats/DataRecord/interface/GEMELMapRcd.h"

#include <iostream>
#include <fstream>

//using namespace std;
//using namespace edm;

// class declaration
class GEMELMapRcdReader : public edm::EDAnalyzer {
public:
  explicit GEMELMapRcdReader( const edm::ParameterSet& );
  ~GEMELMapRcdReader();
  void beginJob();
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
  void endJob();

private:
  std::string dumpFileName_;
  std::ofstream dumpFout_;
  //bool m_flag;
};

GEMELMapRcdReader::GEMELMapRcdReader( const edm::ParameterSet& iConfig ) :
  dumpFileName_(iConfig.getUntrackedParameter<std::string>("dumpFileName","")),
  dumpFout_()
{
}

GEMELMapRcdReader::~GEMELMapRcdReader(){}


void GEMELMapRcdReader::beginJob()
{
  if (dumpFileName_.size()) {
    dumpFout_.open(dumpFileName_.c_str());
  }
}


void GEMELMapRcdReader::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  std::cout << "====== GEMELMapRcdReader" << std::endl;

  edm::ESHandle<GEMELMap> hInfo;
  iSetup.get<GEMELMapRcd>().get(hInfo);
  if (!hInfo.isValid()) std::cout << "invalid handle hInfo" << std::endl;
  const GEMELMap* info=hInfo.product();
  if(!info) { std::cout << "info is null" << std::endl; return; }

  info->print(std::cout);
  std::cout << std::endl;

  if (dumpFout_.is_open()) info->print(dumpFout_,1);
}

void GEMELMapRcdReader::endJob()
{
  if (dumpFout_.is_open()) {
    dumpFout_ << "GEMELMapRcdReader: dump done\n";
    dumpFout_.close();
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(GEMELMapRcdReader);
