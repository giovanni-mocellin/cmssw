#include "CondTools/GEM/plugins/GEMQC8ConfESSource.h"
#include "DetectorDescription/Parser/interface/DDLParser.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "DetectorDescription/Core/interface/DDRoot.h"

#include <algorithm>


GEMQC8ConfESSource::GEMQC8ConfESSource(const edm::ParameterSet &p) :
  srcHandler(p),
  geoConfig_(p),
  rootNodeName_(p.getParameter<std::string>("rootNodeName"))
{
  setWhatProduced(this,&GEMQC8ConfESSource::produce_QC8Conf);
  findingRecord<GEMQC8ConfRcd>();
  setWhatProduced(this,&GEMQC8ConfESSource::produce_ELMap);
  findingRecord<GEMELMapRcd>();
  setWhatProduced(this,&GEMQC8ConfESSource::produce_Geom);
  findingRecord<IdealGeometryRecord>();

  edm::LogInfo("GEMQC8ConfESProducer") << "constructor" << std::endl;
}


std::unique_ptr<GEMQC8Conf> GEMQC8ConfESSource::produce_QC8Conf(const GEMQC8ConfRcd&)
{
  if (!srcHandler.getQC8conf()) {
    std::cout << "GEMQC8ConfESSource::produce_QC8Conf calling getNewObjects" << std::endl;
    srcHandler.setOnlyConfDef(1);
    srcHandler.getNewObjects();
    std::cout << "GEMQC8ConfESSource::produce_QC8Conf calling getNewObjects done" << std::endl;
    if (!srcHandler.getQC8conf()) {
      std::cout << "GEMQC8ConfESProducer::produce_QC8Conf failed to get a ptr" << std::endl;
    }
  }
  else {
    std::cout << "GEMQC8ConfESSource::produce_QC8Conf srchandler is set" << std::endl;
  }

  std::unique_ptr<GEMQC8Conf> qc8conf(new GEMQC8Conf(srcHandler.getQC8conf()));
  //(*qc8conf).assign(*srcHandler.getQC8conf());
  if (!qc8conf) std::cout << "qc8conf is null" << std::endl;
  else { std::cout << "constructed qc8conf "; qc8conf->print(std::cout); }

  std::cout << "quitting the produce_QC8Conf !!!" << std::endl;
  return qc8conf;
}


std::unique_ptr<GEMELMap> GEMQC8ConfESSource::produce_ELMap(const GEMELMapRcd&)
{
  if (!srcHandler.getQC8elMap()) {
    std::cout << "GEMQC8ConfESSource::produce_ELMap calling getNewObjects" << std::endl;
    srcHandler.setOnlyConfDef(0);
    srcHandler.getNewObjects();
    std::cout << "GEMQC8ConfESSource::produce_ELMap calling getNewObjects done" << std::endl;
    if (!srcHandler.getQC8elMap()) {
      std::cout << "GEMQC8ConfESProducer::produce_ELMap failed to get a ptr" << std::endl;
    }
  }
  else {
    std::cout << "GEMQC8ConfESSource::produce_ELMap srcHandler is set" << std::endl;
  }

  std::unique_ptr<GEMELMap> qc8elMap(new GEMELMap(srcHandler.getQC8elMap()));
  if (!qc8elMap) std::cout << "qc8elMap is null" << std::endl;
  else { std::cout << "constructed qc8elmap "; qc8elMap->print(std::cout); }

  std::cout << "quitting the produce_ELMap !!!" << std::endl;
  return qc8elMap;
}


std::unique_ptr<DDCompactView> GEMQC8ConfESSource::produce_Geom(const IdealGeometryRecord&)
{
  if (!srcHandler.getQC8conf()) {
    std::cout << "GEMQC8ConfESSource::produce_Geom calling getNewObjects" << std::endl;
    srcHandler.setOnlyConfDef(1);
    srcHandler.getNewObjects();
    std::cout << "GEMQC8ConfESSource::produce_Geom calling getNewObjects done" << std::endl;
    if (!srcHandler.getQC8conf()) {
      std::cout << "GEMQC8ConfESProducer::produce_Geom failed to get a ptr" << std::endl;
    }
  }
  else {
    std::cout << "GEMQC8ConfESSource::produce_Geom srcHandler is set" << std::endl;
  }

  // add XML files, corresponding to the set-up
  const std::vector<std::string> *chNames= & srcHandler.getQC8conf()->chSerNums();
  const std::vector<std::string> *chPos = & srcHandler.getQC8conf()->chPositions();
  //spChmbrNames_ = *chNames;
  std::vector<int> loaded;
  std::vector<char> chSize;
  loaded.reserve(15); chSize.reserve(15);
  for (unsigned int i=0; i<chNames->size(); i++) {
    char chamberSize='0';
    if (chNames->at(i).find('L') != std::string::npos) chamberSize='L';
    if (chNames->at(i).find('S') != std::string::npos) chamberSize='S';
    if (chamberSize == '0') continue;

    std::stringstream ss(chPos->at(i));
    int ir=0, ic=0;
    char c1;
    ss >> ir >> c1 >> ic;
    std::cout << "from " << chPos->at(i) << " " << ss.str() << ", got ir=" << ir << ", ic=" << ic << "\n";
    std::vector<int>::const_iterator it = std::find(loaded.begin(),loaded.end(),ir*10+ic);
    if (it != loaded.end()) {
      int idx=(it-loaded.begin());
      if (chSize.at(idx)!=chamberSize) {
	std::cout << "different chamber size! at " << ir << "," << ic << std::endl;
      }
      continue;
    }
    else {
      loaded.push_back(ir*10+ic);
      chSize.push_back(chamberSize);
      std::stringstream sout;
      sout << "Geometry/MuonCommonData/data/GEMQC8/gem11" << chamberSize
	   << "_c" << ic << "_r" << ir << ".xml";
      std::cout << "GEMQC8ConfESSource::produce_Geom: for chamber "
		<< chNames->at(i) << " adding file " << sout.str() << std::endl;
      geoConfig_.addFile(sout.str());
    }
  }

  // produce the result
  //std::cout << "\n\nproducing the result" << std::endl;

  //from XMLIdealGeometryESSource::produce();
  DDName ddName(rootNodeName_);
  DDLogicalPart rootNode(ddName);
  DDRootDef::instance().set(rootNode);
  std::unique_ptr<DDCompactView> returnValue(new DDCompactView(rootNode));
  DDLParser parser(*returnValue); //* parser = DDLParser::instance();
  parser.getDDLSAX2FileHandler()->setUserNS(false);
  int result2 = parser.parse(geoConfig_);
  if (result2 != 0) throw cms::Exception("DDException") << "DDD-Parser: parsing failed!";

  // after parsing the root node should be valid!
  if( !rootNode.isValid() ){
    throw cms::Exception("Geometry")<<"GEMQC8ConfESSource::produce_Geom: "
				    << "There is no valid node named \""
                                    <<rootNodeName_<<"\"";
  }
  returnValue->lockdown();
  std::cout << "quitting the produce_Geom !!!" << std::endl;
  return returnValue;
}


void GEMQC8ConfESSource::setIntervalFor(const edm::eventsetup::EventSetupRecordKey &,
					const edm::IOVSyncValue & iosv,
					edm::ValidityInterval & oValidity)
{
  edm::ValidityInterval infinity(iosv.beginOfTime(), iosv.endOfTime());
  oValidity = infinity;
}


#include "FWCore/Framework/interface/SourceFactory.h"
DEFINE_FWK_EVENTSETUP_SOURCE(GEMQC8ConfESSource);
