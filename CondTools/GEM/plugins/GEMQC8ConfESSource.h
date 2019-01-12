#ifndef CondTools_GEM_plugins_GEMQC8ConfESSource_h
#define CondTools_GEM_plugins_GEMQC8ConfESSource_h

/** \class GEMQC8ConfESSource
 *
 *  A class to produce 3 records required for GEMQC8 stand:
 *  - chamber setup (GEMQC8ConfRcd)
 *  - geometrical definition (IdealGeometryRecord)
 *  - electronics map (GEMELMapRcd)
 *
 *  \author A. Juodagalvis - Vilnius University
 *  Oct 2018
 */

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ESProducts.h"
#include "FWCore/Framework/interface/EventSetupRecordIntervalFinder.h"
#include "CondTools/GEM/interface/GEMQC8ConfSourceHandler.h"
#include "CondFormats/GEMObjects/interface/GEMQC8Conf.h"
#include "CondFormats/DataRecord/interface/GEMQC8ConfRcd.h"
#include "CondFormats/GEMObjects/interface/GEMELMap.h"
#include "CondFormats/DataRecord/interface/GEMELMapRcd.h"

#include "GeometryReaders/XMLIdealGeometryESSource/interface/GeometryConfiguration.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"


class GEMQC8ConfESSource : public edm::ESProducer, public edm::EventSetupRecordIntervalFinder {
 public:
  GEMQC8ConfESSource(const edm::ParameterSet&);
  ~GEMQC8ConfESSource() override {}

  std::unique_ptr<GEMQC8Conf> produce_QC8Conf(const GEMQC8ConfRcd&);
  std::unique_ptr<GEMELMap>   produce_ELMap(const GEMELMapRcd&);
  std::unique_ptr<DDCompactView> produce_Geom(const IdealGeometryRecord&);

 protected:
  void setIntervalFor(const edm::eventsetup::EventSetupRecordKey &,
		      const edm::IOVSyncValue&, edm::ValidityInterval &) override;

 private:
  popcon::GEMQC8ConfSourceHandler srcHandler;

  // geometrical description of superchambers
  GeometryConfiguration geoConfig_;
  std::string rootNodeName_;
};



#endif
