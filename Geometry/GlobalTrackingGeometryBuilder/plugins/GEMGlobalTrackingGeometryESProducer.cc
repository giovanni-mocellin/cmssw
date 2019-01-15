/** \file GEMGlobalTrackingGeometryESProducer.cc
 *
 *  \author Matteo Sani
 *  Jan 12, 2019 adapted for GEM QC8 by A.Juodagalvis
 */

#include "Geometry/GlobalTrackingGeometryBuilder/plugins/GEMGlobalTrackingGeometryESProducer.h"
#include "Geometry/GlobalTrackingGeometryBuilder/plugins/GlobalTrackingGeometryBuilder.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"

#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/NoProxyException.h"
#include "FWCore/Framework/interface/NoRecordException.h"

#include <memory>
#include <iostream>

using namespace edm;

GEMGlobalTrackingGeometryESProducer::GEMGlobalTrackingGeometryESProducer(const edm::ParameterSet & p){
  setWhatProduced(this);
}

GEMGlobalTrackingGeometryESProducer::~GEMGlobalTrackingGeometryESProducer(){}

std::shared_ptr<GlobalTrackingGeometry>
GEMGlobalTrackingGeometryESProducer::produce(const GlobalTrackingGeometryRecord& record) {

  std::cout << "GEMGlobalTrackingGeometryESProducer::produce" << std::endl;

  GEMGeometry const* gem = nullptr;
  ME0Geometry const *me0 = nullptr;

  try {
    edm::ESHandle<GEMGeometry> gemH;
    record.getRecord<MuonGeometryRecord>().get(gemH);
    if(gemH.isValid()) {
      gem = gemH.product();
      std::cout << "GEMGeometry ok" << std::endl;
    } else {
      std::cout << "GEMGeometry missing" << std::endl;
      LogInfo("GeometryGlobalTrackingGeometryBuilder") << "No GEM geometry is available.";
    }
  } catch (edm::eventsetup::NoRecordException<MuonGeometryRecord>& e){
    LogWarning("GeometryGlobalTrackingGeometryBuilder") << "No MuonGeometryRecord is available.";    
  }

  std::cout << "start working" << std::endl;

  GlobalTrackingGeometryBuilder builder;
  std::cout << "builder.build and return" << std::endl;
  return std::shared_ptr<GlobalTrackingGeometry>(builder.build(nullptr, nullptr, nullptr, nullptr, nullptr, gem, me0));
}

DEFINE_FWK_EVENTSETUP_MODULE(GEMGlobalTrackingGeometryESProducer);
