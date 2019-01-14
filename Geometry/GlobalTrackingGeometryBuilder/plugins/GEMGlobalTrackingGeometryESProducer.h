#ifndef Geometry_GEMGeometryBuilder_GEMGlobalTrackingGeometryESProducer_h
#define Geometry_GEMGeometryBuilder_GEMGlobalTrackingGeometryESProducer_h

/** \class GlobalTrackingGeometry
 * 
 *  ESProducer for GlobalTrackingGeometry in MuonGeometryRecord
 *
 *  \author Matteo Sani
 *  Jan 12, 2019 adapted for GEM QC8 by Andrius Juodagalvis
 */

#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

#include <memory>
#include <string>

class GlobalTrackingGeometry;

class GEMGlobalTrackingGeometryESProducer : public edm::ESProducer {
public:

  /// Constructor
  GEMGlobalTrackingGeometryESProducer(const edm::ParameterSet & p);

  /// Destructor
  ~GEMGlobalTrackingGeometryESProducer() override;

  /// Produce GlobalTrackingGeometry
  std::shared_ptr<GlobalTrackingGeometry> produce(const GlobalTrackingGeometryRecord& record);

private:  

};
#endif






