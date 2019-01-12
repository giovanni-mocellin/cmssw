#ifndef GEMGeometry_GEMGeometryPrint_h
#define GEMGeometry_GEMGeometryPrint_h

/** A helper unit to print GEMGeometry classes
 *
 *  \author A. Juodagalvis - Vilnius University, Sep 2018
 */

#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/GEMGeometry/interface/GEMChamber.h"
#include "Geometry/GEMGeometry/interface/GEMSuperChamber.h"
#include "Geometry/GEMGeometry/interface/GEMRing.h"
#include "Geometry/GEMGeometry/interface/GEMStation.h"
#include "Geometry/GEMGeometry/interface/GEMRegion.h"
#include <vector>
#include <map>
#include <iostream>


typedef TrackingGeometry::DetTypeContainer DetTypeContainer; // v<const GeomDetType*>
typedef TrackingGeometry::DetContainer     DetContainer; // v<const GeomDet*>
typedef TrackingGeometry::DetIdContainer   DetIdContainer; // v<DetId>
typedef TrackingGeometry::mapIdToDet       mapIdToDet; // std::unordered_map<unsigned int, const GeomDet*>


void printGEMGeometry(const GEMGeometry &gg, std::ostream &out = std::cout);
void printPtrAddr(const void *ptr, std::ostream &out = std::cout, const char *preMsg=NULL, const char *postMsg=NULL);

//void printDetContainer(const DetContainer &dc, std::ostream &out = std::cout, int detailed=1);
//void printDetContainer(const std::vector<const GeomDet*> &dc, std::ostream &out = std::cout, int detailed=1);
void printAllGeomDets(const std::vector<const GeomDet*> &dets, std::ostream &out = std::cout, int detailed=1);
void printGeomDet(const GeomDet &gd, std::ostream &out = std::cout);
void printDetTypeContainer(const DetTypeContainer &dtc, std::ostream &out = std::cout, int detailed=1);
void printAllGeomDetTypes(const std::vector<const GeomDetType*> &gdts, std::ostream &out = std::cout, int detailed=1);
void printGeomDetType(const GeomDetType &gdt, std::ostream &out = std::cout);
void printDetIdContainer(const DetIdContainer &ids, std::ostream &out = std::cout, int detailed=1);
void printAllDetIds(const std::vector<DetId> &detIds, std::ostream &out = std::cout, int detailed=1);
void printMapIdToDet(const mapIdToDet &mp, std::ostream &out = std::cout, int detailed=1);
void printAllMapUIntToGeomDet(const std::unordered_map<unsigned int, const GeomDet*> &mp, std::ostream &out = std::cout, int detailed=1);
void printAllGEMEtaPartitions(const std::vector<const GEMEtaPartition*> &etaPartitions, std::ostream &out = std::cout, int detailed=1);
void printGEMEtaPartition(const GEMEtaPartition &etaPartition, std::ostream &out = std::cout, int detailed=1);
void printGEMEtaPartitionSpecs(const GEMEtaPartitionSpecs &eps, std::ostream &out = std::cout, int detailed=1);
void printAllGEMChambers(const std::vector<const GEMChamber*> &chambers, std::ostream &out = std::cout, int detailed=1);
void printGEMChamber(const GEMChamber &chamber, std::ostream &out = std::cout, int detailed=1);
void printAllGEMSuperChambers(const std::vector<const GEMSuperChamber*> &spchambers, std::ostream &out = std::cout, int detailed=1);
void printGEMSuperChamber(const GEMSuperChamber &spch, std::ostream &out = std::cout, int detailed=1);
void printAllGEMRings(const std::vector<const GEMRing*> &rings, std::ostream &out = std::cout, int detailed=1);
void printGEMRing(const GEMRing &ring, std::ostream &out = std::cout, int detailed=1);
void printGEMDetIds(const std::vector<GEMDetId> &ids, std::ostream &out = std::cout, int detailed=1);
void printAllGEMStations(const std::vector<const GEMStation*> &stations, std::ostream &out = std::cout, int detailed=1);
void printGEMStation(const GEMStation &st, std::ostream &out = std::cout, int detailed=1);
void printAllGEMRegions(const std::vector<const GEMRegion*> &regions, std::ostream &out = std::cout, int detailed=1);
void printGEMRegion(const GEMRegion &reg, std::ostream &out = std::cout, int detailed=1);

#endif
