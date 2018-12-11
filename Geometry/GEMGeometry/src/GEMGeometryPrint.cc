/** Implementation of print GEMGeometry
 *
 *  \author A. Juodagalvis - Vilnius University, Sep 2018
 */

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/GEMGeometry/interface/GEMGeometryPrint.h"

void printGEMGeometry(const GEMGeometry &gg, std::ostream &out)
{
  out << "GEMGeometry" << std::endl;
  out << " theEtaPartitions: ";
  printAllGeomDets(gg.detUnits(),out,0);
  out << " theDets: ";
  printAllGeomDets(gg.dets(),out,0);
  out << " theEtaPartitionTypes: ";
  printDetTypeContainer(gg.detTypes(),out,0);
  out << " theEtaPartitionIds: ";
  printDetIdContainer(gg.detUnitIds(),out);
  out << " theDetIds: ";
  printDetIdContainer(gg.detIds(),out);
  out << " theMap: ";
  printMapIdToDet(gg.getIdToDetMap(),out,0);
  out << " allEtaPartitions: ";
  printAllGEMEtaPartitions(gg.etaPartitions(),out);
  out << " allChambers: ";
  printAllGEMChambers(gg.chambers(),out);
  out << " allSuperChambers: ";
  printAllGEMSuperChambers(gg.superChambers(),out);
  out << " allRings: ";
  printAllGEMRings(gg.rings(),out);
  out << " allStations: ";
  printAllGEMStations(gg.stations(),out);
  out << " allRegions: ";
  printAllGEMRegions(gg.regions(),out);
  out << "GEMGeometry END" << std::endl;
}

void printPtrAddr(const void *ptr, std::ostream &out, const char *preMsg, const char *postMsg)
{
  if (preMsg) out << preMsg;
  out << std::hex << ptr << std::dec;
  if (postMsg) out << postMsg;
}

/*  
void printDetContainter(const DetContainer &dc, std::ostream &out, int detailed)
{
  out << " [DetContainer] ";
  printAllGeomDets(dc,out,detailed);
}

void printDetContainter(const std::vector<const GeomDet*> &dets, std::ostream &out, int detailed)
{
  out << " [DetContainer] ";
  printAllGeomDets(dets,out,detailed);
}
*/

void printAllGeomDets(const std::vector<const GeomDet*> &dets, std::ostream &out, int detailed)
{
  out << " printAllGeomDets size[" << dets.size() << "]\n";
  if (detailed) out << " (verb=" << detailed << ") ";
  for (unsigned int i=0; i<dets.size(); i++) {
    out << "(" << i << ") ";
    if (dets[i]==NULL) out << "NULL ";
    else {
      if (detailed==2) printPtrAddr(dets[i],out,NULL," ");
      if (detailed) printGeomDet(*dets[i],out);
      else out << "gid=" << dets[i]->geographicalId().rawId() << "\n";
    }
    out << "\n";
  }
}

void printGeomDet(const GeomDet &gd, std::ostream &out)
{
  out << "printGeomDet not ready\n";
}

void printDetTypeContainer(const DetTypeContainer &dtc, std::ostream &out, int detailed)
{
  out << " [DetTypeContainer] ";
  printAllGeomDetTypes(dtc,out,detailed);
}

void printAllGeomDetTypes(const std::vector<const GeomDetType*> &gdts, std::ostream &out, int detailed)
{
  out << " printAllGeomDetTypes size[" << gdts.size() << "]\n";
  if (detailed) out << " (verb=" << detailed << ") ";
  for (unsigned int i=0; i<gdts.size(); i++) {
    out << "(" << i << ") ";
    if (gdts[i]==NULL) out << "NULL ";
    else {
      if (detailed==2) printPtrAddr(gdts[i],out,NULL," ");
      if (detailed) printGeomDetType(*gdts[i],out);
      else out << gdts[i]->name() << " " << "\n";
    }
  }
  out << "\n";
}

void printGeomDetType(const GeomDetType &gdt, std::ostream &out)
{
  out << "printGeomDetType not ready\n";
}

void printDetIdContainer(const DetIdContainer &ids, std::ostream &out, int detailed)
{
  out << " [DetIdContainer] ";
  printAllDetIds(ids,out,detailed);
}

void printAllDetIds(const std::vector<DetId> &detIds, std::ostream &out, int detailed)
{
  out << " printAllDetIds size[" << detIds.size() << "]: ";
  if (detailed) out << " (verb=" << detailed << ") ";
  for (unsigned int i=0; i<detIds.size(); i++) {
    out << " " << detIds[i].rawId();
  }
  out << "\n";
}

void printMapIdToDet(const mapIdToDet &mp, std::ostream &out, int detailed)
{
  out << " [mapIdToDet] ";
  printAllMapUIntToGeomDet(mp,out,detailed);
}

void printAllMapUIntToGeomDet(const std::unordered_map<unsigned int, const GeomDet*> &mp, std::ostream &out, int detailed)
{
  out << " printAllMapUIntToGeomDet size[" << mp.size() << "]: ";
  if (detailed) out << " (verb=" << detailed << ") ";
  int idx=0;
  for (std::unordered_map<unsigned int, const GeomDet*>::const_iterator it=mp.begin(); it!=mp.end(); it++, idx++) {
    out << "(" << idx << ") ";
    out << it->first << ", ";
    if (detailed==2) printPtrAddr(it->second,out,NULL," ");
    if (detailed) printGeomDet(*(it->second),out);
    else out << " rawId=" << it->second->geographicalId().rawId() << "\n";
  }
  out << "\n";
}

void printAllGEMEtaPartitions(const std::vector<const GEMEtaPartition*> &etaPartitions, std::ostream &out, int detailed)
{
  out << " printAllGEMEtaPartitions size[" << etaPartitions.size() << "]\n";
  for (unsigned int i=0; i<etaPartitions.size(); i++) {
    out << "(" << i << ") ";
    printGEMEtaPartition(*etaPartitions[i],out,detailed);
  }
  out << "\n";
}

void printGEMEtaPartition(const GEMEtaPartition &etaPartition, std::ostream &out, int detailed)
{
  out << " GEMEtaPartition ";
  if (detailed) out << " (verb=" << detailed << ") ";
  out << " GEMDetId=" << etaPartition.id() << " ";
  printGeomDet(*(GeomDet*)&etaPartition,out);
  printGEMEtaPartitionSpecs(*etaPartition.specs(),out,detailed);
}

void printGEMEtaPartitionSpecs(const GEMEtaPartitionSpecs &eps, std::ostream &out, int detailed)
{
  out << " GEMEtaPartitionSpecs ";
  for (unsigned int i=0; i<eps.parameters().size(); i++) {
    out << " " << eps.parameters().at(i);
  }
  out << "\n";
}

void printAllGEMChambers(const std::vector<const GEMChamber*> &chambers, std::ostream &out, int detailed)
{
  out << " printAllGEMChambers size[" << chambers.size() << "]";
  if (detailed) out << " (verb=" << detailed << ")";
  out << "\n";
  for (unsigned int i=0; i<chambers.size(); i++) {
    out << "(" << i << ") ";
    if (detailed==2) printPtrAddr(chambers.at(i),out,NULL," ");
    if (detailed) printGEMChamber(*chambers.at(i),out,detailed);
    else out << "id=" << chambers.at(i)->id();
  }
  out << "\n";
}

void printGEMChamber(const GEMChamber &chamber, std::ostream &out, int detailed)
{
  out << "GEMChamber id=" << chamber.id() << " ";
  out << "  has " << chamber.etaPartitions().size() << " partitions\n";
  printAllGEMEtaPartitions(chamber.etaPartitions(),out,detailed);
}

void printAllGEMSuperChambers(const std::vector<const GEMSuperChamber*> &spchambers, std::ostream &out, int detailed)
{
  out << " printAllGEMSuperChambers size[" << spchambers.size() << "]\n";
  if (detailed) out << " (verb=" << detailed << ") ";
  for (unsigned int i=0; i<spchambers.size(); i++) {
    out << "(" << i << ") ";
    if (spchambers[i]==NULL) out << "NULL ";
    else {
      if (detailed==2) printPtrAddr(spchambers[i],out,NULL, " ");
      if (detailed) printGEMSuperChamber(*spchambers[i],out,detailed);
      else {
	out << spchambers[i]->geographicalId().rawId() << " ";
	printGEMDetIds(spchambers[i]->ids(),out,detailed);
      }
    }
  }
  out << "\n";
}

void printGEMSuperChamber(const GEMSuperChamber &spch, std::ostream &out, int detailed)
{
  out << "printGEMSuperChamber id=" << spch.id()
      << ", geoId=" << spch.geographicalId().rawId()
      << ", with " << spch.chambers().size() << " chambers\n";
  if (detailed) {
    printAllGEMChambers(spch.chambers(),out,detailed);
  }
}

void printAllGEMRings(const std::vector<const GEMRing*> &rings, std::ostream &out, int detailed)
{
  out << "printAllGEMRings size[" << rings.size() << "]\n";
  if (detailed) out << " (verb=" << detailed << ")\n";
  for (unsigned int i=0; i<rings.size(); i++) {
    out << "(" << i << ") ";
    if (rings[i]==NULL) out << "NULL ";
    else {
      if (detailed==2) printPtrAddr(rings[i],out,NULL," ");
      if (detailed) printGEMRing(*rings[i],out,detailed);
      else {
	out << "region=" << rings[i]->region() << " station="
	    << rings[i]->station() << " ring=" << rings[i]->ring()
	    << "  nSuperCh=" << rings[i]->superChambers().size() << "\n";
	printGEMDetIds(rings[i]->ids(),out,0);
      }
    }
  }
  out << "\n";
}

void printGEMRing(const GEMRing &ring, std::ostream &out, int detailed)
{
  out << "GEMRing ";
  out << "region=" << ring.region() << " station="
      << ring.station() << " ring=" << ring.ring()
      << "  nSuperCh=" << ring.superChambers().size() << "\n";
  printGEMDetIds(ring.ids(),out,0);
  out << " superChambers:\n";
  for (unsigned int i=0; i<ring.superChambers().size(); i++) {
    out << "(" << i << ") ";
    if (detailed==2) printPtrAddr(ring.superChambers()[i],out,NULL," ");
    if (detailed) printGEMSuperChamber(*ring.superChambers()[i],out,detailed);
    else {
      out << " " << ring.superChambers()[i]->id() << "\n";
    }
  }
  out << "\n";
}

void printGEMDetIds(const std::vector<GEMDetId> &ids, std::ostream &out, int detailed)
{
  out << "printGEMDetIds size[" << ids.size() << "]: ";
  for (unsigned int i=0; i<ids.size(); i++) {
    out << " " << ids[i].rawId();
  }
  out << "\n";
}

void printAllGEMStations(const std::vector<const GEMStation*> &stations, std::ostream &out, int detailed)
{
  out << "printAllGEMStations size[" << stations.size() << "]\n";
  if (detailed) out << " (verb=" << detailed << ") ";
  for (unsigned int i=0; i<stations.size(); i++) {
    out << "(" << i << ") ";
    if (stations[i]==NULL) out << "NULL ";
    else {
      if (detailed==2) printPtrAddr(stations[i],out,NULL," ");
      if (detailed) printGEMStation(*stations[i],out,detailed);
      else {
	out << "name=" << stations[i]->getName()
	    << " region=" << stations[i]->region()
	    << " station=" << stations[i]->station()
	    << " nRings=" << stations[i]->rings().size()
	    << "\n";
      }
    }
  }
  out << "\n";
}

void printGEMStation(const GEMStation &station, std::ostream &out, int detailed)
{
  out << "printGEMStation ";
  if (detailed) out << " (verb=" << detailed << ") ";
  out << "name=" << station.getName()
      << " region=" << station.region()
      << " station=" << station.station()
      << " nRings=" << station.rings().size()
      << "\n";
  printAllGEMRings(station.rings(),out,detailed);
}

void printAllGEMRegions(const std::vector<const GEMRegion*> &regions, std::ostream &out, int detailed)
{
  out << "printAllGEMRegions size[" << regions.size() << "] ";
  if (detailed) out << " (verb=" << detailed << ")\n";
  for (unsigned int i=0; i<regions.size(); i++) {
    out << "(" << i << ") ";
    if (detailed==2) printPtrAddr(regions[i],out,NULL," ");
    if (detailed) printGEMRegion(*regions[i],out,detailed);
    else {
      out << " region=" << regions[i]->region() << ", nStations="
	  << regions[i]->stations().size() << "\n";
    }
    out << "\n";
  }
}

void printGEMRegion(const GEMRegion &reg, std::ostream &out, int detailed)
{
  out << "printGEMRegion region=" << reg.region() << ", nStations="
      << reg.stations().size() << "\n";
  if (detailed) out << " (verb=" << detailed << ")\n";
  printAllGEMStations(reg.stations(),out,detailed);
}
