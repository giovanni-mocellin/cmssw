#include "CondFormats/GEMObjects/interface/GEMROmap.h"

void GEMROmap::printElDetMap(std::ostream &out) const
{
  out << "(amcId,gebId,vfatId) -> (vfatType,gemDetId,locPhi)\n";
  int idx=0;
  for (auto it= roMapED_.begin(); it!=roMapED_.end(); it++) {
    out << "idx=" << idx << " (" << it->first.amcId << "," << it->first.gebId
	<< "," << it->first.vfatId << ") -> (" << it->second.vfatType
	<< "," << it->second.gemDetId << "," << it->second.locPhi << ")\n";
    idx++;
  }
  return;
}


void GEMROmap::printDetElMap(std::ostream &out) const
{
  out << "(vfatType,gemDetId,locPhi) -> (amcId,gebId,vfatId)\n";
  int idx=0;
  for (auto it= roMapDE_.begin(); it!=roMapDE_.end(); it++) {
    out << "idx=" << idx << " (" << it->first.vfatType
	<< "," << it->first.gemDetId << "," << it->first.locPhi << ") -> ("
	<< it->second.amcId << "," << it->second.gebId
	<< "," << it->second.vfatId << ")\n";
    idx++;
  }
  return;
}
