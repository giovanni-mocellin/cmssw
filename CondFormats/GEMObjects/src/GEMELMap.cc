#ifndef noFileInPath_H
#include "CondFormats/GEMObjects/interface/GEMELMap.h"
#include "CondFormats/GEMObjects/interface/GEMROmap.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#endif

GEMELMap::GEMELMap():
  theVFatMap_(), theStripMap_(),
  theVersion("") {}

GEMELMap::GEMELMap(const GEMELMap *ptr) :
  theVFatMap_(), theStripMap_(),
  theVersion("")
{
  if (!ptr) return;
  theVFatMap_ = ptr->theVFatMap_;
  theStripMap_ = ptr->theStripMap_;
  theVersion = ptr->theVersion;
}

GEMELMap::GEMELMap(const std::string & version):
  theVFatMap_(), theStripMap_(),
  theVersion(version) {}

GEMELMap::~GEMELMap() {}

const std::string & GEMELMap::version() const {
  return theVersion;
}

#ifndef noFileInPath_H
void GEMELMap::convert(GEMROmap & romap) {
  // amc->geb->vfat mapping to GEMDetId
  for (auto imap : theVFatMap_) {
    for (unsigned int ix=0;ix<imap.vfatId.size();ix++) {
      // vfat v2
      if (imap.vfatType[ix] < vfatTypeV3_) {
        GEMROmap::eCoord ec;
        ec.vfatId = imap.vfatId[ix] & chipIdMask_;
        ec.gebId = imap.gebId[ix];
        ec.amcId = imap.amcId[ix];
        
        int st = std::abs(imap.z_direction[ix]);
        GEMROmap::dCoord dc;
        dc.gemDetId = GEMDetId(imap.z_direction[ix], 1, st, imap.depth[ix], imap.sec[ix], imap.iEta[ix]);
        dc.vfatType = imap.vfatType[ix];
        dc.locPhi = (imap.iPhi[ix]-1)%maxVFatGE11_;
        romap.add(ec,dc);
        romap.add(dc,ec);
        
        GEMROmap::eCoord ecGEB;
        ecGEB.vfatId = 0;
        ecGEB.gebId = ec.gebId;
        ecGEB.amcId = ec.amcId;
        // this is for backward compatablity for v2
        if (!romap.isValidChipID(ecGEB)) {
          GEMROmap::dCoord dcGEB;
          dcGEB.gemDetId = dc.gemDetId.chamberId();
          dcGEB.vfatType = dc.vfatType;
          romap.add(ecGEB,dcGEB);
          romap.add(dcGEB,ecGEB);
        }
      }
      // vfat v3
      else {
        GEMROmap::eCoord ec;
        ec.vfatId = imap.vfatId[ix]; // vfatId is chip position in v3
        ec.gebId = imap.gebId[ix];
        ec.amcId = imap.amcId[ix];
        
        int st = std::abs(imap.z_direction[ix]);
        GEMROmap::dCoord dc;
        dc.gemDetId = GEMDetId(imap.z_direction[ix], 1, st, imap.depth[ix], imap.sec[ix], imap.iEta[ix]);
        dc.vfatType = imap.vfatType[ix];
        int maxVFat = maxVFatGE11_;
        if (st == 2) maxVFat = maxVFatGE21_;
        if (st == 0) maxVFat = maxVFatGE0_;
        dc.locPhi = (imap.iPhi[ix]-1)%maxVFat;
        romap.add(ec,dc);
        romap.add(dc,ec);
      }
    }
  }

  // channel mapping
  for (auto imap : theStripMap_) {
    for (unsigned int ix=0;ix<imap.vfatType.size();ix++) {
      GEMROmap::channelNum cMap;
      cMap.vfatType = imap.vfatType[ix];
      cMap.chNum = imap.vfatCh[ix];

      GEMROmap::stripNum sMap;
      sMap.vfatType = imap.vfatType[ix];
      sMap.stNum = imap.vfatStrip[ix];

      romap.add(cMap, sMap);
      romap.add(sMap, cMap);
    }
  }

  romap.printElDetMap(std::cout);
  romap.printDetElMap(std::cout);
}

void GEMELMap::convertDummy(GEMROmap & romap) {
  // 12 bits for vfat, 5 bits for geb, 8 bit long GLIB serial number
  uint16_t amcId = 1; //amc
  uint16_t gebId = 0; 

  for (int re = -1; re <= 1; re = re+2) {
    for (int st = GEMDetId::minStationId; st<=GEMDetId::maxStationId; ++st) {
      int maxVFat = maxVFatGE11_;
      if (st == 2) maxVFat = maxVFatGE21_;
      if (st == 0) maxVFat = maxVFatGE0_;

      for (int ch = 1; ch<=GEMDetId::maxChamberId; ++ch) {
        for (int ly = 1; ly<=GEMDetId::maxLayerId; ++ly) {
          // 1 geb per chamber
          gebId++;
          uint16_t chipPos = 0;
          for (int lphi = 0; lphi < maxVFat; ++lphi) {
            for (int roll = 1; roll<=maxEtaPartition_; ++roll) {

              GEMDetId gemId(re, 1, st, ly, ch, roll);

              GEMROmap::eCoord ec;
              ec.vfatId = chipPos;
              ec.gebId = gebId;
              ec.amcId = amcId;

              GEMROmap::dCoord dc;
              dc.gemDetId = gemId;
              dc.vfatType = vfatTypeV3_;// > 10 is vfat v3
              dc.locPhi = lphi;

              romap.add(ec,dc);
              romap.add(dc,ec);
              chipPos++;
            }
          }

          // 5 bits for gebId 
          if (gebId == maxGEBs_) {
            // 24 gebs per amc
            gebId = 0;
            amcId++;
          }
        }
      }
    }
  }

  for (int i = 0; i < maxChan_; ++i) {
    // only 1 vfat type for dummy map
    GEMROmap::channelNum cMap;
    cMap.vfatType = vfatTypeV3_;
    cMap.chNum = i;

    GEMROmap::stripNum sMap;
    sMap.vfatType = vfatTypeV3_;
    sMap.stNum = i;

    romap.add(cMap, sMap);
    romap.add(sMap, cMap);
  }
}
#endif


void GEMELMap::print(std::ostream &out, int detailed, int checkArrays) const
{
  out << "GEMELMap: theVFATMap[" << theVFatMap_.size() << "], theStripMap["
      << theStripMap_.size() << "]\n";
  if (detailed) {
    out << "  theVFATMap[" << theVFatMap_.size() << "]:\n";
    for (unsigned int i=0; i<theVFatMap_.size(); i++) {
      const GEMVFatMap *ptr= & theVFatMap_.at(i);
      unsigned int sz= ptr->vfat_position.size();
      out << " i=" << i << ", VFATmapTypeId=" << ptr->VFATmapTypeId << ", sz="
	  << sz << "\n";
      int ok=1;
      if (checkArrays) {
	if ((sz!=ptr->z_direction.size()) || (sz!=ptr->iEta.size()) || (sz!=ptr->iPhi.size()) || (sz!=ptr->depth.size()) || (sz!=ptr->vfatType.size()) || (sz!=ptr->vfatId.size()) || (sz!=ptr->amcId.size()) || (sz!=ptr->gebId.size()) || (sz!=ptr->sec.size())) {
	  out << "GEMVFatMap size error : z_direction[" << ptr->z_direction.size() << "], iEta[" << ptr->iEta.size() << "], iPhi[" << ptr->iPhi.size() << "], depth[" << ptr->depth.size() << "], vfatType[" << ptr->vfatType.size() << "], vfatId[" << ptr->vfatId.size() << "], amcId[" << ptr->amcId.size() << "], gebId[" << ptr->gebId.size() << "], sec[" << ptr->sec.size() << "\n";
	  ok=0;
	}
      }
      if (ok) {
	for (unsigned int ii=0; ii<sz; ii++) {
	  out << " " << ii << " vfat_pos=" << ptr->vfat_position[ii];
	  if (ii < ptr->z_direction.size()) out << " z_dir=" << ptr->z_direction[ii];
	  if (ii < ptr->iEta.size()) out << " iEta=" << ptr->iEta[ii];
	  if (ii < ptr->iPhi.size()) out << " iPhi=" << ptr->iPhi[ii];
	  if (ii < ptr->depth.size()) out << " depth=" << ptr->depth[ii];
	  if (ii < ptr->vfatType.size()) out << " vfatType=" << ptr->vfatType[ii];
	  if (ii < ptr->vfatId.size()) out << " vfatId=0x" << std::hex << ptr->vfatId[ii] << std::dec;
	  if (ii < ptr->amcId.size()) out << " amcId=" << ptr->amcId[ii];
	  if (ii < ptr->gebId.size()) out << " gebId=" << ptr->gebId[ii];
	  if (ii < ptr->sec.size()) out << " sect=" << ptr->sec[ii];
	  out << "\n";
	}
      }
    } // for theVFatMap_

    out << " theStripMap[" << theStripMap_.size() << "]:\n";
    for (unsigned int i=0; i<theStripMap_.size(); i++) {
      const GEMStripMap *ptr = & theStripMap_.at(i);
      unsigned int sz= ptr->vfatType.size();
      out << " i=" << i << ", sz=" << sz << "\n";
      int ok=1;
      if (checkArrays) {
	if ((sz!=ptr->vfatCh.size()) || (sz!=ptr->vfatStrip.size())) {
	  out << "GEMStripMap size error : vfatType[" << ptr->vfatType.size() << "], vfatCh[" << ptr->vfatCh.size() << "], vfatStrip[" << ptr->vfatStrip.size() << "]\n";
	  ok=0;
	}
      }
      if (ok) {
	for (unsigned int ii=0; ii<sz; ii++) {
	  out << " " << ii << " vfatType=" << ptr->vfatType[ii]
	      << " vfatCh=" << ptr->vfatCh[ii]
	      << " vfatStrip=" << ptr->vfatStrip[ii] << "\n";
	}
      }
    } // for theStripMap
  } // if (detailed)
}


void GEMELMap::GEMVFatMap::printLast(std::ostream &out, int printEOL) const
{
  out << "GEMVFatMap::printLast ";
  if (vfat_position.size()==0) {
    out << "empty";
    return;
  }
  out << "vfat_pos=" << vfat_position.back() << " z_dir=" << z_direction.back()
      << " iEta=" << iEta.back() << " iPhi=" << iPhi.back()
      << " depth=" << depth.back() << " vfatType=" << vfatType.back()
      << " vfatId=0x" << std::hex << vfatId.back() << std::dec
      << " amcId=" << amcId.back()
      << " gebId=" << gebId.back() << " sec=" << sec.back();
  if (printEOL) out << "\n";
}

int GEMELMap::GEMStripMap::areIdentical(const GEMELMap::GEMStripMap &mp, int printDiff) const
{
  int res = ( (vfatType == mp.vfatType) &&
	      (vfatCh == mp.vfatCh) &&
	      (vfatStrip == mp.vfatStrip) ) ? 1:0;
  if (!res && printDiff) {
    if (vfatType.size()!=mp.vfatType.size()) {
      std::cout << "areIdentical=false: sizes are different: "
		<< vfatType.size() << " " << mp.vfatType.size() << "\n";
    }
    else {
      for (unsigned int i=0; i<vfatType.size(); i++) {
	std::cout << " i=" << i << "  " << vfatType[i] << " " << vfatCh[i]
		  << " " << vfatStrip[i] << "  "
		  << mp.vfatType[i] << " " << mp.vfatCh[i] << " "
		  << mp.vfatStrip[i] << "\n";
      }
    }
  }
  return res;
}

void GEMELMap::GEMStripMap::printLast(std::ostream &out, int printEOL) const
{
  out << "GEMStripMap::printLast ";
  if (vfatType.size()==0) {
    out << "empty";
    return;
  }
  out << "vfatType=" << vfatType.back() << " vfatCh=" << vfatCh.back()
      << " vfatStrip=" << vfatStrip.back();
  if (printEOL) out << "\n";
}
