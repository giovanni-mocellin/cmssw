#ifndef CondFormats_GEMObjects_GEMELMapHelper_h
#define CondFormats_GEMObjects_GEMELMapHelper_h

/** \class GEMELMapHelper
 *
 *  An auxiliary class to read GEMELMap vfatpos->strip2chan info from a file
 *
 *  \author A. Juodagalvis - Vilnius University
 *  Oct 2018
 */


#ifndef noFileInPath_H
#include "CondFormats/GEMObjects/interface/GEMELMap.h"
#endif

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>


// -------------------------------------------------------

class GEMELMapHelper {
public:
  typedef enum { _vfatMap_none=0, _vfatMap_conf0, _vfatMap_conf1,
		 _vfatMap_conf2, _vfatMap_conf3,
		 _vfatMap_last
  } TVFatMaps;
  typedef enum { _chamber_none=0,
		 _chamber_longVFat2, _chamber_shortVFat2,
		 _chamber_longVFat3bV2, _chamber_shortVFat3bV2,
		 _chamber_last
  } TChamberId;


  GEMELMapHelper();
  GEMELMapHelper(const GEMELMapHelper &h);
  int assign(const GEMELMap &mp);
  int addStrip2ChanInfo(GEMELMap &mp);
  int createStrip2ChanInfo(GEMELMap::GEMStripMap &mp);
  int load(TChamberId, int autoAdjustVFat2ChNum=1); // load definitions

protected:
  int loadVFat2(const std::string &fname, TChamberId);
  int loadVFat3bV2(const std::string &fname, TChamberId);

public:
  void print(std::ostream &out = std::cout, int asList=0) const;

  // access functions
  TChamberId chamberId() const { return chamberId_; }
  const std::vector<int>& vfat2mapId() const { return vfat2mapId_; }
  const std::map<int,std::vector<int> > & mapId2stripCh() const
  { return mapId2stripCh_; }

  int strip2channelKey_atIdx(int idx) const;
  const std::vector<int>* strip2channel_ptr(int iMap) const;
  const std::vector<int>* strip2channel_atKeyIdx_ptr(int iKey) const;

  void modifyChanNum(int addVal);
  void reverseStripOrdering();
  int swapStrip2ChanMap(int iMap1, int iMap2);

public:
  static std::string vfatMapName(TVFatMaps);
  static std::string chamberIdName(TChamberId);

protected:
  TChamberId chamberId_;
  std::vector<int> vfat2mapId_;   // vfat position to mapId
  std::map<int,std::vector<int> > mapId2stripCh_;
};

// -------------------------------------------------------

#ifndef noFileInPath_H
int read_GEMELMap(std::ifstream &fin, GEMELMap &elMap, int onlyGebNum=-1);
int gemelmap_areIdentical(const GEMELMap &elMap1, const GEMELMap &elMap2, int manipulate=0);
#endif

// -------------------------------------------------------

#endif // GEMELMapHelper_h
