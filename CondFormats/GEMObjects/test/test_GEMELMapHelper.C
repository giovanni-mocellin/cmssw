#include <string>
#include <iostream>

// remove dependency on FWCore ... FileInPath
#define FWCore_ParameterSet_FileInPath_h
#define test_GEMELMapHelper_C
#include "noFileInPath.h"
#include "../src/GEMELMapHelper.cc"

void test_statics();

// ---------------------------------------------------------

void test_GEMELMapHelper(int chId_user=0)
{
  if (0) { test_statics(); return; }

  GEMELMapHelper::TChamberId chId=(chId_user==0) ? GEMELMapHelper::_chamber_longVFat2 : GEMELMapHelper::_chamber_shortVFat2;
  if (chId_user==2) chId= GEMELMapHelper::_chamber_longVFat3bV2;
  else if (chId_user==3) chId= GEMELMapHelper::_chamber_shortVFat3bV2;
  std::cout << "chamberId=" << GEMELMapHelper::chamberIdName(chId) << "\n";

  GEMELMapHelper h;
  if (!h.load(chId)) {
    std::cout << "failed\n";
    return;
  }
  h.print();
}

// ---------------------------------------------------------

void test_statics()
{
  for (int i=int(GEMELMapHelper::_vfatMap_none); i<int(GEMELMapHelper::_vfatMap_last); i++) {
    GEMELMapHelper::TVFatMaps m= GEMELMapHelper::TVFatMaps(i);
    std::cout << "vfatMap i=" << i << " " << GEMELMapHelper::vfatMapName(m) << "\n";
  }
  for (int i=int(GEMELMapHelper::_chamber_none); i<int(GEMELMapHelper::_chamber_last); i++) {
    GEMELMapHelper::TChamberId id= GEMELMapHelper::TChamberId(i);
    std::cout << "chamber i=" << i << " " << GEMELMapHelper::chamberIdName(id) << "\n";
  }
}

// ---------------------------------------------------------
