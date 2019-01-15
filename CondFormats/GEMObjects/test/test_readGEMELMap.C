#include <string>
#include <iostream>

// remove dependency on FWCore ... FileInPath
#define FWCore_ParameterSet_FileInPath_h
#define test_GEMELMapHelper_C
#include "noFileInPath.h"
#include "../interface/GEMELMap.h"
#include "../src/GEMELMap.cc"
#include "../interface/GEMELMapHelper.h"
#include "../src/GEMELMapHelper.cc"

int test_mapId2stripCh_matching(const GEMELMapHelper &h1, std::string name1,
				const GEMELMapHelper &h2, std::string name2,
				int detailedPrint=0);
int test_vfat2mapId_matching(const GEMELMapHelper &h1, std::string name1,
			     const GEMELMapHelper &h2, std::string name2);

void testLongVsShort();
void testVFat2VsVFat3b2v(int theCase);
void testDuplicate(const GEMELMap &elMap);

// ---------------------------------------------------------

void test_readGEMELMap(int chId_user=0, int onlyGeb=-1)
{

  //testLongVsShort(); return;
  //testVFat2VsVFat3b2v(chId_user); return;

  if (onlyGeb>=0) {
    if ((onlyGeb/2)%2==0) chId_user=0; else chId_user=1;
    if (onlyGeb<2) chId_user=1;
  }

  std::ifstream fin("dumpELMap-from-sqlite.out");
  GEMELMap elMap;
  int res=read_GEMELMap(fin,elMap,onlyGeb);
  fin.close();

  //elMap.print(std::cout,1);
  testDuplicate(elMap); return;

  GEMELMapHelper hSqlite;
  if (!hSqlite.assign(elMap)) {
    std::cout << "conversion failed\n";
    return;
  }
  hSqlite.print();
  //return;

  GEMELMapHelper::TChamberId chId=(chId_user==0) ? GEMELMapHelper::_chamber_longVFat2 : GEMELMapHelper::_chamber_shortVFat2;
  if (chId_user==2) chId= GEMELMapHelper::_chamber_longVFat3bV2;
  else if (chId_user==3) chId= GEMELMapHelper::_chamber_shortVFat3bV2;
  std::cout << "chamberId=" << GEMELMapHelper::chamberIdName(chId) << "\n";

  GEMELMapHelper h;
  int autoAdjustVFat2ChNum=1;
  if (!h.load(chId,autoAdjustVFat2ChNum)) {
    std::cout << "failed\n";
    return;
  }
  h.print();

  // hypotheses
  //h.modifyChanNum(-1);  // taken care by autoAdjustVFat2ChNum
  //h.swapStrip2ChanMap(1,3); // obsolete w.r.t reverseStripOrdering
  //h.swapStrip2ChanMap(2,4); // obsolete w.r.t reverseStripOrdering

  // take care of incorrect definition order
  h.reverseStripOrdering();
  // different order of mapId definitions as compared to Hyunyong
  h.swapStrip2ChanMap(1,2);
  h.swapStrip2ChanMap(3,4);

  std::string label2="hExcel_" + GEMELMapHelper::chamberIdName(chId);
  int detailedPrint=0;
  test_mapId2stripCh_matching(hSqlite,"hSqlite",h,label2,detailedPrint);
  //test_vfat2mapId_matching(hSqlite,"hSqlite",h,label2);
}

// ---------------------------------------------------------
// ---------------------------------------------------------

int test_mapId2stripCh_matching(const GEMELMapHelper &h1, std::string name1,
				const GEMELMapHelper &h2, std::string name2,
				int detailedPrint)
{
  std::cout << "test_mapId2stripCh_matching: " << name1 << " vs " << name2 << "\n";

  std::stringstream sout;
  for (unsigned int imap1=0; imap1<h1.mapId2stripCh().size(); imap1++) {
    const std::vector<int> *s2c_1 = h1.strip2channel_atKeyIdx_ptr(imap1);
    if (!s2c_1) { std::cout << "s2c_1 is null\n"; return 0; }
    for (unsigned int imap2=0; imap2<h2.mapId2stripCh().size(); imap2++) {
      const std::vector<int> *s2c_2 = h2.strip2channel_atKeyIdx_ptr(imap2);
      if (!s2c_2) { std::cout << "s2c_2 is null\n"; return 0; }
      std::cout << "compare imaps " << imap1 << " to " << imap2 << "  : "
		<< (((*s2c_1)==(*s2c_2)) ? 1 : 0) << "\n";
      if (detailedPrint) {
	double dist=0;
	for (unsigned int ii=0; ii<s2c_1->size(); ii++) {
	  if (detailedPrint==1) {
	    std::cout << " ii " << ii << "  " << s2c_1->at(ii)
		      << " " << s2c_2->at(ii) << "\n";
	  }
	  dist+= (s2c_1->at(ii) - s2c_2->at(ii))*(s2c_1->at(ii) - s2c_2->at(ii));
	}
	sout << " dist imaps " << imap1 << " to " << imap2 << "  = "
	     << sqrt(dist) << "\n";
      }
    }
  }

  if (sout.str().size()) std::cout << "\n" << sout.str() << "\n";

  return 1;
}

// ---------------------------------------------------------

int test_vfat2mapId_matching(const GEMELMapHelper &h1, std::string name1,
			     const GEMELMapHelper &h2, std::string name2)
{
  std::cout << "test_vfat2mapId_matching: " << name1 << " vs " << name2 << "\n";

  if (h1.vfat2mapId() == h2.vfat2mapId()) {
    std::cout << "maps are equal\n";
  }
  else {
    std::cout << "maps are different\n";
    for (unsigned int i=0; i<h1.vfat2mapId().size(); i++) {
      std::cout << "i=" << i << "  " << h1.vfat2mapId()[i] << " "
		<< h2.vfat2mapId()[i] << "   "
		<< ((h1.vfat2mapId()[i]==h2.vfat2mapId()[i]) ? "ok" : "diff")
		<< "\n";
    }
  }

  return 1;
}

// ---------------------------------------------------------

void testLongVsShort()
{

  GEMELMapHelper hShort, hLong;
  int autoAdjustVFat2ChNum=1;
  if (!hShort.load(GEMELMapHelper::_chamber_shortVFat2,autoAdjustVFat2ChNum) ||
      !hLong.load(GEMELMapHelper::_chamber_longVFat2,autoAdjustVFat2ChNum)) {
    std::cout << "failed\n";
    return;
  }
  hShort.print();
  hLong.print();

  test_mapId2stripCh_matching(hLong,"vfat2_long",hShort,"vfat2_short");

}

// ---------------------------------------------------------

void testVFat2VsVFat3b2v(int theCase)
{
  GEMELMapHelper::TChamberId chId2=(theCase%2==0) ? GEMELMapHelper::_chamber_longVFat2 : GEMELMapHelper::_chamber_shortVFat2;
  GEMELMapHelper::TChamberId chId3=(theCase/2==0) ? GEMELMapHelper::_chamber_longVFat3bV2 : GEMELMapHelper::_chamber_shortVFat3bV2;

  GEMELMapHelper hVFat2, hVFat3bV2;
  int autoAdjustVFat2ChNum=1;

  if (!hVFat2.load(chId2,autoAdjustVFat2ChNum) ||
      !hVFat3bV2.load(chId3,autoAdjustVFat2ChNum)) {
    std::cout << "failed\n";
    return;
  }
  hVFat2.print();
  hVFat3bV2.print();

  int detailedPrint=-1;
  test_mapId2stripCh_matching(hVFat2,GEMELMapHelper::chamberIdName(chId2),
			      hVFat3bV2,GEMELMapHelper::chamberIdName(chId3),
			      detailedPrint);

}

// ---------------------------------------------------------

void testDuplicate(const GEMELMap &elMap)
{

  GEMELMapHelper hSqlite;
  if (!hSqlite.assign(elMap)) {
    std::cout << "conversion failed\n";
    return;
  }
  hSqlite.print();

  GEMELMap elMapCopy(elMap);
  for (unsigned int i=0; i<elMapCopy.theVFatMap_.size(); i++) {
    elMapCopy.theVFatMap_[i].vfatType.clear();
  }
  //std::cout << "check " << gemelmap_areIdentical(elMap,elMapCopy) << "\n";

  elMapCopy.theStripMap_.clear();

  hSqlite.addStrip2ChanInfo(elMapCopy);

  if (gemelmap_areIdentical(elMap,elMapCopy,1,-1)) {
    std::cout << "map and its re-engineered copy are equal\n";
  }
  else {
    std::cout << "map and its re-engineered copy are different\n";
    elMap.print(std::cout,1);
    elMapCopy.print(std::cout,1);
  }
}
