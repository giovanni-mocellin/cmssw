#ifndef noFileInPath_H
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "CondFormats/GEMObjects/interface/GEMELMap.h"
#include "CondFormats/GEMObjects/interface/GEMELMapHelper.h"
#else
#include "../interface/GEMELMap.h"
#include "../interface/GEMELMapHelper.h"
#endif

#include <map>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>


// -------------------------------------------------------

GEMELMapHelper::GEMELMapHelper() :
  chamberId_(_chamber_none),
  vfat2mapId_(), mapId2stripCh_()
{}

// -------------------------------------------------------

GEMELMapHelper::GEMELMapHelper(const GEMELMapHelper &h) :
  chamberId_(h.chamberId_),
  vfat2mapId_(h.vfat2mapId_), mapId2stripCh_(h.mapId2stripCh_)
{}

// -------------------------------------------------------

std::string GEMELMapHelper::vfatMapName(GEMELMapHelper::TVFatMaps v)
{
  std::string name="vfat_unknown";
  switch(v) {
  case _vfatMap_none: name="vfatMap_none"; break;
  case _vfatMap_conf0: name="vfatMap_conf0"; break;
  case _vfatMap_conf1: name="vfatMap_conf1"; break;
  case _vfatMap_conf2: name="vfatMap_conf2"; break;
  case _vfatMap_conf3: name="vfatMap_conf3"; break;
  case _vfatMap_last: name="vfatMap_last"; break;
  default: ;
  }
  return name;
}

// -------------------------------------------------------

std::string GEMELMapHelper::chamberIdName(GEMELMapHelper::TChamberId id)
{
  std::string name="chamber_unknown";
  switch(id) {
  case _chamber_none: name="chamber_none"; break;
  case _chamber_longVFat2: name="chamber_longVFat2"; break;
  case _chamber_shortVFat2: name="chamber_shortVFat2"; break;
  case _chamber_longVFat3bV2: name="chamber_longVFat3bV2"; break;
  case _chamber_shortVFat3bV2: name="chamber_shortVat3bV2"; break;
  case _chamber_last: name="chamber_last"; break;
  default: ;
  }
  return name;
}

// -------------------------------------------------------

int GEMELMapHelper::assign(const GEMELMap &elMap)
{
  chamberId_= _chamber_none;

  vfat2mapId_ = std::vector<int>(24,-1);
  for (unsigned int i=0; i<elMap.theVFatMap_.size(); i++) {
    const GEMELMap::GEMVFatMap * p = & elMap.theVFatMap_.at(i);
    for (unsigned int ii=0; ii<p->size(); ii++) {
      int vfatPos= p->vfat_position[ii];
      int vfatType= p->vfatType[ii];
      if (vfat2mapId_[vfatPos]==-1) vfat2mapId_[vfatPos]=vfatType;
      else if (vfat2mapId_[vfatPos]!=vfatType) {
	std::cout << "mixed GEMELMap\n";
	return 0;
      }
    }
  }
  // check the assignment
  if (std::find(vfat2mapId_.begin(),vfat2mapId_.end(),-1)!=vfat2mapId_.end()) {
    std::cout << "vfat2mapId_ contains -1\n";
    return 0;
  }


  int stripMapCount=0;
  for (unsigned int i=0; i < elMap.theStripMap_.size(); i++) {
    if (elMap.theStripMap_[i].size()==0) {
      std::cout << "elMap.theStripMap_[" << i << "] size is 0\n";
      return 0;
    }
    std::vector<int>::const_iterator it = std::max_element(elMap.theStripMap_[i].vfatType.begin(),elMap.theStripMap_[i].vfatType.end());
    const int test= *it;
    std::cout << "test=" << test << "\n";
    if (test>stripMapCount) stripMapCount=test;
    if (test==0) {
      std::cout << "stripMapType cannot be 0\n";
      return 0;
    }
  }

  for (int imap=0; imap<stripMapCount; imap++) {
    //mapId2stripCh_[imap+1];
    std::vector<int> str2ch = std::vector<int>(128,-1);
    for (unsigned int i=0; i<elMap.theStripMap_.size(); i++) {
      const GEMELMap::GEMStripMap *p = & elMap.theStripMap_.at(i);
      for (unsigned int ii=0; ii<p->size(); ii++) {
	if (p->vfatType[ii] == imap+1) {
	  //std::cout << "strip->ch " << p->vfatStrip[ii] << " " << p->vfatCh[ii] << "\n";
	  if (str2ch[p->vfatStrip[ii]]!=-1) {
	    std::cout << "duplicated map " << imap+1 << "\n";
	    return 0;
	  }
	  str2ch[p->vfatStrip[ii]] = p->vfatCh[ii];
	  //std::cout << "strip->ch " << p->vfatStrip[ii] << " " << p->vfatCh[ii] << "  " << str2ch[p->vfatStrip[ii]] << "\n";
	}
      }
    }
      std::cout << "imap=" << imap << "\n";
    // check validity
    if (std::find(str2ch.begin(),str2ch.end(),-1)!=str2ch.end()) {
      std::cout << "there was a strip without a channel\n";
      std::vector<int>::iterator it= std::find(str2ch.begin(),str2ch.end(),-1);
      std::cout << "  stripNo=" << int(it-str2ch.begin()) << "\n";
      for (unsigned int i=0; i<str2ch.size(); i++) {
	std::cout << " i = " << i << " val=" << str2ch[i] << "\n";
      }
      return 0;
    }
    // assign
    mapId2stripCh_[imap+1] = str2ch;
  }

  return 1;
}

// -------------------------------------------------------

int GEMELMapHelper::addStrip2ChanInfo(GEMELMap &elMap)
{

  // create a self-copy to be able to modify
  GEMELMapHelper h(*this);

  // if the map is not empty, adjust strip2channel information
  if ((elMap.theVFatMap_.size()!=0) &&
      (elMap.theStripMap_.size()!=0)) {

    GEMELMapHelper hInp;
    if (!hInp.assign(elMap)) {
      std::cout << "failed to make a copy of the input elMap\n";
      return 0;
    }

    // match strip2chan
    std::vector<int> match; // indices, not keys!
    int changeNeeded=0;

    for (unsigned int iOur=0; iOur<h.mapId2stripCh().size(); iOur++) {
      if ((*h.strip2channel_atKeyIdx_ptr(iOur)) == (*hInp.strip2channel_atKeyIdx_ptr(iOur))) {
	match.push_back(iOur);
      }
      else {
	changeNeeded=1;
	int idx=-1;
	for (unsigned int iInp=0; iInp<hInp.mapId2stripCh().size(); iInp++) {
	  if ((*h.strip2channel_atKeyIdx_ptr(iOur)) == (*hInp.strip2channel_atKeyIdx_ptr(iInp))) {
	    idx=int(iInp);
	    break;
	  }
	}
	if (idx<0) {
	  std::cout << "Failed to find the corresponding strip2channel map\n";
	  return 0;
	}
	match.push_back(idx);
      }
    }

    if (changeNeeded) {
      for (unsigned int i=0; i<match.size(); i++) {
	if (int(i)==match[i]) continue;
	int map1= h.strip2channelKey_atIdx(i);
	int map2= h.strip2channelKey_atIdx(match[i]);
	if (map1>=map2) continue;
	h.swapStrip2ChanMap(map1,map2);
      }
    }
  }

  // if the strip2channel map is empty, just create it
  if (elMap.theStripMap_.size()==0) {
    GEMELMap::GEMStripMap mp;
    createStrip2ChanInfo(mp);
    elMap.theStripMap_.push_back(mp);
  }

  // update vfat map with the strip2channel map info
  for (unsigned int ivfm=0; ivfm<elMap.theVFatMap_.size(); ivfm++) {
    GEMELMap::GEMVFatMap *p= & elMap.theVFatMap_[ivfm];
    //std::cout << "p->size=" << p->size() << ", p->vfatType.size()=" << p->vfatType.size() << "\n";
    p->vfatType.clear();
    p->vfatType = std::vector<int>(p->size(),-1);
    for (unsigned int i=0; i<p->size(); i++) {
      p->vfatType[i] = h.vfat2mapId_[p->vfat_position.at(i)];
    }
  }

  return 1;
}

// -------------------------------------------------------

int GEMELMapHelper::createStrip2ChanInfo(GEMELMap::GEMStripMap &mp) {
  mp.vfatType.clear();
  mp.vfatCh.clear();
  mp.vfatStrip.clear();
  mp.vfatType.reserve(4*128);
  mp.vfatCh.reserve(4*128);
  mp.vfatStrip.reserve(4*128);
  for (auto it : mapId2stripCh_) {
    for (unsigned int i=0; i<it.second.size(); i++) {
      mp.vfatType.push_back( it.first );
      mp.vfatCh  .push_back( it.second.at(i) );
      mp.vfatStrip.push_back( int(i) );
    }
  }
  return 1;
}

// -------------------------------------------------------

int GEMELMapHelper::load(GEMELMapHelper::TChamberId set_chamberId, int autoAdjustVFat2ChNum)
{
  if (chamberId_ == set_chamberId) return 1;
  chamberId_ = set_chamberId;
  vfat2mapId_.clear();
  mapId2stripCh_.clear();
  if (set_chamberId == _chamber_none) return 1;

  // Prepare the file name
  std::string fname;
  switch(chamberId_) {
  case _chamber_none: return 1;
  case _chamber_longVFat2: fname="vfat2_longChannelMap.txt"; break;
  case _chamber_shortVFat2: fname="vfat2_shortChannelMap.txt"; break;
  case _chamber_longVFat3bV2: fname="vfat3bV2_long.csv"; break;
  case _chamber_shortVFat3bV2: fname="vfat3bV2_short.csv"; break;
  default:
    std::cout << "GEMELMapHelper::load is not ready for this chamber type" << std::endl;
    return 0;
  }

  fname.insert(0, "CondFormats/GEMObjects/data/");
  edm::FileInPath fp(fname);
  fname= fp.fullPath();

  int res=0;
  switch(chamberId_) {
  case _chamber_longVFat2:
  case _chamber_shortVFat2:
    res= loadVFat2(fname,chamberId_);
    if (res && autoAdjustVFat2ChNum) modifyChanNum(-1);
    break;
  case _chamber_longVFat3bV2:
  case _chamber_shortVFat3bV2: res= loadVFat3bV2(fname,chamberId_); break;
  default:
    std::cout << "GEMELMapHelper::load code error" << std::endl;
  }

  return res;
}

// -------------------------------------------------------

int GEMELMapHelper::loadVFat2(const std::string &fname, GEMELMapHelper::TChamberId set_chamberId)
{
  std::cout << "GEMELMapHelper::loadVFat2\n";

  // load the file
  std::ifstream fin(fname.c_str());
  if (!fin.is_open()) {
    std::cout << "Failed to open the file " << fname << std::endl;
    return 0;
  }

  std::string line;
  getline(fin,line);
  //std::cout << "header =" << line << "\n";

  int iVfat_old=-1, iVfat,iStr,iCh;
  int iLine=0;
  std::vector<int> str2ch;
  while (!fin.eof()) {
    iLine++;
    getline(fin,line);
    std::stringstream ss(line);
    ss >> iVfat >> iStr >> iCh;
    //std::cout << "iLine=" << iLine << ", iVfat=" << iVfat << " iStrip=" << iStr
    //<< " iChan=" << iCh << "   " << "Vfat_old=" << iVfat_old << "\n";
    //iCh += iChNoCorr;
    if (((iVfat!=iVfat_old) || fin.eof()) && str2ch.size()) {
      int idx=-1;
      int i=0;
      for (auto it : mapId2stripCh_) {
	if (it.second == str2ch) {
	  idx= i + 1; // 0th map is 'none'
	  break;
	}
	i++;
      }
      //std::cout << "stripCh map idx=" << idx << ", vfat2mapId.size=" << vfat2mapId_.size() << ", mapId2stripCh_.size=" << mapId2stripCh_.size() << "\n";
      // this is a new map
      if (idx==-1) {
	idx= mapId2stripCh_.size() + 1;
	mapId2stripCh_[idx] = str2ch;
      }
      if ((iVfat>=0) && (iVfat_old!=int(vfat2mapId_.size()))) {
	std::cout << "file format assumption error (vfat sequence)" << std::endl;
	return 0;
      }
      vfat2mapId_.push_back(idx);
      str2ch.clear();
    }

    if ((iStr!=int(str2ch.size())) && !fin.eof()) {
      std::cout << "file format assumption error (strip-channel)" << std::endl;
      return 0;
    }
    str2ch.push_back(iCh);
    iVfat_old = iVfat;
  }
  fin.close();
  std::cout << "file <" << fname << "> successfully loaded\n";
  return 1;
}

// -------------------------------------------------------

int GEMELMapHelper::loadVFat3bV2(const std::string &fname, GEMELMapHelper::TChamberId set_chamberId)
{
  std::cout << "GEMELMapHelper::loadVFat3bV2\n";

  // load the file
  std::ifstream fin(fname.c_str());
  if (!fin.is_open()) {
    std::cout << "Failed to open the file " << fname << std::endl;
    return 0;
  }

  // Skip 3 lines and get hybrid pos definitions
  std::string line;
  for (int i=0; i<4; i++) {
    getline(fin,line);
    //std::cout << "header =" << line << "\n";
  }

  // read hybrid positions
  {
    vfat2mapId_ = std::vector<int>(24,-1);
    //for (unsigned int i=0; i<vfat2mapId_.size(); i++) { std::cout << "i=" << i << " val=" << vfat2mapId_[i] << "\n"; }
    size_t p=0;
    for (int iMap=0; iMap<4; iMap++) {
      p= line.find("Positions",p+1);
      std::stringstream ss(line.c_str()+p+strlen("Positions")+1);
      //std::cout << "read from stringstream " << ss.str() << "\n";
      char sep=' '; // separator
      int vfatPos=0;
      while (sep!='"') {
	ss >> vfatPos >> sep;
	vfat2mapId_[vfatPos] = iMap + int(_vfatMap_conf0) + GEMELMap::vfatTypeV3_;
	//std::cout << "vfatPos=" << vfatPos << " map=" << vfat2mapId_[vfatPos] << "\n";
      }
    }
  }

  // check for missing values
  for (unsigned int i=0; i<vfat2mapId_.size(); i++) {
    if (vfat2mapId_[i]==-1) {
      std::cout << "detected bad entry in vfat2mapId_\n";
      return 0;
    }
  }

  // skip last header line
  getline(fin,line);
  //std::cout << "last line <" << line << ">\n";

  // load the strip-to-channel map
  for (int mp= int(_vfatMap_conf0); mp<int(_vfatMap_last); mp++) {
    mapId2stripCh_[mp + GEMELMap::vfatTypeV3_] = std::vector<int>(128,-1);
  }


  int iPin=0, iCh=0, iStrip=0;
  char sep1,sep2,sep3;
  while (!fin.eof() && getline(fin,line)) {
    if (line.size()<7) continue;
    std::stringstream ss(line.c_str());
    ss >> iPin >> sep1 >> sep2 >> sep3 >> iCh;
    for (int iMap=0; iMap<4; iMap++) {
      ss >> sep1 >> iStrip;
      mapId2stripCh_[iMap + int(_vfatMap_conf0) + GEMELMap::vfatTypeV3_][iStrip] = iCh;
    }
  }

  // check for missing entries
  for (auto it : mapId2stripCh_) {
    for (unsigned int i=0; i<it.second.size(); i++) {
      if (it.second[i]==-1) {
	std::cout << "detected undefined strip2channel entry\n";
	return 0;
      }
    }
  }

  fin.close();
  std::cout << "file <" << fname << "> successfully loaded\n";
  return 1;
}

// -------------------------------------------------------

void GEMELMapHelper::print(std::ostream &out, int asList) const
{
  out << "GEMELMapHelper::print\n";
  out << "chamberId=" << GEMELMapHelper::chamberIdName(chamberId_) << "\n";
  out << "vfatPos --> mapId [" << vfat2mapId_.size() << "]:\n";
  for (unsigned int i=0; i<vfat2mapId_.size(); i++) {
    out << "vfatPos=" << i << " mapId=" << vfat2mapId_[i] << " "
	<< GEMELMapHelper::vfatMapName(TVFatMaps(vfat2mapId_[i])) << "\n";
  }
  out << "mapId2stripCh[" << mapId2stripCh_.size() << "]:\n";
  for (auto it : mapId2stripCh_) {
    out << "mapId=" << it.first << " "
	<< GEMELMapHelper::vfatMapName(TVFatMaps(it.first)) << ":\n";
    for (unsigned int i=0; i<it.second.size(); i++) {
      out << " " << it.second.at(i);
    }
    out << "\n";
  }
  if (asList) {
    out << "as list:\n";
    out << "vfatPos Strip Channel\n";
    char buf1[10], buf2[10], buf3[10];
    for (unsigned int i=0; i<vfat2mapId_.size(); i++) {
      std::map<int,std::vector<int> >::const_iterator it= mapId2stripCh_.find(i);
      if (it==mapId2stripCh_.end()) {
	std::cout << "code error\n";
	break;
      }
      const std::vector<int> * ptr= & it->second;
      for (unsigned int ii=0; ii<ptr->size(); ii++) {
	snprintf(buf1,10,"%2d",i);
	snprintf(buf2,10,"%3d",ii);
	snprintf(buf3,10,"%3d",ptr->at(ii));
	out << " " << buf1 << " " << buf2 << " " << buf3 << "\n";
      }
    }
  }
}

// -------------------------------------------------------

int GEMELMapHelper::strip2channelKey_atIdx(int idx) const {
  auto it = mapId2stripCh_.begin();
  while ((idx>0) && (it!=mapId2stripCh_.end())) { it++; idx--; }
  if (it == mapId2stripCh_.end()) return -1;
  else return it->first;
}

// -------------------------------------------------------

const std::vector<int>* GEMELMapHelper::strip2channel_ptr(int iMap) const
{
  auto it = mapId2stripCh_.find(iMap);
  if (it == mapId2stripCh_.end()) return NULL;
  else return & it->second;
}

// -------------------------------------------------------

const std::vector<int>* GEMELMapHelper::strip2channel_atKeyIdx_ptr(int iKey) const
{
  auto it = mapId2stripCh_.begin();
  while ((iKey>0) && (it!=mapId2stripCh_.end())) { it++; iKey--; }
  if (it == mapId2stripCh_.end()) return NULL;
  else return & it->second;
}

// -------------------------------------------------------

void GEMELMapHelper::modifyChanNum(int addVal)
{
  std::cout << "GEMELMapHelper::modifyChanNum(addVal=" << addVal << ")\n";
  for (std::map<int,std::vector<int> >::iterator it = mapId2stripCh_.begin(); it!=mapId2stripCh_.end(); it++) {
    std::vector<int> *p = & it->second;
    for (unsigned int i=0; i<p->size(); i++) {
      (*p)[i] += addVal;
    }
  }
}

// -------------------------------------------------------

void GEMELMapHelper::reverseStripOrdering()
{
  int first=0;
  for (std::map<int,std::vector<int> >::iterator it = mapId2stripCh_.begin(); it!=mapId2stripCh_.end(); it++) {
    std::vector<int> *p = & it->second;
    for (unsigned int i=0; i<p->size()/2; i++) {
      int t=(*p)[i];
      int idx= p->size()-1-i;
      if (first) std::cout << "reverse " << i << " and " << idx << "\n";
      (*p)[i] = (*p)[idx];
      (*p)[idx] = t;
    }
    first=0;
  }
}

// -------------------------------------------------------

int GEMELMapHelper::swapStrip2ChanMap(int iMap1, int iMap2)
{
  std::map<int,std::vector<int> >::iterator it1 = mapId2stripCh_.find(iMap1);
  std::map<int,std::vector<int> >::iterator it2 = mapId2stripCh_.find(iMap2);
  if ((it1==mapId2stripCh_.end()) ||
      (it2==mapId2stripCh_.end())) {
    std::cout << "swapStrip2ChanMap: at least one iMap value is not correct: "
	      << iMap1 << " and " << iMap2 << "\n";
    return 0;
  }
  it1->second.swap(it2->second);

  for (unsigned int i=0; i<vfat2mapId_.size(); i++) {
    if (vfat2mapId_[i] == iMap2) vfat2mapId_[i]=iMap1;
    else if (vfat2mapId_[i] == iMap1) vfat2mapId_[i]=iMap2;
  }
  return 1;
}

// -------------------------------------------------------
// -------------------------------------------------------
// -------------------------------------------------------

void select_target_lines(std::vector<std::string> &lines, std::vector<std::string> &v, const std::string &target)
{
  for (std::vector<std::string>::iterator it=lines.begin(); it<lines.end(); it++) {
    if (it->find(target)!=std::string::npos) {
      v.push_back(*it);
      lines.erase(it);
      it--;
    }
  }
}

// -------------------------------------------------------

int read_GEMELMap(std::ifstream &fin, GEMELMap &elMap, int onlyGeb) {
  if (!fin.is_open()) {
    std::cout << "read_GEMELMap fin is not open\n";
    return 0;
  }

  std::string line;

  // find the "header" line
  do {
    getline(fin,line);
  } while (!fin.eof() && (line.find("GEMELMap:")==std::string::npos));

  getline(fin,line);
  std::cout << "header=" << line << "\n";
  size_t p= line.find('[');
  int cntVFATMap=atoi(line.c_str()+p+1);
  std::cout << "cntVFATMap=" << cntVFATMap << "\n";

  GEMELMap::GEMVFatMap vfatMap;
  for (int icnt=0; icnt<cntVFATMap; icnt++) {
    getline(fin,line);
    //std::cout << "decode " << line << "\n";
    p= line.find('=',line.find('=')+1);
    vfatMap.VFATmapTypeId= atoi(line.c_str()+p+1);
    p= line.find('=', p+1);
    int sz= atoi(line.c_str()+p+1);
    //std::cout << " VFATmapTypeId=" << vfatMap.VFATmapTypeId << ", size=" << sz << "\n";
    // collect into buffer
    std::vector<std::string> lines;
    for (int i=0; i<sz; i++) {
      getline(fin,line);
      lines.push_back(line);
    }

    // sort by sector
    //for (int iSect=1; (iSect<31) && lines.size(); iSect++) {
    for (int iGeb=0; (iGeb<10) && lines.size(); iGeb++) {
      if ((onlyGeb!=-1) && (iGeb!=onlyGeb)) continue;
      std::vector<std::string> v;
      std::stringstream sout;
      sout << "gebId=" << iGeb << " ";
      std::string target= sout.str();
      select_target_lines(lines,v,target);
      if (v.size()==0) continue;
      //std::cout << "gebId=" << iGeb << " has " << v.size() << " lines\n";
      //for (unsigned int i=0; i<v.size(); i++) {
      //std::cout << "i=" << i << " | " << v[i] << "\n";
      //}
      if (v.size()!=24) {
	std::cout << "each GEB should have 24 vfats\n";
	return 0;
      }
      std::vector<std::string> vfatInfo;
      for (unsigned int vfatP=0; vfatP<24; vfatP++) {
	sout.str(std::string());
	sout << "vfat_pos=" << vfatP << " ";
	vfatInfo.clear();
	select_target_lines(v,vfatInfo,sout.str());
	if (vfatInfo.size()==0) {
	  std::cout << "vfatInfo.size=0 for target=" << sout.str() << "\n";
	  return 0;
	}
	line= vfatInfo[0];
	//std::cout << "line=" << line << "\n";
	size_t p= line.find('=');
	vfatMap.vfat_position.push_back(atoi(line.c_str()+p+1));
	p= line.find('=',p+1);
	vfatMap.z_direction.push_back(atoi(line.c_str()+p+1));
	p= line.find('=',p+1);
	vfatMap.iEta.push_back(atoi(line.c_str()+p+1));
	p= line.find('=',p+1);
	vfatMap.iPhi.push_back(atoi(line.c_str()+p+1));
	p= line.find('=',p+1);
	vfatMap.depth.push_back(atoi(line.c_str()+p+1));
	p= line.find('=',p+1);
	vfatMap.vfatType.push_back(atoi(line.c_str()+p+1));
	p= line.find('=',p+1);
	int vfatId=0;
	{
	  std::stringstream sinp(line.c_str()+p+1);
	  sinp >> std::hex >> vfatId;
	}
	vfatMap.vfatId.push_back(vfatId);
	p= line.find('=',p+1);
	vfatMap.amcId.push_back(atoi(line.c_str()+p+1));
	p= line.find('=',p+1);
	vfatMap.gebId.push_back(atoi(line.c_str()+p+1));
	p= line.find('=',p+1);
	vfatMap.sec.push_back(atoi(line.c_str()+p+1));
	//vfatMap.printLast();
      }
    }
    elMap.theVFatMap_.push_back(vfatMap);
  } // for cntVFATMap

  getline(fin,line);
  if (line.find("theStripMap[")==std::string::npos) {
    std::cout << "'theStripMap[' was not found where expected\n";
    return 0;
  }
  p= line.find('[');
  int cntStripMap= atoi(line.c_str()+p+1);
  std::cout << "cntStripMap=" << cntStripMap << "\n";

  GEMELMap::GEMStripMap stripMap;
  for (int icnt=0; icnt<cntStripMap; icnt++) {
    getline(fin,line);
    //std::cout << "decode " << line << "\n";
    p= line.find('=',line.find('=')+1);
    int sz= atoi(line.c_str()+p+1);
    //std::cout << " stripMap size=" << sz << "\n";

    // collect into buffer
    std::vector<std::string> lines;
    for (int i=0; i<sz; i++) {
      getline(fin,line);
      lines.push_back(line);
    }

    // read the buffer
    for (int vfatType=1; (vfatType<10) && lines.size(); vfatType++) {
      std::vector<std::string> v;
      std::stringstream sout;
      sout << "vfatType=" << vfatType << " ";
      std::string target= sout.str();
      select_target_lines(lines,v,target);
      if (v.size()==0) continue;
      //std::cout << "vfatType=" << vfatType << " has " << v.size() << " lines\n";
      //for (unsigned int i=0; i<v.size(); i++) {
      //std::cout << "i=" << i << " | " << v[i] << "\n";
      //}
      if (v.size()!=128) {
	std::cout << "each vfatType should have 128 channels\n";
	return 0;
      }
      for (unsigned int i=0; i<v.size(); i++) {
	line= v[i];
	//std::cout << "line=" << line << "\n";
	size_t p= line.find('=');
	stripMap.vfatType.push_back(atoi(line.c_str()+p+1));
	p= line.find('=',p+1);
	stripMap.vfatCh.push_back(atoi(line.c_str()+p+1));
	p= line.find('=',p+1);
	stripMap.vfatStrip.push_back(atoi(line.c_str()+p+1));
	//stripMap.printLast();
      }
    }
    elMap.theStripMap_.push_back(stripMap);
  } // for cntStripMap

  return 1;
}

// -------------------------------------------------------

int gemelmap_areIdentical(const GEMELMap &em1, const GEMELMap &em2, int manipulate=0,
			  int partOnly=-1)
{
  if ((em1.theVFatMap_.size()!=em2.theVFatMap_.size()) ||
      (em1.theStripMap_.size()!=em2.theStripMap_.size())) return 0;

  int equal=1;
  if ((partOnly==-1) || (partOnly==1)) {
    for (unsigned int i=0; (i<em1.theVFatMap_.size()) && equal; i++) {
      const GEMELMap::GEMVFatMap *m1 = & em1.theVFatMap_.at(i);
      const GEMELMap::GEMVFatMap *m2 = & em2.theVFatMap_.at(i);
      if ( (m1->VFATmapTypeId != m2->VFATmapTypeId) ||
	   (m1->vfat_position != m2->vfat_position) ||
	   (m1->z_direction != m2->z_direction) ||
	   (m1->iEta != m2->iEta) ||
	   (m1->iPhi != m2->iPhi) ||
	   (m1->depth != m2->depth) ||
	   (m1->vfatType != m2->vfatType) ||
	   (m1->vfatId != m2->vfatId) ||
	   (m1->amcId != m2->amcId) ||
	   (m1->gebId != m2->gebId) ||
	   (m1->sec != m2->sec) ) equal=0;
    }
  }

  if ((partOnly==-1) || (partOnly==2)) {
    int equal2=1;
    for (unsigned int i=0; (i<em1.theStripMap_.size()) && equal2; i++) {
      const GEMELMap::GEMStripMap *m1 = & em1.theStripMap_.at(i);
      const GEMELMap::GEMStripMap *m2 = & em2.theStripMap_.at(i);
      if ( (m1->vfatType != m2->vfatType) ||
	   (m1->vfatCh != m2->vfatCh) ||
	   (m1->vfatStrip != m2->vfatStrip) ) equal2=0;
    }

    if (!equal2 && manipulate) {
      equal2=1;
      for (unsigned int i=0; (i<em1.theStripMap_.size()) && equal2; i++) {
	const GEMELMap::GEMStripMap *m1 = & em1.theStripMap_.at(i);
	const GEMELMap::GEMStripMap *m2 = & em2.theStripMap_.at(i);
	for (unsigned int ii=0; (ii<m1->size()) && equal2; ii++) {
	  int ok=0;
	  for (unsigned int jj=0; jj<m2->size(); jj++) {
	    if ( (m1->vfatType[ii] == m2->vfatType[jj]) &&
		 (m1->vfatCh[ii] == m2->vfatCh[jj]) ) {
	      if (m1->vfatStrip[ii]==m2->vfatStrip[jj]) ok=1;
	      break;
	    }
	  }
	  if (!ok) equal2=0;
	}
      }
    }

    if (!equal2) equal=0;
  }

  return equal;
}

// -------------------------------------------------------

