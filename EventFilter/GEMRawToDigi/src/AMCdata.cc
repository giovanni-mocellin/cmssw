#include "EventFilter/GEMRawToDigi/interface/AMCdata.h"
#include <string>
#include <sstream>

std::string gem::AMCdata::getAMCheader1_str() const {
  AMCheader1 h1(amch1_);
  std::stringstream ss;
  ss << "AMCheader1: reserved=" << (uint32_t)(h1.reserved) << " AMCnum=" << (uint32_t)(h1.AMCnum) << " l1AID=" << (uint32_t)(h1.l1AID) << " bxID=" << (uint32_t)(h1.bxID) << " dataLen=" << (uint32_t)(h1.dataLength);
  return ss.str();
}

std::string gem::AMCdata::getAMCheader2_str() const {
  AMCheader2 h2(amch2_);
  std::stringstream ss;
  ss << "AMCheader2: boardID=" << (uint32_t)(h2.boardID) << " orbitNum=" << (uint32_t)(h2.orbitNum) << " param3=" << (uint32_t)(h2.param3) << " param2=" << (uint32_t)(h2.param2) << " param1=" << (uint32_t)(h2.param1) << " runType=" << (uint32_t)(h2.runType) << " formatVer=" << (uint32_t)(h2.formatVer);
  return ss.str();
}

std::string gem::AMCdata::getAMCtrailer_str() const {
  AMCTrailer t(amct_);
  std::stringstream ss;
  ss << "AMCtrailer: dataLenT=" << (uint32_t)(t.dataLengthT) << " cb0=" << (uint32_t)(t.cb0) << " l1AIDT=" << (uint32_t)(t.l1AIDT) << " crc=" << (uint32_t)(t.crc);
  return ss.str();
}

std::string gem::AMCdata::getGEMeventHeader_str() const {
  EventHeader eh(eh_);
  std::stringstream ss;
  ss << "GEMeventHeader: ttsState=" << (uint32_t)(eh.ttsState) << " reserved=" << (uint32_t)(eh.reserved) << " davCnt=" << (uint32_t)(eh.davCnt) << " buffState=" << (uint32_t)(eh.buffState) << " davList=" << (uint32_t)(eh.davList);
  return ss.str();
}

std::string gem::AMCdata::getGEMeventTrailer_str() const {
  EventTrailer et(et_);
  std::stringstream ss;
  ss << "GEMeventTrailer: unused1=" << (uint32_t)(et.unused1) << " ctrl5=" << (uint32_t)(et.ctrl5) << " unused2=" << (uint32_t)(et.unused2) << " oosGlib=" << (uint32_t)(et.oosGlib) << " chTimeOut=" << (uint32_t)(et.chTimeOut);
  return ss.str();
}
