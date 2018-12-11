#include "EventFilter/GEMRawToDigi/interface/GEBdata.h"
#include <string>
#include <sstream>

std::string gem::GEBdata::getChamberHeader_str() const {
  GEBchamberHeader ch(ch_);
  std::stringstream ss;
  ss << "ChamberHeader: unused1=" << (uint32_t)(ch.unused1) << " inputStatus=" << (uint32_t)(ch.inputStatus) << " [flag_noVFATmarker=" << flag_noVFATmarker() << "]" << " vfatWordCnt=" << (uint32_t)(ch.vfatWordCnt) << " inputID=" << (uint32_t)(ch.inputID) << " zeroSupWordsCnt=" << (uint32_t)(ch.zeroSupWordsCnt) << " unused2=" << (uint32_t)(ch.unused2);
  return ss.str();
}

std::string gem::GEBdata::getChamberTrailer_str() const {
  GEBchamberTrailer ct(ct_);
  std::stringstream ss;
  ss << "ChamberTrailer: unused=" << (uint32_t)(ct.unused) << " inUfw=" << (uint32_t)(ct.inUfw) << " stuckData=" << (uint32_t)(ct.stuckData) << " inFIFOund=" << (uint32_t)(ct.inFIFOund) << " vfatWordCntT=" << (uint32_t)(ct.vfatWordCntT) << " ohcrc=" << (uint32_t)(ct.ohcrc);
  return ss.str();
}
