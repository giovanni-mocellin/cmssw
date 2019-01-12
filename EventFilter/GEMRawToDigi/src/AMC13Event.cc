#include <cstdint>
#include "EventFilter/GEMRawToDigi/interface/AMC13Event.h"
#include <sstream>

using namespace gem;

void AMC13Event::setCDFHeader(uint8_t Evt_ty, uint32_t LV1_id, uint16_t BX_id, uint16_t Source_id)
{
  CDFHeader u;
  u.cb5 = 0x5;
  u.eventType = Evt_ty;
  u.lv1Id = LV1_id;
  u.bxId = BX_id;
  u.sourceId = Source_id;
  u.fov = 0; // not used
  cdfh_ = u.word;
}

void AMC13Event::setAMC13Header(uint8_t CalTyp, uint8_t nAMC, uint32_t OrN)
{
  AMC13Header u;
  u.uFov = 0; // format version # FIX ME
  u.calType = CalTyp;
  u.nAMC = nAMC;
  u.reserved0 = 0;
  u.orbitN = OrN;
  u.cb0 = 0x0;
  amc13h_ = u.word;
}

void AMC13Event::setAMC13Trailer(uint8_t Blk_NoT, uint8_t LV1_idT, uint16_t BX_idT)
{
  AMC13Trailer u;
  u.crc32 = 0; // FIX ME
  u.blkN = Blk_NoT;
  u.lv1IdT = LV1_idT;
  u.bxIdT = BX_idT;
  amc13t_ = u.word;
}

void AMC13Event::setCDFTrailer(uint32_t EvtLength)
{
  CDFTrailer u;
  u.cbA = 0xA;
  u.evtType = 0; // FIX ME
  u.evtLength = EvtLength;
  u.crcCDF = 0; // FIX ME
  u.cfxx = 0; // FIX ME
  u.evtStat = 0; // event status # FIX ME
  u.tts = 0;
  u.trdd = 0;
  cdft_ = u.word;
}

void AMC13Event::addAMCheader(uint64_t word)
{
  amcHeaders_.push_back(word);
}

// FIX ME: not verified wrt to AMC13DataDef
void AMC13Event::addAMCheader(uint32_t AMC_size, uint8_t Blk_No, uint8_t AMC_No, uint16_t BoardID)
{
  // AMC Header word 
  // 55 - 32  | 27 - 20 | 19 - 16 | 15 - 0  |
  // AMC_size | Blk_No  | AMC_No  | BoardID |  
  //uint64_t word =
  //  (static_cast<uint64_t>(AMC_size & 0x00ffffff) << 32) |
  //  (static_cast<uint64_t>(Blk_No & 0xff) << 20) |
  //  (static_cast<uint64_t>(AMC_No & 0x0f) << 16) |
  //  (static_cast<uint64_t>(BoardID & 0xffff));
  AMC13DataDef u(0);
  u.dataSize = AMC_size;
  u.blkSeqNo = Blk_No;
  u.amcNo = AMC_No;
  u.boardId = BoardID;
  amcHeaders_.push_back(u.word);
}

std::string AMC13Event::getCDFHeader_str() const
{
  CDFHeader u;
  u.word = cdfh_;
  std::stringstream ss;
  ss << "CDFHeader: fov=" << (uint32_t)(u.fov) << " sourceId=" << (uint32_t)(u.sourceId) << " bxId=" << (uint32_t)(u.bxId) << " lv1Id=" << (uint32_t)(u.lv1Id) << " eventType=" << (uint32_t)(u.eventType) << " cb5=" << (uint32_t)(u.cb5);
  return ss.str();
}

std::string AMC13Event::getCDFTrailer_str() const
{
  CDFTrailer u;
  u.word = cdft_;
  std::stringstream ss;
  ss << "CDFTrailer: trdd=" << (uint32_t)(u.trdd) << " tts=" << (uint32_t)(u.tts) << " evtStat=" << (uint32_t)(u.evtStat) << " cfxx=" << (uint32_t)(u.cfxx) << " crcCDF=" << (uint32_t)(u.crcCDF) << " evtLength=" << (uint32_t)(u.evtLength) << " evtType=" << (uint32_t)(u.evtType) << " cbA=" << (uint32_t)(u.cbA);
  return ss.str();
}

std::string AMC13Event::getAMC13Header_str() const
{
  AMC13Header u;
  u.word = amc13h_;
  std::stringstream ss;
  ss << "AMC13Header: cb0=" << (uint32_t)(u.cb0) << " orbitN=" << (uint32_t)(u.orbitN) << " reserved0=" << (uint32_t)(u.reserved0) << " nAMC=" << (uint32_t)(u.nAMC) << " calType=" << (uint32_t)(u.calType) << " uFov=" << (uint32_t)(u.uFov);
  return ss.str();
}

std::string AMC13Event::getAMC13Trailer_str() const
{
  AMC13Trailer u;
  u.word = amc13t_;
  std::stringstream ss;
  ss << "AMC13Trailer: bxIdT=" << (uint32_t)(u.bxIdT) << " lv1IdT=" << (uint32_t)(u.lv1IdT) << " blkN=" << (uint32_t)(u.blkN) << " crc32=" << (uint32_t)(u.crc32);
  return ss.str();
}

std::string AMC13Event::getAMC13DataDef_str(int i) const
{
  AMC13DataDef d = amcHeaders_.at(i);
  std::stringstream ss;
  ss << "AMC13DataDef: boardId=" << (uint32_t)(d.boardId) << " amcNo=" << (uint32_t)(d.amcNo) << " blkSeqNo=" << (uint32_t)(d.blkSeqNo) << " cb0=" << (uint32_t)(d.cb0) << " dataSize=" << (uint32_t)(d.dataSize) << " endBits8" << (uint32_t)(d.endBits8);
  return ss.str();
}
