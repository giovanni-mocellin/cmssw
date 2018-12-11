#ifndef EventFilter_GEMRawToDigi_GEBdata_h
#define EventFilter_GEMRawToDigi_GEBdata_h
#include "VFATdata.h"
#include <vector>
#include <string>

namespace gem {
  // Input status 1 bit for each
  // BX mismatch GLIB OH / BX mismatch GLIB VFAT / OOS GLIB OH / OOS GLIB VFAT / No VFAT marker
  // Event size warn / L1AFIFO near full / InFIFO near full / EvtFIFO near full / Event size overflow
  // L1AFIFO full / InFIFO full / EvtFIFO full
  union GEBchamberHeader {
    GEBchamberHeader(uint64_t x=0) : word(x) {}

    uint64_t word;
    struct {
      uint64_t unused1         : 10;
      uint64_t inputStatus     : 13; // Input status (critical) only highest 13 bits used
      uint64_t vfatWordCnt     : 12; // Size of VFAT payload in 64bit words expected to send to AMC13
      uint64_t inputID         : 5 ; // GLIB input ID (starting at 0)
      uint64_t zeroSupWordsCnt : 12; // Bitmask indicating if certain VFAT blocks have been zero suppressed
      uint64_t unused2         : 12;
    };
  };

  union GEBchamberTrailer {      
    GEBchamberTrailer(uint64_t x=0) : word(x) {}

    uint64_t word;
    struct {
      uint64_t unused       : 33;
      uint64_t inUfw        : 1;
      uint64_t stuckData    : 1;  // Input status (warning): There was data in InFIFO or EvtFIFO when L1A FIFO was empty
      uint64_t inFIFOund    : 1;  // Input status (critical): Input FIFO underflow occurred while sending this event
      uint64_t vfatWordCntT : 12; // Size of VFAT payload in 64bit words sent to AMC13 
      uint64_t ohcrc        : 16; // CRC of OH data (currently not available – filled with 0)
    };
  };

  class GEBdata
  {
  public:
    
    GEBdata() {};
    ~GEBdata() {vfatd_.clear();}

    //!Read chamberHeader from the block.
    void setChamberHeader(uint64_t word) { ch_ = word;}
    uint64_t getChamberHeader() const { return ch_;}

    //!Read chamberTrailer from the block.
    void setChamberTrailer(uint64_t word) { ct_ = word;}
    uint64_t getChamberTrailer() const { return ct_;}

    std::string getChamberHeader_str() const;
    std::string getChamberTrailer_str() const;

    uint16_t unused1()         const {return GEBchamberHeader{ch_}.unused1;}
    uint16_t inputStatus()     const {return GEBchamberHeader{ch_}.inputStatus;}
    int      flag_noVFATmarker()    const {return (((uint16_t)(GEBchamberHeader{ch_}.inputStatus) & ((uint16_t)(1)<<5))!=0) ? 1:0;}
    uint16_t vfatWordCnt()     const {return GEBchamberHeader{ch_}.vfatWordCnt;}
    uint8_t  inputID()         const {return GEBchamberHeader{ch_}.inputID;}
    uint16_t zeroSupWordsCnt() const {return GEBchamberHeader{ch_}.zeroSupWordsCnt;}
    uint16_t unused2()         const {return GEBchamberHeader{ch_}.unused2;}

    uint32_t get_inputID()     const {return (uint32_t)(GEBchamberHeader{ch_}.inputID);}

    //uint32_t ecOH()            const {return ct_.ecOH;}
    //uint16_t bcOH()            const {return ct_.bcOH;}
    uint64_t unused()          const {return GEBchamberTrailer{ct_}.unused;}
    uint8_t  inUfw()           const {return GEBchamberTrailer{ct_}.inUfw;}
    uint8_t  stuckData()       const {return GEBchamberTrailer{ct_}.stuckData;}
    uint8_t  inFIFOund()       const {return GEBchamberTrailer{ct_}.inFIFOund;}
    uint16_t vfatWordCntT()    const {return GEBchamberTrailer{ct_}.vfatWordCntT;}
    uint16_t ohcrc()           const {return GEBchamberTrailer{ct_}.ohcrc;}

    uint32_t get_vfatWordCnt() const {return (uint32_t)(GEBchamberHeader{ch_}.vfatWordCnt);}
    uint32_t get_vfatWordCntT()const {return (uint32_t)(GEBchamberTrailer{ct_}.vfatWordCntT);}

    void setVfatWordCnt(uint16_t n)
    { GEBchamberHeader h(ch_); h.vfatWordCnt=n; ch_=h.word;
      GEBchamberTrailer t(ct_); t.vfatWordCntT=n; ct_=t.word;}

    void setInputID(uint8_t n)
    {GEBchamberHeader h; h.word=ch_; h.inputID=n; ch_=h.word;}

    //!Adds VFAT data to the vector
    void addVFAT(VFATdata v) {vfatd_.push_back(v);}
    //!Returns the vector of FVAT data
    const std::vector<VFATdata> * vFATs() const {return &vfatd_;}  

    static const int sizeGebID = 5;
    
  private:
    //GEBchamberHeader ch_;
    //GEBchamberTrailer ct_;
    uint64_t ch_;
    uint64_t ct_;

    std::vector<VFATdata> vfatd_;     ///<Vector of VFAT data
  };
}
#endif
