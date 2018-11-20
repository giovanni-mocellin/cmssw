/** \unpacker for gem
 *  \author J. Lee - UoS
 */
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"

#include "EventFilter/GEMRawToDigi/plugins/GEMRawToDigiModule.h"

#include <bitset>

using namespace gem;

GEMRawToDigiModule::GEMRawToDigiModule(const edm::ParameterSet & pset) :
  fed_token(consumes<FEDRawDataCollection>( pset.getParameter<edm::InputTag>("InputLabel") )),
  useDBEMap_(pset.getParameter<bool>("useDBEMap")),
  unPackStatusDigis_(pset.getParameter<bool>("unPackStatusDigis"))
{
  edm::InputTag inp= pset.getParameter<edm::InputTag>("InputLabel");
  //std::cout << "GEMRawToDigiModule: inp = " << inp << std::endl;

  produces<GEMDigiCollection>(); 
  if (unPackStatusDigis_) {
    produces<GEMVfatStatusDigiCollection>("vfatStatus");
    produces<GEMGEBdataCollection>("gebStatus");
    produces<GEMAMCdataCollection>("AMCdata"); 
    produces<GEMAMC13EventCollection>("AMC13Event"); 
  }
}

void GEMRawToDigiModule::fillDescriptions(edm::ConfigurationDescriptions & descriptions)
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("InputLabel", edm::InputTag("rawDataCollector")); 
  desc.add<bool>("useDBEMap", false);
  desc.add<bool>("unPackStatusDigis", false);
  descriptions.add("muonGEMDigisDefault", desc);  
}

std::shared_ptr<GEMROmap> GEMRawToDigiModule::globalBeginRun(edm::Run const&, edm::EventSetup const& iSetup) const
{
  auto gemORmap = std::make_shared<GEMROmap>();
  if (useDBEMap_) {
    edm::ESHandle<GEMELMap> gemEMapRcd;
    iSetup.get<GEMELMapRcd>().get(gemEMapRcd);
    auto gemEMap = std::make_unique<GEMELMap>(*(gemEMapRcd.product()));
    gemEMap->convert(*gemORmap);
    gemEMap.reset();    
  }
  else {
    // no EMap in DB, using dummy
    auto gemEMap = std::make_unique<GEMELMap>();
    gemEMap->convertDummy(*gemORmap);
    gemEMap.reset();    
  }
  return gemORmap;
}

void GEMRawToDigiModule::produce(edm::StreamID iID, edm::Event & iEvent, edm::EventSetup const&) const
{
  auto outGEMDigis = std::make_unique<GEMDigiCollection>();
  auto outVFATStatus = std::make_unique<GEMVfatStatusDigiCollection>();
  auto outGEBStatus = std::make_unique<GEMGEBdataCollection>();
  auto outAMCdata = std::make_unique<GEMAMCdataCollection>();
  auto outAMC13Event = std::make_unique<GEMAMC13EventCollection>();

  // Take raw from the event
  edm::Handle<FEDRawDataCollection> fed_buffers;
  iEvent.getByToken( fed_token, fed_buffers );
  
  auto gemROMap = runCache(iEvent.getRun().index());
  
  for (unsigned int id=FEDNumbering::MINGEMFEDID; id<=FEDNumbering::MAXGEMFEDID; ++id) { 
    //if (id!=1472) continue;
    const FEDRawData& fedData = fed_buffers->FEDData(id);
    
    int nWords = fedData.size()/sizeof(uint64_t);
    LogDebug("GEMRawToDigiModule") <<" words " << nWords;

    if (nWords<5) continue;
    //std::cout << "GEMRawToDigiModule: nWords= " << nWords << std::endl;
    const unsigned char * data = fedData.data();
    
    auto amc13Event = std::make_unique<AMC13Event>();
    
    const uint64_t* word = reinterpret_cast<const uint64_t* >(data);
    
    amc13Event->setCDFHeader(*word);
    amc13Event->setAMC13Header(*(++word));

    std::cout << "GEMRawToDigi: cb5=" << amc13Event->get_cb5() << std::endl;
    std::cout << "GEMRawToDigi: cb0=" << amc13Event->get_cb0() << std::endl;
    std::cout << "nAMC=" << amc13Event->get_nAMC() << std::endl;
    
    // Readout out AMC headers
    for (uint8_t i = 0; i < amc13Event->nAMC(); ++i) {
      amc13Event->addAMCheader(*(++word));
      std::cout << "i=" << i << " cb0=" << amc13Event->get_cb0(i) << " amcNr=" << amc13Event->get_amcNr(i) << " dataSize=" << amc13Event->get_dataSize(i) << std::endl;
    }
    
    // Readout out AMC payloads
    for (uint8_t i = 0; i < amc13Event->nAMC(); ++i) {
      auto amcData = std::make_unique<AMCdata>();
      amcData->setAMCheader1(*(++word));      
      amcData->setAMCheader2(*(++word));
      amcData->setGEMeventHeader(*(++word));
      uint16_t amcId = amcData->boardId();
      uint16_t amcBx = amcData->bx();

      std::cout << "iamc=i=" << i << " amcData->amcNr=" << (uint32_t)(amcData->amcNum()) << " amcData->davCnt=" << (uint32_t)(amcData->davCnt()) << std::endl;

      // Fill GEB
      std::cout << "amcData->davCnt=" << amcData->get_davCnt() << std::endl;
      for (uint8_t j = 0; j < amcData->davCnt(); ++j) {
	auto gebData = std::make_unique<GEBdata>();
	gebData->setChamberHeader(*(++word));

	std::cout << "igeb=j=" << j << " vfatWordCnt=" << gebData->get_vfatWordCnt() << std::endl;
	uint16_t gebId = gebData->inputID();
	uint16_t vfatId=0;
	GEMROmap::eCoord geb_ec = {amcId, gebId, vfatId};
	std::cout << "amcId=" << amcId << " gebId=" << gebId << " vfatId=" << vfatId << std::endl;
	if (!gemROMap->isValidChipID(geb_ec)) {
	  std::cout << "is not valid geb_ec chipID" << std::endl;
	  gemROMap->printElDetMap(std::cout);
	}
	GEMROmap::dCoord geb_dc = gemROMap->hitPosition(geb_ec);
	GEMDetId gemId = geb_dc.gemDetId;

	std::cout << "gebData->vfatWordCnt()= " << gebData->get_vfatWordCnt() << std::endl;
	for (uint16_t k = 0; k < gebData->vfatWordCnt()/3; k++) {
	  auto vfatData = std::make_unique<VFATdata>();
	  vfatData->read_fw(*(++word));
	  vfatData->read_sw(*(++word));
	  vfatData->read_tw(*(++word));
	  
	  if (geb_dc.vfatType < 10) {
	    // vfat v2
	    vfatId = vfatData->chipID();
	    vfatData->setVersion(2);
	  }
	  else {
	    // vfat v3
	    vfatId = vfatData->position();
	    vfatData->setVersion(3);
	  }
	  
	  uint16_t bc=vfatData->bc();
	  // strip bx = vfat bx - amc bx
	  int bx = bc-amcBx;
	  
	  if (vfatData->quality()) {
	    edm::LogWarning("GEMRawToDigiModule") << "Quality "<< vfatData->quality()
						  << " b1010 "<< int(vfatData->b1010())
						  << " b1100 "<< int(vfatData->b1100())
						  << " b1110 "<< int(vfatData->b1110());
	    //if (vfatData->crc() != vfatData->checkCRC() ) {
	    //  edm::LogWarning("GEMRawToDigiModule") << "DIFFERENT CRC :"
	    //					    <<vfatData->crc()<<"   "<<vfatData->checkCRC();
	    //}
	  }
	  
	  //check if ChipID exists.
	  GEMROmap::eCoord ec = {amcId, gebId, vfatId};
	  if (!gemROMap->isValidChipID(ec)) {
	    edm::LogWarning("GEMRawToDigiModule") << "InValid: amcId "<<ec.amcId
						  << " gebId "<< ec.gebId
						  << " vfatId "<< ec.vfatId;
	    continue;
	  }

	  GEMROmap::dCoord dc = gemROMap->hitPosition(ec);
	  gemId = dc.gemDetId;
	  vfatData->setPhi(dc.locPhi);

	  for (int chan = 0; chan < VFATdata::nChannels; ++chan) {
	    uint8_t chan0xf = 0;
	    if (chan < 64) chan0xf = ((vfatData->lsData() >> chan) & 0x1);
	    else chan0xf = ((vfatData->msData() >> (chan-64)) & 0x1);

	    // no hits
	    if (chan0xf==0) continue;
	    	             
            GEMROmap::channelNum chMap = {dc.vfatType, chan};
            GEMROmap::stripNum stMap = gemROMap->hitPosition(chMap);

            int stripId = stMap.stNum + vfatData->phi()*GEMELMap::maxChan_;    

	    GEMDigi digi(stripId,bx);
	    LogDebug("GEMRawToDigiModule") <<" vfatId "<<ec.vfatId
					   <<" gemDetId "<< gemId
					   <<" chan "<< chMap.chNum
					   <<" strip "<< stripId
					   <<" bx "<< digi.bx();
	    
	    outGEMDigis.get()->insertDigi(gemId,digi);
	    
	  }// end of channel loop
	  
	  if (unPackStatusDigis_) {
            outVFATStatus.get()->insertDigi(gemId, GEMVfatStatusDigi(*vfatData));
	  }
	  
	} // end of vfat loop
	//std::cout << "end of vfat loop" << std::endl;
	
	gebData->setChamberTrailer(*(++word));
	
        if (unPackStatusDigis_) {
	  outGEBStatus.get()->insertDigi(gemId.chamberId(), (*gebData)); 
        }
	
      } // end of geb loop
      
      amcData->setGEMeventTrailer(*(++word));
      amcData->setAMCTrailer(*(++word));
      
      if (unPackStatusDigis_) {
        outAMCdata.get()->insertDigi(amcData->boardId(), (*amcData));
      }

    } // end of amc loop
    
    amc13Event->setAMC13Trailer(*(++word));
    amc13Event->setCDFTrailer(*(++word));
     
    if (unPackStatusDigis_) {
      outAMC13Event.get()->insertDigi(amc13Event->bxId(), AMC13Event(*amc13Event));
    }
    
  } // end of amc13Event
  
  iEvent.put(std::move(outGEMDigis));
  
  if (unPackStatusDigis_) {
    iEvent.put(std::move(outVFATStatus), "vfatStatus");
    iEvent.put(std::move(outGEBStatus), "gebStatus");
    iEvent.put(std::move(outAMCdata), "AMCdata");
    iEvent.put(std::move(outAMC13Event), "AMC13Event");
  }
  
}
