#ifndef CondTools_GEM_GEMQC8ConfSourceHandler_h
#define CondTools_GEM_GEMQC8ConfSourceHandler_h

/** \class GEMQC8ConfSourceHandler
 *
 *  Used for interaction with DB to get GEMQC8Conf and GEMELMap records
 *  Relies on GEMELMapHelper to fill vfatpos->strip2chan mapping
 *
 *  \author A. Juodagalvis - Vilnius University
 *  Oct 2018
 */


#include "CondCore/PopCon/interface/PopConSourceHandler.h"
#include "CondFormats/GEMObjects/interface/GEMQC8Conf.h"
#include "CondFormats/GEMObjects/interface/GEMELMap.h"
#include "CondFormats/GEMObjects/interface/GEMELMapHelper.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>

namespace popcon
{
  class GEMQC8ConfSourceHandler : public popcon::PopConSourceHandler<GEMQC8Conf>
  {

  public:

    GEMQC8ConfSourceHandler( const edm::ParameterSet& ps );
    ~GEMQC8ConfSourceHandler();
    void getNewObjects();
    void ConnectOnlineDB( const std::string& connect, const edm::ParameterSet& connectionPset );
    void DisconnectOnlineDB();
    void readGEMQC8Conf();
    void readGEMQC8EMap();
    void readGEMQC8EMap_specDummy();
    std::string id() const { return m_name; }
    const GEMQC8Conf* getQC8conf() const { return qc8conf; }
    const GEMELMap* getQC8elMap() const { return qc8elMap; }
    void setOnlyConfDef(int val) { m_onlyConfDef = val; }

    // temporary function
    int gemELMap_vfat_autoFill(const std::string &chamberName, GEMELMap::GEMVFatMap &vfats, GEMELMap::GEMStripMap &stripMap, GEMELMapHelper &elMapAid);

  private:
    GEMQC8Conf * qc8conf;
    GEMELMap * qc8elMap;
    cond::persistency::Session session;
    std::string m_name;
    int m_dummy;
    int m_debugMode;
    std::string m_connect;
    std::string m_authpath;
    edm::ParameterSet m_connectionPset;
    int m_runNumber;
    int m_allowRollBack;
    int m_noDBOutput;   // whether DBOutput module is expected
    int m_onlyConfDef;
    bool m_printValues;
  };
}
#endif
