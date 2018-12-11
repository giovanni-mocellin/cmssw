#ifndef CondFormats_GEMObjects_GEMQC8Conf_h
#define CondFormats_GEMObjects_GEMQC8Conf_h

/** \class GEMQC8Conf
 *
 *  A class contains definition of the QC8 (cosmic stand) configuration
 *  It is used for GEMQC8ConfRcd
 *
 *  \author A. Juodagalvis - Vilnius University
 *  Oct 2018
 */


#include "CondFormats/Serialization/interface/Serializable.h"
#include "CondFormats/GEMObjects/interface/GEMELMap.h"
#include <vector>
#include <iostream>

class GEMQC8Conf {

 public:
  GEMQC8Conf();
  GEMQC8Conf(const GEMQC8Conf *);
  ~GEMQC8Conf(){}

  int runNo() const { return run_number_; }
  const std::vector<std::string>& chSerNums() const { return chSerialNums_; }
  const std::string& chSerNum(int idx) const { return chSerialNums_.at(idx); }
  const std::vector<std::string>& chPositions() const { return chPositions_; }
  const std::string& chPos(int idx) const { return chPositions_.at(idx); }
  const std::vector<float> chGasFlow() const { return chGasFlow_; }
  float chGasFlow(int idx) const { return chGasFlow_.at(idx); }
  int hasELMap() const { return hasELMap_; }
  const std::vector<GEMELMap>& elMap() const { return elMap_; }
  const GEMELMap& elMap(int idx) const { return elMap_.at(idx); }

  void assign(const GEMQC8Conf&);
  void print(std::ostream &out = std::cout, int detailed=0) const;

 public:
  int run_number_;
  std::vector<std::string> chSerialNums_;
  std::vector<std::string> chPositions_;
  std::vector<float> chGasFlow_;

  // optional electronics map
  int hasELMap_;
  std::vector<GEMELMap> elMap_;

  COND_SERIALIZABLE;
};
#endif
