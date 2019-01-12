#include "CondFormats/GEMObjects/interface/GEMQC8Conf.h"
#include "FWCore/Utilities/interface/Exception.h"

GEMQC8Conf::GEMQC8Conf() :
  run_number_(0), chSerialNums_(),
  chPositions_(), chGasFlow_(),
  hasELMap_(0), elMap_()
{
}

GEMQC8Conf::GEMQC8Conf(const GEMQC8Conf *g) :
  run_number_(0), chSerialNums_(),
  chPositions_(), chGasFlow_(),
  hasELMap_(0), elMap_()
{
  if (!g) return;
  run_number_ = g->run_number_;
  chSerialNums_ = g->chSerialNums_;
  chPositions_ = g->chPositions_;
  chGasFlow_ = g->chGasFlow_;
  hasELMap_ = g->hasELMap_;
  if (hasELMap_) elMap_ = g->elMap_;
}

void GEMQC8Conf::assign(const GEMQC8Conf &g)
{
  run_number_ = g.run_number_;
  chSerialNums_ = g.chSerialNums_;
  chPositions_ = g.chPositions_;
  chGasFlow_ = g.chGasFlow_;
  hasELMap_ = g.hasELMap_;
  elMap_ = g.elMap_;
}


void GEMQC8Conf::print(std::ostream &out, int detailed) const
{
  out << "GEMQC8Conf for run_number=" << run_number_ << " (hasElMap="
      << hasELMap_ << ")\n";
  unsigned int sz= chSerialNums_.size();
  if ((sz!=chPositions_.size()) || (sz!=chGasFlow_.size())) {
    out << " different sizes of arrays: chSerialNums[" << sz
	<< "], chPositions[" << chPositions_.size() << "], chGasFlow["
	<< chGasFlow_.size() << "]\n";
  }
  else {
    for (unsigned int i=0; i<chSerialNums_.size(); i++) {
      out << " " << chPositions_.at(i) << " " << chSerialNums_.at(i) << "  "
	  << chGasFlow_.at(i) << "\n";
    }
  }

  if (hasELMap_) {
    out << " electronics map [" << elMap_.size() << "]\n";
    if (elMap_.size()!=chSerialNums_.size()) {
      out << " ERROR: there are " << chSerialNums_.size() << " chambers and "
	  << elMap_.size() << " elMap elements\n";
    }
    else {
      for (unsigned int i=0; i<elMap_.size(); i++) {
	out << " elMap@" << i << " (chamber=" << chSerialNums_[i] <<") ";
	elMap_[i].print(out,detailed);
      }
    }
  }
}
