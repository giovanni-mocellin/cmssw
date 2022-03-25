#include "L1Trigger/CSCTriggerPrimitives/interface/CSCGEMMatcher.h"
#include "L1Trigger/CSCTriggerPrimitives/interface/GEMInternalCluster.h"
#include "DataFormats/CSCDigi/interface/CSCConstants.h"
#include "DataFormats/CSCDigi/interface/CSCALCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCLCTDigi.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <algorithm>
#include <cmath>

CSCGEMMatcher::CSCGEMMatcher(
    int endcap, unsigned station, unsigned chamber, const edm::ParameterSet& tmbParams, const edm::ParameterSet& conf)
    : endcap_(endcap), station_(station), chamber_(chamber) {
  isEven_ = (chamber_ % 2 == 0);

  maxDeltaBXALCTGEM_ = tmbParams.getParameter<unsigned>("maxDeltaBXALCTGEM");
  maxDeltaBXCLCTGEM_ = tmbParams.getParameter<unsigned>("maxDeltaBXCLCTGEM");

  maxDeltaHsEven_ = tmbParams.getParameter<unsigned>("maxDeltaHsEven");
  maxDeltaHsOdd_ = tmbParams.getParameter<unsigned>("maxDeltaHsOdd");

  matchCLCTpropagation_ = tmbParams.getParameter<bool>("matchCLCTpropagation");

  if (station_ == 1) {
    maxDeltaHsEvenME1a_ = tmbParams.getParameter<unsigned>("maxDeltaHsEvenME1a");
    maxDeltaHsOddME1a_ = tmbParams.getParameter<unsigned>("maxDeltaHsOddME1a");
  }

  mitigateSlopeByCosi_ = tmbParams.getParameter<bool>("mitigateSlopeByCosi");
  assign_gem_csc_bending_ = tmbParams.getParameter<bool>("assignGEMCSCBending");
}

void CSCGEMMatcher::setESLookupTables(const CSCL1TPLookupTableME11ILT* conf) { lookupTableME11ILT_ = conf; }

void CSCGEMMatcher::setESLookupTables(const CSCL1TPLookupTableME21ILT* conf) { lookupTableME21ILT_ = conf; }

//function to replace the CLCT slope by the slope indicated by the strip difference between the CLCT and its matching GEM internal cluster
int CSCGEMMatcher::calculateGEMCSCBending(const CSCCLCTDigi& clct, const GEMInternalCluster& cluster) const {
  const bool isME1a(station_ == 1 and clct.getKeyStrip() > CSCConstants::MAX_HALF_STRIP_ME1B);

  //ME1a necessitates a different treatment because of a different strip numbering scheme and strip width
  const int SignedEighthStripDiff = matchedClusterDistES(clct, cluster);
  const unsigned eighthStripDiff = abs(SignedEighthStripDiff); //LUTs consider only absolute change

  //use LUTs to determine absolute slope, default 0
  int slopeShift = 0;
  if (station_ == 2) {
    if (isEven_) slopeShift = lookupTableME21ILT_->es_diff_slope_L1_ME21_even(eighthStripDiff);
    else         slopeShift = lookupTableME21ILT_->es_diff_slope_L1_ME21_odd(eighthStripDiff);
  } else if (station_ == 1) {
    if (isME1a) {  //is in ME1a
      if (isEven_) slopeShift = lookupTableME11ILT_->es_diff_slope_L1_ME11a_even(eighthStripDiff);
      else 	   slopeShift = lookupTableME11ILT_->es_diff_slope_L1_ME11a_odd(eighthStripDiff);
    } else {
      if (isEven_) slopeShift = lookupTableME11ILT_->es_diff_slope_L1_ME11b_even(eighthStripDiff);
      else         slopeShift = lookupTableME11ILT_->es_diff_slope_L1_ME11b_odd(eighthStripDiff);
    }
  }

  //account for the sign of the difference and take into account whether CLCT slope propagation is on or not
  slopeShift *= pow(-1, std::signbit(SignedEighthStripDiff));
  int NewSlope = matchCLCTpropagation_ ? clct.getSlope() * pow(-1, clct.getBend()) + slopeShift : slopeShift;
  int NewSlopeSign = pow(-1, std::signbit(NewSlope));
  NewSlope = std::min(15, abs(NewSlope)) * NewSlopeSign;

  //Debugging
  //std::cout<<"old slope "<<clct.getSlope() * pow(-1, clct.getBend())<<" vs new slope "<<NewSlope<<std::endl;

  return NewSlope;
}



// match an ALCT to GEMInternalCluster by bunch-crossing
void CSCGEMMatcher::matchingClustersBX(const CSCALCTDigi& alct,
                                       const GEMInternalClusters& clusters,
                                       GEMInternalClusters& output) const {
  if (!alct.isValid() or clusters.empty())
    return;

  // select clusters matched in time
  for (const auto& cl : clusters) {
    const unsigned diff = std::abs(int(alct.getBX()) - cl.bx());
    if (diff <= maxDeltaBXALCTGEM_)
      output.push_back(cl);
  }
}

// match a CLCT to GEMInternalCluster by bunch-crossing
void CSCGEMMatcher::matchingClustersBX(const CSCCLCTDigi& clct,
                                       const GEMInternalClusters& clusters,
                                       GEMInternalClusters& output) const {
  if (!clct.isValid() or clusters.empty())
    return;

  // select clusters matched in time
  for (const auto& cl : clusters) {
    const unsigned diff = std::abs(int(clct.getBX()) - cl.bx());
    if (diff <= maxDeltaBXCLCTGEM_)
      output.push_back(cl);
  }
}

// match an ALCT and CLCT to GEMInternalCluster by bunch-crossing
void CSCGEMMatcher::matchingClustersBX(const CSCALCTDigi& alct,
                                       const CSCCLCTDigi& clct,
                                       const GEMInternalClusters& clusters,
                                       GEMInternalClusters& output) const {
  // both need to be valid
  if (!alct.isValid() or !clct.isValid() or clusters.empty())
    return;

  // get the single matches
  GEMInternalClusters alctClusters, clctClusters;
  matchingClustersBX(alct, clusters, alctClusters);
  matchingClustersBX(clct, clusters, clctClusters);

  // get the intersection
  for (const auto& p : alctClusters) {
    for (const auto& q : clctClusters) {
      if (p == q) {
        output.push_back(p);
      }
    }
  }
}

void CSCGEMMatcher::matchingClustersLoc(const CSCALCTDigi& alct,
                                        const GEMInternalClusters& clusters,
                                        GEMInternalClusters& output) const {
  if (!alct.isValid() or clusters.empty())
    return;

  // select clusters matched in wiregroup
  for (const auto& cl : clusters) {
    // for now add 10 wiregroups to make sure the matching can be done
    // this should be quite generous
    unsigned deltaWG(station_ == 1 ? 10 : 20);
    if (cl.min_wg() <= alct.getKeyWG() and alct.getKeyWG() <= cl.max_wg() + deltaWG) {
      output.push_back(cl);
    }
  }
}

void CSCGEMMatcher::matchingClustersLoc(const CSCCLCTDigi& clct,
                                        const GEMInternalClusters& clusters,
                                        GEMInternalClusters& output) const {
  if (!clct.isValid() or clusters.empty())
    return;

  const bool isME1a(station_ == 1 and clct.getKeyStrip() > CSCConstants::MAX_HALF_STRIP_ME1B);

  //determine window size
  unsigned eighthStripCut;
  if (isEven_) {
    eighthStripCut = 4 * (isME1a ? maxDeltaHsEvenME1a_ : maxDeltaHsEven_); // Cut in 1/8 = 4 * cut in 1/2
  } else {
    eighthStripCut = 4 * (isME1a ? maxDeltaHsOddME1a_ : maxDeltaHsOdd_); // Cut in 1/8 = 4 * cut in 1/2
  }

  // select clusters matched by 1/2-strip or 1/8-strip, picking closest option in eighth-strips
  std::vector<int> distances;
  for (const auto& cl : clusters) {
    const unsigned distanceES = abs(matchedClusterDistES(clct, cl));
    if (distanceES <= eighthStripCut){ //only accept clusters in the window around the CLCT

      //assure building an ordered in increasing absolute distance list of GEM matches
      unsigned ListPosition = output.size();
      for (unsigned i = 0; i < output.size(); ++i){
	if (distanceES < abs(distances[i])) ListPosition = i;
	break;
      }
      distances.insert(distances.begin() + ListPosition, distanceES);
      output.insert(output.begin() + ListPosition, cl);
    }
  }
}

// calculate distance in eighth-strip units between CLCT and GEM
int CSCGEMMatcher::matchedClusterDistES(const CSCCLCTDigi& clct, const GEMInternalCluster& cl) const {
  const bool isME1a(station_ == 1 and clct.getKeyStrip() > CSCConstants::MAX_HALF_STRIP_ME1B);

  int cl_es = isME1a ? cl.getKeyStripME1a(8) : cl.getKeyStrip(8);

  int eighthStripDiff = cl_es - clct.getKeyStrip(8);

  //Debugging
  //std::cout<<"Diff GEM "<<cl_es<<" - CSC "<<clct.getKeyStrip(8)<<" = "<<eighthStripDiff<<std::endl;

  if (matchCLCTpropagation_) { //modification of DeltaStrip by CLCT slope
    int SlopeShift = 0;
    uint16_t baseSlope = -1;
    baseSlope =  mitigateSlopeByCosi_ ? mitigatedSlopeByConsistency(clct) : clct.getSlope();

    int clctSlope = pow(-1, clct.getBend()) * baseSlope;

    SlopeShift = CSCGEMSlopeCorrector(isME1a, clctSlope);
    //Debugging
    //std::cout<<"SlopeShift = "<<SlopeShift<<" at slope "<<clctSlope<<std::endl;
    eighthStripDiff -= SlopeShift;
  }

  //Debugging
  //uint16_t mCOSIslope = mitigatedSlopeByConsistency(clct);
  //std::cout<<"COSI slope = "<<pow(-1, clct.getBend()) * mCOSIslope<<std::endl;

  return eighthStripDiff;
}

void CSCGEMMatcher::matchingClustersLoc(const CSCALCTDigi& alct,
                                        const CSCCLCTDigi& clct,
                                        const GEMInternalClusters& clusters,
                                        GEMInternalClusters& output) const {
  // both need to be valid
  if (!alct.isValid() or !clct.isValid() or clusters.empty())
    return;

  // get the single matches
  GEMInternalClusters alctClusters, clctClusters;
  matchingClustersLoc(alct, clusters, alctClusters);
  matchingClustersLoc(clct, clusters, clctClusters);

  // get the intersection
  for (const auto& p : clctClusters) { //start with clct clusters, as this list is ordered in rising distance
    for (const auto& q : alctClusters) {
      if (p == q) {
        output.push_back(p);
      }
    }
  }
}

void CSCGEMMatcher::matchingClustersBXLoc(const CSCALCTDigi& alct,
                                          const GEMInternalClusters& clusters,
                                          GEMInternalClusters& output) const {
  if (!alct.isValid() or clusters.empty()){
    edm::LogError("CSCGEMMatcher") << "either clusters empty or alct invalid \n";
    return;
  }

  // match by BX
  GEMInternalClusters clustersBX;
  matchingClustersBX(alct, clusters, clustersBX);

  // match spatially
  matchingClustersLoc(alct, clustersBX, output);
}

void CSCGEMMatcher::matchingClustersBXLoc(const CSCCLCTDigi& clct,
                                          const GEMInternalClusters& clusters,
                                          GEMInternalClusters& output) const {
  if (!clct.isValid() or clusters.empty())
    return;

  // match by BX
  GEMInternalClusters clustersBX;
  matchingClustersBX(clct, clusters, clustersBX);

  // match spatially
  matchingClustersLoc(clct, clustersBX, output);
}

void CSCGEMMatcher::matchingClustersBXLoc(const CSCALCTDigi& alct,
                                          const CSCCLCTDigi& clct,
                                          const GEMInternalClusters& clusters,
                                          GEMInternalClusters& selected) const {
  // both need to be valid
  if (!alct.isValid() or !clct.isValid() or clusters.empty()){
    edm::LogError("CSCGEMMatcher") << "either clusters empty or alct invalid \n";
    return;
  }

  // match by BX
  GEMInternalClusters clustersBX;
  matchingClustersBX(alct, clct, clusters, clustersBX);

  // match spatially
  matchingClustersLoc(alct, clct, clustersBX, selected);
}

void CSCGEMMatcher::bestClusterBXLoc(const CSCALCTDigi& alct,
                                     const GEMInternalClusters& clusters,
                                     GEMInternalCluster& best) const {
  if (!alct.isValid() or clusters.empty())
    return;

  GEMInternalClusters clustersBXLoc;
  matchingClustersBXLoc(alct, clusters, clustersBXLoc);

  // simply pick the first matching one
  if (!clustersBXLoc.empty())
    best = clustersBXLoc[0];
}

void CSCGEMMatcher::bestClusterBXLoc(const CSCCLCTDigi& clct,
                                     const GEMInternalClusters& clusters,
                                     GEMInternalCluster& best) const {
  if (!clct.isValid() or clusters.empty())
    return;

  // match by BX
  GEMInternalClusters clustersBXLoc;
  matchingClustersBXLoc(clct, clusters, clustersBXLoc);

  // the first matching one is also the closest in phi distance (to expected position, if extrapolating), by ordered list in CLCT matching
  if (!clustersBXLoc.empty())
    best = clustersBXLoc[0];
}

void CSCGEMMatcher::bestClusterBXLoc(const CSCALCTDigi& alct,
                                     const CSCCLCTDigi& clct,
                                     const GEMInternalClusters& clusters,
                                     GEMInternalCluster& best) const {
  // match by BX
  GEMInternalClusters clustersBXLoc;
  matchingClustersBXLoc(alct, clct, clusters, clustersBXLoc);

  // the first matching one is also the closest in phi distance (to expected position, if extrapolating), by ordered list in CLCT matching
  if (!clustersBXLoc.empty())
    best = clustersBXLoc[0];
}

// function to determine CLCT consistency of slope indicator (COSI) and use it to mitigate slope according to LUT
uint16_t CSCGEMMatcher::mitigatedSlopeByConsistency(const CSCCLCTDigi& clct) const {
  
  const bool isME1a(station_ == 1 and clct.getKeyStrip() > CSCConstants::MAX_HALF_STRIP_ME1B);

  // extract hit values from CLCT hit matrix
  std::vector<std::vector<uint16_t>> CLCTHitMatrix = clct.getHits();
  int CLCTHits[6] = {-1, -1, -1, -1, -1, -1};

  for (unsigned layer = 0; layer < CLCTHitMatrix.size(); ++layer) {
    for (unsigned position = 0; position < CLCTHitMatrix.at(layer).size(); ++position) {
      const uint16_t value = CLCTHitMatrix.at(layer).at(position);
      if (value != 0 && value != 65535) {
        CLCTHits[layer] = (int)value;
        break;
      }
    }
  }

  //Debugging
  //std::cout<<"CLCT Hits = "<<CLCTHits[0]<<", "<<CLCTHits[1]<<", "<<CLCTHits[2]<<", "<<CLCTHits[3]<<", "<<CLCTHits[4]<<", "<<CLCTHits[5]<<std::endl;

  //calculate slope consistency
  float MinMaxPairDifferences[2] = {999., -999.};
  for (unsigned First = 0; First < 5; ++First) {
    //skip empty layers
    if (CLCTHits[First] == -1)
      continue;
    for (unsigned Second = First + 1; Second < 6; ++Second) {
      //skip empty layers
      if (CLCTHits[Second] == -1)
        continue;
      float PairDifference = (CLCTHits[First] - CLCTHits[Second]) / (float)(Second - First);
      if (PairDifference < MinMaxPairDifferences[0])
        MinMaxPairDifferences[0] = PairDifference;
      if (PairDifference > MinMaxPairDifferences[1])
        MinMaxPairDifferences[1] = PairDifference;
    }
  }

  //calculate consistency of slope indicator: cosi
  uint16_t cosi = std::ceil(std::abs(MinMaxPairDifferences[1] - MinMaxPairDifferences[0]));
  //Debugging
  //std::cout<<"COSI = "<<cosi<<std::endl;

  //disambiguate cosi cases

  //extremely inconsistent track, deprecate slope
  if (cosi > 3)
    return 0;
  //consistent slope, do not change
  else if (cosi < 2)
    return clct.getSlope();
  //need to look up in table 2->1
  else if (cosi == 2) {
    if (station_ == 1){
      if (isME1a){
        if (chamber_ % 2 == 0)
          return lookupTableME11ILT_->CSC_slope_cosi_2to1_L1_ME11a_even(clct.getSlope());
        else
          return lookupTableME11ILT_->CSC_slope_cosi_2to1_L1_ME11a_odd(clct.getSlope());
      } else {
        if (chamber_ % 2 == 0)
          return lookupTableME11ILT_->CSC_slope_cosi_2to1_L1_ME11b_even(clct.getSlope());
        else
          return lookupTableME11ILT_->CSC_slope_cosi_2to1_L1_ME11b_odd(clct.getSlope());
      }
    } else {
      if (chamber_ % 2 == 0)
        return lookupTableME21ILT_->CSC_slope_cosi_2to1_L1_ME21_even(clct.getSlope());
      else
        return lookupTableME21ILT_->CSC_slope_cosi_2to1_L1_ME21_odd(clct.getSlope());
    }
  }
  //need to look up in table 3->1
  else if (cosi == 3) {
    if (station_ == 1){
      if (isME1a) {
        if (chamber_ % 2 == 0)
          return lookupTableME11ILT_->CSC_slope_cosi_3to1_L1_ME11a_even(clct.getSlope());
        else
          return lookupTableME11ILT_->CSC_slope_cosi_3to1_L1_ME11a_odd(clct.getSlope());
        } else {
        if (chamber_ % 2 == 0)
          return lookupTableME11ILT_->CSC_slope_cosi_3to1_L1_ME11b_even(clct.getSlope());
        else
          return lookupTableME11ILT_->CSC_slope_cosi_3to1_L1_ME11b_odd(clct.getSlope());
      }
    } else {
      if (chamber_ % 2 == 0)
        return lookupTableME21ILT_->CSC_slope_cosi_3to1_L1_ME21_even(clct.getSlope());
      else
        return lookupTableME21ILT_->CSC_slope_cosi_3to1_L1_ME21_odd(clct.getSlope());
    }
  }
  //just to avoid compiler errors an error code
  else {
    return 999;
  }
}

//function to correct expected GEM position in phi by CSC slope measurement
int CSCGEMMatcher::CSCGEMSlopeCorrector(bool isME1a, int cscSlope) const {
  int SlopeShift = 0;
  int SlopeSign = pow(-1, std::signbit(cscSlope));
  //account for slope mitigation by cosi, if opted-in
  if (mitigateSlopeByCosi_) {
    if (station_ == 1) {
      if (chamber_ % 2 == 0) SlopeShift = isME1a ? lookupTableME11ILT_->CSC_slope_cosi_corr_L1_ME11a_even(std::abs(cscSlope))
		 				 : lookupTableME11ILT_->CSC_slope_cosi_corr_L1_ME11b_even(std::abs(cscSlope));
      else		     SlopeShift = isME1a ? lookupTableME11ILT_->CSC_slope_cosi_corr_L1_ME11a_odd(std::abs(cscSlope))
 						 : lookupTableME11ILT_->CSC_slope_cosi_corr_L1_ME11b_odd(std::abs(cscSlope));
    }
    else if (station_ == 2) {
      if (chamber_ % 2 == 0) SlopeShift = lookupTableME21ILT_->CSC_slope_cosi_corr_L1_ME21_even(std::abs(cscSlope));
      else                   SlopeShift = lookupTableME21ILT_->CSC_slope_cosi_corr_L1_ME21_odd(std::abs(cscSlope));
    }
  }
  else{ //account for slope without mitigation, if opted out
    if (station_ == 1) {
      if (chamber_ % 2 == 0) SlopeShift = isME1a ? lookupTableME11ILT_->CSC_slope_corr_L1_ME11a_even(std::abs(cscSlope))
                                                 : lookupTableME11ILT_->CSC_slope_corr_L1_ME11b_even(std::abs(cscSlope));
      else                   SlopeShift = isME1a ? lookupTableME11ILT_->CSC_slope_corr_L1_ME11a_odd(std::abs(cscSlope))
                                                 : lookupTableME11ILT_->CSC_slope_corr_L1_ME11b_odd(std::abs(cscSlope));
    }
    else if (station_ == 2) {
      if (chamber_ % 2 == 0) SlopeShift = lookupTableME21ILT_->CSC_slope_corr_L1_ME21_even(std::abs(cscSlope));
      else                   SlopeShift = lookupTableME21ILT_->CSC_slope_corr_L1_ME21_odd(std::abs(cscSlope));
    }     
  }
  return std::round(SlopeShift * SlopeSign);
}
