#ifndef gemcrValidation_H
#define gemcrValidation_H

#include "Validation/MuonGEMHits/interface/GEMBaseValidation.h"
#include <DataFormats/GEMRecHit/interface/GEMRecHit.h>
#include <DataFormats/GEMRecHit/interface/GEMRecHitCollection.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackingRecHit/interface/TrackingRecHit.h>
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMChamber.h"

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "RecoMuon/CosmicMuonProducer/interface/CosmicMuonSmoother.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"

#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "RecoMuon/StandAloneTrackFinder/interface/StandAloneMuonSmoother.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "SimDataFormats/Track/interface/SimTrack.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TFile.h>
#include <TTree.h>

typedef struct tagGPSeed {
  GlobalPoint P1;
  GlobalPoint P2;
} GPSeed;


class gemcrValidation : public GEMBaseValidation
{
public:
  explicit gemcrValidation( const edm::ParameterSet& );
  ~gemcrValidation();
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  void analyze(const edm::Event& e, const edm::EventSetup&) override;
  int findIndex(GEMDetId id_);
  int findvfat(float x, float a, float b);
  const GEMGeometry* initGeometry(edm::EventSetup const & iSetup);
  double maxCLS, minCLS,maxRes, trackChi2, trackResY, trackResX, MulSigmaOnWindow;
  bool makeTrack, isMC;

private:
  const GEMGeometry* GEMGeometry_;
  std::vector<MonitorElement*> gem_chamber_x_y;
  std::vector<MonitorElement*> gem_chamber_cl_size;
  std::vector<MonitorElement*> gem_chamber_bx;
  std::vector<MonitorElement*> gem_chamber_pad_vfat;
  std::vector<MonitorElement*> gem_chamber_copad_vfat;
  std::vector<MonitorElement*> gem_chamber_pad_vfat_withmul;
  std::vector<MonitorElement*> gem_chamber_copad_vfat_withmul;
  std::vector<MonitorElement*> gem_chamber_tr2D_eff;
  std::vector<MonitorElement*> gem_chamber_th2D_eff;
  std::vector<MonitorElement*> gem_chamber_trxroll_eff;
  std::vector<MonitorElement*> gem_chamber_thxroll_eff;
  std::vector<MonitorElement*> gem_chamber_trxy_eff;
  std::vector<MonitorElement*> gem_chamber_thxy_eff;
  std::vector<MonitorElement*> gem_chamber_residual;
  std::vector<MonitorElement*> gem_chamber_residualX1DSim;
  std::vector<MonitorElement*> gem_chamber_residualY1DSim;
  std::vector<MonitorElement*> gem_chamber_local_x;
  std::vector<MonitorElement*> gem_chamber_digi_digi;
  std::vector<MonitorElement*> gem_chamber_digi_recHit;
  std::vector<MonitorElement*> gem_chamber_digi_CLS;
  std::vector<MonitorElement*> gem_chamber_hitMul;
  std::vector<MonitorElement*> gem_chamber_vfatHitMul;
  std::vector<MonitorElement*> gem_chamber_stripHitMul;
  std::vector<MonitorElement*> gem_chamber_bestChi2;
  std::vector<MonitorElement*> gem_chamber_track;
  
  std::vector<MonitorElement*> gem_chamber_th2D_eff_scint;
  std::vector<MonitorElement*> gem_chamber_thxroll_eff_scint;
  std::vector<MonitorElement*> gem_chamber_thxy_eff_scint;
  
  std::vector<MonitorElement*> gem_chamber_tr2D_eff_scint;
  std::vector<MonitorElement*> gem_chamber_trxroll_eff_scint;
  std::vector<MonitorElement*> gem_chamber_trxy_eff_scint;
  std::vector<MonitorElement*> gem_chamber_local_x_scint;
  std::vector<MonitorElement*> gem_chamber_residual_scint;
  MonitorElement* rh3_chamber_scint;


  MonitorElement* gemcr_g;
  MonitorElement* gemcrTr_g;
  MonitorElement* gemcrCf_g;
  MonitorElement* gemcrTrScint_g;
  MonitorElement* gemcrCfScint_g;
  MonitorElement* gem_cls_tot;
  MonitorElement* gem_bx_tot;
  MonitorElement* tr_size;
  MonitorElement* tr_hit_size;
  MonitorElement* tr_chamber;
  MonitorElement* th_chamber;
  MonitorElement* rh_chamber;
  MonitorElement* rh1_chamber;
  MonitorElement* rh2_chamber;
  MonitorElement* rh3_chamber;
  MonitorElement* trajectoryh;
  MonitorElement* events_withtraj;
  MonitorElement* firedMul;
  MonitorElement* firedChamber;
  MonitorElement* scinUpperHit;
  MonitorElement* scinLowerHit;
  MonitorElement* scinUpperRecHit;
  MonitorElement* scinLowerRecHit;
  MonitorElement* resXSim;
  MonitorElement* resXByErrSim;
  MonitorElement* resYByErrSim;
  MonitorElement* hitXErr;
  MonitorElement* hitYErr;
  MonitorElement* aftershots;
  MonitorElement* projtheta_dist_sim;
  MonitorElement* projtheta_dist_edge_sim;
  //MonitorElement* diffTrajGenRec;
  

  std::vector<GEMChamber> gemChambers;
  int n_ch;
  MuonServiceProxy* theService;
  CosmicMuonSmoother* theSmoother;
  KFUpdator* theUpdator;
  edm::EDGetToken InputTagToken_, InputTagToken_RH, InputTagToken_TR, InputTagToken_TS, InputTagToken_TJ, InputTagToken_TI, InputTagToken_TT, InputTagToken_DG, InputTagToken_US;
  
  
  bool isPassedScintillators(GlobalPoint p1, GlobalPoint p2);
  float CalcWindowWidthX(GPSeed *pVecSeed, GlobalPoint *pPCurr);
  float CalcWindowWidthY(GPSeed *pVecSeed, GlobalPoint *pPCurr);
  //int CalcDiffGenRec(GPSeed *pVecSeed, GlobalPoint *pPCurr);


  TTree *tree;
  int run;
  int lumi;
  int nev;
  float genMuPt;
  float genMuTheta;
  float genMuPhi;
  float genMuX;
  float genMuY;
  float genMuZ;
  const static int maxNlayer = 30;
  const static int maxNphi = 3;
  const static int maxNeta = 8;
  int vfatI[maxNlayer][maxNphi][maxNeta];
  int vfatF[maxNlayer][maxNphi][maxNeta];

  const static int maxNRecHit = 100;
  int nglobalRecHit;
  float globalRecHitTheta[maxNRecHit];
  float globalRecHitPhi[maxNRecHit];
  float globalRecHitX[maxNRecHit];
  float globalRecHitY[maxNRecHit];
  float globalRecHitZ[maxNRecHit];

  int nTotTrajRecHit;
  float trajRecHitTheta[maxNRecHit];
  float trajRecHitPhi[maxNRecHit];
  float trajRecHitX[maxNRecHit];
  float trajRecHitY[maxNRecHit];
  float trajRecHitZ[maxNRecHit];

};

#endif
