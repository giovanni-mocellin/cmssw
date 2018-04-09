#include "Validation/GEMCR/interface/gemcrValidation.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/GEMGeometry/interface/GEMSuperChamber.h"
#include <Geometry/GEMGeometry/interface/GEMGeometry.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/GeometrySurface/interface/Bounds.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/Math/interface/Vector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Common/interface/RefToBase.h" 
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoMuon/CosmicMuonProducer/interface/HeaderForQC8.h"

#include <iomanip>

#include <TCanvas.h>
using namespace std;
using namespace edm;

gemcrValidation::gemcrValidation(const edm::ParameterSet& cfg): GEMBaseValidation(cfg)
{
  time_t rawTime;
  time(&rawTime);
  printf("Begin of gemcrValidation::gemcrValidation() at %s\n", asctime(localtime(&rawTime)));
  isMC = cfg.getParameter<bool>("isMC");
  InputTagToken_RH = consumes<GEMRecHitCollection>(cfg.getParameter<edm::InputTag>("recHitsInputLabel"));
  InputTagToken_TR = consumes<vector<reco::Track>>(cfg.getParameter<edm::InputTag>("tracksInputLabel"));
  InputTagToken_TS = consumes<vector<TrajectorySeed>>(cfg.getParameter<edm::InputTag>("seedInputLabel"));
  InputTagToken_TJ = consumes<vector<Trajectory>>(cfg.getParameter<edm::InputTag>("trajInputLabel"));
  InputTagToken_TI = consumes<vector<int>>(cfg.getParameter<edm::InputTag>("chNoInputLabel"));
  InputTagToken_TT = consumes<vector<unsigned int>>(cfg.getParameter<edm::InputTag>("seedTypeInputLabel"));
  InputTagToken_DG = consumes<GEMDigiCollection>(cfg.getParameter<edm::InputTag>("gemDigiLabel"));
  if ( isMC ) InputTagToken_US = consumes<edm::HepMCProduct>(cfg.getParameter<edm::InputTag>("genVtx"));
  edm::ParameterSet serviceParameters = cfg.getParameter<edm::ParameterSet>("ServiceParameters");
  theService = new MuonServiceProxy(serviceParameters);
  minCLS = cfg.getParameter<double>("minClusterSize"); 
  maxCLS = cfg.getParameter<double>("maxClusterSize");
  maxRes = cfg.getParameter<double>("maxResidual");
  makeTrack = cfg.getParameter<bool>("makeTrack");
  trackChi2 = cfg.getParameter<double>("trackChi2");
  trackResY = cfg.getParameter<double>("trackResY"); 
  trackResX = cfg.getParameter<double>("trackResX");
  MulSigmaOnWindow = cfg.getParameter<double>("MulSigmaOnWindow");
  SuperChamType = cfg.getParameter<vector<string>>("SuperChamberType");
  vecChamType = cfg.getParameter<vector<double>>("SuperChamberSeedingLayers");
  edm::ParameterSet smootherPSet = cfg.getParameter<edm::ParameterSet>("MuonSmootherParameters");
  theSmoother = new CosmicMuonSmoother(smootherPSet, theService);
  theUpdator = new KFUpdator();
  time(&rawTime);

  //created by Jongseok Lee ==>
  edm::Service<TFileService> fs;
  hev = fs->make<TH1D>("hev","EventSummary",3,1,3);

  genTree = fs->make<TTree>("genTree", "gen info for QC8");
  genTree->Branch("genMuPt",&genMuPt,"genMuPt/F");
  genTree->Branch("genMuTheta",&genMuTheta,"genMuTheta/F");
  genTree->Branch("genMuPhi",&genMuPhi,"genMuPhi/F");
  genTree->Branch("genMuX",&genMuX,"genMuX/F");
  genTree->Branch("genMuY",&genMuY,"genMuY/F");
  genTree->Branch("genMuZ",&genMuZ,"genMuZ/F");
  genTree->Branch("nTraj",&nTraj,"nTraj/I");

  genTree->Branch("nrecHit",&nrecHit,"nrecHit/I");
  genTree->Branch("grecHitX",grecHitX,"grecHitX[nrecHit]/F");
  genTree->Branch("grecHitY",grecHitY,"grecHitY[nrecHit]/F");
  genTree->Branch("grecHitZ",grecHitZ,"grecHitZ[nrecHit]/F");

  hvfatHit_numerator = fs->make<TH3D>("hvfatHit_numerator","vfat hit (numerator of efficiency)",8,0,8,9,0,9,10,0,10);
  hvfatHit_denominator = fs->make<TH3D>("hvfatHit_denominator","vfat hit (denominator of efficienct)",8,0,8,9,0,9,10,0,10);
  tree = fs->make<TTree>("tree", "Tree for QC8");
  tree->Branch("run",&run,"run/I");
  tree->Branch("lumi",&lumi,"lumi/I");
  tree->Branch("ev",&nev,"ev/I");

  tree->Branch("genMuPt",&genMuPt,"genMuPt/F");
  tree->Branch("genMuTheta",&genMuTheta,"genMuTheta/F");
  tree->Branch("genMuPhi",&genMuPhi,"genMuPhi/F");
  tree->Branch("genMuX",&genMuX,"genMuX/F");
  tree->Branch("genMuY",&genMuY,"genMuY/F");
  tree->Branch("genMuZ",&genMuZ,"genMuZ/F");

  ////tree->Branch("nrecHitGP",&nrecHitGP,"nrecHitGP/I");
  //tree->Branch("nrecHit",&nrecHit,"nrecHit/I");
  //tree->Branch("grecHitX",grecHitX,"grecHitX[nrecHit]/F");
  //tree->Branch("grecHitY",grecHitY,"grecHitY[nrecHit]/F");
  //tree->Branch("grecHitZ",grecHitZ,"grecHitZ[nrecHit]/F");

  tree->Branch("trajTheta",&trajTheta,"trajTheta/F");
  tree->Branch("trajPhi",&trajPhi,"trajPhi/F");
  tree->Branch("trajX",&trajX,"trajX/F");
  tree->Branch("trajY",&trajY,"trajY/F");
  tree->Branch("trajZ",&trajZ,"trajZ/F");
  tree->Branch("trajPx",&trajPx,"trajPx/F");
  tree->Branch("trajPy",&trajPy,"trajPy/F");
  tree->Branch("trajPz",&trajPz,"trajPz/F");
  tree->Branch("nTrajHit",&ntrajHit,"nTrajHit/I");
  tree->Branch("nTrajRecHit",&ntrajRecHit,"nTrajRecHit/I");


  tree->Branch("vfatI",vfatI,"vfatI[30][3][8]/I");
  tree->Branch("vfatF",vfatF,"vfatF[30][3][8]/I");

  tree->Branch("trajHitX",trajHitX,"trajHitX[30][3][8]/F");
  tree->Branch("trajHitY",trajHitY,"trajHitY[30][3][8]/F");
  tree->Branch("trajHitZ",trajHitZ,"trajHitZ[30][3][8]/F");
  tree->Branch("recHitX",recHitX,"recHitX[30][3][8]/F");
  tree->Branch("recHitY",recHitY,"recHitY[30][3][8]/F");
  tree->Branch("recHitZ",recHitZ,"recHitZ[30][3][8]/F");
  tree->Branch("genHitX",genHitX,"genHitX[30][3][8]/F");
  tree->Branch("genHitY",genHitY,"genHitY[30][3][8]/F");
  tree->Branch("genHitZ",genHitZ,"genHitZ[30][3][8]/F");

  tree->Branch("floorHitX",floorHitX,"floorHitX[10]/F");
  tree->Branch("floorHitY",floorHitY,"floorHitY[10]/F");
  tree->Branch("floorHitZ",floorHitZ,"floorHitZ[10]/F");
  //created by Jongseok Lee <==

  printf("End of gemcrValidation::gemcrValidation() at %s\n", asctime(localtime(&rawTime)));
}

MonitorElement* g_resXRTSim;
MonitorElement* g_resYRTByErrSim;

void gemcrValidation::bookHistograms(DQMStore::IBooker & ibooker, edm::Run const & Run, edm::EventSetup const & iSetup ) {
  time_t rawTime;
  time(&rawTime);
  printf("Begin of gemcrValidation::bookHistograms() at %s\n", asctime(localtime(&rawTime)));
  GEMGeometry_ = initGeometry(iSetup);
  if ( GEMGeometry_ == nullptr) return ;  

  const std::vector<const GEMSuperChamber*>& superChambers_ = GEMGeometry_->superChambers();   
  for (auto sch : superChambers_)
  {
    int n_lay = sch->nChambers();
    for (int l=0;l<n_lay;l++)
    {
   	  gemChambers.push_back(*sch->chamber(l+1));
    }
  }
  n_ch = gemChambers.size();
  
  ibooker.setCurrentFolder("MuonGEMRecHitsV/GEMRecHitsTask");

  gemcr_g = ibooker.book3D("gemcr_g","GEMCR GLOBAL RECHITS", 200,-100,100,156,-61,95,167,-12,155);
  gemcrGen_g = ibooker.book3D("gemcrGen_g","GEMCR GLOBAL GEN HITS", 200,-100,100,156,-61,95,167,-12,155);
  gemcrTr_g = ibooker.book3D("gemcrTr_g","GEMCR GLOBAL RECHITS", 200,-100,100,156,-61,95,167,-12,155);
  gemcrCf_g = ibooker.book3D("gemcrCf_g","GEMCR GLOBAL RECHITS CONFIRMED", 200,-100,100,156,-61,95,167,-12,155);
  gemcrTrScint_g = ibooker.book3D("gemcrTrScint_g","GEMCR GLOBAL RECHITS", 200,-100,100,156,-61,95,167,-12,155);
  gemcrCfScint_g = ibooker.book3D("gemcrCfScint_g","GEMCR GLOBAL RECHITS SCINTILLATED", 200,-100,100,156,-61,95,167,-12,155);
  gem_cls_tot = ibooker.book1D("cluster_size","Cluseter size",20,0,20);  
  gem_bx_tot = ibooker.book1D("bx", "BX" , 30, -15,15);
  tr_size = ibooker.book1D("tr_size", "track size",10,0,10);
  tr_hit_size = ibooker.book1D("tr_hit_size", "hit size in track",15,0,15); 
  trajectoryh = ibooker.book1D("trajectory","trajectory", 4,0,4);
  trajectoryh->setBinLabel(1, "total seeds");
  trajectoryh->setBinLabel(2, "unvalid");
  trajectoryh->setBinLabel(3, "valid");
  trajectoryh->setBinLabel(4, "passed chi2");
  firedMul = ibooker.book1D("firedMul","fired chamber multiplicity",n_ch+1,0,n_ch+1);
  firedChamber = ibooker.book1D("firedChamber", "fired chamber",n_ch,0,n_ch);
  scinUpperHit = ibooker.book3D("scinUpperHit","UPPER SCINTILLATOR GLOBAL", 200,-100,100,156,-61,95,167,-12,155);
  scinLowerHit = ibooker.book3D("scinUpperHit","LOWER SCINTILLATOR GLOBAL", 200,-100,100,156,-61,95,167,-12,155);
  scinUpperRecHit = ibooker.book3D("scinLowerRecHit","UPPER SCINTILLATOR GLOBAL RECHITS", 200,-100,100,156,-61,95,167,-12,155);
  scinLowerRecHit = ibooker.book3D("scinLowerRecHit","LOWER SCINTILLATOR GLOBAL RECHITS", 200,-100,100,156,-61,95,167,-12,155);
  
  resXSim = ibooker.book1D("residualx_sim", " residual x (sim)",200,-3,3);
  resXByErrSim = ibooker.book1D("residualxbyeff_sim", " residual x / x_err (sim)",200,-10,10);
  resYByErrSim = ibooker.book1D("residualy_sim", " residual y / y_err (sim)",200,-3,3);
  hitXErr = ibooker.book1D("x_err", "x_err",200,0,5);
  hitYErr = ibooker.book1D("y_err", "y_err",200,0,50);
  
  g_resXRTSim = ibooker.book1D("residualx_reco_track", " residual x (reco traj vs track)",200,-3,3);
  g_resYRTByErrSim = ibooker.book1D("residualy_reco_track", " residual y (reco traj vs track)",200,-3,3);

  tr_chamber = ibooker.book1D("tr_eff_ch", "tr rec /chamber",n_ch,1,n_ch); 
  th_chamber = ibooker.book1D("th_eff_ch", "tr hit/chamber",n_ch,0,n_ch); 
  rh_chamber = ibooker.book1D("rh_eff_ch", "rec hit/chamber",n_ch,0,n_ch); 
  rh1_chamber = ibooker.book1D("rh1_chamber", "all recHits",n_ch,0,n_ch); 
  rh2_chamber = ibooker.book1D("rh2_chamber", "cut passed recHits ",n_ch,0,n_ch); 
  rh3_chamber = ibooker.book1D("rh3_chamber", "tracking recHits",n_ch,0,n_ch); 
  
  rh3_chamber_scint = ibooker.book1D("rh3_chamber_scint", "tracking recHits, scintillated",n_ch,0,n_ch); 
  
  events_withtraj = ibooker.book1D("events_withtraj", "# of events with a reconstructed trajectory",8,0,8); 
  
  events_withtraj->setBinLabel(1, "events with firedCh >= 4");
  events_withtraj->setBinLabel(2, "events with firedCh >= 5");
  events_withtraj->setBinLabel(3, "events with firedCh >= 6");
  events_withtraj->setBinLabel(4, "events with top-bottom seed");
  events_withtraj->setBinLabel(5, "events with tr, firedCh >= 4");
  events_withtraj->setBinLabel(6, "events with tr, firedCh >= 5");
  events_withtraj->setBinLabel(7, "events with tr, firedCh >= 6");
  events_withtraj->setBinLabel(8, "events with tr, top-bottom seed");
  
  aftershots = ibooker.book2D("plain_aftershots", "aftershots", 25, -130, 130, 25, -82.5, 82.5); 
  
  projtheta_dist_sim = ibooker.book1D("projtheta_dist_sim", "theta distribution (SIM)", 300, 0, 15);
  projtheta_dist_edge_sim = ibooker.book1D("projtheta_dist_edge_sim", "theta distribution (SIM)", 300, 0, 15);
  
  
  for(int c = 0; c<n_ch; c++){
   //cout << gemChambers[c].id() << endl;
   GEMDetId gid = gemChambers[c].id();
   string b_name = "chamber_"+to_string(gid.chamber())+"_layer_"+to_string(gid.layer());
   tr_chamber->setBinLabel(c+1,b_name);
   th_chamber->setBinLabel(c+1,b_name);
   rh_chamber->setBinLabel(c+1,b_name);
   rh1_chamber->setBinLabel(c+1,b_name);
   rh2_chamber->setBinLabel(c+1,b_name);
   rh3_chamber->setBinLabel(c+1,b_name);
   firedChamber->setBinLabel(c+1,b_name);
   rh3_chamber_scint->setBinLabel(c+1,b_name);
  }
  for(int c = 0; c<n_ch;c++){
     GEMDetId gid = gemChambers[c].id();
     string h_name = "chamber_"+to_string(gid.chamber())+"_layer_"+to_string(gid.layer());
     gem_chamber_x_y.push_back(ibooker.book2D(h_name+"_recHit",h_name+" recHit", 500,-25,25,8,1,9));
     gem_chamber_cl_size.push_back(ibooker.book2D(h_name+"_recHit_size", h_name+" recHit size", 50,0,50,24,0,24));
     gem_chamber_bx.push_back(ibooker.book2D(h_name+"_bx", h_name+" BX", 30,-15,15,10,0,10));
     gem_chamber_pad_vfat.push_back(ibooker.book2D(h_name+"_pad_vfat", h_name+" pad (vfat)", 3,1,4,8,1,9));
     gem_chamber_copad_vfat.push_back(ibooker.book2D(h_name+"_copad_vfat", h_name+" copad (vfat)", 3,1,4,8,1,9));
     gem_chamber_pad_vfat_withmul.push_back(ibooker.book2D(h_name+"_pad_vfat_withmul", h_name+" pad (vfat), with multiplicity", 3,1,4,8,1,9));
     gem_chamber_copad_vfat_withmul.push_back(ibooker.book2D(h_name+"_copad_vfat_withmul", h_name+" copad (vfat), with multiplicity", 3,1,4,8,1,9));
     gem_chamber_tr2D_eff.push_back(ibooker.book2D(h_name+"_recHit_efficiency", h_name+" recHit efficiency", 3,1,4,8,1,9));
     gem_chamber_th2D_eff.push_back(ibooker.book2D(h_name+"_th2D_eff", h_name+"_th2D_eff", 3,1,4,8,1,9));
     gem_chamber_trxroll_eff.push_back(ibooker.book2D(h_name+"_trxroll_eff", h_name+" recHit efficiency", 50,-25,25,8,1,9));
     gem_chamber_trxy_eff.push_back(ibooker.book2D(h_name+"_trxy_eff", h_name+" recHit efficiency", 50,-25,25,120,-60,60));
     gem_chamber_thxroll_eff.push_back(ibooker.book2D(h_name+"_thxroll_eff", h_name+"_th2D_eff", 50,-25,25,8,1,9));
     gem_chamber_thxy_eff.push_back(ibooker.book2D(h_name+"_thxy_eff", h_name+"_th2D_eff", 50,-25,25,120,-60,60));
     gem_chamber_residual.push_back(ibooker.book2D(h_name+"_residual", h_name+" residual", 500,-25,25,100,-5,5));
     gem_chamber_residualX1DSim.push_back(ibooker.book1D(h_name+"_residualX1DSim", h_name+" residual x 1D (sim)", 200,-3,3));
     gem_chamber_residualY1DSim.push_back(ibooker.book1D(h_name+"_residualY1DSim", h_name+" residual y 1D (sim)", 200,-3,3));
     gem_chamber_local_x.push_back(ibooker.book2D(h_name+"_local_x", h_name+" local x",500,-25,25,500,-25,25));
     gem_chamber_digi_digi.push_back(ibooker.book2D(h_name+"_digi_gemDigi", h_name+" gemDigi (DIGI)", 384,0,384,8,1,9));
     gem_chamber_digi_recHit.push_back(ibooker.book2D(h_name+"_recHit_gemDigi", h_name+" gemDigi (recHit)", 384,0,384,8,1,9));
     gem_chamber_digi_CLS.push_back(ibooker.book2D(h_name+"_CLS_gemDigi", h_name+" gemDigi (CLS)", 384,0,384,8,1,9));
     gem_chamber_hitMul.push_back(ibooker.book1D(h_name+"_hit_mul", h_name+" hit multiplicity",25,0,25 ));
     gem_chamber_vfatHitMul.push_back(ibooker.book2D(h_name+"_vfatHit_mul", h_name+" vfat hit multiplicity",25,0,25, 24,0,24));
     gem_chamber_stripHitMul.push_back(ibooker.book2D(h_name+"_stripHit_mul", h_name+" strip hit multiplicity", 150,0,150,9,0,9));
     gem_chamber_bestChi2.push_back(ibooker.book1D(h_name+"_bestChi2", h_name+" #chi^{2} distribution", trackChi2*10,0,trackChi2));
     gem_chamber_track.push_back(ibooker.book1D(h_name+"_track", h_name+" track",7,0,7));
     
     gem_chamber_th2D_eff_scint.push_back(ibooker.book2D(h_name+"_th2D_eff_scint", h_name+"_th2D_eff, scintillated", 3,1,4,8,1,9));
     gem_chamber_thxroll_eff_scint.push_back(ibooker.book2D(h_name+"_thxroll_eff_scint", h_name+"_th2D_eff, scintillated", 50,-25,25,8,1,9));
     gem_chamber_thxy_eff_scint.push_back(ibooker.book2D(h_name+"_thxy_eff_scint", h_name+"_th2D_eff, scintillated", 50,-25,25,120,-60,60));
     
     gem_chamber_tr2D_eff_scint.push_back(ibooker.book2D(h_name+"_recHit_efficiency_scint", h_name+" recHit efficiency, scintillated", 3,1,4,8,1,9));
     gem_chamber_trxroll_eff_scint.push_back(ibooker.book2D(h_name+"_trxroll_eff_scint", h_name+" recHit efficiency, scintillated", 50,-25,25,8,1,9));
     gem_chamber_trxy_eff_scint.push_back(ibooker.book2D(h_name+"_trxy_eff_scint", h_name+" recHit efficiency, scintillated", 50,-25,25,120,-60,60));
     gem_chamber_local_x_scint.push_back(ibooker.book2D(h_name+"_local_x_scint", h_name+" local x, scintillated",500,-25,25,500,-25,25));
     gem_chamber_residual_scint.push_back(ibooker.book2D(h_name+"_residual_scint", h_name+" residual, scintillated", 500,-25,25,100,-5,5));
  }
  time(&rawTime);
  printf("End of gemcrValidation::bookHistograms() at %s\n", asctime(localtime(&rawTime)));
}

int gemcrValidation::findIndex(GEMDetId id_, bool bIsFindCopad = false) {
  int index=-1;
  for(int c =0;c<n_ch;c++){
    if((gemChambers[c].id().chamber() == id_.chamber())&&(bIsFindCopad ^ ( gemChambers[c].id().layer() == id_.layer() )) ){index = c;}
  }
  return index;
}

int gemcrValidation::findvfat(float x, float a, float b) {
  float step = abs(b-a)/3.0;
  if ( x < (min(a,b)+step) ) { return 1;}
  else if ( x < (min(a,b)+2.0*step) ) { return 2;}
  else { return 3;}
}

const GEMGeometry* gemcrValidation::initGeometry(edm::EventSetup const & iSetup) {
  const GEMGeometry* GEMGeometry_ = nullptr;
  try {
    edm::ESHandle<GEMGeometry> hGeom;
    iSetup.get<MuonGeometryRecord>().get(hGeom);
    GEMGeometry_ = &*hGeom;
  }
  catch( edm::eventsetup::NoProxyException<GEMGeometry>& e) {
    edm::LogError("MuonGEMBaseValidation") << "+++ Error : GEM geometry is unavailable on event loop. +++\n";
    return nullptr;
  }
  return GEMGeometry_;
}

////////////////////////////////////////////////////////////////////////////////
// To turn on or turn off methods for bad trajectory veto, seek them: 
// CONE_WINDOW
// VETO_TOO_FEW_RECHITS
// VETO_TOO_SHORT_SEED
////////////////////////////////////////////////////////////////////////////////

int g_nEvt = 0;
int g_nNumRecHit = 0;
int g_nNumFiredCh = 0;
int g_nNumFiredChValid = 0;
int g_nNumTrajHit = 0;
int g_nNumTrajHit6 = 0;
int g_nNumMatched = 0;

double g_dMinX = 100000.0, g_dMaxX = -10000000.0;
double g_dMinY = 100000.0, g_dMaxY = -10000000.0;
double g_dMinZ = 100000.0, g_dMaxZ = -10000000.0;

gemcrValidation::~gemcrValidation() {
  printf("res1 : %i\n", g_nEvt);
  printf("res2 : %i\n", g_nNumRecHit);
  printf("res3 : %i\n", g_nNumFiredCh);
  printf("res4 : %i\n", g_nNumFiredChValid);
  printf("res5 : %i\n", g_nNumTrajHit6);
  printf("res6 : %i (g_nNumTrajHit)\n", g_nNumTrajHit);
  printf("res7 : %i (g_nNumMatched)\n", g_nNumMatched);
  if(g_nNumTrajHit>0) printf("eff : %f\n", double(g_nNumMatched)/g_nNumTrajHit);
  
  printf("X range : %lf, %lf\n", g_dMinX, g_dMaxX);
  printf("Y range : %lf, %lf\n", g_dMinY, g_dMaxY);
  printf("Z range : %lf, %lf\n", g_dMinZ, g_dMaxZ);
}


bool gemcrValidation::isPassedScintillators(GlobalPoint p1, GlobalPoint p2)
{
  bool isPassedScint = false;
  
  double ScintilLowerZ =  -11.485;
  double ScintilUpperZ =  154.015;
  double ScintilXMin   = -100.0;
  double ScintilXMax   =  100.0;
  double ScintilYMin   =  -60.56;
  double ScintilYMax   =   63.0;
  //https://brilliant.org/wiki/3d-coordinate-geometry-equation-of-a-line
  //3D line equation for genMuon (generated track) : (x - genMuX)/gen_px = (y - genMuY)/gen_py = (z - genMuZ)/gen_pz;
  //if y=y0 then x = genMuX + (y0-genMuY)*(gen_px/gen_py),  z = genMuZ + (y0-genMuY)*(gen_pz/gen_py)
  //if z=z0 then x = genMuX + (z0-genMuZ)*(gen_px/gen_pz),  y = genMuY + (z0-genMuZ)*(gen_py/gen_pz)
  
  // Finding the X coordinate of the track in the upper and lower scintillator planes : (x)=(m_xy)*(y)+(q_xy), (m_xy)=(x2-x1)/(y2-y1), (q_xy)=x2-(m_xy)*y2
  double m_xz = ( p2.x() - p1.x() ) / ( p2.z() - p1.z() ) ;
  double q_xz = p2.x() - m_xz * p2.z() ;
  double HitXupScint  = m_xz * ScintilUpperZ + q_xz ;
  double HitXlowScint = m_xz * ScintilLowerZ + q_xz ;
  //double HitXupScint  = p1.x() + (ScintilUpperZ - p1.z())*m_xz;
  //double HitXlowScint = p1.x() + (ScintilLowerZ - p1.z())*m_xz;
  
  // Finding the Z coordinate of the track in the upper and lower scintillator planes : (z)=(m_zy)*(y)+(q_zy), (m_zy)=(z2-z1)/(y2-y1), (q_zy)=z2-(m_zy)*y2
  double m_yz = ( p2.y() - p1.y() ) / ( p2.z() - p1.z() ) ;
  double q_yz = p2.y() - m_yz * p2.z() ;
  double HitYupScint  = m_yz * ScintilUpperZ + q_yz ;
  double HitYlowScint = m_yz * ScintilLowerZ + q_yz ;
  
  if (( ScintilXMin <= HitXupScint  && HitXupScint  <= ScintilXMax ) && ( ScintilYMin  <= HitYupScint  && HitYupScint  <= ScintilYMax ) &&
      ( ScintilXMin <= HitXlowScint && HitXlowScint <= ScintilXMax ) && ( ScintilYMin  <= HitYlowScint && HitYlowScint <= ScintilYMax ) )
  {
    isPassedScint = true;
  }
  
  scinUpperHit -> Fill(HitXupScint, HitYupScint, ScintilUpperZ);
  scinLowerHit -> Fill(HitXlowScint, HitYlowScint, ScintilLowerZ);
  
  if  (isPassedScint == true)
  {
      scinUpperRecHit->Fill(HitXupScint, HitYupScint, ScintilUpperZ);
      scinLowerRecHit->Fill(HitXlowScint, HitYlowScint, ScintilLowerZ);
  }
  
  return isPassedScint;
}


int g_nNumTest = 0;


void gemcrValidation::analyze(const edm::Event& e, const edm::EventSetup& iSetup){
  g_nEvt++;

  //created by Jongseok Lee ==>
  run = e.id().run();
  lumi = e.id().luminosityBlock();
  nev = e.id().event();
  //if(nev%1000==0) cout<<"nev "<<nev<<endl;
  //cout<<"nev "<<nev<<endl;
  hev->Fill(1);

  genMuPt = -10;
  genMuTheta = -10;
  genMuPhi = -10;
  genMuX = 0;
  genMuY = 0;
  genMuZ = 0;

  nTraj = 0;

  for(int i=0;i<maxNlayer;i++)
  {
    for(int j=0;j<maxNphi;j++)
    {
      for(int k=0;k<maxNeta;k++)
      {
        vfatI[i][j][k] = 0;
        vfatF[i][j][k] = 0;
        trajHitX[i][j][k] = 0;
        trajHitY[i][j][k] = 0;
        trajHitZ[i][j][k] = 0;
        recHitX[i][j][k] = 0;
        recHitY[i][j][k] = 0;
        recHitZ[i][j][k] = 0;
        genHitX[i][j][k] = 0;
        genHitY[i][j][k] = 0;
        genHitZ[i][j][k] = 0;
      }
    }
  }
  trajTheta = -10;
  trajPhi = -10;
  trajX = 0;
  trajY = 0;
  trajZ = 0;
  trajPx = 0;
  trajPy = 0;
  trajPz = 0;
  ntrajHit = 0;
  ntrajRecHit = 0;

  for(int i=0;i<maxNfloor;i++)
  {
    floorHitX[i] = 0;
    floorHitY[i] = 0;
    floorHitZ[i] = 0;
  }

  nrecHit = 0;
  //created by Jongseok Lee <==

  theService->update(iSetup);

  edm::Handle<GEMRecHitCollection> gemRecHits;
  if (!isMC)
  {
    edm::Handle<GEMDigiCollection> digis;
    e.getByToken( this->InputTagToken_DG, digis);
    int nNumDigi = 0;
    for (GEMDigiCollection::DigiRangeIterator gemdgIt = digis->begin(); gemdgIt != digis->end(); ++gemdgIt)
    {
      nNumDigi++;
      const GEMDetId& gemId = (*gemdgIt).first;
      int index = findIndex(gemId);
      const GEMDigiCollection::Range& range = (*gemdgIt).second;
      for ( auto digi = range.first; digi != range.second; ++digi )
      {
        gem_chamber_digi_digi[index]->Fill(digi->strip(),gemId.roll());
      }
    }
  }
  e.getByToken( this->InputTagToken_RH, gemRecHits);
  if (!gemRecHits.isValid())
  {
    edm::LogError("gemcrValidation") << "Cannot get strips by Token RecHits Token.\n";
    return ;
  }
  
  if ( g_nEvt == -1 )
  {
    int nIdxChTmp = 0;
    for (auto tch : gemChambers)
    {
      nIdxChTmp++;
      int nIdxEP = 0;
      
      for (auto etaPart : tch.etaPartitions())
      {
        nIdxEP++;
        printf("(%2i, %i) : %i (%X) - ", nIdxChTmp, nIdxEP, etaPart->id().rawId(), etaPart->id().rawId());
        std::cout << etaPart->id() << std::endl;
      }
    }
    
    double dXTest = gemChambers[ 0 ].etaPartition(8)->centreOfStrip(0).x() + 0.01;
    Local3DPoint tlproll1(dXTest, 0.0, 0.0);
    
    printf("Interval between two strips : %lf\n", gemChambers[ 0 ].etaPartition(8)->centreOfStrip(1).x() - gemChambers[ 0 ].etaPartition(8)->centreOfStrip(0).x()); fflush(stdout);
        
    for ( int i = 1 ; i <= 8 ; i++ ) {
      const BoundPlane& bprollCurr = GEMGeometry_->idToDet(gemChambers[ 0 ].etaPartition(i)->id())->surface();
      printf("Roll area %i test : %s ; %lf, %lf, %i ; %lf\n", i, 
        ( bprollCurr.bounds().inside(tlproll1) ? "inside" : "outside" ), 
        gemChambers[ 0 ].etaPartition(i)->centreOfStrip(1).x(), 
        gemChambers[ 0 ].etaPartition(i)->centreOfStrip(gemChambers[ 0 ].etaPartition(i)->nstrips()).x(), 
        gemChambers[ 0 ].etaPartition(i)->nstrips(), 
        ( i <= 1 ? 0.0 : 
          bprollCurr.toLocal(GEMGeometry_->idToDet(gemChambers[ 0 ].etaPartition(i-1)->id())->surface().position()).y() ));
    }
    
    for ( int i = 1 ; i <= 8 ; i++ ) {
      const BoundPlane& bprollCurr = GEMGeometry_->idToDet(gemChambers[ 0 ].etaPartition(i)->id())->surface();
      printf("Roll area %i dimension : %lf, %lf, %lf ; (%lf, %lf, %lf)\n", i, 
        bprollCurr.bounds().length(), 
        bprollCurr.bounds().width(), 
        bprollCurr.bounds().thickness(), 
        bprollCurr.position().x(), 
        bprollCurr.position().z(), 
        bprollCurr.position().y());
    }
    
    for ( int i = 10 ; i <= 19 ; i++ ) {
      const BoundPlane& bprollCurr = GEMGeometry_->idToDet(gemChambers[ i ].etaPartition(1)->id())->surface();
      printf("Chamber %i Roll area 0 dimension : %lf, %lf, %lf ; (%lf, %lf, %lf)\n", i, 
        bprollCurr.bounds().length(), 
        bprollCurr.bounds().width(), 
        bprollCurr.bounds().thickness(), 
        bprollCurr.position().x(), 
        bprollCurr.position().z(), 
        bprollCurr.position().y());
    }
    
    for ( int i = -2500 ; i <= 2500 ; i++ ) {
      Local3DPoint tltest1(i * 0.01, -9.69, 0);
      Local3DPoint tltest2(i * 0.01,  9.69, 0);
      printf("x test - %5.2f : %s (-), %s (+)\n", 0.01 * i, 
        GEMGeometry_->idToDet(gemChambers[ 0 ].etaPartition(1)->id())->surface().bounds().inside(tltest1) ? "in " : "out", 
        GEMGeometry_->idToDet(gemChambers[ 0 ].etaPartition(1)->id())->surface().bounds().inside(tltest2) ? "in " : "out");
    }
    
    for ( int i = -2 ; i <= 50 ; i++ ) {
      Local3DPoint tltestJump(0, -9.71 - 0.01 * i, 0);
      Global3DPoint gptestJump = GEMGeometry_->idToDet(gemChambers[ 0 ].etaPartition(1)->id())->surface().toGlobal(tltestJump);
      Local3DPoint tltestJump2 = GEMGeometry_->idToDet(gemChambers[ 0 ].etaPartition(2)->id())->surface().toLocal(gptestJump);
      printf("jump test - %5.2f : %s\n", -9.71 - 0.01 * i, 
        ( GEMGeometry_->idToDet(gemChambers[ 0 ].etaPartition(2)->id())->surface().bounds().inside(tltestJump2) ? "in " : "out" ));
    }
  }
  
  TString strKeep("");
  
  float fXGenGP1x = 0.0, fXGenGP1y = 0.0, fXGenGP1z = 0.0;
  float fXGenGP2x = 0.0, fXGenGP2y = 0.0, fXGenGP2z = 0.0;
  
  int nNumCurrFiredCh = 0;
  HepMC::GenParticle *genMuon = NULL;
  
  //created by Jongseok Lee ==>
  double gen_px = 0;
  double gen_py = 0;
  double gen_pz = 0;
  double gen_pt = 0;
  double gen_theta = 0;
  double gen_phi = 0;

  //nrecHit = gemRecHits->size();
  //created by Jongseok Lee <==

  if ( isMC )
  {
    edm::Handle<edm::HepMCProduct> genVtx;
    e.getByToken( this->InputTagToken_US, genVtx);
    genMuon = genVtx->GetEvent()->barcode_to_particle(1);
    
    double dUnitGen = 0.1;
    
    fXGenGP1x = dUnitGen * genMuon->production_vertex()->position().x();
    fXGenGP1y = dUnitGen * genMuon->production_vertex()->position().y();
    fXGenGP1z = dUnitGen * genMuon->production_vertex()->position().z();
    
    fXGenGP2x = fXGenGP1x + dUnitGen * genMuon->momentum().x();
    fXGenGP2y = fXGenGP1y + dUnitGen * genMuon->momentum().y();
    fXGenGP2z = fXGenGP1z + dUnitGen * genMuon->momentum().z();
    
    Float_t fVecX, fVecY;
    int arrnFired[ 32 ] = {0, };
    
    fVecX = genMuon->momentum().x() / genMuon->momentum().z();
    fVecY = genMuon->momentum().y() / genMuon->momentum().z();

    //created by Jongseok Lee ==>
    gen_px = genMuon->momentum().x();
    gen_py = genMuon->momentum().y();
    gen_pz = genMuon->momentum().z();
    gen_pt = genMuon->momentum().perp();
    gen_theta = genMuon->momentum().theta();
    gen_phi = genMuon->momentum().phi();

    genMuPt = gen_pt;
    genMuTheta = gen_theta;
    genMuPhi = gen_phi;

    genMuX = fXGenGP1x;
    genMuY = fXGenGP1y;
    genMuZ = fXGenGP1z;

    //https://brilliant.org/wiki/3d-coordinate-geometry-equation-of-a-line
    //3D line equation for genMuon (generated track) : (x - genMuX)/gen_px = (y - genMuY)/gen_py = (z - genMuZ)/gen_pz;
    //if y=y0 then x = genMuX + (y0-genMuY)*(gen_px/gen_py),  z = genMuZ + (y0-genMuY)*(gen_pz/gen_py)
    //if z=z0 then x = genMuX + (z0-genMuZ)*(gen_px/gen_pz),  y = genMuY + (z0-genMuZ)*(gen_py/gen_pz)
    for ( GEMRecHitCollection::const_iterator rechit = gemRecHits->begin(); rechit != gemRecHits->end(); ++rechit )
    {
      //int nIdxCh = findIndex((*rechit).gemId());
      if ((*rechit).clusterSize()<minCLS) continue;
      if ((*rechit).clusterSize()>maxCLS) continue;

      GlobalPoint recHitGP = GEMGeometry_->idToDet((*rechit).gemId())->surface().toGlobal(rechit->localPosition());
      grecHitX[nrecHit] = recHitGP.x();
      grecHitY[nrecHit] = recHitGP.y();
      grecHitZ[nrecHit] = recHitGP.z();
      nrecHit++;
    }
    //created by Jongseok Lee <==

    for ( GEMRecHitCollection::const_iterator recHit = gemRecHits->begin(); recHit != gemRecHits->end(); ++recHit )
    {
      int nIdxCh = findIndex((*recHit).gemId());
      GlobalPoint recHitGP = GEMGeometry_->idToDet((*recHit).gemId())->surface().toGlobal(recHit->localPosition());
      
      if ( arrnFired[ nIdxCh ] == 0 )
      {
        arrnFired[ nIdxCh ] = 1;
        nNumCurrFiredCh++;
      }
      
      Float_t fDiffZ = recHitGP.z() - fXGenGP1z;
      
      Float_t fXGenHitX = fXGenGP1x + fDiffZ * fVecX;
      Float_t fXGenHitY = fXGenGP1y + fDiffZ * fVecY;

      strKeep += TString::Format("  recHit : %i, RECO : (%0.5f, %0.5f, %0.5f) <... GEN : (%0.5f, %0.5f, %0.5f)\n", nIdxCh + 1, recHitGP.x(), recHitGP.z(), recHitGP.y(), fXGenHitX, fXGenHitY, recHitGP.z());
      
      resYByErrSim->Fill(( recHitGP.y() - fXGenHitY ) / recHit->localPositionError().yy());
      gem_chamber_residualY1DSim[ nIdxCh ]->Fill(( recHitGP.y() - fXGenHitY ) / recHit->localPositionError().yy());
      hitYErr->Fill(recHit->localPositionError().yy());
      
      g_nNumRecHit++;
    }
    
    int nIdxCh = 0;
    for (auto tch : gemChambers){
      nIdxCh++;

      for (auto etaPart : tch.etaPartitions()){
        GEMDetId etaPartID = etaPart->id();
        GEMRecHitCollection::range range = gemRecHits->get(etaPartID);
        for (GEMRecHitCollection::const_iterator rechit = range.first; rechit!=range.second; ++rechit){
            int nContinueCLS = 0;
            if ((*rechit).clusterSize()<minCLS) nContinueCLS |= 0x1;
            if ((*rechit).clusterSize()>maxCLS) nContinueCLS |= 0x2;
            
            GlobalPoint recHitGP = GEMGeometry_->idToDet((*rechit).gemId())->surface().toGlobal(rechit->localPosition());
            
            Float_t fDiffZ = recHitGP.z() - fXGenGP1z;
            
            Float_t fXGenHitX = fXGenGP1x + fDiffZ * fVecX;
            Float_t fXGenHitY = fXGenGP1y + fDiffZ * fVecY;
            
            resXSim->Fill(recHitGP.x() - fXGenHitX);
            resXByErrSim->Fill(( recHitGP.x() - fXGenHitX ) / rechit->localPositionError().xx());
            gem_chamber_residualX1DSim[ nIdxCh - 1 ]->Fill(recHitGP.x() - fXGenHitX);
            hitXErr->Fill(rechit->localPositionError().xx());
            
            resYByErrSim->Fill(( recHitGP.y() - fXGenHitY ) / rechit->localPositionError().yy());
            gem_chamber_residualY1DSim[ nIdxCh - 1 ]->Fill(( recHitGP.y() - fXGenHitY ) / rechit->localPositionError().yy());
            hitYErr->Fill(rechit->localPositionError().yy());
            
            strKeep += TString::Format("  recHit (TCH) : %i (%i), RECO : (%0.5f, %0.5f, %0.5f) <... GEN : (%0.5f, %0.5f, %0.5f)\n", 
              nIdxCh, nContinueCLS, recHitGP.x(), recHitGP.y(), recHitGP.z(), fXGenHitX, fXGenHitY, recHitGP.z());
            
            strKeep += TString::Format("                 LOCAL : (%0.5f, %0.5f, %0.5f) ... (%0.5f, %0.5f, %0.5f)\n",
              tch.surface().position().x(), tch.surface().position().y(), tch.surface().position().z(), 
              GEMGeometry_->idToDet((*rechit).gemId())->surface().position().x(), 
              GEMGeometry_->idToDet((*rechit).gemId())->surface().position().z(), 
              GEMGeometry_->idToDet((*rechit).gemId())->surface().position().y());
        }
      }
    }
    
    aftershots->Fill(fXGenGP1x + 400.0 * fVecX, fXGenGP1y + 400.0 * fVecY);
    
    g_nNumFiredCh += nNumCurrFiredCh;
    if ( nNumCurrFiredCh > 6 ) g_nNumFiredChValid += nNumCurrFiredCh;
    
    projtheta_dist_sim->Fill(180.0 / 3.141592 * atan2(genMuon->momentum().z(), genMuon->momentum().x()));
  }
  
  GlobalPoint genGPos1(fXGenGP1x, fXGenGP1y, fXGenGP1z);
  GlobalPoint genGPos2(fXGenGP2x, fXGenGP2y, fXGenGP2z);
  
  vector<bool> firedCh;
  vector<int> rMul;
  vector<vector<int>> vMul(n_ch, vector<int>(24, 0));
  vector<vector<int>> sMul(n_ch, vector<int>(9, 0));
  for (int c=0;c<n_ch;c++){
    firedCh.push_back(0);
    rMul.push_back(0);
  }
  TString strListRecHit("");
  for (GEMRecHitCollection::const_iterator recHit = gemRecHits->begin(); recHit != gemRecHits->end(); ++recHit){

    Float_t  rh_l_x = recHit->localPosition().x();
    Int_t  bx = recHit->BunchX();
    Int_t  clusterSize = recHit->clusterSize();
    Int_t  firstClusterStrip = recHit->firstClusterStrip();

    GEMDetId id((*recHit).gemId());
    int index = findIndex(id);
    firedCh[index] = 1;
    rMul[index] += 1;
    Short_t rh_roll = (Short_t) id.roll();
    LocalPoint recHitLP = recHit->localPosition();

    if ( GEMGeometry_->idToDet((*recHit).gemId()) == nullptr) {
      std::cout<<"This gem recHit did not matched with GEMGeometry."<<std::endl;
      continue;
    }
    GlobalPoint recHitGP = GEMGeometry_->idToDet((*recHit).gemId())->surface().toGlobal(recHitLP);
    Float_t     rh_g_X = recHitGP.x();
    Float_t     rh_g_Y = recHitGP.y();
    Float_t     rh_g_Z = recHitGP.z();
    int nVfat = 8*(findvfat(firstClusterStrip+clusterSize*0.5, 0, 128*3)-1) + (8-rh_roll);
    vMul[index][nVfat] += 1;
    gem_chamber_x_y[index]->Fill(rh_l_x, rh_roll);
    gem_chamber_cl_size[index]->Fill(clusterSize, nVfat);
    gem_chamber_bx[index]->Fill(bx,rh_roll);
    gemcr_g->Fill(-rh_g_X,rh_g_Y,rh_g_Z);
    gem_cls_tot->Fill(clusterSize);
    gem_bx_tot->Fill(bx);
    rh1_chamber->Fill(index);
    for(int i = firstClusterStrip; i < (firstClusterStrip + clusterSize); i++){
      gem_chamber_digi_recHit[index]->Fill(i,rh_roll);
    }
    strListRecHit += TString::Format("recHit : %lf, %lf, %lf (%i, %i) (%s)\n", rh_g_X, rh_g_Z, rh_g_Y, index + 1 - ( index % 2 ), ( index % 2 ) + 1, ( clusterSize >= minCLS && clusterSize <= maxCLS ? "okay" : "cut" ));
    if (clusterSize < minCLS) continue;
    if (clusterSize > maxCLS) continue;
    rh2_chamber->Fill(index);
    if ( rMul[ index ] <= 1 ) gem_chamber_track[ index ]->Fill(3.5);
    for(int i = firstClusterStrip; i < (firstClusterStrip + clusterSize); i++){
      sMul[index][rh_roll] +=1;
      sMul[index][0] +=1;
      gem_chamber_digi_CLS[index]->Fill(i,rh_roll);
    }
    
    float fPGenRecX = fXGenGP1x + ( rh_g_Z - fXGenGP1z ) * genMuon->momentum().x() / genMuon->momentum().z();
    float fPGenRecY = fXGenGP1z + ( rh_g_Z - fXGenGP1z ) * genMuon->momentum().y() / genMuon->momentum().z();
    gemcrGen_g->Fill(-fPGenRecX, fPGenRecY, rh_g_Z);
  }
  
  for ( int ich = 0 ; ich < n_ch ; ich++ ) {
    for ( int ivfat = 0 ; ivfat < 24 ; ivfat++ ) {
      if ( vMul[ ich ][ ivfat ] > 0 ) {
        int nRoll = 8 - ivfat % 8;
        int nVFat = ivfat / 8 + 1;
        
        gem_chamber_pad_vfat[ ich ]->Fill(nVFat, nRoll);
        for ( int mulidx = 0 ; mulidx < vMul[ ich ][ ivfat ] ; mulidx++ ) gem_chamber_pad_vfat_withmul[ ich ]->Fill(nVFat, nRoll);
        
        if ( vMul[ ich ^ 0x1 ][ ivfat ] > 0 ) {
          gem_chamber_copad_vfat[ ich ]->Fill(nVFat, nRoll);
          for ( int mulidx = 0 ; mulidx < vMul[ ich ^ 0x1 ][ ivfat ] ; mulidx++ ) gem_chamber_copad_vfat_withmul[ ich ]->Fill(nVFat, nRoll);
        }
      }
    }
  }
  
  int fChMul = 0;
  for(int c=0;c<n_ch;c++) {
    gem_chamber_hitMul[c]->Fill(rMul[c]);
    for(int v=0; v<24;v++){
      gem_chamber_vfatHitMul[c]->Fill(vMul[c][v],v);    
    }
    for(int r=0; r<9;r++) {
      gem_chamber_stripHitMul[c]->Fill(sMul[c][r],r);
    } 
    if (firedCh[c]) { 
      firedChamber->Fill(c+0.5);
      fChMul += 1;
    }
  }
  firedMul->Fill(fChMul);
    
  bool bIsScint = isPassedScintillators(genGPos1, genGPos2);
  
  /// Tracking start
  
  int fCha = 10;
  int lCha = 19;
  
  edm::Handle<std::vector<int>> idxChTraj;
  e.getByToken( this->InputTagToken_TI, idxChTraj);
  
  edm::Handle<std::vector<TrajectorySeed>> seedGCM;
  e.getByToken( this->InputTagToken_TS, seedGCM);
  
  edm::Handle<std::vector<Trajectory>> trajGCM;
  e.getByToken( this->InputTagToken_TJ, trajGCM);
  
  edm::Handle<vector<reco::Track>> trackCollection;
  e.getByToken( this->InputTagToken_TR, trackCollection);
  
  edm::Handle<std::vector<unsigned int>> seedTypes;
  e.getByToken( this->InputTagToken_TT, seedTypes);

  //created by Jongseok Lee ==>  
  //ntraj = idxChTraj->size();
  genTree->Fill();
/*
  float bestNChi2 = 10000;

  int nTRAJHIT = 0;
  const int maxNtraj = 30;
  int Bidx[maxNtraj];
  int Bivfat[maxNtraj];
  int BimRoll[maxNtraj];
  for(int i=0;i<maxNtraj;i++)
  {
    Bidx[i] = -1;
    Bivfat[i] = -1;
    BimRoll[i] = -1;
  }


  float BgenHitX[maxNtraj] = {0,};
  float BgenHitY[maxNtraj] = {0,};
  float BgenHitZ[maxNtraj] = {0,};
  float BtrajHitX[maxNtraj] = {0,};
  float BtrajHitY[maxNtraj] = {0,};
  float BtrajHitZ[maxNtraj] = {0,};
  float BrecHitX[maxNtraj] = {0,};
  float BrecHitY[maxNtraj] = {0,};
  float BrecHitZ[maxNtraj] = {0,};
  float BfloorHitX[maxNtraj] = {0,};
  float BfloorHitY[maxNtraj] = {0,};
  float BfloorHitZ[maxNtraj] = {0,};
  
  int ntch = 0;
  int bestTchIdx = -1;
  MuonTransientTrackingRecHit::MuonRecHitContainer testRecHits;
  //Finding a best trajectory per event
  for (auto tch : gemChambers)
  {
    ntch += 1;
    for (auto etaPart : tch.etaPartitions()){
      GEMDetId etaPartID = etaPart->id();
      GEMRecHitCollection::range range = gemRecHits->get(etaPartID);
      for (GEMRecHitCollection::const_iterator rechit = range.first; rechit!=range.second; ++rechit){
          const GeomDet* geomDet(etaPart);
          if ((*rechit).clusterSize()<minCLS) continue;
          if ((*rechit).clusterSize()>maxCLS) continue;
          testRecHits.push_back(MuonTransientTrackingRecHit::specificBuild(geomDet,&*rechit));        
      }
    }

    std::vector<int>::const_iterator it1 = idxChTraj->begin();
    std::vector<Trajectory>::const_iterator it3 = trajGCM->begin();
    
    TrajectorySeed bestSeed;
    Trajectory bestTraj;
    
    for ( ; it1 != idxChTraj->end() ; )
    {
      if ( *it1 == ntch )
      {
        bestTraj = *it3;
        bestSeed = (*it3).seed();
        break;
      }
      it1++;
      it3++;
    }
    
    if ( it1 == idxChTraj->end() ) continue;

    float NormalizedChi2 = bestTraj.chiSquared()/float(bestTraj.ndof());
    //cout<<"nev "<<nev<<", ntch "<<ntch-1<<", NormChi2 "<<NormalizedChi2<<endl;
    if(fabs(bestNChi2-1) > fabs(NormalizedChi2-1)) 
    {
      bestNChi2 = NormalizedChi2;
      bestTchIdx = ntch-1;
    }
  }
  //cout<<"nev "<<nev<<", bestTchIdx "<<bestTchIdx<<endl;

  if(bestTchIdx==-1) return;

  const int MaxNSeed = 2;
  int chSeedIdx[MaxNSeed] = {0,};
*/
  //created by Jongseok Lee <==  

  if ( idxChTraj->size() == 0 ) return;

  if (!makeTrack) return; 
  int countTC = 0;
  int nIsTraceGet = 0;
  int nIsLongSeed = 0;
  if ( fChMul == 4 ) events_withtraj->Fill(0.5);
  if ( fChMul == 5 ) events_withtraj->Fill(1.5);
  if ( fChMul >= 6 ) events_withtraj->Fill(2.5);

  //modified by Jongseok Lee ==>
  //GEMChamber tch = gemChambers[bestTchIdx];
  //countTC = bestTchIdx+1;
  for (auto tch : gemChambers)
  {
    countTC += 1;
    MuonTransientTrackingRecHit::MuonRecHitContainer testRecHits;
    for (auto etaPart : tch.etaPartitions())
    {
      GEMDetId etaPartID = etaPart->id();
      GEMRecHitCollection::range range = gemRecHits->get(etaPartID);
      for (GEMRecHitCollection::const_iterator rechit = range.first; rechit!=range.second; ++rechit)
      {
        const GeomDet* geomDet(etaPart);
        if ((*rechit).clusterSize()<minCLS) continue;
        if ((*rechit).clusterSize()>maxCLS) continue;
        testRecHits.push_back(MuonTransientTrackingRecHit::specificBuild(geomDet,&*rechit));        

        //GlobalPoint recHitGP = GEMGeometry_->idToDet((*rechit).gemId())->surface().toGlobal(rechit->localPosition());
        //grecHitX[nrecHit] = recHitGP.x();
        //grecHitY[nrecHit] = recHitGP.y();
        //grecHitZ[nrecHit] = recHitGP.z();
        //nrecHit++;
      }
    }
    //modified by Jongseok Lee <==

    std::vector<int>::const_iterator it1 = idxChTraj->begin();
    std::vector<TrajectorySeed>::const_iterator it2 = seedGCM->begin();
    std::vector<Trajectory>::const_iterator it3 = trajGCM->begin();
    std::vector<reco::Track>::const_iterator it4 = trackCollection->begin();
    std::vector<unsigned int>::const_iterator it5 = seedTypes->begin();
    
    TrajectorySeed bestSeed;
    Trajectory bestTraj;
    reco::Track bestTrack;
    unsigned int unTypeSeed = 0;
    
    for ( ; it1 != idxChTraj->end() ; )
    {
      if ( *it1 == countTC )
      {
        bestTraj = *it3;
        bestSeed = (*it3).seed();
        bestTrack = *it4;
        unTypeSeed = *it5;
        break;
      }
      
      it1++;
      it2++;
      it3++;
      it4++;
      it5++;
    }
    
    if ( it1 == idxChTraj->end() ) continue;
    
    const FreeTrajectoryState* ftsAtVtx = bestTraj.geometricalInnermostState().freeState();
    
    GlobalPoint trackPCA = ftsAtVtx->position();
    GlobalVector gvecTrack = ftsAtVtx->momentum();
    
    Float_t fTrackVelX = gvecTrack.x() / gvecTrack.z();
    Float_t fTrackVelY = gvecTrack.y() / gvecTrack.z();
    
    Float_t fSeedP1x = 0.0, fSeedP1y = 0.0, fSeedP1z = 0.0;
    Float_t fSeedP2x = 0.0, fSeedP2y = 0.0, fSeedP2z = 0.0;
    
    PTrajectoryStateOnDet ptsd1(bestSeed.startingState());
    DetId did(ptsd1.detId());
    const BoundPlane& bp = theService->trackingGeometry()->idToDet(did)->surface();
    TrajectoryStateOnSurface tsos = trajectoryStateTransform::transientState(ptsd1,&bp,&*theService->magneticField());
    TrajectoryStateOnSurface tsosCurrent = tsos;
    
    //created by Jongseok Lee ==>
    //int nch=0;
    //TrajectorySeed::range range = bestSeed.recHits();
    //for (OwnVector<TrackingRecHit>::const_iterator rechit = range.first; rechit!=range.second; ++rechit){
    //  GEMDetId hitID(rechit->rawId());
    //  chSeedIdx[nch] = hitID.chamber()+hitID.layer()-2;
    //  //cout<<"nev "<<nev<<", chSeedIdx "<<chSeedIdx[nch]<<endl;
    //  nch++;
    //}

    nTraj++;
    //float NormalizedChi2 = bestTraj.chiSquared()/float(bestTraj.ndof());
    //if(fabs(bestNChi2-1) <= fabs(NormalizedChi2-1)) continue;
    //bestNChi2 = NormalizedChi2;
    trajTheta = gvecTrack.theta();
    trajPhi = gvecTrack.phi();
    trajX = trackPCA.x();
    trajY = trackPCA.y();
    trajZ = trackPCA.z();
    trajPx = gvecTrack.x();
    trajPy = gvecTrack.y();
    trajPz = gvecTrack.z();
    //ntrajHit = 0;
    //for(int i=0;i<maxNtraj;i++)
    //{
    //  Bidx[i] = -1;
    //  Bivfat[i] = -1;
    //  BimRoll[i] = -1;
    //}
//cout<<"nev "<<nev<<", nTraj "<<nTraj<<", Theta "<<trajTheta<<", Phi "<<trajPhi<<", X "<<trajX<<", Y "<<trajY<<", Z "<<trajZ<<", Px "<<trajPx<<", Py "<<trajPy<<", Pz "<<trajPz<<endl;
//cout<<"nev "<<nev<<", nTraj "<<nTraj<<", Theta "<<trajTheta<<", Phi "<<trajPhi<<", X "<<trajX<<", Y "<<trajY<<", Z "<<trajZ<<", NormChi2 "<<NormalizedChi2<<endl;
//cout<<"nev "<<nev<<", nTraj "<<nTraj<<", Theta "<<trajTheta<<", Phi "<<trajPhi<<", X "<<trajX<<", Y "<<trajY<<", Z "<<trajZ<<", NormChi2 "<<NormalizedChi2<<endl;
    //created by Jongseok Lee <==

    int nTrajHit = 0, nTrajRecHit = 0, nTestHit = 0;
    //nTrajHit = number of TrajHit in a chamber(tch)
    //nTrajRecHit = number of RecHit in a chamber(tch)
    for(int c=0; c<n_ch;c++)
    {
      GEMChamber ch = gemChambers[c];
      const BoundPlane& bpch = GEMGeometry_->idToDet(ch.id())->surface();
      tsosCurrent = theService->propagator("SteppingHelixPropagatorAny")->propagate(tsosCurrent, bpch);
      if (!tsosCurrent.isValid()) continue;
      Global3DPoint gtrp = tsosCurrent.freeTrajectoryState()->position();
      Global3DPoint gtrpGEN(
        fXGenGP1x + ( genMuon->momentum().x() / genMuon->momentum().z() ) * ( gtrp.z() - fXGenGP1z ), 
        fXGenGP1y + ( genMuon->momentum().y() / genMuon->momentum().z() ) * ( gtrp.z() - fXGenGP1z ),
        gtrp.z()); 
      Local3DPoint tlp = bpch.toLocal(gtrp);
      if ( c == 10 ) {fSeedP1x = gtrp.x(); fSeedP1y = gtrp.y(); fSeedP1z = gtrp.z();}
      if ( c == 19 ) {fSeedP2x = gtrp.x(); fSeedP2y = gtrp.y(); fSeedP2z = gtrp.z();}
      Global3DPoint gtrp2(trackPCA.x() + fTrackVelX * ( gtrp.z() - trackPCA.z() ), trackPCA.y() + fTrackVelY * ( gtrp.z() - trackPCA.z() ), gtrp.z());
      if (!bpch.bounds().inside(tlp)){continue;}
      // modified by Jongseok Lee ==>
      if (ch==tch)
      {
      // modified by Jongseok Lee <==
        int n_roll = ch.nEtaPartitions();
        double minDely = 50.;
        int mRoll = -1;
        for (int r=0;r<n_roll;r++)
        {
          const BoundPlane& bproll = GEMGeometry_->idToDet(ch.etaPartition(r+1)->id())->surface();
          Local3DPoint rtlp = bproll.toLocal(gtrp);
          if(minDely > abs(rtlp.y())){minDely = abs(rtlp.y()); mRoll = r+1;}
        }

        if(1 == 0 && ( mRoll == 1 || mRoll == 8 )){
          bool tester = 1;
          for (int chId = fCha; chId < lCha+1; chId++){
            if (chId == findIndex(ch.id())) continue;
            if (!firedCh[chId]) tester = 0;
          }
          if (!tester) continue;
        }

        if (mRoll == -1){cout << "no mRoll" << endl;continue;}
        
        if ( mRoll == 1 && ( countTC == 11 || countTC == 20 ) && 
          ( unTypeSeed & QC8FLAG_SEEDINFO_MASK_REFVERTROLL18 ) != 0 )
        {
          projtheta_dist_edge_sim->Fill(180.0 / 3.141592 * atan2(genMuon->momentum().z(), genMuon->momentum().x()));
        }
        
        int n_strip = ch.etaPartition(mRoll)->nstrips();
        double min_x = ch.etaPartition(mRoll)->centreOfStrip(1).x();
        double max_x = ch.etaPartition(mRoll)->centreOfStrip(n_strip).x();
        
        if ( (tlp.x()>(min_x)) & (tlp.x() < (max_x)) )
        {

          // For testing the edge eta partition on the top and bottom layers only vertical seeds are allowed!
          if ( ( vecChamType[ countTC - 1 ] == 2 || vecChamType[ countTC - 1 ] == 1 ) && 
               ( mRoll == 1 || mRoll == 8 ) && 
               ( unTypeSeed & QC8FLAG_SEEDINFO_MASK_REFVERTROLL18 ) == 0 ) continue;
            
          uint32_t topOrBottomEtaIs8 = (unTypeSeed & QC8FLAG_SEEDINFO_MASK_TOPBOTTOMETA8) >> QC8FLAG_SEEDINFO_SHIFT_TOPBOTTOMETA8;       

          if (vecChamType[ countTC - 1 ] == 0 &&
          	  SuperChamType[ int((countTC-1)/2) ] == "S" &&
          	  topOrBottomEtaIs8 != 0)
          {
          	continue;
          }
          
          uint32_t topAndBottomEtaIs7 = (unTypeSeed & QC8FLAG_SEEDINFO_MASK_TOPBOTTOMETA7) >> QC8FLAG_SEEDINFO_SHIFT_TOPBOTTOMETA7;
          
          if (vecChamType[ countTC - 1 ] == 0 &&
          	  SuperChamType[ int((countTC-1)/2) ] == "S" &&
          	  topAndBottomEtaIs7 != 0)
          {
          	continue;
          }

          uint32_t unDiffCol = ( unTypeSeed & QC8FLAG_SEEDINFO_MASK_DIFFCOL ) >> QC8FLAG_SEEDINFO_SHIFT_DIFFCOL;
            
          if ( ! ( (tlp.x()>(min_x + 1.5)) & (tlp.x() < (max_x - 1.5)) ) )
          {
            if ( unDiffCol == 1 ) 
            {
              continue;
            }
            else if ( ( vecChamType[ countTC - 1 ] == 2 || vecChamType[ countTC - 1 ] == 1 ) )
            {
              continue;
            }
          }
          
          gem_chamber_track[findIndex(ch.id())]->Fill(4.5);
          int index = findIndex(ch.id());
          double vfat = findvfat(tlp.x(), min_x, max_x);

          //created by Jongseok Lee ==>

          int idx = index;
          int ivfat = (int)vfat - 1;
          int imRoll = mRoll - 1;
          vfatI[idx][ivfat][imRoll]=1;
          vfatF[idx][ivfat][imRoll]=0;

          int n1 = imRoll;
          int n2 = 2-ivfat + int(2-idx/10)*3;
          int n3 = idx%10;
          hvfatHit_denominator->Fill(n1,n2,n3);

          genHitX[idx][ivfat][imRoll] = genMuX + (gtrp.z()-genMuZ)*(gen_px/gen_pz);
          genHitY[idx][ivfat][imRoll] = genMuY + (gtrp.z()-genMuZ)*(gen_py/gen_pz);
          genHitZ[idx][ivfat][imRoll] = gtrp.z();
          trajHitX[idx][ivfat][imRoll] = gtrp.x();
          trajHitY[idx][ivfat][imRoll] = gtrp.y();
          trajHitZ[idx][ivfat][imRoll] = gtrp.z();
          recHitX[idx][ivfat][imRoll] = 0;
          recHitY[idx][ivfat][imRoll] = 0;
          recHitZ[idx][ivfat][imRoll] = 0;

          int floor = idx%10;
          floorHitX[floor] = gtrp.x();
          floorHitY[floor] = gtrp.y();
          floorHitZ[floor] = gtrp.z();

/*
          Bidx[ntrajHit] = index;
          Bivfat[ntrajHit] = (int)vfat - 1;
          BimRoll[ntrajHit] = mRoll - 1;
          BgenHitX[ntrajHit] = genMuX + (gtrp.y()-genMuY)*(gen_px/gen_py);
          BgenHitY[ntrajHit] = gtrp.y();
          BgenHitZ[ntrajHit] = genMuZ + (gtrp.y()-genMuY)*(gen_pz/gen_py);
          BtrajHitX[ntrajHit] = gtrp.x();
          BtrajHitY[ntrajHit] = gtrp.y();
          BtrajHitZ[ntrajHit] = gtrp.z();
          BfloorHitX[ntrajHit] = gtrp.x();
          BfloorHitY[ntrajHit] = gtrp.y();
          BfloorHitZ[ntrajHit] = gtrp.z();
          BrecHitX[ntrajHit] = 0;
          BrecHitY[ntrajHit] = 0;
          BrecHitZ[ntrajHit] = 0;
          nTRAJHIT++;
*/
          ntrajHit++;
          //created by Jongseok Lee <==

          gem_chamber_th2D_eff[index]->Fill(vfat, mRoll);                
          gem_chamber_thxroll_eff[index]->Fill(tlp.x(), mRoll);
          gem_chamber_thxy_eff[index]->Fill(tlp.x(), gtrp.y());
          gemcrTr_g->Fill(-gtrp.x(), gtrp.y(), gtrp.z());
          g_nNumTrajHit++;
          if ( nNumCurrFiredCh == 6 ) g_nNumTrajHit6++;
          
          if ( g_dMinX > gtrp.x() ) g_dMinX = gtrp.x();
          if ( g_dMaxX < gtrp.x() ) g_dMaxX = gtrp.x();
          if ( g_dMinY > gtrp.y() ) g_dMinY = gtrp.y();
          if ( g_dMaxY < gtrp.y() ) g_dMaxY = gtrp.y();
          if ( g_dMinZ > gtrp.z() ) g_dMinZ = gtrp.z();
          if ( g_dMaxZ < gtrp.z() ) g_dMaxZ = gtrp.z();
          
          nTrajHit++;
          
          if ( bIsScint )
          {
            gem_chamber_th2D_eff_scint[index]->Fill(vfat, mRoll);                
            gem_chamber_thxroll_eff_scint[index]->Fill(tlp.x(), mRoll);
            gem_chamber_thxy_eff_scint[index]->Fill(tlp.x(), gtrp.y());
            gemcrTrScint_g->Fill(-gtrp.x(), gtrp.y(), gtrp.z());
          }
          
          double maxR = 99.9;
          shared_ptr<MuonTransientTrackingRecHit> tmpRecHit;
          
          for (auto hit : testRecHits)
          {
            GEMDetId hitID(hit->rawId());
            if (hitID.chamberId() != ch.id()) continue;
            GlobalPoint hitGP = hit->globalPosition();
            if (abs(hitGP.x() - gtrp.x()) > maxRes) continue;
            if (abs(hitID.roll() - mRoll)>1) continue;
            //double diffXZ = sqrt((hitGP.x() - gtrp.x())*(hitGP.x() - gtrp.x())+(hitGP.z() - gtrp.z())*(hitGP.z() - gtrp.z()));
            //if(diffXZ>20) continue;
            double deltaR = (hitGP - gtrp).mag();
            if (maxR > deltaR)
            {
              tmpRecHit = hit;
              maxR = deltaR;
            }
          }
          
          if ( !tmpRecHit && testRecHits.size() > 0 ) gem_chamber_track[findIndex(ch.id())]->Fill(6.5);

          if(tmpRecHit){
            gem_chamber_track[findIndex(ch.id())]->Fill(5.5);
            Local3DPoint hitLP = tmpRecHit->localPosition();
            Global3DPoint recHitGP = tmpRecHit->globalPosition();
            
            gemcrCf_g->Fill(-recHitGP.x(), recHitGP.y(), recHitGP.z());
          
            nTrajRecHit++;

            //created by Jongseok Lee ==>
            recHitX[idx][ivfat][imRoll] = recHitGP.x();
            recHitY[idx][ivfat][imRoll] = recHitGP.y();
            recHitZ[idx][ivfat][imRoll] = recHitGP.z();
            vfatF[idx][ivfat][imRoll]=1;
            hvfatHit_numerator->Fill(n1,n2,n3);
            
            //BrecHitX[ntrajHit-1] = recHitGP.x();
            //BrecHitY[ntrajHit-1] = recHitGP.y();
            //BrecHitZ[ntrajHit-1] = recHitGP.z();
            ntrajRecHit++;
            //created by Jongseok Lee <==
            
            gem_chamber_tr2D_eff[index]->Fill(vfat, mRoll);
            gem_chamber_trxroll_eff[index]->Fill(tlp.x(), mRoll);
            gem_chamber_trxy_eff[index]->Fill(tlp.x(), gtrp.y());
            gem_chamber_local_x[index]->Fill(tlp.x(), hitLP.x());
            gem_chamber_residual[index]->Fill(tlp.x(), hitLP.x() - tlp.x());
            rh3_chamber->Fill(index);
            
            if ( bIsScint )
            {
              gemcrCfScint_g->Fill(-recHitGP.x(), recHitGP.y(), recHitGP.z());
              gem_chamber_tr2D_eff_scint[index]->Fill(vfat, mRoll);
              gem_chamber_trxroll_eff_scint[index]->Fill(tlp.x(), mRoll);
              gem_chamber_trxy_eff_scint[index]->Fill(tlp.x(), gtrp.y());
              gem_chamber_local_x_scint[index]->Fill(tlp.x(), hitLP.x());
              gem_chamber_residual_scint[index]->Fill(tlp.x(), hitLP.x() - tlp.x());
              rh3_chamber_scint->Fill(index);
            }
            
            g_nNumMatched++;
            
          } else {
            
            if ( countTC == 17 && mRoll == 6 && vfat == 1 ) {
              Float_t fVecX, fVecY;
              double dUnitGen = 0.1;
              
              fVecX = genMuon->momentum().x() / genMuon->momentum().z();
              fVecY = genMuon->momentum().y() / genMuon->momentum().z();
              
              fXGenGP1x = dUnitGen * genMuon->production_vertex()->position().x();
              fXGenGP1y = dUnitGen * genMuon->production_vertex()->position().y();
              fXGenGP1z = dUnitGen * genMuon->production_vertex()->position().z();
              
              Float_t fDiffZ = gtrp.z() - fXGenGP1z;
              Float_t fXGenHitX = fXGenGP1x + fDiffZ * fVecX;
              Float_t fXGenHitY = fXGenGP1y + fDiffZ * fVecY;
              
              printf("17_2_roll6_VFAT1 : event no. = %i ; # = %i\n", g_nEvt, (int)testRecHits.size());
              printf("GEN velocity : (%lf, %lf, 1.0)\n", fVecX, fVecY);
              printf(strKeep.Data());
              
              printf("reco trj hit : GEN (%lf, %lf, %lf) vs RECO_TRAJ (%lf, %lf, %lf)\n", 
                fXGenHitX, fXGenHitY, gtrp.z(), gtrp.x(), gtrp.y(), gtrp.z());
              
              for (GEMRecHitCollection::const_iterator hit = gemRecHits->begin(); hit != gemRecHits->end(); ++hit)
              {
                GlobalPoint hitGP = GEMGeometry_->idToDet((*hit).gemId())->surface().toGlobal(hit->localPosition());
                printf("    associated recHits : (%lf, %lf, %lf) ; (%lf, %lf, %lf)\n", 
                  hitGP.x(), hitGP.y(), hitGP.z(), 
                  trackPCA.x() + fTrackVelX * ( hitGP.z() - trackPCA.z() ), 
                  trackPCA.y() + fTrackVelY * ( hitGP.z() - trackPCA.z() ),
                  hitGP.z()); 
              }
            }
          }
        }
        continue;
      }
    }
//cout<<"nev "<<nev<<", nTraj "<<nTraj<<", Theta "<<trajTheta<<", Phi "<<trajPhi<<", X "<<trajX<<", Y "<<trajY<<", Z "<<trajZ<<", NormChi2 "<<NormalizedChi2<<", nTRAJHIT "<<nTRAJHIT<<", nTrajHit "<<nTrajHit<<", nTrajRecHit "<<nTrajRecHit<<", eff "<<double(nTrajRecHit)/nTrajHit;
//if(nTrajHit>0) 
//{
//  cout<<", vfatID ";
//  for(int nth=0;nth<nTrajHit;nth++) cout<<Bidx[nth]<<" "<<Bivfat[nth]<<" "<<BimRoll[nth]<<", ";
//}
//cout<<endl;
    
    if ( 11 <= countTC && countTC <= 20 )
    {
      Float_t fSeedDiffZ = fSeedP2z - fSeedP1z;
      Float_t fSeedVelX = ( fSeedP2x - fSeedP1x ) / fSeedDiffZ;
      Float_t fSeedVelY = ( fSeedP2y - fSeedP1y ) / fSeedDiffZ;
      
      Float_t fSeedVtxX = fSeedP1x + fSeedVelX * ( fXGenGP1z - fSeedP1z );
      Float_t fSeedVtxY = fSeedP1y + fSeedVelY * ( fXGenGP1z - fSeedP1z );
      
      Float_t fTrackVtxX = trackPCA.x() + fTrackVelX * ( fXGenGP1z - trackPCA.z() );
      Float_t fTrackVtxY = trackPCA.y() + fTrackVelY * ( fXGenGP1z - trackPCA.z() );
      
      g_resXRTSim->Fill(fTrackVtxX - fSeedVtxX);
      g_resYRTByErrSim->Fill(fTrackVtxY - fSeedVtxY);
    }
    
    if ( nTestHit != 0 && 1 == 0 )
    {
      printf("### missing hit occurs! ###\n");
    }
  }
  if ( nIsLongSeed != 0 ) events_withtraj->Fill(3.5);
  if ( nIsTraceGet != 0 )
  {
    if ( fChMul == 4 ) events_withtraj->Fill(4.5);
    if ( fChMul == 5 ) events_withtraj->Fill(5.5);
    if ( fChMul >= 6 ) events_withtraj->Fill(6.5);
    events_withtraj->Fill(7.5);
  }
  
  g_nNumTest++;

  //created by Jongseok Lee ==>
  //cout<<"nev "<<nev<<", nTraj "<<nTraj<<endl;

/*
  for(int i=0;i<ntrajHit;i++)
  {
    if(Bidx[i]!=-1 && Bivfat[i]!=-1 && BimRoll[i]!=-1 && Bidx[i]!=chSeedIdx[0] && Bidx[i]!=chSeedIdx[1])
    {
      int n1 = BimRoll[i];
      int n2 = 2-Bivfat[i] + int(2-Bidx[i]/10)*3;
      int n3 = Bidx[i]%10;
      hvfatHit_denominator->Fill(n1,n2,n3);
      vfatI[Bidx[i]][Bivfat[i]][BimRoll[i]] = 1;
      if(BrecHitX[i]!=0 && BrecHitY[i]!=0 && BrecHitZ[i]!=0)
      {
        hvfatHit_numerator->Fill(n1,n2,n3);
        vfatF[Bidx[i]][Bivfat[i]][BimRoll[i]] = 1;
      }
      genHitX[Bidx[i]][Bivfat[i]][BimRoll[i]] = BgenHitX[i];
      genHitY[Bidx[i]][Bivfat[i]][BimRoll[i]] = BgenHitY[i];
      genHitZ[Bidx[i]][Bivfat[i]][BimRoll[i]] = BgenHitZ[i];
      trajHitX[Bidx[i]][Bivfat[i]][BimRoll[i]] = BtrajHitX[i];
      trajHitY[Bidx[i]][Bivfat[i]][BimRoll[i]] = BtrajHitY[i];
      trajHitZ[Bidx[i]][Bivfat[i]][BimRoll[i]] = BtrajHitZ[i];
      floorHitX[Bidx[i]%10] = BfloorHitX[i];
      floorHitY[Bidx[i]%10] = BfloorHitY[i];
      floorHitZ[Bidx[i]%10] = BfloorHitZ[i];
      recHitX[Bidx[i]][Bivfat[i]][BimRoll[i]] = BrecHitX[i];
      recHitY[Bidx[i]][Bivfat[i]][BimRoll[i]] = BrecHitY[i];
      recHitZ[Bidx[i]][Bivfat[i]][BimRoll[i]] = BrecHitZ[i];
    }
  }
*/
  tree->Fill();
  hev->Fill(2);
  //if(ntrajHit>=2) 
  //{
  //  tree->Fill();
  //  hev->Fill(2);
  //}
  //created by Jongseok Lee <==

}
