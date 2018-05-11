// class GEMCosmicMuonForQC8

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/GEMGeometry/interface/GEMGeometry.h>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/CosmicMuonProducer/interface/CosmicMuonSmoother.h"
#include "RecoMuon/StandAloneTrackFinder/interface/StandAloneMuonSmoother.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include <iostream>
#include <iomanip>
#include <string.h>
#include <vector>
#include "TH1.h"
#include "TH2D.h"
#include "TH2.h"
#include "TH3F.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMath.h"

using namespace std;

class GEMCosmicMuonForQC8 : public edm::stream::EDProducer<> {
public:
  /// Constructor
  explicit GEMCosmicMuonForQC8(const edm::ParameterSet&);
  /// Destructor
  virtual ~GEMCosmicMuonForQC8() {}
  /// Produce the GEMSegment collection
  void produce(edm::Event&, const edm::EventSetup&) override;
  std::vector<std::string> g_SuperChamType;
  vector<double> g_vecChamType;
private:
  int iev; // events through
  edm::EDGetTokenT<GEMRecHitCollection> theGEMRecHitToken;
  CosmicMuonSmoother* theSmoother;
  MuonServiceProxy* theService;
  KFUpdator* theUpdator;
  int findSeeds(std::vector<TrajectorySeed> *tmptrajectorySeeds, MuonTransientTrackingRecHit::MuonRecHitContainer &muRecHits);
  Trajectory makeTrajectory(TrajectorySeed seed, MuonTransientTrackingRecHit::MuonRecHitContainer &muRecHits, std::vector<GEMChamber> gemChambers);
  const GEMGeometry* gemGeom;
};

GEMCosmicMuonForQC8::GEMCosmicMuonForQC8(const edm::ParameterSet& ps) : iev(0) {
  time_t rawTime;
  time(&rawTime);
  printf("Begin of GEMCosmicMuonForQC8::GEMCosmicMuonForQC8() at %s\n", asctime(localtime(&rawTime)));
  g_SuperChamType = ps.getParameter<vector<string>>("SuperChamberType");
  g_vecChamType = ps.getParameter<vector<double>>("SuperChamberSeedingLayers");
  theGEMRecHitToken = consumes<GEMRecHitCollection>(ps.getParameter<edm::InputTag>("gemRecHitLabel"));
  edm::ParameterSet serviceParameters = ps.getParameter<edm::ParameterSet>("ServiceParameters");
  theService = new MuonServiceProxy(serviceParameters);
  edm::ParameterSet smootherPSet = ps.getParameter<edm::ParameterSet>("MuonSmootherParameters");
  theSmoother = new CosmicMuonSmoother(smootherPSet,theService);
  theUpdator = new KFUpdator();
  printf("End of GEMCosmicMuonForQC8::GEMCosmicMuonForQC8() at %s\n", asctime(localtime(&rawTime)));
}

void GEMCosmicMuonForQC8::produce(edm::Event& ev, const edm::EventSetup& setup)
{
  theService->update(setup);

  edm::ESHandle<GEMGeometry> gemg;
  setup.get<MuonGeometryRecord>().get(gemg);
  const GEMGeometry* mgeom = &*gemg;
  gemGeom = &*gemg;

  vector<GEMChamber> gemChambers;

  const std::vector<const GEMSuperChamber*>& superChambers_ = mgeom->superChambers();
  for (auto sch : superChambers_)
  {
    int n_lay = sch->nChambers();
    for (int l=0;l<n_lay;l++)
    {
      gemChambers.push_back(*sch->chamber(l+1));
    }
  }

  edm::Handle<GEMRecHitCollection> gemRecHits;
  ev.getByToken(theGEMRecHitToken,gemRecHits);

  if (gemRecHits->size() < 4) return;

  MuonTransientTrackingRecHit::MuonRecHitContainer muRecHits;
  MuonTransientTrackingRecHit::MuonRecHitContainer seedRecHits;

  for (auto ch : gemChambers)
  {
    for (auto etaPart : ch.etaPartitions())
    {
      GEMDetId etaPartID = etaPart->id();
      GEMRecHitCollection::range range = gemRecHits->get(etaPartID);

      for (GEMRecHitCollection::const_iterator rechit = range.first; rechit!=range.second; ++rechit)
      {
        const GeomDet* geomDet(etaPart);
        if ((*rechit).clusterSize()<1 || (*rechit).clusterSize()>10) continue;
        muRecHits.push_back(MuonTransientTrackingRecHit::specificBuild(geomDet,&*rechit));
      }
    }
  }

  if (muRecHits.size() < 4) return;

  vector<TrajectorySeed> trajSeedsBody;
  std::vector<TrajectorySeed> *trajSeeds = &trajSeedsBody;
  findSeeds(trajSeeds, muRecHits);
  Trajectory bestTrajectory;
  TrajectorySeed bestSeed;

  float maxChi2 = 10000000.0;

  for (auto seed : *trajSeeds)
  {
    Trajectory smoothed = makeTrajectory(seed, muRecHits, gemChambers);

    if (smoothed.isValid())
    {
      float dProbChiNDF = smoothed.chiSquared()/float(smoothed.ndof());

      if (fabs(maxChi2-1) > fabs(dProbChiNDF-1))
      {
        maxChi2 = dProbChiNDF;
        bestTrajectory = smoothed;
        bestSeed = seed;
      }
    }
  }

  if (!bestTrajectory.isValid()) return;
  if (maxChi2 > 3) return;

  /// Alignment section STARTS HERE!!!!!

  TH1F * resX[30][8];
  char * histname = new char[10];
  for (Int_t i = 0 ; i < 30 ; i++)
  {
    for (Int_t j = 0 ; j < 8 ; i++)
    {
      sprintf(histname, "h_resX_%d_%u",i,j);
      resX[i][j] = new TH1F(histname,"",100,-5,5);
    }
  }
  float dx = 0.0;
  float rz = 0.0;

  PTrajectoryStateOnDet ptsd(bestSeed.startingState());
  DetId did(ptsd.detId());
  const BoundPlane& bp = theService->trackingGeometry()->idToDet(did)->surface();
  TrajectoryStateOnSurface tsos = trajectoryStateTransform::transientState(ptsd,&bp,&*theService->magneticField());
  TrajectoryStateOnSurface tsosCurrent = tsos;

  for (auto ch : gemChambers)
  {
    tsosCurrent = theService->propagator("SteppingHelixPropagatorAny")->propagate(tsosCurrent, theService->trackingGeometry()->idToDet(ch.id())->surface());
    if (!tsosCurrent.isValid()) continue;
    GlobalPoint tsosGP = tsosCurrent.freeTrajectoryState()->position();

    for (auto hit : bestTrajectory.recHits())
    {
      GEMDetId hitID(hit->rawId());
      if (hitID.chamberId() == ch.id() )
      {
        GlobalPoint hitGP = hit->globalPosition();

        int chNumber = ch.id().chamber() + ch.id().layer() - 1;

        int n_roll = ch.nEtaPartitions();
        double minDely = 50.;
        int Roll = 0;
        for (int r = 1 ; r <= n_roll ; r++)
        {
          const BoundPlane& bproll = theService->trackingGeometry()->idToDet(ch.etaPartition(r)->id())->surface();
          Local3DPoint rtlp = bproll.toLocal(hitGP);

          cout << Roll << " " << rtlp.y() << " " << ch.etaPartition(Roll)->centreOfStrip(1).y() << endl;

          if(minDely > abs(rtlp.y()))
          {
            minDely = abs(rtlp.y());
            Roll = r;
          }
        }
        if (Roll == 0) continue;

        resX[chNumber-1][Roll-1] -> Fill(hitGP.x() - tsosGP.x());
      }
    }
  }

  TCanvas *cnv = new TCanvas("cnv", "Plot",950,730);
  resX[0][0] -> Draw();
  cnv -> Update();
}


int GEMCosmicMuonForQC8::findSeeds(std::vector<TrajectorySeed> *tmptrajectorySeeds, MuonTransientTrackingRecHit::MuonRecHitContainer &muRecHits)
{
  for (auto hit1 : muRecHits)
  {
    for (auto hit2 : muRecHits)
    {
      if (hit1->globalPosition().z() < hit2->globalPosition().z())
      {
        LocalPoint segPos = hit1->localPosition();
        GlobalVector segDirGV(hit2->globalPosition().x() - hit1->globalPosition().x(),
                              hit2->globalPosition().y() - hit1->globalPosition().y(),
                              hit2->globalPosition().z() - hit1->globalPosition().z());

        segDirGV *=10;
        LocalVector segDir = hit1->det()->toLocal(segDirGV);

        int charge= 1;
        LocalTrajectoryParameters param(segPos, segDir, charge);

        AlgebraicSymMatrix mat(5,0);
        mat = hit1->parametersError().similarityT( hit1->projectionMatrix() );
        LocalTrajectoryError error(asSMatrix<5>(mat));

        TrajectoryStateOnSurface tsos(param, error, hit1->det()->surface(), &*theService->magneticField());
        uint32_t id = hit1->rawId();
        PTrajectoryStateOnDet const & seedTSOS = trajectoryStateTransform::persistentState(tsos, id);

        edm::OwnVector<TrackingRecHit> seedHits;
        seedHits.push_back(hit1->hit()->clone());
        seedHits.push_back(hit2->hit()->clone());

        TrajectorySeed seed(seedTSOS,seedHits,alongMomentum);
        tmptrajectorySeeds->push_back(seed);
      }
    }
  }
  return 0;
}


Trajectory GEMCosmicMuonForQC8::makeTrajectory(TrajectorySeed seed, MuonTransientTrackingRecHit::MuonRecHitContainer &muRecHits, std::vector<GEMChamber> gemChambers)
{
  PTrajectoryStateOnDet ptsd1(seed.startingState());
  DetId did(ptsd1.detId());
  const BoundPlane& bp = theService->trackingGeometry()->idToDet(did)->surface();
  TrajectoryStateOnSurface tsos = trajectoryStateTransform::transientState(ptsd1,&bp,&*theService->magneticField());
  TrajectoryStateOnSurface tsosCurrent = tsos;
  TransientTrackingRecHit::ConstRecHitContainer consRecHits;

  TrajectorySeed::range range = seed.recHits();
  int nseed = 0;
  GlobalPoint seedGP[2];
  for (edm::OwnVector<TrackingRecHit>::const_iterator rechit = range.first; rechit!=range.second; ++rechit)
  {
    GEMDetId hitID(rechit->rawId());
    seedGP[nseed] = gemGeom->idToDet((*rechit).rawId())->surface().toGlobal(rechit->localPosition());
    nseed++;
  }
  std::map<double,int> rAndhit;

  for (auto ch : gemChambers)
  {
    std::shared_ptr<MuonTransientTrackingRecHit> tmpRecHit;
    tsosCurrent = theService->propagator("SteppingHelixPropagatorAny")->propagate(tsosCurrent, theService->trackingGeometry()->idToDet(ch.id())->surface());
    if (!tsosCurrent.isValid()) return Trajectory();
    GlobalPoint tsosGP = tsosCurrent.freeTrajectoryState()->position();

    float maxR = 9999.;
    int nhit=-1;
    int tmpNhit=-1;
    double tmpR=-1.;
    for (auto hit : muRecHits)
    {
      nhit++;
      GEMDetId hitID(hit->rawId());
      if (hitID.chamberId() == ch.id() )
      {
        GlobalPoint hitGP = hit->globalPosition();
        if (fabs(hitGP.x() - tsosGP.x()) > 3) continue;
        if (fabs(hitGP.y() - tsosGP.y()) > 30) continue;
        float deltaR = (hitGP - tsosGP).mag();
        if (maxR > deltaR)
        {
          tmpRecHit = hit;
          maxR = deltaR;
          tmpNhit = nhit;
          tmpR = (hitGP - seedGP[0]).mag();
        }
      }
    }
    if (tmpRecHit)
    {
      rAndhit[tmpR] = tmpNhit;
    }
  }

  if (rAndhit.size() < 5) return Trajectory();
  vector<pair<double,int>> rAndhitV;
  copy(rAndhit.begin(), rAndhit.end(), back_inserter<vector<pair<double,int>>>(rAndhitV));
  for(unsigned int i=0;i<rAndhitV.size();i++)
  {
    consRecHits.push_back(muRecHits[rAndhitV[i].second]);
  }

  if (consRecHits.size() < 5) return Trajectory();
  vector<Trajectory> fitted = theSmoother->trajectories(seed, consRecHits, tsos);
  if ( fitted.size() <= 0 ) return Trajectory();

  Trajectory smoothed = fitted.front();
  return fitted.front();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GEMCosmicMuonForQC8);
