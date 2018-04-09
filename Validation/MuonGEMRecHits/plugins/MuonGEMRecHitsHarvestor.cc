// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"

///Data Format
#include <DataFormats/GEMRecHit/interface/GEMRecHit.h>
#include <DataFormats/GEMRecHit/interface/GEMRecHitCollection.h>
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/GeometrySurface/interface/LocalError.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

///Geometry
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "DQMServices/Core/interface/DQMStore.h"

///Log messages
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Validation/MuonGEMRecHits/plugins/MuonGEMRecHitsHarvestor.h"


MuonGEMRecHitsHarvestor::MuonGEMRecHitsHarvestor(const edm::ParameterSet& ps)
{
  dbe_path_ = std::string("MuonGEMRecHitsV/GEMRecHitsTask/");
  outputFile_ = ps.getUntrackedParameter<std::string>("outputFile", "myfile.root");
}


MuonGEMRecHitsHarvestor::~MuonGEMRecHitsHarvestor()
{
}
TProfile* MuonGEMRecHitsHarvestor::ComputeEff(TH1F* num, TH1F* denum )
{
  std::string name  = "eff_"+std::string(num->GetName());
  std::string title = "Eff. "+std::string(num->GetTitle());
  TProfile * efficHist = new TProfile(name.c_str(), title.c_str(),denum->GetXaxis()->GetNbins(), denum->GetXaxis()->GetXmin(),denum->GetXaxis()->GetXmax());

  for (int i=1; i <= denum->GetNbinsX(); i++) {

    double nNum = num->GetBinContent(i);
    double nDenum = denum->GetBinContent(i);
    if ( nDenum == 0 || nNum ==0  ) {
      continue;
    }
    if ( nNum > nDenum ) {
      double temp = nDenum;
      nDenum = nNum;
      nNum = temp;
      std::cout<<"Alert! specific bin's num is bigger than denum"<<std::endl;
    }
    const double effVal = nNum/nDenum;
    efficHist->SetBinContent(i, effVal);
    efficHist->SetBinEntries(i,1);
    efficHist->SetBinError(i,0);
    const double errLo = TEfficiency::ClopperPearson((int)nDenum,(int)nNum,0.683,false);
    const double errUp = TEfficiency::ClopperPearson((int)nDenum,(int)nNum,0.683,true);
    const double errVal = (effVal - errLo > errUp - effVal) ? effVal - errLo : errUp - effVal;
    efficHist->SetBinError(i, sqrt(effVal * effVal + errVal * errVal));
  }
  return efficHist;
}


TProfile* MuonGEMRecHitsHarvestor::ComputeEff2DForHitEff(TH2F* num, TH2F* denum )
{
  std::string name  = "eff_"+std::string(num->GetName());
  std::string title = "Eff. "+std::string(num->GetTitle());
  //TProfile * efficHist = new TProfile(name.c_str(), title.c_str(),denum->GetXaxis()->GetNbins(), denum->GetXaxis()->GetXmin(),denum->GetXaxis()->GetXmax());
  TProfile * efficHist = new TProfile(name.c_str(), title.c_str(), 24, 0 - 0.5, 24 - 0.5);

  for (int i=1; i <= denum->GetNbinsX(); i++) {
    for (int j=1; j <= denum->GetNbinsY(); j++) {

      int nIdx = 8 * ( i - 1 ) + j;
      
      double nNum = num->GetBinContent(i, j);
      double nDenum = denum->GetBinContent(i, j);
      if ( nDenum == 0 || nNum ==0  ) {
        continue;
      }
      if ( nNum > nDenum ) {
        double temp = nDenum;
        nDenum = nNum;
        nNum = temp;
        std::cout<<"Alert! specific bin's num is bigger than denum"<<std::endl;
      }
      const double effVal = nNum/nDenum;
      efficHist->SetBinContent(nIdx, effVal);
      efficHist->SetBinEntries(nIdx,1);
      efficHist->SetBinError(nIdx,0);
      const double errLo = TEfficiency::ClopperPearson((int)nDenum,(int)nNum,0.683,false);
      const double errUp = TEfficiency::ClopperPearson((int)nDenum,(int)nNum,0.683,true);
      const double errVal = (effVal - errLo > errUp - effVal) ? effVal - errLo : errUp - effVal;
      efficHist->SetBinError(nIdx, sqrt(effVal * effVal + errVal * errVal));
    }
  }
  return efficHist;
}


TH2F* MuonGEMRecHitsHarvestor::ComputeEff2DForCopad(TH2F* num, TH2F* denum )
{
  std::string name  = "eff_"+std::string(num->GetName());
  std::string title = "Eff. "+std::string(num->GetTitle());
  //TProfile * efficHist = new TProfile(name.c_str(), title.c_str(),denum->GetXaxis()->GetNbins(), denum->GetXaxis()->GetXmin(),denum->GetXaxis()->GetXmax());
  TH2F * efficHist = new TH2F(name.c_str(), title.c_str(), 3, 1, 4, 8, 1, 9);

  for (int i=1; i <= denum->GetNbinsX(); i++) {
    for (int j=1; j <= denum->GetNbinsY(); j++) {

      //int nIdx = 8 * ( i - 1 ) + j;
      
      double nNum = num->GetBinContent(i, j);
      double nDenum = denum->GetBinContent(i, j);
      if ( nDenum == 0 || nNum ==0  ) {
        continue;
      }
      if ( nNum > nDenum ) {
        double temp = nDenum;
        nDenum = nNum;
        nNum = temp;
        std::cout<<"Alert! specific bin's num is bigger than denum"<<std::endl;
      }
      const double effVal = nNum/nDenum;
      efficHist->SetBinContent(i, j, effVal);
      //efficHist->SetBinEntries(i, j,1);
      efficHist->SetBinError(i, j,0);
      const double errLo = TEfficiency::ClopperPearson((int)nDenum,(int)nNum,0.683,false);
      const double errUp = TEfficiency::ClopperPearson((int)nDenum,(int)nNum,0.683,true);
      const double errVal = (effVal - errLo > errUp - effVal) ? effVal - errLo : errUp - effVal;
      efficHist->SetBinError(i, j, sqrt(effVal * effVal + errVal * errVal));
    }
  }
  return efficHist;
}

void MuonGEMRecHitsHarvestor::ProcessBooking( DQMStore::IBooker& ibooker, DQMStore::IGetter& ig, const char* label, TString suffix, TH1F* track_hist, TH1F* sh_hist )
{
  TString dbe_label = TString(dbe_path_)+label+suffix;
  if( ig.get(dbe_label.Data()) != nullptr && sh_hist !=nullptr && track_hist !=nullptr ) {
    TH1F* hist =    (TH1F*)ig.get( dbe_label.Data() )->getTH1F()->Clone();
    TProfile* profile = ComputeEff( hist, track_hist);
    TProfile* profile_sh = ComputeEff( hist, sh_hist );
    profile_sh->SetName( (profile->GetName()+std::string("_sh")).c_str());
    TString x_axis_title = TString(hist->GetXaxis()->GetTitle());
    TString title  = TString::Format("Eff. for a SimTrack to have an associated GEM RecHits in %s;%s;Eff.",suffix.Data(),x_axis_title.Data());
    TString title2 = TString::Format("Eff. for a SimTrack to have an associated GEM RecHits in %s with a matched SimHit;%s;Eff.",suffix.Data(),x_axis_title.Data() );
    profile->SetTitle( title.Data());
    profile_sh->SetTitle( title2.Data() );
    ibooker.bookProfile( profile->GetName(),profile); 
    ibooker.bookProfile( profile_sh->GetName(),profile_sh); 
  }
  else {
    std::cout<<"Can not found histogram of "<<dbe_label<<std::endl; 
    if ( track_hist == nullptr) std::cout<<"track not found"<<std::endl;
    if ( sh_hist    == nullptr) std::cout<<"sh_hist not found"<<std::endl;
  }
  return;
}


void 
MuonGEMRecHitsHarvestor::dqmEndJob(DQMStore::IBooker& ibooker, DQMStore::IGetter& ig)
{
  ig.setCurrentFolder(dbe_path_.c_str());
 
  const char* l_suffix[4] = {"_l1","_l2","_l1or2","_l1and2"};
  const char* s_suffix[2] = {"_st1","_st2"};   
  const char* c_suffix[3] = {"_even","_odd","_all"};   

  TH1F* gem_trk_eta[2];
  TH1F* gem_trk_phi[2][2];  

  TH1F* sh_eta[2][4];
  TH1F* sh_phi[2][4][3];
  
  std::cout << dbe_path_ << std::endl;
  
  for( int i = 0 ; i < 2 ; i++) {
    TString eta_label = TString(dbe_path_)+"track_eta"+s_suffix[i];
    TString phi_label;
    if ( ig.get(eta_label.Data()) != nullptr ) {
      gem_trk_eta[i] = (TH1F*)ig.get(eta_label.Data())->getTH1F()->Clone();
      gem_trk_eta[i]->Sumw2();
    }
    else std::cout<<"Can not found track_eta"<<std::endl;
    for ( int k=0 ; k <3 ; k++) {
      phi_label = TString(dbe_path_.c_str())+"track_phi"+s_suffix[i]+c_suffix[k];
      if ( ig.get(phi_label.Data()) !=nullptr ) {
        gem_trk_phi[i][k] = (TH1F*)ig.get(phi_label.Data())->getTH1F()->Clone();
        gem_trk_phi[i][k]->Sumw2();
      }
      else std::cout<<"Can not found track_phi"<<std::endl;
    }
    
    if ( ig.get(eta_label.Data()) != nullptr && ig.get(phi_label.Data()) !=nullptr ) {
      for( int j = 0; j < 4 ; j++) { 
        TString suffix = TString( s_suffix[i] )+TString( l_suffix[j]);
        TString eta_label = TString(dbe_path_)+"rh_sh_eta"+suffix;
        if( ig.get(eta_label.Data()) !=nullptr ) {
          sh_eta[i][j] = (TH1F*)ig.get(eta_label.Data())->getTH1F()->Clone();
          sh_eta[i][j]->Sumw2();
        }
        else std::cout<<"Can not found eta histogram : "<<eta_label<<std::endl;
        ProcessBooking( ibooker, ig, "rh_eta", suffix, gem_trk_eta[i], sh_eta[i][j]); 
        for ( int k= 0 ; k< 3 ; k++) {
          suffix = TString( s_suffix[i])+TString( l_suffix[j]) +TString(c_suffix[k]);
          TString phi_label = TString(dbe_path_)+"rh_sh_phi"+suffix;
          if( ig.get(phi_label.Data()) !=nullptr ) {
           sh_phi[i][j][k] = (TH1F*)ig.get(phi_label.Data())->getTH1F()->Clone();
           sh_phi[i][j][k]->Sumw2();
          }
          else { std::cout<<"Can not found phi plots : "<<phi_label<<std::endl; continue; }
          ProcessBooking( ibooker, ig, "rh_phi",suffix, gem_trk_phi[i][k], sh_phi[i][j][k]);
        }
      }
    }
    else std::cout<<"Can not find eta or phi of all track"<<std::endl;
  }
  
  int i, j;
  
  // Making the efficiency plots of trajectory-hit
  for ( i = 1 ; i <= 29 ; i += 2 ) {
    for ( j = 1 ; j <= 2 ; j++ ) {
      TString strChamber = TString::Format("MuonGEMRecHitsV/GEMRecHitsTask/chamber_%i_layer_%i_", i, j);
      TString strNameDen = strChamber + "th2D_eff";
      TString strNameNum = strChamber + "recHit_efficiency";
      
      TString strNameEff = TString::Format("chamber_%i_layer_%i_efficiency", i, j);
      
      TH2F *hDen = (TH2F *)ig.get(strNameDen.Data())->getTH2F()->Clone();
      TH2F *hNum = (TH2F *)ig.get(strNameNum.Data())->getTH2F()->Clone();
      
      TProfile *profileTest = ComputeEff2DForHitEff(hNum, hDen);
      profileTest->SetName(strNameEff.Data());
      ibooker.bookProfile(profileTest->GetName(), profileTest);
    }
  }
  
  for ( i = 1 ; i <= 29 ; i += 2 ) {
    for ( j = 1 ; j <= 2 ; j++ ) {
      TString strChamber = TString::Format("MuonGEMRecHitsV/GEMRecHitsTask/chamber_%i_layer_%i_", i, j);
      TString strNameDen = strChamber + "copad_vfat";
      TString strNameNum = strChamber + "pad_vfat";
      
      TString strNameEff = TString::Format("chamber_%i_layer_%i_copad_efficiency", i, j);
      
      TH2F *hDen = (TH2F *)ig.get(strNameDen.Data())->getTH2F()->Clone();
      TH2F *hNum = (TH2F *)ig.get(strNameNum.Data())->getTH2F()->Clone();
      
      TH2F *profileTest = ComputeEff2DForCopad(hNum, hDen);
      profileTest->SetName(strNameEff.Data());
      ibooker.book2D(profileTest->GetName(), profileTest);
    }
  }
  
  for ( i = 1 ; i <= 29 ; i += 2 ) {
    for ( j = 1 ; j <= 2 ; j++ ) {
      TString strChamber = TString::Format("MuonGEMRecHitsV/GEMRecHitsTask/chamber_%i_layer_%i_", i, j);
      TString strNameDen = strChamber + "pad_vfat";
      TString strNameNum = strChamber + "recHit_efficiency";
      
      TString strNameEff = TString::Format("chamber_%i_layer_%i_matched_efficiency", i, j);
      
      TH2F *hDen = (TH2F *)ig.get(strNameDen.Data())->getTH2F()->Clone();
      TH2F *hNum = (TH2F *)ig.get(strNameNum.Data())->getTH2F()->Clone();
      
      TH2F *profileTest = ComputeEff2DForCopad(hNum, hDen);
      profileTest->SetName(strNameEff.Data());
      ibooker.book2D(profileTest->GetName(), profileTest);
    }
  }
  
  TH1F *histRH2 = ig.get("MuonGEMRecHitsV/GEMRecHitsTask/rh2_chamber")->getTH1F();
  TH1F *histRH3 = ig.get("MuonGEMRecHitsV/GEMRecHitsTask/rh3_chamber")->getTH1F();
  TH1F *histRHeff = ig.get("MuonGEMRecHitsV/GEMRecHitsTask/rh_eff_ch")->getTH1F();
  
  for ( i = 1 ; i <= 30 ; i++ ) {
    float fDen = histRH2->GetBinContent(i);
    float fNum = histRH3->GetBinContent(i);
    histRHeff->SetBinContent(i, fNum / fDen);
  }
}


//define this as a plug-in
DEFINE_FWK_MODULE(MuonGEMRecHitsHarvestor);
